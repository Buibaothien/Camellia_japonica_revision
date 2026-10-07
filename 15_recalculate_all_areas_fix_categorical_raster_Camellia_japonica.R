# =====================================================================
# 15_RECALCULATE_ALL_AREAS_FIX_CATEGORICAL_RASTER_CAMELLIA_JAPONICA.R
# =====================================================================
# PURPOSE
# -------
# Repair the area tables from Scripts 13 and 14 WITHOUT rerunning MaxEnt,
# future projections, GCM ensemble, 4-class classification, or transition maps.
#
# WHY THIS REPAIR IS REQUIRED
# ---------------------------
# The first area-calculation scripts assigned category labels to raster classes
# (for example, "Stable suitable"). terra::expanse(..., byValue=TRUE) then
# returned labels rather than numeric class codes. The old code attempted to
# convert those text labels with as.integer(), resulting in unmatched values
# and therefore zero areas in the CSV outputs.
#
# This repair script uses a robust geodesic calculation:
#   1. Compute the true geographic area of each grid cell once with
#      terra::cellSize(..., unit = "km", transform = TRUE).
#   2. Sum these cell areas by raster class using terra::zonal(..., "sum").
#   3. Match either numeric codes OR text labels to the known class legend.
#
# It overwrites ONLY area tables in:
#   - 13_MTSPS_4class_future_ensemble_area/03_area_tables/
#   - 14_binary_distribution_change_stability/03_area_tables/
#
# Existing classified maps and transition maps are NOT changed.
# =====================================================================

# ---------------------------------------------------------------------
# 0. PACKAGE
# ---------------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_if_missing("terra")

suppressPackageStartupMessages({
  library(terra)
})

# ---------------------------------------------------------------------
# 1. USER SETTINGS
# ---------------------------------------------------------------------

species_name <- "Camellia_japonica"
project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

current_model_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

classification_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "13_MTSPS_4class_future_ensemble_area"
)

transition_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "14_binary_distribution_change_stability"
)

repair_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "15_area_recalculation_repair"
)

# FALSE = reuse an existing geodesic cell-area raster when present.
# TRUE = regenerate only the repair support files.
clean_repair_support_files <- FALSE

# ---------------------------------------------------------------------
# 2. INPUT PATHS
# ---------------------------------------------------------------------

threshold_csv <- file.path(
  current_model_root,
  "04_CV_metrics_MTSPS_and_areas",
  "03_FINAL_single_MTSPS_threshold.csv"
)

classification_manifest_csv <- file.path(
  classification_root,
  "04_manifests",
  "classification_status_manifest.csv"
)

transition_manifest_csv <- file.path(
  transition_root,
  "04_manifests",
  "transition_status_manifest.csv"
)

# ---------------------------------------------------------------------
# 3. OVERWRITTEN OUTPUT TABLES
# ---------------------------------------------------------------------

classification_area_dir <- file.path(
  classification_root,
  "03_area_tables"
)

classification_long_csv <- file.path(
  classification_area_dir,
  "02_four_class_area_all_current_and_future_maps.csv"
)

classification_future_wide_csv <- file.path(
  classification_area_dir,
  "03_future_ensemble_four_class_area_wide.csv"
)

classification_current_csv <- file.path(
  classification_area_dir,
  "04_current_four_class_area.csv"
)

classification_suitable_summary_csv <- file.path(
  classification_area_dir,
  "05_total_suitable_area_summary.csv"
)

transition_area_dir <- file.path(
  transition_root,
  "03_area_tables"
)

transition_long_csv <- file.path(
  transition_area_dir,
  "02_transition_area_long.csv"
)

transition_wide_csv <- file.path(
  transition_area_dir,
  "03_transition_area_wide.csv"
)

distribution_change_summary_csv <- file.path(
  transition_area_dir,
  "04_distribution_change_summary.csv"
)

balance_check_csv <- file.path(
  transition_area_dir,
  "05_area_balance_checks.csv"
)

# Repair provenance outputs.
repair_support_dir <- file.path(
  repair_root,
  "01_supporting_files"
)

repair_manifest_dir <- file.path(
  repair_root,
  "02_manifests"
)

cell_area_tif <- file.path(
  repair_support_dir,
  "global_geodesic_cell_area_km2.tif"
)

four_class_repair_manifest_csv <- file.path(
  repair_manifest_dir,
  "four_class_area_repair_manifest.csv"
)

transition_repair_manifest_csv <- file.path(
  repair_manifest_dir,
  "transition_area_repair_manifest.csv"
)

repair_summary_csv <- file.path(
  repair_manifest_dir,
  "area_repair_run_summary.csv"
)

# ---------------------------------------------------------------------
# 4. HELPERS
# ---------------------------------------------------------------------

safe_remove <- function(path) {
  if (file.exists(path) || dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
}

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

find_column <- function(data, candidates, label) {
  names_original <- names(data)
  names_normalized <- normalise_name(names_original)
  candidates_normalized <- normalise_name(candidates)

  direct_hits <- match(
    candidates_normalized,
    names_normalized,
    nomatch = 0
  )

  direct_hits <- direct_hits[direct_hits > 0]

  if (length(direct_hits) > 0) {
    return(names_original[direct_hits[1]])
  }

  for (candidate in candidates_normalized) {
    partial_hits <- grep(candidate, names_normalized, fixed = TRUE)

    if (length(partial_hits) > 0) {
      return(names_original[partial_hits[1]])
    }
  }

  stop(
    "\nCould not identify required column in ", label, ".\n",
    "Available columns: ", paste(names_original, collapse = ", "),
    call. = FALSE
  )
}

read_mtsps_threshold <- function(path) {
  threshold_data <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  threshold_col <- find_column(
    threshold_data,
    c("MTSPS_threshold", "threshold"),
    "MTSPS threshold CSV"
  )

  values <- suppressWarnings(as.numeric(threshold_data[[threshold_col]]))
  values <- values[is.finite(values)]

  if (length(values) != 1) {
    stop(
      "\nExpected exactly one finite MTSPS threshold in:\n",
      path,
      call. = FALSE
    )
  }

  values[1]
}

assert_same_grid <- function(reference, candidate, label) {
  is_same <- terra::compareGeom(
    reference,
    candidate,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )

  if (!is_same) {
    stop(
      "\nRaster geometry mismatch: ", label,
      "\nAll classified and transition rasters must be exactly aligned with ",
      "the global template before geodesic area calculation.",
      call. = FALSE
    )
  }
}

find_zonal_columns <- function(zonal_output, label) {
  if (ncol(zonal_output) < 2) {
    stop(
      "\nUnexpected zonal output with fewer than two columns for ",
      label,
      call. = FALSE
    )
  }

  # terra::zonal returns zone values in the first column and the requested
  # statistic (sum) in the final column.
  list(
    zone = names(zonal_output)[1],
    sum = names(zonal_output)[ncol(zonal_output)]
  )
}

convert_zone_values_to_codes <- function(
  zones,
  lookup,
  code_column,
  label_columns,
  label
) {
  zone_text <- trimws(as.character(zones))

  # First handle zones returned directly as numeric class codes.
  numeric_codes <- suppressWarnings(as.integer(zone_text))

  numeric_valid <- !is.na(numeric_codes) &
    numeric_codes %in% lookup[[code_column]]

  if (all(numeric_valid)) {
    return(numeric_codes)
  }

  # Then handle categorical rasters, for which terra may return the active
  # category label rather than the stored integer code.
  for (label_column in label_columns) {
    if (!label_column %in% names(lookup)) {
      next
    }

    matched <- match(zone_text, as.character(lookup[[label_column]]))

    if (all(!is.na(matched))) {
      return(as.integer(lookup[[code_column]][matched]))
    }
  }

  stop(
    "\nCould not match zonal category values to known codes for ",
    label, ".\nReturned zones: ", paste(unique(zone_text), collapse = ", "),
    call. = FALSE
  )
}

calculate_geodesic_area_by_class <- function(
  class_raster,
  cell_area_raster,
  lookup,
  code_column,
  label_columns,
  label
) {
  assert_same_grid(cell_area_raster, class_raster, label)

  # Sum the precomputed geodesic area raster independently for every class.
  zonal_output <- terra::zonal(
    cell_area_raster,
    class_raster,
    fun = "sum",
    na.rm = TRUE
  )

  columns <- find_zonal_columns(zonal_output, label)

  class_codes <- convert_zone_values_to_codes(
    zones = zonal_output[[columns$zone]],
    lookup = lookup,
    code_column = code_column,
    label_columns = label_columns,
    label = label
  )

  observed <- data.frame(
    class_code = class_codes,
    area_km2 = as.numeric(zonal_output[[columns$sum]]),
    stringsAsFactors = FALSE
  )

  observed <- stats::aggregate(
    area_km2 ~ class_code,
    data = observed,
    FUN = sum
  )

  output <- merge(
    lookup,
    observed,
    by.x = code_column,
    by.y = "class_code",
    all.x = TRUE,
    sort = TRUE
  )

  output$area_km2[is.na(output$area_km2)] <- 0
  output$area_10_4_km2 <- output$area_km2 / 10000

  output
}

write_manifest <- function(records, path) {
  if (length(records) > 0) {
    write.csv(do.call(rbind, records), path, row.names = FALSE)
  }
}

# ---------------------------------------------------------------------
# 5. LOOKUP TABLES
# ---------------------------------------------------------------------

four_class_lookup <- data.frame(
  class_code = 1:4,
  class_name = c(
    "Unsuitable",
    "Low suitability",
    "Moderate suitability",
    "High suitability"
  ),
  class_definition = c(
    "Suitability < MTSPS",
    "MTSPS <= Suitability < 0.4",
    "0.4 <= Suitability < 0.6",
    "Suitability >= 0.6"
  ),
  stringsAsFactors = FALSE
)

transition_lookup <- data.frame(
  transition_code = 1:4,
  transition_name = c(
    "Stable unsuitable",
    "Expansion (gain)",
    "Contraction (loss)",
    "Stable suitable"
  ),
  transition_definition = c(
    "Current unsuitable -> Future unsuitable",
    "Current unsuitable -> Future suitable",
    "Current suitable -> Future unsuitable",
    "Current suitable -> Future suitable"
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------
# 6. INPUT VALIDATION
# ---------------------------------------------------------------------

required_paths <- c(
  MTSPS_threshold = threshold_csv,
  classification_manifest = classification_manifest_csv,
  transition_manifest = transition_manifest_csv
)

missing_paths <- required_paths[!file.exists(required_paths)]

if (length(missing_paths) > 0) {
  stop(
    "\nMissing required input(s):\n",
    paste(names(missing_paths), ": ", missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

mtsps_threshold <- read_mtsps_threshold(threshold_csv)

classification_manifest <- read.csv(
  classification_manifest_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

transition_manifest <- read.csv(
  transition_manifest_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Script 13 manifest columns.
class_map_id_col <- find_column(
  classification_manifest,
  c("map_id", "mapid"),
  "Script 13 classification manifest"
)

class_map_type_col <- find_column(
  classification_manifest,
  c("map_type", "maptype"),
  "Script 13 classification manifest"
)

class_scenario_col <- find_column(
  classification_manifest,
  c("scenario"),
  "Script 13 classification manifest"
)

class_ssp_col <- find_column(
  classification_manifest,
  c("ssp"),
  "Script 13 classification manifest"
)

class_period_start_col <- find_column(
  classification_manifest,
  c("period_start", "periodstart"),
  "Script 13 classification manifest"
)

class_period_end_col <- find_column(
  classification_manifest,
  c("period_end", "periodend"),
  "Script 13 classification manifest"
)

class_status_col <- find_column(
  classification_manifest,
  c("status"),
  "Script 13 classification manifest"
)

class_output_col <- find_column(
  classification_manifest,
  c("output_4class_tif", "output4classtif"),
  "Script 13 classification manifest"
)

# Script 14 manifest columns.
trans_map_id_col <- find_column(
  transition_manifest,
  c("map_id", "mapid"),
  "Script 14 transition manifest"
)

trans_ssp_col <- find_column(
  transition_manifest,
  c("ssp"),
  "Script 14 transition manifest"
)

trans_period_start_col <- find_column(
  transition_manifest,
  c("period_start", "periodstart"),
  "Script 14 transition manifest"
)

trans_period_end_col <- find_column(
  transition_manifest,
  c("period_end", "periodend"),
  "Script 14 transition manifest"
)

trans_status_col <- find_column(
  transition_manifest,
  c("status"),
  "Script 14 transition manifest"
)

trans_output_col <- find_column(
  transition_manifest,
  c("transition_tif", "transitiontif"),
  "Script 14 transition manifest"
)

completed_statuses <- c("completed", "skipped_existing_output")

class_rows <- classification_manifest[
  classification_manifest[[class_status_col]] %in% completed_statuses &
    !is.na(classification_manifest[[class_output_col]]) &
    nzchar(classification_manifest[[class_output_col]]),
  ,
  drop = FALSE
]

if (nrow(class_rows) != 18L) {
  stop(
    "\nExpected 18 completed Script 13 classified maps (2 current + 16 future), ",
    "but found ", nrow(class_rows), ".",
    call. = FALSE
  )
}

if (any(!file.exists(class_rows[[class_output_col]]))) {
  print(
    class_rows[
      !file.exists(class_rows[[class_output_col]]),
      c(class_map_id_col, class_output_col),
      drop = FALSE
    ]
  )

  stop(
    "\nOne or more Script 13 classified maps are missing.",
    call. = FALSE
  )
}

trans_rows <- transition_manifest[
  transition_manifest[[trans_status_col]] %in% completed_statuses &
    !is.na(transition_manifest[[trans_output_col]]) &
    nzchar(transition_manifest[[trans_output_col]]),
  ,
  drop = FALSE
]

if (nrow(trans_rows) != 16L) {
  stop(
    "\nExpected 16 completed Script 14 transition maps, but found ",
    nrow(trans_rows), ".",
    call. = FALSE
  )
}

if (any(!file.exists(trans_rows[[trans_output_col]]))) {
  print(
    trans_rows[
      !file.exists(trans_rows[[trans_output_col]]),
      c(trans_map_id_col, trans_output_col),
      drop = FALSE
    ]
  )

  stop(
    "\nOne or more Script 14 transition maps are missing.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 7. PREPARE REPAIR SUPPORT FILES
# ---------------------------------------------------------------------

if (clean_repair_support_files) {
  safe_remove(repair_root)
}

for (folder in c(
  classification_area_dir,
  transition_area_dir,
  repair_support_dir,
  repair_manifest_dir
)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

# Use the final full-data current map as the common global grid template.
template_index <- which(
  class_rows[[class_map_id_col]] == "Current_final_full_data"
)

if (length(template_index) != 1) {
  stop(
    "\nCould not identify exactly one Current_final_full_data classified map ",
    "in Script 13 manifest.",
    call. = FALSE
  )
}

template_raster <- terra::rast(
  class_rows[[class_output_col]][template_index]
)

if (
  file.exists(cell_area_tif) &&
  !clean_repair_support_files
) {
  cell_area <- terra::rast(cell_area_tif)
  assert_same_grid(template_raster, cell_area, "existing geodesic cell-area raster")
} else {
  cat("\nCreating one global geodesic cell-area raster in km²...\n")

  # cellSize computes latitude-corrected area for every grid cell.
  # mask = FALSE is deliberate: zonal() later masks via each class raster.
  cell_area <- terra::cellSize(
    template_raster,
    unit = "km",
    transform = TRUE,
    mask = FALSE
  )

  names(cell_area) <- "geodesic_cell_area_km2"

  terra::writeRaster(
    cell_area,
    cell_area_tif,
    overwrite = TRUE,
    wopt = list(gdal = "COMPRESS=LZW")
  )
}

cell_area_total_km2 <- as.numeric(
  terra::global(cell_area, "sum", na.rm = TRUE)[1, 1]
)

cat("\n======================================================================\n")
cat("AREA-TABLE REPAIR — Camellia japonica\n")
cat("======================================================================\n")
cat(
  "MTSPS threshold retained: ",
  format(mtsps_threshold, scientific = TRUE, digits = 8),
  "\n",
  sep = ""
)
cat(
  "Common global grid cell-area total (km²): ",
  format(cell_area_total_km2, big.mark = ",", scientific = FALSE),
  "\n",
  sep = ""
)
cat("Classified maps to recalculate: ", nrow(class_rows), "\n", sep = "")
cat("Transition maps to recalculate: ", nrow(trans_rows), "\n", sep = "")

# ---------------------------------------------------------------------
# 8. RECALCULATE SCRIPT 13 FOUR-CLASS AREAS
# ---------------------------------------------------------------------

class_rows <- class_rows[
  order(
    class_rows[[class_map_type_col]],
    class_rows[[class_ssp_col]],
    suppressWarnings(as.integer(class_rows[[class_period_start_col]]))
  ),
  ,
  drop = FALSE
]

four_class_area_records <- list()
four_class_manifest_records <- list()

for (index in seq_len(nrow(class_rows))) {
  row <- class_rows[index, , drop = FALSE]

  map_id <- as.character(row[[class_map_id_col]])
  map_type <- as.character(row[[class_map_type_col]])
  scenario <- as.character(row[[class_scenario_col]])
  ssp <- as.character(row[[class_ssp_col]])
  period_start <- suppressWarnings(as.integer(row[[class_period_start_col]]))
  period_end <- suppressWarnings(as.integer(row[[class_period_end_col]]))
  map_file <- as.character(row[[class_output_col]])

  cat("\n[Four-class ", index, "/", nrow(class_rows), "] ", map_id, "\n", sep = "")

  map_start <- Sys.time()

  result <- tryCatch(
    {
      map_raster <- terra::rast(map_file)
      assert_same_grid(template_raster, map_raster, map_id)

      area_data <- calculate_geodesic_area_by_class(
        class_raster = map_raster,
        cell_area_raster = cell_area,
        lookup = four_class_lookup,
        code_column = "class_code",
        label_columns = c("class_name"),
        label = map_id
      )

      area_data$map_id <- map_id
      area_data$map_type <- map_type
      area_data$source_model <- if (
        identical(map_id, "Current_CV100_average")
      ) {
        "Mean of 100 cross-validation replicates"
      } else if (identical(map_id, "Current_final_full_data")) {
        "Final full-data MaxEnt model"
      } else {
        "Mean of ACCESS-CM2, BCC-CSM2-MR, and CMCC-ESM2"
      }
      area_data$scenario <- scenario
      area_data$ssp <- ssp
      area_data$period_start <- period_start
      area_data$period_end <- period_end
      area_data$MTSPS_threshold <- mtsps_threshold

      area_data <- area_data[
        c(
          "map_id",
          "map_type",
          "source_model",
          "scenario",
          "ssp",
          "period_start",
          "period_end",
          "class_code",
          "class_name",
          "class_definition",
          "MTSPS_threshold",
          "area_km2",
          "area_10_4_km2"
        )
      ]

      total_valid_area <- sum(area_data$area_km2)

      area_data$percent_of_valid_environment <- if (total_valid_area > 0) {
        100 * area_data$area_km2 / total_valid_area
      } else {
        NA_real_
      }

      if (total_valid_area <= 0) {
        stop(
          "Calculated four-class total area is zero for ",
          map_id,
          call. = FALSE
        )
      }

      list(
        status = "completed",
        area_data = area_data,
        total_area = total_valid_area,
        message = ""
      )
    },
    error = function(error) {
      list(
        status = "failed",
        area_data = NULL,
        total_area = NA_real_,
        message = conditionMessage(error)
      )
    }
  )

  elapsed_minutes <- as.numeric(
    difftime(Sys.time(), map_start, units = "mins")
  )

  four_class_manifest_records[[length(four_class_manifest_records) + 1]] <- data.frame(
    map_id = map_id,
    map_type = map_type,
    input_classified_tif = normalizePath(map_file, winslash = "/", mustWork = TRUE),
    status = result$status,
    total_valid_environment_km2 = result$total_area,
    elapsed_minutes = elapsed_minutes,
    message = result$message,
    stringsAsFactors = FALSE
  )

  if (identical(result$status, "failed")) {
    write_manifest(
      four_class_manifest_records,
      four_class_repair_manifest_csv
    )

    stop(
      "\nFour-class area repair failed for ", map_id, ".\n",
      "Read the repair manifest:\n",
      four_class_repair_manifest_csv,
      call. = FALSE
    )
  }

  four_class_area_records[[length(four_class_area_records) + 1]] <- result$area_data

  cat(
    "  Valid domain area (km²): ",
    format(round(result$total_area, 2), big.mark = ","),
    "\n",
    sep = ""
  )

  write_manifest(
    four_class_manifest_records,
    four_class_repair_manifest_csv
  )
}

four_class_long <- do.call(rbind, four_class_area_records)

four_class_long <- four_class_long[
  order(
    four_class_long$map_type,
    four_class_long$ssp,
    four_class_long$period_start,
    four_class_long$class_code
  ),
  ,
  drop = FALSE
]

write.csv(
  four_class_long,
  classification_long_csv,
  row.names = FALSE
)

future_four_class_long <- four_class_long[
  four_class_long$map_type == "Future three-GCM ensemble mean",
  ,
  drop = FALSE
]

current_four_class_long <- four_class_long[
  grepl("^Current", four_class_long$map_type),
  ,
  drop = FALSE
]

write.csv(
  current_four_class_long,
  classification_current_csv,
  row.names = FALSE
)

future_keys <- unique(
  future_four_class_long[
    c(
      "map_id",
      "scenario",
      "ssp",
      "period_start",
      "period_end",
      "MTSPS_threshold"
    )
  ]
)

future_keys <- future_keys[
  order(future_keys$ssp, future_keys$period_start),
  ,
  drop = FALSE
]

future_wide_rows <- lapply(
  seq_len(nrow(future_keys)),
  function(index) {
    key <- future_keys[index, , drop = FALSE]

    subset <- future_four_class_long[
      future_four_class_long$map_id == key$map_id,
      ,
      drop = FALSE
    ]

    get_area <- function(class_name) {
      value <- subset$area_km2[subset$class_name == class_name]
      if (length(value) == 0) 0 else value[1]
    }

    unsuitable <- get_area("Unsuitable")
    low <- get_area("Low suitability")
    moderate <- get_area("Moderate suitability")
    high <- get_area("High suitability")
    total_valid <- unsuitable + low + moderate + high
    total_suitable <- low + moderate + high

    data.frame(
      scenario = key$scenario,
      ssp = key$ssp,
      period_start = key$period_start,
      period_end = key$period_end,
      MTSPS_threshold = key$MTSPS_threshold,
      unsuitable_km2 = unsuitable,
      low_suitability_km2 = low,
      moderate_suitability_km2 = moderate,
      high_suitability_km2 = high,
      total_suitable_km2 = total_suitable,
      total_suitable_10_4_km2 = total_suitable / 10000,
      total_valid_environment_km2 = total_valid,
      percent_suitable = if (total_valid > 0) {
        100 * total_suitable / total_valid
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
)

future_four_class_wide <- do.call(rbind, future_wide_rows)

write.csv(
  future_four_class_wide,
  classification_future_wide_csv,
  row.names = FALSE
)

suitable_summary <- do.call(
  rbind,
  lapply(
    split(four_class_long, four_class_long$map_id),
    function(table) {
      total_valid <- sum(table$area_km2)
      total_suitable <- sum(
        table$area_km2[table$class_code %in% c(2, 3, 4)]
      )

      data.frame(
        map_id = unique(table$map_id),
        map_type = unique(table$map_type),
        scenario = unique(table$scenario),
        ssp = unique(table$ssp),
        period_start = unique(table$period_start),
        period_end = unique(table$period_end),
        MTSPS_threshold = unique(table$MTSPS_threshold),
        total_valid_environment_km2 = total_valid,
        total_suitable_km2 = total_suitable,
        total_suitable_10_4_km2 = total_suitable / 10000,
        percent_suitable = if (total_valid > 0) {
          100 * total_suitable / total_valid
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  suitable_summary,
  classification_suitable_summary_csv,
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 9. RECALCULATE SCRIPT 14 TRANSITION AREAS
# ---------------------------------------------------------------------

trans_rows <- trans_rows[
  order(
    trans_rows[[trans_ssp_col]],
    suppressWarnings(as.integer(trans_rows[[trans_period_start_col]]))
  ),
  ,
  drop = FALSE
]

transition_area_records <- list()
transition_manifest_records <- list()

for (index in seq_len(nrow(trans_rows))) {
  row <- trans_rows[index, , drop = FALSE]

  map_id <- as.character(row[[trans_map_id_col]])
  ssp <- as.character(row[[trans_ssp_col]])
  period_start <- suppressWarnings(
    as.integer(row[[trans_period_start_col]])
  )
  period_end <- suppressWarnings(
    as.integer(row[[trans_period_end_col]])
  )
  transition_file <- as.character(row[[trans_output_col]])

  cat("\n[Transition ", index, "/", nrow(trans_rows), "] ", map_id, "\n", sep = "")

  map_start <- Sys.time()

  result <- tryCatch(
    {
      transition_raster <- terra::rast(transition_file)

      assert_same_grid(
        template_raster,
        transition_raster,
        map_id
      )

      area_data <- calculate_geodesic_area_by_class(
        class_raster = transition_raster,
        cell_area_raster = cell_area,
        lookup = transition_lookup,
        code_column = "transition_code",
        label_columns = c("transition_name"),
        label = map_id
      )

      area_data$map_id <- map_id
      area_data$ssp <- ssp
      area_data$period_start <- period_start
      area_data$period_end <- period_end
      area_data$MTSPS_threshold <- mtsps_threshold

      area_data <- area_data[
        c(
          "map_id",
          "ssp",
          "period_start",
          "period_end",
          "transition_code",
          "transition_name",
          "transition_definition",
          "MTSPS_threshold",
          "area_km2",
          "area_10_4_km2"
        )
      ]

      shared_valid_domain <- sum(area_data$area_km2)

      area_data$percent_of_shared_valid_domain <- if (
        shared_valid_domain > 0
      ) {
        100 * area_data$area_km2 / shared_valid_domain
      } else {
        NA_real_
      }

      if (shared_valid_domain <= 0) {
        stop(
          "Calculated transition total area is zero for ",
          map_id,
          call. = FALSE
        )
      }

      list(
        status = "completed",
        area_data = area_data,
        total_area = shared_valid_domain,
        message = ""
      )
    },
    error = function(error) {
      list(
        status = "failed",
        area_data = NULL,
        total_area = NA_real_,
        message = conditionMessage(error)
      )
    }
  )

  elapsed_minutes <- as.numeric(
    difftime(Sys.time(), map_start, units = "mins")
  )

  transition_manifest_records[[length(transition_manifest_records) + 1]] <- data.frame(
    map_id = map_id,
    input_transition_tif = normalizePath(
      transition_file,
      winslash = "/",
      mustWork = TRUE
    ),
    status = result$status,
    shared_valid_domain_km2 = result$total_area,
    elapsed_minutes = elapsed_minutes,
    message = result$message,
    stringsAsFactors = FALSE
  )

  if (identical(result$status, "failed")) {
    write_manifest(
      transition_manifest_records,
      transition_repair_manifest_csv
    )

    stop(
      "\nTransition area repair failed for ", map_id, ".\n",
      "Read the repair manifest:\n",
      transition_repair_manifest_csv,
      call. = FALSE
    )
  }

  transition_area_records[[length(transition_area_records) + 1]] <- result$area_data

  cat(
    "  Shared valid domain area (km²): ",
    format(round(result$total_area, 2), big.mark = ","),
    "\n",
    sep = ""
  )

  write_manifest(
    transition_manifest_records,
    transition_repair_manifest_csv
  )
}

transition_long <- do.call(rbind, transition_area_records)

transition_long <- transition_long[
  order(
    transition_long$ssp,
    transition_long$period_start,
    transition_long$transition_code
  ),
  ,
  drop = FALSE
]

write.csv(
  transition_long,
  transition_long_csv,
  row.names = FALSE
)

transition_keys <- unique(
  transition_long[
    c(
      "map_id",
      "ssp",
      "period_start",
      "period_end",
      "MTSPS_threshold"
    )
  ]
)

transition_keys <- transition_keys[
  order(transition_keys$ssp, transition_keys$period_start),
  ,
  drop = FALSE
]

transition_wide_rows <- lapply(
  seq_len(nrow(transition_keys)),
  function(index) {
    key <- transition_keys[index, , drop = FALSE]

    subset <- transition_long[
      transition_long$map_id == key$map_id,
      ,
      drop = FALSE
    ]

    get_area <- function(transition_name) {
      value <- subset$area_km2[
        subset$transition_name == transition_name
      ]

      if (length(value) == 0) 0 else value[1]
    }

    stable_unsuitable <- get_area("Stable unsuitable")
    expansion <- get_area("Expansion (gain)")
    contraction <- get_area("Contraction (loss)")
    stable_suitable <- get_area("Stable suitable")

    shared_valid_domain <- sum(subset$area_km2)
    current_suitable <- contraction + stable_suitable
    future_suitable <- expansion + stable_suitable
    net_change <- future_suitable - current_suitable

    data.frame(
      scenario = paste0(
        key$ssp,
        " ",
        key$period_start,
        "-",
        key$period_end
      ),
      ssp = key$ssp,
      period_start = key$period_start,
      period_end = key$period_end,
      MTSPS_threshold = key$MTSPS_threshold,
      stable_unsuitable_km2 = stable_unsuitable,
      expansion_gain_km2 = expansion,
      contraction_loss_km2 = contraction,
      stable_suitable_km2 = stable_suitable,
      current_suitable_km2 = current_suitable,
      future_suitable_km2 = future_suitable,
      net_suitable_change_km2 = net_change,
      net_suitable_change_10_4_km2 = net_change / 10000,
      percent_change_relative_to_current_suitable = if (
        current_suitable > 0
      ) {
        100 * net_change / current_suitable
      } else {
        NA_real_
      },
      shared_valid_domain_km2 = shared_valid_domain,
      expansion_percent_of_current_unsuitable = if (
        (stable_unsuitable + expansion) > 0
      ) {
        100 * expansion / (stable_unsuitable + expansion)
      } else {
        NA_real_
      },
      contraction_percent_of_current_suitable = if (
        current_suitable > 0
      ) {
        100 * contraction / current_suitable
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
)

transition_wide <- do.call(rbind, transition_wide_rows)

write.csv(
  transition_wide,
  transition_wide_csv,
  row.names = FALSE
)

distribution_change_summary <- transition_wide[
  c(
    "scenario",
    "ssp",
    "period_start",
    "period_end",
    "MTSPS_threshold",
    "current_suitable_km2",
    "future_suitable_km2",
    "net_suitable_change_km2",
    "net_suitable_change_10_4_km2",
    "percent_change_relative_to_current_suitable",
    "expansion_gain_km2",
    "contraction_loss_km2",
    "stable_suitable_km2",
    "stable_unsuitable_km2",
    "expansion_percent_of_current_unsuitable",
    "contraction_percent_of_current_suitable",
    "shared_valid_domain_km2"
  )
]

write.csv(
  distribution_change_summary,
  distribution_change_summary_csv,
  row.names = FALSE
)

balance_checks <- data.frame(
  scenario = transition_wide$scenario,
  ssp = transition_wide$ssp,
  period_start = transition_wide$period_start,
  period_end = transition_wide$period_end,
  current_suitable_balance_error_km2 =
    transition_wide$current_suitable_km2 -
    (
      transition_wide$contraction_loss_km2 +
        transition_wide$stable_suitable_km2
    ),
  future_suitable_balance_error_km2 =
    transition_wide$future_suitable_km2 -
    (
      transition_wide$expansion_gain_km2 +
        transition_wide$stable_suitable_km2
    ),
  net_change_balance_error_km2 =
    transition_wide$net_suitable_change_km2 -
    (
      transition_wide$expansion_gain_km2 -
        transition_wide$contraction_loss_km2
    ),
  stringsAsFactors = FALSE
)

write.csv(
  balance_checks,
  balance_check_csv,
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 10. FINAL QUALITY CHECKS AND REPORT
# ---------------------------------------------------------------------

four_class_totals <- aggregate(
  area_km2 ~ map_id,
  data = four_class_long,
  FUN = sum
)

transition_totals <- aggregate(
  area_km2 ~ map_id,
  data = transition_long,
  FUN = sum
)

if (any(four_class_totals$area_km2 <= 0)) {
  stop(
    "\nAt least one corrected four-class map still has zero total area.",
    call. = FALSE
  )
}

if (any(transition_totals$area_km2 <= 0)) {
  stop(
    "\nAt least one corrected transition map still has zero total area.",
    call. = FALSE
  )
}

maximum_balance_error <- max(
  abs(
    unlist(
      balance_checks[
        c(
          "current_suitable_balance_error_km2",
          "future_suitable_balance_error_km2",
          "net_change_balance_error_km2"
        )
      ]
    )
  ),
  na.rm = TRUE
)

run_summary <- data.frame(
  MTSPS_threshold = mtsps_threshold,
  corrected_four_class_maps = nrow(four_class_totals),
  corrected_transition_maps = nrow(transition_totals),
  global_cell_area_km2 = cell_area_total_km2,
  min_four_class_valid_domain_km2 = min(four_class_totals$area_km2),
  max_four_class_valid_domain_km2 = max(four_class_totals$area_km2),
  min_transition_shared_domain_km2 = min(transition_totals$area_km2),
  max_transition_shared_domain_km2 = max(transition_totals$area_km2),
  maximum_area_balance_error_km2 = maximum_balance_error,
  stringsAsFactors = FALSE
)

write.csv(
  run_summary,
  repair_summary_csv,
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("AREA-TABLE REPAIR COMPLETED — Camellia japonica\n")
cat("======================================================================\n")
cat("Corrected four-class maps: ", nrow(four_class_totals), "\n", sep = "")
cat("Corrected transition maps: ", nrow(transition_totals), "\n", sep = "")
cat(
  "Maximum area-balance error (km²): ",
  format(maximum_balance_error, scientific = TRUE, digits = 8),
  "\n",
  sep = ""
)
cat("\nCorrected Script 13 tables:\n")
cat(normalizePath(classification_long_csv, winslash = "/"), "\n")
cat(normalizePath(classification_future_wide_csv, winslash = "/"), "\n")
cat(normalizePath(classification_suitable_summary_csv, winslash = "/"), "\n")
cat("\nCorrected Script 14 tables:\n")
cat(normalizePath(transition_long_csv, winslash = "/"), "\n")
cat(normalizePath(transition_wide_csv, winslash = "/"), "\n")
cat(normalizePath(distribution_change_summary_csv, winslash = "/"), "\n")
cat("\nRepair summary:\n")
cat(normalizePath(repair_summary_csv, winslash = "/"), "\n")
