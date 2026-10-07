# =====================================================================
# 14_BINARY_DISTRIBUTION_CHANGE_STABILITY_CAMELLIA_JAPONICA.R
# =====================================================================
# PURPOSE
# -------
# Compare the final full-data CURRENT MaxEnt distribution with each of the
# 16 three-GCM future ensemble distributions using the SAME fixed MTSPS
# threshold already applied in Script 13.
#
# BASELINE
# --------
# Current final full-data MaxEnt map only:
#   Camellia_japonica_current_final_full_data_MTSPS_4classes.tif
#
# BINARY SUITABILITY DEFINITION
# -----------------------------
# Unsuitable = class 1 = suitability < MTSPS
# Suitable   = classes 2, 3, or 4 = suitability >= MTSPS
#
# TRANSITION CODES
# ----------------
# 1 = Stable unsuitable: Current unsuitable -> Future unsuitable
# 2 = Expansion (gain):  Current unsuitable -> Future suitable
# 3 = Contraction (loss):Current suitable   -> Future unsuitable
# 4 = Stable suitable:   Current suitable   -> Future suitable
#
# OUTPUTS
# -------
# 16 transition GeoTIFFs, one per future SSP x period
# Global geodesic area tables in km² and 10^4 km²
# Consistency diagnostics:
#   current suitable = contraction + stable suitable
#   future suitable  = expansion + stable suitable
#   net change       = expansion - contraction
#
# IMPORTANT
# ---------
# This script does NOT use the CV100 average current map as the change baseline.
# The baseline is the final full-data MaxEnt current map, while the future map
# is the continuous three-GCM ensemble classified in Script 13.
# =====================================================================

# ---------------------------------------------------------------------
# 0. PACKAGE
# ---------------------------------------------------------------------

install_if_missing <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cloud.r-project.org")
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

# Completed Script 13 output.
classification_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "13_MTSPS_4class_future_ensemble_area"
)

# Final current model folder, used only to record and verify the fixed MTSPS.
current_model_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

# New output folder created by Script 14.
output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "14_binary_distribution_change_stability"
)

# Set TRUE only for a completely fresh rerun. This deletes output_root only.
clean_output_root <- FALSE

# Safe resume setting. Existing transition GeoTIFFs are retained and their
# area statistics are recalculated.
skip_existing_transition_maps <- TRUE

# Stop at the first failed scenario to avoid wasting time when a shared raster
# geometry/path problem exists.
stop_after_first_failure <- TRUE

# ---------------------------------------------------------------------
# 2. INPUT AND OUTPUT PATHS
# ---------------------------------------------------------------------

classification_manifest_csv <- file.path(
  classification_root,
  "04_manifests",
  "classification_status_manifest.csv"
)

threshold_csv <- file.path(
  current_model_root,
  "04_CV_metrics_MTSPS_and_areas",
  "03_FINAL_single_MTSPS_threshold.csv"
)

transition_map_root <- file.path(
  output_root,
  "01_binary_transition_maps"
)

binary_map_root <- file.path(
  output_root,
  "02_binary_suitability_maps"
)

area_table_root <- file.path(
  output_root,
  "03_area_tables"
)

manifest_root <- file.path(
  output_root,
  "04_manifests"
)

transition_manifest_csv <- file.path(
  manifest_root,
  "transition_status_manifest.csv"
)

transition_legend_csv <- file.path(
  area_table_root,
  "01_binary_transition_legend.csv"
)

transition_area_long_csv <- file.path(
  area_table_root,
  "02_transition_area_long.csv"
)

transition_area_wide_csv <- file.path(
  area_table_root,
  "03_transition_area_wide.csv"
)

distribution_change_summary_csv <- file.path(
  area_table_root,
  "04_distribution_change_summary.csv"
)

balance_check_csv <- file.path(
  area_table_root,
  "05_area_balance_checks.csv"
)

run_summary_csv <- file.path(
  manifest_root,
  "transition_run_summary.csv"
)

# ---------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ---------------------------------------------------------------------

safe_remove <- function(path) {
  if (file.exists(path) || dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
}

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

find_column <- function(data, candidates, object_label) {
  original_names <- names(data)
  normalized_names <- normalise_name(original_names)
  normalized_candidates <- normalise_name(candidates)

  direct_matches <- match(
    normalized_candidates,
    normalized_names,
    nomatch = 0
  )

  direct_matches <- direct_matches[direct_matches > 0]

  if (length(direct_matches) > 0) {
    return(original_names[direct_matches[1]])
  }

  for (candidate in normalized_candidates) {
    partial_matches <- grep(
      candidate,
      normalized_names,
      fixed = TRUE
    )

    if (length(partial_matches) > 0) {
      return(original_names[partial_matches[1]])
    }
  }

  stop(
    "\nCould not identify required column in ", object_label, ".\n",
    "Available columns: ", paste(original_names, collapse = ", "),
    call. = FALSE
  )
}

read_single_MTSPS_threshold <- function(path) {
  threshold_table <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  threshold_column <- find_column(
    threshold_table,
    c("MTSPS_threshold", "threshold"),
    "MTSPS threshold CSV"
  )

  threshold_values <- suppressWarnings(
    as.numeric(threshold_table[[threshold_column]])
  )

  threshold_values <- threshold_values[
    is.finite(threshold_values)
  ]

  if (length(threshold_values) != 1) {
    stop(
      "\nThe threshold CSV must contain exactly one finite MTSPS value.",
      call. = FALSE
    )
  }

  threshold_values[1]
}

assert_same_grid <- function(reference_raster, candidate_raster, label) {
  same_grid <- terra::compareGeom(
    reference_raster,
    candidate_raster,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )

  if (!same_grid) {
    stop(
      "\nCurrent and future classified maps do not have identical geometry for ",
      label, ".\n",
      "No automatic resampling is performed because binary distribution-change ",
      "analysis requires exact cell-to-cell correspondence.",
      call. = FALSE
    )
  }
}

validate_four_class_raster <- function(class_raster, label) {
  range_values <- terra::global(
    class_raster,
    fun = c("min", "max"),
    na.rm = TRUE
  )

  minimum <- as.numeric(range_values[1, "min"])
  maximum <- as.numeric(range_values[1, "max"])

  if (!is.finite(minimum) || !is.finite(maximum)) {
    stop(
      "\nNo finite class values are present in ",
      label,
      call. = FALSE
    )
  }

  if (minimum < 1 || maximum > 4) {
    stop(
      "\nExpected four-class integer values 1-4 in ",
      label,
      ". Observed range: ", minimum, " to ", maximum,
      call. = FALSE
    )
  }
}

build_binary_suitability <- function(four_class_raster) {
  # Class 1 is unsuitable. Classes 2, 3, and 4 are suitable.
  binary <- terra::ifel(
    is.na(four_class_raster),
    NA,
    terra::ifel(four_class_raster >= 2, 1, 0)
  )

  names(binary) <- "suitable_MTSPS_binary"

  binary
}

build_transition_raster <- function(
  current_binary,
  future_binary
) {
  # Use NA wherever either current or future prediction is NA.
  transition <- terra::ifel(
    is.na(current_binary) | is.na(future_binary),
    NA,
    terra::ifel(
      current_binary == 0 & future_binary == 0,
      1,
      terra::ifel(
        current_binary == 0 & future_binary == 1,
        2,
        terra::ifel(
          current_binary == 1 & future_binary == 0,
          3,
          4
        )
      )
    )
  )

  names(transition) <- "binary_distribution_transition"

  transition
}

set_transition_categories <- function(transition_raster) {
  levels(transition_raster) <- data.frame(
    ID = 1:4,
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

  transition_raster
}

calculate_transition_areas <- function(
  transition_raster,
  map_id,
  ssp,
  period_start,
  period_end,
  mtsps_threshold
) {
  observed <- terra::expanse(
    transition_raster,
    unit = "km",
    transform = TRUE,
    byValue = TRUE
  )

  value_column <- find_column(
    observed,
    c("value"),
    paste0(map_id, " area output")
  )

  area_column <- find_column(
    observed,
    c("area"),
    paste0(map_id, " area output")
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

  observed_table <- data.frame(
    transition_code = as.integer(observed[[value_column]]),
    area_km2 = as.numeric(observed[[area_column]]),
    stringsAsFactors = FALSE
  )

  output <- merge(
    transition_lookup,
    observed_table,
    by = "transition_code",
    all.x = TRUE,
    sort = TRUE
  )

  output$area_km2[is.na(output$area_km2)] <- 0
  output$area_10_4_km2 <- output$area_km2 / 10000

  total_valid_area <- sum(output$area_km2)

  output$percent_of_shared_valid_domain <- if (total_valid_area > 0) {
    100 * output$area_km2 / total_valid_area
  } else {
    NA_real_
  }

  output$map_id <- map_id
  output$ssp <- ssp
  output$period_start <- period_start
  output$period_end <- period_end
  output$MTSPS_threshold <- mtsps_threshold

  output <- output[
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
      "area_10_4_km2",
      "percent_of_shared_valid_domain"
    )
  ]

  output
}

build_wide_transition_table <- function(long_table) {
  keys <- unique(
    long_table[
      c("map_id", "ssp", "period_start", "period_end", "MTSPS_threshold")
    ]
  )

  keys <- keys[
    order(keys$ssp, keys$period_start),
    ,
    drop = FALSE
  ]

  output_rows <- lapply(
    seq_len(nrow(keys)),
    function(index) {
      key <- keys[index, , drop = FALSE]

      subset <- long_table[
        long_table$map_id == key$map_id,
        ,
        drop = FALSE
      ]

      get_area <- function(transition_name) {
        values <- subset$area_km2[
          subset$transition_name == transition_name
        ]

        if (length(values) == 0) 0 else values[1]
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
        percent_change_relative_to_current_suitable = if (current_suitable > 0) {
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

  do.call(rbind, output_rows)
}

build_balance_checks <- function(wide_table) {
  # These identities must equal zero, apart from negligible numeric rounding:
  # Current suitable = contraction + stable suitable
  # Future suitable = expansion + stable suitable
  # Net change = expansion - contraction
  data.frame(
    scenario = wide_table$scenario,
    ssp = wide_table$ssp,
    period_start = wide_table$period_start,
    period_end = wide_table$period_end,
    current_suitable_balance_error_km2 =
      wide_table$current_suitable_km2 -
      (
        wide_table$contraction_loss_km2 +
          wide_table$stable_suitable_km2
      ),
    future_suitable_balance_error_km2 =
      wide_table$future_suitable_km2 -
      (
        wide_table$expansion_gain_km2 +
          wide_table$stable_suitable_km2
      ),
    net_change_balance_error_km2 =
      wide_table$net_suitable_change_km2 -
      (
        wide_table$expansion_gain_km2 -
          wide_table$contraction_loss_km2
      ),
    stringsAsFactors = FALSE
  )
}

write_progress_manifest <- function(records) {
  if (length(records) > 0) {
    write.csv(
      do.call(rbind, records),
      transition_manifest_csv,
      row.names = FALSE
    )
  }
}

# ---------------------------------------------------------------------
# 4. VERIFY INPUTS
# ---------------------------------------------------------------------

required_paths <- c(
  classification_manifest = classification_manifest_csv,
  MTSPS_threshold = threshold_csv
)

missing_paths <- required_paths[
  !file.exists(required_paths)
]

if (length(missing_paths) > 0) {
  stop(
    "\nRequired input file(s) are missing:\n",
    paste(names(missing_paths), ": ", missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

mtsps_threshold <- read_single_MTSPS_threshold(threshold_csv)

classification_manifest <- read.csv(
  classification_manifest_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

map_id_column <- find_column(
  classification_manifest,
  c("map_id", "mapid"),
  "Script 13 classification manifest"
)

status_column <- find_column(
  classification_manifest,
  c("status"),
  "Script 13 classification manifest"
)

output_tif_column <- find_column(
  classification_manifest,
  c("output_4class_tif", "output4classtif"),
  "Script 13 classification manifest"
)

ssp_column <- find_column(
  classification_manifest,
  c("ssp"),
  "Script 13 classification manifest"
)

period_start_column <- find_column(
  classification_manifest,
  c("period_start", "periodstart"),
  "Script 13 classification manifest"
)

period_end_column <- find_column(
  classification_manifest,
  c("period_end", "periodend"),
  "Script 13 classification manifest"
)

usable_statuses <- c("completed", "skipped_existing_output")

current_baseline_rows <- classification_manifest[
  classification_manifest[[map_id_column]] == "Current_final_full_data" &
    classification_manifest[[status_column]] %in% usable_statuses,
  ,
  drop = FALSE
]

if (nrow(current_baseline_rows) != 1) {
  stop(
    "\nThe Script 13 manifest must contain exactly one completed current ",
    "final full-data baseline row with map_id = Current_final_full_data.",
    call. = FALSE
  )
}

current_baseline_tif <- as.character(
  current_baseline_rows[[output_tif_column]][1]
)

if (!file.exists(current_baseline_tif)) {
  stop(
    "\nCurrent final full-data classified baseline map is missing:\n",
    current_baseline_tif,
    call. = FALSE
  )
}

future_rows <- classification_manifest[
  grepl(
    "^Future_ensemble_",
    classification_manifest[[map_id_column]]
  ) &
    classification_manifest[[status_column]] %in% usable_statuses &
    !is.na(classification_manifest[[output_tif_column]]) &
    nzchar(classification_manifest[[output_tif_column]]),
  ,
  drop = FALSE
]

if (nrow(future_rows) != 16L) {
  stop(
    "\nExpected 16 completed classified future ensemble maps from Script 13, ",
    "but found ", nrow(future_rows), ".",
    call. = FALSE
  )
}

if (any(!file.exists(future_rows[[output_tif_column]]))) {
  missing_future <- future_rows[
    !file.exists(future_rows[[output_tif_column]]),
    c(map_id_column, output_tif_column),
    drop = FALSE
  ]

  print(missing_future)

  stop(
    "\nAt least one classified future ensemble map is missing.",
    call. = FALSE
  )
}

if (anyDuplicated(future_rows[[map_id_column]])) {
  stop(
    "\nDuplicate future map IDs were found in the Script 13 manifest.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 5. PREPARE OUTPUT FOLDERS
# ---------------------------------------------------------------------

if (clean_output_root) {
  safe_remove(output_root)
}

for (folder in c(
  transition_map_root,
  binary_map_root,
  area_table_root,
  manifest_root
)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

transition_legend <- data.frame(
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
  binary_suitability_definition = c(
    "Unsuitable = suitability < MTSPS",
    "Suitable = suitability >= MTSPS",
    "Suitable = suitability >= MTSPS",
    "Suitable = suitability >= MTSPS"
  ),
  MTSPS_threshold = rep(mtsps_threshold, 4),
  stringsAsFactors = FALSE
)

write.csv(
  transition_legend,
  transition_legend_csv,
  row.names = FALSE
)

cat("\n======================================================================\n")
cat("BINARY DISTRIBUTION CHANGE + STABILITY — Camellia japonica\n")
cat("======================================================================\n")
cat(
  "Fixed current-model MTSPS threshold: ",
  format(mtsps_threshold, scientific = TRUE, digits = 8),
  "\n",
  sep = ""
)
cat("Current baseline: final full-data MaxEnt map\n")
cat("Future maps: 16 three-GCM ensemble maps\n")
cat("Output root:\n", normalizePath(output_root, winslash = "/"), "\n", sep = "")

# ---------------------------------------------------------------------
# 6. PREPARE CURRENT BINARY SUITABILITY BASELINE
# ---------------------------------------------------------------------

current_four_class <- terra::rast(current_baseline_tif)

validate_four_class_raster(
  current_four_class,
  "current final full-data classified baseline"
)

current_binary_tif <- file.path(
  binary_map_root,
  "Camellia_japonica_current_final_full_data_MTSPS_binary_suitable.tif"
)

if (
  !skip_existing_transition_maps ||
  !file.exists(current_binary_tif)
) {
  current_binary <- build_binary_suitability(current_four_class)

  terra::writeRaster(
    current_binary,
    current_binary_tif,
    overwrite = TRUE,
    wopt = list(
      gdal = "COMPRESS=LZW",
      datatype = "INT1U"
    )
  )
} else {
  current_binary <- terra::rast(current_binary_tif)
}

# ---------------------------------------------------------------------
# 7. CREATE THE 16 FUTURE TRANSITION MAPS
# ---------------------------------------------------------------------

future_rows <- future_rows[
  order(
    future_rows[[ssp_column]],
    as.integer(future_rows[[period_start_column]])
  ),
  ,
  drop = FALSE
]

status_records <- list()
area_records <- list()
run_start_time <- Sys.time()

for (index in seq_len(nrow(future_rows))) {
  row <- future_rows[index, , drop = FALSE]

  future_map_id <- as.character(row[[map_id_column]])
  future_four_class_tif <- as.character(row[[output_tif_column]])
  ssp <- as.character(row[[ssp_column]])
  period_start <- as.integer(row[[period_start_column]])
  period_end <- as.integer(row[[period_end_column]])

  ensemble_id <- sub(
    "^Future_ensemble_",
    "",
    future_map_id
  )

  scenario_label <- paste0(
    ssp,
    " ",
    period_start,
    "-",
    period_end
  )

  scenario_transition_dir <- file.path(
    transition_map_root,
    ensemble_id
  )

  scenario_binary_dir <- file.path(
    binary_map_root,
    ensemble_id
  )

  dir.create(
    scenario_transition_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  dir.create(
    scenario_binary_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  future_binary_tif <- file.path(
    scenario_binary_dir,
    paste0(
      species_name,
      "_",
      ensemble_id,
      "_MTSPS_binary_suitable.tif"
    )
  )

  transition_tif <- file.path(
    scenario_transition_dir,
    paste0(
      species_name,
      "_",
      ensemble_id,
      "_binary_distribution_transition.tif"
    )
  )

  map_start_time <- Sys.time()

  cat(
    "\n[", index, "/16] ", scenario_label, "\n",
    sep = ""
  )

  result <- tryCatch(
    {
      future_four_class <- terra::rast(future_four_class_tif)

      validate_four_class_raster(
        future_four_class,
        scenario_label
      )

      assert_same_grid(
        current_four_class,
        future_four_class,
        scenario_label
      )

      if (
        !skip_existing_transition_maps ||
        !file.exists(future_binary_tif)
      ) {
        future_binary <- build_binary_suitability(future_four_class)

        terra::writeRaster(
          future_binary,
          future_binary_tif,
          overwrite = TRUE,
          wopt = list(
            gdal = "COMPRESS=LZW",
            datatype = "INT1U"
          )
        )
      } else {
        future_binary <- terra::rast(future_binary_tif)
      }

      if (
        !skip_existing_transition_maps ||
        !file.exists(transition_tif)
      ) {
        transition <- build_transition_raster(
          current_binary,
          future_binary
        )

        transition <- set_transition_categories(transition)

        terra::writeRaster(
          transition,
          transition_tif,
          overwrite = TRUE,
          wopt = list(
            gdal = "COMPRESS=LZW",
            datatype = "INT1U"
          )
        )

        output_status <- "completed"
      } else {
        transition <- terra::rast(transition_tif)
        output_status <- "skipped_existing_output"
      }

      transition_area_table <- calculate_transition_areas(
        transition_raster = transition,
        map_id = future_map_id,
        ssp = ssp,
        period_start = period_start,
        period_end = period_end,
        mtsps_threshold = mtsps_threshold
      )

      list(
        status = output_status,
        area_table = transition_area_table,
        message = ""
      )
    },
    error = function(error) {
      list(
        status = "failed",
        area_table = NULL,
        message = conditionMessage(error)
      )
    }
  )

  elapsed_minutes <- as.numeric(
    difftime(Sys.time(), map_start_time, units = "mins")
  )

  status_records[[length(status_records) + 1]] <- data.frame(
    map_id = future_map_id,
    scenario = scenario_label,
    ssp = ssp,
    period_start = period_start,
    period_end = period_end,
    status = result$status,
    current_baseline_4class_tif = normalizePath(
      current_baseline_tif,
      winslash = "/",
      mustWork = TRUE
    ),
    future_4class_tif = normalizePath(
      future_four_class_tif,
      winslash = "/",
      mustWork = TRUE
    ),
    current_binary_tif = normalizePath(
      current_binary_tif,
      winslash = "/",
      mustWork = TRUE
    ),
    future_binary_tif = if (file.exists(future_binary_tif)) {
      normalizePath(future_binary_tif, winslash = "/", mustWork = TRUE)
    } else {
      NA_character_
    },
    transition_tif = if (file.exists(transition_tif)) {
      normalizePath(transition_tif, winslash = "/", mustWork = TRUE)
    } else {
      NA_character_
    },
    MTSPS_threshold = mtsps_threshold,
    elapsed_minutes = elapsed_minutes,
    message = result$message,
    stringsAsFactors = FALSE
  )

  if (!is.null(result$area_table)) {
    area_records[[length(area_records) + 1]] <- result$area_table
  }

  write_progress_manifest(status_records)

  if (identical(result$status, "failed")) {
    cat("  Status: FAILED.\n")

    if (stop_after_first_failure) {
      stop(
        "\nTransition analysis stopped after the first failure.\n",
        "Read the status manifest:\n",
        transition_manifest_csv,
        call. = FALSE
      )
    }
  } else {
    cat(
      "  Status: ",
      result$status,
      " in ",
      format(round(elapsed_minutes, 2), nsmall = 2),
      " minutes.\n",
      sep = ""
    )
  }
}

# ---------------------------------------------------------------------
# 8. WRITE TABLES AND QUALITY-CONTROL SUMMARIES
# ---------------------------------------------------------------------

transition_long_table <- do.call(rbind, area_records)

transition_long_table <- transition_long_table[
  order(
    transition_long_table$ssp,
    transition_long_table$period_start,
    transition_long_table$transition_code
  ),
  ,
  drop = FALSE
]

write.csv(
  transition_long_table,
  transition_area_long_csv,
  row.names = FALSE
)

transition_wide_table <- build_wide_transition_table(
  transition_long_table
)

write.csv(
  transition_wide_table,
  transition_area_wide_csv,
  row.names = FALSE
)

distribution_change_summary <- transition_wide_table[
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

balance_check_table <- build_balance_checks(
  transition_wide_table
)

write.csv(
  balance_check_table,
  balance_check_csv,
  row.names = FALSE
)

status_table <- do.call(rbind, status_records)

write_progress_manifest(status_records)

completed_n <- sum(status_table$status == "completed")
skipped_n <- sum(status_table$status == "skipped_existing_output")
failed_n <- sum(status_table$status == "failed")

maximum_balance_error <- max(
  abs(
    unlist(
      balance_check_table[
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

total_runtime_hours <- as.numeric(
  difftime(Sys.time(), run_start_time, units = "hours")
)

run_summary <- data.frame(
  MTSPS_threshold = mtsps_threshold,
  expected_future_transition_maps = 16,
  completed_this_run = completed_n,
  skipped_existing_outputs = skipped_n,
  failed = failed_n,
  maximum_area_balance_error_km2 = maximum_balance_error,
  total_runtime_hours = total_runtime_hours,
  stringsAsFactors = FALSE
)

write.csv(
  run_summary,
  run_summary_csv,
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 9. FINAL REPORT
# ---------------------------------------------------------------------

cat("\n======================================================================\n")
cat("BINARY DISTRIBUTION CHANGE + STABILITY ANALYSIS COMPLETED\n")
cat("======================================================================\n")
cat(
  "MTSPS threshold used unchanged: ",
  format(mtsps_threshold, scientific = TRUE, digits = 8),
  "\n",
  sep = ""
)
cat("Baseline: Current final full-data MaxEnt model\n")
cat("Completed maps this run: ", completed_n, "\n", sep = "")
cat("Skipped existing transition maps: ", skipped_n, "\n", sep = "")
cat("Failed maps: ", failed_n, "\n", sep = "")
cat(
  "Maximum area-balance error (km²): ",
  format(maximum_balance_error, scientific = TRUE, digits = 8),
  "\n",
  sep = ""
)
cat("\nTransition area table (long):\n")
cat(normalizePath(transition_area_long_csv, winslash = "/"), "\n")
cat("\nTransition area table (wide):\n")
cat(normalizePath(transition_area_wide_csv, winslash = "/"), "\n")
cat("\nDistribution-change summary:\n")
cat(normalizePath(distribution_change_summary_csv, winslash = "/"), "\n")
cat("\nArea balance checks:\n")
cat(normalizePath(balance_check_csv, winslash = "/"), "\n")
cat("\nTransition status manifest:\n")
cat(normalizePath(transition_manifest_csv, winslash = "/"), "\n")

if (failed_n == 0 && (completed_n + skipped_n) == 16L) {
  cat(
    "\nAll 16 future ensemble scenarios have stable unsuitable, expansion, ",
    "contraction, and stable suitable maps and areas.\n"
  )
} else {
  cat(
    "\nResolve all failed transition rows before interpreting distribution change.\n"
  )
}
