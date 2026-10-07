# =============================================================================
# 17A_REPAIR_AREA_TABLES_COMMON_VALID_DOMAIN_CAMELLIA_JAPONICA.R
# =============================================================================
# PURPOSE
# -------
# Repair Script 17 area tables using ONE common valid environmental domain:
# cells with finite continuous suitability values in the current full-data model
# AND in all 16 future three-GCM ensemble-mean rasters.
#
# WHY THIS REPAIR IS REQUIRED
# ---------------------------
# Script 17 correctly created all binary and transition maps, but its initial
# area-balance check compared:
#   - current/future binary areas calculated on each raster's own valid domain
# with
#   - transition areas calculated only where BOTH current and future rasters
#     are valid.
#
# When a small number of cells are valid in one map but NA in another, the
# domains differ, and balance errors can appear (e.g., ~3,334 km²) even though
# the binary and transition map logic is correct.
#
# THIS SCRIPT
# -----------
# 1. Does NOT rerun MaxEnt, future projections, ensemble averaging, binary
#    maps, or transition maps.
# 2. Builds a strict common valid domain across current + all 16 future
#    continuous ensemble-mean suitability rasters.
# 3. Recalculates every Script 17 area table within that common domain.
# 4. Overwrites only Script 17 area-table CSVs with corrected values.
#
# SCIENTIFIC IMPLICATION
# ----------------------
# All current-vs-future comparisons become spatially comparable because every
# scenario uses exactly the same valid grid cells. The MTSPS threshold remains
# the primary threshold; MTP and T10 remain supplementary sensitivity checks.
# =============================================================================

options(warn = 1)

cat("\n")
cat("======================================================================\n")
cat("STARTING SCRIPT 17A — COMMON-DOMAIN AREA-TABLE REPAIR\n")
cat("Species: Camellia japonica\n")
cat("======================================================================\n")

script_started_at <- Sys.time()

# =============================================================================
# A. USER SETTINGS
# =============================================================================

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

threshold_order <- c("MTSPS", "MTP", "T10")

# Continuous future suitability path column in the Script 12 ensemble manifest.
future_ensemble_mean_raster_column <- "ensemble_mean_tif"

# =============================================================================
# B. FIXED PATHS
# =============================================================================

final_model_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

future_ensemble_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "12_future_ensemble_3GCM"
)

script17_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "17_threshold_sensitivity_MTSPS_MTP_T10"
)

current_continuous_candidates <- c(
  file.path(
    final_model_root,
    "05_current_suitability_maps",
    "Camellia_japonica_current_final_full_data_cloglog.tif"
  ),
  file.path(
    final_model_root,
    "05_current_suitability_maps",
    "Camellia_japonica_current_final_full_data_cloglog.tiff"
  )
)

ensemble_manifest_candidates <- c(
  file.path(
    future_ensemble_root,
    "04_manifests",
    "ensemble_status_manifest.csv"
  ),
  file.path(
    future_ensemble_root,
    "04_manifests",
    "future_ensemble_status_manifest.csv"
  )
)

threshold_table_csv <- file.path(
  script17_root,
  "02_thresholds_and_current_scores",
  "01_threshold_values.csv"
)

binary_manifest_csv <- file.path(
  script17_root,
  "07_supporting_files",
  "binary_map_status_manifest.csv"
)

transition_manifest_csv <- file.path(
  script17_root,
  "07_supporting_files",
  "transition_map_status_manifest.csv"
)

support_dir <- file.path(script17_root, "07_supporting_files")
area_dir <- file.path(script17_root, "05_area_tables")
status_dir <- file.path(script17_root, "00_run_status")
log_dir <- file.path(script17_root, "08_logs")

common_domain_tif <- file.path(
  support_dir,
  "common_valid_domain_current_plus_all_16_future_ensembles.tif"
)

geodesic_area_tif <- file.path(
  support_dir,
  "global_geodesic_cell_area_km2.tif"
)

repair_status_csv <- file.path(
  support_dir,
  "17A_area_repair_status.csv"
)

common_domain_diagnostics_csv <- file.path(
  support_dir,
  "08_common_valid_domain_diagnostics.csv"
)

threshold_utility_csv <- file.path(
  support_dir,
  "09_threshold_utility_diagnostics.csv"
)

repair_summary_csv <- file.path(
  support_dir,
  "10_area_table_repair_summary.csv"
)

repair_started_marker <- file.path(
  status_dir,
  "17A_AREA_REPAIR_STARTED.txt"
)

repair_completed_marker <- file.path(
  status_dir,
  "17A_AREA_REPAIR_COMPLETED.txt"
)

repair_failed_marker <- file.path(
  status_dir,
  "17A_AREA_REPAIR_FAILED.txt"
)

repair_log <- file.path(
  log_dir,
  "17A_common_domain_area_repair_log.txt"
)

# Script 17 tables overwritten by this repair:
current_area_csv <- file.path(
  area_dir,
  "01_current_suitable_area_by_threshold.csv"
)

future_area_long_csv <- file.path(
  area_dir,
  "02_future_suitable_area_by_threshold_long.csv"
)

future_area_wide_csv <- file.path(
  area_dir,
  "03_future_suitable_area_by_threshold_wide.csv"
)

change_area_long_csv <- file.path(
  area_dir,
  "04_distribution_change_by_threshold_long.csv"
)

change_area_wide_csv <- file.path(
  area_dir,
  "05_distribution_change_by_threshold_wide.csv"
)

area_balance_checks_csv <- file.path(
  area_dir,
  "06_area_balance_checks.csv"
)

threshold_comparison_summary_csv <- file.path(
  area_dir,
  "07_threshold_comparison_summary.csv"
)

# =============================================================================
# C. HELPERS
# =============================================================================

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_log <- function(...) {
  line <- paste0(
    "[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] ",
    paste0(..., collapse = "")
  )
  cat(line, "\n")
  cat(line, "\n", file = repair_log, append = TRUE)
}

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

safe_number <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

find_column <- function(data, candidates, object_label) {
  raw_names <- names(data)
  normalized_names <- normalise_name(raw_names)
  normalized_candidates <- normalise_name(candidates)

  direct <- match(normalized_candidates, normalized_names, nomatch = 0)
  direct <- direct[direct > 0]

  if (length(direct) > 0) return(raw_names[direct[1]])

  for (candidate in normalized_candidates) {
    hit <- grep(candidate, normalized_names, fixed = TRUE)
    if (length(hit) > 0) return(raw_names[hit[1]])
  }

  stop(
    "\nCannot identify a required column in ", object_label, ".\n",
    "Available columns: ", paste(raw_names, collapse = ", "),
    call. = FALSE
  )
}

find_first_existing_file <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop(
      "\nNo ", label, " found. Checked:\n",
      paste(candidates, collapse = "\n"),
      call. = FALSE
    )
  }

  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

find_existing_path_column <- function(data, candidate_columns, label) {
  for (candidate in candidate_columns) {
    if (candidate %in% names(data)) {
      values <- as.character(data[[candidate]])
      valid <- !is.na(values) & nzchar(values) & file.exists(values)

      if (sum(valid) > 0) return(candidate)
    }
  }

  # Fallback: search all columns for the most existing TIFF paths.
  scores <- vapply(names(data), function(column) {
    values <- as.character(data[[column]])
    sum(
      !is.na(values) &
        nzchar(values) &
        file.exists(values) &
        grepl("\\.(tif|tiff)$", values, ignore.case = TRUE)
    )
  }, numeric(1))

  if (max(scores) <= 0) {
    stop(
      "\nNo existing raster-path column found in ", label, ".",
      call. = FALSE
    )
  }

  best <- names(scores)[scores == max(scores)]

  if (length(best) != 1) {
    stop(
      "\nMore than one possible raster path column in ", label, ":\n",
      paste(best, collapse = ", "),
      call. = FALSE
    )
  }

  best
}

assert_same_grid <- function(reference, candidate, label) {
  same_grid <- terra::compareGeom(
    reference,
    candidate,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )

  if (!same_grid) {
    stop(
      "\nRaster geometry mismatch: ", label, ".\n",
      "All continuous, binary, transition, and area rasters must align.",
      call. = FALSE
    )
  }
}

calculate_area_by_code <- function(class_raster, cell_area_raster, domain_raster,
                                   codes, class_names) {
  assert_same_grid(cell_area_raster, class_raster, "class raster")
  assert_same_grid(domain_raster, class_raster, "common analysis domain")

  class_in_domain <- terra::mask(class_raster, domain_raster)

  zonal_output <- terra::zonal(
    cell_area_raster,
    class_in_domain,
    fun = "sum",
    na.rm = TRUE
  )

  output <- data.frame(
    code = codes,
    class_name = class_names,
    area_km2 = 0,
    stringsAsFactors = FALSE
  )

  if (!is.null(zonal_output) && nrow(zonal_output) > 0) {
    zone_codes <- suppressWarnings(
      as.integer(round(as.numeric(zonal_output[[1]])))
    )
    zone_area <- as.numeric(zonal_output[[ncol(zonal_output)]])

    for (i in seq_along(codes)) {
      matched <- zone_area[zone_codes == codes[i]]
      if (length(matched) > 0) {
        output$area_km2[i] <- sum(matched, na.rm = TRUE)
      }
    }
  }

  output$area_10_4_km2 <- output$area_km2 / 10000
  output
}

get_area <- function(area_table, code) {
  value <- area_table$area_km2[area_table$code == code]

  if (length(value) != 1) {
    stop("Could not retrieve one area for code ", code, call. = FALSE)
  }

  value
}

common_domain_area <- function(domain_raster, cell_area_raster) {
  area_table <- calculate_area_by_code(
    class_raster = domain_raster,
    cell_area_raster = cell_area_raster,
    domain_raster = domain_raster,
    codes = 1,
    class_names = "Common valid environmental domain"
  )

  area_table$area_km2[1]
}

# =============================================================================
# D. MAIN REPAIR
# =============================================================================

run_common_domain_area_repair <- function() {

  for (folder in c(status_dir, support_dir, area_dir, log_dir)) {
    ensure_dir(folder)
  }

  unlink(repair_completed_marker, force = TRUE)
  unlink(repair_failed_marker, force = TRUE)

  writeLines(
    c(
      "Script 17A area-table repair started",
      paste("Started:", format(script_started_at, "%Y-%m-%d %H:%M:%S")),
      "Purpose: calculate all areas in a strict common valid environmental domain"
    ),
    repair_started_marker
  )

  if (!requireNamespace("terra", quietly = TRUE)) {
    install.packages("terra", repos = "https://cloud.r-project.org")
  }

  suppressPackageStartupMessages(library(terra))

  write_log("Package terra loaded.")

  # -------------------------------------------------------------------
  # D1. Read input files
  # -------------------------------------------------------------------

  current_continuous_path <- find_first_existing_file(
    current_continuous_candidates,
    "current full-data continuous suitability raster"
  )

  ensemble_manifest_path <- find_first_existing_file(
    ensemble_manifest_candidates,
    "future ensemble manifest"
  )

  required_files <- c(
    threshold_table_csv,
    binary_manifest_csv,
    transition_manifest_csv
  )

  missing_required <- required_files[!file.exists(required_files)]

  if (length(missing_required) > 0) {
    stop(
      "\nRequired Script 17 output file(s) missing:\n",
      paste(missing_required, collapse = "\n"),
      call. = FALSE
    )
  }

  threshold_table <- read.csv(
    threshold_table_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  threshold_id_col <- find_column(
    threshold_table,
    c("threshold_id"),
    "Script 17 threshold table"
  )

  threshold_value_col <- find_column(
    threshold_table,
    c("threshold_cloglog"),
    "Script 17 threshold table"
  )

  threshold_table$threshold_id_analysis <- as.character(
    threshold_table[[threshold_id_col]]
  )
  threshold_table$threshold_value_analysis <- safe_number(
    threshold_table[[threshold_value_col]]
  )

  threshold_table <- threshold_table[
    match(threshold_order, threshold_table$threshold_id_analysis),
    ,
    drop = FALSE
  ]

  if (
    nrow(threshold_table) != length(threshold_order) ||
    any(is.na(threshold_table$threshold_id_analysis)) ||
    any(!is.finite(threshold_table$threshold_value_analysis))
  ) {
    stop(
      "\nThreshold table does not contain valid MTSPS, MTP, and T10 rows.",
      call. = FALSE
    )
  }

  binary_manifest <- read.csv(
    binary_manifest_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  transition_manifest <- read.csv(
    transition_manifest_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  ensemble_manifest <- read.csv(
    ensemble_manifest_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # -------------------------------------------------------------------
  # D2. Identify manifest fields
  # -------------------------------------------------------------------

  ensemble_status_col <- find_column(
    ensemble_manifest,
    c("status"),
    "future ensemble manifest"
  )

  ensemble_ssp_col <- find_column(
    ensemble_manifest,
    c("ssp"),
    "future ensemble manifest"
  )

  ensemble_period_start_col <- find_column(
    ensemble_manifest,
    c("period_start", "periodstart"),
    "future ensemble manifest"
  )

  ensemble_period_end_col <- find_column(
    ensemble_manifest,
    c("period_end", "periodend"),
    "future ensemble manifest"
  )

  if (!future_ensemble_mean_raster_column %in% names(ensemble_manifest)) {
    stop(
      "\nRequired ensemble-mean raster column not found:\n",
      future_ensemble_mean_raster_column,
      "\nAvailable columns:\n",
      paste(names(ensemble_manifest), collapse = ", "),
      call. = FALSE
    )
  }

  ensemble_path_col <- future_ensemble_mean_raster_column

  binary_threshold_col <- find_column(
    binary_manifest,
    c("threshold_id"),
    "binary map manifest"
  )

  binary_map_type_col <- find_column(
    binary_manifest,
    c("map_type", "maptype"),
    "binary map manifest"
  )

  binary_ssp_col <- find_column(
    binary_manifest,
    c("ssp"),
    "binary map manifest"
  )

  binary_period_start_col <- find_column(
    binary_manifest,
    c("period_start", "periodstart"),
    "binary map manifest"
  )

  binary_period_end_col <- find_column(
    binary_manifest,
    c("period_end", "periodend"),
    "binary map manifest"
  )

  binary_path_col <- find_existing_path_column(
    binary_manifest,
    c("binary_tif"),
    "binary map manifest"
  )

  transition_threshold_col <- find_column(
    transition_manifest,
    c("threshold_id"),
    "transition map manifest"
  )

  transition_ssp_col <- find_column(
    transition_manifest,
    c("ssp"),
    "transition map manifest"
  )

  transition_period_start_col <- find_column(
    transition_manifest,
    c("period_start", "periodstart"),
    "transition map manifest"
  )

  transition_period_end_col <- find_column(
    transition_manifest,
    c("period_end", "periodend"),
    "transition map manifest"
  )

  transition_path_col <- find_existing_path_column(
    transition_manifest,
    c("transition_tif"),
    "transition map manifest"
  )

  # -------------------------------------------------------------------
  # D3. Strict common valid domain: current + all 16 future ensembles
  # -------------------------------------------------------------------

  current_continuous <- terra::rast(current_continuous_path)

  valid_ensemble_status <- c("completed", "skipped_existing_output")

  future_ensemble_rows <- ensemble_manifest[
    ensemble_manifest[[ensemble_status_col]] %in% valid_ensemble_status &
      !is.na(ensemble_manifest[[ensemble_path_col]]) &
      nzchar(ensemble_manifest[[ensemble_path_col]]) &
      file.exists(ensemble_manifest[[ensemble_path_col]]),
    ,
    drop = FALSE
  ]

  if (nrow(future_ensemble_rows) != 16L) {
    stop(
      "\nExpected exactly 16 valid future ensemble mean rasters; found ",
      nrow(future_ensemble_rows), ".",
      call. = FALSE
    )
  }

  future_ensemble_rows$ssp_analysis <- as.character(
    future_ensemble_rows[[ensemble_ssp_col]]
  )
  future_ensemble_rows$period_start_analysis <- safe_number(
    future_ensemble_rows[[ensemble_period_start_col]]
  )
  future_ensemble_rows$period_end_analysis <- safe_number(
    future_ensemble_rows[[ensemble_period_end_col]]
  )
  future_ensemble_rows$continuous_path_analysis <- normalizePath(
    as.character(future_ensemble_rows[[ensemble_path_col]]),
    winslash = "/",
    mustWork = TRUE
  )

  future_ensemble_rows <- future_ensemble_rows[
    order(
      future_ensemble_rows$ssp_analysis,
      future_ensemble_rows$period_start_analysis
    ),
    ,
    drop = FALSE
  ]

  # Start with cells valid in the current continuous map.
  common_domain <- terra::ifel(is.na(current_continuous), NA, 1)

  write_log(
    "Building strict common valid domain across current + ",
    nrow(future_ensemble_rows),
    " future ensemble rasters."
  )

  for (i in seq_len(nrow(future_ensemble_rows))) {
    future_continuous <- terra::rast(
      future_ensemble_rows$continuous_path_analysis[i]
    )

    assert_same_grid(
      current_continuous,
      future_continuous,
      paste0(
        future_ensemble_rows$ssp_analysis[i],
        " ",
        future_ensemble_rows$period_start_analysis[i],
        "-",
        future_ensemble_rows$period_end_analysis[i]
      )
    )

    # mask() retains cells only when the mask has non-NA values.
    common_domain <- terra::mask(common_domain, future_continuous)
  }

  names(common_domain) <- "common_valid_domain"

  terra::writeRaster(
    common_domain,
    common_domain_tif,
    overwrite = TRUE,
    datatype = "INT1U",
    wopt = list(gdal = "COMPRESS=LZW")
  )

  # Geodesic cell area.
  if (file.exists(geodesic_area_tif)) {
    cell_area <- terra::rast(geodesic_area_tif)
    assert_same_grid(current_continuous, cell_area, "existing geodesic area raster")
  } else {
    write_log("Creating geodesic cell-area raster (km²).")

    cell_area <- terra::cellSize(
      current_continuous,
      unit = "km",
      transform = TRUE,
      mask = FALSE
    )

    names(cell_area) <- "geodesic_cell_area_km2"

    terra::writeRaster(
      cell_area,
      geodesic_area_tif,
      overwrite = TRUE,
      wopt = list(gdal = "COMPRESS=LZW")
    )
  }

  common_valid_domain_km2 <- common_domain_area(
    common_domain,
    cell_area
  )

  # -------------------------------------------------------------------
  # D4. Verify binary and transition maps
  # -------------------------------------------------------------------

  binary_manifest$threshold_id_analysis <- as.character(
    binary_manifest[[binary_threshold_col]]
  )
  binary_manifest$map_type_analysis <- as.character(
    binary_manifest[[binary_map_type_col]]
  )
  binary_manifest$ssp_analysis <- as.character(
    binary_manifest[[binary_ssp_col]]
  )
  binary_manifest$period_start_analysis <- safe_number(
    binary_manifest[[binary_period_start_col]]
  )
  binary_manifest$period_end_analysis <- safe_number(
    binary_manifest[[binary_period_end_col]]
  )
  binary_manifest$binary_path_analysis <- as.character(
    binary_manifest[[binary_path_col]]
  )

  transition_manifest$threshold_id_analysis <- as.character(
    transition_manifest[[transition_threshold_col]]
  )
  transition_manifest$ssp_analysis <- as.character(
    transition_manifest[[transition_ssp_col]]
  )
  transition_manifest$period_start_analysis <- safe_number(
    transition_manifest[[transition_period_start_col]]
  )
  transition_manifest$period_end_analysis <- safe_number(
    transition_manifest[[transition_period_end_col]]
  )
  transition_manifest$transition_path_analysis <- as.character(
    transition_manifest[[transition_path_col]]
  )

  binary_expected <- length(threshold_order) * 17L
  transition_expected <- length(threshold_order) * 16L

  valid_binary_rows <- binary_manifest[
    binary_manifest$threshold_id_analysis %in% threshold_order &
      file.exists(binary_manifest$binary_path_analysis),
    ,
    drop = FALSE
  ]

  valid_transition_rows <- transition_manifest[
    transition_manifest$threshold_id_analysis %in% threshold_order &
      file.exists(transition_manifest$transition_path_analysis),
    ,
    drop = FALSE
  ]

  if (nrow(valid_binary_rows) != binary_expected) {
    stop(
      "\nExpected ", binary_expected,
      " Script 17 binary maps, but found ", nrow(valid_binary_rows), ".",
      call. = FALSE
    )
  }

  if (nrow(valid_transition_rows) != transition_expected) {
    stop(
      "\nExpected ", transition_expected,
      " Script 17 transition maps, but found ", nrow(valid_transition_rows), ".",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------
  # D5. Current and future binary suitable areas in common domain
  # -------------------------------------------------------------------

  current_area_records <- list()
  future_area_records <- list()
  repair_records <- list()

  for (threshold_id in threshold_order) {

    threshold_value <- threshold_table$threshold_value_analysis[
      threshold_table$threshold_id_analysis == threshold_id
    ]

    current_row <- valid_binary_rows[
      valid_binary_rows$threshold_id_analysis == threshold_id &
        grepl("^current", tolower(valid_binary_rows$map_type_analysis)),
      ,
      drop = FALSE
    ]

    if (nrow(current_row) != 1L) {
      stop(
        "\nExpected exactly one current binary map for threshold ",
        threshold_id, ".",
        call. = FALSE
      )
    }

    current_binary <- terra::rast(current_row$binary_path_analysis)
    assert_same_grid(current_continuous, current_binary, paste0(threshold_id, " current binary"))

    current_binary_area <- calculate_area_by_code(
      class_raster = current_binary,
      cell_area_raster = cell_area,
      domain_raster = common_domain,
      codes = c(0, 1),
      class_names = c("Unsuitable", "Suitable")
    )

    current_suitable <- get_area(current_binary_area, 1)
    current_unsuitable <- get_area(current_binary_area, 0)

    current_area_records[[length(current_area_records) + 1]] <- data.frame(
      threshold_id = threshold_id,
      threshold_cloglog = threshold_value,
      analysis_domain = "Common valid domain: current + all 16 future ensemble rasters",
      current_unsuitable_km2 = current_unsuitable,
      current_suitable_km2 = current_suitable,
      current_suitable_10_4_km2 = current_suitable / 10000,
      current_valid_environment_km2 = current_unsuitable + current_suitable,
      current_percent_suitable = 100 * current_suitable /
        (current_unsuitable + current_suitable),
      current_binary_tif = normalizePath(
        current_row$binary_path_analysis,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )

    future_binary_rows <- valid_binary_rows[
      valid_binary_rows$threshold_id_analysis == threshold_id &
        grepl("^future", tolower(valid_binary_rows$map_type_analysis)),
      ,
      drop = FALSE
    ]

    if (nrow(future_binary_rows) != 16L) {
      stop(
        "\nExpected 16 future binary maps for threshold ",
        threshold_id, ", but found ", nrow(future_binary_rows), ".",
        call. = FALSE
      )
    }

    for (i in seq_len(nrow(future_ensemble_rows))) {
      ssp <- future_ensemble_rows$ssp_analysis[i]
      period_start <- as.integer(future_ensemble_rows$period_start_analysis[i])
      period_end <- as.integer(future_ensemble_rows$period_end_analysis[i])

      future_row <- future_binary_rows[
        future_binary_rows$ssp_analysis == ssp &
          future_binary_rows$period_start_analysis == period_start &
          future_binary_rows$period_end_analysis == period_end,
        ,
        drop = FALSE
      ]

      if (nrow(future_row) != 1L) {
        stop(
          "\nCannot uniquely match a future binary map for ",
          threshold_id, " | ", ssp, " ", period_start, "-", period_end,
          call. = FALSE
        )
      }

      future_binary <- terra::rast(future_row$binary_path_analysis)
      assert_same_grid(
        current_continuous,
        future_binary,
        paste0(threshold_id, " ", ssp, " ", period_start, "-", period_end)
      )

      future_binary_area <- calculate_area_by_code(
        class_raster = future_binary,
        cell_area_raster = cell_area,
        domain_raster = common_domain,
        codes = c(0, 1),
        class_names = c("Unsuitable", "Suitable")
      )

      future_suitable <- get_area(future_binary_area, 1)
      future_unsuitable <- get_area(future_binary_area, 0)

      future_area_records[[length(future_area_records) + 1]] <- data.frame(
        threshold_id = threshold_id,
        threshold_cloglog = threshold_value,
        analysis_domain = "Common valid domain: current + all 16 future ensemble rasters",
        ssp = ssp,
        period_start = period_start,
        period_end = period_end,
        scenario = paste0(ssp, " ", period_start, "-", period_end),
        current_suitable_km2 = current_suitable,
        future_unsuitable_km2 = future_unsuitable,
        future_suitable_km2 = future_suitable,
        future_suitable_10_4_km2 = future_suitable / 10000,
        future_valid_environment_km2 = future_unsuitable + future_suitable,
        future_percent_suitable = 100 * future_suitable /
          (future_unsuitable + future_suitable),
        net_suitable_change_km2 = future_suitable - current_suitable,
        percent_change_relative_to_current = 100 *
          (future_suitable - current_suitable) / current_suitable,
        future_binary_tif = normalizePath(
          future_row$binary_path_analysis,
          winslash = "/",
          mustWork = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  current_area <- do.call(rbind, current_area_records)
  future_area_long <- do.call(rbind, future_area_records)

  current_area <- current_area[
    match(threshold_order, current_area$threshold_id),
    ,
    drop = FALSE
  ]

  future_area_long <- future_area_long[
    order(
      match(future_area_long$threshold_id, threshold_order),
      future_area_long$ssp,
      future_area_long$period_start
    ),
    ,
    drop = FALSE
  ]

  # -------------------------------------------------------------------
  # D6. Transition areas in the same common domain
  # -------------------------------------------------------------------

  change_area_records <- list()

  for (threshold_id in threshold_order) {

    threshold_value <- threshold_table$threshold_value_analysis[
      threshold_table$threshold_id_analysis == threshold_id
    ]

    current_suitable <- current_area$current_suitable_km2[
      current_area$threshold_id == threshold_id
    ]

    threshold_transition_rows <- valid_transition_rows[
      valid_transition_rows$threshold_id_analysis == threshold_id,
      ,
      drop = FALSE
    ]

    if (nrow(threshold_transition_rows) != 16L) {
      stop(
        "\nExpected 16 transition maps for threshold ",
        threshold_id, ", but found ", nrow(threshold_transition_rows), ".",
        call. = FALSE
      )
    }

    for (i in seq_len(nrow(future_ensemble_rows))) {
      ssp <- future_ensemble_rows$ssp_analysis[i]
      period_start <- as.integer(future_ensemble_rows$period_start_analysis[i])
      period_end <- as.integer(future_ensemble_rows$period_end_analysis[i])

      transition_row <- threshold_transition_rows[
        threshold_transition_rows$ssp_analysis == ssp &
          threshold_transition_rows$period_start_analysis == period_start &
          threshold_transition_rows$period_end_analysis == period_end,
        ,
        drop = FALSE
      ]

      if (nrow(transition_row) != 1L) {
        stop(
          "\nCannot uniquely match a transition map for ",
          threshold_id, " | ", ssp, " ", period_start, "-", period_end,
          call. = FALSE
        )
      }

      transition_raster <- terra::rast(
        transition_row$transition_path_analysis
      )

      assert_same_grid(
        current_continuous,
        transition_raster,
        paste0(
          "transition ", threshold_id,
          " ", ssp, " ", period_start, "-", period_end
        )
      )

      transition_area <- calculate_area_by_code(
        class_raster = transition_raster,
        cell_area_raster = cell_area,
        domain_raster = common_domain,
        codes = 1:4,
        class_names = c(
          "Stable unsuitable",
          "Expansion (gain)",
          "Contraction (loss)",
          "Stable suitable"
        )
      )

      stable_unsuitable <- get_area(transition_area, 1)
      expansion <- get_area(transition_area, 2)
      contraction <- get_area(transition_area, 3)
      stable_suitable <- get_area(transition_area, 4)

      future_suitable <- expansion + stable_suitable
      net_change <- future_suitable - current_suitable

      change_area_records[[length(change_area_records) + 1]] <- data.frame(
        threshold_id = threshold_id,
        threshold_cloglog = threshold_value,
        analysis_domain = "Common valid domain: current + all 16 future ensemble rasters",
        ssp = ssp,
        period_start = period_start,
        period_end = period_end,
        scenario = paste0(ssp, " ", period_start, "-", period_end),
        stable_unsuitable_km2 = stable_unsuitable,
        expansion_gain_km2 = expansion,
        contraction_loss_km2 = contraction,
        stable_suitable_km2 = stable_suitable,
        current_suitable_km2 = current_suitable,
        future_suitable_km2 = future_suitable,
        net_suitable_change_km2 = net_change,
        net_suitable_change_10_4_km2 = net_change / 10000,
        percent_change_relative_to_current_suitable = 100 *
          net_change / current_suitable,
        expansion_percent_of_current_unsuitable = 100 * expansion /
          (stable_unsuitable + expansion),
        contraction_percent_of_current_suitable = 100 * contraction /
          current_suitable,
        shared_valid_domain_km2 = stable_unsuitable + expansion +
          contraction + stable_suitable,
        transition_tif = normalizePath(
          transition_row$transition_path_analysis,
          winslash = "/",
          mustWork = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  change_area_long <- do.call(rbind, change_area_records)

  change_area_long <- change_area_long[
    order(
      match(change_area_long$threshold_id, threshold_order),
      change_area_long$ssp,
      change_area_long$period_start
    ),
    ,
    drop = FALSE
  ]

  # -------------------------------------------------------------------
  # D7. Wide tables
  # -------------------------------------------------------------------

  scenario_keys <- unique(
    future_area_long[
      c("ssp", "period_start", "period_end", "scenario"),
      drop = FALSE
    ]
  )

  scenario_keys <- scenario_keys[
    order(scenario_keys$ssp, scenario_keys$period_start),
    ,
    drop = FALSE
  ]

  future_area_wide <- do.call(
    rbind,
    lapply(seq_len(nrow(scenario_keys)), function(i) {
      key <- scenario_keys[i, , drop = FALSE]
      row_output <- key

      for (threshold_id in threshold_order) {
        subset <- future_area_long[
          future_area_long$threshold_id == threshold_id &
            future_area_long$ssp == key$ssp &
            future_area_long$period_start == key$period_start,
          ,
          drop = FALSE
        ]

        row_output[[paste0(threshold_id, "_future_suitable_km2")]] <-
          subset$future_suitable_km2
        row_output[[paste0(threshold_id, "_net_change_km2")]] <-
          subset$net_suitable_change_km2
        row_output[[paste0(threshold_id, "_percent_change")]] <-
          subset$percent_change_relative_to_current
      }

      row_output
    })
  )

  change_area_wide <- do.call(
    rbind,
    lapply(seq_len(nrow(scenario_keys)), function(i) {
      key <- scenario_keys[i, , drop = FALSE]
      row_output <- key

      for (threshold_id in threshold_order) {
        subset <- change_area_long[
          change_area_long$threshold_id == threshold_id &
            change_area_long$ssp == key$ssp &
            change_area_long$period_start == key$period_start,
          ,
          drop = FALSE
        ]

        metrics <- c(
          "expansion_gain_km2",
          "contraction_loss_km2",
          "stable_suitable_km2",
          "net_suitable_change_km2",
          "percent_change_relative_to_current_suitable"
        )

        for (metric in metrics) {
          row_output[[paste0(threshold_id, "_", metric)]] <- subset[[metric]]
        }
      }

      row_output
    })
  )

  # -------------------------------------------------------------------
  # D8. Area-balance checks and diagnostics
  # -------------------------------------------------------------------

  area_balance_checks <- do.call(
    rbind,
    lapply(seq_len(nrow(change_area_long)), function(i) {
      row <- change_area_long[i, , drop = FALSE]

      future_row <- future_area_long[
        future_area_long$threshold_id == row$threshold_id &
          future_area_long$ssp == row$ssp &
          future_area_long$period_start == row$period_start,
        ,
        drop = FALSE
      ]

      current_row <- current_area[
        current_area$threshold_id == row$threshold_id,
        ,
        drop = FALSE
      ]

      data.frame(
        threshold_id = row$threshold_id,
        ssp = row$ssp,
        period_start = row$period_start,
        period_end = row$period_end,
        current_suitable_balance_error_km2 =
          current_row$current_suitable_km2 -
          (row$contraction_loss_km2 + row$stable_suitable_km2),
        future_suitable_balance_error_km2 =
          future_row$future_suitable_km2 -
          (row$expansion_gain_km2 + row$stable_suitable_km2),
        net_change_balance_error_km2 =
          row$net_suitable_change_km2 -
          (row$expansion_gain_km2 - row$contraction_loss_km2),
        stringsAsFactors = FALSE
      )
    })
  )

  max_balance_error <- max(
    abs(
      unlist(
        area_balance_checks[
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

  threshold_utility <- do.call(
    rbind,
    lapply(threshold_order, function(threshold_id) {
      row <- current_area[
        current_area$threshold_id == threshold_id,
        ,
        drop = FALSE
      ]

      p <- row$current_percent_suitable

      data.frame(
        threshold_id = threshold_id,
        threshold_cloglog = row$threshold_cloglog,
        current_percent_suitable = p,
        is_degenerate_global_binary_threshold = p <= 0.001 || p >= 99.999,
        interpretation = if (p >= 99.999) {
          "All cells in the common valid domain are classified as suitable; report only as an uninformative lower-bound sensitivity threshold."
        } else if (p <= 0.001) {
          "Almost no cells in the common valid domain are classified as suitable; interpret cautiously as an extreme threshold."
        } else {
          "Non-degenerate binary threshold suitable for supplementary comparison."
        },
        stringsAsFactors = FALSE
      )
    })
  )

  common_domain_diagnostics <- data.frame(
    analysis_domain = "Cells finite in current full-data suitability and all 16 future ensemble-mean suitability rasters",
    current_continuous_tif = current_continuous_path,
    future_ensemble_count = nrow(future_ensemble_rows),
    common_valid_domain_km2 = common_valid_domain_km2,
    common_valid_domain_10_4_km2 = common_valid_domain_km2 / 10000,
    stringsAsFactors = FALSE
  )

  repair_summary <- data.frame(
    repair = "Script 17A common valid domain area-table recalculation",
    thresholds_recalculated = length(threshold_order),
    current_binary_maps_reused = length(threshold_order),
    future_binary_maps_reused = nrow(future_area_long),
    transition_maps_reused = nrow(change_area_long),
    common_valid_domain_km2 = common_valid_domain_km2,
    maximum_area_balance_error_km2 = max_balance_error,
    stringsAsFactors = FALSE
  )

  # -------------------------------------------------------------------
  # D9. Overwrite Script 17 area tables and write diagnostics
  # -------------------------------------------------------------------

  write.csv(current_area, current_area_csv, row.names = FALSE)
  write.csv(future_area_long, future_area_long_csv, row.names = FALSE)
  write.csv(future_area_wide, future_area_wide_csv, row.names = FALSE)
  write.csv(change_area_long, change_area_long_csv, row.names = FALSE)
  write.csv(change_area_wide, change_area_wide_csv, row.names = FALSE)
  write.csv(area_balance_checks, area_balance_checks_csv, row.names = FALSE)
  write.csv(
    threshold_utility,
    threshold_utility_csv,
    row.names = FALSE
  )

  threshold_comparison_summary <- do.call(
    rbind,
    lapply(threshold_order, function(threshold_id) {
      current_row <- current_area[
        current_area$threshold_id == threshold_id,
        ,
        drop = FALSE
      ]

      change_rows <- change_area_long[
        change_area_long$threshold_id == threshold_id,
        ,
        drop = FALSE
      ]

      utility_row <- threshold_utility[
        threshold_utility$threshold_id == threshold_id,
        ,
        drop = FALSE
      ]

      data.frame(
        threshold_id = threshold_id,
        threshold_cloglog = current_row$threshold_cloglog,
        analysis_domain = current_row$analysis_domain,
        current_suitable_km2 = current_row$current_suitable_km2,
        current_percent_suitable = current_row$current_percent_suitable,
        minimum_future_suitable_km2 = min(change_rows$future_suitable_km2),
        maximum_future_suitable_km2 = max(change_rows$future_suitable_km2),
        minimum_net_change_km2 = min(change_rows$net_suitable_change_km2),
        maximum_net_change_km2 = max(change_rows$net_suitable_change_km2),
        scenarios_with_net_gain = sum(change_rows$net_suitable_change_km2 > 0),
        scenarios_with_net_loss = sum(change_rows$net_suitable_change_km2 < 0),
        scenarios_with_expansion_larger_than_contraction = sum(
          change_rows$expansion_gain_km2 >
            change_rows$contraction_loss_km2
        ),
        scenarios_with_contraction_larger_than_expansion = sum(
          change_rows$contraction_loss_km2 >
            change_rows$expansion_gain_km2
        ),
        mean_stable_suitable_percent_of_current = mean(
          100 * change_rows$stable_suitable_km2 /
            current_row$current_suitable_km2
        ),
        is_degenerate_global_binary_threshold =
          utility_row$is_degenerate_global_binary_threshold,
        threshold_interpretation = utility_row$interpretation,
        stringsAsFactors = FALSE
      )
    })
  )

  write.csv(
    threshold_comparison_summary,
    threshold_comparison_summary_csv,
    row.names = FALSE
  )

  write.csv(
    common_domain_diagnostics,
    common_domain_diagnostics_csv,
    row.names = FALSE
  )

  write.csv(repair_summary, repair_summary_csv, row.names = FALSE)

  repair_status <- data.frame(
    status = "completed",
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    message = "All Script 17 area tables recalculated in strict common valid environmental domain.",
    stringsAsFactors = FALSE
  )

  write.csv(repair_status, repair_status_csv, row.names = FALSE)

  writeLines(
    c(
      "SCRIPT 17A AREA-TABLE REPAIR COMPLETED SUCCESSFULLY",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("Common valid domain (km2):", format(common_valid_domain_km2, scientific = FALSE)),
      paste("Maximum area-balance error (km2):", format(max_balance_error, scientific = TRUE, digits = 8)),
      paste("Corrected future-area table:", future_area_long_csv),
      paste("Corrected transition-area table:", change_area_long_csv),
      paste("Corrected area-balance table:", area_balance_checks_csv)
    ),
    repair_completed_marker
  )

  write_log("Area-table repair completed successfully.")

  cat("\n======================================================================\n")
  cat("SCRIPT 17A COMPLETED — COMMON-DOMAIN AREA-TABLE REPAIR\n")
  cat("======================================================================\n")
  cat(
    "Common valid domain (km²): ",
    format(common_valid_domain_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat(
    "Maximum area-balance error (km²): ",
    format(max_balance_error, scientific = TRUE, digits = 8),
    "\n",
    sep = ""
  )
  cat("\nCorrected tables:\n")
  cat(current_area_csv, "\n")
  cat(future_area_long_csv, "\n")
  cat(change_area_long_csv, "\n")
  cat(area_balance_checks_csv, "\n")
  cat("\nThreshold utility diagnostic:\n")
  cat(threshold_utility_csv, "\n")
  cat("\nDo not interpret the MTP maps as substantive if the diagnostic flags MTP as degenerate.\n")

  invisible(list(
    current_area = current_area,
    future_area = future_area_long,
    change_area = change_area_long,
    threshold_utility = threshold_utility,
    repair_summary = repair_summary
  ))
}

# =============================================================================
# E. EXECUTE WITH EXPLICIT FAILURE REPORTING
# =============================================================================

tryCatch(
  {
    run_common_domain_area_repair()
  },
  error = function(error) {

    dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    failure_message <- paste0(
      "SCRIPT 17A AREA-TABLE REPAIR FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(error), "\n\n",
      "Repair log: ", repair_log, "\n"
    )

    writeLines(failure_message, repair_failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 17A AREA-TABLE REPAIR FAILED\n")
    cat("======================================================================\n")
    cat(failure_message, "\n")

    stop(error)
  }
)
