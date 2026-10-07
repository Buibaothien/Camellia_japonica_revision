# =============================================================================
# 18_INTER_GCM_UNCERTAINTY_AGREEMENT_EXTRAPOLATION_CAMELLIA_JAPONICA.R
# =============================================================================
# PURPOSE
# -------
# Quantify future-projection uncertainty for Camellia japonica across the three
# CMIP6 GCMs used in the final workflow:
#
#   ACCESS-CM2, BCC-CSM2-MR, CMCC-ESM2
#
# The primary suitability threshold is the previously fixed MTSPS threshold.
#
# CORE OUTPUTS FOR EACH OF 16 SSP × PERIOD SCENARIOS
# --------------------------------------------------
# 1. GCM suitable-agreement raster (0, 1, 2, or 3 GCMs predict suitability)
# 2. Ensemble-mean suitable raster at MTSPS
# 3. High-confidence suitable raster (3/3 GCM agreement)
# 4. High-confidence transition map relative to the current final full-data
#    MTSPS prediction
# 5. Geodesic area tables for agreement classes, high-confidence habitat,
#    ensemble suitability, and high-confidence gains/losses
# 6. Inter-GCM standard-deviation summaries using Script 12 ensemble SD rasters
# 7. Inventory of existing MaxEnt MESS/clamping/extrapolation rasters, if any
#
# COMMON ANALYSIS DOMAIN
# ----------------------
# To make all scenario-level areas directly comparable, all calculations are
# restricted to one strict common valid domain: grid cells with finite
# continuous suitability in the current full-data map AND every one of the
# 48 individual future GCM projections.
#
# IMPORTANT INTERPRETATION
# ------------------------
# - 3/3 GCM agreement = high-confidence future climatic suitability
# - 2/3 GCM agreement = moderate support, not high-confidence habitat
# - 1/3 GCM agreement = low support
# - 0/3 GCM agreement = consistently unsuitable under the three GCMs
# - This script does NOT rerun MaxEnt, GCM projections, or ensemble averaging.
# =============================================================================

options(warn = 1)

cat("\n")
cat("======================================================================\n")
cat("STARTING SCRIPT 18 — INTER-GCM UNCERTAINTY AND AGREEMENT\n")
cat("Species: Camellia japonica\n")
cat("Primary threshold: MTSPS\n")
cat("======================================================================\n")

script_started_at <- Sys.time()

# =============================================================================
# A. USER SETTINGS
# =============================================================================

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# MTSPS remains the primary threshold throughout this uncertainty analysis.
threshold_id <- "MTSPS"

# The exact continuous ensemble-mean and SD columns in Script 12 manifest.
ensemble_mean_column <- "ensemble_mean_tif"
ensemble_sd_column <- "ensemble_sd_tif"

# The GCM names expected from the completed Script 11 individual projections.
expected_gcms <- c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2")

# FALSE means existing valid agreement/transition maps will be reused.
overwrite_existing_maps <- FALSE

# Produce global uncertainty figures for the two end-century focal scenarios.
focal_figure_scenarios <- data.frame(
  ssp = c("SSP245", "SSP585"),
  period_start = c(2081L, 2081L),
  period_end = c(2100L, 2100L),
  stringsAsFactors = FALSE
)

# =============================================================================
# B. INPUT PATHS
# =============================================================================

final_model_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

individual_projection_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "11_future_projection_individual_GCM_v2_plotsfix"
)

ensemble_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "12_future_ensemble_3GCM"
)

fixed_mtsps_csv <- file.path(
  final_model_root,
  "04_CV_metrics_MTSPS_and_areas",
  "03_FINAL_single_MTSPS_threshold.csv"
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

individual_manifest_candidates <- c(
  file.path(
    individual_projection_root,
    "04_projection_manifests",
    "future_projection_status_manifest.csv"
  ),
  file.path(
    individual_projection_root,
    "04_manifests",
    "future_projection_status_manifest.csv"
  )
)

ensemble_manifest_candidates <- c(
  file.path(
    ensemble_root,
    "04_manifests",
    "ensemble_status_manifest.csv"
  ),
  file.path(
    ensemble_root,
    "04_manifests",
    "future_ensemble_status_manifest.csv"
  )
)

# =============================================================================
# C. OUTPUT PATHS
# =============================================================================

output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "18_interGCM_uncertainty_agreement_extrapolation"
)

status_dir <- file.path(output_root, "00_run_status")
setup_dir <- file.path(output_root, "01_setup")
agreement_dir <- file.path(output_root, "02_GCM_agreement_maps")
high_confidence_dir <- file.path(output_root, "03_high_confidence_maps")
transition_dir <- file.path(output_root, "04_high_confidence_transitions")
area_dir <- file.path(output_root, "05_area_tables")
figure_dir <- file.path(output_root, "06_figures")
support_dir <- file.path(output_root, "07_supporting_files")
log_dir <- file.path(output_root, "08_logs")

started_marker <- file.path(status_dir, "18_STARTED.txt")
completed_marker <- file.path(status_dir, "18_COMPLETED.txt")
failed_marker <- file.path(status_dir, "18_FAILED.txt")
execution_log <- file.path(log_dir, "18_execution_log.txt")

input_manifest_csv <- file.path(setup_dir, "01_input_manifest.csv")
individual_projection_resolved_csv <- file.path(
  setup_dir,
  "02_resolved_individual_projection_rasters.csv"
)
scenario_inventory_csv <- file.path(
  setup_dir,
  "03_scenario_inventory.csv"
)
extrapolation_inventory_csv <- file.path(
  setup_dir,
  "04_existing_MESS_clamping_inventory.csv"
)

common_domain_tif <- file.path(
  support_dir,
  "common_valid_domain_current_plus_all_48_individual_GCM_projections.tif"
)
geodesic_area_tif <- file.path(
  support_dir,
  "global_geodesic_cell_area_km2.tif"
)
agreement_legend_csv <- file.path(
  support_dir,
  "GCM_agreement_legend.csv"
)
high_confidence_transition_legend_csv <- file.path(
  support_dir,
  "high_confidence_transition_legend.csv"
)
agreement_status_manifest_csv <- file.path(
  support_dir,
  "agreement_map_status_manifest.csv"
)
run_summary_csv <- file.path(
  support_dir,
  "uncertainty_analysis_run_summary.csv"
)

agreement_area_long_csv <- file.path(
  area_dir,
  "01_GCM_agreement_area_by_scenario_long.csv"
)
agreement_area_wide_csv <- file.path(
  area_dir,
  "02_GCM_agreement_area_by_scenario_wide.csv"
)
uncertainty_sd_summary_csv <- file.path(
  area_dir,
  "03_inter_GCM_SD_summary_by_scenario.csv"
)
high_confidence_transition_csv <- file.path(
  area_dir,
  "04_high_confidence_transition_area_by_scenario.csv"
)
area_balance_checks_csv <- file.path(
  area_dir,
  "05_high_confidence_area_balance_checks.csv"
)
scenario_confidence_summary_csv <- file.path(
  area_dir,
  "06_scenario_confidence_summary.csv"
)

# =============================================================================
# D. HELPERS
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
  cat(line, "\n", file = execution_log, append = TRUE)
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
    match_index <- grep(candidate, normalized_names, fixed = TRUE)
    if (length(match_index) > 0) return(raw_names[match_index[1]])
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

read_single_numeric_value <- function(path, candidates, label) {
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  column <- find_column(data, candidates, label)

  values <- safe_number(data[[column]])
  values <- values[is.finite(values)]

  if (length(values) != 1L) {
    stop(
      "\nExpected exactly one finite value in ", label, ".\n",
      "Path: ", path,
      call. = FALSE
    )
  }

  values[1]
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
      "All current, individual-GCM, and ensemble rasters must use the same grid.",
      call. = FALSE
    )
  }
}

find_existing_raster_path_in_row <- function(row_data, label) {
  # 1. Direct file fields in manifest row.
  candidate_values <- unlist(row_data, use.names = FALSE)
  candidate_values <- as.character(candidate_values)
  candidate_values <- candidate_values[
    !is.na(candidate_values) &
      nzchar(candidate_values) &
      file.exists(candidate_values)
  ]

  direct_rasters <- candidate_values[
    grepl("\\.(tif|tiff|asc)$", candidate_values, ignore.case = TRUE) &
      grepl("cloglog", basename(candidate_values), ignore.case = TRUE) &
      !grepl("mess|clamp", basename(candidate_values), ignore.case = TRUE)
  ]

  direct_rasters <- unique(direct_rasters)

  if (length(direct_rasters) == 1L) {
    return(normalizePath(direct_rasters[1], winslash = "/", mustWork = TRUE))
  }

  # 2. Search existing directory fields in manifest row.
  candidate_dirs <- candidate_values[dir.exists(candidate_values)]
  candidate_dirs <- unique(candidate_dirs)

  searched_rasters <- character(0)

  for (directory in candidate_dirs) {
    files <- list.files(
      directory,
      pattern = "\\.(tif|tiff|asc)$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    files <- files[
      grepl("cloglog", basename(files), ignore.case = TRUE) &
        !grepl("mess|clamp", basename(files), ignore.case = TRUE)
    ]

    searched_rasters <- c(searched_rasters, files)
  }

  searched_rasters <- unique(searched_rasters)

  if (length(searched_rasters) == 1L) {
    return(
      normalizePath(searched_rasters[1], winslash = "/", mustWork = TRUE)
    )
  }

  stop(
    "\nCould not uniquely resolve the individual-GCM cloglog raster for:\n",
    label,
    "\nDirect cloglog candidates found: ",
    length(direct_rasters),
    "\nDirectory-scanned cloglog candidates found: ",
    length(searched_rasters),
    "\nInspect the individual projection manifest.",
    call. = FALSE
  )
}

calculate_area_by_code <- function(class_raster, cell_area_raster, domain_raster,
                                   codes, class_names) {
  assert_same_grid(cell_area_raster, class_raster, "cell-area raster")
  assert_same_grid(domain_raster, class_raster, "common domain")

  class_in_domain <- terra::mask(class_raster, domain_raster)

  zonal_table <- terra::zonal(
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

  if (!is.null(zonal_table) && nrow(zonal_table) > 0) {
    zone_codes <- suppressWarnings(
      as.integer(round(as.numeric(zonal_table[[1]])))
    )
    zone_areas <- as.numeric(zonal_table[[ncol(zonal_table)]])

    for (i in seq_along(codes)) {
      matching_values <- zone_areas[zone_codes == codes[i]]
      if (length(matching_values) > 0) {
        output$area_km2[i] <- sum(matching_values, na.rm = TRUE)
      }
    }
  }

  output$area_10_4_km2 <- output$area_km2 / 10000
  output
}

get_code_area <- function(area_table, code) {
  output <- area_table$area_km2[area_table$code == code]

  if (length(output) != 1L) {
    stop("Cannot retrieve unique area for code ", code, call. = FALSE)
  }

  output
}

get_global_stat <- function(raster_object, statistic_name) {
  value <- terra::global(
    raster_object,
    fun = statistic_name,
    na.rm = TRUE
  )

  as.numeric(value[1, 1])
}

mean_sd_by_agreement <- function(sd_raster, agreement_raster, domain_raster) {
  sd_in_domain <- terra::mask(sd_raster, domain_raster)
  agreement_in_domain <- terra::mask(agreement_raster, domain_raster)

  zonal_table <- terra::zonal(
    sd_in_domain,
    agreement_in_domain,
    fun = "mean",
    na.rm = TRUE
  )

  output <- data.frame(
    agreement_count = 0:3,
    mean_inter_GCM_SD = NA_real_,
    stringsAsFactors = FALSE
  )

  if (!is.null(zonal_table) && nrow(zonal_table) > 0) {
    agreement_codes <- suppressWarnings(
      as.integer(round(as.numeric(zonal_table[[1]])))
    )
    mean_values <- as.numeric(zonal_table[[ncol(zonal_table)]])

    for (i in seq_len(nrow(output))) {
      hit <- mean_values[agreement_codes == output$agreement_count[i]]
      if (length(hit) > 0) output$mean_inter_GCM_SD[i] <- hit[1]
    }
  }

  output
}

create_global_uncertainty_figure <- function(
  agreement_raster,
  sd_raster,
  high_confidence_raster,
  scenario_label,
  output_pdf
) {
  grDevices::pdf(output_pdf, width = 13, height = 5)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(1, 3), mar = c(2.8, 2.8, 2.6, 4.3))

  plot(
    agreement_raster,
    col = c("#F0F0F0", "#FEE08B", "#FDAE61", "#1A9850"),
    breaks = seq(-0.5, 3.5, by = 1),
    main = paste0(scenario_label, "\nGCM suitable agreement"),
    axes = TRUE,
    plg = list(
      title = "Suitable\nGCMs",
      at = 0:3,
      labels = c("0", "1", "2", "3")
    )
  )

  plot(
    sd_raster,
    main = paste0(scenario_label, "\nInter-GCM SD"),
    axes = TRUE,
    plg = list(title = "SD")
  )

  plot(
    high_confidence_raster,
    col = c("#F0F0F0", "#2166AC"),
    breaks = c(-0.5, 0.5, 1.5),
    main = paste0(scenario_label, "\nHigh-confidence suitability (3/3)"),
    axes = TRUE,
    plg = list(
      title = "Class",
      at = c(0, 1),
      labels = c("No", "Yes")
    )
  )
}

# =============================================================================
# E. MAIN ANALYSIS
# =============================================================================

run_inter_gcm_uncertainty_analysis <- function() {

  for (folder in c(
    output_root, status_dir, setup_dir, agreement_dir, high_confidence_dir,
    transition_dir, area_dir, figure_dir, support_dir, log_dir
  )) {
    ensure_dir(folder)
  }

  unlink(completed_marker, force = TRUE)
  unlink(failed_marker, force = TRUE)

  writeLines(
    c(
      "Script 18 started",
      paste("Started:", format(script_started_at, "%Y-%m-%d %H:%M:%S")),
      "Purpose: inter-GCM agreement and uncertainty analysis at MTSPS"
    ),
    started_marker
  )

  if (!requireNamespace("terra", quietly = TRUE)) {
    install.packages("terra", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(terra))

  write_log("Package terra loaded.")

  # -------------------------------------------------------------------
  # E1. Inputs and MTSPS
  # -------------------------------------------------------------------

  current_continuous_path <- find_first_existing_file(
    current_continuous_candidates,
    "current final full-data continuous suitability raster"
  )

  individual_manifest_path <- find_first_existing_file(
    individual_manifest_candidates,
    "individual-GCM projection manifest"
  )

  ensemble_manifest_path <- find_first_existing_file(
    ensemble_manifest_candidates,
    "three-GCM ensemble manifest"
  )

  if (!file.exists(fixed_mtsps_csv)) {
    stop(
      "\nMissing fixed MTSPS threshold CSV:\n",
      fixed_mtsps_csv,
      call. = FALSE
    )
  }

  mtsps_threshold <- read_single_numeric_value(
    fixed_mtsps_csv,
    c("MTSPS_threshold", "threshold"),
    "fixed MTSPS threshold CSV"
  )

  current_continuous <- terra::rast(current_continuous_path)

  individual_manifest <- read.csv(
    individual_manifest_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  ensemble_manifest <- read.csv(
    ensemble_manifest_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(
      input_type = c(
        "Current final full-data continuous suitability",
        "Fixed primary MTSPS threshold",
        "Individual-GCM projection manifest",
        "Three-GCM ensemble manifest"
      ),
      path = c(
        current_continuous_path,
        normalizePath(fixed_mtsps_csv, winslash = "/", mustWork = TRUE),
        individual_manifest_path,
        ensemble_manifest_path
      ),
      stringsAsFactors = FALSE
    ),
    input_manifest_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E2. Identify fields in individual-GCM and ensemble manifests
  # -------------------------------------------------------------------

  ind_status_col <- find_column(
    individual_manifest,
    c("status"),
    "individual-GCM projection manifest"
  )
  ind_gcm_col <- find_column(
    individual_manifest,
    c("gcm", "gcm_name", "model"),
    "individual-GCM projection manifest"
  )
  ind_ssp_col <- find_column(
    individual_manifest,
    c("ssp"),
    "individual-GCM projection manifest"
  )
  ind_start_col <- find_column(
    individual_manifest,
    c("period_start", "periodstart"),
    "individual-GCM projection manifest"
  )
  ind_end_col <- find_column(
    individual_manifest,
    c("period_end", "periodend"),
    "individual-GCM projection manifest"
  )

  ens_status_col <- find_column(
    ensemble_manifest,
    c("status"),
    "three-GCM ensemble manifest"
  )
  ens_ssp_col <- find_column(
    ensemble_manifest,
    c("ssp"),
    "three-GCM ensemble manifest"
  )
  ens_start_col <- find_column(
    ensemble_manifest,
    c("period_start", "periodstart"),
    "three-GCM ensemble manifest"
  )
  ens_end_col <- find_column(
    ensemble_manifest,
    c("period_end", "periodend"),
    "three-GCM ensemble manifest"
  )

  if (!ensemble_mean_column %in% names(ensemble_manifest)) {
    stop(
      "\nMissing expected ensemble mean column: ",
      ensemble_mean_column,
      "\nAvailable columns: ", paste(names(ensemble_manifest), collapse = ", "),
      call. = FALSE
    )
  }

  if (!ensemble_sd_column %in% names(ensemble_manifest)) {
    stop(
      "\nMissing expected ensemble SD column: ",
      ensemble_sd_column,
      "\nAvailable columns: ", paste(names(ensemble_manifest), collapse = ", "),
      call. = FALSE
    )
  }

  valid_status <- c("completed", "skipped_existing_output")

  # -------------------------------------------------------------------
  # E3. Resolve 48 individual GCM cloglog rasters
  # -------------------------------------------------------------------

  individual_rows <- individual_manifest[
    individual_manifest[[ind_status_col]] %in% valid_status,
    ,
    drop = FALSE
  ]

  individual_rows$gcm_analysis <- as.character(individual_rows[[ind_gcm_col]])
  individual_rows$ssp_analysis <- as.character(individual_rows[[ind_ssp_col]])
  individual_rows$period_start_analysis <- safe_number(individual_rows[[ind_start_col]])
  individual_rows$period_end_analysis <- safe_number(individual_rows[[ind_end_col]])

  individual_rows <- individual_rows[
    individual_rows$gcm_analysis %in% expected_gcms &
      is.finite(individual_rows$period_start_analysis) &
      is.finite(individual_rows$period_end_analysis),
    ,
    drop = FALSE
  ]

  if (nrow(individual_rows) != 48L) {
    stop(
      "\nExpected 48 completed individual GCM projections (3 GCM × 16 scenarios), ",
      "but found ", nrow(individual_rows), ".\n",
      "Check: ", individual_manifest_path,
      call. = FALSE
    )
  }

  write_log("Resolving 48 individual cloglog suitability rasters.")

  individual_rows$cloglog_path_analysis <- vapply(
    seq_len(nrow(individual_rows)),
    function(i) {
      row <- individual_rows[i, , drop = FALSE]

      find_existing_raster_path_in_row(
        row,
        paste0(
          row$gcm_analysis,
          " | ",
          row$ssp_analysis,
          " ",
          row$period_start_analysis,
          "-",
          row$period_end_analysis
        )
      )
    },
    character(1)
  )

  duplicate_keys <- duplicated(
    paste(
      individual_rows$gcm_analysis,
      individual_rows$ssp_analysis,
      individual_rows$period_start_analysis,
      individual_rows$period_end_analysis,
      sep = "__"
    )
  )

  if (any(duplicate_keys)) {
    stop(
      "\nDuplicate GCM × scenario projection records found after resolving rasters.",
      call. = FALSE
    )
  }

  individual_rows <- individual_rows[
    order(
      individual_rows$ssp_analysis,
      individual_rows$period_start_analysis,
      individual_rows$gcm_analysis
    ),
    ,
    drop = FALSE
  ]

  write.csv(
    individual_rows[
      c(
        "gcm_analysis",
        "ssp_analysis",
        "period_start_analysis",
        "period_end_analysis",
        "cloglog_path_analysis"
      )
    ],
    individual_projection_resolved_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E4. Resolve 16 ensemble mean and SD rasters
  # -------------------------------------------------------------------

  ensemble_rows <- ensemble_manifest[
    ensemble_manifest[[ens_status_col]] %in% valid_status &
      !is.na(ensemble_manifest[[ensemble_mean_column]]) &
      nzchar(ensemble_manifest[[ensemble_mean_column]]) &
      file.exists(ensemble_manifest[[ensemble_mean_column]]) &
      !is.na(ensemble_manifest[[ensemble_sd_column]]) &
      nzchar(ensemble_manifest[[ensemble_sd_column]]) &
      file.exists(ensemble_manifest[[ensemble_sd_column]]),
    ,
    drop = FALSE
  ]

  ensemble_rows$ssp_analysis <- as.character(ensemble_rows[[ens_ssp_col]])
  ensemble_rows$period_start_analysis <- safe_number(ensemble_rows[[ens_start_col]])
  ensemble_rows$period_end_analysis <- safe_number(ensemble_rows[[ens_end_col]])
  ensemble_rows$ensemble_mean_path_analysis <- normalizePath(
    as.character(ensemble_rows[[ensemble_mean_column]]),
    winslash = "/",
    mustWork = TRUE
  )
  ensemble_rows$ensemble_sd_path_analysis <- normalizePath(
    as.character(ensemble_rows[[ensemble_sd_column]]),
    winslash = "/",
    mustWork = TRUE
  )

  if (nrow(ensemble_rows) != 16L) {
    stop(
      "\nExpected 16 completed ensemble mean + SD rasters, but found ",
      nrow(ensemble_rows), ".",
      call. = FALSE
    )
  }

  ensemble_rows <- ensemble_rows[
    order(
      ensemble_rows$ssp_analysis,
      ensemble_rows$period_start_analysis
    ),
    ,
    drop = FALSE
  ]

  # Confirm every scenario has exactly 3 individual GCM maps and one ensemble.
  scenario_keys <- unique(
    ensemble_rows[
      c("ssp_analysis", "period_start_analysis", "period_end_analysis"),
      drop = FALSE
    ]
  )

  scenario_keys <- scenario_keys[
    order(scenario_keys$ssp_analysis, scenario_keys$period_start_analysis),
    ,
    drop = FALSE
  ]

  if (nrow(scenario_keys) != 16L) {
    stop("\nEnsemble manifest does not contain 16 unique scenarios.", call. = FALSE)
  }

  scenario_inventory <- do.call(
    rbind,
    lapply(seq_len(nrow(scenario_keys)), function(i) {
      key <- scenario_keys[i, , drop = FALSE]

      individual_count <- sum(
        individual_rows$ssp_analysis == key$ssp_analysis &
          individual_rows$period_start_analysis == key$period_start_analysis &
          individual_rows$period_end_analysis == key$period_end_analysis
      )

      data.frame(
        ssp = key$ssp_analysis,
        period_start = key$period_start_analysis,
        period_end = key$period_end_analysis,
        individual_GCM_count = individual_count,
        stringsAsFactors = FALSE
      )
    })
  )

  if (any(scenario_inventory$individual_GCM_count != 3L)) {
    print(scenario_inventory)
    stop(
      "\nAt least one future scenario does not have exactly 3 individual GCM projections.",
      call. = FALSE
    )
  }

  write.csv(scenario_inventory, scenario_inventory_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # E5. Build strict common valid domain across current + all 48 individual
  #     GCM future suitability maps
  # -------------------------------------------------------------------

  write_log(
    "Building strict common valid domain across current + all 48 individual GCM projections."
  )

  common_domain <- terra::ifel(is.na(current_continuous), NA, 1)

  for (i in seq_len(nrow(individual_rows))) {
    gcm_raster <- terra::rast(individual_rows$cloglog_path_analysis[i])

    assert_same_grid(
      current_continuous,
      gcm_raster,
      paste0(
        individual_rows$gcm_analysis[i], " | ",
        individual_rows$ssp_analysis[i], " ",
        individual_rows$period_start_analysis[i], "-",
        individual_rows$period_end_analysis[i]
      )
    )

    common_domain <- terra::mask(common_domain, gcm_raster)
  }

  names(common_domain) <- "strict_common_valid_domain"

  terra::writeRaster(
    common_domain,
    common_domain_tif,
    overwrite = TRUE,
    datatype = "INT1U",
    wopt = list(gdal = "COMPRESS=LZW")
  )

  # Cell area on the full grid; all later area calculations are masked to domain.
  if (file.exists(geodesic_area_tif)) {
    cell_area <- terra::rast(geodesic_area_tif)
    assert_same_grid(current_continuous, cell_area, "existing geodesic area raster")
  } else {
    write_log("Creating latitude-corrected geodesic cell-area raster (km²).")

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

  common_domain_area_table <- calculate_area_by_code(
    class_raster = common_domain,
    cell_area_raster = cell_area,
    domain_raster = common_domain,
    codes = 1,
    class_names = "Strict common valid domain"
  )

  strict_common_domain_km2 <- get_code_area(common_domain_area_table, 1)

  current_binary <- terra::ifel(
    is.na(common_domain),
    NA,
    terra::ifel(current_continuous >= mtsps_threshold, 1, 0)
  )

  # -------------------------------------------------------------------
  # E6. Legends and extrapolation inventory
  # -------------------------------------------------------------------

  write.csv(
    data.frame(
      agreement_count = 0:3,
      agreement_class = c(
        "Consistently unsuitable",
        "Low support",
        "Moderate support",
        "High-confidence suitable"
      ),
      definition = c(
        "0 of 3 GCMs predict suitability at MTSPS",
        "1 of 3 GCMs predicts suitability at MTSPS",
        "2 of 3 GCMs predict suitability at MTSPS",
        "3 of 3 GCMs predict suitability at MTSPS"
      ),
      stringsAsFactors = FALSE
    ),
    agreement_legend_csv,
    row.names = FALSE
  )

  write.csv(
    data.frame(
      transition_code = 1:4,
      transition_class = c(
        "Stable not high-confidence suitable",
        "High-confidence gain",
        "Loss of high-confidence suitability",
        "Stable high-confidence suitable"
      ),
      definition = c(
        "Current unsuitable -> future not 3/3 GCM suitable",
        "Current unsuitable -> future 3/3 GCM suitable",
        "Current suitable -> future not 3/3 GCM suitable",
        "Current suitable -> future 3/3 GCM suitable"
      ),
      stringsAsFactors = FALSE
    ),
    high_confidence_transition_legend_csv,
    row.names = FALSE
  )

  extrapolation_files <- list.files(
    individual_projection_root,
    pattern = "mess|clamp|extrapol",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  extrapolation_files <- extrapolation_files[
    grepl("\\.(tif|tiff|asc)$", extrapolation_files, ignore.case = TRUE)
  ]

  extrapolation_inventory <- if (length(extrapolation_files) > 0) {
    data.frame(
      file_path = normalizePath(
        extrapolation_files,
        winslash = "/",
        mustWork = TRUE
      ),
      file_name = basename(extrapolation_files),
      detected_type = ifelse(
        grepl("mess", basename(extrapolation_files), ignore.case = TRUE),
        "MESS",
        ifelse(
          grepl("clamp", basename(extrapolation_files), ignore.case = TRUE),
          "Clamping",
          "Other extrapolation-related output"
        )
      ),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      file_path = NA_character_,
      file_name = NA_character_,
      detected_type = "No existing MESS/clamping/extrapolation raster found in Script 11 output tree",
      stringsAsFactors = FALSE
    )
  }

  write.csv(
    extrapolation_inventory,
    extrapolation_inventory_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E7. Create agreement, high-confidence and transition maps
  # -------------------------------------------------------------------

  agreement_area_records <- list()
  sd_summary_records <- list()
  transition_area_records <- list()
  status_records <- list()

  for (i in seq_len(nrow(scenario_keys))) {
    key <- scenario_keys[i, , drop = FALSE]

    ssp <- as.character(key$ssp_analysis)
    period_start <- as.integer(key$period_start_analysis)
    period_end <- as.integer(key$period_end_analysis)
    scenario_label <- paste0(ssp, " ", period_start, "-", period_end)

    cat(
      "\n[Scenario ", i, "/", nrow(scenario_keys), "] ",
      scenario_label, "\n",
      sep = ""
    )

    scenario_individual <- individual_rows[
      individual_rows$ssp_analysis == ssp &
        individual_rows$period_start_analysis == period_start &
        individual_rows$period_end_analysis == period_end,
      ,
      drop = FALSE
    ]

    scenario_individual <- scenario_individual[
      match(expected_gcms, scenario_individual$gcm_analysis),
      ,
      drop = FALSE
    ]

    if (
      nrow(scenario_individual) != 3L ||
      any(is.na(scenario_individual$gcm_analysis))
    ) {
      stop(
        "\nCould not order exactly three expected GCMs for ",
        scenario_label,
        call. = FALSE
      )
    }

    ensemble_row <- ensemble_rows[
      ensemble_rows$ssp_analysis == ssp &
        ensemble_rows$period_start_analysis == period_start &
        ensemble_rows$period_end_analysis == period_end,
      ,
      drop = FALSE
    ]

    if (nrow(ensemble_row) != 1L) {
      stop(
        "\nCould not uniquely match ensemble rasters for ",
        scenario_label,
        call. = FALSE
      )
    }

    agreement_path <- file.path(
      agreement_dir,
      paste0(
        species_name,
        "_GCM_suitable_agreement_",
        ssp, "_", period_start, "_", period_end,
        "_MTSPS.tif"
      )
    )

    ensemble_binary_path <- file.path(
      agreement_dir,
      paste0(
        species_name,
        "_ensemble_mean_suitable_",
        ssp, "_", period_start, "_", period_end,
        "_MTSPS.tif"
      )
    )

    high_confidence_path <- file.path(
      high_confidence_dir,
      paste0(
        species_name,
        "_high_confidence_3of3_GCM_suitable_",
        ssp, "_", period_start, "_", period_end,
        "_MTSPS.tif"
      )
    )

    transition_path <- file.path(
      transition_dir,
      paste0(
        species_name,
        "_high_confidence_transition_current_to_",
        ssp, "_", period_start, "_", period_end,
        "_MTSPS.tif"
      )
    )

    output_exists <- all(file.exists(
      c(
        agreement_path,
        ensemble_binary_path,
        high_confidence_path,
        transition_path
      )
    ))

    if (!output_exists || overwrite_existing_maps) {
      write_log("Building agreement and confidence maps: ", scenario_label)

      gcm_rasters <- lapply(
        scenario_individual$cloglog_path_analysis,
        terra::rast
      )

      for (j in seq_along(gcm_rasters)) {
        assert_same_grid(
          current_continuous,
          gcm_rasters[[j]],
          paste0(
            scenario_individual$gcm_analysis[j],
            " | ", scenario_label
          )
        )
      }

      # At the strict common domain, all 3 future GCM rasters are finite.
      gcm_binary <- lapply(
        gcm_rasters,
        function(r) {
          terra::ifel(
            is.na(common_domain),
            NA,
            terra::ifel(r >= mtsps_threshold, 1, 0)
          )
        }
      )

      agreement <- gcm_binary[[1]] + gcm_binary[[2]] + gcm_binary[[3]]
      agreement <- terra::mask(agreement, common_domain)
      names(agreement) <- "suitable_GCM_agreement_count"

      ensemble_mean <- terra::rast(ensemble_row$ensemble_mean_path_analysis)
      ensemble_sd <- terra::rast(ensemble_row$ensemble_sd_path_analysis)

      assert_same_grid(current_continuous, ensemble_mean, paste0("Ensemble mean ", scenario_label))
      assert_same_grid(current_continuous, ensemble_sd, paste0("Ensemble SD ", scenario_label))

      ensemble_binary <- terra::ifel(
        is.na(common_domain),
        NA,
        terra::ifel(ensemble_mean >= mtsps_threshold, 1, 0)
      )
      names(ensemble_binary) <- "ensemble_mean_suitable_MTSPS"

      high_confidence <- terra::ifel(
        is.na(common_domain),
        NA,
        terra::ifel(agreement == 3, 1, 0)
      )
      names(high_confidence) <- "high_confidence_3of3_GCM_suitable"

      high_conf_transition <- 1 + high_confidence + 2 * current_binary
      high_conf_transition <- terra::mask(
        high_conf_transition,
        common_domain
      )
      names(high_conf_transition) <- "high_confidence_transition"

      terra::writeRaster(
        agreement,
        agreement_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      terra::writeRaster(
        ensemble_binary,
        ensemble_binary_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      terra::writeRaster(
        high_confidence,
        high_confidence_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      terra::writeRaster(
        high_conf_transition,
        transition_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      map_status <- "completed"
    } else {
      write_log("Reusing existing agreement/confidence maps: ", scenario_label)
      map_status <- "skipped_existing_output"

      agreement <- terra::rast(agreement_path)
      ensemble_binary <- terra::rast(ensemble_binary_path)
      high_confidence <- terra::rast(high_confidence_path)
      high_conf_transition <- terra::rast(transition_path)

      ensemble_sd <- terra::rast(ensemble_row$ensemble_sd_path_analysis)
    }

    # Ensure objects are available when maps were created in this iteration.
    if (!exists("agreement") || !inherits(agreement, "SpatRaster")) {
      agreement <- terra::rast(agreement_path)
    }
    if (!exists("ensemble_binary") || !inherits(ensemble_binary, "SpatRaster")) {
      ensemble_binary <- terra::rast(ensemble_binary_path)
    }
    if (!exists("high_confidence") || !inherits(high_confidence, "SpatRaster")) {
      high_confidence <- terra::rast(high_confidence_path)
    }
    if (!exists("high_conf_transition") || !inherits(high_conf_transition, "SpatRaster")) {
      high_conf_transition <- terra::rast(transition_path)
    }
    if (!exists("ensemble_sd") || !inherits(ensemble_sd, "SpatRaster")) {
      ensemble_sd <- terra::rast(ensemble_row$ensemble_sd_path_analysis)
    }

    # Agreement areas.
    agreement_area <- calculate_area_by_code(
      class_raster = agreement,
      cell_area_raster = cell_area,
      domain_raster = common_domain,
      codes = 0:3,
      class_names = c(
        "Consistently unsuitable (0/3)",
        "Low support (1/3)",
        "Moderate support (2/3)",
        "High-confidence suitable (3/3)"
      )
    )

    ensemble_binary_area <- calculate_area_by_code(
      class_raster = ensemble_binary,
      cell_area_raster = cell_area,
      domain_raster = common_domain,
      codes = c(0, 1),
      class_names = c("Ensemble unsuitable", "Ensemble suitable")
    )

    high_confidence_area <- get_code_area(agreement_area, 3)
    ensemble_suitable_area <- get_code_area(ensemble_binary_area, 1)

    agreement_area_records[[length(agreement_area_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = period_start,
      period_end = period_end,
      scenario = scenario_label,
      MTSPS_threshold = mtsps_threshold,
      strict_common_valid_domain_km2 = strict_common_domain_km2,
      agreement_0of3_km2 = get_code_area(agreement_area, 0),
      agreement_1of3_km2 = get_code_area(agreement_area, 1),
      agreement_2of3_km2 = get_code_area(agreement_area, 2),
      agreement_3of3_km2 = high_confidence_area,
      agreement_3of3_10_4_km2 = high_confidence_area / 10000,
      ensemble_mean_suitable_km2 = ensemble_suitable_area,
      ensemble_suitable_but_not_3of3_km2 =
        max(0, ensemble_suitable_area - high_confidence_area),
      percent_common_domain_3of3 =
        100 * high_confidence_area / strict_common_domain_km2,
      percent_ensemble_suitable_with_3of3_agreement =
        if (ensemble_suitable_area > 0) {
          100 * high_confidence_area / ensemble_suitable_area
        } else {
          NA_real_
        },
      agreement_tif = normalizePath(agreement_path, winslash = "/", mustWork = TRUE),
      high_confidence_tif = normalizePath(high_confidence_path, winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )

    # Inter-GCM SD summaries.
    sd_in_domain <- terra::mask(ensemble_sd, common_domain)

    mean_sd_classes <- mean_sd_by_agreement(
      sd_raster = ensemble_sd,
      agreement_raster = agreement,
      domain_raster = common_domain
    )

    sd_summary_records[[length(sd_summary_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = period_start,
      period_end = period_end,
      scenario = scenario_label,
      MTSPS_threshold = mtsps_threshold,
      inter_GCM_SD_mean = get_global_stat(sd_in_domain, "mean"),
      inter_GCM_SD_sd = get_global_stat(sd_in_domain, "sd"),
      inter_GCM_SD_min = get_global_stat(sd_in_domain, "min"),
      inter_GCM_SD_max = get_global_stat(sd_in_domain, "max"),
      mean_SD_agreement_0of3 = mean_sd_classes$mean_inter_GCM_SD[
        mean_sd_classes$agreement_count == 0
      ],
      mean_SD_agreement_1of3 = mean_sd_classes$mean_inter_GCM_SD[
        mean_sd_classes$agreement_count == 1
      ],
      mean_SD_agreement_2of3 = mean_sd_classes$mean_inter_GCM_SD[
        mean_sd_classes$agreement_count == 2
      ],
      mean_SD_agreement_3of3 = mean_sd_classes$mean_inter_GCM_SD[
        mean_sd_classes$agreement_count == 3
      ],
      ensemble_SD_tif = normalizePath(
        ensemble_row$ensemble_sd_path_analysis,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )

    # High-confidence transition areas.
    transition_area <- calculate_area_by_code(
      class_raster = high_conf_transition,
      cell_area_raster = cell_area,
      domain_raster = common_domain,
      codes = 1:4,
      class_names = c(
        "Stable not high-confidence suitable",
        "High-confidence gain",
        "Loss of high-confidence suitability",
        "Stable high-confidence suitable"
      )
    )

    stable_not_high_conf <- get_code_area(transition_area, 1)
    high_conf_gain <- get_code_area(transition_area, 2)
    loss_high_conf <- get_code_area(transition_area, 3)
    stable_high_conf <- get_code_area(transition_area, 4)

    current_suitable_area <- loss_high_conf + stable_high_conf
    future_high_conf_area <- high_conf_gain + stable_high_conf

    transition_area_records[[length(transition_area_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = period_start,
      period_end = period_end,
      scenario = scenario_label,
      MTSPS_threshold = mtsps_threshold,
      strict_common_valid_domain_km2 = strict_common_domain_km2,
      stable_not_high_confidence_suitable_km2 = stable_not_high_conf,
      high_confidence_gain_km2 = high_conf_gain,
      loss_of_high_confidence_suitability_km2 = loss_high_conf,
      stable_high_confidence_suitable_km2 = stable_high_conf,
      current_MTSPS_suitable_km2 = current_suitable_area,
      future_high_confidence_suitable_km2 = future_high_conf_area,
      net_high_confidence_change_km2 = future_high_conf_area - current_suitable_area,
      percent_future_high_confidence_relative_to_current_suitable =
        100 * future_high_conf_area / current_suitable_area,
      transition_tif = normalizePath(
        transition_path,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )

    status_records[[length(status_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = period_start,
      period_end = period_end,
      scenario = scenario_label,
      status = map_status,
      agreement_tif = normalizePath(agreement_path, winslash = "/", mustWork = TRUE),
      ensemble_binary_tif = normalizePath(ensemble_binary_path, winslash = "/", mustWork = TRUE),
      high_confidence_tif = normalizePath(high_confidence_path, winslash = "/", mustWork = TRUE),
      high_confidence_transition_tif = normalizePath(transition_path, winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )

    # Focal scenario figures.
    is_focal <- any(
      focal_figure_scenarios$ssp == ssp &
        focal_figure_scenarios$period_start == period_start &
        focal_figure_scenarios$period_end == period_end
    )

    if (is_focal) {
      figure_path <- file.path(
        figure_dir,
        paste0(
          "Figure_S_inter_GCM_agreement_SD_",
          ssp, "_", period_start, "_", period_end,
          ".pdf"
        )
      )

      tryCatch(
        {
          create_global_uncertainty_figure(
            agreement_raster = agreement,
            sd_raster = sd_in_domain,
            high_confidence_raster = high_confidence,
            scenario_label = scenario_label,
            output_pdf = figure_path
          )
          write_log("Created focal uncertainty figure: ", figure_path)
        },
        error = function(e) {
          write_log(
            "Focal figure skipped for ", scenario_label,
            " because: ", conditionMessage(e)
          )
        }
      )
    }

    # Remove large objects explicitly before next global scenario.
    rm(
      agreement,
      ensemble_binary,
      high_confidence,
      high_conf_transition,
      ensemble_sd,
      sd_in_domain
    )
    gc()
  }

  # -------------------------------------------------------------------
  # E8. Tables, balances, and final summaries
  # -------------------------------------------------------------------

  agreement_area_long <- do.call(rbind, agreement_area_records)
  uncertainty_sd_summary <- do.call(rbind, sd_summary_records)
  high_conf_transition <- do.call(rbind, transition_area_records)
  agreement_status_manifest <- do.call(rbind, status_records)

  agreement_area_long <- agreement_area_long[
    order(agreement_area_long$ssp, agreement_area_long$period_start),
    ,
    drop = FALSE
  ]

  uncertainty_sd_summary <- uncertainty_sd_summary[
    order(uncertainty_sd_summary$ssp, uncertainty_sd_summary$period_start),
    ,
    drop = FALSE
  ]

  high_conf_transition <- high_conf_transition[
    order(high_conf_transition$ssp, high_conf_transition$period_start),
    ,
    drop = FALSE
  ]

  agreement_status_manifest <- agreement_status_manifest[
    order(agreement_status_manifest$ssp, agreement_status_manifest$period_start),
    ,
    drop = FALSE
  ]

  agreement_area_wide <- agreement_area_long[
    c(
      "scenario",
      "ssp",
      "period_start",
      "period_end",
      "MTSPS_threshold",
      "agreement_0of3_km2",
      "agreement_1of3_km2",
      "agreement_2of3_km2",
      "agreement_3of3_km2",
      "ensemble_mean_suitable_km2",
      "ensemble_suitable_but_not_3of3_km2",
      "percent_common_domain_3of3",
      "percent_ensemble_suitable_with_3of3_agreement"
    )
  ]

  area_balance_checks <- data.frame(
    scenario = high_conf_transition$scenario,
    ssp = high_conf_transition$ssp,
    period_start = high_conf_transition$period_start,
    period_end = high_conf_transition$period_end,
    current_suitable_balance_error_km2 =
      high_conf_transition$current_MTSPS_suitable_km2 -
      (
        high_conf_transition$loss_of_high_confidence_suitability_km2 +
          high_conf_transition$stable_high_confidence_suitable_km2
      ),
    future_high_confidence_balance_error_km2 =
      high_conf_transition$future_high_confidence_suitable_km2 -
      (
        high_conf_transition$high_confidence_gain_km2 +
          high_conf_transition$stable_high_confidence_suitable_km2
      ),
    stringsAsFactors = FALSE
  )

  max_balance_error <- max(
    abs(
      unlist(
        area_balance_checks[
          c(
            "current_suitable_balance_error_km2",
            "future_high_confidence_balance_error_km2"
          )
        ]
      )
    ),
    na.rm = TRUE
  )

  scenario_confidence_summary <- merge(
    agreement_area_long,
    uncertainty_sd_summary,
    by = c("scenario", "ssp", "period_start", "period_end", "MTSPS_threshold"),
    all.x = TRUE,
    sort = FALSE
  )

  scenario_confidence_summary <- merge(
    scenario_confidence_summary,
    high_conf_transition[
      c(
        "scenario",
        "ssp",
        "period_start",
        "period_end",
        "MTSPS_threshold",
        "high_confidence_gain_km2",
        "loss_of_high_confidence_suitability_km2",
        "stable_high_confidence_suitable_km2",
        "future_high_confidence_suitable_km2"
      )
    ],
    by = c("scenario", "ssp", "period_start", "period_end", "MTSPS_threshold"),
    all.x = TRUE,
    sort = FALSE
  )

  write.csv(agreement_area_long, agreement_area_long_csv, row.names = FALSE)
  write.csv(agreement_area_wide, agreement_area_wide_csv, row.names = FALSE)
  write.csv(uncertainty_sd_summary, uncertainty_sd_summary_csv, row.names = FALSE)
  write.csv(high_conf_transition, high_confidence_transition_csv, row.names = FALSE)
  write.csv(area_balance_checks, area_balance_checks_csv, row.names = FALSE)
  write.csv(
    scenario_confidence_summary,
    scenario_confidence_summary_csv,
    row.names = FALSE
  )
  write.csv(
    agreement_status_manifest,
    agreement_status_manifest_csv,
    row.names = FALSE
  )

  run_summary <- data.frame(
    species = species_name,
    MTSPS_threshold = mtsps_threshold,
    individual_GCMs = paste(expected_gcms, collapse = "; "),
    individual_future_rasters = nrow(individual_rows),
    ensemble_scenarios = nrow(scenario_keys),
    strict_common_valid_domain_km2 = strict_common_domain_km2,
    agreement_maps_processed = nrow(agreement_status_manifest),
    maximum_high_confidence_area_balance_error_km2 = max_balance_error,
    existing_MESS_or_clamping_rasters_detected = sum(
      !is.na(extrapolation_inventory$file_path)
    ),
    stringsAsFactors = FALSE
  )

  write.csv(run_summary, run_summary_csv, row.names = FALSE)

  writeLines(
    c(
      "SCRIPT 18 COMPLETED SUCCESSFULLY",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("MTSPS threshold:", format(mtsps_threshold, scientific = TRUE, digits = 8)),
      paste("Individual GCM projections resolved:", nrow(individual_rows)),
      paste("Scenarios processed:", nrow(agreement_status_manifest)),
      paste("Strict common valid domain (km2):", format(strict_common_domain_km2, scientific = FALSE)),
      paste("Maximum high-confidence area-balance error (km2):", format(max_balance_error, scientific = TRUE, digits = 8)),
      paste("Existing MESS/clamping rasters detected:", sum(!is.na(extrapolation_inventory$file_path)))
    ),
    completed_marker
  )

  write_log("Script completed successfully.")

  cat("\n======================================================================\n")
  cat("SCRIPT 18 COMPLETED — INTER-GCM UNCERTAINTY AND AGREEMENT\n")
  cat("======================================================================\n")
  cat(
    "MTSPS threshold: ",
    format(mtsps_threshold, scientific = TRUE, digits = 8),
    "\n",
    sep = ""
  )
  cat("Individual GCM projections resolved: ", nrow(individual_rows), "\n", sep = "")
  cat("Scenarios processed: ", nrow(agreement_status_manifest), "\n", sep = "")
  cat(
    "Strict common valid domain (km²): ",
    format(strict_common_domain_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat(
    "Maximum high-confidence area-balance error (km²): ",
    format(max_balance_error, scientific = TRUE, digits = 8),
    "\n",
    sep = ""
  )
  cat(
    "Existing MESS/clamping rasters detected: ",
    sum(!is.na(extrapolation_inventory$file_path)),
    "\n",
    sep = ""
  )
  cat("\nKey tables:\n")
  cat(agreement_area_long_csv, "\n")
  cat(uncertainty_sd_summary_csv, "\n")
  cat(high_confidence_transition_csv, "\n")
  cat(scenario_confidence_summary_csv, "\n")
  cat("\nExtrapolation inventory:\n")
  cat(extrapolation_inventory_csv, "\n")
  cat("\nOutput folder:\n")
  cat(output_root, "\n")

  invisible(list(
    agreement_area = agreement_area_long,
    SD_summary = uncertainty_sd_summary,
    high_confidence_transition = high_conf_transition,
    scenario_confidence = scenario_confidence_summary
  ))
}

# =============================================================================
# F. EXECUTE WITH EXPLICIT FAILURE REPORTING
# =============================================================================

tryCatch(
  {
    run_inter_gcm_uncertainty_analysis()
  },
  error = function(error) {

    dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    failure_message <- paste0(
      "SCRIPT 18 FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(error), "\n\n",
      "Execution log: ", execution_log, "\n"
    )

    writeLines(failure_message, failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 18 FAILED\n")
    cat("======================================================================\n")
    cat(failure_message, "\n")

    stop(error)
  }
)
