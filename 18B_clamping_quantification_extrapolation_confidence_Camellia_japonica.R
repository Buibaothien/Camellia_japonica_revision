# =============================================================================
# 18B_CLAMPING_QUANTIFICATION_EXTRAPOLATION_CONFIDENCE_CAMELLIA_JAPONICA.R
# =============================================================================
# PURPOSE
# -------
# Quantify MaxEnt clamping across the 48 completed individual future GCM
# projections and create extrapolation-aware high-confidence suitability masks.
#
# FINAL MODEL CONTEXT
# -------------------
# Species: Camellia japonica
# Primary binary threshold: MTSPS
# GCMs: ACCESS-CM2, BCC-CSM2-MR, CMCC-ESM2
# Scenarios: 4 SSPs × 4 time periods = 16 scenarios
#
# CLAMPING INTERPRETATION
# -----------------------
# MaxEnt clamping layers are interpreted as follows:
#   - 0       = no predictor clamping
#   - non-zero = one or more predictor values were clamped
#
# This script does NOT treat clamping as model failure. Instead, it identifies
# where future suitability predictions rely on environmental extrapolation
# relative to calibration conditions and separates:
#
# A. Strict high-confidence future suitability:
#       3/3 GCM agreement AND no clamping in any GCM (clamping count = 0)
#
# B. Relaxed high-confidence future suitability:
#       3/3 GCM agreement AND clamping in no more than one GCM
#       (clamping count <= 1)
#
# The strict layer is intended for Script 19 climate-stable refugia and
# highest-confidence in situ conservation priorities.
#
# INPUTS REUSED
# -------------
# - Script 18: agreement maps, strict common valid domain, and resolved
#   individual GCM cloglog raster paths.
# - Script 11: MaxEnt clamping rasters.
#
# OUTPUTS
# -------
# 1. Individual GCM clamping extent by scenario
# 2. Number of GCMs with clamping (0-3) by scenario
# 3. Strict / relaxed high-confidence suitable masks
# 4. Strict-clamping-aware transition maps relative to current MTSPS suitability
# 5. Tables required to define a defensible Script 19 priority rule
# =============================================================================

options(warn = 1)

cat("\n")
cat("======================================================================\n")
cat("STARTING SCRIPT 18B — CLAMPING QUANTIFICATION AND CONFIDENCE MASKS\n")
cat("Species: Camellia japonica\n")
cat("Primary threshold: MTSPS\n")
cat("======================================================================\n")

script_started_at <- Sys.time()

# =============================================================================
# A. USER SETTINGS
# =============================================================================

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# The primary threshold is intentionally inherited from the final MTSPS analysis.
threshold_id <- "MTSPS"

# Clamping is normally exactly zero when no environmental predictor is clamped.
# A tiny tolerance avoids treating numerical round-off as clamping.
clamping_nonzero_tolerance <- 1e-12

# Expected GCM ordering in every scenario.
expected_gcms <- c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2")

# FALSE: reuse existing valid maps when present.
overwrite_existing_maps <- FALSE

# Create focal uncertainty figures for these two end-century scenarios.
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

script11_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "11_future_projection_individual_GCM_v2_plotsfix"
)

script18_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "18_interGCM_uncertainty_agreement_extrapolation"
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

resolved_individual_projection_csv <- file.path(
  script18_root,
  "01_setup",
  "02_resolved_individual_projection_rasters.csv"
)

script18_agreement_manifest_csv <- file.path(
  script18_root,
  "07_supporting_files",
  "agreement_map_status_manifest.csv"
)

script18_clamping_inventory_csv <- file.path(
  script18_root,
  "01_setup",
  "04_existing_MESS_clamping_inventory.csv"
)

strict_common_domain_tif <- file.path(
  script18_root,
  "07_supporting_files",
  "common_valid_domain_current_plus_all_48_individual_GCM_projections.tif"
)

script18_cell_area_tif <- file.path(
  script18_root,
  "07_supporting_files",
  "global_geodesic_cell_area_km2.tif"
)

# =============================================================================
# C. OUTPUT PATHS
# =============================================================================

output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "18B_clamping_quantification_extrapolation_confidence"
)

status_dir <- file.path(output_root, "00_run_status")
setup_dir <- file.path(output_root, "01_setup")
clamp_count_dir <- file.path(output_root, "02_clamping_GCM_count_maps")
strict_mask_dir <- file.path(output_root, "03_strict_unclamped_high_confidence_maps")
relaxed_mask_dir <- file.path(output_root, "04_relaxed_high_confidence_maps")
transition_dir <- file.path(output_root, "05_strict_unclamped_transitions")
area_dir <- file.path(output_root, "06_area_tables")
figure_dir <- file.path(output_root, "07_figures")
support_dir <- file.path(output_root, "08_supporting_files")
log_dir <- file.path(output_root, "09_logs")

started_marker <- file.path(status_dir, "18B_STARTED.txt")
completed_marker <- file.path(status_dir, "18B_COMPLETED.txt")
failed_marker <- file.path(status_dir, "18B_FAILED.txt")
execution_log <- file.path(log_dir, "18B_execution_log.txt")

input_manifest_csv <- file.path(setup_dir, "01_input_manifest.csv")
resolved_clamping_manifest_csv <- file.path(
  setup_dir,
  "02_resolved_clamping_raster_manifest.csv"
)
clamping_value_diagnostics_csv <- file.path(
  setup_dir,
  "03_clamping_value_and_coverage_diagnostics.csv"
)

clamping_count_legend_csv <- file.path(
  support_dir,
  "clamping_GCM_count_legend.csv"
)
strict_transition_legend_csv <- file.path(
  support_dir,
  "strict_unclamped_high_confidence_transition_legend.csv"
)
map_status_manifest_csv <- file.path(
  support_dir,
  "clamping_confidence_map_status_manifest.csv"
)
run_summary_csv <- file.path(
  support_dir,
  "clamping_analysis_run_summary.csv"
)

individual_clamping_area_csv <- file.path(
  area_dir,
  "01_individual_GCM_clamping_area_by_scenario_long.csv"
)
scenario_clamping_agreement_csv <- file.path(
  area_dir,
  "02_clamping_GCM_count_and_high_confidence_area_by_scenario.csv"
)
strict_transition_area_csv <- file.path(
  area_dir,
  "03_strict_unclamped_high_confidence_transition_area_by_scenario.csv"
)
area_balance_checks_csv <- file.path(
  area_dir,
  "04_strict_unclamped_area_balance_checks.csv"
)
priority_readiness_summary_csv <- file.path(
  area_dir,
  "05_clamping_priority_readiness_summary.csv"
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

read_single_numeric_value <- function(path, candidates, label) {
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  column <- find_column(data, candidates, label)

  values <- safe_number(data[[column]])
  values <- values[is.finite(values)]

  if (length(values) != 1L) {
    stop(
      "\nExpected exactly one finite numeric value in ", label, ".\n",
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
      "All relevant raster inputs must be aligned.",
      call. = FALSE
    )
  }
}

split_path_components <- function(path) {
  unlist(strsplit(
    gsub("\\\\", "/", normalizePath(path, winslash = "/", mustWork = FALSE)),
    "/",
    fixed = TRUE
  ))
}

common_prefix_depth <- function(path_a, path_b) {
  a <- split_path_components(path_a)
  b <- split_path_components(path_b)
  n <- min(length(a), length(b))
  depth <- 0L

  if (n == 0L) return(depth)

  for (i in seq_len(n)) {
    if (tolower(a[i]) == tolower(b[i])) {
      depth <- depth + 1L
    } else {
      break
    }
  }

  depth
}

extension_rank <- function(path) {
  extension <- tolower(tools::file_ext(path))

  if (extension %in% c("tif", "tiff")) return(3L)
  if (extension == "asc") return(2L)
  1L
}

calculate_area_by_code <- function(class_raster, cell_area_raster, domain_raster,
                                   codes, class_names) {
  assert_same_grid(cell_area_raster, class_raster, "cell-area raster")
  assert_same_grid(domain_raster, class_raster, "common valid-domain raster")

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
    zone_areas <- as.numeric(zonal_output[[ncol(zonal_output)]])

    for (i in seq_along(codes)) {
      matched <- zone_areas[zone_codes == codes[i]]
      if (length(matched) > 0) {
        output$area_km2[i] <- sum(matched, na.rm = TRUE)
      }
    }
  }

  output$area_10_4_km2 <- output$area_km2 / 10000
  output
}

get_code_area <- function(area_table, code) {
  result <- area_table$area_km2[area_table$code == code]

  if (length(result) != 1L) {
    stop("Could not retrieve a unique area for code ", code, call. = FALSE)
  }

  result
}

count_non_na_cells <- function(raster_object) {
  result <- terra::global(
    terra::ifel(is.na(raster_object), NA, 1),
    fun = "sum",
    na.rm = TRUE
  )

  as.numeric(result[1, 1])
}

make_binary_clamping_raster <- function(clamping_raster, domain_raster,
                                        tolerance) {
  # 0 = no clamping; 1 = non-zero clamping.
  output <- terra::ifel(
    is.na(domain_raster),
    NA,
    terra::ifel(abs(clamping_raster) > tolerance, 1, 0)
  )

  names(output) <- "clamping_present"
  output
}

resolve_clamping_raster <- function(
  projection_path,
  gcm,
  ssp,
  period_start,
  period_end,
  clamping_inventory_paths,
  reference_raster
) {
  valid_paths <- clamping_inventory_paths[
    !is.na(clamping_inventory_paths) &
      nzchar(clamping_inventory_paths) &
      file.exists(clamping_inventory_paths)
  ]

  valid_paths <- unique(valid_paths[
    grepl("\\.(tif|tiff|asc)$", valid_paths, ignore.case = TRUE)
  ])

  valid_paths <- valid_paths[
    grepl("clamp", basename(valid_paths), ignore.case = TRUE)
  ]

  if (length(valid_paths) == 0L) {
    stop("\nNo readable clamping rasters found in inventory.", call. = FALSE)
  }

  # A score uses the matching projection directory as the strongest signal.
  # Scenario strings are additional safeguards where those tokens occur in paths.
  projection_norm <- normalise_name(projection_path)
  gcm_norm <- normalise_name(gcm)
  ssp_norm <- normalise_name(ssp)
  start_text <- as.character(period_start)
  end_text <- as.character(period_end)

  candidate_scores <- vapply(valid_paths, function(candidate) {
    candidate_norm <- normalise_name(candidate)

    prefix_score <- 1000 * common_prefix_depth(
      dirname(projection_path),
      dirname(candidate)
    )

    token_score <- 0
    if (grepl(gcm_norm, candidate_norm, fixed = TRUE)) token_score <- token_score + 250
    if (grepl(ssp_norm, candidate_norm, fixed = TRUE)) token_score <- token_score + 150
    if (grepl(start_text, candidate, fixed = TRUE)) token_score <- token_score + 75
    if (grepl(end_text, candidate, fixed = TRUE)) token_score <- token_score + 75

    # Minor preference for TIFF when ASC and TIFF represent duplicate exports.
    prefix_score + token_score + extension_rank(candidate)
  }, numeric(1))

  candidate_table <- data.frame(
    path = valid_paths,
    score = candidate_scores,
    stringsAsFactors = FALSE
  )

  candidate_table <- candidate_table[
    order(candidate_table$score, decreasing = TRUE),
    ,
    drop = FALSE
  ]

  # Try candidate rasters in score order and retain the first one matching
  # the individual cloglog raster geometry.
  matching_paths <- character(0)

  for (candidate in candidate_table$path) {
    candidate_raster <- try(terra::rast(candidate), silent = TRUE)

    if (inherits(candidate_raster, "try-error")) next

    same_grid <- terra::compareGeom(
      reference_raster,
      candidate_raster,
      crs = TRUE,
      ext = TRUE,
      rowcol = TRUE,
      res = TRUE,
      stopOnError = FALSE
    )

    if (same_grid) matching_paths <- c(matching_paths, candidate)
  }

  matching_paths <- unique(matching_paths)

  if (length(matching_paths) == 0L) {
    stop(
      "\nCould not resolve a grid-matched clamping raster for:\n",
      gcm, " | ", ssp, " ", period_start, "-", period_end,
      "\nProjection: ", projection_path,
      call. = FALSE
    )
  }

  # Select highest-scored matching candidate; tif/tiff wins when duplicated.
  matching_scores <- candidate_scores[
    match(matching_paths, valid_paths)
  ]

  matching_order <- order(
    matching_scores,
    vapply(matching_paths, extension_rank, integer(1)),
    decreasing = TRUE
  )

  selected <- matching_paths[matching_order[1]]

  list(
    selected_path = normalizePath(selected, winslash = "/", mustWork = TRUE),
    matching_grid_candidate_count = length(matching_paths),
    candidates_considered = nrow(candidate_table)
  )
}

get_global_stat <- function(raster_object, statistic_name) {
  result <- terra::global(raster_object, fun = statistic_name, na.rm = TRUE)
  as.numeric(result[1, 1])
}

create_focal_clamping_figure <- function(
  clamping_count,
  strict_mask,
  strict_transition,
  scenario_label,
  output_pdf
) {
  grDevices::pdf(output_pdf, width = 14, height = 5.5, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(1, 3), mar = c(2.7, 2.8, 3.0, 1.0))

  plot(
    clamping_count,
    col = c("#F0F0F0", "#FEE08B", "#FDAE61", "#D73027"),
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
    plg = FALSE,
    main = "A. GCM clamping count",
    axes = TRUE
  )
  legend(
    "bottomleft",
    legend = c("0", "1", "2", "3"),
    fill = c("#F0F0F0", "#FEE08B", "#FDAE61", "#D73027"),
    title = "GCMs with clamping",
    bty = "n",
    cex = 0.78,
    inset = 0.02
  )

  plot(
    strict_mask,
    col = c("#F0F0F0", "#2166AC"),
    breaks = c(-0.5, 0.5, 1.5),
    plg = FALSE,
    main = "B. Strict high-confidence habitat",
    axes = TRUE
  )
  legend(
    "bottomleft",
    legend = c("No", "Yes"),
    fill = c("#F0F0F0", "#2166AC"),
    title = "3/3 GCM suitable\nand no clamping",
    bty = "n",
    cex = 0.78,
    inset = 0.02
  )

  plot(
    strict_transition,
    col = c("#F0F0F0", "#92C5DE", "#F4A582", "#0571B0"),
    breaks = c(0.5, 1.5, 2.5, 3.5, 4.5),
    plg = FALSE,
    main = "C. Strict confidence transition",
    axes = TRUE
  )
  legend(
    "bottomleft",
    legend = c(
      "Stable not strict",
      "Strict gain",
      "Loss of strict suitability",
      "Stable strict"
    ),
    fill = c("#F0F0F0", "#92C5DE", "#F4A582", "#0571B0"),
    bty = "n",
    cex = 0.70,
    inset = 0.02
  )

  mtext(
    scenario_label,
    side = 3,
    outer = FALSE,
    line = 0.6,
    font = 2,
    cex = 1.15
  )
}

# =============================================================================
# E. MAIN ANALYSIS
# =============================================================================

run_clamping_analysis <- function() {

  for (folder in c(
    output_root, status_dir, setup_dir, clamp_count_dir, strict_mask_dir,
    relaxed_mask_dir, transition_dir, area_dir, figure_dir, support_dir,
    log_dir
  )) {
    ensure_dir(folder)
  }

  unlink(completed_marker, force = TRUE)
  unlink(failed_marker, force = TRUE)

  writeLines(
    c(
      "Script 18B started",
      paste("Started:", format(script_started_at, "%Y-%m-%d %H:%M:%S")),
      "Purpose: clamping quantification and extrapolation-aware confidence masks"
    ),
    started_marker
  )

  if (!requireNamespace("terra", quietly = TRUE)) {
    install.packages("terra", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(terra))

  write_log("Package terra loaded.")

  # -------------------------------------------------------------------
  # E1. Read all required inputs
  # -------------------------------------------------------------------

  current_continuous_path <- find_first_existing_file(
    current_continuous_candidates,
    "current full-data continuous suitability raster"
  )

  required_files <- c(
    fixed_mtsps_csv,
    resolved_individual_projection_csv,
    script18_agreement_manifest_csv,
    script18_clamping_inventory_csv,
    strict_common_domain_tif,
    script18_cell_area_tif
  )

  missing_required <- required_files[!file.exists(required_files)]

  if (length(missing_required) > 0L) {
    stop(
      "\nRequired input file(s) missing:\n",
      paste(missing_required, collapse = "\n"),
      call. = FALSE
    )
  }

  mtsps_threshold <- read_single_numeric_value(
    fixed_mtsps_csv,
    c("MTSPS_threshold", "threshold"),
    "fixed MTSPS threshold CSV"
  )

  current_continuous <- terra::rast(current_continuous_path)
  strict_domain <- terra::rast(strict_common_domain_tif)
  cell_area <- terra::rast(script18_cell_area_tif)

  assert_same_grid(current_continuous, strict_domain, "Script 18 strict common domain")
  assert_same_grid(current_continuous, cell_area, "Script 18 cell-area raster")

  current_binary <- terra::ifel(
    is.na(strict_domain),
    NA,
    terra::ifel(current_continuous >= mtsps_threshold, 1, 0)
  )
  names(current_binary) <- "current_MTSPS_suitable"

  strict_domain_cells <- count_non_na_cells(strict_domain)

  strict_domain_area <- calculate_area_by_code(
    class_raster = strict_domain,
    cell_area_raster = cell_area,
    domain_raster = strict_domain,
    codes = 1,
    class_names = "Strict common valid domain"
  )

  strict_domain_km2 <- get_code_area(strict_domain_area, 1)

  resolved_projection <- read.csv(
    resolved_individual_projection_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  agreement_manifest <- read.csv(
    script18_agreement_manifest_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  clamping_inventory <- read.csv(
    script18_clamping_inventory_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(
      input_type = c(
        "Current final full-data continuous suitability",
        "Fixed MTSPS threshold",
        "Script 18 strict common domain",
        "Script 18 geodesic cell-area raster",
        "Script 18 resolved individual projections",
        "Script 18 agreement map manifest",
        "Script 18 clamping inventory"
      ),
      path = c(
        current_continuous_path,
        normalizePath(fixed_mtsps_csv, winslash = "/", mustWork = TRUE),
        normalizePath(strict_common_domain_tif, winslash = "/", mustWork = TRUE),
        normalizePath(script18_cell_area_tif, winslash = "/", mustWork = TRUE),
        normalizePath(resolved_individual_projection_csv, winslash = "/", mustWork = TRUE),
        normalizePath(script18_agreement_manifest_csv, winslash = "/", mustWork = TRUE),
        normalizePath(script18_clamping_inventory_csv, winslash = "/", mustWork = TRUE)
      ),
      stringsAsFactors = FALSE
    ),
    input_manifest_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E2. Identify table fields and validate 48 individual projections
  # -------------------------------------------------------------------

  projection_gcm_col <- find_column(
    resolved_projection,
    c("gcm_analysis", "gcm", "gcm_name"),
    "resolved individual projection table"
  )
  projection_ssp_col <- find_column(
    resolved_projection,
    c("ssp_analysis", "ssp"),
    "resolved individual projection table"
  )
  projection_start_col <- find_column(
    resolved_projection,
    c("period_start_analysis", "period_start", "periodstart"),
    "resolved individual projection table"
  )
  projection_end_col <- find_column(
    resolved_projection,
    c("period_end_analysis", "period_end", "periodend"),
    "resolved individual projection table"
  )
  projection_path_col <- find_column(
    resolved_projection,
    c("cloglog_path_analysis", "cloglog_path", "projection_path"),
    "resolved individual projection table"
  )

  resolved_projection$gcm_analysis <- as.character(
    resolved_projection[[projection_gcm_col]]
  )
  resolved_projection$ssp_analysis <- as.character(
    resolved_projection[[projection_ssp_col]]
  )
  resolved_projection$period_start_analysis <- safe_number(
    resolved_projection[[projection_start_col]]
  )
  resolved_projection$period_end_analysis <- safe_number(
    resolved_projection[[projection_end_col]]
  )
  resolved_projection$projection_path_analysis <- as.character(
    resolved_projection[[projection_path_col]]
  )

  resolved_projection <- resolved_projection[
    resolved_projection$gcm_analysis %in% expected_gcms &
      !is.na(resolved_projection$ssp_analysis) &
      is.finite(resolved_projection$period_start_analysis) &
      is.finite(resolved_projection$period_end_analysis) &
      file.exists(resolved_projection$projection_path_analysis),
    ,
    drop = FALSE
  ]

  if (nrow(resolved_projection) != 48L) {
    stop(
      "\nExpected exactly 48 resolved individual projection rasters, but found ",
      nrow(resolved_projection), ".",
      call. = FALSE
    )
  }

  agreement_ssp_col <- find_column(
    agreement_manifest,
    c("ssp"),
    "Script 18 agreement map manifest"
  )
  agreement_start_col <- find_column(
    agreement_manifest,
    c("period_start", "periodstart"),
    "Script 18 agreement map manifest"
  )
  agreement_end_col <- find_column(
    agreement_manifest,
    c("period_end", "periodend"),
    "Script 18 agreement map manifest"
  )
  agreement_path_col <- find_column(
    agreement_manifest,
    c("agreement_tif"),
    "Script 18 agreement map manifest"
  )

  agreement_manifest$ssp_analysis <- as.character(
    agreement_manifest[[agreement_ssp_col]]
  )
  agreement_manifest$period_start_analysis <- safe_number(
    agreement_manifest[[agreement_start_col]]
  )
  agreement_manifest$period_end_analysis <- safe_number(
    agreement_manifest[[agreement_end_col]]
  )
  agreement_manifest$agreement_path_analysis <- as.character(
    agreement_manifest[[agreement_path_col]]
  )

  agreement_manifest <- agreement_manifest[
    file.exists(agreement_manifest$agreement_path_analysis),
    ,
    drop = FALSE
  ]

  if (nrow(agreement_manifest) != 16L) {
    stop(
      "\nExpected 16 agreement maps from Script 18, but found ",
      nrow(agreement_manifest), ".",
      call. = FALSE
    )
  }

  inventory_path_col <- find_column(
    clamping_inventory,
    c("file_path", "path"),
    "Script 18 clamping inventory"
  )

  clamping_paths <- as.character(clamping_inventory[[inventory_path_col]])
  clamping_paths <- clamping_paths[
    !is.na(clamping_paths) &
      nzchar(clamping_paths) &
      file.exists(clamping_paths)
  ]

  clamping_paths <- unique(clamping_paths[
    grepl("\\.(tif|tiff|asc)$", clamping_paths, ignore.case = TRUE) &
      grepl("clamp", basename(clamping_paths), ignore.case = TRUE)
  ])

  if (length(clamping_paths) < 48L) {
    stop(
      "\nExpected at least 48 readable clamping files in the inventory, but found ",
      length(clamping_paths), ".",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------
  # E3. Resolve exactly one grid-matched clamping raster for each GCM map
  # -------------------------------------------------------------------

  write_log("Resolving one grid-matched clamping raster for each of 48 projections.")

  clamping_resolution_records <- list()

  for (i in seq_len(nrow(resolved_projection))) {
    reference_projection <- terra::rast(
      resolved_projection$projection_path_analysis[i]
    )

    assert_same_grid(
      current_continuous,
      reference_projection,
      paste0(
        resolved_projection$gcm_analysis[i],
        " | ",
        resolved_projection$ssp_analysis[i],
        " ",
        resolved_projection$period_start_analysis[i], "-",
        resolved_projection$period_end_analysis[i]
      )
    )

    resolution <- resolve_clamping_raster(
      projection_path = resolved_projection$projection_path_analysis[i],
      gcm = resolved_projection$gcm_analysis[i],
      ssp = resolved_projection$ssp_analysis[i],
      period_start = resolved_projection$period_start_analysis[i],
      period_end = resolved_projection$period_end_analysis[i],
      clamping_inventory_paths = clamping_paths,
      reference_raster = reference_projection
    )

    clamping_resolution_records[[length(clamping_resolution_records) + 1]] <- data.frame(
      gcm = resolved_projection$gcm_analysis[i],
      ssp = resolved_projection$ssp_analysis[i],
      period_start = resolved_projection$period_start_analysis[i],
      period_end = resolved_projection$period_end_analysis[i],
      individual_cloglog_tif = normalizePath(
        resolved_projection$projection_path_analysis[i],
        winslash = "/",
        mustWork = TRUE
      ),
      selected_clamping_raster = resolution$selected_path,
      selected_clamping_extension = tolower(
        tools::file_ext(resolution$selected_path)
      ),
      matching_grid_candidate_count = resolution$matching_grid_candidate_count,
      inventory_candidates_considered = resolution$candidates_considered,
      stringsAsFactors = FALSE
    )
  }

  resolved_clamping <- do.call(rbind, clamping_resolution_records)

  duplicate_selected_paths <- duplicated(resolved_clamping$selected_clamping_raster)

  if (any(duplicate_selected_paths)) {
    print(resolved_clamping[duplicate_selected_paths, , drop = FALSE])
    stop(
      "\nThe same clamping raster was assigned to more than one individual projection. ",
      "Automatic matching is ambiguous; do not proceed.",
      call. = FALSE
    )
  }

  resolved_clamping <- resolved_clamping[
    order(
      resolved_clamping$ssp,
      resolved_clamping$period_start,
      resolved_clamping$gcm
    ),
    ,
    drop = FALSE
  ]

  write.csv(
    resolved_clamping,
    resolved_clamping_manifest_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E4. Legends
  # -------------------------------------------------------------------

  write.csv(
    data.frame(
      clamping_GCM_count = 0:3,
      clamping_class = c(
        "No clamping",
        "Clamping in one GCM",
        "Clamping in two GCMs",
        "Clamping in all three GCMs"
      ),
      interpretation = c(
        "No individual GCM projection requires MaxEnt clamping at this cell.",
        "One of three individual GCM projections requires clamping.",
        "Two of three individual GCM projections require clamping.",
        "All three individual GCM projections require clamping."
      ),
      stringsAsFactors = FALSE
    ),
    clamping_count_legend_csv,
    row.names = FALSE
  )

  write.csv(
    data.frame(
      transition_code = 1:4,
      transition_class = c(
        "Stable not strict-unclamped high-confidence suitable",
        "Strict-unclamped high-confidence gain",
        "Loss of strict-unclamped high-confidence suitability",
        "Stable strict-unclamped high-confidence suitable"
      ),
      definition = c(
        "Current MTSPS unsuitable -> future not (3/3 GCM suitable and unclamped)",
        "Current MTSPS unsuitable -> future 3/3 GCM suitable and unclamped",
        "Current MTSPS suitable -> future not 3/3 GCM suitable and unclamped",
        "Current MTSPS suitable -> future 3/3 GCM suitable and unclamped"
      ),
      stringsAsFactors = FALSE
    ),
    strict_transition_legend_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E5. Process 16 scenarios
  # -------------------------------------------------------------------

  scenario_keys <- unique(
    resolved_clamping[
      c("ssp", "period_start", "period_end"),
      drop = FALSE
    ]
  )

  scenario_keys <- scenario_keys[
    order(scenario_keys$ssp, scenario_keys$period_start),
    ,
    drop = FALSE
  ]

  if (nrow(scenario_keys) != 16L) {
    stop(
      "\nExpected 16 unique SSP × period scenarios, but found ",
      nrow(scenario_keys), ".",
      call. = FALSE
    )
  }

  individual_area_records <- list()
  scenario_area_records <- list()
  transition_area_records <- list()
  diagnostic_records <- list()
  map_status_records <- list()

  for (i in seq_len(nrow(scenario_keys))) {
    ssp <- as.character(scenario_keys$ssp[i])
    period_start <- as.integer(scenario_keys$period_start[i])
    period_end <- as.integer(scenario_keys$period_end[i])
    scenario_label <- paste0(ssp, " ", period_start, "-", period_end)

    cat(
      "\n[Scenario ", i, "/", nrow(scenario_keys), "] ",
      scenario_label, "\n",
      sep = ""
    )

    scenario_rows <- resolved_clamping[
      resolved_clamping$ssp == ssp &
        resolved_clamping$period_start == period_start &
        resolved_clamping$period_end == period_end,
      ,
      drop = FALSE
    ]

    scenario_rows <- scenario_rows[
      match(expected_gcms, scenario_rows$gcm),
      ,
      drop = FALSE
    ]

    if (
      nrow(scenario_rows) != 3L ||
      any(is.na(scenario_rows$gcm))
    ) {
      stop(
        "\nCannot obtain exactly the expected 3 GCM rows for ",
        scenario_label,
        call. = FALSE
      )
    }

    agreement_row <- agreement_manifest[
      agreement_manifest$ssp_analysis == ssp &
        agreement_manifest$period_start_analysis == period_start &
        agreement_manifest$period_end_analysis == period_end,
      ,
      drop = FALSE
    ]

    if (nrow(agreement_row) != 1L) {
      stop(
        "\nCannot uniquely match Script 18 agreement raster for ",
        scenario_label,
        call. = FALSE
      )
    }

    agreement_raster <- terra::rast(agreement_row$agreement_path_analysis)
    assert_same_grid(
      current_continuous,
      agreement_raster,
      paste0("Script 18 agreement raster ", scenario_label)
    )

    clamp_count_path <- file.path(
      clamp_count_dir,
      paste0(
        species_name,
        "_clamping_GCM_count_",
        ssp, "_", period_start, "_", period_end,
        ".tif"
      )
    )

    strict_mask_path <- file.path(
      strict_mask_dir,
      paste0(
        species_name,
        "_strict_3of3_GCM_suitable_no_clamping_",
        ssp, "_", period_start, "_", period_end,
        "_MTSPS.tif"
      )
    )

    relaxed_mask_path <- file.path(
      relaxed_mask_dir,
      paste0(
        species_name,
        "_relaxed_3of3_GCM_suitable_clamping_at_most_1_",
        ssp, "_", period_start, "_", period_end,
        "_MTSPS.tif"
      )
    )

    strict_transition_path <- file.path(
      transition_dir,
      paste0(
        species_name,
        "_strict_unclamped_transition_current_to_",
        ssp, "_", period_start, "_", period_end,
        "_MTSPS.tif"
      )
    )

    all_maps_exist <- all(file.exists(c(
      clamp_count_path,
      strict_mask_path,
      relaxed_mask_path,
      strict_transition_path
    )))

    clamping_binary_list <- vector("list", length(expected_gcms))

    # Read each clamping raster, compute diagnostics, and calculate individual
    # clamped/unclamped areas for every GCM.
    for (j in seq_len(nrow(scenario_rows))) {
      gcm <- scenario_rows$gcm[j]
      clamping_raster <- terra::rast(
        scenario_rows$selected_clamping_raster[j]
      )

      assert_same_grid(
        current_continuous,
        clamping_raster,
        paste0("clamping raster ", gcm, " | ", scenario_label)
      )

      clamping_in_domain <- terra::mask(clamping_raster, strict_domain)

      covered_cells <- count_non_na_cells(clamping_in_domain)

      if (covered_cells != strict_domain_cells) {
        stop(
          "\nClamping raster does not cover the full strict common valid domain:\n",
          gcm, " | ", scenario_label,
          "\nCovered cells: ", covered_cells,
          "\nExpected cells: ", strict_domain_cells,
          call. = FALSE
        )
      }

      clamping_binary <- make_binary_clamping_raster(
        clamping_raster = clamping_raster,
        domain_raster = strict_domain,
        tolerance = clamping_nonzero_tolerance
      )

      clamping_binary_list[[j]] <- clamping_binary

      clamp_area <- calculate_area_by_code(
        class_raster = clamping_binary,
        cell_area_raster = cell_area,
        domain_raster = strict_domain,
        codes = c(0, 1),
        class_names = c("No clamping", "Clamping")
      )

      unclamped_area <- get_code_area(clamp_area, 0)
      clamped_area <- get_code_area(clamp_area, 1)

      individual_area_records[[length(individual_area_records) + 1]] <- data.frame(
        gcm = gcm,
        ssp = ssp,
        period_start = period_start,
        period_end = period_end,
        scenario = scenario_label,
        clamping_nonzero_tolerance = clamping_nonzero_tolerance,
        strict_common_valid_domain_km2 = strict_domain_km2,
        unclamped_km2 = unclamped_area,
        clamped_km2 = clamped_area,
        clamped_10_4_km2 = clamped_area / 10000,
        percent_common_domain_clamped = 100 * clamped_area / strict_domain_km2,
        clamping_raster = scenario_rows$selected_clamping_raster[j],
        stringsAsFactors = FALSE
      )

      diagnostic_records[[length(diagnostic_records) + 1]] <- data.frame(
        gcm = gcm,
        ssp = ssp,
        period_start = period_start,
        period_end = period_end,
        scenario = scenario_label,
        clamping_raster = scenario_rows$selected_clamping_raster[j],
        clamping_min = get_global_stat(clamping_in_domain, "min"),
        clamping_max = get_global_stat(clamping_in_domain, "max"),
        clamping_mean = get_global_stat(clamping_in_domain, "mean"),
        strict_domain_cells = strict_domain_cells,
        covered_domain_cells = covered_cells,
        matching_grid_candidate_count =
          scenario_rows$matching_grid_candidate_count[j],
        stringsAsFactors = FALSE
      )
    }

    if (!all_maps_exist || overwrite_existing_maps) {
      write_log("Building clamping-confidence maps: ", scenario_label)

      clamping_count <- clamping_binary_list[[1]] +
        clamping_binary_list[[2]] +
        clamping_binary_list[[3]]

      clamping_count <- terra::mask(clamping_count, strict_domain)
      names(clamping_count) <- "clamping_GCM_count"

      strict_mask <- terra::ifel(
        is.na(strict_domain),
        NA,
        terra::ifel(agreement_raster == 3 & clamping_count == 0, 1, 0)
      )
      names(strict_mask) <- "strict_3of3_GCM_suitable_no_clamping"

      relaxed_mask <- terra::ifel(
        is.na(strict_domain),
        NA,
        terra::ifel(agreement_raster == 3 & clamping_count <= 1, 1, 0)
      )
      names(relaxed_mask) <- "relaxed_3of3_GCM_suitable_clamping_at_most_1"

      strict_transition <- 1 + strict_mask + 2 * current_binary
      strict_transition <- terra::mask(strict_transition, strict_domain)
      names(strict_transition) <- "strict_unclamped_high_confidence_transition"

      terra::writeRaster(
        clamping_count,
        clamp_count_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      terra::writeRaster(
        strict_mask,
        strict_mask_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      terra::writeRaster(
        relaxed_mask,
        relaxed_mask_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      terra::writeRaster(
        strict_transition,
        strict_transition_path,
        overwrite = TRUE,
        datatype = "INT1U",
        wopt = list(gdal = "COMPRESS=LZW")
      )

      map_status <- "completed"
    } else {
      write_log("Reusing existing clamping-confidence maps: ", scenario_label)
      map_status <- "skipped_existing_output"

      clamping_count <- terra::rast(clamp_count_path)
      strict_mask <- terra::rast(strict_mask_path)
      relaxed_mask <- terra::rast(relaxed_mask_path)
      strict_transition <- terra::rast(strict_transition_path)
    }

    # Area of each clamping-GCM-count class over the strict common domain.
    clamping_count_area <- calculate_area_by_code(
      class_raster = clamping_count,
      cell_area_raster = cell_area,
      domain_raster = strict_domain,
      codes = 0:3,
      class_names = c(
        "No clamping in any GCM",
        "Clamping in 1 GCM",
        "Clamping in 2 GCMs",
        "Clamping in all 3 GCMs"
      )
    )

    # Intersect the Script 18 high-confidence (3/3 agreement) class with
    # clamping-count classes to show how much high-agreement habitat is
    # extrapolation-free versus clamped.
    agreement_clamp_class <- terra::ifel(
      is.na(strict_domain),
      NA,
      terra::ifel(agreement_raster == 3, clamping_count + 1, 0)
    )
    names(agreement_clamp_class) <- "high_confidence_suitability_by_clamping"

    high_conf_clamp_area <- calculate_area_by_code(
      class_raster = agreement_clamp_class,
      cell_area_raster = cell_area,
      domain_raster = strict_domain,
      codes = 0:4,
      class_names = c(
        "Not high-confidence suitable (not 3/3 agreement)",
        "3/3 suitable, 0 GCM clamped",
        "3/3 suitable, 1 GCM clamped",
        "3/3 suitable, 2 GCMs clamped",
        "3/3 suitable, 3 GCMs clamped"
      )
    )

    strict_area <- get_code_area(high_conf_clamp_area, 1)
    high_conf_clamp1_area <- get_code_area(high_conf_clamp_area, 2)
    high_conf_clamp2_area <- get_code_area(high_conf_clamp_area, 3)
    high_conf_clamp3_area <- get_code_area(high_conf_clamp_area, 4)

    total_high_conf_area <- strict_area +
      high_conf_clamp1_area +
      high_conf_clamp2_area +
      high_conf_clamp3_area

    relaxed_area <- strict_area + high_conf_clamp1_area

    scenario_area_records[[length(scenario_area_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = period_start,
      period_end = period_end,
      scenario = scenario_label,
      MTSPS_threshold = mtsps_threshold,
      strict_common_valid_domain_km2 = strict_domain_km2,
      clamping_0of3_GCM_km2 = get_code_area(clamping_count_area, 0),
      clamping_1of3_GCM_km2 = get_code_area(clamping_count_area, 1),
      clamping_2of3_GCM_km2 = get_code_area(clamping_count_area, 2),
      clamping_3of3_GCM_km2 = get_code_area(clamping_count_area, 3),
      high_confidence_3of3_GCM_suitable_km2 = total_high_conf_area,
      high_confidence_suitable_unclamped_km2 = strict_area,
      high_confidence_suitable_clamped_in_1GCM_km2 = high_conf_clamp1_area,
      high_confidence_suitable_clamped_in_2GCMs_km2 = high_conf_clamp2_area,
      high_confidence_suitable_clamped_in_3GCMs_km2 = high_conf_clamp3_area,
      strict_high_confidence_unclamped_km2 = strict_area,
      relaxed_high_confidence_clamp_at_most_1_km2 = relaxed_area,
      strict_percent_of_3of3_agreement_suitable = if (total_high_conf_area > 0) {
        100 * strict_area / total_high_conf_area
      } else {
        NA_real_
      },
      relaxed_percent_of_3of3_agreement_suitable = if (total_high_conf_area > 0) {
        100 * relaxed_area / total_high_conf_area
      } else {
        NA_real_
      },
      clamping_GCM_count_tif = normalizePath(
        clamp_count_path,
        winslash = "/",
        mustWork = TRUE
      ),
      strict_unclamped_mask_tif = normalizePath(
        strict_mask_path,
        winslash = "/",
        mustWork = TRUE
      ),
      relaxed_mask_tif = normalizePath(
        relaxed_mask_path,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )

    # Strict extrapolation-aware transition areas.
    strict_transition_area <- calculate_area_by_code(
      class_raster = strict_transition,
      cell_area_raster = cell_area,
      domain_raster = strict_domain,
      codes = 1:4,
      class_names = c(
        "Stable not strict-unclamped high-confidence suitable",
        "Strict-unclamped high-confidence gain",
        "Loss of strict-unclamped high-confidence suitability",
        "Stable strict-unclamped high-confidence suitable"
      )
    )

    stable_not_strict <- get_code_area(strict_transition_area, 1)
    strict_gain <- get_code_area(strict_transition_area, 2)
    strict_loss <- get_code_area(strict_transition_area, 3)
    stable_strict <- get_code_area(strict_transition_area, 4)

    current_suitable <- strict_loss + stable_strict
    future_strict <- strict_gain + stable_strict

    transition_area_records[[length(transition_area_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = period_start,
      period_end = period_end,
      scenario = scenario_label,
      MTSPS_threshold = mtsps_threshold,
      strict_common_valid_domain_km2 = strict_domain_km2,
      stable_not_strict_unclamped_high_confidence_km2 = stable_not_strict,
      strict_unclamped_high_confidence_gain_km2 = strict_gain,
      loss_of_strict_unclamped_high_confidence_km2 = strict_loss,
      stable_strict_unclamped_high_confidence_km2 = stable_strict,
      current_MTSPS_suitable_km2 = current_suitable,
      future_strict_unclamped_high_confidence_km2 = future_strict,
      net_strict_unclamped_high_confidence_change_km2 =
        future_strict - current_suitable,
      transition_tif = normalizePath(
        strict_transition_path,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )

    map_status_records[[length(map_status_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = period_start,
      period_end = period_end,
      scenario = scenario_label,
      status = map_status,
      clamping_GCM_count_tif = normalizePath(
        clamp_count_path,
        winslash = "/",
        mustWork = TRUE
      ),
      strict_unclamped_mask_tif = normalizePath(
        strict_mask_path,
        winslash = "/",
        mustWork = TRUE
      ),
      relaxed_mask_tif = normalizePath(
        relaxed_mask_path,
        winslash = "/",
        mustWork = TRUE
      ),
      strict_transition_tif = normalizePath(
        strict_transition_path,
        winslash = "/",
        mustWork = TRUE
      ),
      stringsAsFactors = FALSE
    )

    # Create figures for pre-defined focal scenarios.
    is_focal <- any(
      focal_figure_scenarios$ssp == ssp &
        focal_figure_scenarios$period_start == period_start &
        focal_figure_scenarios$period_end == period_end
    )

    if (is_focal) {
      figure_path <- file.path(
        figure_dir,
        paste0(
          "Figure_S_clamping_confidence_",
          ssp, "_", period_start, "_", period_end,
          ".pdf"
        )
      )

      tryCatch(
        {
          create_focal_clamping_figure(
            clamping_count = clamping_count,
            strict_mask = strict_mask,
            strict_transition = strict_transition,
            scenario_label = scenario_label,
            output_pdf = figure_path
          )
          write_log("Created focal clamping figure: ", figure_path)
        },
        error = function(e) {
          write_log(
            "Focal clamping figure skipped for ", scenario_label,
            " because: ", conditionMessage(e)
          )
        }
      )
    }

    rm(
      clamping_binary_list,
      agreement_raster,
      clamping_count,
      strict_mask,
      relaxed_mask,
      strict_transition,
      agreement_clamp_class
    )
    gc()
  }

  # -------------------------------------------------------------------
  # E6. Final tables, quality checks, and completion markers
  # -------------------------------------------------------------------

  individual_clamping_area <- do.call(rbind, individual_area_records)
  scenario_clamping_area <- do.call(rbind, scenario_area_records)
  strict_transition_area <- do.call(rbind, transition_area_records)
  clamping_diagnostics <- do.call(rbind, diagnostic_records)
  map_status_manifest <- do.call(rbind, map_status_records)

  individual_clamping_area <- individual_clamping_area[
    order(
      individual_clamping_area$ssp,
      individual_clamping_area$period_start,
      individual_clamping_area$gcm
    ),
    ,
    drop = FALSE
  ]

  scenario_clamping_area <- scenario_clamping_area[
    order(
      scenario_clamping_area$ssp,
      scenario_clamping_area$period_start
    ),
    ,
    drop = FALSE
  ]

  strict_transition_area <- strict_transition_area[
    order(
      strict_transition_area$ssp,
      strict_transition_area$period_start
    ),
    ,
    drop = FALSE
  ]

  clamping_diagnostics <- clamping_diagnostics[
    order(
      clamping_diagnostics$ssp,
      clamping_diagnostics$period_start,
      clamping_diagnostics$gcm
    ),
    ,
    drop = FALSE
  ]

  map_status_manifest <- map_status_manifest[
    order(
      map_status_manifest$ssp,
      map_status_manifest$period_start
    ),
    ,
    drop = FALSE
  ]

  area_balance_checks <- data.frame(
    scenario = strict_transition_area$scenario,
    ssp = strict_transition_area$ssp,
    period_start = strict_transition_area$period_start,
    period_end = strict_transition_area$period_end,
    current_suitable_balance_error_km2 =
      strict_transition_area$current_MTSPS_suitable_km2 -
      (
        strict_transition_area$loss_of_strict_unclamped_high_confidence_km2 +
          strict_transition_area$stable_strict_unclamped_high_confidence_km2
      ),
    future_strict_balance_error_km2 =
      strict_transition_area$future_strict_unclamped_high_confidence_km2 -
      (
        strict_transition_area$strict_unclamped_high_confidence_gain_km2 +
          strict_transition_area$stable_strict_unclamped_high_confidence_km2
      ),
    stringsAsFactors = FALSE
  )

  max_balance_error <- max(
    abs(
      unlist(
        area_balance_checks[
          c(
            "current_suitable_balance_error_km2",
            "future_strict_balance_error_km2"
          )
        ]
      )
    ),
    na.rm = TRUE
  )

  priority_readiness_summary <- scenario_clamping_area[
    c(
      "scenario",
      "ssp",
      "period_start",
      "period_end",
      "MTSPS_threshold",
      "high_confidence_3of3_GCM_suitable_km2",
      "strict_high_confidence_unclamped_km2",
      "relaxed_high_confidence_clamp_at_most_1_km2",
      "strict_percent_of_3of3_agreement_suitable",
      "relaxed_percent_of_3of3_agreement_suitable"
    )
  ]

  priority_readiness_summary$priority_rule_strict <- paste0(
    "Future 3/3 GCM suitable AND clamping GCM count = 0; MTSPS = ",
    format(mtsps_threshold, scientific = TRUE, digits = 8)
  )

  priority_readiness_summary$priority_rule_relaxed <- paste0(
    "Future 3/3 GCM suitable AND clamping GCM count <= 1; MTSPS = ",
    format(mtsps_threshold, scientific = TRUE, digits = 8)
  )

  write.csv(
    individual_clamping_area,
    individual_clamping_area_csv,
    row.names = FALSE
  )
  write.csv(
    scenario_clamping_area,
    scenario_clamping_agreement_csv,
    row.names = FALSE
  )
  write.csv(
    strict_transition_area,
    strict_transition_area_csv,
    row.names = FALSE
  )
  write.csv(
    area_balance_checks,
    area_balance_checks_csv,
    row.names = FALSE
  )
  write.csv(
    priority_readiness_summary,
    priority_readiness_summary_csv,
    row.names = FALSE
  )
  write.csv(
    clamping_diagnostics,
    clamping_value_diagnostics_csv,
    row.names = FALSE
  )
  write.csv(
    map_status_manifest,
    map_status_manifest_csv,
    row.names = FALSE
  )

  run_summary <- data.frame(
    species = species_name,
    primary_threshold_id = threshold_id,
    MTSPS_threshold = mtsps_threshold,
    expected_GCMs = paste(expected_gcms, collapse = "; "),
    individual_projection_count = nrow(resolved_clamping),
    unique_clamping_rasters_resolved = length(
      unique(resolved_clamping$selected_clamping_raster)
    ),
    scenario_count = nrow(scenario_keys),
    strict_common_valid_domain_km2 = strict_domain_km2,
    clamping_nonzero_tolerance = clamping_nonzero_tolerance,
    maximum_strict_transition_area_balance_error_km2 = max_balance_error,
    stringsAsFactors = FALSE
  )

  write.csv(run_summary, run_summary_csv, row.names = FALSE)

  writeLines(
    c(
      "SCRIPT 18B COMPLETED SUCCESSFULLY",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("MTSPS threshold:", format(mtsps_threshold, scientific = TRUE, digits = 8)),
      paste("Individual clamping rasters resolved:", nrow(resolved_clamping)),
      paste("Unique clamping rasters:", length(unique(resolved_clamping$selected_clamping_raster))),
      paste("Scenarios processed:", nrow(scenario_keys)),
      paste("Strict common valid domain (km2):", format(strict_domain_km2, scientific = FALSE)),
      paste("Maximum strict transition balance error (km2):", format(max_balance_error, scientific = TRUE, digits = 8))
    ),
    completed_marker
  )

  write_log("Script completed successfully.")

  cat("\n======================================================================\n")
  cat("SCRIPT 18B COMPLETED — CLAMPING QUANTIFICATION AND CONFIDENCE MASKS\n")
  cat("======================================================================\n")
  cat(
    "MTSPS threshold: ",
    format(mtsps_threshold, scientific = TRUE, digits = 8),
    "\n",
    sep = ""
  )
  cat(
    "Individual clamping rasters resolved: ",
    nrow(resolved_clamping),
    "\n",
    sep = ""
  )
  cat(
    "Unique clamping rasters: ",
    length(unique(resolved_clamping$selected_clamping_raster)),
    "\n",
    sep = ""
  )
  cat("Scenarios processed: ", nrow(scenario_keys), "\n", sep = "")
  cat(
    "Strict common valid domain (km²): ",
    format(strict_domain_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat(
    "Maximum strict transition balance error (km²): ",
    format(max_balance_error, scientific = TRUE, digits = 8),
    "\n",
    sep = ""
  )
  cat("\nKey tables:\n")
  cat(individual_clamping_area_csv, "\n")
  cat(scenario_clamping_agreement_csv, "\n")
  cat(strict_transition_area_csv, "\n")
  cat(priority_readiness_summary_csv, "\n")
  cat("\nResolved clamping manifest:\n")
  cat(resolved_clamping_manifest_csv, "\n")
  cat("\nOutput folder:\n")
  cat(output_root, "\n")

  invisible(list(
    resolved_clamping = resolved_clamping,
    individual_clamping_area = individual_clamping_area,
    scenario_clamping_area = scenario_clamping_area,
    strict_transition_area = strict_transition_area,
    priority_readiness = priority_readiness_summary
  ))
}

# =============================================================================
# F. EXECUTE WITH EXPLICIT FAILURE REPORTING
# =============================================================================

tryCatch(
  {
    run_clamping_analysis()
  },
  error = function(error) {

    dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    failure_message <- paste0(
      "SCRIPT 18B FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(error), "\n\n",
      "Execution log: ", execution_log, "\n"
    )

    writeLines(failure_message, failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 18B FAILED\n")
    cat("======================================================================\n")
    cat(failure_message, "\n")

    stop(error)
  }
)
