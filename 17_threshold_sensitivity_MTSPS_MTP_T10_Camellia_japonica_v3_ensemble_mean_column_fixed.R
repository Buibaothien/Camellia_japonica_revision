# =============================================================================
# 17_THRESHOLD_SENSITIVITY_MTSPS_MTP_T10_CAMELLIA_JAPONICA_V3_ENSEMBLE_MEAN_COLUMN_FIXED.R
# =============================================================================
# PURPOSE
# -------
# Threshold-sensitivity analysis for the final calibrated Camellia japonica
# MaxEnt model. MTSPS remains the primary threshold. MTP and T10 are
# supplementary thresholds used to test whether conclusions about suitable
# area, expansion, contraction and stable suitable habitat are robust.
#
# Outputs are written ONLY to:
# E:/PREDICTION TREE/Camellia_japonica/global_reanalysis_outputs/
# 17_threshold_sensitivity_MTSPS_MTP_T10/
#
# This script does NOT rerun MaxEnt, individual future projections, GCM
# ensemble prediction, or the existing MTSPS analyses.
# =============================================================================

options(warn = 1)

cat("\n======================================================================\n")
cat("STARTING SCRIPT 17 — THRESHOLD SENSITIVITY ANALYSIS\n")
cat("Species: Camellia japonica\n")
cat("Primary threshold: MTSPS | Supplementary: MTP and T10\n")
cat("======================================================================\n")

# =============================================================================
# 1. USER SETTINGS
# =============================================================================

species_name <- "Camellia_japonica"
project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# Keep MTSPS as primary; MTP and T10 are supplementary only.
threshold_order <- c("MTSPS", "MTP", "T10")

# Number of valid current-prediction cells sampled only for the diagnostic
# threshold-context figure and table. It does not affect area calculations.
background_sample_size <- 50000L
random_seed <- 20260707L

# FALSE retains valid Script 17 rasters if this script is rerun.
overwrite_existing_maps <- FALSE

# Leave NA to locate the final full-data maxentResults.csv automatically.
final_maxent_results_csv <- NA_character_

# Exact continuous ensemble-mean suitability column in the Script 12 manifest.
# This must be the mean suitability raster, not ensemble_sd_tif or
# valid_GCM_count_tif.
future_ensemble_mean_raster_column <- "ensemble_mean_tif"

# =============================================================================
# 2. FIXED PATHS
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

fixed_mtsps_csv <- file.path(
  final_model_root,
  "04_CV_metrics_MTSPS_and_areas",
  "03_FINAL_single_MTSPS_threshold.csv"
)

final_occurrence_candidates <- c(
  file.path(
    project_dir,
    "global_reanalysis_outputs",
    "09_kuenm_global_OR10",
    "02_occurrences",
    "Camellia_japonica_joint.csv"
  ),
  file.path(
    project_dir,
    "global_reanalysis_outputs",
    "09_kuenm_global_OR10",
    "02_occurrences",
    "Camellia_japonica_final_occurrences.csv"
  )
)

current_final_map_candidates <- c(
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

# =============================================================================
# 3. OUTPUT PATHS
# =============================================================================

output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "17_threshold_sensitivity_MTSPS_MTP_T10"
)

status_dir <- file.path(output_root, "00_run_status")
setup_dir <- file.path(output_root, "01_setup")
threshold_dir <- file.path(output_root, "02_thresholds_and_current_scores")
binary_root <- file.path(output_root, "03_binary_maps")
transition_root <- file.path(output_root, "04_transition_maps")
area_dir <- file.path(output_root, "05_area_tables")
figure_dir <- file.path(output_root, "06_figures")
support_dir <- file.path(output_root, "07_supporting_files")
log_dir <- file.path(output_root, "08_logs")

started_marker <- file.path(status_dir, "17_STARTED.txt")
completed_marker <- file.path(status_dir, "17_COMPLETED.txt")
failed_marker <- file.path(status_dir, "17_FAILED.txt")
execution_log <- file.path(log_dir, "17_execution_log.txt")

maxent_candidates_csv <- file.path(
  setup_dir,
  "final_maxentResults_candidates.csv"
)

input_manifest_csv <- file.path(setup_dir, "input_manifest.csv")
thresholds_csv <- file.path(threshold_dir, "01_threshold_values.csv")
presence_summary_csv <- file.path(
  threshold_dir,
  "02_current_prediction_summary_at_occurrences.csv"
)
presence_scores_csv <- file.path(
  threshold_dir,
  "03_current_prediction_values_at_final_occurrences.csv"
)
background_scores_csv <- file.path(
  threshold_dir,
  "04_random_global_background_current_predictions.csv"
)
threshold_context_pdf <- file.path(
  figure_dir,
  "Figure_S_threshold_context_current_predictions.pdf"
)

binary_legend_csv <- file.path(support_dir, "binary_suitability_legend.csv")
transition_legend_csv <- file.path(support_dir, "transition_legend.csv")
cell_area_tif <- file.path(support_dir, "global_geodesic_cell_area_km2.tif")
binary_manifest_csv <- file.path(support_dir, "binary_map_status_manifest.csv")
transition_manifest_csv <- file.path(support_dir, "transition_map_status_manifest.csv")
run_summary_csv <- file.path(support_dir, "threshold_sensitivity_run_summary.csv")

current_area_csv <- file.path(area_dir, "01_current_suitable_area_by_threshold.csv")
future_area_long_csv <- file.path(area_dir, "02_future_suitable_area_by_threshold_long.csv")
future_area_wide_csv <- file.path(area_dir, "03_future_suitable_area_by_threshold_wide.csv")
change_area_long_csv <- file.path(area_dir, "04_distribution_change_by_threshold_long.csv")
change_area_wide_csv <- file.path(area_dir, "05_distribution_change_by_threshold_wide.csv")
balance_csv <- file.path(area_dir, "06_area_balance_checks.csv")
comparison_csv <- file.path(area_dir, "07_threshold_comparison_summary.csv")

# =============================================================================
# 4. FUNCTIONS
# =============================================================================

ensure_dir <- function(x) {
  dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

log_message <- function(...) {
  txt <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  out <- paste0("[", stamp, "] ", txt)
  cat(out, "\n")
  cat(out, "\n", file = execution_log, append = TRUE)
}

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

find_column <- function(dat, candidates, label) {
  nms <- names(dat)
  norm <- normalise_name(nms)
  cand <- normalise_name(candidates)
  hit <- match(cand, norm, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit)) return(nms[hit[1]])
  for (x in cand) {
    p <- grep(x, norm, fixed = TRUE)
    if (length(p)) return(nms[p[1]])
  }
  stop(
    "\nCannot identify required column in ", label, ".\nAvailable columns:\n",
    paste(nms, collapse = ", "),
    call. = FALSE
  )
}

first_existing <- function(paths, label) {
  found <- paths[file.exists(paths)]
  if (!length(found)) {
    stop(
      "\nNo ", label, " found. Checked:\n",
      paste(paths, collapse = "\n"),
      call. = FALSE
    )
  }
  normalizePath(found[1], winslash = "/", mustWork = TRUE)
}

same_grid <- function(reference, candidate, label) {
  ok <- terra::compareGeom(
    reference, candidate, crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
    stopOnError = FALSE
  )
  if (!ok) {
    stop(
      "\nRaster geometry mismatch: ", label,
      "\nAll current and future suitability rasters must be exactly aligned.",
      call. = FALSE
    )
  }
}

read_scalar <- function(path, candidates, label) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  col <- find_column(x, candidates, label)
  val <- as_num(x[[col]])
  val <- val[is.finite(val)]
  if (length(val) != 1) {
    stop("\nExpected exactly one finite value in ", label, ".", call. = FALSE)
  }
  val[1]
}

discover_final_maxent_results <- function(root, manual_path, candidate_file) {
  if (!is.na(manual_path) && nzchar(manual_path)) {
    if (!file.exists(manual_path)) {
      stop("\nManual maxentResults.csv path does not exist:\n", manual_path, call. = FALSE)
    }
    return(normalizePath(manual_path, winslash = "/", mustWork = TRUE))
  }

  files <- list.files(
    root, pattern = "^maxentResults\\.csv$", full.names = TRUE,
    recursive = TRUE, ignore.case = TRUE
  )

  if (!length(files)) {
    stop("\nNo maxentResults.csv found under:\n", root, call. = FALSE)
  }

  n <- normalise_name(files)
  score <- 4 * grepl("final", n) +
    3 * grepl("fulldata|full", n) +
    2 * grepl("model", n) -
    4 * grepl("cv|fold|replicate", n)

  tab <- data.frame(
    path = normalizePath(files, winslash = "/", mustWork = TRUE),
    score = score,
    stringsAsFactors = FALSE
  )
  write.csv(tab, candidate_file, row.names = FALSE)

  best <- files[score == max(score)]
  if (length(best) != 1) {
    stop(
      "\nMore than one candidate final maxentResults.csv was found.\n",
      "Inspect:\n", candidate_file,
      "\nThen set final_maxent_results_csv manually at the top of this script.",
      call. = FALSE
    )
  }
  normalizePath(best[1], winslash = "/", mustWork = TRUE)
}

read_threshold_from_maxent <- function(path, type) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  nms <- names(x)
  nn <- normalise_name(nms)

  needed <- switch(
    type,
    MTP = c("minimumtrainingpresence", "cloglogthreshold"),
    T10 = c("10percentiletrainingpresence", "cloglogthreshold"),
    stop("Unknown threshold type.", call. = FALSE)
  )

  hit <- rep(TRUE, length(nms))
  for (term in needed) hit <- hit & grepl(term, nn, fixed = TRUE)
  cols <- nms[hit]

  if (length(cols) != 1) {
    threshold_cols <- nms[grepl("threshold", nn, fixed = TRUE)]
    stop(
      "\nCould not uniquely identify ", type, " cloglog threshold.\n",
      "Threshold-related columns in maxentResults.csv:\n",
      paste(threshold_cols, collapse = "\n"),
      call. = FALSE
    )
  }

  value <- as_num(x[[cols]][1])
  if (!is.finite(value)) {
    stop("\n", type, " in maxentResults.csv is not finite.", call. = FALSE)
  }
  list(value = value, column = cols)
}

find_raster_path_column <- function(dat, preferred, label) {
  nms <- names(dat)
  nn <- normalise_name(nms)
  pref <- normalise_name(preferred)

  direct <- match(pref, nn, nomatch = 0)
  direct <- direct[direct > 0]

  for (i in direct) {
    vals <- as.character(dat[[nms[i]]])
    good <- !is.na(vals) & nzchar(vals) &
      file.exists(vals) & grepl("\\.(tif|tiff)$", vals, ignore.case = TRUE)
    if (sum(good) > 0) return(nms[i])
  }

  score <- vapply(nms, function(nm) {
    vals <- as.character(dat[[nm]])
    sum(!is.na(vals) & nzchar(vals) &
      file.exists(vals) & grepl("\\.(tif|tiff)$", vals, ignore.case = TRUE))
  }, numeric(1))

  if (max(score) == 0) {
    stop(
      "\nNo existing TIFF path found in ", label, ".\nColumns:\n",
      paste(nms, collapse = ", "),
      call. = FALSE
    )
  }

  top <- nms[score == max(score)]
  if (length(top) != 1) {
    stop(
      "\nMore than one possible raster path column in ", label, ":\n",
      paste(top, collapse = ", "),
      call. = FALSE
    )
  }
  top[1]
}

area_by_code <- function(class_raster, cell_area, codes, labels) {
  same_grid(cell_area, class_raster, "area calculation raster")
  z <- terra::zonal(cell_area, class_raster, fun = "sum", na.rm = TRUE)
  out <- data.frame(
    code = codes, class_name = labels, area_km2 = 0,
    stringsAsFactors = FALSE
  )
  if (!is.null(z) && nrow(z)) {
    zid <- suppressWarnings(as.integer(round(as.numeric(z[[1]]))))
    a <- as.numeric(z[[ncol(z)]])
    for (i in seq_along(codes)) {
      out$area_km2[i] <- sum(a[zid == codes[i]], na.rm = TRUE)
    }
  }
  out$area_10_4_km2 <- out$area_km2 / 10000
  out
}

make_binary <- function(suitability, threshold, outpath) {
  x <- terra::ifel(
    is.na(suitability), NA,
    terra::ifel(suitability >= threshold, 1, 0)
  )
  terra::writeRaster(
    x, outpath, overwrite = TRUE, datatype = "INT1U",
    wopt = list(gdal = "COMPRESS=LZW")
  )
}

make_transition <- function(current_binary, future_binary, outpath) {
  # 1 = stable unsuitable; 2 = gain; 3 = loss; 4 = stable suitable.
  x <- 1 + future_binary + 2 * current_binary
  terra::writeRaster(
    x, outpath, overwrite = TRUE, datatype = "INT1U",
    wopt = list(gdal = "COMPRESS=LZW")
  )
}

make_context_figure <- function(pres, bg, thresholds, outpath) {
  grDevices::pdf(outpath, width = 10, height = 6.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  hist(
    bg, breaks = 100, freq = FALSE, col = "grey90", border = "grey70",
    main = expression(
      "Current cloglog predictions: global background and final " *
        italic("C. japonica") * " occurrences"
    ),
    xlab = "Current cloglog suitability"
  )
  lines(stats::density(pres), lwd = 2)

  lty <- c(MTSPS = 1, MTP = 2, T10 = 3)
  for (i in seq_len(nrow(thresholds))) {
    abline(
      v = thresholds$threshold_cloglog[i],
      lty = lty[thresholds$threshold_id[i]],
      lwd = 2
    )
  }
  legend(
    "topright",
    legend = c(
      "Random valid global background cells",
      "Final occurrence predictions",
      paste0(
        thresholds$threshold_id, " = ",
        format(thresholds$threshold_cloglog, scientific = TRUE, digits = 5)
      )
    ),
    pch = c(15, NA, NA, NA, NA),
    lty = c(NA, 1, lty[thresholds$threshold_id]),
    lwd = c(NA, 2, rep(2, nrow(thresholds))),
    col = c("grey70", "black", rep("black", nrow(thresholds))),
    bty = "n"
  )
}

# =============================================================================
# 5. MAIN WORKFLOW
# =============================================================================

run_script17 <- function() {

  for (d in c(
    output_root, status_dir, setup_dir, threshold_dir, binary_root,
    transition_root, area_dir, figure_dir, support_dir, log_dir
  )) ensure_dir(d)

  unlink(completed_marker, force = TRUE)
  unlink(failed_marker, force = TRUE)

  writeLines(
    c(
      "SCRIPT 17 STARTED",
      paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
    ),
    started_marker
  )

  if (!requireNamespace("terra", quietly = TRUE)) {
    install.packages("terra", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(terra))

  log_message("Package terra loaded.")

  # -------------------------------------------------------------------
  # 5.1 Read all required inputs
  # -------------------------------------------------------------------

  if (!file.exists(fixed_mtsps_csv)) {
    stop("\nMissing MTSPS CSV:\n", fixed_mtsps_csv, call. = FALSE)
  }

  current_map <- first_existing(
    current_final_map_candidates,
    "final full-data current cloglog raster"
  )
  occurrences_csv <- first_existing(
    final_occurrence_candidates,
    "final 828 occurrence CSV"
  )
  ensemble_manifest_path <- first_existing(
    ensemble_manifest_candidates,
    "future ensemble manifest"
  )

  maxent_csv <- discover_final_maxent_results(
    final_model_root, final_maxent_results_csv, maxent_candidates_csv
  )

  mtsps <- read_scalar(
    fixed_mtsps_csv, c("MTSPS_threshold", "threshold"),
    "fixed MTSPS CSV"
  )
  mtp <- read_threshold_from_maxent(maxent_csv, "MTP")
  t10 <- read_threshold_from_maxent(maxent_csv, "T10")

  threshold_table <- data.frame(
    threshold_id = c("MTSPS", "MTP", "T10"),
    threshold_name = c(
      "Maximum training sensitivity plus specificity",
      "Minimum training presence",
      "10th percentile training presence"
    ),
    role = c("Primary", "Supplementary", "Supplementary"),
    threshold_cloglog = c(mtsps, mtp$value, t10$value),
    source_file = c(fixed_mtsps_csv, maxent_csv, maxent_csv),
    source_column = c("MTSPS_threshold", mtp$column, t10$column),
    stringsAsFactors = FALSE
  )

  threshold_table <- threshold_table[
    match(threshold_order, threshold_table$threshold_id), ,
    drop = FALSE
  ]

  if (any(!is.finite(threshold_table$threshold_cloglog)) ||
      any(threshold_table$threshold_cloglog < 0) ||
      any(threshold_table$threshold_cloglog > 1)) {
    stop(
      "\nAt least one threshold is invalid. Expected cloglog values in [0,1].",
      call. = FALSE
    )
  }

  write.csv(threshold_table, thresholds_csv, row.names = FALSE)

  log_message(
    "Thresholds: ",
    paste(
      threshold_table$threshold_id,
      format(threshold_table$threshold_cloglog, scientific = TRUE, digits = 8),
      sep = "=",
      collapse = "; "
    )
  )

  current <- terra::rast(current_map)

  em <- read.csv(
    ensemble_manifest_path, check.names = FALSE, stringsAsFactors = FALSE
  )
  status_col <- find_column(em, c("status"), "future ensemble manifest")
  ssp_col <- find_column(em, c("ssp"), "future ensemble manifest")
  ps_col <- find_column(
    em, c("period_start", "periodstart"), "future ensemble manifest"
  )
  pe_col <- find_column(
    em, c("period_end", "periodend"), "future ensemble manifest"
  )
  if (!future_ensemble_mean_raster_column %in% names(em)) {
    stop(
      "\nThe required ensemble-mean raster column was not found in the future ensemble manifest:\n",
      future_ensemble_mean_raster_column,
      "\nAvailable columns:\n",
      paste(names(em), collapse = ", "),
      "\n\nSet future_ensemble_mean_raster_column at the top of this script only if your manifest uses a different mean-suitability column.",
      call. = FALSE
    )
  }

  raster_col <- future_ensemble_mean_raster_column

  ok_status <- c("completed", "skipped_existing_output")
  future <- em[
    em[[status_col]] %in% ok_status &
      !is.na(em[[raster_col]]) & nzchar(em[[raster_col]]) &
      file.exists(em[[raster_col]]),
    ,
    drop = FALSE
  ]

  if (nrow(future) != 16L) {
    stop(
      "\nExpected 16 completed ensemble rasters, found ", nrow(future),
      ". Inspect:\n", ensemble_manifest_path,
      call. = FALSE
    )
  }

  future$path <- normalizePath(
    as.character(future[[raster_col]]), winslash = "/", mustWork = TRUE
  )
  future$ssp <- as.character(future[[ssp_col]])
  future$period_start <- as.integer(as_num(future[[ps_col]]))
  future$period_end <- as.integer(as_num(future[[pe_col]]))
  future <- future[order(future$ssp, future$period_start), , drop = FALSE]

  for (i in seq_len(nrow(future))) {
    same_grid(
      current, terra::rast(future$path[i]),
      paste0(future$ssp[i], " ", future$period_start[i], "-", future$period_end[i])
    )
  }

  input_manifest <- data.frame(
    input = c(
      "Fixed MTSPS CSV", "Final maxentResults.csv",
      "Current final full-data cloglog raster",
      "Final occurrence CSV", "Future ensemble manifest"
    ),
    path = c(
      normalizePath(fixed_mtsps_csv, winslash = "/", mustWork = TRUE),
      maxent_csv, current_map, occurrences_csv, ensemble_manifest_path
    ),
    stringsAsFactors = FALSE
  )
  write.csv(input_manifest, input_manifest_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # 5.2 Cell-area layer
  # -------------------------------------------------------------------

  if (file.exists(cell_area_tif)) {
    cell_area <- terra::rast(cell_area_tif)
    same_grid(current, cell_area, "existing geodesic area raster")
  } else {
    log_message("Creating geodesic cell-area raster (km²).")
    cell_area <- terra::cellSize(
      current, unit = "km", transform = TRUE, mask = FALSE
    )
    names(cell_area) <- "geodesic_cell_area_km2"
    terra::writeRaster(
      cell_area, cell_area_tif, overwrite = TRUE,
      wopt = list(gdal = "COMPRESS=LZW")
    )
  }

  # -------------------------------------------------------------------
  # 5.3 Current prediction at final occurrences and random background
  # -------------------------------------------------------------------

  occ <- read.csv(occurrences_csv, check.names = FALSE, stringsAsFactors = FALSE)
  lon_col <- find_column(
    occ, c("longitude", "decimalLongitude", "lon", "x"), "final occurrence CSV"
  )
  lat_col <- find_column(
    occ, c("latitude", "decimalLatitude", "lat", "y"), "final occurrence CSV"
  )

  occ$longitude_analysis <- as_num(occ[[lon_col]])
  occ$latitude_analysis <- as_num(occ[[lat_col]])
  occ <- occ[
    is.finite(occ$longitude_analysis) &
      is.finite(occ$latitude_analysis),
    ,
    drop = FALSE
  ]
  occ <- occ[
    !duplicated(paste(occ$longitude_analysis, occ$latitude_analysis, sep = "_")),
    ,
    drop = FALSE
  ]

  if (nrow(occ) != 828L) {
    stop(
      "\nExpected 828 unique final occurrence coordinates, found ", nrow(occ),
      ".", call. = FALSE
    )
  }

  occ_vect <- terra::vect(
    occ,
    geom = c("longitude_analysis", "latitude_analysis"),
    crs = "EPSG:4326"
  )
  p <- terra::extract(current, occ_vect, ID = FALSE)
  p <- as_num(p[[1]])

  if (sum(is.finite(p)) != 828L) {
    stop(
      "\nNot all 828 final occurrences have valid current suitability values.",
      call. = FALSE
    )
  }

  presence_scores <- data.frame(
    record_id = seq_len(828),
    longitude = occ$longitude_analysis,
    latitude = occ$latitude_analysis,
    current_cloglog_suitability = p,
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(threshold_table))) {
    id <- threshold_table$threshold_id[i]
    val <- threshold_table$threshold_cloglog[i]
    presence_scores[[paste0("suitable_", id)]] <- as.integer(p >= val)
  }
  write.csv(presence_scores, presence_scores_csv, row.names = FALSE)

  set.seed(random_seed)
  bg <- terra::spatSample(
    current, size = background_sample_size, method = "random",
    as.df = TRUE, values = TRUE, na.rm = TRUE
  )
  bg <- as_num(bg[[1]])
  bg <- bg[is.finite(bg)]
  if (length(bg) < 1000L) {
    stop("\nInsufficient valid global background sample.", call. = FALSE)
  }
  write.csv(
    data.frame(current_cloglog_suitability = bg),
    background_scores_csv, row.names = FALSE
  )

  psum <- do.call(rbind, lapply(seq_len(nrow(threshold_table)), function(i) {
    val <- threshold_table$threshold_cloglog[i]
    data.frame(
      threshold_id = threshold_table$threshold_id[i],
      threshold_cloglog = val,
      final_occurrence_count = length(p),
      occurrence_min = min(p),
      occurrence_mean = mean(p),
      occurrence_median = stats::median(p),
      occurrence_10th_percentile = as.numeric(stats::quantile(p, 0.10, names = FALSE)),
      occurrence_max = max(p),
      occurrences_at_or_above_threshold = sum(p >= val),
      percent_occurrences_at_or_above_threshold = 100 * mean(p >= val),
      background_sample_count = length(bg),
      background_mean = mean(bg),
      background_median = stats::median(bg),
      background_95th_percentile = as.numeric(stats::quantile(bg, 0.95, names = FALSE)),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(psum, presence_summary_csv, row.names = FALSE)
  make_context_figure(p, bg, threshold_table, threshold_context_pdf)

  # -------------------------------------------------------------------
  # 5.4 Legends
  # -------------------------------------------------------------------

  write.csv(
    data.frame(
      binary_code = c(0, 1),
      binary_class = c("Unsuitable", "Suitable"),
      rule = c(
        "Cloglog suitability < selected threshold",
        "Cloglog suitability >= selected threshold"
      ),
      stringsAsFactors = FALSE
    ),
    binary_legend_csv, row.names = FALSE
  )

  write.csv(
    data.frame(
      transition_code = 1:4,
      transition_class = c(
        "Stable unsuitable", "Expansion (gain)",
        "Contraction (loss)", "Stable suitable"
      ),
      definition = c(
        "Current unsuitable -> Future unsuitable",
        "Current unsuitable -> Future suitable",
        "Current suitable -> Future unsuitable",
        "Current suitable -> Future suitable"
      ),
      stringsAsFactors = FALSE
    ),
    transition_legend_csv, row.names = FALSE
  )

  # -------------------------------------------------------------------
  # 5.5 Binary maps and suitable area
  # -------------------------------------------------------------------

  current_area_list <- list()
  future_area_list <- list()
  binary_manifest_list <- list()
  current_binary_paths <- list()

  for (t in seq_len(nrow(threshold_table))) {
    id <- threshold_table$threshold_id[t]
    val <- threshold_table$threshold_cloglog[t]
    dir_t <- file.path(binary_root, id)
    ensure_dir(dir_t)

    current_bin_path <- file.path(
      dir_t,
      paste0(species_name, "_current_final_full_data_", id, "_binary.tif")
    )

    if (!file.exists(current_bin_path) || overwrite_existing_maps) {
      log_message("Writing current binary map: ", id)
      make_binary(current, val, current_bin_path)
      map_status <- "completed"
    } else {
      map_status <- "skipped_existing_output"
    }

    cbin <- terra::rast(current_bin_path)
    ca <- area_by_code(cbin, cell_area, c(0, 1), c("Unsuitable", "Suitable"))
    curr_suitable <- ca$area_km2[ca$code == 1]
    curr_valid <- sum(ca$area_km2)

    current_area_list[[length(current_area_list) + 1]] <- data.frame(
      threshold_id = id,
      threshold_cloglog = val,
      current_unsuitable_km2 = ca$area_km2[ca$code == 0],
      current_suitable_km2 = curr_suitable,
      current_suitable_10_4_km2 = curr_suitable / 10000,
      current_valid_environment_km2 = curr_valid,
      current_percent_suitable = 100 * curr_suitable / curr_valid,
      current_binary_tif = normalizePath(current_bin_path, winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )

    binary_manifest_list[[length(binary_manifest_list) + 1]] <- data.frame(
      threshold_id = id, map_type = "Current final full-data",
      ssp = NA_character_, period_start = NA_integer_, period_end = NA_integer_,
      binary_tif = normalizePath(current_bin_path, winslash = "/", mustWork = TRUE),
      status = map_status, stringsAsFactors = FALSE
    )
    current_binary_paths[[id]] <- current_bin_path

    for (i in seq_len(nrow(future))) {
      ssp <- future$ssp[i]
      start <- future$period_start[i]
      end <- future$period_end[i]
      future_bin_path <- file.path(
        dir_t,
        paste0(
          species_name, "_future_ensemble_", ssp, "_",
          start, "_", end, "_", id, "_binary.tif"
        )
      )

      if (!file.exists(future_bin_path) || overwrite_existing_maps) {
        log_message("Writing future binary map: ", id, " | ", ssp, " ", start, "-", end)
        make_binary(terra::rast(future$path[i]), val, future_bin_path)
        fstatus <- "completed"
      } else {
        fstatus <- "skipped_existing_output"
      }

      fbin <- terra::rast(future_bin_path)
      fa <- area_by_code(fbin, cell_area, c(0, 1), c("Unsuitable", "Suitable"))
      fut_suitable <- fa$area_km2[fa$code == 1]
      fut_valid <- sum(fa$area_km2)

      future_area_list[[length(future_area_list) + 1]] <- data.frame(
        threshold_id = id,
        threshold_cloglog = val,
        ssp = ssp,
        period_start = start,
        period_end = end,
        scenario = paste0(ssp, " ", start, "-", end),
        current_suitable_km2 = curr_suitable,
        future_unsuitable_km2 = fa$area_km2[fa$code == 0],
        future_suitable_km2 = fut_suitable,
        future_suitable_10_4_km2 = fut_suitable / 10000,
        future_valid_environment_km2 = fut_valid,
        future_percent_suitable = 100 * fut_suitable / fut_valid,
        net_suitable_change_km2 = fut_suitable - curr_suitable,
        percent_change_relative_to_current = 100 * (fut_suitable - curr_suitable) / curr_suitable,
        future_binary_tif = normalizePath(future_bin_path, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
      )

      binary_manifest_list[[length(binary_manifest_list) + 1]] <- data.frame(
        threshold_id = id, map_type = "Future three-GCM ensemble mean",
        ssp = ssp, period_start = start, period_end = end,
        binary_tif = normalizePath(future_bin_path, winslash = "/", mustWork = TRUE),
        status = fstatus, stringsAsFactors = FALSE
      )
    }
  }

  current_area <- do.call(rbind, current_area_list)
  future_area <- do.call(rbind, future_area_list)
  binary_manifest <- do.call(rbind, binary_manifest_list)

  current_area <- current_area[
    match(threshold_order, current_area$threshold_id), , drop = FALSE
  ]
  future_area <- future_area[
    order(match(future_area$threshold_id, threshold_order),
          future_area$ssp, future_area$period_start), ,
    drop = FALSE
  ]

  write.csv(current_area, current_area_csv, row.names = FALSE)
  write.csv(future_area, future_area_long_csv, row.names = FALSE)
  write.csv(binary_manifest, binary_manifest_csv, row.names = FALSE)

  scenarios <- unique(future_area[, c("ssp", "period_start", "period_end", "scenario")])
  scenarios <- scenarios[order(scenarios$ssp, scenarios$period_start), , drop = FALSE]

  wide_future <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    out <- scenarios[i, , drop = FALSE]
    for (id in threshold_order) {
      z <- future_area[
        future_area$threshold_id == id &
          future_area$ssp == out$ssp &
          future_area$period_start == out$period_start,
        ,
        drop = FALSE
      ]
      out[[paste0(id, "_future_suitable_km2")]] <- z$future_suitable_km2
      out[[paste0(id, "_net_change_km2")]] <- z$net_suitable_change_km2
      out[[paste0(id, "_percent_change")]] <- z$percent_change_relative_to_current
    }
    out
  }))
  write.csv(wide_future, future_area_wide_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # 5.6 Transition maps and transition areas
  # -------------------------------------------------------------------

  change_list <- list()
  transition_manifest_list <- list()

  for (t in seq_len(nrow(threshold_table))) {
    id <- threshold_table$threshold_id[t]
    val <- threshold_table$threshold_cloglog[t]
    dir_t <- file.path(transition_root, id)
    ensure_dir(dir_t)
    cbin <- terra::rast(current_binary_paths[[id]])

    for (i in seq_len(nrow(future))) {
      ssp <- future$ssp[i]
      start <- future$period_start[i]
      end <- future$period_end[i]

      fbin_path <- file.path(
        binary_root, id,
        paste0(
          species_name, "_future_ensemble_", ssp, "_",
          start, "_", end, "_", id, "_binary.tif"
        )
      )

      trans_path <- file.path(
        dir_t,
        paste0(
          species_name, "_transition_current_to_", ssp, "_",
          start, "_", end, "_", id, ".tif"
        )
      )

      if (!file.exists(trans_path) || overwrite_existing_maps) {
        log_message("Writing transition map: ", id, " | ", ssp, " ", start, "-", end)
        fbin <- terra::rast(fbin_path)
        same_grid(cbin, fbin, paste(id, ssp, start, end))
        make_transition(cbin, fbin, trans_path)
        tstatus <- "completed"
      } else {
        tstatus <- "skipped_existing_output"
      }

      tr <- terra::rast(trans_path)
      ta <- area_by_code(
        tr, cell_area, 1:4,
        c("Stable unsuitable", "Expansion (gain)",
          "Contraction (loss)", "Stable suitable")
      )

      get_area <- function(code) ta$area_km2[ta$code == code]
      stable_unsuitable <- get_area(1)
      expansion <- get_area(2)
      contraction <- get_area(3)
      stable_suitable <- get_area(4)
      current_suitable <- contraction + stable_suitable
      future_suitable <- expansion + stable_suitable

      change_list[[length(change_list) + 1]] <- data.frame(
        threshold_id = id,
        threshold_cloglog = val,
        ssp = ssp,
        period_start = start,
        period_end = end,
        scenario = paste0(ssp, " ", start, "-", end),
        stable_unsuitable_km2 = stable_unsuitable,
        expansion_gain_km2 = expansion,
        contraction_loss_km2 = contraction,
        stable_suitable_km2 = stable_suitable,
        current_suitable_km2 = current_suitable,
        future_suitable_km2 = future_suitable,
        net_suitable_change_km2 = future_suitable - current_suitable,
        net_suitable_change_10_4_km2 = (future_suitable - current_suitable) / 10000,
        percent_change_relative_to_current_suitable =
          100 * (future_suitable - current_suitable) / current_suitable,
        expansion_percent_of_current_unsuitable =
          100 * expansion / (stable_unsuitable + expansion),
        contraction_percent_of_current_suitable =
          100 * contraction / current_suitable,
        shared_valid_domain_km2 = sum(ta$area_km2),
        transition_tif = normalizePath(trans_path, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
      )

      transition_manifest_list[[length(transition_manifest_list) + 1]] <- data.frame(
        threshold_id = id,
        ssp = ssp,
        period_start = start,
        period_end = end,
        transition_tif = normalizePath(trans_path, winslash = "/", mustWork = TRUE),
        status = tstatus,
        stringsAsFactors = FALSE
      )
    }
  }

  change_area <- do.call(rbind, change_list)
  transition_manifest <- do.call(rbind, transition_manifest_list)

  change_area <- change_area[
    order(match(change_area$threshold_id, threshold_order),
          change_area$ssp, change_area$period_start),
    ,
    drop = FALSE
  ]

  write.csv(change_area, change_area_long_csv, row.names = FALSE)
  write.csv(transition_manifest, transition_manifest_csv, row.names = FALSE)

  wide_change <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
    out <- scenarios[i, , drop = FALSE]
    for (id in threshold_order) {
      z <- change_area[
        change_area$threshold_id == id &
          change_area$ssp == out$ssp &
          change_area$period_start == out$period_start,
        ,
        drop = FALSE
      ]
      for (metric in c(
        "expansion_gain_km2", "contraction_loss_km2",
        "stable_suitable_km2", "net_suitable_change_km2",
        "percent_change_relative_to_current_suitable"
      )) {
        out[[paste0(id, "_", metric)]] <- z[[metric]]
      }
    }
    out
  }))
  write.csv(wide_change, change_area_wide_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # 5.7 Balance checks and threshold summary
  # -------------------------------------------------------------------

  balance <- do.call(rbind, lapply(seq_len(nrow(change_area)), function(i) {
    z <- change_area[i, , drop = FALSE]
    curr <- current_area[current_area$threshold_id == z$threshold_id, , drop = FALSE]
    fut <- future_area[
      future_area$threshold_id == z$threshold_id &
        future_area$ssp == z$ssp &
        future_area$period_start == z$period_start,
      ,
      drop = FALSE
    ]
    data.frame(
      threshold_id = z$threshold_id,
      ssp = z$ssp,
      period_start = z$period_start,
      period_end = z$period_end,
      current_suitable_balance_error_km2 =
        curr$current_suitable_km2 -
        (z$contraction_loss_km2 + z$stable_suitable_km2),
      future_suitable_balance_error_km2 =
        fut$future_suitable_km2 -
        (z$expansion_gain_km2 + z$stable_suitable_km2),
      net_change_balance_error_km2 =
        z$net_suitable_change_km2 -
        (z$expansion_gain_km2 - z$contraction_loss_km2),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(balance, balance_csv, row.names = FALSE)

  comparison <- do.call(rbind, lapply(threshold_order, function(id) {
    z <- change_area[change_area$threshold_id == id, , drop = FALSE]
    curr <- current_area[current_area$threshold_id == id, , drop = FALSE]
    data.frame(
      threshold_id = id,
      threshold_cloglog = curr$threshold_cloglog,
      current_suitable_km2 = curr$current_suitable_km2,
      current_percent_suitable = curr$current_percent_suitable,
      minimum_future_suitable_km2 = min(z$future_suitable_km2),
      maximum_future_suitable_km2 = max(z$future_suitable_km2),
      minimum_net_change_km2 = min(z$net_suitable_change_km2),
      maximum_net_change_km2 = max(z$net_suitable_change_km2),
      scenarios_with_net_gain = sum(z$net_suitable_change_km2 > 0),
      scenarios_with_net_loss = sum(z$net_suitable_change_km2 < 0),
      scenarios_expansion_gt_contraction = sum(z$expansion_gain_km2 > z$contraction_loss_km2),
      scenarios_contraction_gt_expansion = sum(z$contraction_loss_km2 > z$expansion_gain_km2),
      mean_stable_suitable_percent_of_current = mean(
        100 * z$stable_suitable_km2 / curr$current_suitable_km2
      ),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(comparison, comparison_csv, row.names = FALSE)

  max_balance <- max(
    abs(unlist(balance[
      c(
        "current_suitable_balance_error_km2",
        "future_suitable_balance_error_km2",
        "net_change_balance_error_km2"
      )
    ])),
    na.rm = TRUE
  )

  run_summary <- data.frame(
    species = species_name,
    baseline = "Current final full-data MaxEnt model",
    primary_threshold = "MTSPS",
    MTSPS_cloglog = threshold_table$threshold_cloglog[threshold_table$threshold_id == "MTSPS"],
    MTP_cloglog = threshold_table$threshold_cloglog[threshold_table$threshold_id == "MTP"],
    T10_cloglog = threshold_table$threshold_cloglog[threshold_table$threshold_id == "T10"],
    current_binary_maps = nrow(current_area),
    future_binary_maps = nrow(future_area),
    transition_maps = nrow(change_area),
    maximum_area_balance_error_km2 = max_balance,
    stringsAsFactors = FALSE
  )
  write.csv(run_summary, run_summary_csv, row.names = FALSE)

  writeLines(
    c(
      "SCRIPT 17 COMPLETED SUCCESSFULLY",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("MTSPS:", format(run_summary$MTSPS_cloglog, scientific = TRUE, digits = 8)),
      paste("MTP:", format(run_summary$MTP_cloglog, scientific = TRUE, digits = 8)),
      paste("T10:", format(run_summary$T10_cloglog, scientific = TRUE, digits = 8)),
      paste("Maximum area-balance error (km²):", format(max_balance, scientific = TRUE, digits = 8))
    ),
    completed_marker
  )

  log_message("Script completed successfully.")

  cat("\n======================================================================\n")
  cat("SCRIPT 17 COMPLETED — THRESHOLD SENSITIVITY ANALYSIS\n")
  cat("======================================================================\n")
  cat(
    "MTSPS: ",
    format(run_summary$MTSPS_cloglog, scientific = TRUE, digits = 8),
    "\n", sep = ""
  )
  cat(
    "MTP: ",
    format(run_summary$MTP_cloglog, scientific = TRUE, digits = 8),
    "\n", sep = ""
  )
  cat(
    "T10: ",
    format(run_summary$T10_cloglog, scientific = TRUE, digits = 8),
    "\n", sep = ""
  )
  cat("Current binary maps: ", nrow(current_area), "\n", sep = "")
  cat("Future binary maps: ", nrow(future_area), "\n", sep = "")
  cat("Transition maps: ", nrow(change_area), "\n", sep = "")
  cat(
    "Maximum area-balance error (km²): ",
    format(max_balance, scientific = TRUE, digits = 8),
    "\n", sep = ""
  )
  cat("\nKey tables:\n")
  cat(thresholds_csv, "\n")
  cat(presence_summary_csv, "\n")
  cat(future_area_long_csv, "\n")
  cat(change_area_long_csv, "\n")
  cat(comparison_csv, "\n")
}

tryCatch(
  run_script17(),
  error = function(e) {
    dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    message_text <- paste0(
      "SCRIPT 17 FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(e), "\n\n",
      "Execution log: ", execution_log, "\n"
    )
    writeLines(message_text, failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 17 FAILED\n")
    cat("======================================================================\n")
    cat(message_text, "\n")

    stop(e)
  }
)
