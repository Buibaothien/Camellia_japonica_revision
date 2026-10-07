# =====================================================================
# 12_ENSEMBLE_3GCM_FUTURE_CAMELLIA_JAPONICA.R
# =====================================================================
# Create 16 continuous ensemble future suitability rasters from the
# successfully completed 48 individual GCM projections:
# 3 GCMs x 4 SSPs x 4 time periods.
#
# The ensemble is calculated as the arithmetic mean of continuous cloglog
# suitability values from ACCESS-CM2, BCC-CSM2-MR, and CMCC-ESM2.
#
# Only cells with finite predictions from ALL three GCMs are retained.
# This script does NOT classify suitability or calculate areas.
# Run Script 13 after this script completes all 16 ensembles.
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
suppressPackageStartupMessages(library(terra))

# ---------------------------------------------------------------------
# 1. USER SETTINGS
# ---------------------------------------------------------------------

species_name <- "Camellia_japonica"
project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# Successful Script 11 output folder.
individual_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "11_future_projection_individual_GCM_v2_plotsfix"
)

# New output folder created by this script.
ensemble_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "12_future_ensemble_3GCM"
)

# FALSE retains any completed ensemble outputs and safely resumes.
clean_ensemble_root <- FALSE
skip_existing_ensemble <- TRUE

expected_gcms <- c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2")

# ---------------------------------------------------------------------
# 2. PATHS
# ---------------------------------------------------------------------

projection_manifest_csv <- file.path(
  individual_root,
  "04_projection_manifests",
  "future_projection_status_manifest.csv"
)

mean_root <- file.path(ensemble_root, "01_ensemble_mean_cloglog")
sd_root <- file.path(ensemble_root, "02_interGCM_standard_deviation")
count_root <- file.path(ensemble_root, "03_valid_GCM_count")
manifest_root <- file.path(ensemble_root, "04_manifests")
log_root <- file.path(ensemble_root, "05_logs")

ensemble_manifest_csv <- file.path(
  manifest_root,
  "ensemble_status_manifest.csv"
)
source_manifest_csv <- file.path(
  manifest_root,
  "ensemble_source_manifest.csv"
)
summary_csv <- file.path(
  manifest_root,
  "ensemble_run_summary.csv"
)

# ---------------------------------------------------------------------
# 3. HELPERS
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
  names_norm <- normalise_name(names_original)
  candidates_norm <- normalise_name(candidates)

  hit <- match(candidates_norm, names_norm, nomatch = 0)
  hit <- hit[hit > 0]

  if (length(hit) > 0) {
    return(names_original[hit[1]])
  }

  for (candidate in candidates_norm) {
    partial <- grep(candidate, names_norm, fixed = TRUE)
    if (length(partial) > 0) return(names_original[partial[1]])
  }

  stop(
    "\nRequired column missing from ", label, ".\nAvailable columns: ",
    paste(names_original, collapse = ", "),
    call. = FALSE
  )
}

assert_same_grid <- function(reference, candidate, label) {
  ok <- terra::compareGeom(
    reference, candidate,
    crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
    stopOnError = FALSE
  )

  if (!ok) {
    stop(
      "\nGrid mismatch between future GCM rasters for ", label, ".\n",
      "Do not continue: all 3 projections should have the same global grid.",
      call. = FALSE
    )
  }
}

raster_summary <- function(r) {
  rng <- terra::global(r, fun = c("min", "max"), na.rm = TRUE)
  valid <- terra::global(!is.na(r), fun = "sum", na.rm = TRUE)

  list(
    min = as.numeric(rng[1, "min"]),
    max = as.numeric(rng[1, "max"]),
    valid_cells = as.numeric(valid[1, 1])
  )
}

strict_mean <- function(x) {
  if (all(is.finite(x))) mean(x) else NA_real_
}

strict_sd <- function(x) {
  if (all(is.finite(x))) stats::sd(x) else NA_real_
}

valid_count <- function(x) sum(is.finite(x))

write_progress <- function(records, source_records) {
  if (length(records) > 0) {
    write.csv(
      do.call(rbind, records),
      ensemble_manifest_csv,
      row.names = FALSE
    )
  }

  if (length(source_records) > 0) {
    write.csv(
      do.call(rbind, source_records),
      source_manifest_csv,
      row.names = FALSE
    )
  }
}

# ---------------------------------------------------------------------
# 4. VERIFY SCRIPT 11 OUTPUT
# ---------------------------------------------------------------------

if (!file.exists(projection_manifest_csv)) {
  stop(
    "\nScript 11 manifest was not found:\n",
    projection_manifest_csv,
    "\nRun Script 11 successfully before Script 12.",
    call. = FALSE
  )
}

projection_manifest <- read.csv(
  projection_manifest_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

gcm_col <- find_column(projection_manifest, c("gcm"), "Script 11 manifest")
ssp_col <- find_column(projection_manifest, c("ssp"), "Script 11 manifest")
start_col <- find_column(
  projection_manifest,
  c("period_start", "periodstart"),
  "Script 11 manifest"
)
end_col <- find_column(
  projection_manifest,
  c("period_end", "periodend"),
  "Script 11 manifest"
)
status_col <- find_column(projection_manifest, c("status"), "Script 11 manifest")
output_col <- find_column(
  projection_manifest,
  c("output_continuous_tif", "outputcontinuoustif"),
  "Script 11 manifest"
)

usable <- projection_manifest[
  projection_manifest[[status_col]] %in% c("completed", "skipped_existing_output") &
    !is.na(projection_manifest[[output_col]]) &
    nzchar(projection_manifest[[output_col]]),
  ,
  drop = FALSE
]

if (nrow(usable) != 48) {
  stop(
    "\nExpected 48 completed individual GCM projections, but found ",
    nrow(usable), ".\nResolve Script 11 outputs before creating the ensemble.",
    call. = FALSE
  )
}

if (any(!file.exists(usable[[output_col]]))) {
  missing <- usable[!file.exists(usable[[output_col]]), c(
    gcm_col, ssp_col, start_col, end_col, output_col
  ), drop = FALSE]
  print(missing)
  stop("\nSome completed Script 11 output GeoTIFFs are missing.", call. = FALSE)
}

# ---------------------------------------------------------------------
# 5. PREPARE OUTPUT FOLDERS
# ---------------------------------------------------------------------

if (clean_ensemble_root) safe_remove(ensemble_root)

for (folder in c(mean_root, sd_root, count_root, manifest_root, log_root)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

cat("\n======================================================================\n")
cat("THREE-GCM CONTINUOUS ENSEMBLE — Camellia japonica\n")
cat("======================================================================\n")
cat("Input GCM rasters: ", nrow(usable), "\n", sep = "")
cat("Expected ensemble rasters: 16\n")
cat("Output root:\n", normalizePath(ensemble_root, winslash = "/"), "\n", sep = "")

# ---------------------------------------------------------------------
# 6. BUILD 16 ENSEMBLES
# ---------------------------------------------------------------------

keys <- unique(data.frame(
  ssp = as.character(usable[[ssp_col]]),
  period_start = as.integer(usable[[start_col]]),
  period_end = as.integer(usable[[end_col]]),
  stringsAsFactors = FALSE
))

keys <- keys[order(keys$ssp, keys$period_start), , drop = FALSE]

if (nrow(keys) != 16L) {
  stop(
    "\nExpected 16 SSP x period combinations but found ", nrow(keys), ".",
    call. = FALSE
  )
}

status_records <- list()
source_records <- list()
run_start <- Sys.time()

for (i in seq_len(nrow(keys))) {
  key <- keys[i, , drop = FALSE]
  ssp <- key$ssp
  p_start <- key$period_start
  p_end <- key$period_end
  scenario_id <- paste0(tolower(ssp), "_", p_start, "_", p_end)

  cat("\n[", i, "/16] ", ssp, " ", p_start, "-", p_end, "\n", sep = "")

  rows <- usable[
    as.character(usable[[ssp_col]]) == ssp &
      as.integer(usable[[start_col]]) == p_start &
      as.integer(usable[[end_col]]) == p_end,
    ,
    drop = FALSE
  ]

  present <- as.character(rows[[gcm_col]])
  missing_gcms <- setdiff(expected_gcms, present)

  if (nrow(rows) != 3L || length(missing_gcms) > 0) {
    stop(
      "\nIncomplete GCM set for ", ssp, " ", p_start, "-", p_end, ".\n",
      "Present: ", paste(present, collapse = ", "),
      "\nMissing: ", paste(missing_gcms, collapse = ", "),
      call. = FALSE
    )
  }

  rows <- rows[match(expected_gcms, as.character(rows[[gcm_col]])), , drop = FALSE]
  source_files <- as.character(rows[[output_col]])

  out_mean_dir <- file.path(mean_root, scenario_id)
  out_sd_dir <- file.path(sd_root, scenario_id)
  out_count_dir <- file.path(count_root, scenario_id)

  for (folder in c(out_mean_dir, out_sd_dir, out_count_dir)) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  }

  mean_tif <- file.path(
    out_mean_dir,
    paste0(species_name, "_ensemble_", scenario_id, "_cloglog_mean.tif")
  )
  sd_tif <- file.path(
    out_sd_dir,
    paste0(species_name, "_ensemble_", scenario_id, "_cloglog_sd.tif")
  )
  count_tif <- file.path(
    out_count_dir,
    paste0(species_name, "_ensemble_", scenario_id, "_valid_GCM_count.tif")
  )
  source_csv <- file.path(out_mean_dir, "source_GCM_rasters.csv")
  scenario_log <- file.path(log_root, paste0("ensemble_", scenario_id, ".txt"))

  scenario_start <- Sys.time()

  if (
    skip_existing_ensemble &&
    file.exists(mean_tif) &&
    file.exists(sd_tif) &&
    file.exists(count_tif)
  ) {
    cat("  Status: existing ensemble outputs skipped.\n")

    status_records[[length(status_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = p_start,
      period_end = p_end,
      ensemble_id = scenario_id,
      status = "skipped_existing_output",
      ensemble_mean_tif = normalizePath(mean_tif, winslash = "/", mustWork = TRUE),
      ensemble_sd_tif = normalizePath(sd_tif, winslash = "/", mustWork = TRUE),
      valid_GCM_count_tif = normalizePath(count_tif, winslash = "/", mustWork = TRUE),
      min_mean_suitability = NA_real_,
      max_mean_suitability = NA_real_,
      mean_inter_GCM_sd = NA_real_,
      n_cells_with_all_3_GCMs = NA_real_,
      elapsed_minutes = 0,
      message = "",
      stringsAsFactors = FALSE
    )

    source_records[[length(source_records) + 1]] <- data.frame(
      ssp = ssp,
      period_start = p_start,
      period_end = p_end,
      gcm = expected_gcms,
      source_continuous_tif = normalizePath(source_files, winslash = "/", mustWork = TRUE),
      stringsAsFactors = FALSE
    )

    write_progress(status_records, source_records)
    next
  }

  result <- tryCatch(
    {
      gcm_rasters <- lapply(source_files, terra::rast)

      for (j in seq_along(gcm_rasters)) {
        assert_same_grid(
          gcm_rasters[[1]],
          gcm_rasters[[j]],
          expected_gcms[j]
        )
      }

      gcm_summaries <- lapply(gcm_rasters, raster_summary)

      for (j in seq_along(gcm_summaries)) {
        info <- gcm_summaries[[j]]
        if (!is.finite(info$min) || !is.finite(info$max)) {
          stop("No finite prediction values in ", expected_gcms[j], ".", call. = FALSE)
        }
        if (info$min < -1e-8 || info$max > 1 + 1e-8) {
          warning(
            "Suitability values outside [0,1] detected in ",
            expected_gcms[j], " for ", scenario_id
          )
        }
      }

      x <- terra::rast(gcm_rasters)
      names(x) <- expected_gcms

      # Strict ensemble mean and SD: only use cells with all 3 finite GCM values.
      mean_raster <- terra::app(x, fun = strict_mean)
      sd_raster <- terra::app(x, fun = strict_sd)
      count_raster <- terra::app(x, fun = valid_count)

      names(mean_raster) <- "ensemble_mean_cloglog"
      names(sd_raster) <- "inter_GCM_standard_deviation"
      names(count_raster) <- "valid_GCM_count"

      terra::writeRaster(
        mean_raster,
        mean_tif,
        overwrite = TRUE,
        wopt = list(gdal = "COMPRESS=LZW")
      )
      terra::writeRaster(
        sd_raster,
        sd_tif,
        overwrite = TRUE,
        wopt = list(gdal = "COMPRESS=LZW")
      )
      terra::writeRaster(
        count_raster,
        count_tif,
        overwrite = TRUE,
        wopt = list(gdal = "COMPRESS=LZW")
      )

      mean_info <- raster_summary(mean_raster)
      mean_sd <- terra::global(sd_raster, "mean", na.rm = TRUE)
      all_three_cells <- terra::global(count_raster == 3, "sum", na.rm = TRUE)

      sources <- data.frame(
        ssp = ssp,
        period_start = p_start,
        period_end = p_end,
        gcm = expected_gcms,
        source_continuous_tif = normalizePath(source_files, winslash = "/", mustWork = TRUE),
        source_min_suitability = vapply(gcm_summaries, function(z) z$min, numeric(1)),
        source_max_suitability = vapply(gcm_summaries, function(z) z$max, numeric(1)),
        source_valid_cells = vapply(gcm_summaries, function(z) z$valid_cells, numeric(1)),
        stringsAsFactors = FALSE
      )

      write.csv(sources, source_csv, row.names = FALSE)

      writeLines(
        c(
          paste0("Species: ", species_name),
          paste0("Scenario: ", ssp, " ", p_start, "-", p_end),
          paste0("GCMs: ", paste(expected_gcms, collapse = "; ")),
          "Ensemble: strict arithmetic mean across all three continuous cloglog rasters",
          "NoData rule: NoData if one or more GCM predictions are NoData",
          paste0("Mean output: ", normalizePath(mean_tif, winslash = "/")),
          paste0("SD output: ", normalizePath(sd_tif, winslash = "/")),
          paste0("Count output: ", normalizePath(count_tif, winslash = "/"))
        ),
        scenario_log
      )

      list(
        status = "completed",
        mean_tif = normalizePath(mean_tif, winslash = "/", mustWork = TRUE),
        sd_tif = normalizePath(sd_tif, winslash = "/", mustWork = TRUE),
        count_tif = normalizePath(count_tif, winslash = "/", mustWork = TRUE),
        min_mean = mean_info$min,
        max_mean = mean_info$max,
        mean_sd = as.numeric(mean_sd[1, 1]),
        all_three_cells = as.numeric(all_three_cells[1, 1]),
        sources = sources,
        message = ""
      )
    },
    error = function(e) {
      list(
        status = "failed",
        mean_tif = NA_character_,
        sd_tif = NA_character_,
        count_tif = NA_character_,
        min_mean = NA_real_,
        max_mean = NA_real_,
        mean_sd = NA_real_,
        all_three_cells = NA_real_,
        sources = data.frame(),
        message = conditionMessage(e)
      )
    }
  )

  elapsed <- as.numeric(difftime(Sys.time(), scenario_start, units = "mins"))

  if (nrow(result$sources) > 0) {
    source_records[[length(source_records) + 1]] <- result$sources
  }

  status_records[[length(status_records) + 1]] <- data.frame(
    ssp = ssp,
    period_start = p_start,
    period_end = p_end,
    ensemble_id = scenario_id,
    status = result$status,
    ensemble_mean_tif = result$mean_tif,
    ensemble_sd_tif = result$sd_tif,
    valid_GCM_count_tif = result$count_tif,
    min_mean_suitability = result$min_mean,
    max_mean_suitability = result$max_mean,
    mean_inter_GCM_sd = result$mean_sd,
    n_cells_with_all_3_GCMs = result$all_three_cells,
    elapsed_minutes = elapsed,
    message = result$message,
    stringsAsFactors = FALSE
  )

  write_progress(status_records, source_records)

  if (result$status == "completed") {
    cat("  Status: completed in ", format(round(elapsed, 2), nsmall = 2), " minutes.\n", sep = "")
  } else {
    stop(
      "\nEnsemble failed for ", scenario_id, ".\nRead the manifest:\n",
      ensemble_manifest_csv,
      call. = FALSE
    )
  }
}

# ---------------------------------------------------------------------
# 7. FINAL REPORT
# ---------------------------------------------------------------------

status <- do.call(rbind, status_records)
write_progress(status_records, source_records)

completed <- sum(status$status == "completed")
skipped <- sum(status$status == "skipped_existing_output")
failed <- sum(status$status == "failed")
runtime_hours <- as.numeric(difftime(Sys.time(), run_start, units = "hours"))

summary <- data.frame(
  expected_ensemble_rasters = 16,
  completed_this_run = completed,
  skipped_existing = skipped,
  failed = failed,
  total_completed_or_existing = completed + skipped,
  runtime_hours = runtime_hours,
  stringsAsFactors = FALSE
)

write.csv(summary, summary_csv, row.names = FALSE)

cat("\n======================================================================\n")
cat("THREE-GCM ENSEMBLE FINISHED — Camellia japonica\n")
cat("======================================================================\n")
cat("Completed this run: ", completed, "\n", sep = "")
cat("Skipped existing outputs: ", skipped, "\n", sep = "")
cat("Failed ensembles: ", failed, "\n", sep = "")
cat("Total completed/existing: ", completed + skipped, " of 16\n", sep = "")
cat("Status manifest:\n", normalizePath(ensemble_manifest_csv, winslash = "/"), "\n", sep = "")
cat("Run summary:\n", normalizePath(summary_csv, winslash = "/"), "\n", sep = "")

if (completed + skipped == 16 && failed == 0) {
  cat("\nAll 16 continuous ensemble rasters are ready for Script 13.\n")
} else {
  cat("\nDo not run Script 13 until all 16 ensemble means are completed.\n")
}
