# =====================================================================
# 13_CLASSIFY_MTSPS_4CLASSES_AND_AREA_CAMELLIA_JAPONICA.R
# =====================================================================
# Apply the single MTSPS threshold from the final current model to:
#   - 16 three-GCM continuous future ensemble rasters;
#   - current CV100 average and current final full-data continuous maps.
#
# Four classes:
#   1 = Unsuitable:           p < MTSPS
#   2 = Low suitability:      MTSPS <= p < 0.4
#   3 = Moderate suitability: 0.4 <= p < 0.6
#   4 = High suitability:     p >= 0.6
#
# Geographic areas are calculated geodesically in km² and 10^4 km².
# =====================================================================

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

current_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

ensemble_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "12_future_ensemble_3GCM"
)

output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "13_MTSPS_4class_future_ensemble_area"
)

clean_output_root <- FALSE
skip_existing_classified_maps <- TRUE
stop_after_first_failure <- TRUE

# ---------------------------------------------------------------------
# 2. INPUT / OUTPUT PATHS
# ---------------------------------------------------------------------

threshold_csv <- file.path(
  current_root,
  "04_CV_metrics_MTSPS_and_areas",
  "03_FINAL_single_MTSPS_threshold.csv"
)

current_cv100_tif <- file.path(
  current_root,
  "05_current_suitability_maps",
  "Camellia_japonica_current_CV100_average_cloglog.tif"
)

current_final_tif <- file.path(
  current_root,
  "05_current_suitability_maps",
  "Camellia_japonica_current_final_full_data_cloglog.tif"
)

ensemble_manifest_csv <- file.path(
  ensemble_root,
  "04_manifests",
  "ensemble_status_manifest.csv"
)

current_output_dir <- file.path(output_root, "01_current_4class_maps")
future_output_root <- file.path(output_root, "02_future_ensemble_4class_maps")
area_dir <- file.path(output_root, "03_area_tables")
manifest_dir <- file.path(output_root, "04_manifests")

status_csv <- file.path(manifest_dir, "classification_status_manifest.csv")
legend_csv <- file.path(area_dir, "01_suitability_4class_legend.csv")
all_areas_csv <- file.path(area_dir, "02_four_class_area_all_current_and_future_maps.csv")
future_wide_csv <- file.path(area_dir, "03_future_ensemble_four_class_area_wide.csv")
current_areas_csv <- file.path(area_dir, "04_current_four_class_area.csv")
suitable_summary_csv <- file.path(area_dir, "05_total_suitable_area_summary.csv")
run_summary_csv <- file.path(manifest_dir, "classification_run_summary.csv")

# ---------------------------------------------------------------------
# 3. FUNCTIONS
# ---------------------------------------------------------------------

safe_remove <- function(path) {
  if (file.exists(path) || dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
}

normalise_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

find_column <- function(data, candidates, label) {
  nm <- names(data)
  nm_norm <- normalise_name(nm)
  cand_norm <- normalise_name(candidates)

  hit <- match(cand_norm, nm_norm, nomatch = 0)
  hit <- hit[hit > 0]
  if (length(hit) > 0) return(nm[hit[1]])

  for (candidate in cand_norm) {
    partial <- grep(candidate, nm_norm, fixed = TRUE)
    if (length(partial) > 0) return(nm[partial[1]])
  }

  stop(
    "\nCannot find required column in ", label,
    ". Available columns: ", paste(nm, collapse = ", "),
    call. = FALSE
  )
}

read_mtsps <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  col <- find_column(d, c("MTSPS_threshold", "threshold"), "MTSPS threshold CSV")
  val <- suppressWarnings(as.numeric(d[[col]]))
  val <- val[is.finite(val)]

  if (length(val) != 1) {
    stop("MTSPS threshold CSV must contain exactly one finite threshold.", call. = FALSE)
  }
  if (val[1] < 0 || val[1] >= 0.4) {
    stop(
      "MTSPS threshold must be >= 0 and < 0.4 for the requested four-class scheme. Found: ",
      val[1],
      call. = FALSE
    )
  }
  val[1]
}

check_continuous <- function(r, label) {
  rng <- terra::global(r, c("min", "max"), na.rm = TRUE)
  mn <- as.numeric(rng[1, "min"])
  mx <- as.numeric(rng[1, "max"])

  if (!is.finite(mn) || !is.finite(mx)) {
    stop("No finite suitability values in ", label, call. = FALSE)
  }
  if (mn < -1e-8 || mx > 1 + 1e-8) {
    warning("Suitability values outside [0,1] in ", label, ": ", mn, " to ", mx)
  }
  invisible(list(min = mn, max = mx))
}

classify_4 <- function(r, mtsps) {
  out <- terra::ifel(
    is.na(r), NA,
    terra::ifel(
      r < mtsps, 1,
      terra::ifel(r < 0.4, 2, terra::ifel(r < 0.6, 3, 4))
    )
  )
  names(out) <- "suitability_class"
  out
}

set_categories <- function(r, mtsps) {
  levels(r) <- data.frame(
    ID = 1:4,
    class_name = c(
      "Unsuitable",
      "Low suitability",
      "Moderate suitability",
      "High suitability"
    ),
    class_definition = c(
      paste0("Suitability < ", format(mtsps, scientific = FALSE, trim = TRUE)),
      paste0(format(mtsps, scientific = FALSE, trim = TRUE), " <= Suitability < 0.4"),
      "0.4 <= Suitability < 0.6",
      "Suitability >= 0.6"
    ),
    stringsAsFactors = FALSE
  )
  r
}

area_table <- function(
  classified, map_id, map_type, source_model, scenario,
  ssp, period_start, period_end, mtsps
) {
  observed <- terra::expanse(
    classified,
    unit = "km",
    transform = TRUE,
    byValue = TRUE
  )

  value_col <- find_column(observed, c("value"), paste0(map_id, " area output"))
  area_col <- find_column(observed, c("area"), paste0(map_id, " area output"))

  lookup <- data.frame(
    class_code = 1:4,
    class_name = c(
      "Unsuitable",
      "Low suitability",
      "Moderate suitability",
      "High suitability"
    ),
    class_definition = c(
      paste0("Suitability < ", format(mtsps, scientific = FALSE, trim = TRUE)),
      paste0(format(mtsps, scientific = FALSE, trim = TRUE), " <= Suitability < 0.4"),
      "0.4 <= Suitability < 0.6",
      "Suitability >= 0.6"
    ),
    stringsAsFactors = FALSE
  )

  observed <- data.frame(
    class_code = as.integer(observed[[value_col]]),
    area_km2 = as.numeric(observed[[area_col]]),
    stringsAsFactors = FALSE
  )

  out <- merge(lookup, observed, by = "class_code", all.x = TRUE, sort = TRUE)
  out$area_km2[is.na(out$area_km2)] <- 0
  out$area_10_4_km2 <- out$area_km2 / 10000
  total <- sum(out$area_km2)
  out$percent_of_valid_environment <- if (total > 0) 100 * out$area_km2 / total else NA_real_

  out$map_id <- map_id
  out$map_type <- map_type
  out$source_model <- source_model
  out$scenario <- scenario
  out$ssp <- ssp
  out$period_start <- period_start
  out$period_end <- period_end
  out$MTSPS_threshold <- mtsps

  out[, c(
    "map_id", "map_type", "source_model", "scenario", "ssp",
    "period_start", "period_end", "class_code", "class_name",
    "class_definition", "MTSPS_threshold", "area_km2", "area_10_4_km2",
    "percent_of_valid_environment"
  )]
}

write_status <- function(records) {
  if (length(records) > 0) write.csv(do.call(rbind, records), status_csv, row.names = FALSE)
}

make_future_wide <- function(long_data) {
  keys <- unique(long_data[, c("map_id", "scenario", "ssp", "period_start", "period_end", "MTSPS_threshold")])
  keys <- keys[order(keys$ssp, keys$period_start), , drop = FALSE]

  out <- lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, , drop = FALSE]
    x <- long_data[long_data$map_id == k$map_id, , drop = FALSE]

    get_area <- function(class_name) {
      z <- x$area_km2[x$class_name == class_name]
      if (length(z) == 0) 0 else z[1]
    }

    unsuitable <- get_area("Unsuitable")
    low <- get_area("Low suitability")
    moderate <- get_area("Moderate suitability")
    high <- get_area("High suitability")
    total <- unsuitable + low + moderate + high
    suitable <- low + moderate + high

    data.frame(
      scenario = k$scenario,
      ssp = k$ssp,
      period_start = k$period_start,
      period_end = k$period_end,
      MTSPS_threshold = k$MTSPS_threshold,
      unsuitable_km2 = unsuitable,
      low_suitability_km2 = low,
      moderate_suitability_km2 = moderate,
      high_suitability_km2 = high,
      total_suitable_km2 = suitable,
      total_suitable_10_4_km2 = suitable / 10000,
      total_valid_environment_km2 = total,
      percent_suitable = if (total > 0) 100 * suitable / total else NA_real_,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, out)
}

# ---------------------------------------------------------------------
# 4. VALIDATE INPUTS
# ---------------------------------------------------------------------

needed <- c(
  MTSPS_threshold = threshold_csv,
  current_CV100_average = current_cv100_tif,
  current_final_full_data = current_final_tif,
  ensemble_manifest = ensemble_manifest_csv
)

missing <- needed[!file.exists(needed)]
if (length(missing) > 0) {
  stop("\nMissing input(s):\n", paste(names(missing), ": ", missing, collapse = "\n"), call. = FALSE)
}

mtsps <- read_mtsps(threshold_csv)

manifest <- read.csv(ensemble_manifest_csv, check.names = FALSE, stringsAsFactors = FALSE)
status_col <- find_column(manifest, c("status"), "Script 12 ensemble manifest")
mean_col <- find_column(manifest, c("ensemble_mean_tif", "ensemblemeantif"), "Script 12 ensemble manifest")
ssp_col <- find_column(manifest, c("ssp"), "Script 12 ensemble manifest")
start_col <- find_column(manifest, c("period_start", "periodstart"), "Script 12 ensemble manifest")
end_col <- find_column(manifest, c("period_end", "periodend"), "Script 12 ensemble manifest")
id_col <- find_column(manifest, c("ensemble_id", "ensembleid"), "Script 12 ensemble manifest")

ensembles <- manifest[
  manifest[[status_col]] %in% c("completed", "skipped_existing_output") &
    !is.na(manifest[[mean_col]]) &
    nzchar(manifest[[mean_col]]),
  ,
  drop = FALSE
]

if (nrow(ensembles) != 16L) {
  stop("Expected 16 completed ensemble means, found ", nrow(ensembles), ".", call. = FALSE)
}
if (any(!file.exists(ensembles[[mean_col]]))) {
  print(ensembles[!file.exists(ensembles[[mean_col]]), c(id_col, mean_col), drop = FALSE])
  stop("At least one ensemble mean GeoTIFF is missing.", call. = FALSE)
}
if (anyDuplicated(ensembles[[id_col]])) stop("Duplicate ensemble IDs in manifest.", call. = FALSE)

# ---------------------------------------------------------------------
# 5. PREPARE OUTPUTS
# ---------------------------------------------------------------------

if (clean_output_root) safe_remove(output_root)
for (folder in c(current_output_dir, future_output_root, area_dir, manifest_dir)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

legend <- data.frame(
  class_code = 1:4,
  class_name = c("Unsuitable", "Low suitability", "Moderate suitability", "High suitability"),
  class_definition = c(
    paste0("Suitability < ", format(mtsps, scientific = FALSE, trim = TRUE)),
    paste0(format(mtsps, scientific = FALSE, trim = TRUE), " <= Suitability < 0.4"),
    "0.4 <= Suitability < 0.6",
    "Suitability >= 0.6"
  ),
  stringsAsFactors = FALSE
)
write.csv(legend, legend_csv, row.names = FALSE)

cat("\n======================================================================\n")
cat("MTSPS FOUR-CLASS CLASSIFICATION + AREA — Camellia japonica\n")
cat("======================================================================\n")
cat("Fixed MTSPS threshold: ", format(mtsps, scientific = TRUE, digits = 8), "\n", sep = "")
cat("Future ensemble means: ", nrow(ensembles), "\n", sep = "")
cat("Output root:\n", normalizePath(output_root, winslash = "/"), "\n", sep = "")

# ---------------------------------------------------------------------
# 6. CLASSIFY CURRENT MAPS
# ---------------------------------------------------------------------

current_inputs <- data.frame(
  map_id = c("Current_CV100_average", "Current_final_full_data"),
  map_type = c("Current CV100 average", "Current final full-data model"),
  source_model = c("Mean of 100 cross-validation replicates", "Final full-data MaxEnt model"),
  input_tif = c(current_cv100_tif, current_final_tif),
  output_tif = c(
    file.path(current_output_dir, "Camellia_japonica_current_CV100_average_MTSPS_4classes.tif"),
    file.path(current_output_dir, "Camellia_japonica_current_final_full_data_MTSPS_4classes.tif")
  ),
  stringsAsFactors = FALSE
)

status_records <- list()
area_records <- list()
start_time <- Sys.time()

for (i in seq_len(nrow(current_inputs))) {
  x <- current_inputs[i, , drop = FALSE]
  map_start <- Sys.time()

  cat("\nCurrent map: ", x$map_id, "\n", sep = "")

  if (!skip_existing_classified_maps || !file.exists(x$output_tif)) {
    continuous <- terra::rast(x$input_tif)
    check_continuous(continuous, x$map_id)
    classified <- set_categories(classify_4(continuous, mtsps), mtsps)

    terra::writeRaster(
      classified,
      x$output_tif,
      overwrite = TRUE,
      wopt = list(gdal = "COMPRESS=LZW", datatype = "INT1U")
    )
  }

  classified <- terra::rast(x$output_tif)
  area_records[[length(area_records) + 1]] <- area_table(
    classified,
    map_id = x$map_id,
    map_type = x$map_type,
    source_model = x$source_model,
    scenario = "Current",
    ssp = NA_character_,
    period_start = NA_integer_,
    period_end = NA_integer_,
    mtsps = mtsps
  )

  elapsed <- as.numeric(difftime(Sys.time(), map_start, units = "mins"))
  status_records[[length(status_records) + 1]] <- data.frame(
    map_id = x$map_id,
    map_type = x$map_type,
    scenario = "Current",
    ssp = NA_character_,
    period_start = NA_integer_,
    period_end = NA_integer_,
    status = "completed",
    input_continuous_tif = normalizePath(x$input_tif, winslash = "/", mustWork = TRUE),
    output_4class_tif = normalizePath(x$output_tif, winslash = "/", mustWork = TRUE),
    elapsed_minutes = elapsed,
    message = "",
    stringsAsFactors = FALSE
  )
  write_status(status_records)
}

# ---------------------------------------------------------------------
# 7. CLASSIFY 16 FUTURE THREE-GCM ENSEMBLES
# ---------------------------------------------------------------------

ensembles <- ensembles[
  order(ensembles[[ssp_col]], as.integer(ensembles[[start_col]])),
  ,
  drop = FALSE
]

for (i in seq_len(nrow(ensembles))) {
  x <- ensembles[i, , drop = FALSE]

  ensemble_id <- as.character(x[[id_col]])
  ssp <- as.character(x[[ssp_col]])
  p_start <- as.integer(x[[start_col]])
  p_end <- as.integer(x[[end_col]])
  input_tif <- as.character(x[[mean_col]])

  output_dir <- file.path(future_output_root, ensemble_id)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  output_tif <- file.path(
    output_dir,
    paste0(species_name, "_", ensemble_id, "_MTSPS_4classes.tif")
  )

  cat("\n[", i, "/16] Future ensemble: ", ssp, " ", p_start, "-", p_end, "\n", sep = "")
  map_start <- Sys.time()

  result <- tryCatch(
    {
      status <- "completed"

      if (skip_existing_classified_maps && file.exists(output_tif)) {
        status <- "skipped_existing_output"
      } else {
        continuous <- terra::rast(input_tif)
        check_continuous(continuous, paste(ssp, p_start, p_end))
        classified <- set_categories(classify_4(continuous, mtsps), mtsps)

        terra::writeRaster(
          classified,
          output_tif,
          overwrite = TRUE,
          wopt = list(gdal = "COMPRESS=LZW", datatype = "INT1U")
        )
      }

      classified <- terra::rast(output_tif)
      a <- area_table(
        classified,
        map_id = paste0("Future_ensemble_", ensemble_id),
        map_type = "Future three-GCM ensemble mean",
        source_model = "Mean of ACCESS-CM2, BCC-CSM2-MR, and CMCC-ESM2",
        scenario = paste0(ssp, " ", p_start, "-", p_end),
        ssp = ssp,
        period_start = p_start,
        period_end = p_end,
        mtsps = mtsps
      )

      list(status = status, areas = a, message = "")
    },
    error = function(e) list(status = "failed", areas = NULL, message = conditionMessage(e))
  )

  elapsed <- as.numeric(difftime(Sys.time(), map_start, units = "mins"))

  status_records[[length(status_records) + 1]] <- data.frame(
    map_id = paste0("Future_ensemble_", ensemble_id),
    map_type = "Future three-GCM ensemble mean",
    scenario = paste0(ssp, " ", p_start, "-", p_end),
    ssp = ssp,
    period_start = p_start,
    period_end = p_end,
    status = result$status,
    input_continuous_tif = normalizePath(input_tif, winslash = "/", mustWork = TRUE),
    output_4class_tif = if (file.exists(output_tif)) normalizePath(output_tif, winslash = "/", mustWork = TRUE) else NA_character_,
    elapsed_minutes = elapsed,
    message = result$message,
    stringsAsFactors = FALSE
  )

  if (!is.null(result$areas)) area_records[[length(area_records) + 1]] <- result$areas
  write_status(status_records)

  if (result$status == "failed") {
    cat("  Status: FAILED.\n")
    if (stop_after_first_failure) {
      stop("\nStopped after failed classification. Read:\n", status_csv, call. = FALSE)
    }
  } else {
    cat("  Status: ", result$status, " in ", format(round(elapsed, 2), nsmall = 2), " minutes.\n", sep = "")
  }
}

# ---------------------------------------------------------------------
# 8. WRITE AREA TABLES
# ---------------------------------------------------------------------

all_areas <- do.call(rbind, area_records)
all_areas <- all_areas[
  order(all_areas$map_type, all_areas$ssp, all_areas$period_start, all_areas$class_code),
  ,
  drop = FALSE
]

write.csv(all_areas, all_areas_csv, row.names = FALSE)

future_areas <- all_areas[
  all_areas$map_type == "Future three-GCM ensemble mean",
  ,
  drop = FALSE
]

current_areas <- all_areas[
  grepl("^Current", all_areas$map_type),
  ,
  drop = FALSE
]

write.csv(current_areas, current_areas_csv, row.names = FALSE)
write.csv(make_future_wide(future_areas), future_wide_csv, row.names = FALSE)

suitable_summary <- do.call(rbind, lapply(split(all_areas, all_areas$map_id), function(x) {
  total_valid <- sum(x$area_km2)
  suitable <- sum(x$area_km2[x$class_code %in% c(2, 3, 4)])

  data.frame(
    map_id = unique(x$map_id),
    map_type = unique(x$map_type),
    scenario = unique(x$scenario),
    ssp = unique(x$ssp),
    period_start = unique(x$period_start),
    period_end = unique(x$period_end),
    MTSPS_threshold = unique(x$MTSPS_threshold),
    total_valid_environment_km2 = total_valid,
    total_suitable_km2 = suitable,
    total_suitable_10_4_km2 = suitable / 10000,
    percent_suitable = if (total_valid > 0) 100 * suitable / total_valid else NA_real_,
    stringsAsFactors = FALSE
  )
}))

write.csv(suitable_summary, suitable_summary_csv, row.names = FALSE)

status_table <- do.call(rbind, status_records)
write_status(status_records)

completed <- sum(status_table$status == "completed")
skipped <- sum(status_table$status == "skipped_existing_output")
failed <- sum(status_table$status == "failed")
runtime_hours <- as.numeric(difftime(Sys.time(), start_time, units = "hours"))

write.csv(
  data.frame(
    MTSPS_threshold = mtsps,
    expected_future_ensemble_maps = 16,
    current_maps = 2,
    completed_this_run = completed,
    skipped_existing = skipped,
    failed = failed,
    runtime_hours = runtime_hours,
    stringsAsFactors = FALSE
  ),
  run_summary_csv,
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 9. REPORT
# ---------------------------------------------------------------------

cat("\n======================================================================\n")
cat("MTSPS FOUR-CLASS CLASSIFICATION + AREA COMPLETED\n")
cat("======================================================================\n")
cat("MTSPS threshold used unchanged: ", format(mtsps, scientific = TRUE, digits = 8), "\n", sep = "")
cat("Completed maps this run: ", completed, "\n", sep = "")
cat("Skipped existing maps: ", skipped, "\n", sep = "")
cat("Failed maps: ", failed, "\n", sep = "")
cat("\nAll classes and area table:\n", normalizePath(all_areas_csv, winslash = "/"), "\n", sep = "")
cat("Wide future ensemble area table:\n", normalizePath(future_wide_csv, winslash = "/"), "\n", sep = "")
cat("Total suitable area summary:\n", normalizePath(suitable_summary_csv, winslash = "/"), "\n", sep = "")
cat("Status manifest:\n", normalizePath(status_csv, winslash = "/"), "\n", sep = "")

if (failed == 0) {
  cat("\nAll 16 future ensemble maps are classified. Next: Script 14 for stable habitat, expansion, and contraction.\n")
} else {
  cat("\nResolve failed rows before Script 14.\n")
}
