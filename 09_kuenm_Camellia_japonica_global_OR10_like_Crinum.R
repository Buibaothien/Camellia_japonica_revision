# =====================================================================
# KUENM GLOBAL FC/RM CALIBRATION — Camellia japonica
# Uses the completed global pipeline (Steps 05, 07, and 08)
# =====================================================================
#
# INPUTS ALREADY CREATED BY YOUR PREVIOUS SCRIPTS
# ------------------------------------------------
# Occurrences (Step 05):
#   global strict natural-land-cover candidates after 5 km thinning
#   and one-record-per-2.5 arc-minute raster cell deduplication.
#
# Predictors (Step 07):
#   Bio2, Bio6, Bio14, Bio18, Elv
#
# Background (Step 08):
#   10,000 area-weighted global valid-environmental background cells.
#
# CALIBRATION DESIGN
# ------------------
# - Global model: no accessible-area M restriction.
# - FC: all 31 valid L/Q/P/T/H feature-class combinations.
# - RM: 0.5 to 4.0 by 0.5.
# - 31 FC x 8 RM = 248 parameter settings.
# - kuenm generates an "_all" and a "_cal" model per setting:
#   496 MaxEnt runs in total.
# - OR10 selection:
#   significant partial ROC + omission rate <= 10% + Delta AICc <= 2.
#
# IMPORTANT
# ---------
# The 10,000 background points from Step 08 are converted to a binary
# MaxEnt bias/mask grid. The grid contains all 10,000 background cells
# plus occurrence cells so MaxEnt can access environmental values at each
# presence locality. All valid cells have value 1; all other cells are NoData.
# Therefore, background selection is restricted to the Step 08 global
# area-weighted background design rather than an unrestricted global raster.
# =====================================================================

# ---------------------------------------------------------------------
# 0. INSTALL AND LOAD REQUIRED PACKAGES
# ---------------------------------------------------------------------

install_if_missing <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cloud.r-project.org")
  }
}

install_if_missing("terra")
install_if_missing("remotes")

if (!requireNamespace("kuenm", quietly = TRUE)) {
  remotes::install_github("marlonecobos/kuenm")
}

suppressPackageStartupMessages({
  library(terra)
  library(kuenm)
})

# ---------------------------------------------------------------------
# 1. USER SETTINGS
# ---------------------------------------------------------------------

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# Outputs from the completed global occurrence and predictor workflow.
occ_csv <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "05_spatial_thinning_5km",
  "05_global_strict_natural_5km_thinned_one_per_raster_cell.csv"
)

background_csv <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "08_global_background",
  "08_global_area_weighted_background_n10000.csv"
)

current_predictor_dir <- file.path(
  project_dir,
  "model_inputs",
  "Primary_Global_ClimateTopo_Current"
)

# Folder containing maxent.jar.
maxent_dir <- "D:/Prediction_tree/maxent"

# A new workspace. It is safe to keep it on E:, because the script
# writes fresh ASCII layers here rather than copying pre-existing ASC files.
workdir <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "09_kuenm_global_OR10"
)

# Set TRUE only when a complete fresh rerun is intended.
# This deletes old contents inside workdir but never changes Steps 05, 07, or 08.
clean_previous_run <- TRUE

# Expected outputs from the completed global workflow.
expected_occurrence_n <- 828L
expected_background_n <- 10000L

# Reproducible 75/25 calibration and independent-test split.
train_fraction <- 0.75
random_seed <- 20260705

# Full FC and RM search, exactly as requested.
reg_mult <- seq(0.5, 4.0, by = 0.5)
feature_classes <- "all"

# OR10 evaluation.
omission_threshold <- 10
partial_roc_random_percent <- 50
partial_roc_iterations <- 500

# MaxEnt.
max_memory_mb <- 4096

# ---------------------------------------------------------------------
# 2. HELPER FUNCTIONS
# ---------------------------------------------------------------------

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

find_column <- function(data, candidates, object_label) {
  available <- names(data)
  available_norm <- normalise_name(available)
  candidate_norm <- normalise_name(candidates)

  matched <- which(available_norm %in% candidate_norm)

  if (length(matched) == 0) {
    stop(
      "\nCould not identify a required column in ", object_label, ".\n",
      "Expected one of: ", paste(candidates, collapse = ", "),
      "\nAvailable columns: ", paste(available, collapse = ", "),
      call. = FALSE
    )
  }

  available[matched[1]]
}

safe_remove <- function(path) {
  if (file.exists(path) || dir.exists(path)) {
    unlink(path, recursive = TRUE, force = TRUE)
  }
}

check_grid_alignment <- function(reference_raster, candidate_raster, label) {
  aligned <- terra::compareGeom(
    reference_raster,
    candidate_raster,
    stopOnError = FALSE,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE
  )

  if (!aligned) {
    stop(
      "\nGrid mismatch detected for: ", label,
      "\nAll five selected predictors must have identical CRS, extent, ",
      "resolution, rows, columns, and cell alignment.",
      call. = FALSE
    )
  }
}

require_valid_values <- function(raster_stack, points, label) {
  extracted <- terra::extract(raster_stack, points, ID = FALSE)

  # Convert all extracted predictor columns to a matrix.
  extracted_matrix <- as.matrix(extracted)
  valid_rows <- apply(
    extracted_matrix,
    1,
    function(x) all(is.finite(x) & x != -9999)
  )

  if (!all(valid_rows)) {
    stop(
      "\n", sum(!valid_rows), " record(s) in ", label,
      " fall on NoData cells for one or more selected predictors.\n",
      "Do not continue until the occurrence/background input is reconciled ",
      "with the selected current predictor rasters.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

# ---------------------------------------------------------------------
# 3. VERIFY DIRECT INPUTS
# ---------------------------------------------------------------------

required_paths <- c(
  occurrence_csv = occ_csv,
  background_csv = background_csv,
  current_predictor_directory = current_predictor_dir,
  maxent_directory = maxent_dir,
  maxent_jar = file.path(maxent_dir, "maxent.jar")
)

missing_paths <- required_paths[
  !file.exists(required_paths) & !dir.exists(required_paths)
]

if (length(missing_paths) > 0) {
  stop(
    "\nMissing required input path(s):\n",
    paste(names(missing_paths), ": ", missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

if (Sys.which("java") == "") {
  warning(
    "Java was not found on the system PATH. MaxEnt cannot run until a compatible ",
    "Java installation is available on PATH."
  )
}

# ---------------------------------------------------------------------
# 4. DEFINE THE FINAL PRIMARY PREDICTOR SET
# ---------------------------------------------------------------------

selected_predictors <- c("Bio2", "Bio6", "Bio14", "Bio18", "Elv")

predictor_files <- stats::setNames(
  file.path(current_predictor_dir, paste0(selected_predictors, ".tif")),
  selected_predictors
)

missing_predictors <- predictor_files[!file.exists(predictor_files)]

if (length(missing_predictors) > 0) {
  stop(
    "\nThe following final primary predictor raster(s) are missing:\n",
    paste(missing_predictors, collapse = "\n"),
    "\n\nRun Step 07 again or check the file extension/capitalization.",
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 5. CREATE A CLEAN KUENM WORKSPACE
# ---------------------------------------------------------------------

dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
setwd(workdir)

variables_root <- "01_global_variables"
set1_dir <- file.path(variables_root, "set1_primary5")
occ_dir <- "02_occurrences"
background_dir <- "03_background_mask"
candidate_dir <- "04_candidate_models"
evaluation_dir <- "05_calibration_results"
log_dir <- "06_logs"

if (clean_previous_run) {
  for (item in c(
    variables_root,
    occ_dir,
    background_dir,
    candidate_dir,
    evaluation_dir,
    log_dir
  )) {
    safe_remove(item)
  }
}

dir.create(set1_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(occ_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(background_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(evaluation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("\n======================================================================\n")
cat("KUENM GLOBAL FC/RM CALIBRATION — Camellia japonica\n")
cat("======================================================================\n")
cat("Occurrence input:\n", normalizePath(occ_csv, winslash = "/"), "\n", sep = "")
cat("Background input:\n", normalizePath(background_csv, winslash = "/"), "\n", sep = "")
cat("Predictor directory:\n", normalizePath(current_predictor_dir, winslash = "/"), "\n", sep = "")
cat("MaxEnt JAR:\n", normalizePath(file.path(maxent_dir, "maxent.jar"), winslash = "/"), "\n", sep = "")
cat("Workspace:\n", normalizePath(workdir, winslash = "/"), "\n\n", sep = "")

# ---------------------------------------------------------------------
# 6. READ AND VERIFY FINAL OCCURRENCE DATA
# ---------------------------------------------------------------------

occ_raw <- read.csv(occ_csv, check.names = FALSE, stringsAsFactors = FALSE)

occ_lon_column <- find_column(
  occ_raw,
  c("decimalLongitude", "longitude", "lon", "x"),
  "final occurrence CSV"
)

occ_lat_column <- find_column(
  occ_raw,
  c("decimalLatitude", "latitude", "lat", "y"),
  "final occurrence CSV"
)

occ_final <- data.frame(
  species = species_name,
  longitude = suppressWarnings(as.numeric(occ_raw[[occ_lon_column]])),
  latitude = suppressWarnings(as.numeric(occ_raw[[occ_lat_column]])),
  stringsAsFactors = FALSE
)

occ_final <- occ_final[
  is.finite(occ_final$longitude) &
    is.finite(occ_final$latitude) &
    occ_final$longitude >= -180 &
    occ_final$longitude <= 180 &
    occ_final$latitude >= -90 &
    occ_final$latitude <= 90,
  ,
  drop = FALSE
]

occ_final <- unique(occ_final)

if (nrow(occ_final) != expected_occurrence_n) {
  stop(
    "\nThe Step 05 occurrence input contains ", nrow(occ_final),
    " valid unique coordinates, but the completed current workflow is expected ",
    "to contain ", expected_occurrence_n, ".\n",
    "Check that occ_csv points to:\n",
    occ_csv,
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 7. READ AND VERIFY THE STEP 08 GLOBAL BACKGROUND SAMPLE
# ---------------------------------------------------------------------

background_raw <- read.csv(
  background_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

bg_lon_column <- find_column(
  background_raw,
  c("decimalLongitude", "longitude", "lon", "x"),
  "Step 08 background CSV"
)

bg_lat_column <- find_column(
  background_raw,
  c("decimalLatitude", "latitude", "lat", "y"),
  "Step 08 background CSV"
)

background_final <- data.frame(
  longitude = suppressWarnings(as.numeric(background_raw[[bg_lon_column]])),
  latitude = suppressWarnings(as.numeric(background_raw[[bg_lat_column]])),
  stringsAsFactors = FALSE
)

background_final <- background_final[
  is.finite(background_final$longitude) &
    is.finite(background_final$latitude) &
    background_final$longitude >= -180 &
    background_final$longitude <= 180 &
    background_final$latitude >= -90 &
    background_final$latitude <= 90,
  ,
  drop = FALSE
]

background_final <- unique(background_final)

if (nrow(background_final) != expected_background_n) {
  stop(
    "\nThe Step 08 background CSV contains ", nrow(background_final),
    " valid unique coordinates, not ", expected_background_n, ".\n",
    "Check that background_csv points to the final output of Step 08:\n",
    background_csv,
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 8. PREPARE GLOBAL ASCII PREDICTORS FOR KUENM/MAXENT
# ---------------------------------------------------------------------

reference_raster <- terra::rast(predictor_files["Bio2"])

reference_extent <- terra::ext(reference_raster)

if (
  terra::xmin(reference_extent) < -181 ||
  terra::xmax(reference_extent) > 181 ||
  terra::ymin(reference_extent) < -91 ||
  terra::ymax(reference_extent) > 91
) {
  stop(
    "\nCurrent predictors are not in a geographic longitude-latitude grid.\n",
    "This script expects the 2.5 arc-minute global WGS84 layers produced by Step 07.",
    call. = FALSE
  )
}

cat("Converting the five final current GeoTIFFs to ASCII grids for MaxEnt ...\n")

ascii_files <- character(length(predictor_files))
names(ascii_files) <- names(predictor_files)

for (i in seq_along(predictor_files)) {
  predictor_name <- names(predictor_files)[i]
  source_file <- predictor_files[i]
  output_file <- file.path(set1_dir, paste0(predictor_name, ".asc"))

  cat("  [", i, "/", length(predictor_files), "] ", basename(source_file), "\n", sep = "")

  current_raster <- terra::rast(source_file)

  check_grid_alignment(
    reference_raster,
    current_raster,
    basename(source_file)
  )

  terra::writeRaster(
    current_raster,
    filename = output_file,
    filetype = "AAIGrid",
    NAflag = -9999,
    overwrite = TRUE
  )

  if (!file.exists(output_file) || file.info(output_file)$size <= 0) {
    stop(
      "\nFailed to create a valid ASCII predictor file:\n",
      output_file,
      call. = FALSE
    )
  }

  ascii_files[i] <- output_file
}

write.csv(
  data.frame(
    predictor = names(predictor_files),
    source_geotiff = unname(predictor_files),
    kuenm_ascii = unname(ascii_files),
    stringsAsFactors = FALSE
  ),
  file.path(log_dir, "predictor_manifest.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 9. CHECK OCCURRENCES AND BACKGROUND AGAINST ALL FIVE PREDICTORS
# ---------------------------------------------------------------------

predictor_stack <- terra::rast(unname(predictor_files))

occ_points <- terra::vect(
  occ_final,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

background_points <- terra::vect(
  background_final,
  geom = c("longitude", "latitude"),
  crs = "EPSG:4326"
)

if (!terra::same.crs(occ_points, reference_raster) ||
    !terra::same.crs(background_points, reference_raster)) {
  stop(
    "\nThe current predictor CRS is not EPSG:4326 despite being a geographic grid.\n",
    "Revise the coordinate transformation before continuing.",
    call. = FALSE
  )
}

cat("Checking 828 occurrence records against all five selected predictors ...\n")
require_valid_values(predictor_stack, occ_points, "final occurrence input")

cat("Checking 10,000 global background records against all five selected predictors ...\n")
require_valid_values(predictor_stack, background_points, "Step 08 background input")

# ---------------------------------------------------------------------
# 10. CREATE A BINARY BACKGROUND BIAS/MASK GRID
# ---------------------------------------------------------------------
# The global Step 08 background cells are retained as the background domain.
# Occurrence cells are additionally included so the bias grid is defined at
# all presence localities. All valid cells have the same positive value (= 1).

background_cells <- terra::cellFromXY(
  reference_raster,
  as.matrix(background_final[, c("longitude", "latitude")])
)

occurrence_cells <- terra::cellFromXY(
  reference_raster,
  as.matrix(occ_final[, c("longitude", "latitude")])
)

background_cells <- background_cells[is.finite(background_cells)]
occurrence_cells <- occurrence_cells[is.finite(occurrence_cells)]

allowed_bias_cells <- sort(unique(c(background_cells, occurrence_cells)))

if (length(allowed_bias_cells) < expected_background_n) {
  stop(
    "\nThe bias/mask grid has fewer than ", expected_background_n,
    " valid cells after mapping Step 08 coordinates to the reference raster.\n",
    "Inspect the background CSV and the Bio2 reference raster alignment.",
    call. = FALSE
  )
}

# Create an empty raster with NoData everywhere, then assign 1 only to
# the Step 08 background cells and occurrence cells.
bias_raster <- terra::rast(
  nrows = terra::nrow(reference_raster),
  ncols = terra::ncol(reference_raster),
  ext = terra::ext(reference_raster),
  crs = terra::crs(reference_raster)
)

bias_raster[allowed_bias_cells] <- 1

bias_asc <- file.path(
  background_dir,
  "Camellia_japonica_global_background_mask.asc"
)

terra::writeRaster(
  bias_raster,
  filename = bias_asc,
  filetype = "AAIGrid",
  NAflag = -9999,
  overwrite = TRUE
)

if (!file.exists(bias_asc) || file.info(bias_asc)$size <= 0) {
  stop("\nFailed to create the global background bias/mask grid.", call. = FALSE)
}

write.csv(
  data.frame(
    source_background_points = nrow(background_final),
    unique_background_cells = length(unique(background_cells)),
    unique_occurrence_cells = length(unique(occurrence_cells)),
    total_positive_bias_cells = length(allowed_bias_cells),
    bias_value_at_positive_cells = 1,
    stringsAsFactors = FALSE
  ),
  file.path(log_dir, "background_mask_summary.csv"),
  row.names = FALSE
)

# maximumbackground is set larger than the valid bias-mask domain so MaxEnt
# retains all cells in the specified background domain rather than subsampling.
maximum_background <- length(allowed_bias_cells) + 1L

bias_path_windows <- gsub(
  "/",
  "\\\\",
  normalizePath(bias_asc, winslash = "/", mustWork = TRUE),
  fixed = TRUE
)

maxent_args <- paste0(
  'biasfile="', bias_path_windows,
  '" biastype=3 maximumbackground=', maximum_background,
  ' cache=true'
)

cat("\nBackground-mask summary:\n")
cat("  Step 08 unique background cells: ", length(unique(background_cells)), "\n", sep = "")
cat("  Occurrence cells added to mask: ", length(unique(occurrence_cells)), "\n", sep = "")
cat("  Total positive bias-mask cells: ", length(allowed_bias_cells), "\n", sep = "")

# ---------------------------------------------------------------------
# 11. CREATE THE REPRODUCIBLE 75/25 TRAINING AND TEST SPLIT
# ---------------------------------------------------------------------

set.seed(random_seed)

n_occ <- nrow(occ_final)
n_train <- floor(train_fraction * n_occ)

training_index <- sample(seq_len(n_occ), size = n_train, replace = FALSE)

occ_train <- occ_final[training_index, , drop = FALSE]
occ_test <- occ_final[-training_index, , drop = FALSE]

occ_joint_file <- file.path(occ_dir, "Camellia_japonica_joint.csv")
occ_train_file <- file.path(occ_dir, "Camellia_japonica_train.csv")
occ_test_file <- file.path(occ_dir, "Camellia_japonica_test.csv")

write.csv(occ_final, occ_joint_file, row.names = FALSE)
write.csv(occ_train, occ_train_file, row.names = FALSE)
write.csv(occ_test, occ_test_file, row.names = FALSE)

cat("\nOccurrence partition:\n")
cat("  Joint records: ", nrow(occ_final), "\n", sep = "")
cat("  Training records: ", nrow(occ_train), "\n", sep = "")
cat("  Testing records: ", nrow(occ_test), "\n", sep = "")

# ---------------------------------------------------------------------
# 12. RUN FULL KUENM CALIBRATION: 31 FC × 8 RM
# ---------------------------------------------------------------------

batch_name <- "Camellia_japonica_global_OR10_allFC"
n_feature_classes <- 31L
n_rm <- length(reg_mult)
n_parameter_settings <- n_feature_classes * n_rm
n_expected_maxent_runs <- 2L * n_parameter_settings

cat(
  "\nCandidate-model configuration:\n",
  "  Predictor set: Bio2, Bio6, Bio14, Bio18, Elv\n",
  "  Background: Step 08 global area-weighted 10,000-cell sample\n",
  "  Feature classes: ALL 31 valid L/Q/P/T/H combinations\n",
  "  RM values: ", paste(reg_mult, collapse = ", "), "\n",
  "  Parameter settings: 31 x ", n_rm, " = ", n_parameter_settings, "\n",
  "  Expected MaxEnt runs: ", n_expected_maxent_runs, "\n",
  "  Selection: partial ROC + omission rate <= ", omission_threshold,
  "% + Delta AICc <= 2\n\n",
  sep = ""
)

kuenm_cal(
  occ.joint = occ_joint_file,
  occ.tra = occ_train_file,
  M.var.dir = variables_root,
  batch = batch_name,
  out.dir = candidate_dir,
  max.memory = max_memory_mb,
  reg.mult = reg_mult,
  f.clas = feature_classes,
  args = maxent_args,
  maxent.path = maxent_dir,
  wait = TRUE,
  run = TRUE
)

# ---------------------------------------------------------------------
# 13. VERIFY THAT ALL MAXENT CANDIDATE RUNS FINISHED
# ---------------------------------------------------------------------

result_files <- list.files(
  candidate_dir,
  pattern = "^maxentResults\\.csv$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

lambda_files <- list.files(
  candidate_dir,
  pattern = "\\.lambdas$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)

completion_check <- data.frame(
  expected_parameter_settings = n_parameter_settings,
  expected_total_MaxEnt_runs = n_expected_maxent_runs,
  maxentResults_csv_found = length(result_files),
  lambdas_found = length(lambda_files),
  stringsAsFactors = FALSE
)

write.csv(
  completion_check,
  file.path(log_dir, "candidate_model_completion_check.csv"),
  row.names = FALSE
)

cat("\nCandidate-model completion check:\n")
print(completion_check)

if (
  length(result_files) < n_expected_maxent_runs ||
  length(lambda_files) < n_expected_maxent_runs
) {
  stop(
    "\nAt least one MaxEnt candidate run failed or remained incomplete.\n",
    "kuenm_ceval was NOT started.\n",
    "Read the terminal output and inspect:\n",
    normalizePath(
      file.path(log_dir, "candidate_model_completion_check.csv"),
      winslash = "/",
      mustWork = TRUE
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 14. EVALUATE: PARTIAL ROC + OR10 + AICc
# ---------------------------------------------------------------------

cat("\nAll candidate outputs are present. Starting kuenm_ceval.\n")

evaluation_results <- kuenm_ceval(
  path = candidate_dir,
  occ.joint = occ_joint_file,
  occ.tra = occ_train_file,
  occ.test = occ_test_file,
  batch = batch_name,
  out.eval = evaluation_dir,
  threshold = omission_threshold,
  rand.percent = partial_roc_random_percent,
  iterations = partial_roc_iterations,
  kept = TRUE,
  selection = "OR_AICc",
  parallel.proc = FALSE
)

saveRDS(
  evaluation_results,
  file.path(evaluation_dir, "Camellia_japonica_global_OR10_kuenm_results.rds")
)

writeLines(
  list.files(evaluation_dir, recursive = TRUE, full.names = TRUE),
  file.path(log_dir, "evaluation_output_files.txt")
)

cat("\n======================================================================\n")
cat("KUENM CALIBRATION AND EVALUATION FINISHED — Camellia japonica\n")
cat("======================================================================\n")
cat("Candidate models:\n", normalizePath(candidate_dir, winslash = "/"), "\n", sep = "")
cat("Evaluation results:\n", normalizePath(evaluation_dir, winslash = "/"), "\n", sep = "")
cat("Open selected_models.csv and calibration_results.csv to identify the chosen FC/RM.\n")
