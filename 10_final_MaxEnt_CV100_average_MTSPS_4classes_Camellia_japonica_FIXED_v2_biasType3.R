# =====================================================================
# FINAL MAXENT CV100 + 4-CLASS AREA ANALYSIS — FIXED VERSION
# Camellia japonica
# =====================================================================
#
# FIXES APPLIED
# -------------
# 1) Uses system2() to write logs. The prior version incorrectly passed
#    "> log_file 2>&1" into MaxEnt as initialization flags.
#
# 2) Rebuilds the bias raster on the EXACT Bio2.ASC grid used by MaxEnt.
#    The new bias raster has:
#       - value 1.0 at the 10,000 Step-08 background cells;
#       - value 1.0 at occurrence cells;
#       - value 1e-6 at all other valid environmental cells;
#       - NoData outside the common valid environmental domain.
#    Thus, every valid occurrence has a positive bias value while the
#    Step-08 background cells receive overwhelmingly higher sampling weight.
#
# 3) Explicitly sets biasType=3 because MaxEnt otherwise assigns the
#    bias layer type Unknown (0), which causes an Unknown feature type error.
#
# MODEL DESIGN
# ------------
# Tuned settings from kuenm:
#   RM = 1
#   FC = lph
#
# Outputs:
#   A. Mean prediction across 100 MaxEnt cross-validation replicates.
#   B. Standard deviation across CV replicates, when MaxEnt writes it.
#   C. Final full-data model for the one MTSPS threshold.
#   D. Four suitability classes and geographic areas for both maps.
# =====================================================================

# ---------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------

install_if_missing <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    install.packages(package_name, repos = "https://cloud.r-project.org")
  }
}

install_if_missing("terra")
suppressPackageStartupMessages(library(terra))

# ---------------------------------------------------------------------
# 1. SETTINGS
# ---------------------------------------------------------------------

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

kuenm_workdir <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "09_kuenm_global_OR10"
)

maxent_dir <- "D:/Prediction_tree/maxent"
maxent_jar <- file.path(maxent_dir, "maxent.jar")

output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

# TRUE deletes only this script's own output folder.
clean_previous_run <- TRUE

# Selected kuenm settings.
regularization_multiplier <- 1
feature_class_label <- "lph"
selected_predictors <- c("Bio2", "Bio6", "Bio14", "Bio18", "Elv")

# Cross-validation.
n_cross_validation_replicates <- 100L
max_memory_mb <- 4096L

# A strictly positive low weight for all valid environmental cells not in the
# Step-08 10,000-point background. Do not use zero in a MaxEnt bias raster.
bias_floor <- 1e-6

# ---------------------------------------------------------------------
# 2. INPUT PATHS
# ---------------------------------------------------------------------

occurrence_csv <- file.path(
  kuenm_workdir,
  "02_occurrences",
  "Camellia_japonica_joint.csv"
)

background_csv <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "08_global_background",
  "08_global_area_weighted_background_n10000.csv"
)

environmental_asc_dir <- file.path(
  kuenm_workdir,
  "01_global_variables",
  "set1_primary5"
)

current_predictor_dir <- file.path(
  project_dir,
  "model_inputs",
  "Primary_Global_ClimateTopo_Current"
)

selected_models_csv <- file.path(
  kuenm_workdir,
  "05_calibration_results",
  "selected_models.csv"
)

expected_selected_model <- "M_1_F_lph_set1_primary5"

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

find_column <- function(data, candidates, object_label, required = TRUE) {
  original_names <- names(data)
  normal_names <- normalise_name(original_names)
  normal_candidates <- normalise_name(candidates)

  direct <- match(normal_candidates, normal_names, nomatch = 0)
  direct <- direct[direct > 0]

  if (length(direct) > 0) {
    return(original_names[direct[1]])
  }

  for (candidate in normal_candidates) {
    partial <- grep(candidate, normal_names, fixed = TRUE)
    if (length(partial) > 0) {
      return(original_names[partial[1]])
    }
  }

  if (required) {
    stop(
      "\nCould not find a required column in ", object_label, ".\n",
      "Available columns: ", paste(original_names, collapse = ", "),
      call. = FALSE
    )
  }

  NA_character_
}

read_xy_csv <- function(path, object_label, add_species = FALSE) {
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  longitude_col <- find_column(
    data,
    c("longitude", "decimalLongitude", "decimalLon", "lon", "x"),
    object_label
  )

  latitude_col <- find_column(
    data,
    c("latitude", "decimalLatitude", "decimalLat", "lat", "y"),
    object_label
  )

  output <- data.frame(
    longitude = suppressWarnings(as.numeric(data[[longitude_col]])),
    latitude = suppressWarnings(as.numeric(data[[latitude_col]])),
    stringsAsFactors = FALSE
  )

  output <- output[
    is.finite(output$longitude) &
      is.finite(output$latitude) &
      output$longitude >= -180 &
      output$longitude <= 180 &
      output$latitude >= -90 &
      output$latitude <= 90,
    ,
    drop = FALSE
  ]

  output <- unique(output)

  if (nrow(output) == 0) {
    stop("No valid coordinates remain in ", object_label, ".", call. = FALSE)
  }

  if (add_species) {
    output <- data.frame(
      species = rep(species_name, nrow(output)),
      output,
      stringsAsFactors = FALSE
    )
  }

  output
}

xy_to_vect <- function(data) {
  terra::vect(
    data,
    geom = c("longitude", "latitude"),
    crs = "EPSG:4326"
  )
}

quoted_parameter_path <- function(parameter_name, path, must_work = TRUE) {
  normalized_path <- normalizePath(
    path,
    winslash = "\\",
    mustWork = must_work
  )

  paste0(parameter_name, "=", shQuote(normalized_path))
}

assert_grid_alignment <- function(reference, candidate, label) {
  same <- terra::compareGeom(
    reference,
    candidate,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )

  if (!same) {
    stop(
      "\nRaster-grid mismatch: ", label,
      "\nAll final predictors must have identical geometry.",
      call. = FALSE
    )
  }
}

build_weighted_bias_raster <- function(
  env_asc_paths,
  background_xy,
  occurrence_xy,
  output_asc,
  output_log_csv
) {
  # IMPORTANT: template is Bio2.ASC, which is exactly what MaxEnt reads.
  template <- terra::rast(env_asc_paths["Bio2"])

  # Create a low, strictly positive floor over the common valid domain of
  # all five environmental variables.
  valid_floor <- terra::ifel(!is.na(template), bias_floor, NA)

  for (predictor_name in names(env_asc_paths)[-1]) {
    predictor <- terra::rast(env_asc_paths[predictor_name])

    assert_grid_alignment(
      template,
      predictor,
      paste0(predictor_name, ".asc")
    )

    valid_floor <- terra::mask(valid_floor, predictor)
  }

  # Selected 10,000 background cells and all occurrence cells receive weight 1.
  # This keeps the Step-08 design dominant while avoiding missing bias values
  # at sample locations.
  weighted_points <- rbind(
    data.frame(
      longitude = background_xy$longitude,
      latitude = background_xy$latitude,
      group = "background"
    ),
    data.frame(
      longitude = occurrence_xy$longitude,
      latitude = occurrence_xy$latitude,
      group = "occurrence"
    )
  )

  high_weight <- terra::rasterize(
    xy_to_vect(weighted_points),
    template,
    field = 1,
    background = NA
  )

  bias_raster <- terra::cover(high_weight, valid_floor)
  names(bias_raster) <- "global_background_weight"

  dir.create(dirname(output_asc), recursive = TRUE, showWarnings = FALSE)

  terra::writeRaster(
    bias_raster,
    output_asc,
    filetype = "AAIGrid",
    NAflag = -9999,
    overwrite = TRUE
  )

  # Re-read exactly as MaxEnt will read it, then validate every occurrence.
  written_bias <- terra::rast(output_asc)

  background_values <- terra::extract(
    written_bias,
    xy_to_vect(background_xy),
    ID = FALSE
  )[[1]]

  occurrence_values <- terra::extract(
    written_bias,
    xy_to_vect(occurrence_xy),
    ID = FALSE
  )[[1]]

  if (
    any(!is.finite(occurrence_values)) ||
    any(occurrence_values <= 0)
  ) {
    stop(
      "\nThe rebuilt bias raster still has missing/non-positive values at ",
      sum(!is.finite(occurrence_values) | occurrence_values <= 0),
      " occurrence location(s).",
      call. = FALSE
    )
  }

  if (
    any(!is.finite(background_values)) ||
    any(background_values <= 0)
  ) {
    stop(
      "\nThe rebuilt bias raster has missing/non-positive values at ",
      sum(!is.finite(background_values) | background_values <= 0),
      " Step-08 background location(s).",
      call. = FALSE
    )
  }

  all_values <- terra::global(
    written_bias,
    fun = c("min", "max"),
    na.rm = TRUE
  )

  log_table <- data.frame(
    item = c(
      "bias_floor",
      "n_background_input_points",
      "n_occurrence_input_points",
      "n_background_cells_weight_1",
      "n_occurrence_cells_weight_1",
      "bias_min_positive_value",
      "bias_max_value"
    ),
    value = c(
      bias_floor,
      nrow(background_xy),
      nrow(occurrence_xy),
      sum(background_values >= 0.999999),
      sum(occurrence_values >= 0.999999),
      all_values[1, "min"],
      all_values[1, "max"]
    ),
    stringsAsFactors = FALSE
  )

  write.csv(log_table, output_log_csv, row.names = FALSE)

  written_bias
}

run_maxent <- function(
  samples_csv,
  output_dir,
  environmental_dir,
  bias_asc,
  log_file,
  model_type = c("cv100", "full_data")
) {
  model_type <- match.arg(model_type)

  safe_remove(output_dir)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  java_exe <- Sys.which("java")

  if (identical(java_exe, "")) {
    stop("Java is not available on PATH.", call. = FALSE)
  }

  base_parameters <- c(
    quoted_parameter_path("samplesfile", samples_csv),
    quoted_parameter_path("environmentallayers", environmental_dir),
    quoted_parameter_path("outputdirectory", output_dir, must_work = FALSE),
    quoted_parameter_path("biasfile", bias_asc),

    # MaxEnt's hidden biasType default is 0 (= Unknown). Set 3 explicitly:
    # Layer type 3 is the internal "Bias" layer type required for biasfile.
    "biasType=3",

    # Exactly 10,000 MaxEnt background points are sampled from the
    # weighted bias surface.
    "maximumbackground=10000",

    # Tuned MaxEnt settings.
    "linear=true",
    "quadratic=false",
    "product=true",
    "threshold=false",
    "hinge=true",
    "autofeature=false",
    paste0("betamultiplier=", regularization_multiplier),

    # Do not alter the externally specified weighted background scheme.
    "addsamplestobackground=false",
    "addallsamplestobackground=false",

    "outputformat=cloglog",
    "outputfiletype=asc",
    "removeduplicates=false",
    "randomseed=false",
    "askoverwrite=false",
    "visible=false",
    "cache=true",
    "appendtoresultsfile=false"
  )

  if (model_type == "cv100") {
    parameters <- c(
      base_parameters,
      paste0("replicates=", n_cross_validation_replicates),
      "replicatetype=Crossvalidate",
      "randomtestpoints=0",

      # Summary maps remain enabled; individual replicate rasters are omitted.
      "outputgrids=false",
      "plots=true",
      "responsecurves=false",
      "jackknife=false",
      "writeplotdata=false",
      "writebackgroundpredictions=false"
    )
  } else {
    parameters <- c(
      base_parameters,
      "replicates=1",
      "randomtestpoints=0",
      "outputgrids=true",
      "plots=true",
      "responsecurves=true",
      "jackknife=true",
      "writeplotdata=true",
      "writebackgroundpredictions=true",
      "writeclampgrid=true",
      "writemess=true"
    )
  }

  args <- c(
    paste0("-mx", max_memory_mb, "m"),
    "-jar",
    shQuote(normalizePath(maxent_jar, winslash = "\\", mustWork = TRUE)),
    parameters
  )

  command_preview <- paste(
    shQuote(java_exe),
    paste(args, collapse = " ")
  )

  writeLines(command_preview, file.path(output_dir, "maxent_command.txt"))

  # CRITICAL: stdout/stderr are handled by system2(), not appended as
  # command-line arguments. This prevents ">" and "2>&1" being read by MaxEnt.
  exit_status <- system2(
    command = java_exe,
    args = args,
    stdout = log_file,
    stderr = log_file,
    wait = TRUE
  )

  result_csv <- file.path(output_dir, "maxentResults.csv")

  if (exit_status != 0 || !file.exists(result_csv)) {
    stop(
      "\nMaxEnt failed for ", model_type, ".\n",
      "Read the log file:\n", log_file,
      call. = FALSE
    )
  }

  result_csv
}

find_prediction_grid <- function(output_dir, mode = c("cv_mean", "cv_sd", "final")) {
  mode <- match.arg(mode)

  asc_files <- list.files(
    output_dir,
    pattern = "\\.asc$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(asc_files) == 0) {
    stop("No ASC grids were written by MaxEnt in:\n", output_dir, call. = FALSE)
  }

  names_only <- basename(asc_files)

  if (mode == "cv_mean") {
    candidate <- asc_files[
      grepl("(_avg|_average)\\.asc$", names_only, ignore.case = TRUE) &
        !grepl("std|stdev|clamp|mess", names_only, ignore.case = TRUE)
    ]

    if (length(candidate) == 1) {
      return(candidate)
    }
  }

  if (mode == "cv_sd") {
    candidate <- asc_files[
      grepl("(stddev|stdev|_std)", names_only, ignore.case = TRUE) &
        !grepl("clamp|mess", names_only, ignore.case = TRUE)
    ]

    if (length(candidate) == 1) {
      return(candidate)
    }

    if (length(candidate) == 0) {
      return(NA_character_)
    }
  }

  if (mode == "final") {
    expected <- file.path(output_dir, paste0(species_name, ".asc"))

    if (file.exists(expected)) {
      return(expected)
    }

    candidate <- asc_files[
      !grepl("avg|average|std|stdev|clamp|mess", names_only, ignore.case = TRUE)
    ]

    if (length(candidate) == 1) {
      return(candidate)
    }
  }

  stop(
    "\nCould not identify the requested MaxEnt prediction grid: ", mode,
    "\nAvailable ASC files:\n",
    paste(names_only, collapse = "\n"),
    call. = FALSE
  )
}

restore_crs <- function(prediction_asc, reference_tif) {
  prediction <- terra::rast(prediction_asc)

  if (is.na(terra::crs(prediction)) || terra::crs(prediction) == "") {
    terra::crs(prediction) <- terra::crs(reference_tif)
  }

  prediction
}

extract_scores <- function(prediction, xy_data, label) {
  scores <- terra::extract(
    prediction,
    xy_to_vect(xy_data),
    ID = FALSE
  )[[1]]

  scores <- suppressWarnings(as.numeric(scores))

  if (any(!is.finite(scores))) {
    stop(
      sum(!is.finite(scores)),
      " non-finite prediction(s) for ", label, ".",
      call. = FALSE
    )
  }

  scores
}

mtsps_threshold <- function(presence_scores, background_scores) {
  presence_scores <- sort(as.numeric(presence_scores))
  background_scores <- sort(as.numeric(background_scores))

  thresholds <- sort(unique(c(presence_scores, background_scores)))

  tolerance <- pmax(
    .Machine$double.eps * pmax(1, abs(thresholds)) * 16,
    1e-15
  )

  presence_below <- findInterval(thresholds - tolerance, presence_scores)
  background_below <- findInterval(thresholds - tolerance, background_scores)

  sensitivity <- (length(presence_scores) - presence_below) / length(presence_scores)
  specificity <- background_below / length(background_scores)

  table <- data.frame(
    threshold = thresholds,
    sensitivity = sensitivity,
    specificity = specificity,
    sensitivity_plus_specificity = sensitivity + specificity,
    TSS = sensitivity + specificity - 1,
    stringsAsFactors = FALSE
  )

  best <- max(table$sensitivity_plus_specificity)

  best_rows <- table[
    abs(table$sensitivity_plus_specificity - best) < 1e-12,
    ,
    drop = FALSE
  ]

  list(
    selected = best_rows[which.min(best_rows$threshold), , drop = FALSE],
    curve = table
  )
}

binary_metrics <- function(presence_scores, background_scores, threshold) {
  tp <- sum(presence_scores >= threshold)
  fn <- sum(presence_scores < threshold)
  fp <- sum(background_scores >= threshold)
  tn <- sum(background_scores < threshold)

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)
  tss <- sensitivity + specificity - 1
  accuracy <- (tp + tn) / (tp + fn + fp + tn)

  total <- tp + fn + fp + tn
  chance_accuracy <- (
    ((tp + fp) * (tp + fn)) +
      ((tn + fn) * (tn + fp))
  ) / total^2

  kappa <- if (abs(1 - chance_accuracy) < .Machine$double.eps) {
    NA_real_
  } else {
    (accuracy - chance_accuracy) / (1 - chance_accuracy)
  }

  data.frame(
    threshold = threshold,
    sensitivity = sensitivity,
    specificity = specificity,
    TSS = tss,
    accuracy = accuracy,
    Kappa = kappa,
    TP = tp,
    FN = fn,
    FP = fp,
    TN = tn,
    stringsAsFactors = FALSE
  )
}

extract_cv_summary <- function(maxent_results_csv, expected_replicates) {
  result <- read.csv(
    maxent_results_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  species_col <- find_column(
    result,
    c("Species", "species"),
    "CV100 maxentResults.csv"
  )

  ids <- as.character(result[[species_col]])

  replicate_rows <- grepl(
    paste0("^", species_name, "_[0-9]+$"),
    ids,
    ignore.case = TRUE
  )

  replicate_data <- result[replicate_rows, , drop = FALSE]

  if (nrow(replicate_data) != expected_replicates) {
    stop(
      "\nExpected ", expected_replicates, " cross-validation replicate rows, found ",
      nrow(replicate_data), ".\nInspect:\n", maxent_results_csv,
      call. = FALSE
    )
  }

  numeric_columns <- character(0)

  for (name in names(replicate_data)) {
    if (name == species_col) next

    converted <- suppressWarnings(as.numeric(replicate_data[[name]]))

    if (sum(is.finite(converted)) > 0) {
      replicate_data[[name]] <- converted
      numeric_columns <- c(numeric_columns, name)
    }
  }

  summary <- do.call(
    rbind,
    lapply(numeric_columns, function(name) {
      values <- replicate_data[[name]]
      values <- values[is.finite(values)]

      data.frame(
        metric = name,
        n_replicates = length(values),
        mean = mean(values),
        sd = stats::sd(values),
        min = min(values),
        max = max(values),
        stringsAsFactors = FALSE
      )
    })
  )

  list(replicates = replicate_data, summary = summary)
}

classify_four_classes <- function(suitability, mtsps) {
  # 1: unsuitable; 2: low; 3: moderate; 4: high
  terra::ifel(
    suitability < mtsps,
    1,
    terra::ifel(
      suitability < 0.4,
      2,
      terra::ifel(suitability < 0.6, 3, 4)
    )
  )
}

calculate_class_areas <- function(class_raster, map_name, mtsps) {
  class_area <- terra::expanse(
    class_raster,
    unit = "km",
    transform = TRUE,
    byValue = TRUE
  )

  value_col <- find_column(class_area, c("value"), "class area table")
  area_col <- find_column(class_area, c("area"), "class area table")

  class_lookup <- data.frame(
    class_code = 1:4,
    class_name = c(
      "Unsuitable",
      "Low suitability",
      "Moderate suitability",
      "High suitability"
    ),
    class_definition = c(
      paste0("Suitability < ", format(mtsps, trim = TRUE, scientific = FALSE)),
      paste0(
        format(mtsps, trim = TRUE, scientific = FALSE),
        " <= Suitability < 0.4"
      ),
      "0.4 <= Suitability < 0.6",
      "Suitability >= 0.6"
    ),
    stringsAsFactors = FALSE
  )

  observed <- data.frame(
    class_code = as.integer(class_area[[value_col]]),
    area_km2 = as.numeric(class_area[[area_col]]),
    stringsAsFactors = FALSE
  )

  output <- merge(class_lookup, observed, by = "class_code", all.x = TRUE)
  output$area_km2[is.na(output$area_km2)] <- 0
  output$area_10_4_km2 <- output$area_km2 / 10000

  total <- sum(output$area_km2)

  output$percent_of_valid_environment <- if (total > 0) {
    100 * output$area_km2 / total
  } else {
    NA_real_
  }

  output$map_name <- map_name
  output$MTSPS_threshold <- mtsps

  output <- output[
    c(
      "map_name", "class_code", "class_name", "class_definition",
      "MTSPS_threshold", "area_km2", "area_10_4_km2",
      "percent_of_valid_environment"
    )
  ]

  output
}

# ---------------------------------------------------------------------
# 4. VERIFY INPUTS
# ---------------------------------------------------------------------

required_paths <- c(
  occurrence_csv = occurrence_csv,
  background_csv = background_csv,
  environmental_ASC_directory = environmental_asc_dir,
  current_predictor_directory = current_predictor_dir,
  selected_models_csv = selected_models_csv,
  maxent_jar = maxent_jar
)

missing <- required_paths[
  !file.exists(required_paths) & !dir.exists(required_paths)
]

if (length(missing) > 0) {
  stop(
    "\nMissing input(s):\n",
    paste(names(missing), ": ", missing, collapse = "\n"),
    call. = FALSE
  )
}

if (identical(Sys.which("java"), "")) {
  stop("\nJava is not available on system PATH.", call. = FALSE)
}

selected_models <- read.csv(
  selected_models_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

selected_model_col <- find_column(
  selected_models,
  c("Model", "model_name"),
  "selected_models.csv"
)

if (!expected_selected_model %in% selected_models[[selected_model_col]]) {
  stop(
    "\nThe selected kuenm model was not found:\n",
    expected_selected_model,
    call. = FALSE
  )
}

env_asc_paths <- stats::setNames(
  file.path(environmental_asc_dir, paste0(selected_predictors, ".asc")),
  selected_predictors
)

if (any(!file.exists(env_asc_paths))) {
  stop(
    "\nMissing predictor ASC file(s):\n",
    paste(env_asc_paths[!file.exists(env_asc_paths)], collapse = "\n"),
    call. = FALSE
  )
}

current_tif_paths <- file.path(
  current_predictor_dir,
  paste0(selected_predictors, ".tif")
)

if (any(!file.exists(current_tif_paths))) {
  stop(
    "\nMissing predictor GeoTIFF(s):\n",
    paste(current_tif_paths[!file.exists(current_tif_paths)], collapse = "\n"),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 5. PREPARE OUTPUT FOLDERS
# ---------------------------------------------------------------------

if (clean_previous_run) {
  safe_remove(output_root)
}

cv_dir <- file.path(output_root, "01_CV100_MaxEnt")
final_dir <- file.path(output_root, "02_final_full_data_MaxEnt")
bias_dir <- file.path(output_root, "03_weighted_background_bias")
metrics_dir <- file.path(output_root, "04_CV_metrics_MTSPS_and_areas")
maps_dir <- file.path(output_root, "05_current_suitability_maps")
logs_dir <- file.path(output_root, "06_logs")

for (folder in c(cv_dir, final_dir, bias_dir, metrics_dir, maps_dir, logs_dir)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

cat("\n======================================================================\n")
cat("FINAL MAXENT CV100 + 4-CLASS AREA ANALYSIS — FIXED VERSION\n")
cat("======================================================================\n")
cat("Tuned configuration: RM = 1; FC = lph\n")
cat("CV replicates: ", n_cross_validation_replicates, "\n", sep = "")
cat("Output folder:\n", normalizePath(output_root, winslash = "/"), "\n", sep = "")

# ---------------------------------------------------------------------
# 6. READ OCCURRENCE AND BACKGROUND DATA
# ---------------------------------------------------------------------

occurrences <- read_xy_csv(
  occurrence_csv,
  "final occurrence CSV",
  add_species = TRUE
)

background <- read_xy_csv(
  background_csv,
  "Step-08 background CSV",
  add_species = FALSE
)

if (nrow(occurrences) != 828L) {
  stop(
    "\nExpected 828 final occurrences but found ", nrow(occurrences),
    ". Check the Step-05 final input.",
    call. = FALSE
  )
}

if (nrow(background) != 10000L) {
  stop(
    "\nExpected 10,000 Step-08 background points but found ", nrow(background),
    ".",
    call. = FALSE
  )
}

maxent_samples_csv <- file.path(
  logs_dir,
  "Camellia_japonica_final_occurrences_828.csv"
)

write.csv(occurrences, maxent_samples_csv, row.names = FALSE)

write.csv(
  data.frame(
    item = c("Final occurrences", "Step-08 background points"),
    n = c(nrow(occurrences), nrow(background)),
    stringsAsFactors = FALSE
  ),
  file.path(logs_dir, "input_record_counts.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 7. REBUILD THE MAXENT-WEIGHTED BIAS FILE
# ---------------------------------------------------------------------

weighted_bias_asc <- file.path(
  bias_dir,
  "Camellia_japonica_global_weighted_background_bias.asc"
)

weighted_bias_log <- file.path(
  bias_dir,
  "weighted_bias_validation_summary.csv"
)

cat("\nBuilding a validated weighted global background bias raster ...\n")

weighted_bias <- build_weighted_bias_raster(
  env_asc_paths = env_asc_paths,
  background_xy = background,
  occurrence_xy = occurrences,
  output_asc = weighted_bias_asc,
  output_log_csv = weighted_bias_log
)

# ---------------------------------------------------------------------
# 8. RUN MAXENT AVERAGED ACROSS 100 CROSS-VALIDATION REPLICATES
# ---------------------------------------------------------------------

cv_log <- file.path(logs_dir, "01_CV100_MaxEnt.log")

cat("\nRunning MaxEnt averaged across 100 cross-validation replicates ...\n")

cv_results_csv <- run_maxent(
  samples_csv = maxent_samples_csv,
  output_dir = cv_dir,
  environmental_dir = environmental_asc_dir,
  bias_asc = weighted_bias_asc,
  log_file = cv_log,
  model_type = "cv100"
)

cv_metrics <- extract_cv_summary(
  cv_results_csv,
  expected_replicates = n_cross_validation_replicates
)

write.csv(
  cv_metrics$replicates,
  file.path(metrics_dir, "01_CV100_replicate_metrics_raw.csv"),
  row.names = FALSE
)

write.csv(
  cv_metrics$summary,
  file.path(metrics_dir, "02_CV100_metrics_mean_SD_min_max.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 9. RUN THE FINAL FULL-DATA MODEL FOR THE SINGLE MTSPS THRESHOLD
# ---------------------------------------------------------------------

final_log <- file.path(logs_dir, "02_final_full_data_MaxEnt.log")

cat("\nRunning final full-data MaxEnt model ...\n")

final_results_csv <- run_maxent(
  samples_csv = maxent_samples_csv,
  output_dir = final_dir,
  environmental_dir = environmental_asc_dir,
  bias_asc = weighted_bias_asc,
  log_file = final_log,
  model_type = "full_data"
)

# ---------------------------------------------------------------------
# 10. LOAD CONTINUOUS PREDICTION GRIDS
# ---------------------------------------------------------------------

reference_tif <- terra::rast(current_tif_paths[1])

cv_average_asc <- find_prediction_grid(cv_dir, "cv_mean")
cv_average <- restore_crs(cv_average_asc, reference_tif)

cv_sd_asc <- find_prediction_grid(cv_dir, "cv_sd")
cv_sd <- if (is.na(cv_sd_asc)) NULL else restore_crs(cv_sd_asc, reference_tif)

final_prediction_asc <- find_prediction_grid(final_dir, "final")
final_prediction <- restore_crs(final_prediction_asc, reference_tif)

# ---------------------------------------------------------------------
# 11. DERIVE ONE FINAL MTSPS THRESHOLD
# ---------------------------------------------------------------------

background_cells <- terra::cellFromXY(
  reference_tif,
  as.matrix(background[, c("longitude", "latitude")])
)

occurrence_cells <- terra::cellFromXY(
  reference_tif,
  as.matrix(occurrences[, c("longitude", "latitude")])
)

negative_background <- background[
  !background_cells %in% occurrence_cells,
  ,
  drop = FALSE
]

if (nrow(negative_background) < 9000L) {
  stop(
    "\nFewer than 9,000 independent background points remain after removing ",
    "cells shared with occurrence points.",
    call. = FALSE
  )
}

presence_scores <- extract_scores(
  final_prediction,
  occurrences,
  "final occurrences"
)

background_scores <- extract_scores(
  final_prediction,
  negative_background,
  "independent Step-08 background"
)

mtsps <- mtsps_threshold(presence_scores, background_scores)
final_threshold <- mtsps$selected$threshold
threshold_metrics <- binary_metrics(
  presence_scores,
  background_scores,
  final_threshold
)

threshold_summary <- data.frame(
  species = species_name,
  model_configuration = "RM = 1; FC = lph",
  output_scale = "cloglog",
  threshold_method = "Maximum training sensitivity plus specificity (MTSPS)",
  MTSPS_threshold = final_threshold,
  sensitivity = threshold_metrics$sensitivity,
  specificity = threshold_metrics$specificity,
  TSS = threshold_metrics$TSS,
  Kappa = threshold_metrics$Kappa,
  n_occurrences = nrow(occurrences),
  n_negative_background_points = nrow(negative_background),
  stringsAsFactors = FALSE
)

write.csv(
  threshold_summary,
  file.path(metrics_dir, "03_FINAL_single_MTSPS_threshold.csv"),
  row.names = FALSE
)

write.csv(
  mtsps$curve,
  file.path(metrics_dir, "04_FINAL_MTSPS_threshold_curve.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 12. EXPORT CONTINUOUS MAPS
# ---------------------------------------------------------------------

cv_average_tif <- file.path(
  maps_dir,
  "Camellia_japonica_current_CV100_average_cloglog.tif"
)

final_prediction_tif <- file.path(
  maps_dir,
  "Camellia_japonica_current_final_full_data_cloglog.tif"
)

cv_sd_tif <- file.path(
  maps_dir,
  "Camellia_japonica_current_CV100_standard_deviation.tif"
)

terra::writeRaster(
  cv_average,
  cv_average_tif,
  overwrite = TRUE,
  wopt = list(gdal = "COMPRESS=LZW")
)

terra::writeRaster(
  final_prediction,
  final_prediction_tif,
  overwrite = TRUE,
  wopt = list(gdal = "COMPRESS=LZW")
)

if (!is.null(cv_sd)) {
  terra::writeRaster(
    cv_sd,
    cv_sd_tif,
    overwrite = TRUE,
    wopt = list(gdal = "COMPRESS=LZW")
  )
}

# ---------------------------------------------------------------------
# 13. CREATE FOUR SUITABILITY CLASSES
# ---------------------------------------------------------------------

cv_classes <- classify_four_classes(cv_average, final_threshold)
final_classes <- classify_four_classes(final_prediction, final_threshold)

names(cv_classes) <- "suitability_class"
names(final_classes) <- "suitability_class"

cv_class_tif <- file.path(
  maps_dir,
  "Camellia_japonica_current_CV100_average_MTSPS_4classes.tif"
)

final_class_tif <- file.path(
  maps_dir,
  "Camellia_japonica_current_final_full_data_MTSPS_4classes.tif"
)

terra::writeRaster(
  cv_classes,
  cv_class_tif,
  overwrite = TRUE,
  wopt = list(gdal = "COMPRESS=LZW")
)

terra::writeRaster(
  final_classes,
  final_class_tif,
  overwrite = TRUE,
  wopt = list(gdal = "COMPRESS=LZW")
)

class_legend <- data.frame(
  class_code = 1:4,
  class_name = c(
    "Unsuitable",
    "Low suitability",
    "Moderate suitability",
    "High suitability"
  ),
  class_definition = c(
    paste0("Suitability < ", format(final_threshold, scientific = FALSE, trim = TRUE)),
    paste0(
      format(final_threshold, scientific = FALSE, trim = TRUE),
      " <= Suitability < 0.4"
    ),
    "0.4 <= Suitability < 0.6",
    "Suitability >= 0.6"
  ),
  stringsAsFactors = FALSE
)

write.csv(
  class_legend,
  file.path(metrics_dir, "05_suitability_4class_legend.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 14. CALCULATE GEOGRAPHIC AREA OF FOUR CLASSES
# ---------------------------------------------------------------------

cv_area <- calculate_class_areas(
  cv_classes,
  "CV100 average suitability",
  final_threshold
)

final_area <- calculate_class_areas(
  final_classes,
  "Final full-data suitability",
  final_threshold
)

four_class_areas <- rbind(cv_area, final_area)

write.csv(
  four_class_areas,
  file.path(metrics_dir, "06_four_class_area_by_map.csv"),
  row.names = FALSE
)

suitable_area_summary <- do.call(
  rbind,
  lapply(
    split(four_class_areas, four_class_areas$map_name),
    function(table) {
      total_area <- sum(table$area_km2)
      suitable_area <- sum(table$area_km2[table$class_code >= 2])

      data.frame(
        map_name = unique(table$map_name),
        total_valid_environment_km2 = total_area,
        total_suitable_km2 = suitable_area,
        total_suitable_10_4_km2 = suitable_area / 10000,
        percent_suitable = 100 * suitable_area / total_area,
        stringsAsFactors = FALSE
      )
    }
  )
)

write.csv(
  suitable_area_summary,
  file.path(metrics_dir, "07_total_suitable_area_summary.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 15. SAVE FINAL MAXENT RESULT TABLE + VARIABLE IMPORTANCE
# ---------------------------------------------------------------------

file.copy(
  final_results_csv,
  file.path(metrics_dir, "08_final_full_data_maxentResults.csv"),
  overwrite = TRUE
)

final_results <- read.csv(
  final_results_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

result_names <- names(final_results)

importance <- do.call(
  rbind,
  lapply(selected_predictors, function(predictor) {
    contribution_col <- result_names[
      grepl(
        paste0("^", predictor, "\\s+contribution$"),
        result_names,
        ignore.case = TRUE
      )
    ]

    permutation_col <- result_names[
      grepl(
        paste0("^", predictor, "\\s+permutation\\s+importance$"),
        result_names,
        ignore.case = TRUE
      )
    ]

    data.frame(
      predictor = predictor,
      percent_contribution = if (length(contribution_col) == 1) {
        suppressWarnings(as.numeric(final_results[[contribution_col]][1]))
      } else NA_real_,
      permutation_importance = if (length(permutation_col) == 1) {
        suppressWarnings(as.numeric(final_results[[permutation_col]][1]))
      } else NA_real_,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  importance,
  file.path(metrics_dir, "09_final_variable_importance.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 16. REPORT
# ---------------------------------------------------------------------

cat("\n======================================================================\n")
cat("FIXED CV100 MAXENT + MTSPS + 4-CLASS ANALYSIS COMPLETED\n")
cat("======================================================================\n")
cat(
  "Final MTSPS threshold (cloglog): ",
  format(round(final_threshold, 8), nsmall = 8),
  "\n",
  sep = ""
)

cat("\nKey files:\n")
cat(
  "CV100 mean/SD metrics:\n",
  normalizePath(
    file.path(metrics_dir, "02_CV100_metrics_mean_SD_min_max.csv"),
    winslash = "/"
  ),
  "\n",
  sep = ""
)

cat(
  "Four-class area table:\n",
  normalizePath(
    file.path(metrics_dir, "06_four_class_area_by_map.csv"),
    winslash = "/"
  ),
  "\n",
  sep = ""
)

cat(
  "CV100 average four-class map:\n",
  normalizePath(cv_class_tif, winslash = "/"),
  "\n",
  sep = ""
)

cat(
  "Final full-data four-class map:\n",
  normalizePath(final_class_tif, winslash = "/"),
  "\n",
  sep = ""
)
