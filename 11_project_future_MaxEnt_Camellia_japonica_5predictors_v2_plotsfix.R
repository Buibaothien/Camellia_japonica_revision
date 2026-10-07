# =====================================================================
# 11_PROJECT_FUTURE_MAXENT_CAMELLIA_JAPONICA_V2_PLOTSFIX.R
# =====================================================================
# PURPOSE
# -------
# Project the final calibrated MaxEnt model of Camellia japonica onto all
# future climate scenarios:
#
#   GCMs: ACCESS-CM2, BCC-CSM2-MR, CMCC-ESM2
#   SSPs: SSP126, SSP245, SSP370, SSP585
#   Periods: 2021-2040, 2041-2060, 2061-2080, 2081-2100
#
# Total expected continuous projections: 3 x 4 x 4 = 48.
#
# FINAL MODEL CONFIGURATION
# -------------------------
# RM = 1
# FC = lph
# Predictors = Bio2, Bio6, Bio14, Bio18, Elv
# Output scale = cloglog
#
# IMPORTANT SCIENTIFIC DESIGN
# ---------------------------
# 1. This script does NOT derive new thresholds and does NOT classify maps.
# 2. It writes continuous suitability rasters for each individual GCM.
# 3. The next script will create a three-GCM continuous ensemble for every
#    SSP x period. Classification and four-class area calculation occur only
#    after ensemble creation.
# 4. The single current MTSPS threshold will be applied later to the ensemble
#    rasters, not recalculated for future climates.
#
# FUTURE INPUT STRUCTURE EXPECTED
# -------------------------------
# E:/PREDICTION TREE/Camellia_japonica/future_environmental/
#   ACCESS-CM2/
#     ssp126_2021_2040/
#       wc2.1_2.5m_bioc_ACCESS-CM2_ssp126_2021-2040.tif  [19 BIO bands]
#       elv.tif
#     ...
#   BCC-CSM2-MR/
#   CMCC-ESM2/
#
# Each scenario directory must contain:
# - one 19-band WorldClim "bioc" GeoTIFF, where bands 1-19 = Bio1-Bio19;
# - one elevation GeoTIFF.
#
# OUTPUT
# ------
# global_reanalysis_outputs/11_future_projection_individual_GCM/
#   01_individual_GCM_continuous/
#     ACCESS-CM2/ssp126_2021_2040/
#       Camellia_japonica_ACCESS-CM2_ssp126_2021_2040_cloglog.tif
#       01_maxent_raw_outputs/
#       02_metadata/
#   02_projection_logs/
#   03_projection_manifests/
#
# =====================================================================
# PATCH IN V2
# -----------
# MaxEnt attempted to write its automatic prediction PNG to a missing
# output_directory/plots/ folder, causing all 48 scenarios to fail.
# This version pre-creates plots/ and html/ before every MaxEnt call and
# stops immediately after any first failure.
# =====================================================================

# ---------------------------------------------------------------------
# 0. INSTALL AND LOAD REQUIRED PACKAGE
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

# Root folder shown in your File Explorer screenshots.
future_root_dir <- file.path(
  project_dir,
  "future_environmental"
)

# Folder that contains the successfully completed final current MaxEnt model.
# This is the output folder printed in your completed CV100 console.
final_current_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

# Workspace that contains the five current ASCII predictors used by kuenm.
kuenm_workdir <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "09_kuenm_global_OR10"
)

# Folder containing maxent.jar.
maxent_dir <- "D:/Prediction_tree/maxent"
maxent_jar <- file.path(maxent_dir, "maxent.jar")

# New output folder created by this script.
output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "11_future_projection_individual_GCM_v2_plotsfix"
)

# Use TRUE only for a completely fresh rerun of all future projections.
# This deletes only output_root; it never changes current model outputs or
# source future environmental data.
clean_output_root <- FALSE

# If TRUE, a completed suitability GeoTIFF is not regenerated.
# This permits safe resume after interruption.
skip_completed_scenarios <- TRUE

# Stop immediately after the first failed scenario, preventing a repeated
# systemic failure from consuming hours across all remaining scenarios.
stop_immediately_after_first_failure <- TRUE

# Delete temporary future ASCII predictor folders after a successful scenario.
# Keeping them is unnecessary and consumes substantial disk space.
delete_temporary_future_ASC_after_success <- TRUE

# Preserve raw MaxEnt output folders for traceability, including MESS/clamp.
keep_raw_maxent_outputs <- TRUE

# Optional MaxEnt threads. Set to 1 if memory is limited.
maxent_threads <- 4L

# Memory allocated to Java/MaxEnt in MB.
max_memory_mb <- 4096L

# FINAL tuned configuration.
regularization_multiplier <- 1
feature_class_label <- "lph"
selected_predictors <- c("Bio2", "Bio6", "Bio14", "Bio18", "Elv")

# Expected folder names.
expected_gcms <- c("ACCESS-CM2", "BCC-CSM2-MR", "CMCC-ESM2")
expected_ssps <- c("ssp126", "ssp245", "ssp370", "ssp585")
expected_periods <- c(
  "2021_2040",
  "2041_2060",
  "2061_2080",
  "2081_2100"
)

# ---------------------------------------------------------------------
# 2. INPUT PATHS FROM THE COMPLETED CURRENT-MODEL WORKFLOW
# ---------------------------------------------------------------------

# All 828 final global natural-land-cover occurrence records.
occurrence_csv <- file.path(
  kuenm_workdir,
  "02_occurrences",
  "Camellia_japonica_joint.csv"
)

# Current environmental layers used for model training.
current_ASC_dir <- file.path(
  kuenm_workdir,
  "01_global_variables",
  "set1_primary5"
)

# Validated weighted bias raster used in the final current model.
weighted_bias_ASC <- file.path(
  final_current_root,
  "03_weighted_background_bias",
  "Camellia_japonica_global_weighted_background_bias.asc"
)

# Final full-data model outputs; existence is checked to verify that the
# current model actually completed before future projection begins.
final_model_dir <- file.path(
  final_current_root,
  "02_final_full_data_MaxEnt"
)

final_model_results_csv <- file.path(
  final_model_dir,
  "maxentResults.csv"
)

final_model_lambdas <- file.path(
  final_model_dir,
  paste0(species_name, ".lambdas")
)

# Kuenm model-selection output used as an explicit safeguard.
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

  direct_hits <- match(normal_candidates, normal_names, nomatch = 0)
  direct_hits <- direct_hits[direct_hits > 0]

  if (length(direct_hits) > 0) {
    return(original_names[direct_hits[1]])
  }

  for (candidate in normal_candidates) {
    partial_hits <- grep(candidate, normal_names, fixed = TRUE)

    if (length(partial_hits) > 0) {
      return(original_names[partial_hits[1]])
    }
  }

  if (required) {
    stop(
      "\nCould not identify a required column in ", object_label, ".\n",
      "Available columns: ", paste(original_names, collapse = ", "),
      call. = FALSE
    )
  }

  NA_character_
}

read_occurrences_for_maxent <- function(path) {
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

  lon_col <- find_column(
    data,
    c("longitude", "decimalLongitude", "decimalLon", "lon", "x"),
    "occurrence CSV"
  )

  lat_col <- find_column(
    data,
    c("latitude", "decimalLatitude", "decimalLat", "lat", "y"),
    "occurrence CSV"
  )

  out <- data.frame(
    species = rep(species_name, nrow(data)),
    longitude = suppressWarnings(as.numeric(data[[lon_col]])),
    latitude = suppressWarnings(as.numeric(data[[lat_col]])),
    stringsAsFactors = FALSE
  )

  out <- out[
    is.finite(out$longitude) &
      is.finite(out$latitude) &
      out$longitude >= -180 &
      out$longitude <= 180 &
      out$latitude >= -90 &
      out$latitude <= 90,
    ,
    drop = FALSE
  ]

  out <- unique(out)

  if (nrow(out) == 0) {
    stop("No valid occurrence coordinates remain in:\n", path, call. = FALSE)
  }

  out
}

quote_path_parameter <- function(parameter_name, path, must_work = TRUE) {
  normalized_path <- normalizePath(
    path,
    winslash = "\\",
    mustWork = must_work
  )

  paste0(parameter_name, "=", shQuote(normalized_path))
}

assert_aligned <- function(reference_raster, candidate_raster, label) {
  is_aligned <- terra::compareGeom(
    reference_raster,
    candidate_raster,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )

  if (!is_aligned) {
    stop(
      "\nFuture predictor grid mismatch: ", label,
      "\nBio2, Bio6, Bio14, Bio18, and Elv must have the same CRS, extent, ",
      "resolution, rows, columns, and origin within every scenario.",
      call. = FALSE
    )
  }
}

list_tif_files <- function(directory) {
  list.files(
    directory,
    pattern = "\\.(tif|tiff)$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
}

detect_elevation_file <- function(tif_files, scenario_label) {
  stems <- normalise_name(tools::file_path_sans_ext(basename(tif_files)))

  elevation_indices <- which(
    stems %in% c("elv", "elev", "elevation", "dem") |
      grepl("elev|elevation|(^|[^a-z])elv", stems)
  )

  if (length(elevation_indices) == 1) {
    return(tif_files[elevation_indices])
  }

  if (length(elevation_indices) == 0) {
    stop(
      "\nNo elevation GeoTIFF was identified in scenario:\n",
      scenario_label,
      "\nExpected a filename such as elv.tif or elevation.tif.",
      call. = FALSE
    )
  }

  stop(
    "\nMore than one elevation-like GeoTIFF was found in scenario:\n",
    scenario_label,
    "\nCandidates:\n",
    paste(tif_files[elevation_indices], collapse = "\n"),
    call. = FALSE
  )
}

detect_multiband_bioc_file <- function(tif_files, elevation_file, scenario_label) {
  climate_files <- tif_files[
    normalizePath(tif_files, winslash = "/", mustWork = TRUE) !=
      normalizePath(elevation_file, winslash = "/", mustWork = TRUE)
  ]

  if (length(climate_files) == 0) {
    stop(
      "\nNo climate GeoTIFF remains after identifying elevation in:\n",
      scenario_label,
      call. = FALSE
    )
  }

  details <- lapply(climate_files, function(path) {
    raster <- terra::rast(path)

    data.frame(
      file = path,
      n_layers = terra::nlyr(raster),
      looks_like_bioc = grepl(
        "bioc|bioclim|bio",
        basename(path),
        ignore.case = TRUE
      ),
      stringsAsFactors = FALSE
    )
  })

  details <- do.call(rbind, details)

  candidates <- details[
    details$n_layers >= 19 & details$looks_like_bioc,
    ,
    drop = FALSE
  ]

  if (nrow(candidates) == 1) {
    return(candidates$file[1])
  }

  # Fallback: accept exactly one 19+ band raster even if its name lacks "bioc".
  fallback <- details[details$n_layers >= 19, , drop = FALSE]

  if (nrow(fallback) == 1) {
    return(fallback$file[1])
  }

  if (nrow(candidates) == 0 && nrow(fallback) == 0) {
    stop(
      "\nNo 19-band bioclimatic GeoTIFF was found in scenario:\n",
      scenario_label,
      "\nThe expected file should contain Bio1-Bio19 as bands 1-19.",
      call. = FALSE
    )
  }

  stop(
    "\nMore than one possible 19-band bioclimatic GeoTIFF was found in:\n",
    scenario_label,
    "\nCandidates:\n",
    paste(unique(c(candidates$file, fallback$file)), collapse = "\n"),
    call. = FALSE
  )
}

parse_scenario_folder <- function(folder_name) {
  standardized <- tolower(folder_name)

  match <- regexec(
    "^ssp(126|245|370|585)_([0-9]{4})_([0-9]{4})$",
    standardized
  )

  values <- regmatches(standardized, match)[[1]]

  if (length(values) != 4) {
    return(NULL)
  }

  list(
    scenario_folder = folder_name,
    ssp = paste0("SSP", values[2]),
    period_start = as.integer(values[3]),
    period_end = as.integer(values[4]),
    period_label = paste0(values[3], "-", values[4]),
    scenario_id = paste0(
      "SSP", values[2], "_", values[3], "_", values[4]
    )
  )
}

find_future_scenarios <- function(root_dir, gcms) {
  rows <- list()

  for (gcm in gcms) {
    gcm_dir <- file.path(root_dir, gcm)

    if (!dir.exists(gcm_dir)) {
      warning("Expected GCM folder is missing and will be skipped: ", gcm_dir)
      next
    }

    child_dirs <- list.dirs(
      gcm_dir,
      recursive = FALSE,
      full.names = TRUE
    )

    for (scenario_dir in child_dirs) {
      parsed <- parse_scenario_folder(basename(scenario_dir))

      if (is.null(parsed)) {
        next
      }

      rows[[length(rows) + 1]] <- data.frame(
        gcm = gcm,
        scenario_folder = parsed$scenario_folder,
        ssp = parsed$ssp,
        period_start = parsed$period_start,
        period_end = parsed$period_end,
        period_label = parsed$period_label,
        scenario_id = parsed$scenario_id,
        input_directory = normalizePath(
          scenario_dir,
          winslash = "/",
          mustWork = TRUE
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    stop(
      "\nNo scenario directories matching ssp###_YYYY_YYYY were found under:\n",
      root_dir,
      call. = FALSE
    )
  }

  scenarios <- do.call(rbind, rows)

  scenarios <- scenarios[
    order(
      match(scenarios$gcm, gcms),
      match(tolower(scenarios$ssp), tolower(expected_ssps)),
      scenarios$period_start
    ),
    ,
    drop = FALSE
  ]

  rownames(scenarios) <- NULL
  scenarios
}

build_future_ASC_predictors <- function(
  scenario_directory,
  scenario_label,
  temporary_ASC_directory,
  metadata_csv
) {
  tif_files <- list_tif_files(scenario_directory)

  if (length(tif_files) < 2) {
    stop(
      "\nScenario contains fewer than two GeoTIFFs:\n",
      scenario_directory,
      "\nExpected one 19-band bioclimatic file and one elevation file.",
      call. = FALSE
    )
  }

  elevation_file <- detect_elevation_file(tif_files, scenario_label)
  bioc_file <- detect_multiband_bioc_file(
    tif_files,
    elevation_file,
    scenario_label
  )

  bioc_stack <- terra::rast(bioc_file)

  if (terra::nlyr(bioc_stack) < 19) {
    stop(
      "\nThe selected bioclimatic file has fewer than 19 bands:\n",
      bioc_file,
      call. = FALSE
    )
  }

  # WorldClim bioc layer ordering is Bio1 through Bio19.
  future_bio2 <- bioc_stack[[2]]
  future_bio6 <- bioc_stack[[6]]
  future_bio14 <- bioc_stack[[14]]
  future_bio18 <- bioc_stack[[18]]
  future_elv <- terra::rast(elevation_file)

  assert_aligned(future_bio2, future_bio6, "Bio6")
  assert_aligned(future_bio2, future_bio14, "Bio14")
  assert_aligned(future_bio2, future_bio18, "Bio18")
  assert_aligned(future_bio2, future_elv, "Elv")

  future_stack <- c(
    future_bio2,
    future_bio6,
    future_bio14,
    future_bio18,
    future_elv
  )

  names(future_stack) <- selected_predictors

  if (
    is.na(terra::crs(future_stack)) ||
    terra::crs(future_stack) == ""
  ) {
    stop(
      "\nFuture predictor stack has no CRS in scenario:\n",
      scenario_label,
      call. = FALSE
    )
  }

  # Confirm that every selected layer contains at least some valid cells.
  valid_counts <- vapply(
    seq_len(terra::nlyr(future_stack)),
    function(index) {
      summary <- terra::global(
        !is.na(future_stack[[index]]),
        fun = "sum",
        na.rm = TRUE
      )
      as.numeric(summary[1, 1])
    },
    numeric(1)
  )

  if (any(!is.finite(valid_counts) | valid_counts <= 0)) {
    stop(
      "\nAt least one selected future predictor is entirely NoData in scenario:\n",
      scenario_label,
      call. = FALSE
    )
  }

  safe_remove(temporary_ASC_directory)
  dir.create(temporary_ASC_directory, recursive = TRUE, showWarnings = FALSE)

  asc_paths <- file.path(
    temporary_ASC_directory,
    paste0(selected_predictors, ".asc")
  )

  for (index in seq_along(selected_predictors)) {
    terra::writeRaster(
      future_stack[[index]],
      filename = asc_paths[index],
      filetype = "AAIGrid",
      NAflag = -9999,
      overwrite = TRUE
    )

    if (
      !file.exists(asc_paths[index]) ||
      file.info(asc_paths[index])$size <= 0
    ) {
      stop(
        "\nFailed to create future ASCII predictor:\n",
        asc_paths[index],
        call. = FALSE
      )
    }
  }

  metadata <- data.frame(
    predictor = selected_predictors,
    source_file = c(
      rep(normalizePath(bioc_file, winslash = "/", mustWork = TRUE), 4),
      normalizePath(elevation_file, winslash = "/", mustWork = TRUE)
    ),
    source_band = c(2, 6, 14, 18, 1),
    source_description = c(
      "WorldClim bioc band 2 = Bio2",
      "WorldClim bioc band 6 = Bio6",
      "WorldClim bioc band 14 = Bio14",
      "WorldClim bioc band 18 = Bio18",
      "Elevation GeoTIFF"
    ),
    generated_ASC = normalizePath(
      asc_paths,
      winslash = "/",
      mustWork = TRUE
    ),
    ncol = terra::ncol(future_stack),
    nrow = terra::nrow(future_stack),
    resolution_x = terra::res(future_stack)[1],
    resolution_y = terra::res(future_stack)[2],
    CRS = terra::crs(future_stack),
    valid_cell_count = valid_counts,
    stringsAsFactors = FALSE
  )

  write.csv(metadata, metadata_csv, row.names = FALSE)

  list(
    temporary_ASC_directory = temporary_ASC_directory,
    future_ASC_paths = asc_paths,
    bioc_file = bioc_file,
    elevation_file = elevation_file,
    reference_future_raster = future_bio2,
    metadata = metadata
  )
}

run_maxent_future_projection <- function(
  samples_csv,
  current_environmental_ASC_directory,
  future_projection_ASC_directory,
  output_directory,
  log_file
) {
  safe_remove(output_directory)
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  # MaxEnt automatically writes a main prediction PNG to output_directory/plots.
  # Create it explicitly because otherwise Windows can raise FileNotFoundException.
  dir.create(
    file.path(output_directory, "plots"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Create the HTML folder as well to prevent a similar output-path issue.
  dir.create(
    file.path(output_directory, "html"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  java_executable <- Sys.which("java")

  if (identical(java_executable, "")) {
    stop(
      "Java is not available on the system PATH. MaxEnt cannot be executed.",
      call. = FALSE
    )
  }

  # MaxEnt trains the deterministic final current model and projects that
  # calibrated model onto the five future ASCII layers specified here.
  parameters <- c(
    quote_path_parameter("samplesfile", samples_csv),
    quote_path_parameter(
      "environmentallayers",
      current_environmental_ASC_directory
    ),
    quote_path_parameter(
      "projectionlayers",
      future_projection_ASC_directory
    ),
    quote_path_parameter(
      "outputdirectory",
      output_directory,
      must_work = FALSE
    ),
    quote_path_parameter("biasfile", weighted_bias_ASC),

    # Required for the validated weighted bias surface used in the current run.
    "biasType=3",
    "maximumbackground=10000",

    # Fixed tuned configuration selected by kuenm.
    "linear=true",
    "quadratic=false",
    "product=true",
    "threshold=false",
    "hinge=true",
    "autofeature=false",
    paste0("betamultiplier=", regularization_multiplier),

    # Preserve the same background design used for final current modelling.
    "addsamplestobackground=false",
    "addallsamplestobackground=false",

    # Fixed deterministic full-data model.
    "replicates=1",
    "randomtestpoints=0",
    "randomseed=false",
    "removeduplicates=false",

    # Projection settings.
    "outputformat=cloglog",
    "outputfiletype=asc",
    "outputgrids=true",
    "writeclampgrid=true",
    "writemess=true",
    "doclamp=true",
    "extrapolate=true",
    "fadebyclamping=false",

    # Keep the standard prediction figure inside the pre-created plots folder.
    # Response curves and jackknife outputs remain disabled.
    "responsecurves=false",
    "jackknife=false",
    "plots=true",
    "writeplotdata=false",
    "writebackgroundpredictions=false",

    "askoverwrite=false",
    "visible=false",
    "cache=true",
    "appendtoresultsfile=false",
    paste0("threads=", maxent_threads)
  )

  args <- c(
    paste0("-mx", max_memory_mb, "m"),
    "-jar",
    shQuote(normalizePath(maxent_jar, winslash = "\\", mustWork = TRUE)),
    parameters
  )

  command_preview <- paste(
    shQuote(java_executable),
    paste(args, collapse = " ")
  )

  writeLines(
    command_preview,
    file.path(output_directory, "maxent_command.txt")
  )

  exit_status <- system2(
    command = java_executable,
    args = args,
    stdout = log_file,
    stderr = log_file,
    wait = TRUE
  )

  results_csv <- file.path(output_directory, "maxentResults.csv")

  if (exit_status != 0 || !file.exists(results_csv)) {
    stop(
      "\nMaxEnt projection failed.\n",
      "Read the scenario log:\n",
      log_file,
      call. = FALSE
    )
  }

  results_csv
}

find_projected_ASC <- function(maxent_output_directory, scenario_folder_name) {
  asc_files <- list.files(
    maxent_output_directory,
    pattern = "\\.asc$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  if (length(asc_files) == 0) {
    stop(
      "No ASCII grids were written in:\n",
      maxent_output_directory,
      call. = FALSE
    )
  }

  names_only <- basename(asc_files)
  standard_current_name <- paste0(species_name, ".asc")

  # MaxEnt appends the projection directory name to the species output name.
  expected_name <- paste0(species_name, "_", scenario_folder_name, ".asc")

  exact_match <- asc_files[
    tolower(names_only) == tolower(expected_name)
  ]

  if (length(exact_match) == 1) {
    return(exact_match)
  }

  # Robust fallback for minor naming differences.
  projection_candidates <- asc_files[
    grepl(
      scenario_folder_name,
      names_only,
      ignore.case = TRUE,
      fixed = TRUE
    ) &
      !grepl(
        "clamp|mess|explain|std|avg",
        names_only,
        ignore.case = TRUE
      )
  ]

  if (length(projection_candidates) == 1) {
    return(projection_candidates)
  }

  # Final fallback: identify the non-current, non-auxiliary suitability ASC.
  non_auxiliary <- asc_files[
    tolower(names_only) != tolower(standard_current_name) &
      !grepl(
        "clamp|mess|explain|std|avg",
        names_only,
        ignore.case = TRUE
      )
  ]

  if (length(non_auxiliary) == 1) {
    return(non_auxiliary)
  }

  stop(
    "\nCould not uniquely identify the projected suitability ASC for scenario ",
    scenario_folder_name,
    ".\nAvailable ASC files:\n",
    paste(names_only, collapse = "\n"),
    call. = FALSE
  )
}

find_optional_auxiliary_ASC <- function(maxent_output_directory, keyword) {
  asc_files <- list.files(
    maxent_output_directory,
    pattern = "\\.asc$",
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )

  candidates <- asc_files[
    grepl(keyword, basename(asc_files), ignore.case = TRUE)
  ]

  if (length(candidates) == 1) {
    candidates
  } else {
    NA_character_
  }
}

write_projection_geotiff <- function(
  projected_ASC,
  projection_reference_raster,
  output_tif
) {
  prediction <- terra::rast(projected_ASC)

  # ASCII grids commonly lack CRS metadata. Restore the future grid CRS.
  if (is.na(terra::crs(prediction)) || terra::crs(prediction) == "") {
    terra::crs(prediction) <- terra::crs(projection_reference_raster)
  }

  # Ensure that MaxEnt output has the exact projection geometry.
  if (!terra::compareGeom(
    prediction,
    projection_reference_raster,
    crs = TRUE,
    ext = TRUE,
    rowcol = TRUE,
    res = TRUE,
    stopOnError = FALSE
  )) {
    stop(
      "\nProjected MaxEnt raster geometry does not match the future predictor grid:\n",
      projected_ASC,
      call. = FALSE
    )
  }

  values_range <- terra::global(
    prediction,
    fun = c("min", "max"),
    na.rm = TRUE
  )

  if (
    !is.finite(values_range[1, "min"]) ||
    !is.finite(values_range[1, "max"])
  ) {
    stop(
      "\nProjected suitability grid contains no finite values:\n",
      projected_ASC,
      call. = FALSE
    )
  }

  if (
    values_range[1, "min"] < -1e-8 ||
    values_range[1, "max"] > 1 + 1e-8
  ) {
    warning(
      "Projected suitability values fall outside the expected cloglog interval [0, 1]: ",
      basename(projected_ASC)
    )
  }

  terra::writeRaster(
    prediction,
    output_tif,
    overwrite = TRUE,
    wopt = list(gdal = "COMPRESS=LZW")
  )

  list(
    output_tif = output_tif,
    min_suitability = values_range[1, "min"],
    max_suitability = values_range[1, "max"],
    ncell = terra::ncell(prediction),
    n_non_NA = as.numeric(
      terra::global(!is.na(prediction), fun = "sum", na.rm = TRUE)[1, 1]
    )
  )
}

write_status <- function(status_data, output_csv) {
  if (nrow(status_data) > 0) {
    write.csv(status_data, output_csv, row.names = FALSE)
  }
}

# ---------------------------------------------------------------------
# 4. VERIFY REQUIRED INPUTS
# ---------------------------------------------------------------------

required_paths <- c(
  future_root_dir = future_root_dir,
  occurrence_csv = occurrence_csv,
  current_ASC_dir = current_ASC_dir,
  weighted_bias_ASC = weighted_bias_ASC,
  final_model_results_csv = final_model_results_csv,
  final_model_lambdas = final_model_lambdas,
  selected_models_csv = selected_models_csv,
  maxent_jar = maxent_jar
)

missing_paths <- required_paths[
  !file.exists(required_paths) & !dir.exists(required_paths)
]

if (length(missing_paths) > 0) {
  stop(
    "\nMissing required path(s):\n",
    paste(names(missing_paths), ": ", missing_paths, collapse = "\n"),
    call. = FALSE
  )
}

if (identical(Sys.which("java"), "")) {
  stop(
    "\nJava is not available on the system PATH. Configure Java before running MaxEnt.",
    call. = FALSE
  )
}

selected_models <- read.csv(
  selected_models_csv,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

model_col <- find_column(
  selected_models,
  c("Model", "model_name"),
  "selected_models.csv"
)

if (!expected_selected_model %in% selected_models[[model_col]]) {
  stop(
    "\nThe expected kuenm-selected model was not found:\n",
    expected_selected_model,
    call. = FALSE
  )
}

current_ASC_paths <- file.path(
  current_ASC_dir,
  paste0(selected_predictors, ".asc")
)

if (any(!file.exists(current_ASC_paths))) {
  stop(
    "\nThe following current ASCII predictor(s) are missing:\n",
    paste(
      current_ASC_paths[!file.exists(current_ASC_paths)],
      collapse = "\n"
    ),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------
# 5. PREPARE OUTPUT DIRECTORY STRUCTURE
# ---------------------------------------------------------------------

if (clean_output_root) {
  safe_remove(output_root)
}

individual_projection_root <- file.path(
  output_root,
  "01_individual_GCM_continuous"
)

temporary_ASC_root <- file.path(
  output_root,
  "02_temporary_future_ASC"
)

logs_root <- file.path(
  output_root,
  "03_projection_logs"
)

manifest_root <- file.path(
  output_root,
  "04_projection_manifests"
)

for (folder in c(
  individual_projection_root,
  temporary_ASC_root,
  logs_root,
  manifest_root
)) {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
}

manifest_csv <- file.path(
  manifest_root,
  "future_projection_status_manifest.csv"
)

scenario_inventory_csv <- file.path(
  manifest_root,
  "future_scenario_inventory.csv"
)

run_settings_csv <- file.path(
  manifest_root,
  "future_projection_run_settings.csv"
)

# ---------------------------------------------------------------------
# 6. READ OCCURRENCES AND DISCOVER ALL FUTURE SCENARIOS
# ---------------------------------------------------------------------

occurrences <- read_occurrences_for_maxent(occurrence_csv)

if (nrow(occurrences) != 828L) {
  stop(
    "\nExpected 828 final occurrence records but found ", nrow(occurrences),
    ". Check the completed global workflow input.",
    call. = FALSE
  )
}

samples_for_projection_csv <- file.path(
  manifest_root,
  "Camellia_japonica_final_occurrences_828.csv"
)

write.csv(
  occurrences,
  samples_for_projection_csv,
  row.names = FALSE
)

future_scenarios <- find_future_scenarios(
  future_root_dir,
  expected_gcms
)

write.csv(
  future_scenarios,
  scenario_inventory_csv,
  row.names = FALSE
)

expected_total <- length(expected_gcms) *
  length(expected_ssps) *
  length(expected_periods)

cat("\n======================================================================\n")
cat("FUTURE MAXENT PROJECTION — Camellia japonica\n")
cat("======================================================================\n")
cat("Tuned MaxEnt configuration: RM = 1; FC = lph\n")
cat("Predictors: ", paste(selected_predictors, collapse = ", "), "\n", sep = "")
cat("Detected future scenarios: ", nrow(future_scenarios), "\n", sep = "")
cat("Expected complete design: ", expected_total, " (3 GCM x 4 SSP x 4 periods)\n", sep = "")
cat("Output folder:\n", normalizePath(output_root, winslash = "/"), "\n", sep = "")

if (nrow(future_scenarios) != expected_total) {
  warning(
    "\nDetected ", nrow(future_scenarios),
    " scenarios rather than the expected ", expected_total,
    ". The script will run all detected valid scenario folders and record any ",
    "missing GCM/SSP/period combinations in the manifest."
  )
}

write.csv(
  data.frame(
    setting = c(
      "Species",
      "Final MaxEnt RM",
      "Final MaxEnt FC",
      "Predictors",
      "Current model source",
      "Current model lambdas source",
      "Future environmental root",
      "Expected GCMs",
      "Expected SSPs",
      "Expected time periods",
      "Output scale",
      "Future map classification",
      "MTSPS handling"
    ),
    value = c(
      species_name,
      regularization_multiplier,
      feature_class_label,
      paste(selected_predictors, collapse = "; "),
      normalizePath(final_model_results_csv, winslash = "/"),
      normalizePath(final_model_lambdas, winslash = "/"),
      normalizePath(future_root_dir, winslash = "/"),
      paste(expected_gcms, collapse = "; "),
      paste(toupper(expected_ssps), collapse = "; "),
      paste(gsub("_", "-", expected_periods), collapse = "; "),
      "cloglog",
      "Not performed in Script 11",
      "Applied after continuous three-GCM ensemble in Script 13"
    ),
    stringsAsFactors = FALSE
  ),
  run_settings_csv,
  row.names = FALSE
)

# ---------------------------------------------------------------------
# 7. PROJECT ALL INDIVIDUAL GCM SCENARIOS
# ---------------------------------------------------------------------

status_records <- list()
start_all <- Sys.time()

for (scenario_index in seq_len(nrow(future_scenarios))) {
  scenario <- future_scenarios[scenario_index, , drop = FALSE]

  gcm_name <- scenario$gcm
  scenario_folder <- scenario$scenario_folder
  scenario_id <- scenario$scenario_id
  scenario_label <- paste(gcm_name, scenario_id, sep = " | ")

  scenario_output_dir <- file.path(
    individual_projection_root,
    gcm_name,
    scenario_folder
  )

  raw_maxent_output_dir <- file.path(
    scenario_output_dir,
    "01_maxent_raw_outputs"
  )

  metadata_dir <- file.path(
    scenario_output_dir,
    "02_metadata"
  )

  output_tif <- file.path(
    scenario_output_dir,
    paste0(
      species_name, "_",
      gcm_name, "_",
      scenario_folder,
      "_cloglog.tif"
    )
  )

  scenario_log <- file.path(
    logs_root,
    paste0(gcm_name, "_", scenario_folder, "_MaxEnt.log")
  )

  temporary_ASC_dir <- file.path(
    temporary_ASC_root,
    gcm_name,
    scenario_folder
  )

  predictor_metadata_csv <- file.path(
    metadata_dir,
    "future_predictor_manifest.csv"
  )

  scenario_start <- Sys.time()

  cat(
    "\n[", scenario_index, "/", nrow(future_scenarios), "] ",
    scenario_label,
    "\n",
    sep = ""
  )

  if (
    skip_completed_scenarios &&
    file.exists(output_tif) &&
    file.info(output_tif)$size > 0
  ) {
    cat("  Status: output GeoTIFF already exists; skipped.\n")

    status_records[[length(status_records) + 1]] <- data.frame(
      gcm = gcm_name,
      ssp = scenario$ssp,
      period_start = scenario$period_start,
      period_end = scenario$period_end,
      scenario_folder = scenario_folder,
      scenario_id = scenario_id,
      input_directory = scenario$input_directory,
      status = "skipped_existing_output",
      output_continuous_tif = normalizePath(
        output_tif,
        winslash = "/",
        mustWork = TRUE
      ),
      projected_ASC = NA_character_,
      source_bioc_tif = NA_character_,
      source_elv_tif = NA_character_,
      min_suitability = NA_real_,
      max_suitability = NA_real_,
      valid_prediction_cells = NA_real_,
      elapsed_minutes = 0,
      message = "",
      stringsAsFactors = FALSE
    )

    write_status(
      do.call(rbind, status_records),
      manifest_csv
    )

    next
  }

  result <- tryCatch(
    {
      dir.create(scenario_output_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

      future_inputs <- build_future_ASC_predictors(
        scenario_directory = scenario$input_directory,
        scenario_label = scenario_label,
        temporary_ASC_directory = temporary_ASC_dir,
        metadata_csv = predictor_metadata_csv
      )

      maxent_results_csv <- run_maxent_future_projection(
        samples_csv = samples_for_projection_csv,
        current_environmental_ASC_directory = current_ASC_dir,
        future_projection_ASC_directory =
          future_inputs$temporary_ASC_directory,
        output_directory = raw_maxent_output_dir,
        log_file = scenario_log
      )

      projected_ASC <- find_projected_ASC(
        raw_maxent_output_dir,
        scenario_folder_name = basename(
          future_inputs$temporary_ASC_directory
        )
      )

      projection_info <- write_projection_geotiff(
        projected_ASC = projected_ASC,
        projection_reference_raster =
          future_inputs$reference_future_raster,
        output_tif = output_tif
      )

      # Export optional MESS and clamp rasters as GeoTIFF whenever MaxEnt
      # creates them. These are useful for novelty/extrapolation assessment.
      mess_ASC <- find_optional_auxiliary_ASC(
        raw_maxent_output_dir,
        "mess"
      )

      clamp_ASC <- find_optional_auxiliary_ASC(
        raw_maxent_output_dir,
        "clamp"
      )

      optional_outputs <- data.frame(
        output_type = c(
          "continuous_cloglog",
          "MESS",
          "clamp"
        ),
        source_ASC = c(
          normalizePath(projected_ASC, winslash = "/", mustWork = TRUE),
          if (is.na(mess_ASC)) NA_character_ else normalizePath(
            mess_ASC,
            winslash = "/",
            mustWork = TRUE
          ),
          if (is.na(clamp_ASC)) NA_character_ else normalizePath(
            clamp_ASC,
            winslash = "/",
            mustWork = TRUE
          )
        ),
        output_GeoTIFF = c(
          normalizePath(output_tif, winslash = "/", mustWork = TRUE),
          NA_character_,
          NA_character_
        ),
        stringsAsFactors = FALSE
      )

      convert_auxiliary <- function(input_ASC, suffix) {
        if (is.na(input_ASC)) {
          return(NA_character_)
        }

        auxiliary <- terra::rast(input_ASC)

        if (is.na(terra::crs(auxiliary)) || terra::crs(auxiliary) == "") {
          terra::crs(auxiliary) <- terra::crs(
            future_inputs$reference_future_raster
          )
        }

        auxiliary_tif <- file.path(
          scenario_output_dir,
          paste0(
            species_name, "_",
            gcm_name, "_",
            scenario_folder,
            "_", suffix, ".tif"
          )
        )

        terra::writeRaster(
          auxiliary,
          auxiliary_tif,
          overwrite = TRUE,
          wopt = list(gdal = "COMPRESS=LZW")
        )

        normalizePath(auxiliary_tif, winslash = "/", mustWork = TRUE)
      }

      optional_outputs$output_GeoTIFF[
        optional_outputs$output_type == "MESS"
      ] <- convert_auxiliary(mess_ASC, "MESS")

      optional_outputs$output_GeoTIFF[
        optional_outputs$output_type == "clamp"
      ] <- convert_auxiliary(clamp_ASC, "clamp")

      write.csv(
        optional_outputs,
        file.path(metadata_dir, "scenario_output_inventory.csv"),
        row.names = FALSE
      )

      # Keep only required provenance files if raw MaxEnt folders are not wanted.
      if (!keep_raw_maxent_outputs) {
        preserve_files <- c(
          "maxentResults.csv",
          paste0(species_name, ".lambdas"),
          "maxent_command.txt"
        )

        preservation_dir <- file.path(
          metadata_dir,
          "preserved_MaxEnt_provenance"
        )

        dir.create(preservation_dir, recursive = TRUE, showWarnings = FALSE)

        for (file_name in preserve_files) {
          source_file <- file.path(raw_maxent_output_dir, file_name)

          if (file.exists(source_file)) {
            file.copy(
              source_file,
              file.path(preservation_dir, file_name),
              overwrite = TRUE
            )
          }
        }

        safe_remove(raw_maxent_output_dir)
      }

      if (delete_temporary_future_ASC_after_success) {
        safe_remove(future_inputs$temporary_ASC_directory)
      }

      list(
        status = "completed",
        output_tif = normalizePath(output_tif, winslash = "/", mustWork = TRUE),
        projected_ASC = normalizePath(
          projected_ASC,
          winslash = "/",
          mustWork = TRUE
        ),
        bioc_file = normalizePath(
          future_inputs$bioc_file,
          winslash = "/",
          mustWork = TRUE
        ),
        elevation_file = normalizePath(
          future_inputs$elevation_file,
          winslash = "/",
          mustWork = TRUE
        ),
        min_suitability = projection_info$min_suitability,
        max_suitability = projection_info$max_suitability,
        valid_prediction_cells = projection_info$n_non_NA,
        message = paste0(
          "MaxEnt result: ",
          normalizePath(maxent_results_csv, winslash = "/")
        )
      )
    },
    error = function(error) {
      # Preserve temporary files after an error for diagnosis.
      list(
        status = "failed",
        output_tif = NA_character_,
        projected_ASC = NA_character_,
        bioc_file = NA_character_,
        elevation_file = NA_character_,
        min_suitability = NA_real_,
        max_suitability = NA_real_,
        valid_prediction_cells = NA_real_,
        message = conditionMessage(error)
      )
    }
  )

  elapsed_minutes <- as.numeric(
    difftime(Sys.time(), scenario_start, units = "mins")
  )

  status_records[[length(status_records) + 1]] <- data.frame(
    gcm = gcm_name,
    ssp = scenario$ssp,
    period_start = scenario$period_start,
    period_end = scenario$period_end,
    scenario_folder = scenario_folder,
    scenario_id = scenario_id,
    input_directory = scenario$input_directory,
    status = result$status,
    output_continuous_tif = result$output_tif,
    projected_ASC = result$projected_ASC,
    source_bioc_tif = result$bioc_file,
    source_elv_tif = result$elevation_file,
    min_suitability = result$min_suitability,
    max_suitability = result$max_suitability,
    valid_prediction_cells = result$valid_prediction_cells,
    elapsed_minutes = elapsed_minutes,
    message = result$message,
    stringsAsFactors = FALSE
  )

  status_table <- do.call(rbind, status_records)

  write_status(status_table, manifest_csv)

  if (identical(result$status, "completed")) {
    cat(
      "  Status: completed in ",
      format(round(elapsed_minutes, 2), nsmall = 2),
      " minutes.\n",
      sep = ""
    )
  } else {
    cat(
      "  Status: FAILED. Read log or manifest message before running Script 12.\n"
    )

    if (stop_immediately_after_first_failure) {
      stop(
        "\nProjection stopped after the first failed scenario to protect runtime.\n",
        "Read the scenario log and manifest before re-running.",
        call. = FALSE
      )
    }
  }
}

# ---------------------------------------------------------------------
# 8. FINAL REPORT
# ---------------------------------------------------------------------

status_table <- do.call(rbind, status_records)

if (is.null(status_table)) {
  stop("No projection records were created.", call. = FALSE)
}

write_status(status_table, manifest_csv)

completed_n <- sum(status_table$status == "completed")
skipped_n <- sum(status_table$status == "skipped_existing_output")
failed_n <- sum(status_table$status == "failed")

elapsed_total_hours <- as.numeric(
  difftime(Sys.time(), start_all, units = "hours")
)

summary_table <- data.frame(
  expected_scenarios = expected_total,
  detected_scenarios = nrow(future_scenarios),
  completed_this_run = completed_n,
  skipped_existing = skipped_n,
  failed = failed_n,
  total_completed_or_existing = completed_n + skipped_n,
  total_elapsed_hours = elapsed_total_hours,
  stringsAsFactors = FALSE
)

summary_csv <- file.path(
  manifest_root,
  "future_projection_run_summary.csv"
)

write.csv(summary_table, summary_csv, row.names = FALSE)

cat("\n======================================================================\n")
cat("FUTURE MAXENT PROJECTION FINISHED — Camellia japonica\n")
cat("======================================================================\n")
cat("Completed this run: ", completed_n, "\n", sep = "")
cat("Skipped existing outputs: ", skipped_n, "\n", sep = "")
cat("Failed scenarios: ", failed_n, "\n", sep = "")
cat(
  "Total completed/existing: ",
  completed_n + skipped_n,
  " of ",
  nrow(future_scenarios),
  "\n",
  sep = ""
)
cat(
  "Total runtime (hours): ",
  format(round(elapsed_total_hours, 2), nsmall = 2),
  "\n",
  sep = ""
)
cat("\nStatus manifest:\n", normalizePath(manifest_csv, winslash = "/"), "\n", sep = "")
cat("Run summary:\n", normalizePath(summary_csv, winslash = "/"), "\n", sep = "")

if (failed_n > 0) {
  cat(
    "\nAt least one scenario failed. Do NOT proceed to GCM ensemble until ",
    "the failed rows in future_projection_status_manifest.csv are resolved.\n"
  )
} else {
  cat(
    "\nAll detected future scenarios have continuous GCM-specific suitability ",
    "rasters. The next step is Script 12: three-GCM continuous ensemble.\n"
  )
}
