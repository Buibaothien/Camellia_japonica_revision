# =============================================================================
# 16A_NATIVE_VS_EXTRA_NATIVE_PCAENV_CAMELLIA_JAPONICA_V3_NATIVE_LAYER_FIXED.R
# =============================================================================
# PURPOSE
# -------
# Robust native-versus-extra-native environmental niche assessment for
# Camellia japonica, based on the FINAL 828 occurrence candidates and the
# five selected model predictors:
#
#   Bio2, Bio6, Bio14, Bio18, Elv
#
# This script addresses reviewer concerns that occurrences outside the native
# range may represent cultivated, managed, naturalized, or otherwise
# non-natural populations. It does NOT claim that a point inside the native
# range is necessarily wild; the native/extra-native label is a spatial
# membership classification based solely on the supplied native-range polygon.
#
# ANALYSES
# --------
# 1. Spatially classify the 828 final records as within vs outside the exact
#    native-range polygon used in the preceding workflow.
# 2. Extract five-predictor environmental values at all records.
# 3. Compare univariate distributions of the five predictors.
# 4. Conduct PCA-env niche comparison.
# 5. Calculate Schoener's D and Hellinger's I.
# 6. Run niche equivalency and bidirectional niche-similarity tests
#    (1,000 permutations by default).
# 7. Export figures, results, R objects, logs, and an explicit recommendation
#    for whether Script 16B (native-only MaxEnt sensitivity model) is needed.
#
# IMPORTANT
# ---------
# This script will NOT silently finish. It writes:
#   00_run_status/16A_STARTED.txt immediately,
#   00_run_status/16A_COMPLETED.txt when successful, or
#   00_run_status/16A_FAILED.txt if any error occurs.
#
# Before running:
# - Keep native_range_vector <- NA to allow safe automatic discovery.
# - If automatic discovery reports more than one candidate, set the exact path
#   to the earlier native-range/Admin-1 polygon below.
# - Do NOT replace it with broad country boundaries.
# =============================================================================

# =============================================================================
# A. STARTUP STATUS — PRINTED BEFORE PACKAGES ARE LOADED
# =============================================================================

options(warn = 1)

cat("\n")
cat("======================================================================\n")
cat("STARTING SCRIPT 16A — NATIVE VS EXTRA-NATIVE PCA-ENV\n")
cat("Species: Camellia japonica\n")
cat("======================================================================\n")

script_started_at <- Sys.time()

# =============================================================================
# B. USER SETTINGS
# =============================================================================

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# ---------------------------------------------------------------------
# Exact native-range vector used in the previous native/Admin-1 workflow.
# ---------------------------------------------------------------------
# Leave as NA to use safe automatic discovery.
#
# If the script stops because it finds zero or multiple possible candidates,
# open:
#   .../16A_native_vs_extra_native_PCAenv_v2/00_setup/
#       native_range_vector_candidates.csv
#
# Then paste the EXACT native-range vector path here, e.g.
# native_range_vector <- "E:/PREDICTION TREE/Camellia_japonica/.../native_range_admin1_40units.shp"
#
native_range_vector <- NA_character_

# Exact layer in the native-range GeoPackage used by the prior workflow.
# "native_wild_range_union" is the dissolved native wild distribution polygon.
native_range_layer <- "native_wild_range_union"

# TRUE: search the project folder for a unique vector with clear native-range
# naming. FALSE: require native_range_vector to be set manually.
auto_discover_native_range_vector <- TRUE

# ---------------------------------------------------------------------
# Environmental data and test settings
# ---------------------------------------------------------------------

selected_predictors <- c("Bio2", "Bio6", "Bio14", "Bio18", "Elv")

# Availability samples used to calibrate PCA-env (native + extra-native).
background_samples_per_context <- 20000L

# Standard robust number for randomization tests.
n_permutations <- 1000L

# PCA-env grid resolution.
pca_grid_resolution <- 100L

# Random seed used for all stochastic operations.
random_seed <- 20260707L

# Lower bound for each occurrence group after spatial classification.
minimum_occurrences_per_group <- 15L

# Quantile used by equivalency/similarity niche-dynamics summaries.
intersection_quantile <- 0.10

# =====================================================================
# C. PROJECT PATHS
# =====================================================================

# Final 828 occurrence file used in the final MaxEnt workflow.
occurrence_candidates <- c(
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

# Candidate folders for the five selected current predictors. The script uses
# the first folder that contains ALL five required layers.
predictor_directory_candidates <- c(
  file.path(
    project_dir,
    "model_inputs",
    "Primary_Global_ClimateTopo_Current"
  ),
  file.path(
    project_dir,
    "global_reanalysis_outputs",
    "09_kuenm_global_OR10",
    "01_global_variables",
    "set1_primary5"
  )
)

output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "16A_native_vs_extra_native_PCAenv_v2"
)

status_dir <- file.path(output_root, "00_run_status")
setup_dir <- file.path(output_root, "00_setup")
occurrence_dir <- file.path(output_root, "01_occurrence_classification")
environment_dir <- file.path(output_root, "02_environmental_data")
result_dir <- file.path(output_root, "03_niche_test_results")
figure_dir <- file.path(output_root, "04_figures")
object_dir <- file.path(output_root, "05_R_objects")
log_dir <- file.path(output_root, "06_logs")

started_marker <- file.path(status_dir, "16A_STARTED.txt")
completed_marker <- file.path(status_dir, "16A_COMPLETED.txt")
failed_marker <- file.path(status_dir, "16A_FAILED.txt")
execution_log <- file.path(log_dir, "16A_execution_log.txt")

vector_candidates_csv <- file.path(
  setup_dir,
  "native_range_vector_candidates.csv"
)

settings_csv <- file.path(setup_dir, "analysis_settings.csv")
native_vector_provenance_csv <- file.path(
  setup_dir,
  "native_range_vector_provenance.csv"
)

occurrence_status_csv <- file.path(
  occurrence_dir,
  "Camellia_japonica_final_828_native_status.csv"
)

within_native_csv <- file.path(
  occurrence_dir,
  "Camellia_japonica_within_native_range.csv"
)

outside_native_csv <- file.path(
  occurrence_dir,
  "Camellia_japonica_outside_native_range.csv"
)

occurrence_count_csv <- file.path(
  occurrence_dir,
  "occurrence_counts_by_native_status.csv"
)

occurrence_environment_csv <- file.path(
  environment_dir,
  "occurrence_environmental_values_5_predictors.csv"
)

native_background_csv <- file.path(
  environment_dir,
  "native_environmental_availability_sample.csv"
)

outside_background_csv <- file.path(
  environment_dir,
  "outside_native_environmental_availability_sample.csv"
)

univariate_comparison_csv <- file.path(
  result_dir,
  "five_predictor_univariate_comparison.csv"
)

pca_variance_csv <- file.path(result_dir, "PCA_variance.csv")
pca_loadings_csv <- file.path(result_dir, "PCA_loadings.csv")
overlap_csv <- file.path(result_dir, "niche_overlap_metrics.csv")
test_results_csv <- file.path(
  result_dir,
  "niche_equivalency_similarity_tests.csv"
)

null_distribution_csv <- file.path(
  result_dir,
  "niche_test_null_distributions.csv"
)

decision_note_txt <- file.path(
  result_dir,
  "decision_note_for_native_only_sensitivity_model.txt"
)

map_png <- file.path(
  figure_dir,
  "Fig_S1_native_vs_extra_native_occurrences.png"
)

pca_png <- file.path(
  figure_dir,
  "Fig_S2_PCAenv_occurrence_space.png"
)

niche_dynamics_pdf <- file.path(
  figure_dir,
  "Fig_S3_PCAenv_niche_dynamics.pdf"
)

equivalency_D_pdf <- file.path(
  figure_dir,
  "Fig_S4_equivalency_test_Schoener_D.pdf"
)

equivalency_I_pdf <- file.path(
  figure_dir,
  "Fig_S5_equivalency_test_Hellinger_I.pdf"
)

# =============================================================================
# D. HELPER FUNCTIONS
# =============================================================================

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

normalise_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}

find_column <- function(data, candidates, object_label) {
  raw_names <- names(data)
  norm_names <- normalise_name(raw_names)
  norm_candidates <- normalise_name(candidates)

  direct <- match(norm_candidates, norm_names, nomatch = 0)
  direct <- direct[direct > 0]

  if (length(direct) > 0) {
    return(raw_names[direct[1]])
  }

  for (candidate in norm_candidates) {
    matches <- grep(candidate, norm_names, fixed = TRUE)
    if (length(matches) > 0) return(raw_names[matches[1]])
  }

  stop(
    "\nCannot identify a required column in ", object_label, ".\n",
    "Available columns: ", paste(raw_names, collapse = ", "),
    call. = FALSE
  )
}

write_log <- function(...) {
  message <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", message)
  cat(line, "\n")
  cat(line, "\n", file = execution_log, append = TRUE)
}

assert_file_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop("\nMissing ", label, ":\n", path, call. = FALSE)
  }
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
      "\nPredictor-grid mismatch for ", label, ".\n",
      "All five predictor rasters must have identical CRS, extent, resolution, ",
      "rows, and columns.",
      call. = FALSE
    )
  }
}

safe_number <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

find_first_existing_file <- function(candidates, label) {
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0) {
    stop(
      "\nNo ", label, " was found. Checked:\n",
      paste(candidates, collapse = "\n"),
      call. = FALSE
    )
  }

  normalizePath(existing[1], winslash = "/", mustWork = TRUE)
}

detect_predictor_set <- function(candidate_directories, variables) {
  candidates <- list()

  for (directory in candidate_directories) {
    if (!dir.exists(directory)) next

    files <- list.files(
      directory,
      pattern = "\\.(tif|tiff|asc)$",
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )

    if (length(files) == 0) next

    stem <- tolower(tools::file_path_sans_ext(basename(files)))
    selected_files <- character(length(variables))

    for (i in seq_along(variables)) {
      selected_files[i] <- files[
        match(tolower(variables[i]), stem, nomatch = 0)
      ][1]
    }

    if (all(!is.na(selected_files)) && all(file.exists(selected_files))) {
      candidates[[length(candidates) + 1]] <- list(
        directory = normalizePath(directory, winslash = "/", mustWork = TRUE),
        files = normalizePath(selected_files, winslash = "/", mustWork = TRUE)
      )
    }
  }

  if (length(candidates) == 0) {
    stop(
      "\nCould not find a folder containing ALL five predictors:\n",
      paste(selected_predictors, collapse = ", "),
      "\nChecked folders:\n",
      paste(candidate_directories, collapse = "\n"),
      call. = FALSE
    )
  }

  # Prefer the first candidate directory in the listed priority order.
  candidates[[1]]
}

find_native_range_vector <- function(
  project_directory,
  manually_set_path,
  manually_set_layer,
  allow_auto_discovery,
  candidate_csv
) {
  if (!is.na(manually_set_path) && nzchar(manually_set_path)) {
    assert_file_exists(manually_set_path, "manually specified native-range vector")

    return(list(
      path = normalizePath(manually_set_path, winslash = "/", mustWork = TRUE),
      layer = manually_set_layer,
      method = "manual path"
    ))
  }

  if (!allow_auto_discovery) {
    stop(
      "\nNo native-range vector was specified and auto-discovery is disabled.",
      call. = FALSE
    )
  }

  all_vectors <- list.files(
    project_directory,
    pattern = "\\.(shp|gpkg|geojson|json)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(all_vectors) == 0) {
    stop(
      "\nNo vector files (*.shp, *.gpkg, *.geojson) were found below:\n",
      project_directory,
      call. = FALSE
    )
  }

  norm_paths <- normalise_name(all_vectors)

  # Strict keyword rule: requires "native" plus one range-related word.
  likely_native <- grepl("native", norm_paths) &
    grepl("range|admin|wild|distribution", norm_paths)

  # Extra conservative fallback for an explicitly named Camellia native file.
  likely_native <- likely_native |
    (
      grepl("camellia", norm_paths) &
        grepl("native", norm_paths)
    )

  candidate_paths <- all_vectors[likely_native]

  inventory <- data.frame(
    vector_path = normalizePath(all_vectors, winslash = "/", mustWork = TRUE),
    file_name = basename(all_vectors),
    selected_by_native_keyword_rule = all_vectors %in% candidate_paths,
    stringsAsFactors = FALSE
  )

  write.csv(inventory, candidate_csv, row.names = FALSE)

  if (length(candidate_paths) == 1) {
    return(list(
      path = normalizePath(candidate_paths[1], winslash = "/", mustWork = TRUE),
      layer = manually_set_layer,
      method = "automatic unique keyword match"
    ))
  }

  if (length(candidate_paths) == 0) {
    stop(
      "\nNo uniquely named native-range vector was found automatically.\n",
      "Candidate inventory written to:\n",
      candidate_csv,
      "\n\nSet native_range_vector manually to the exact native-range/Admin-1 ",
      "polygon used in your prior workflow, then rerun.",
      call. = FALSE
    )
  }

  stop(
    "\nMore than one possible native-range vector was found.\n",
    "Candidate inventory written to:\n",
    candidate_csv,
    "\n\nSet native_range_vector manually to the exact native-range/Admin-1 ",
    "polygon used in your prior workflow, then rerun.",
    call. = FALSE
  )
}

read_polygon_vector <- function(vector_info) {
  extension <- tolower(tools::file_ext(vector_info$path))

  if (extension == "gpkg") {
    layers <- sf::st_layers(vector_info$path)$name
    selected_layer <- vector_info$layer

    if (is.na(selected_layer) || !nzchar(selected_layer)) {
      if (length(layers) != 1) {
        stop(
          "\nThe selected GeoPackage contains more than one layer:\n",
          paste(layers, collapse = ", "),
          "\nSet native_range_layer manually.",
          call. = FALSE
        )
      }
      selected_layer <- layers[1]
    }

    vector <- sf::st_read(vector_info$path, layer = selected_layer, quiet = TRUE)
  } else {
    selected_layer <- NA_character_
    vector <- sf::st_read(vector_info$path, quiet = TRUE)
  }

  if (is.na(sf::st_crs(vector))) {
    stop(
      "\nNative-range vector has no CRS:\n",
      vector_info$path,
      call. = FALSE
    )
  }

  vector <- sf::st_make_valid(vector)

  geom_type <- unique(as.character(sf::st_geometry_type(vector)))
  if (!all(geom_type %in% c("POLYGON", "MULTIPOLYGON"))) {
    stop(
      "\nNative-range vector must contain polygon geometry.\n",
      "Detected: ", paste(geom_type, collapse = ", "),
      call. = FALSE
    )
  }

  list(vector = vector, layer = selected_layer)
}

sample_available_environment <- function(environment_stack, sample_size, label) {
  sampled <- terra::spatSample(
    environment_stack,
    size = sample_size,
    method = "random",
    as.df = TRUE,
    values = TRUE,
    na.rm = TRUE
  )

  sampled <- as.data.frame(sampled)

  # terra can sometimes append coordinate columns under certain versions.
  sampled <- sampled[, intersect(selected_predictors, names(sampled)), drop = FALSE]

  missing <- setdiff(selected_predictors, names(sampled))
  if (length(missing) > 0) {
    stop(
      "\nEnvironmental sampling for ", label, " lacks: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  sampled <- sampled[
    complete.cases(sampled[, selected_predictors, drop = FALSE]),
    selected_predictors,
    drop = FALSE
  ]

  if (nrow(sampled) < 1000L) {
    stop(
      "\nOnly ", nrow(sampled), " complete environmental cells were sampled for ",
      label, ". At least 1,000 are required.",
      call. = FALSE
    )
  }

  sampled
}

get_suprow_scores <- function(pca_object, new_values) {
  result <- ade4::suprow(pca_object, new_values)

  if (!is.null(result$lisup)) {
    scores <- as.data.frame(result$lisup)
  } else if (!is.null(result$li)) {
    scores <- as.data.frame(result$li)
  } else {
    stop("ade4::suprow did not return PCA row scores.", call. = FALSE)
  }

  scores <- scores[, 1:2, drop = FALSE]
  names(scores) <- c("PC1", "PC2")
  scores
}

get_test_scalar <- function(test_object, field, metric) {
  value <- test_object[[field]]

  if (is.null(value)) return(NA_real_)

  if (is.matrix(value) || is.data.frame(value)) {
    if (metric %in% colnames(value)) return(as.numeric(value[1, metric]))
    if (metric %in% rownames(value)) return(as.numeric(value[metric, 1]))
    return(NA_real_)
  }

  if (!is.null(names(value)) && metric %in% names(value)) {
    return(as.numeric(value[[metric]]))
  }

  if (length(value) == 1 && metric %in% c("D", "I")) {
    return(as.numeric(value))
  }

  NA_real_
}

get_test_pvalue <- function(test_object, metric) {
  candidates <- c(
    paste0("p.", metric),
    paste0("p_", metric),
    paste0("p", metric)
  )

  for (candidate in candidates) {
    if (!is.null(test_object[[candidate]])) {
      return(as.numeric(test_object[[candidate]])[1])
    }
  }

  NA_real_
}

write_null_distribution <- function(test_object, test_name) {
  sim <- test_object$sim

  if (is.null(sim)) return(NULL)

  sim <- as.data.frame(sim)

  if (!"D" %in% names(sim) && ncol(sim) >= 1) names(sim)[1] <- "D"
  if (!"I" %in% names(sim) && ncol(sim) >= 2) names(sim)[2] <- "I"

  sim$test <- test_name
  sim$replicate <- seq_len(nrow(sim))
  sim
}

safe_plot <- function(plot_function, output_path, description) {
  tryCatch(
    {
      plot_function()
      write_log("Figure created: ", description)
      TRUE
    },
    error = function(e) {
      # Close any graphics device left open by a failed plotting call.
      if (grDevices::dev.cur() > 1) {
        try(grDevices::dev.off(), silent = TRUE)
      }

      write_log(
        "Figure skipped because plotting failed [",
        description,
        "]: ",
        conditionMessage(e)
      )
      FALSE
    }
  )
}

# =============================================================================
# E. ANALYSIS FUNCTION
# =============================================================================

run_native_vs_extra_native_pcaenv <- function() {

  # -------------------------------------------------------------------
  # E1. Create folders and startup evidence
  # -------------------------------------------------------------------

  for (folder in c(
    output_root, status_dir, setup_dir, occurrence_dir, environment_dir,
    result_dir, figure_dir, object_dir, log_dir
  )) {
    ensure_dir(folder)
  }

  unlink(completed_marker, force = TRUE)
  unlink(failed_marker, force = TRUE)

  writeLines(
    c(
      "Script 16A started",
      paste("Started:", format(script_started_at, "%Y-%m-%d %H:%M:%S")),
      paste("Project:", project_dir),
      "This file proves that the script entered execution."
    ),
    started_marker
  )

  write_log("Entered run_native_vs_extra_native_pcaenv().")

  # -------------------------------------------------------------------
  # E2. Packages
  # -------------------------------------------------------------------

  required_packages <- c("terra", "sf", "ade4", "ecospat", "ggplot2")

  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_packages) > 0) {
    write_log(
      "Installing missing package(s): ",
      paste(missing_packages, collapse = ", ")
    )

    install.packages(
      missing_packages,
      repos = "https://cloud.r-project.org",
      dependencies = TRUE
    )
  }

  for (package_name in required_packages) {
    suppressPackageStartupMessages(
      library(package_name, character.only = TRUE)
    )
  }

  write_log(
    "Loaded packages: ",
    paste(required_packages, collapse = ", ")
  )

  # -------------------------------------------------------------------
  # E3. Input discovery
  # -------------------------------------------------------------------

  occurrence_csv <- find_first_existing_file(
    occurrence_candidates,
    "final 828 occurrence CSV"
  )

  predictor_info <- detect_predictor_set(
    predictor_directory_candidates,
    selected_predictors
  )

  vector_info <- find_native_range_vector(
    project_directory = project_dir,
    manually_set_path = native_range_vector,
    manually_set_layer = native_range_layer,
    allow_auto_discovery = auto_discover_native_range_vector,
    candidate_csv = vector_candidates_csv
  )

  native_vector_info <- read_polygon_vector(vector_info)

  write_log("Final occurrence file: ", occurrence_csv)
  write_log("Predictor folder: ", predictor_info$directory)
  write_log("Native-range vector: ", vector_info$path)
  write_log("Native-range selection method: ", vector_info$method)

  write.csv(
    data.frame(
      setting = c(
        "species",
        "occurrence_csv",
        "expected_occurrence_count",
        "predictor_directory",
        "predictor_files",
        "native_range_vector",
        "native_range_layer",
        "native_range_selection_method",
        "background_samples_per_context",
        "pca_grid_resolution",
        "n_permutations",
        "intersection_quantile",
        "random_seed"
      ),
      value = c(
        species_name,
        occurrence_csv,
        828,
        predictor_info$directory,
        paste(predictor_info$files, collapse = "; "),
        vector_info$path,
        native_vector_info$layer,
        vector_info$method,
        background_samples_per_context,
        pca_grid_resolution,
        n_permutations,
        intersection_quantile,
        random_seed
      ),
      stringsAsFactors = FALSE
    ),
    settings_csv,
    row.names = FALSE
  )

  write.csv(
    data.frame(
      vector_path = vector_info$path,
      layer = native_vector_info$layer,
      selection_method = vector_info$method,
      feature_count = nrow(native_vector_info$vector),
      CRS = sf::st_crs(native_vector_info$vector)$input,
      interpretation = "Spatial membership only; it does not independently verify wild, cultivated, or naturalized status.",
      stringsAsFactors = FALSE
    ),
    native_vector_provenance_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E4. Read final occurrence data and classify spatial status
  # -------------------------------------------------------------------

  occurrences <- read.csv(
    occurrence_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  lon_col <- find_column(
    occurrences,
    c("longitude", "decimalLongitude", "lon", "x"),
    "final occurrence CSV"
  )

  lat_col <- find_column(
    occurrences,
    c("latitude", "decimalLatitude", "lat", "y"),
    "final occurrence CSV"
  )

  occurrences$longitude_analysis <- safe_number(occurrences[[lon_col]])
  occurrences$latitude_analysis <- safe_number(occurrences[[lat_col]])

  occurrences <- occurrences[
    is.finite(occurrences$longitude_analysis) &
      is.finite(occurrences$latitude_analysis) &
      occurrences$longitude_analysis >= -180 &
      occurrences$longitude_analysis <= 180 &
      occurrences$latitude_analysis >= -90 &
      occurrences$latitude_analysis <= 90,
    ,
    drop = FALSE
  ]

  # Preserve all original attributes but remove exact duplicate coordinate pairs.
  occurrences <- occurrences[
    !duplicated(
      paste(
        occurrences$longitude_analysis,
        occurrences$latitude_analysis,
        sep = "_"
      )
    ),
    ,
    drop = FALSE
  ]

  if (nrow(occurrences) != 828L) {
    stop(
      "\nThe final occurrence input should contain 828 unique valid records, ",
      "but it contains ", nrow(occurrences), ".\n",
      "Use the exact file employed by the completed final model.",
      call. = FALSE
    )
  }

  occurrence_sf <- sf::st_as_sf(
    occurrences,
    coords = c("longitude_analysis", "latitude_analysis"),
    crs = 4326,
    remove = FALSE
  )

  native_range_wgs84 <- sf::st_transform(native_vector_info$vector, 4326)
  native_union <- sf::st_union(native_range_wgs84)

  inside_native <- lengths(
    sf::st_intersects(occurrence_sf, native_union)
  ) > 0

  occurrences$native_range_status <- ifelse(
    inside_native,
    "within_native_range",
    "outside_native_range"
  )

  occurrences$native_range_status <- factor(
    occurrences$native_range_status,
    levels = c("within_native_range", "outside_native_range")
  )

  occurrence_status <- sf::st_drop_geometry(
    cbind(occurrence_sf, native_range_status = occurrences$native_range_status)
  )

  write.csv(occurrence_status, occurrence_status_csv, row.names = FALSE)

  within_native <- occurrence_status[
    occurrence_status$native_range_status == "within_native_range",
    ,
    drop = FALSE
  ]

  outside_native <- occurrence_status[
    occurrence_status$native_range_status == "outside_native_range",
    ,
    drop = FALSE
  ]

  write.csv(within_native, within_native_csv, row.names = FALSE)
  write.csv(outside_native, outside_native_csv, row.names = FALSE)

  status_count <- data.frame(
    native_range_status = c("within_native_range", "outside_native_range"),
    n_occurrences = c(nrow(within_native), nrow(outside_native)),
    percent_of_final_occurrences = 100 * c(
      nrow(within_native), nrow(outside_native)
    ) / nrow(occurrence_status),
    stringsAsFactors = FALSE
  )

  write.csv(status_count, occurrence_count_csv, row.names = FALSE)

  write_log(
    "Occurrence groups: within native = ", nrow(within_native),
    "; outside native = ", nrow(outside_native)
  )

  if (
    nrow(within_native) < minimum_occurrences_per_group ||
    nrow(outside_native) < minimum_occurrences_per_group
  ) {
    stop(
      "\nAt least one group has fewer than ", minimum_occurrences_per_group,
      " occurrence records.\nWithin native = ", nrow(within_native),
      "\nOutside native = ", nrow(outside_native),
      "\nVerify the exact native-range polygon.",
      call. = FALSE
    )
  }

  # -------------------------------------------------------------------
  # E5. Read five selected predictors and extract values at occurrences
  # -------------------------------------------------------------------

  predictor_stack <- terra::rast(predictor_info$files)
  names(predictor_stack) <- selected_predictors

  for (i in seq_len(terra::nlyr(predictor_stack))) {
    assert_same_grid(
      predictor_stack[[1]],
      predictor_stack[[i]],
      selected_predictors[i]
    )
  }

  occurrence_vect <- terra::vect(
    occurrence_status,
    geom = c("longitude_analysis", "latitude_analysis"),
    crs = "EPSG:4326"
  )

  extracted <- terra::extract(
    predictor_stack,
    occurrence_vect,
    ID = FALSE
  )

  extracted <- as.data.frame(extracted)

  missing_predictors <- setdiff(selected_predictors, names(extracted))
  if (length(missing_predictors) > 0) {
    stop(
      "\nPredictor extraction failed. Missing: ",
      paste(missing_predictors, collapse = ", "),
      call. = FALSE
    )
  }

  occurrence_environment <- cbind(
    occurrence_status,
    extracted[, selected_predictors, drop = FALSE]
  )

  occurrence_environment <- occurrence_environment[
    complete.cases(
      occurrence_environment[, selected_predictors, drop = FALSE]
    ),
    ,
    drop = FALSE
  ]

  if (nrow(occurrence_environment) != 828L) {
    stop(
      "\nNot all 828 final occurrences have valid values for all five predictors.\n",
      "Retained after extraction: ", nrow(occurrence_environment),
      call. = FALSE
    )
  }

  write.csv(
    occurrence_environment,
    occurrence_environment_csv,
    row.names = FALSE
  )

  native_env <- occurrence_environment[
    occurrence_environment$native_range_status == "within_native_range",
    selected_predictors,
    drop = FALSE
  ]

  extra_env <- occurrence_environment[
    occurrence_environment$native_range_status == "outside_native_range",
    selected_predictors,
    drop = FALSE
  ]

  # -------------------------------------------------------------------
  # E6. Native and extra-native environmental availability
  # -------------------------------------------------------------------

  set.seed(random_seed)

  template <- predictor_stack[[1]]

  native_mask <- terra::rasterize(
    terra::vect(native_range_wgs84),
    template,
    field = 1,
    background = 0,
    touches = TRUE
  )

  native_stack <- terra::mask(
    predictor_stack,
    native_mask,
    maskvalues = 0
  )

  outside_mask <- terra::ifel(native_mask == 0, 1, NA)
  outside_stack <- terra::mask(predictor_stack, outside_mask)

  write_log("Sampling native-range environmental availability.")
  native_background <- sample_available_environment(
    native_stack,
    background_samples_per_context,
    "native range"
  )

  write_log("Sampling outside-native environmental availability.")
  outside_background <- sample_available_environment(
    outside_stack,
    background_samples_per_context,
    "outside native range"
  )

  native_background$availability_context <- "within_native_range"
  outside_background$availability_context <- "outside_native_range"

  write.csv(native_background, native_background_csv, row.names = FALSE)
  write.csv(outside_background, outside_background_csv, row.names = FALSE)

  native_bg_values <- native_background[, selected_predictors, drop = FALSE]
  outside_bg_values <- outside_background[, selected_predictors, drop = FALSE]

  # -------------------------------------------------------------------
  # E7. Univariate environmental comparison
  # -------------------------------------------------------------------

  univariate <- do.call(
    rbind,
    lapply(selected_predictors, function(variable_name) {
      a <- native_env[[variable_name]]
      b <- extra_env[[variable_name]]

      test <- stats::wilcox.test(a, b, exact = FALSE)

      pooled_sd <- sqrt(
        (
          (length(a) - 1) * stats::var(a) +
          (length(b) - 1) * stats::var(b)
        ) /
        (length(a) + length(b) - 2)
      )

      effect_size <- if (is.finite(pooled_sd) && pooled_sd > 0) {
        (mean(b) - mean(a)) / pooled_sd
      } else {
        NA_real_
      }

      data.frame(
        predictor = variable_name,
        within_native_n = length(a),
        within_native_mean = mean(a),
        within_native_sd = stats::sd(a),
        within_native_median = stats::median(a),
        outside_native_n = length(b),
        outside_native_mean = mean(b),
        outside_native_sd = stats::sd(b),
        outside_native_median = stats::median(b),
        standardised_difference_outside_minus_within = effect_size,
        Wilcoxon_W = unname(test$statistic),
        Wilcoxon_p_raw = test$p.value,
        stringsAsFactors = FALSE
      )
    })
  )

  univariate$Wilcoxon_p_BH <- stats::p.adjust(
    univariate$Wilcoxon_p_raw,
    method = "BH"
  )

  write.csv(univariate, univariate_comparison_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # E8. PCA-env
  # -------------------------------------------------------------------

  write_log("Running PCA-env calibration.")

  # PCA calibration space = combined environmental availability in both
  # geographic contexts. It is centered and scaled before ordination.
  pca_background <- rbind(native_bg_values, outside_bg_values)

  pca <- ade4::dudi.pca(
    pca_background,
    center = TRUE,
    scale = TRUE,
    scannf = FALSE,
    nf = 2
  )

  pca_variance <- data.frame(
    axis = paste0("PC", seq_along(pca$eig)),
    eigenvalue = as.numeric(pca$eig),
    percent_variance = 100 * pca$eig / sum(pca$eig),
    cumulative_percent_variance = cumsum(100 * pca$eig / sum(pca$eig)),
    stringsAsFactors = FALSE
  )

  pca_loadings <- data.frame(
    predictor = rownames(pca$c1),
    PC1_loading = as.numeric(pca$c1[, 1]),
    PC2_loading = as.numeric(pca$c1[, 2]),
    stringsAsFactors = FALSE
  )

  write.csv(pca_variance, pca_variance_csv, row.names = FALSE)
  write.csv(pca_loadings, pca_loadings_csv, row.names = FALSE)

  score_global <- as.data.frame(pca$li[, 1:2, drop = FALSE])
  names(score_global) <- c("PC1", "PC2")

  score_native_bg <- get_suprow_scores(pca, native_bg_values)
  score_outside_bg <- get_suprow_scores(pca, outside_bg_values)
  score_native_occ <- get_suprow_scores(pca, native_env)
  score_outside_occ <- get_suprow_scores(pca, extra_env)

  z_native <- ecospat::ecospat.grid.clim.dyn(
    glob = score_global,
    glob1 = score_native_bg,
    sp = score_native_occ,
    R = pca_grid_resolution,
    th.sp = 0,
    th.env = 0
  )

  z_outside <- ecospat::ecospat.grid.clim.dyn(
    glob = score_global,
    glob1 = score_outside_bg,
    sp = score_outside_occ,
    R = pca_grid_resolution,
    th.sp = 0,
    th.env = 0
  )

  overlap <- ecospat::ecospat.niche.overlap(
    z_native,
    z_outside,
    cor = TRUE
  )

  # Equivalency with lower alternative tests whether observed overlap is
  # significantly LOWER than random occurrence reassignment (niche divergence).
  write_log("Running equivalency test: ", n_permutations, " permutations.")

  equivalency <- ecospat::ecospat.niche.equivalency.test(
    z1 = z_native,
    z2 = z_outside,
    rep = n_permutations,
    intersection = intersection_quantile,
    overlap.alternative = "lower",
    expansion.alternative = "higher",
    stability.alternative = "lower",
    unfilling.alternative = "higher"
  )

  # Directional similarity tests. A small p-value for D under the "higher"
  # alternative means overlap is higher than expected from random shifts in
  # the relevant available environmental space.
  write_log("Running similarity test: outside-native relative to native.")

  similarity_outside_to_native <- ecospat::ecospat.niche.similarity.test(
    z1 = z_native,
    z2 = z_outside,
    rep = n_permutations,
    intersection = intersection_quantile,
    rand.type = 2,
    ncores = 1,
    overlap.alternative = "higher",
    expansion.alternative = "lower",
    stability.alternative = "higher",
    unfilling.alternative = "lower"
  )

  write_log("Running similarity test: native relative to outside-native.")

  similarity_native_to_outside <- ecospat::ecospat.niche.similarity.test(
    z1 = z_outside,
    z2 = z_native,
    rep = n_permutations,
    intersection = intersection_quantile,
    rand.type = 2,
    ncores = 1,
    overlap.alternative = "higher",
    expansion.alternative = "lower",
    stability.alternative = "higher",
    unfilling.alternative = "lower"
  )

  # -------------------------------------------------------------------
  # E9. Results tables and serialized objects
  # -------------------------------------------------------------------

  overlap_table <- data.frame(
    metric = c("Schoener_D", "Hellinger_I"),
    value = c(as.numeric(overlap$D), as.numeric(overlap$I)),
    note = "Availability-corrected overlap calculated in PCA-env space",
    stringsAsFactors = FALSE
  )

  write.csv(overlap_table, overlap_csv, row.names = FALSE)

  test_table <- rbind(
    data.frame(
      test = "Equivalency",
      direction = "Within native vs outside native",
      observed_D = get_test_scalar(equivalency, "obs", "D"),
      observed_I = get_test_scalar(equivalency, "obs", "I"),
      p_D = get_test_pvalue(equivalency, "D"),
      p_I = get_test_pvalue(equivalency, "I"),
      overlap_alternative = "lower (tests divergence)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      test = "Similarity",
      direction = "Outside native relative to native",
      observed_D = get_test_scalar(similarity_outside_to_native, "obs", "D"),
      observed_I = get_test_scalar(similarity_outside_to_native, "obs", "I"),
      p_D = get_test_pvalue(similarity_outside_to_native, "D"),
      p_I = get_test_pvalue(similarity_outside_to_native, "I"),
      overlap_alternative = "higher (tests conservatism)",
      stringsAsFactors = FALSE
    ),
    data.frame(
      test = "Similarity",
      direction = "Native relative to outside native",
      observed_D = get_test_scalar(similarity_native_to_outside, "obs", "D"),
      observed_I = get_test_scalar(similarity_native_to_outside, "obs", "I"),
      p_D = get_test_pvalue(similarity_native_to_outside, "D"),
      p_I = get_test_pvalue(similarity_native_to_outside, "I"),
      overlap_alternative = "higher (tests conservatism)",
      stringsAsFactors = FALSE
    )
  )

  write.csv(test_table, test_results_csv, row.names = FALSE)

  null_data <- do.call(
    rbind,
    Filter(
      Negate(is.null),
      list(
        write_null_distribution(equivalency, "Equivalency"),
        write_null_distribution(
          similarity_outside_to_native,
          "Similarity_outside_relative_to_native"
        ),
        write_null_distribution(
          similarity_native_to_outside,
          "Similarity_native_relative_to_outside"
        )
      )
    )
  )

  if (!is.null(null_data)) {
    write.csv(null_data, null_distribution_csv, row.names = FALSE)
  }

  saveRDS(pca, file.path(object_dir, "PCAenv_model.rds"))
  saveRDS(z_native, file.path(object_dir, "PCAenv_grid_native.rds"))
  saveRDS(z_outside, file.path(object_dir, "PCAenv_grid_outside_native.rds"))
  saveRDS(equivalency, file.path(object_dir, "equivalency_test.rds"))
  saveRDS(
    similarity_outside_to_native,
    file.path(object_dir, "similarity_outside_to_native.rds")
  )
  saveRDS(
    similarity_native_to_outside,
    file.path(object_dir, "similarity_native_to_outside.rds")
  )

  # -------------------------------------------------------------------
  # E10. Figures
  # -------------------------------------------------------------------

  write_log("Creating figures.")

  occurrence_sf <- sf::st_as_sf(
    occurrence_environment,
    coords = c("longitude_analysis", "latitude_analysis"),
    crs = 4326,
    remove = FALSE
  )

  map_plot <- ggplot2::ggplot() +
    ggplot2::geom_sf(
      data = native_range_wgs84,
      fill = "grey92",
      colour = "grey50",
      linewidth = 0.2
    ) +
    ggplot2::geom_sf(
      data = occurrence_sf[
        occurrence_sf$native_range_status == "outside_native_range",
        ,
        drop = FALSE
      ],
      colour = "#D55E00",
      size = 1.0,
      alpha = 0.75
    ) +
    ggplot2::geom_sf(
      data = occurrence_sf[
        occurrence_sf$native_range_status == "within_native_range",
        ,
        drop = FALSE
      ],
      colour = "#0072B2",
      size = 1.0,
      alpha = 0.80
    ) +
    ggplot2::coord_sf(
      xlim = c(-180, 180),
      ylim = c(-60, 85),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = expression(italic("Camellia japonica") ~ "final occurrence candidates"),
      subtitle = "Blue: within supplied native range; orange: outside supplied native range",
      x = "Longitude",
      y = "Latitude"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    map_png,
    map_plot,
    width = 12,
    height = 6.5,
    dpi = 400
  )

  set.seed(random_seed)

  max_background_plot <- min(10000L, nrow(score_global))

  plot_data <- rbind(
    data.frame(
      PC1 = score_global[
        sample(seq_len(nrow(score_global)), max_background_plot),
        "PC1"
      ],
      PC2 = score_global[
        sample(seq_len(nrow(score_global)), max_background_plot),
        "PC2"
      ],
      group = "Environmental availability",
      stringsAsFactors = FALSE
    ),
    data.frame(
      PC1 = score_native_occ$PC1,
      PC2 = score_native_occ$PC2,
      group = "Within native range",
      stringsAsFactors = FALSE
    ),
    data.frame(
      PC1 = score_outside_occ$PC1,
      PC2 = score_outside_occ$PC2,
      group = "Outside native range",
      stringsAsFactors = FALSE
    )
  )

  pc1_var <- round(pca_variance$percent_variance[1], 1)
  pc2_var <- round(pca_variance$percent_variance[2], 1)

  pca_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = PC1, y = PC2)
  ) +
    ggplot2::geom_point(
      data = subset(plot_data, group == "Environmental availability"),
      colour = "grey75",
      size = 0.35,
      alpha = 0.28
    ) +
    ggplot2::geom_point(
      data = subset(plot_data, group == "Outside native range"),
      colour = "#D55E00",
      size = 1.2,
      alpha = 0.75
    ) +
    ggplot2::geom_point(
      data = subset(plot_data, group == "Within native range"),
      colour = "#0072B2",
      size = 1.2,
      alpha = 0.80
    ) +
    ggplot2::labs(
      title = "PCA-env of final occurrence candidates",
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    pca_png,
    pca_plot,
    width = 8.5,
    height = 6.5,
    dpi = 400
  )

  safe_plot(function() {
    grDevices::pdf(niche_dynamics_pdf, width = 9, height = 8)
    on.exit(grDevices::dev.off(), add = TRUE)

    ecospat::ecospat.plot.niche.dyn(
      z_native,
      z_outside,
      interest = 2,
      title = "PCA-env niche dynamics: within native vs outside native range",
      name.axis1 = paste0("PC1 (", pc1_var, "%)"),
      name.axis2 = paste0("PC2 (", pc2_var, "%)")
    )
  }, niche_dynamics_pdf, "PCA-env niche-dynamics PDF")

  safe_plot(function() {
    grDevices::pdf(equivalency_D_pdf, width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)

    ecospat::ecospat.plot.overlap.test(
      equivalency,
      type = "D",
      title = "Equivalency test: Schoener's D"
    )
  }, equivalency_D_pdf, "equivalency Schoener D PDF")

  safe_plot(function() {
    grDevices::pdf(equivalency_I_pdf, width = 8, height = 6)
    on.exit(grDevices::dev.off(), add = TRUE)

    ecospat::ecospat.plot.overlap.test(
      equivalency,
      type = "I",
      title = "Equivalency test: Hellinger's I"
    )
  }, equivalency_I_pdf, "equivalency Hellinger I PDF")

  # -------------------------------------------------------------------
  # E11. Decision note
  # -------------------------------------------------------------------

  equivalency_p_D <- get_test_pvalue(equivalency, "D")
  equivalency_p_I <- get_test_pvalue(equivalency, "I")

  has_divergence <- any(
    c(equivalency_p_D, equivalency_p_I) < 0.05,
    na.rm = TRUE
  )

  recommendation <- if (has_divergence) {
    paste(
      "The lower-tail equivalency test indicates observed niche overlap that is",
      "significantly lower than expected under random reassignment for at least",
      "one overlap metric. Run Script 16B: native-only MaxEnt sensitivity model.",
      "Do not treat the 828-record global model as the sole biogeographic",
      "inference until the native-only sensitivity model is compared."
    )
  } else {
    paste(
      "The lower-tail equivalency test did not detect statistically supported",
      "niche divergence at alpha = 0.05. The 828-record model can remain the",
      "primary climatic-permissiveness model, subject to explicit manuscript",
      "language that outside-native records are not independently verified as",
      "wild populations and that this analysis is a spatial sensitivity check."
    )
  }

  writeLines(
    c(
      "NATIVE VS EXTRA-NATIVE PCA-ENV: DECISION NOTE",
      "================================================",
      "",
      paste0("Occurrences analysed: ", nrow(occurrence_environment)),
      paste0("Within supplied native range: ", nrow(native_env)),
      paste0("Outside supplied native range: ", nrow(extra_env)),
      paste0("Schoener's D: ", format(overlap$D, digits = 6)),
      paste0("Hellinger's I: ", format(overlap$I, digits = 6)),
      paste0("Equivalency lower-tail p(D): ", format(equivalency_p_D, digits = 6)),
      paste0("Equivalency lower-tail p(I): ", format(equivalency_p_I, digits = 6)),
      "",
      "Interpretation:",
      "- Spatial native status does not independently prove wild, cultivated, or naturalized status.",
      "- Lower-tail equivalency p < 0.05 means observed overlap is lower than random reallocation, supporting environmental divergence.",
      "- Directional similarity p < 0.05 with alternative = higher means overlap exceeds random shifts in that direction.",
      "",
      "Recommended next step:",
      recommendation
    ),
    decision_note_txt
  )

  # -------------------------------------------------------------------
  # E12. Completion
  # -------------------------------------------------------------------

  completed_text <- c(
    "SCRIPT 16A COMPLETED SUCCESSFULLY",
    paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste("Occurrences analysed:", nrow(occurrence_environment)),
    paste("Within native range:", nrow(native_env)),
    paste("Outside native range:", nrow(extra_env)),
    paste("Schoener D:", format(overlap$D, digits = 6)),
    paste("Hellinger I:", format(overlap$I, digits = 6)),
    paste("Equivalency p(D):", format(equivalency_p_D, digits = 6)),
    paste("Equivalency p(I):", format(equivalency_p_I, digits = 6)),
    paste("Decision note:", decision_note_txt)
  )

  writeLines(completed_text, completed_marker)

  write_log("Script completed successfully.")

  cat("\n======================================================================\n")
  cat("SCRIPT 16A COMPLETED — NATIVE VS EXTRA-NATIVE PCA-ENV\n")
  cat("======================================================================\n")
  cat("Occurrences analysed: ", nrow(occurrence_environment), "\n", sep = "")
  cat("Within supplied native range: ", nrow(native_env), "\n", sep = "")
  cat("Outside supplied native range: ", nrow(extra_env), "\n", sep = "")
  cat("Schoener's D: ", format(overlap$D, digits = 6), "\n", sep = "")
  cat("Hellinger's I: ", format(overlap$I, digits = 6), "\n", sep = "")
  cat(
    "Equivalency lower-tail p(D): ",
    format(equivalency_p_D, digits = 6),
    "\n",
    sep = ""
  )
  cat(
    "Equivalency lower-tail p(I): ",
    format(equivalency_p_I, digits = 6),
    "\n",
    sep = ""
  )
  cat("\nDecision note:\n", decision_note_txt, "\n", sep = "")
  cat("Execution log:\n", execution_log, "\n", sep = "")

  invisible(list(
    occurrence_status = occurrence_status,
    overlap = overlap,
    equivalency = equivalency,
    similarity_outside_to_native = similarity_outside_to_native,
    similarity_native_to_outside = similarity_native_to_outside
  ))
}

# =============================================================================
# F. EXECUTE WITH EXPLICIT FAILURE REPORTING
# =============================================================================

tryCatch(
  {
    run_native_vs_extra_native_pcaenv()
  },
  error = function(error) {

    # Ensure folders exist even when the error occurs early.
    dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    failure_message <- paste0(
      "SCRIPT 16A FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(error), "\n\n",
      "Read the execution log: ", execution_log, "\n"
    )

    writeLines(failure_message, failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 16A FAILED\n")
    cat("======================================================================\n")
    cat(failure_message, "\n")

    stop(error)
  }
)
