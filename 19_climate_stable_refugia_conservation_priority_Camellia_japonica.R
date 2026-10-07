# =============================================================================
# 19_CLIMATE_STABLE_REFUGIA_CONSERVATION_PRIORITY_CAMELLIA_JAPONICA.R
# =============================================================================
# PURPOSE
# -------
# Identify climate-stable refugia, in situ conservation priorities, and ex situ
# germplasm-collection priorities for Camellia japonica within the native-wild
# range, using the final uncertainty-filtered SDM workflow.
#
# PRIMARY RULE USED FOR MAIN RESULTS
# ----------------------------------
# Strict future suitability:
#   future cell is suitable only when:
#     1) all 3 GCMs predict suitability at the MTSPS threshold, and
#     2) no GCM shows MaxEnt clamping at that cell.
#
# Main conservation priority is restricted to:
#   native_wild_range_union
#
# This is intentionally conservative and responds to reviewer concerns about
# uncertainty, extrapolation, and the distinction between native-range
# conservation versus extra-native climatic permissiveness.
#
# PRIORITY DEFINITIONS IN THE NATIVE RANGE
# ----------------------------------------
# Future persistence = number of the 16 SSP × period scenarios in which a cell
# remains strict future suitable.
#
# Integrated strict priority map:
#   0 = No current/future conservation priority
#   1 = Emerging future suitable native area
#       current unsuitable + strict future persistence >= 12/16
#   2 = Ex situ germplasm collection priority
#       current suitable + strict future persistence <= 4/16
#   3 = Secondary monitoring / adaptive management zone
#       current suitable + strict future persistence 5-11/16
#   4 = High-priority in situ conservation zone
#       current suitable + strict future persistence 12-15/16
#   5 = Core climate-stable refugium
#       current suitable + strict future persistence 16/16
#
# Supplementary relaxed rule:
#   relaxed future suitable = 3/3 GCM suitable + clamping in no more than 1 GCM
#
# OUTPUTS
# -------
# 1. Native range mask aligned to the final raster grid
# 2. Current MTSPS suitability within native range
# 3. Strict and relaxed future-persistence rasters (0-16)
# 4. Integrated priority maps within native range
# 5. Binary masks for:
#      - core climate-stable refugia
#      - high-priority in situ conservation zones
#      - ex situ germplasm-collection priority zones
#      - emerging future suitable native areas
# 6. Area tables and a strict-versus-relaxed sensitivity summary
# =============================================================================

options(warn = 1)

cat("\n")
cat("======================================================================\n")
cat("STARTING SCRIPT 19 — CLIMATE-STABLE REFUGIA AND PRIORITY ZONES\n")
cat("Species: Camellia japonica\n")
cat("Primary future rule: 3/3 GCM suitable + no clamping\n")
cat("Native-range restriction: native_wild_range_union\n")
cat("======================================================================\n")

script_started_at <- Sys.time()

# =============================================================================
# A. USER SETTINGS
# =============================================================================

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

native_range_gpkg <- file.path(
  project_dir,
  "data",
  "native_wild_range_Camellia_japonica.gpkg"
)

native_range_layer <- "native_wild_range_union"

# Main thresholds for persistence categories.
core_persistence_n <- 16L
high_insitu_min_n <- 12L
exsitu_max_n <- 4L

# 5-11/16 is retained as a secondary monitoring/adaptive-management category.
secondary_min_n <- 5L
secondary_max_n <- 11L

# TRUE gives slightly more inclusive rasterization of polygon edges.
native_rasterize_touches <- TRUE

# FALSE: reuse existing output rasters if already written.
overwrite_existing_maps <- FALSE

# =============================================================================
# B. INPUT PATHS
# =============================================================================

final_model_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
)

script18_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "18_interGCM_uncertainty_agreement_extrapolation"
)

script18B_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "18B_clamping_quantification_extrapolation_confidence"
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

strict_common_domain_tif <- file.path(
  script18_root,
  "07_supporting_files",
  "common_valid_domain_current_plus_all_48_individual_GCM_projections.tif"
)

geodesic_cell_area_tif <- file.path(
  script18_root,
  "07_supporting_files",
  "global_geodesic_cell_area_km2.tif"
)

script18B_map_status_manifest_csv <- file.path(
  script18B_root,
  "08_supporting_files",
  "clamping_confidence_map_status_manifest.csv"
)

script18B_priority_readiness_csv <- file.path(
  script18B_root,
  "06_area_tables",
  "05_clamping_priority_readiness_summary.csv"
)

# =============================================================================
# C. OUTPUT PATHS
# =============================================================================

output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "19_climate_stable_refugia_conservation_priority"
)

status_dir <- file.path(output_root, "00_run_status")
setup_dir <- file.path(output_root, "01_setup")
native_mask_dir <- file.path(output_root, "02_native_range_mask")
current_dir <- file.path(output_root, "03_current_suitability_native")
persistence_dir <- file.path(output_root, "04_future_persistence_maps")
priority_map_dir <- file.path(output_root, "05_integrated_priority_maps")
binary_mask_dir <- file.path(output_root, "06_binary_priority_masks")
area_dir <- file.path(output_root, "07_area_tables")
figure_dir <- file.path(output_root, "08_figures")
support_dir <- file.path(output_root, "09_supporting_files")
log_dir <- file.path(output_root, "10_logs")

started_marker <- file.path(status_dir, "19_STARTED.txt")
completed_marker <- file.path(status_dir, "19_COMPLETED.txt")
failed_marker <- file.path(status_dir, "19_FAILED.txt")
execution_log <- file.path(log_dir, "19_execution_log.txt")

input_manifest_csv <- file.path(setup_dir, "01_input_manifest.csv")
future_mask_inventory_csv <- file.path(setup_dir, "02_future_mask_inventory.csv")
native_layer_inventory_csv <- file.path(setup_dir, "03_native_range_layer_inventory.csv")

priority_legend_csv <- file.path(support_dir, "integrated_priority_legend.csv")
persistence_category_legend_csv <- file.path(support_dir, "future_persistence_category_legend.csv")
binary_mask_legend_csv <- file.path(support_dir, "binary_priority_mask_legend.csv")
map_status_manifest_csv <- file.path(support_dir, "priority_map_status_manifest.csv")
run_summary_csv <- file.path(support_dir, "script19_run_summary.csv")

native_mask_tif <- file.path(
  native_mask_dir,
  paste0(species_name, "_native_wild_range_mask_strict_common_domain.tif")
)

current_binary_global_tif <- file.path(
  current_dir,
  paste0(species_name, "_current_MTSPS_suitable_strict_common_domain.tif")
)

current_binary_native_tif <- file.path(
  current_dir,
  paste0(species_name, "_current_MTSPS_suitable_native_range.tif")
)

strict_persistence_global_tif <- file.path(
  persistence_dir,
  paste0(species_name, "_strict_future_persistence_0_16_global.tif")
)

strict_persistence_native_tif <- file.path(
  persistence_dir,
  paste0(species_name, "_strict_future_persistence_0_16_native_range.tif")
)

relaxed_persistence_global_tif <- file.path(
  persistence_dir,
  paste0(species_name, "_relaxed_future_persistence_0_16_global.tif")
)

relaxed_persistence_native_tif <- file.path(
  persistence_dir,
  paste0(species_name, "_relaxed_future_persistence_0_16_native_range.tif")
)

strict_persistence_category_native_tif <- file.path(
  persistence_dir,
  paste0(species_name, "_strict_future_persistence_category_native_range.tif")
)

relaxed_persistence_category_native_tif <- file.path(
  persistence_dir,
  paste0(species_name, "_relaxed_future_persistence_category_native_range.tif")
)

strict_priority_native_tif <- file.path(
  priority_map_dir,
  paste0(species_name, "_strict_integrated_conservation_priority_native_range.tif")
)

relaxed_priority_native_tif <- file.path(
  priority_map_dir,
  paste0(species_name, "_relaxed_integrated_conservation_priority_native_range.tif")
)

strict_core_refugia_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_strict_core_climate_stable_refugia_native_range.tif")
)

strict_high_insitu_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_strict_high_priority_in_situ_12plus_native_range.tif")
)

strict_exsitu_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_strict_ex_situ_germplasm_priority_native_range.tif")
)

strict_emerging_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_strict_emerging_future_suitable_native_range.tif")
)

relaxed_core_refugia_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_relaxed_core_climate_stable_refugia_native_range.tif")
)

relaxed_high_insitu_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_relaxed_high_priority_in_situ_12plus_native_range.tif")
)

relaxed_exsitu_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_relaxed_ex_situ_germplasm_priority_native_range.tif")
)

relaxed_emerging_tif <- file.path(
  binary_mask_dir,
  paste0(species_name, "_relaxed_emerging_future_suitable_native_range.tif")
)

native_summary_csv <- file.path(
  area_dir,
  "01_native_range_and_current_suitability_summary.csv"
)

strict_persistence_area_csv <- file.path(
  area_dir,
  "02_strict_future_persistence_area_native_range.csv"
)

relaxed_persistence_area_csv <- file.path(
  area_dir,
  "03_relaxed_future_persistence_area_native_range.csv"
)

strict_priority_area_csv <- file.path(
  area_dir,
  "04_strict_integrated_priority_area_native_range.csv"
)

relaxed_priority_area_csv <- file.path(
  area_dir,
  "05_relaxed_integrated_priority_area_native_range.csv"
)

binary_priority_area_csv <- file.path(
  area_dir,
  "06_binary_priority_mask_area_summary.csv"
)

strict_vs_relaxed_csv <- file.path(
  area_dir,
  "07_strict_vs_relaxed_priority_sensitivity_summary.csv"
)

manuscript_summary_csv <- file.path(
  area_dir,
  "08_manuscript_ready_priority_summary.csv"
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
    hits <- grep(candidate, normalized_names, fixed = TRUE)
    if (length(hits) > 0) return(raw_names[hits[1]])
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
      "All rasters used in Script 19 must align.",
      call. = FALSE
    )
  }
}

calculate_area_by_code <- function(class_raster, cell_area_raster, domain_raster,
                                   codes, class_names) {
  assert_same_grid(cell_area_raster, class_raster, "cell-area raster")
  assert_same_grid(domain_raster, class_raster, "analysis-domain raster")

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
  value <- area_table$area_km2[area_table$code == code]

  if (length(value) != 1L) {
    stop("Cannot retrieve unique area for code ", code, call. = FALSE)
  }

  value
}

write_raster <- function(raster_object, path, datatype = "INT1U") {
  terra::writeRaster(
    raster_object,
    path,
    overwrite = TRUE,
    datatype = datatype,
    wopt = list(gdal = "COMPRESS=LZW")
  )

  invisible(path)
}

create_persistence_category <- function(persistence_raster, native_mask) {
  # 0 = 0/16
  # 1 = 1-4/16
  # 2 = 5-8/16
  # 3 = 9-11/16
  # 4 = 12-15/16
  # 5 = 16/16
  category <- terra::ifel(is.na(native_mask), NA, 0)

  category <- terra::ifel(
    !is.na(native_mask) & persistence_raster >= 1 & persistence_raster <= 4,
    1,
    category
  )
  category <- terra::ifel(
    !is.na(native_mask) & persistence_raster >= 5 & persistence_raster <= 8,
    2,
    category
  )
  category <- terra::ifel(
    !is.na(native_mask) & persistence_raster >= 9 & persistence_raster <= 11,
    3,
    category
  )
  category <- terra::ifel(
    !is.na(native_mask) & persistence_raster >= 12 & persistence_raster <= 15,
    4,
    category
  )
  category <- terra::ifel(
    !is.na(native_mask) & persistence_raster == 16,
    5,
    category
  )

  names(category) <- "future_persistence_category"
  terra::mask(category, native_mask)
}

create_integrated_priority <- function(persistence_raster, current_native,
                                       native_mask) {
  # Codes:
  # 0 = no priority
  # 1 = emerging future suitable native area
  # 2 = ex situ germplasm collection priority
  # 3 = secondary monitoring / adaptive management
  # 4 = high-priority in situ conservation
  # 5 = core climate-stable refugium
  priority <- terra::ifel(is.na(native_mask), NA, 0)

  priority <- terra::ifel(
    !is.na(native_mask) &
      current_native == 0 &
      persistence_raster >= high_insitu_min_n,
    1,
    priority
  )

  priority <- terra::ifel(
    !is.na(native_mask) &
      current_native == 1 &
      persistence_raster <= exsitu_max_n,
    2,
    priority
  )

  priority <- terra::ifel(
    !is.na(native_mask) &
      current_native == 1 &
      persistence_raster >= secondary_min_n &
      persistence_raster <= secondary_max_n,
    3,
    priority
  )

  priority <- terra::ifel(
    !is.na(native_mask) &
      current_native == 1 &
      persistence_raster >= high_insitu_min_n &
      persistence_raster <= 15,
    4,
    priority
  )

  priority <- terra::ifel(
    !is.na(native_mask) &
      current_native == 1 &
      persistence_raster == core_persistence_n,
    5,
    priority
  )

  names(priority) <- "integrated_priority"
  terra::mask(priority, native_mask)
}

create_binary_priority_masks <- function(persistence_raster, current_native,
                                         native_mask) {
  core <- terra::ifel(
    is.na(native_mask),
    NA,
    terra::ifel(
      current_native == 1 & persistence_raster == core_persistence_n,
      1,
      0
    )
  )
  names(core) <- "core_climate_stable_refugia"

  high_insitu <- terra::ifel(
    is.na(native_mask),
    NA,
    terra::ifel(
      current_native == 1 & persistence_raster >= high_insitu_min_n,
      1,
      0
    )
  )
  names(high_insitu) <- "high_priority_in_situ_12plus"

  exsitu <- terra::ifel(
    is.na(native_mask),
    NA,
    terra::ifel(
      current_native == 1 & persistence_raster <= exsitu_max_n,
      1,
      0
    )
  )
  names(exsitu) <- "ex_situ_germplasm_priority"

  emerging <- terra::ifel(
    is.na(native_mask),
    NA,
    terra::ifel(
      current_native == 0 & persistence_raster >= high_insitu_min_n,
      1,
      0
    )
  )
  names(emerging) <- "emerging_future_suitable_native_area"

  list(
    core = terra::mask(core, native_mask),
    high_insitu = terra::mask(high_insitu, native_mask),
    exsitu = terra::mask(exsitu, native_mask),
    emerging = terra::mask(emerging, native_mask)
  )
}

calculate_binary_mask_area <- function(mask_raster, cell_area, native_mask,
                                       mask_id, rule, future_rule) {
  area <- calculate_area_by_code(
    class_raster = mask_raster,
    cell_area_raster = cell_area,
    domain_raster = native_mask,
    codes = c(0, 1),
    class_names = c("Not selected", "Selected")
  )

  selected_area <- get_code_area(area, 1)

  data.frame(
    future_rule = future_rule,
    mask_id = mask_id,
    rule = rule,
    selected_area_km2 = selected_area,
    selected_area_10_4_km2 = selected_area / 10000,
    stringsAsFactors = FALSE
  )
}

create_priority_figure <- function(priority_raster, persistence_raster,
                                   native_current, output_pdf,
                                   title_text) {
  grDevices::pdf(output_pdf, width = 12.5, height = 5.2, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(1, 3), mar = c(2.7, 2.7, 3.0, 1.0))

  plot(
    native_current,
    col = c("#F0F0F0", "#1B7837"),
    breaks = c(-0.5, 0.5, 1.5),
    plg = FALSE,
    axes = TRUE,
    main = "A. Current suitable"
  )
  legend(
    "bottomleft",
    legend = c("Unsuitable", "Suitable"),
    fill = c("#F0F0F0", "#1B7837"),
    title = "MTSPS",
    bty = "n",
    cex = 0.75
  )

  plot(
    persistence_raster,
    col = hcl.colors(17, "YlOrRd", rev = FALSE),
    breaks = seq(-0.5, 16.5, by = 1),
    plg = FALSE,
    axes = TRUE,
    main = "B. Future persistence"
  )
  legend(
    "bottomleft",
    legend = c("0", "4", "8", "12", "16"),
    title = "Scenarios",
    bty = "n",
    cex = 0.75
  )

  plot(
    priority_raster,
    col = c("#F0F0F0", "#A6CEE3", "#E31A1C", "#FDBF6F", "#33A02C", "#006D2C"),
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
    plg = FALSE,
    axes = TRUE,
    main = "C. Priority class"
  )
  legend(
    "bottomleft",
    legend = c(
      "0 No priority",
      "1 Emerging",
      "2 Ex situ",
      "3 Secondary",
      "4 High in situ",
      "5 Core refugium"
    ),
    fill = c("#F0F0F0", "#A6CEE3", "#E31A1C", "#FDBF6F", "#33A02C", "#006D2C"),
    bty = "n",
    cex = 0.62
  )

  mtext(title_text, side = 3, outer = FALSE, line = 0.6, font = 2, cex = 1.05)
}

# =============================================================================
# E. MAIN ANALYSIS
# =============================================================================

run_script19 <- function() {

  for (folder in c(
    output_root, status_dir, setup_dir, native_mask_dir, current_dir,
    persistence_dir, priority_map_dir, binary_mask_dir, area_dir, figure_dir,
    support_dir, log_dir
  )) {
    ensure_dir(folder)
  }

  unlink(completed_marker, force = TRUE)
  unlink(failed_marker, force = TRUE)

  writeLines(
    c(
      "Script 19 started",
      paste("Started:", format(script_started_at, "%Y-%m-%d %H:%M:%S")),
      "Purpose: climate-stable refugia and conservation priority analysis"
    ),
    started_marker
  )

  if (!requireNamespace("terra", quietly = TRUE)) {
    install.packages("terra", repos = "https://cloud.r-project.org")
  }

  if (!requireNamespace("sf", quietly = TRUE)) {
    install.packages("sf", repos = "https://cloud.r-project.org")
  }

  suppressPackageStartupMessages(library(terra))
  suppressPackageStartupMessages(library(sf))

  write_log("Packages terra and sf loaded.")

  # -------------------------------------------------------------------
  # E1. Validate required inputs
  # -------------------------------------------------------------------

  current_continuous_path <- find_first_existing_file(
    current_continuous_candidates,
    "current final full-data continuous suitability raster"
  )

  required_files <- c(
    fixed_mtsps_csv,
    strict_common_domain_tif,
    geodesic_cell_area_tif,
    script18B_map_status_manifest_csv,
    script18B_priority_readiness_csv,
    native_range_gpkg
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
  strict_common_domain <- terra::rast(strict_common_domain_tif)
  cell_area <- terra::rast(geodesic_cell_area_tif)

  assert_same_grid(current_continuous, strict_common_domain, "strict common domain")
  assert_same_grid(current_continuous, cell_area, "geodesic cell area")

  map_manifest <- read.csv(
    script18B_map_status_manifest_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  readiness <- read.csv(
    script18B_priority_readiness_csv,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  write.csv(
    data.frame(
      input_type = c(
        "Current final full-data continuous suitability",
        "Fixed MTSPS threshold",
        "Strict common valid domain from Script 18",
        "Geodesic cell-area raster",
        "Script 18B map status manifest",
        "Script 18B priority-readiness summary",
        "Native range GPKG"
      ),
      path = c(
        current_continuous_path,
        normalizePath(fixed_mtsps_csv, winslash = "/", mustWork = TRUE),
        normalizePath(strict_common_domain_tif, winslash = "/", mustWork = TRUE),
        normalizePath(geodesic_cell_area_tif, winslash = "/", mustWork = TRUE),
        normalizePath(script18B_map_status_manifest_csv, winslash = "/", mustWork = TRUE),
        normalizePath(script18B_priority_readiness_csv, winslash = "/", mustWork = TRUE),
        normalizePath(native_range_gpkg, winslash = "/", mustWork = TRUE)
      ),
      stringsAsFactors = FALSE
    ),
    input_manifest_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E2. Read and rasterize native-wild range
  # -------------------------------------------------------------------

  native_layers <- sf::st_layers(native_range_gpkg)

  write.csv(
    data.frame(
      layer_name = native_layers$name,
      geometry_type = native_layers$geomtype,
      feature_count = native_layers$features,
      stringsAsFactors = FALSE
    ),
    native_layer_inventory_csv,
    row.names = FALSE
  )

  if (!native_range_layer %in% native_layers$name) {
    stop(
      "\nNative range layer not found in GPKG:\n",
      native_range_layer,
      "\nAvailable layers:\n",
      paste(native_layers$name, collapse = ", "),
      call. = FALSE
    )
  }

  native_sf <- sf::st_read(
    native_range_gpkg,
    layer = native_range_layer,
    quiet = TRUE
  )

  if (nrow(native_sf) < 1L) {
    stop("\nNative range layer contains no features.", call. = FALSE)
  }

  native_sf <- sf::st_make_valid(native_sf)
  native_vect <- terra::vect(native_sf)
  native_vect <- terra::project(native_vect, terra::crs(current_continuous))

  write_log("Rasterizing native range layer: ", native_range_layer)

  native_mask <- terra::rasterize(
    native_vect,
    current_continuous,
    field = 1,
    background = NA,
    touches = native_rasterize_touches
  )

  native_mask <- terra::mask(native_mask, strict_common_domain)
  names(native_mask) <- "native_wild_range_mask"

  write_raster(native_mask, native_mask_tif, datatype = "INT1U")

  native_area <- calculate_area_by_code(
    class_raster = native_mask,
    cell_area_raster = cell_area,
    domain_raster = native_mask,
    codes = 1,
    class_names = "Native wild range within strict common domain"
  )

  native_range_area_km2 <- get_code_area(native_area, 1)

  # -------------------------------------------------------------------
  # E3. Current MTSPS binary suitability
  # -------------------------------------------------------------------

  current_binary_global <- terra::ifel(
    is.na(strict_common_domain),
    NA,
    terra::ifel(current_continuous >= mtsps_threshold, 1, 0)
  )
  names(current_binary_global) <- "current_MTSPS_suitable"

  current_binary_native <- terra::mask(current_binary_global, native_mask)
  names(current_binary_native) <- "current_MTSPS_suitable_native"

  write_raster(current_binary_global, current_binary_global_tif, datatype = "INT1U")
  write_raster(current_binary_native, current_binary_native_tif, datatype = "INT1U")

  current_native_area <- calculate_area_by_code(
    class_raster = current_binary_native,
    cell_area_raster = cell_area,
    domain_raster = native_mask,
    codes = c(0, 1),
    class_names = c("Current unsuitable within native range", "Current suitable within native range")
  )

  current_unsuitable_native_km2 <- get_code_area(current_native_area, 0)
  current_suitable_native_km2 <- get_code_area(current_native_area, 1)

  # -------------------------------------------------------------------
  # E4. Resolve the 16 strict and 16 relaxed future masks
  # -------------------------------------------------------------------

  ssp_col <- find_column(map_manifest, c("ssp"), "Script 18B map status manifest")
  period_start_col <- find_column(map_manifest, c("period_start", "periodstart"), "Script 18B map status manifest")
  period_end_col <- find_column(map_manifest, c("period_end", "periodend"), "Script 18B map status manifest")
  strict_mask_col <- find_column(map_manifest, c("strict_unclamped_mask_tif"), "Script 18B map status manifest")
  relaxed_mask_col <- find_column(map_manifest, c("relaxed_mask_tif"), "Script 18B map status manifest")

  map_manifest$ssp_analysis <- as.character(map_manifest[[ssp_col]])
  map_manifest$period_start_analysis <- safe_number(map_manifest[[period_start_col]])
  map_manifest$period_end_analysis <- safe_number(map_manifest[[period_end_col]])
  map_manifest$strict_mask_path_analysis <- as.character(map_manifest[[strict_mask_col]])
  map_manifest$relaxed_mask_path_analysis <- as.character(map_manifest[[relaxed_mask_col]])

  map_manifest <- map_manifest[
    file.exists(map_manifest$strict_mask_path_analysis) &
      file.exists(map_manifest$relaxed_mask_path_analysis),
    ,
    drop = FALSE
  ]

  if (nrow(map_manifest) != 16L) {
    stop(
      "\nExpected 16 Script 18B future strict/relaxed mask records but found ",
      nrow(map_manifest), ".",
      call. = FALSE
    )
  }

  map_manifest <- map_manifest[
    order(map_manifest$ssp_analysis, map_manifest$period_start_analysis),
    ,
    drop = FALSE
  ]

  write.csv(
    map_manifest[
      ,
      c(
        "ssp_analysis",
        "period_start_analysis",
        "period_end_analysis",
        "strict_mask_path_analysis",
        "relaxed_mask_path_analysis"
      ),
      drop = FALSE
    ],
    future_mask_inventory_csv,
    row.names = FALSE
  )

  # Ensure all masks align with the final current grid.
  for (i in seq_len(nrow(map_manifest))) {
    strict_mask_i <- terra::rast(map_manifest$strict_mask_path_analysis[i])
    relaxed_mask_i <- terra::rast(map_manifest$relaxed_mask_path_analysis[i])

    scenario_label <- paste0(
      map_manifest$ssp_analysis[i],
      " ",
      map_manifest$period_start_analysis[i],
      "-",
      map_manifest$period_end_analysis[i]
    )

    assert_same_grid(
      current_continuous,
      strict_mask_i,
      paste0("strict future mask ", scenario_label)
    )
    assert_same_grid(
      current_continuous,
      relaxed_mask_i,
      paste0("relaxed future mask ", scenario_label)
    )
  }

  # -------------------------------------------------------------------
  # E5. Future persistence counts
  # -------------------------------------------------------------------

  write_log("Calculating strict and relaxed future-persistence rasters.")

  strict_stack <- terra::rast(map_manifest$strict_mask_path_analysis)
  relaxed_stack <- terra::rast(map_manifest$relaxed_mask_path_analysis)

  strict_persistence_global <- sum(strict_stack, na.rm = FALSE)
  strict_persistence_global <- terra::mask(strict_persistence_global, strict_common_domain)
  names(strict_persistence_global) <- "strict_future_persistence_0_16"

  relaxed_persistence_global <- sum(relaxed_stack, na.rm = FALSE)
  relaxed_persistence_global <- terra::mask(relaxed_persistence_global, strict_common_domain)
  names(relaxed_persistence_global) <- "relaxed_future_persistence_0_16"

  strict_persistence_native <- terra::mask(strict_persistence_global, native_mask)
  names(strict_persistence_native) <- "strict_future_persistence_0_16_native"

  relaxed_persistence_native <- terra::mask(relaxed_persistence_global, native_mask)
  names(relaxed_persistence_native) <- "relaxed_future_persistence_0_16_native"

  write_raster(strict_persistence_global, strict_persistence_global_tif, datatype = "INT1U")
  write_raster(strict_persistence_native, strict_persistence_native_tif, datatype = "INT1U")
  write_raster(relaxed_persistence_global, relaxed_persistence_global_tif, datatype = "INT1U")
  write_raster(relaxed_persistence_native, relaxed_persistence_native_tif, datatype = "INT1U")

  strict_persistence_category <- create_persistence_category(
    strict_persistence_native,
    native_mask
  )
  relaxed_persistence_category <- create_persistence_category(
    relaxed_persistence_native,
    native_mask
  )

  write_raster(
    strict_persistence_category,
    strict_persistence_category_native_tif,
    datatype = "INT1U"
  )
  write_raster(
    relaxed_persistence_category,
    relaxed_persistence_category_native_tif,
    datatype = "INT1U"
  )

  # -------------------------------------------------------------------
  # E6. Integrated priority maps and binary masks
  # -------------------------------------------------------------------

  write_log("Creating strict and relaxed integrated priority maps.")

  strict_priority <- create_integrated_priority(
    persistence_raster = strict_persistence_native,
    current_native = current_binary_native,
    native_mask = native_mask
  )

  relaxed_priority <- create_integrated_priority(
    persistence_raster = relaxed_persistence_native,
    current_native = current_binary_native,
    native_mask = native_mask
  )

  write_raster(strict_priority, strict_priority_native_tif, datatype = "INT1U")
  write_raster(relaxed_priority, relaxed_priority_native_tif, datatype = "INT1U")

  strict_masks <- create_binary_priority_masks(
    persistence_raster = strict_persistence_native,
    current_native = current_binary_native,
    native_mask = native_mask
  )

  relaxed_masks <- create_binary_priority_masks(
    persistence_raster = relaxed_persistence_native,
    current_native = current_binary_native,
    native_mask = native_mask
  )

  write_raster(strict_masks$core, strict_core_refugia_tif, datatype = "INT1U")
  write_raster(strict_masks$high_insitu, strict_high_insitu_tif, datatype = "INT1U")
  write_raster(strict_masks$exsitu, strict_exsitu_tif, datatype = "INT1U")
  write_raster(strict_masks$emerging, strict_emerging_tif, datatype = "INT1U")

  write_raster(relaxed_masks$core, relaxed_core_refugia_tif, datatype = "INT1U")
  write_raster(relaxed_masks$high_insitu, relaxed_high_insitu_tif, datatype = "INT1U")
  write_raster(relaxed_masks$exsitu, relaxed_exsitu_tif, datatype = "INT1U")
  write_raster(relaxed_masks$emerging, relaxed_emerging_tif, datatype = "INT1U")

  # -------------------------------------------------------------------
  # E7. Legends
  # -------------------------------------------------------------------

  write.csv(
    data.frame(
      priority_code = 0:5,
      priority_class = c(
        "No current/future conservation priority",
        "Emerging future suitable native area",
        "Ex situ germplasm collection priority",
        "Secondary monitoring / adaptive management zone",
        "High-priority in situ conservation zone",
        "Core climate-stable refugium"
      ),
      rule = c(
        "Native-range cell not selected by the other rules",
        "Current MTSPS unsuitable + future persistence >= 12/16",
        "Current MTSPS suitable + future persistence <= 4/16",
        "Current MTSPS suitable + future persistence 5-11/16",
        "Current MTSPS suitable + future persistence 12-15/16",
        "Current MTSPS suitable + future persistence 16/16"
      ),
      stringsAsFactors = FALSE
    ),
    priority_legend_csv,
    row.names = FALSE
  )

  write.csv(
    data.frame(
      persistence_category = 0:5,
      class_name = c(
        "0/16 scenarios",
        "1-4/16 scenarios",
        "5-8/16 scenarios",
        "9-11/16 scenarios",
        "12-15/16 scenarios",
        "16/16 scenarios"
      ),
      interpretation = c(
        "No future persistence",
        "Low future persistence",
        "Moderate-low future persistence",
        "Moderate-high future persistence",
        "High future persistence",
        "Full future persistence / stable across all scenarios"
      ),
      stringsAsFactors = FALSE
    ),
    persistence_category_legend_csv,
    row.names = FALSE
  )

  write.csv(
    data.frame(
      mask_id = c(
        "core_refugia",
        "high_insitu_12plus",
        "exsitu_germplasm",
        "emerging_future_suitable"
      ),
      rule = c(
        "Current suitable + persistence 16/16",
        "Current suitable + persistence >= 12/16",
        "Current suitable + persistence <= 4/16",
        "Current unsuitable + persistence >= 12/16"
      ),
      stringsAsFactors = FALSE
    ),
    binary_mask_legend_csv,
    row.names = FALSE
  )

  # -------------------------------------------------------------------
  # E8. Area tables
  # -------------------------------------------------------------------

  native_summary <- data.frame(
    species = species_name,
    MTSPS_threshold = mtsps_threshold,
    native_range_layer = native_range_layer,
    native_range_area_within_strict_common_domain_km2 = native_range_area_km2,
    native_range_area_10_4_km2 = native_range_area_km2 / 10000,
    current_MTSPS_unsuitable_native_km2 = current_unsuitable_native_km2,
    current_MTSPS_suitable_native_km2 = current_suitable_native_km2,
    current_MTSPS_suitable_native_10_4_km2 = current_suitable_native_km2 / 10000,
    current_percent_native_range_suitable =
      100 * current_suitable_native_km2 / native_range_area_km2,
    stringsAsFactors = FALSE
  )

  write.csv(native_summary, native_summary_csv, row.names = FALSE)

  strict_persistence_area <- calculate_area_by_code(
    class_raster = strict_persistence_category,
    cell_area_raster = cell_area,
    domain_raster = native_mask,
    codes = 0:5,
    class_names = c(
      "0/16 scenarios",
      "1-4/16 scenarios",
      "5-8/16 scenarios",
      "9-11/16 scenarios",
      "12-15/16 scenarios",
      "16/16 scenarios"
    )
  )

  strict_persistence_area$future_rule <- "strict_3of3_GCM_suitable_no_clamping"
  strict_persistence_area$percent_native_range <-
    100 * strict_persistence_area$area_km2 / native_range_area_km2

  relaxed_persistence_area <- calculate_area_by_code(
    class_raster = relaxed_persistence_category,
    cell_area_raster = cell_area,
    domain_raster = native_mask,
    codes = 0:5,
    class_names = c(
      "0/16 scenarios",
      "1-4/16 scenarios",
      "5-8/16 scenarios",
      "9-11/16 scenarios",
      "12-15/16 scenarios",
      "16/16 scenarios"
    )
  )

  relaxed_persistence_area$future_rule <- "relaxed_3of3_GCM_suitable_clamping_at_most_1"
  relaxed_persistence_area$percent_native_range <-
    100 * relaxed_persistence_area$area_km2 / native_range_area_km2

  write.csv(strict_persistence_area, strict_persistence_area_csv, row.names = FALSE)
  write.csv(relaxed_persistence_area, relaxed_persistence_area_csv, row.names = FALSE)

  strict_priority_area <- calculate_area_by_code(
    class_raster = strict_priority,
    cell_area_raster = cell_area,
    domain_raster = native_mask,
    codes = 0:5,
    class_names = c(
      "No current/future conservation priority",
      "Emerging future suitable native area",
      "Ex situ germplasm collection priority",
      "Secondary monitoring / adaptive management zone",
      "High-priority in situ conservation zone",
      "Core climate-stable refugium"
    )
  )

  strict_priority_area$future_rule <- "strict_3of3_GCM_suitable_no_clamping"
  strict_priority_area$percent_native_range <-
    100 * strict_priority_area$area_km2 / native_range_area_km2

  relaxed_priority_area <- calculate_area_by_code(
    class_raster = relaxed_priority,
    cell_area_raster = cell_area,
    domain_raster = native_mask,
    codes = 0:5,
    class_names = c(
      "No current/future conservation priority",
      "Emerging future suitable native area",
      "Ex situ germplasm collection priority",
      "Secondary monitoring / adaptive management zone",
      "High-priority in situ conservation zone",
      "Core climate-stable refugium"
    )
  )

  relaxed_priority_area$future_rule <- "relaxed_3of3_GCM_suitable_clamping_at_most_1"
  relaxed_priority_area$percent_native_range <-
    100 * relaxed_priority_area$area_km2 / native_range_area_km2

  write.csv(strict_priority_area, strict_priority_area_csv, row.names = FALSE)
  write.csv(relaxed_priority_area, relaxed_priority_area_csv, row.names = FALSE)

  binary_area_summary <- rbind(
    calculate_binary_mask_area(
      strict_masks$core,
      cell_area,
      native_mask,
      "core_refugia",
      "Current suitable + strict persistence 16/16",
      "strict_3of3_GCM_suitable_no_clamping"
    ),
    calculate_binary_mask_area(
      strict_masks$high_insitu,
      cell_area,
      native_mask,
      "high_insitu_12plus",
      "Current suitable + strict persistence >= 12/16",
      "strict_3of3_GCM_suitable_no_clamping"
    ),
    calculate_binary_mask_area(
      strict_masks$exsitu,
      cell_area,
      native_mask,
      "exsitu_germplasm",
      "Current suitable + strict persistence <= 4/16",
      "strict_3of3_GCM_suitable_no_clamping"
    ),
    calculate_binary_mask_area(
      strict_masks$emerging,
      cell_area,
      native_mask,
      "emerging_future_suitable",
      "Current unsuitable + strict persistence >= 12/16",
      "strict_3of3_GCM_suitable_no_clamping"
    ),
    calculate_binary_mask_area(
      relaxed_masks$core,
      cell_area,
      native_mask,
      "core_refugia",
      "Current suitable + relaxed persistence 16/16",
      "relaxed_3of3_GCM_suitable_clamping_at_most_1"
    ),
    calculate_binary_mask_area(
      relaxed_masks$high_insitu,
      cell_area,
      native_mask,
      "high_insitu_12plus",
      "Current suitable + relaxed persistence >= 12/16",
      "relaxed_3of3_GCM_suitable_clamping_at_most_1"
    ),
    calculate_binary_mask_area(
      relaxed_masks$exsitu,
      cell_area,
      native_mask,
      "exsitu_germplasm",
      "Current suitable + relaxed persistence <= 4/16",
      "relaxed_3of3_GCM_suitable_clamping_at_most_1"
    ),
    calculate_binary_mask_area(
      relaxed_masks$emerging,
      cell_area,
      native_mask,
      "Current unsuitable + relaxed persistence >= 12/16",
      "emerging_future_suitable",
      "relaxed_3of3_GCM_suitable_clamping_at_most_1"
    )
  )

  binary_area_summary$percent_native_range <-
    100 * binary_area_summary$selected_area_km2 / native_range_area_km2

  write.csv(binary_area_summary, binary_priority_area_csv, row.names = FALSE)

  strict_vs_relaxed <- merge(
    binary_area_summary[
      binary_area_summary$future_rule == "strict_3of3_GCM_suitable_no_clamping",
      c("mask_id", "selected_area_km2")
    ],
    binary_area_summary[
      binary_area_summary$future_rule == "relaxed_3of3_GCM_suitable_clamping_at_most_1",
      c("mask_id", "selected_area_km2")
    ],
    by = "mask_id",
    suffixes = c("_strict", "_relaxed"),
    all = TRUE
  )

  strict_vs_relaxed$area_difference_relaxed_minus_strict_km2 <-
    strict_vs_relaxed$selected_area_km2_relaxed -
    strict_vs_relaxed$selected_area_km2_strict

  strict_vs_relaxed$percent_difference_relative_to_strict <- ifelse(
    strict_vs_relaxed$selected_area_km2_strict > 0,
    100 *
      strict_vs_relaxed$area_difference_relaxed_minus_strict_km2 /
      strict_vs_relaxed$selected_area_km2_strict,
    NA_real_
  )

  write.csv(strict_vs_relaxed, strict_vs_relaxed_csv, row.names = FALSE)

  # Manuscript-ready concise summary.
  get_binary_area <- function(future_rule, mask_id) {
    value <- binary_area_summary$selected_area_km2[
      binary_area_summary$future_rule == future_rule &
        binary_area_summary$mask_id == mask_id
    ]

    if (length(value) != 1L) return(NA_real_)
    value
  }

  manuscript_summary <- data.frame(
    metric = c(
      "Native range area within strict common domain",
      "Current MTSPS suitable native range",
      "Strict core climate-stable refugia",
      "Strict high-priority in situ zones including core refugia",
      "Strict ex situ germplasm collection priority",
      "Strict emerging future suitable native areas",
      "Relaxed core climate-stable refugia",
      "Relaxed high-priority in situ zones including core refugia",
      "Relaxed ex situ germplasm collection priority",
      "Relaxed emerging future suitable native areas"
    ),
    area_km2 = c(
      native_range_area_km2,
      current_suitable_native_km2,
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "core_refugia"),
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "high_insitu_12plus"),
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "exsitu_germplasm"),
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "emerging_future_suitable"),
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "core_refugia"),
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "high_insitu_12plus"),
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "exsitu_germplasm"),
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "emerging_future_suitable")
    ),
    stringsAsFactors = FALSE
  )

  manuscript_summary$area_10_4_km2 <- manuscript_summary$area_km2 / 10000
  manuscript_summary$percent_native_range <- 100 *
    manuscript_summary$area_km2 / native_range_area_km2

  write.csv(manuscript_summary, manuscript_summary_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # E9. Status manifest and optional figures
  # -------------------------------------------------------------------

  map_status <- data.frame(
    map_id = c(
      "native_mask",
      "current_binary_global",
      "current_binary_native",
      "strict_persistence_global",
      "strict_persistence_native",
      "relaxed_persistence_global",
      "relaxed_persistence_native",
      "strict_priority_native",
      "relaxed_priority_native",
      "strict_core_refugia",
      "strict_high_insitu",
      "strict_exsitu",
      "strict_emerging",
      "relaxed_core_refugia",
      "relaxed_high_insitu",
      "relaxed_exsitu",
      "relaxed_emerging"
    ),
    path = c(
      native_mask_tif,
      current_binary_global_tif,
      current_binary_native_tif,
      strict_persistence_global_tif,
      strict_persistence_native_tif,
      relaxed_persistence_global_tif,
      relaxed_persistence_native_tif,
      strict_priority_native_tif,
      relaxed_priority_native_tif,
      strict_core_refugia_tif,
      strict_high_insitu_tif,
      strict_exsitu_tif,
      strict_emerging_tif,
      relaxed_core_refugia_tif,
      relaxed_high_insitu_tif,
      relaxed_exsitu_tif,
      relaxed_emerging_tif
    ),
    status = "completed",
    stringsAsFactors = FALSE
  )

  map_status$path <- normalizePath(
    map_status$path,
    winslash = "/",
    mustWork = TRUE
  )

  write.csv(map_status, map_status_manifest_csv, row.names = FALSE)

  strict_figure_pdf <- file.path(
    figure_dir,
    "Figure_S_strict_native_priority_persistence.pdf"
  )

  relaxed_figure_pdf <- file.path(
    figure_dir,
    "Figure_S_relaxed_native_priority_persistence.pdf"
  )

  tryCatch(
    {
      create_priority_figure(
        priority_raster = strict_priority,
        persistence_raster = strict_persistence_native,
        native_current = current_binary_native,
        output_pdf = strict_figure_pdf,
        title_text = "Strict native-range priority: 3/3 GCM suitable + no clamping"
      )
      write_log("Created strict priority figure: ", strict_figure_pdf)
    },
    error = function(e) {
      write_log("Strict priority figure skipped because: ", conditionMessage(e))
    }
  )

  tryCatch(
    {
      create_priority_figure(
        priority_raster = relaxed_priority,
        persistence_raster = relaxed_persistence_native,
        native_current = current_binary_native,
        output_pdf = relaxed_figure_pdf,
        title_text = "Relaxed native-range priority: 3/3 GCM suitable + clamping <= 1 GCM"
      )
      write_log("Created relaxed priority figure: ", relaxed_figure_pdf)
    },
    error = function(e) {
      write_log("Relaxed priority figure skipped because: ", conditionMessage(e))
    }
  )

  run_summary <- data.frame(
    species = species_name,
    MTSPS_threshold = mtsps_threshold,
    native_range_layer = native_range_layer,
    native_range_area_km2 = native_range_area_km2,
    current_suitable_native_km2 = current_suitable_native_km2,
    strict_core_refugia_km2 =
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "core_refugia"),
    strict_high_insitu_12plus_km2 =
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "high_insitu_12plus"),
    strict_exsitu_germplasm_km2 =
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "exsitu_germplasm"),
    strict_emerging_future_suitable_km2 =
      get_binary_area("strict_3of3_GCM_suitable_no_clamping", "emerging_future_suitable"),
    relaxed_core_refugia_km2 =
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "core_refugia"),
    relaxed_high_insitu_12plus_km2 =
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "high_insitu_12plus"),
    relaxed_exsitu_germplasm_km2 =
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "exsitu_germplasm"),
    relaxed_emerging_future_suitable_km2 =
      get_binary_area("relaxed_3of3_GCM_suitable_clamping_at_most_1", "emerging_future_suitable"),
    stringsAsFactors = FALSE
  )

  write.csv(run_summary, run_summary_csv, row.names = FALSE)

  writeLines(
    c(
      "SCRIPT 19 COMPLETED SUCCESSFULLY",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("MTSPS threshold:", format(mtsps_threshold, scientific = TRUE, digits = 8)),
      paste("Native range area (km2):", format(native_range_area_km2, scientific = FALSE)),
      paste("Current suitable native area (km2):", format(current_suitable_native_km2, scientific = FALSE)),
      paste("Strict core refugia (km2):", format(run_summary$strict_core_refugia_km2, scientific = FALSE)),
      paste("Strict high-priority in situ 12+ (km2):", format(run_summary$strict_high_insitu_12plus_km2, scientific = FALSE)),
      paste("Strict ex situ germplasm priority (km2):", format(run_summary$strict_exsitu_germplasm_km2, scientific = FALSE))
    ),
    completed_marker
  )

  write_log("Script completed successfully.")

  cat("\n======================================================================\n")
  cat("SCRIPT 19 COMPLETED — CLIMATE-STABLE REFUGIA AND PRIORITY ZONES\n")
  cat("======================================================================\n")
  cat(
    "MTSPS threshold: ",
    format(mtsps_threshold, scientific = TRUE, digits = 8),
    "\n",
    sep = ""
  )
  cat(
    "Native range area (km²): ",
    format(native_range_area_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat(
    "Current suitable native area (km²): ",
    format(current_suitable_native_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat(
    "Strict core climate-stable refugia (km²): ",
    format(run_summary$strict_core_refugia_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat(
    "Strict high-priority in situ 12+ (km²): ",
    format(run_summary$strict_high_insitu_12plus_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat(
    "Strict ex situ germplasm priority (km²): ",
    format(run_summary$strict_exsitu_germplasm_km2, big.mark = ",", scientific = FALSE),
    "\n",
    sep = ""
  )
  cat("\nKey area tables:\n")
  cat(native_summary_csv, "\n")
  cat(strict_priority_area_csv, "\n")
  cat(binary_priority_area_csv, "\n")
  cat(strict_vs_relaxed_csv, "\n")
  cat(manuscript_summary_csv, "\n")
  cat("\nOutput folder:\n")
  cat(output_root, "\n")

  invisible(list(
    native_summary = native_summary,
    strict_priority_area = strict_priority_area,
    relaxed_priority_area = relaxed_priority_area,
    binary_area_summary = binary_area_summary,
    manuscript_summary = manuscript_summary
  ))
}

# =============================================================================
# F. EXECUTE WITH EXPLICIT FAILURE REPORTING
# =============================================================================

tryCatch(
  {
    run_script19()
  },
  error = function(error) {

    dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    failure_message <- paste0(
      "SCRIPT 19 FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(error), "\n\n",
      "Execution log: ", execution_log, "\n"
    )

    writeLines(failure_message, failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 19 FAILED\n")
    cat("======================================================================\n")
    cat(failure_message, "\n")

    stop(error)
  }
)
