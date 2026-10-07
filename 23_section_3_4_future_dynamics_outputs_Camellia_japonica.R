# ============================================================
# SECTION 3.4 OUTPUT GENERATOR
# Camellia japonica
#
# Purpose:
#   Create rasters, figures, and tables for:
#   3.4. Future suitable-area dynamics and threshold sensitivity
#
# Inputs:
#   1) Current and future 4-class MTSPS rasters
#   2) World continent shapefile
#   3) Optional threshold-sensitivity table for Fig. 7
#
# Main outputs:
#   - Binary suitable/unsuitable rasters
#   - Transition rasters
#   - Fig. 7: Global MTSPS vs T10 suitable-area linegraph
#   - Fig. 8: Continental future-change linegraph under MTSPS
#   - Fig. 9: Representative spatial transition map
#   - Fig. Sx: Full 16-panel transition map
#   - Table 4: Compact continental future-change summary
#   - Full continental change table
#
# Class definitions:
#   1 = Unsuitable
#   2 = Low suitability
#   3 = Moderate suitability
#   4 = High suitability
#
# Transition definitions:
#   1 = Stable unsuitable
#   2 = Expansion
#   3 = Contraction
#   4 = Stable suitable
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# 1. Packages
# ------------------------------------------------------------
packages <- c("terra", "sf", "ggplot2", "dplyr", "tidyr")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(terra)
library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 2. User paths
# ------------------------------------------------------------
input_dir <- "E:/PREDICTION TREE/Camellia_japonica/global_reanalysis_outputs/13_MTSPS_4class_future_ensemble_area"

continent_shp <- "D:/Prediction_tree/continents/World_Continents (2).shp"

# Optional: used only for Fig. 7 MTSPS vs T10 threshold sensitivity
threshold_sensitivity_csv <- "E:/PREDICTION TREE/Camellia_japonica/global_reanalysis_outputs/20_manuscript_tables_and_figures/03_supplementary_tables/Table_S_threshold_sensitivity_all_scenarios.csv"

output_root <- "E:/PREDICTION TREE/Camellia_japonica/global_reanalysis_outputs/23_section_3_4_future_dynamics_outputs"

# ------------------------------------------------------------
# 3. Output folders
# ------------------------------------------------------------
dir_setup      <- file.path(output_root, "01_setup")
dir_binary     <- file.path(output_root, "02_binary_rasters")
dir_transition <- file.path(output_root, "03_transition_rasters")
dir_tables     <- file.path(output_root, "04_tables")
dir_figures    <- file.path(output_root, "05_figures")

dir.create(dir_setup, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_binary, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_transition, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_figures, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 4. General settings
# ------------------------------------------------------------
if (.Platform$OS.type == "windows") {
  windowsFonts(TNR = windowsFont("Times New Roman"))
}
base_font <- "TNR"

# Map-display sampling only. Area calculations use full rasters.
max_map_cells_representative <- 900000
max_map_cells_full16 <- 220000

ssp_levels <- c("SSP126", "SSP245", "SSP370", "SSP585")
period_levels <- c("2021–2040", "2041–2060", "2061–2080", "2081–2100")

transition_labels <- c(
  "1" = "Stable unsuitable",
  "2" = "Expansion",
  "3" = "Contraction",
  "4" = "Stable suitable"
)

transition_cols <- c(
  "Stable unsuitable" = "grey92",
  "Expansion" = "#D95F02",
  "Contraction" = "#7570B3",
  "Stable suitable" = "#1B9E77"
)

# ------------------------------------------------------------
# 5. Helper functions
# ------------------------------------------------------------
find_continent_col <- function(x) {
  nms <- names(x)
  candidates <- c(
    "continent", "CONTINENT", "Continent",
    "CONTINENT_N", "CONTINENT_NAME",
    "NAME", "Name", "name"
  )
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[1])
  stop("Could not identify a continent-name field in the continent shapefile. Available fields:\n",
       paste(nms, collapse = ", "))
}

parse_raster_info <- function(path) {
  bn <- basename(path)
  bn_upper <- toupper(bn)

  is_current <- grepl("CURRENT", bn_upper)
  ssp_match <- regmatches(bn_upper, regexpr("SSP[0-9]{3}", bn_upper))
  has_ssp <- length(ssp_match) > 0 && nchar(ssp_match) > 0

  period_patterns <- c("2021[_-]?2040", "2041[_-]?2060", "2061[_-]?2080", "2081[_-]?2100")
  period_clean <- NA_character_
  period_start <- NA_integer_
  period_end <- NA_integer_

  for (pat in period_patterns) {
    m <- regmatches(bn_upper, regexpr(pat, bn_upper))
    if (length(m) > 0 && nchar(m) > 0) {
      yrs <- unlist(regmatches(m, gregexpr("[0-9]{4}", m)))
      if (length(yrs) == 2) {
        period_start <- as.integer(yrs[1])
        period_end <- as.integer(yrs[2])
        period_clean <- paste0(yrs[1], "–", yrs[2])
      }
      break
    }
  }

  data.frame(
    path = path,
    filename = bn,
    is_current = is_current,
    ssp = ifelse(has_ssp, ssp_match, NA_character_),
    period_start = period_start,
    period_end = period_end,
    period_label = period_clean,
    stringsAsFactors = FALSE
  )
}

make_binary <- function(r) {
  # 1 = suitable, 0 = unsuitable
  ifel(r >= 2 & r <= 4, 1, 0)
}

make_transition <- function(current_bin, future_bin) {
  # 1 = stable unsuitable: current 0, future 0
  # 2 = expansion:          current 0, future 1
  # 3 = contraction:        current 1, future 0
  # 4 = stable suitable:    current 1, future 1
  out <- ifel(current_bin == 0 & future_bin == 0, 1,
         ifel(current_bin == 0 & future_bin == 1, 2,
         ifel(current_bin == 1 & future_bin == 0, 3,
         ifel(current_bin == 1 & future_bin == 1, 4, NA))))
  out
}

zonal_area_by_code <- function(code_raster, area_raster, zone_raster, code_type = "class") {
  # Combine continent ID and class/transition code into one zone value.
  # zone_value = continent_id * 100 + code
  combo <- zone_raster * 100 + code_raster
  z <- terra::zonal(area_raster, combo, fun = "sum", na.rm = TRUE, as.raster = FALSE)

  if (nrow(z) == 0) {
    return(data.frame())
  }

  names(z) <- c("zone_value", "area_km2")
  z$continent_id <- floor(z$zone_value / 100)
  z$code <- z$zone_value - z$continent_id * 100

  if (code_type == "class") {
    z$class_code <- z$code
  } else {
    z$transition_code <- z$code
  }

  z$zone_value <- NULL
  z$code <- NULL
  z
}

prep_map_df <- function(r, max_cells = 500000) {
  # Downsample for plotting only.
  # Area tables are calculated from full-resolution rasters.
  n <- terra::ncell(r)

  if (n > max_cells) {
    df <- terra::spatSample(
      r,
      size = max_cells,
      method = "regular",
      as.df = TRUE,
      xy = TRUE,
      na.rm = TRUE
    )
  } else {
    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  }

  value_col <- setdiff(names(df), c("x", "y"))[1]
  names(df)[names(df) == value_col] <- "value"
  df$value <- as.integer(round(df$value))
  df
}

# ------------------------------------------------------------
# 6. Check inputs
# ------------------------------------------------------------
if (!dir.exists(input_dir)) {
  stop("Input raster folder not found:\n", input_dir)
}

if (!file.exists(continent_shp)) {
  stop("Continent shapefile not found:\n", continent_shp)
}

all_tifs <- list.files(input_dir, pattern = "\\.tif$|\\.tiff$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

if (length(all_tifs) == 0) {
  stop("No GeoTIFF files found in:\n", input_dir)
}

manifest <- do.call(rbind, lapply(all_tifs, parse_raster_info))

# Select current raster
current_candidates <- manifest[manifest$is_current, ]
if (nrow(current_candidates) == 0) {
  stop("Could not find a current raster. The filename should contain 'current'.")
}
current_path <- current_candidates$path[1]

# Select future rasters
future_manifest <- manifest[!manifest$is_current & !is.na(manifest$ssp) & !is.na(manifest$period_label), ]
future_manifest <- future_manifest[future_manifest$ssp %in% ssp_levels, ]
future_manifest$ssp <- factor(future_manifest$ssp, levels = ssp_levels)
future_manifest$period_label <- factor(future_manifest$period_label, levels = period_levels)
future_manifest <- future_manifest[order(future_manifest$ssp, future_manifest$period_start), ]

if (nrow(future_manifest) != 16) {
  warning("Expected 16 future rasters, but found ", nrow(future_manifest), ". Check raster filenames.")
}

write.csv(manifest, file.path(dir_setup, "01_all_4class_raster_manifest.csv"), row.names = FALSE)
write.csv(future_manifest, file.path(dir_setup, "02_future_4class_raster_manifest_used.csv"), row.names = FALSE)

cat("\nCurrent raster:\n", current_path, "\n", sep = "")
cat("\nFuture rasters found:\n")
print(future_manifest[, c("ssp", "period_label", "filename")])

# ------------------------------------------------------------
# 7. Read current raster and continent polygons
# ------------------------------------------------------------
current_class <- terra::rast(current_path)
current_class <- round(current_class)

continents <- terra::vect(continent_shp)
continents <- terra::project(continents, terra::crs(current_class))

continent_col <- find_continent_col(continents)
continents$continent_name <- as.character(continents[[continent_col]])

# Remove Antarctica if present
continents <- continents[!(tolower(continents$continent_name) %in% c("antarctica", "antarctique", "antarctida")), ]

# Dissolve continents by name if multiple polygons occur
continents_sf <- sf::st_as_sf(continents)
continents_sf <- continents_sf |>
  dplyr::group_by(continent_name) |>
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")

continents <- terra::vect(continents_sf)
continents$continent_id <- seq_len(nrow(continents))

continent_lookup <- data.frame(
  continent_id = continents$continent_id,
  continent = continents$continent_name,
  stringsAsFactors = FALSE
)

write.csv(continent_lookup, file.path(dir_setup, "03_continent_lookup.csv"), row.names = FALSE)

continent_r <- terra::rasterize(continents, current_class, field = "continent_id")
area_r <- terra::cellSize(current_class, unit = "km")

# ------------------------------------------------------------
# 8. Create current binary raster
# ------------------------------------------------------------
current_bin <- make_binary(current_class)
names(current_bin) <- "current_binary_suitable"

current_bin_path <- file.path(dir_binary, "Camellia_japonica_current_MTSPS_binary_suitable.tif")
terra::writeRaster(current_bin, current_bin_path, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))

# ------------------------------------------------------------
# 9. Current area by continent and class
# ------------------------------------------------------------
current_area_raw <- zonal_area_by_code(current_class, area_r, continent_r, code_type = "class")

current_area <- current_area_raw |>
  left_join(continent_lookup, by = "continent_id") |>
  mutate(
    class_name = case_when(
      class_code == 1 ~ "Unsuitable",
      class_code == 2 ~ "Low suitability",
      class_code == 3 ~ "Moderate suitability",
      class_code == 4 ~ "High suitability",
      TRUE ~ "Other"
    ),
    area_10_4_km2 = area_km2 / 10000
  ) |>
  select(continent, class_code, class_name, area_km2, area_10_4_km2)

current_suitable_continent <- current_area |>
  filter(class_code %in% c(2, 3, 4)) |>
  group_by(continent) |>
  summarise(
    current_suitable_km2 = sum(area_km2, na.rm = TRUE),
    current_suitable_10_4_km2 = current_suitable_km2 / 10000,
    .groups = "drop"
  )

write.csv(current_area, file.path(dir_tables, "01_current_4class_area_by_continent_long.csv"), row.names = FALSE)
write.csv(current_suitable_continent, file.path(dir_tables, "02_current_suitable_area_by_continent.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 10. Create future binary rasters, transition rasters, and area tables
# ------------------------------------------------------------
future_area_list <- list()
transition_area_list <- list()
transition_global_list <- list()
transition_manifest <- data.frame()

for (i in seq_len(nrow(future_manifest))) {

  row_i <- future_manifest[i, ]

  ssp_i <- as.character(row_i$ssp)
  period_i <- as.character(row_i$period_label)
  period_file <- gsub("–", "_", period_i)
  scenario_id <- paste0(ssp_i, "_", period_file)

  cat("\nProcessing: ", scenario_id, "\n", sep = "")

  future_class <- terra::rast(row_i$path)
  future_class <- round(future_class)

  # Align future raster to current raster if needed
  if (!terra::compareGeom(current_class, future_class, stopOnError = FALSE)) {
    future_class <- terra::resample(future_class, current_class, method = "near")
  }

  future_bin <- make_binary(future_class)
  names(future_bin) <- paste0("future_binary_", scenario_id)

  future_bin_path <- file.path(
    dir_binary,
    paste0("Camellia_japonica_future_", scenario_id, "_MTSPS_binary_suitable.tif")
  )

  terra::writeRaster(future_bin, future_bin_path, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))

  # Transition raster
  trans <- make_transition(current_bin, future_bin)
  names(trans) <- paste0("transition_", scenario_id)

  trans_path <- file.path(
    dir_transition,
    paste0("Camellia_japonica_transition_", scenario_id, "_MTSPS.tif")
  )

  terra::writeRaster(trans, trans_path, overwrite = TRUE, wopt = list(gdal = c("COMPRESS=LZW")))

  transition_manifest <- rbind(
    transition_manifest,
    data.frame(
      ssp = ssp_i,
      period_label = period_i,
      scenario_id = scenario_id,
      future_4class_tif = row_i$path,
      future_binary_tif = future_bin_path,
      transition_tif = trans_path,
      stringsAsFactors = FALSE
    )
  )

  # Future area by continent and class
  fut_area_raw <- zonal_area_by_code(future_class, area_r, continent_r, code_type = "class")

  fut_area <- fut_area_raw |>
    left_join(continent_lookup, by = "continent_id") |>
    mutate(
      ssp = ssp_i,
      period_label = period_i,
      period_start = row_i$period_start,
      period_end = row_i$period_end,
      scenario_id = scenario_id,
      class_name = case_when(
        class_code == 1 ~ "Unsuitable",
        class_code == 2 ~ "Low suitability",
        class_code == 3 ~ "Moderate suitability",
        class_code == 4 ~ "High suitability",
        TRUE ~ "Other"
      ),
      area_10_4_km2 = area_km2 / 10000
    ) |>
    select(ssp, period_label, period_start, period_end, scenario_id,
           continent, class_code, class_name, area_km2, area_10_4_km2)

  future_area_list[[scenario_id]] <- fut_area

  # Transition area by continent
  trans_area_raw <- zonal_area_by_code(trans, area_r, continent_r, code_type = "transition")

  trans_area <- trans_area_raw |>
    left_join(continent_lookup, by = "continent_id") |>
    mutate(
      ssp = ssp_i,
      period_label = period_i,
      period_start = row_i$period_start,
      period_end = row_i$period_end,
      scenario_id = scenario_id,
      transition_name = case_when(
        transition_code == 1 ~ "Stable unsuitable",
        transition_code == 2 ~ "Expansion",
        transition_code == 3 ~ "Contraction",
        transition_code == 4 ~ "Stable suitable",
        TRUE ~ "Other"
      ),
      area_10_4_km2 = area_km2 / 10000
    ) |>
    select(ssp, period_label, period_start, period_end, scenario_id,
           continent, transition_code, transition_name, area_km2, area_10_4_km2)

  transition_area_list[[scenario_id]] <- trans_area

  # Global transition area
  trans_global <- trans_area |>
    group_by(ssp, period_label, period_start, period_end, scenario_id,
             transition_code, transition_name) |>
    summarise(
      area_km2 = sum(area_km2, na.rm = TRUE),
      area_10_4_km2 = area_km2 / 10000,
      .groups = "drop"
    )

  transition_global_list[[scenario_id]] <- trans_global
}

write.csv(transition_manifest, file.path(dir_setup, "04_transition_raster_manifest.csv"), row.names = FALSE)

future_area_all <- bind_rows(future_area_list)
transition_area_all <- bind_rows(transition_area_list)
transition_global_all <- bind_rows(transition_global_list)

write.csv(future_area_all, file.path(dir_tables, "03_future_4class_area_by_continent_long.csv"), row.names = FALSE)
write.csv(transition_area_all, file.path(dir_tables, "04_transition_area_by_continent_long_MTSPS.csv"), row.names = FALSE)
write.csv(transition_global_all, file.path(dir_tables, "05_transition_area_global_long_MTSPS.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 11. Future suitable area and change by continent
# ------------------------------------------------------------
future_suitable_continent <- future_area_all |>
  filter(class_code %in% c(2, 3, 4)) |>
  group_by(ssp, period_label, period_start, period_end, scenario_id, continent) |>
  summarise(
    future_suitable_km2 = sum(area_km2, na.rm = TRUE),
    future_suitable_10_4_km2 = future_suitable_km2 / 10000,
    .groups = "drop"
  ) |>
  left_join(current_suitable_continent, by = "continent") |>
  mutate(
    change_km2 = future_suitable_km2 - current_suitable_km2,
    change_10_4_km2 = change_km2 / 10000,
    percent_change = ifelse(current_suitable_km2 > 0, 100 * change_km2 / current_suitable_km2, NA_real_)
  ) |>
  arrange(continent, ssp, period_start)

write.csv(future_suitable_continent, file.path(dir_tables, "06_future_suitable_area_change_by_continent_long_MTSPS.csv"), row.names = FALSE)

# Wide table close to old Table 4 style
future_change_wide <- future_suitable_continent |>
  mutate(
    scenario_period = paste(ssp, period_label),
    change_10_4_km2 = round(change_10_4_km2, 3)
  ) |>
  select(continent, scenario_period, change_10_4_km2) |>
  tidyr::pivot_wider(names_from = scenario_period, values_from = change_10_4_km2)

write.csv(future_change_wide, file.path(dir_tables, "07_Table_S_full_future_change_by_continent_wide_MTSPS_10_4km2.csv"), row.names = FALSE)

# Compact manuscript Table 4
table4_summary <- future_suitable_continent |>
  group_by(continent) |>
  summarise(
    current_suitable_10_4_km2 = first(current_suitable_10_4_km2),
    min_future_suitable_10_4_km2 = min(future_suitable_10_4_km2, na.rm = TRUE),
    max_future_suitable_10_4_km2 = max(future_suitable_10_4_km2, na.rm = TRUE),
    min_change_10_4_km2 = min(change_10_4_km2, na.rm = TRUE),
    max_change_10_4_km2 = max(change_10_4_km2, na.rm = TRUE),
    max_percent_change = max(percent_change, na.rm = TRUE),
    scenario_period_of_max_gain = paste(
      ssp[which.max(change_10_4_km2)],
      period_label[which.max(change_10_4_km2)]
    ),
    .groups = "drop"
  ) |>
  mutate(
    across(
      c(current_suitable_10_4_km2,
        min_future_suitable_10_4_km2,
        max_future_suitable_10_4_km2,
        min_change_10_4_km2,
        max_change_10_4_km2,
        max_percent_change),
      ~ round(.x, 3)
    )
  ) |>
  arrange(desc(max_change_10_4_km2))

write.csv(table4_summary, file.path(dir_tables, "08_Table_4_continental_summary_future_change_MTSPS.csv"), row.names = FALSE)

# Transition wide table: expansion/contraction/stable under MTSPS
transition_wide <- transition_area_all |>
  select(ssp, period_label, period_start, period_end, scenario_id, continent,
         transition_name, area_10_4_km2) |>
  mutate(area_10_4_km2 = round(area_10_4_km2, 3)) |>
  tidyr::pivot_wider(names_from = transition_name, values_from = area_10_4_km2, values_fill = 0) |>
  mutate(
    Net_change = Expansion - Contraction
  ) |>
  arrange(continent, ssp, period_start)

write.csv(transition_wide, file.path(dir_tables, "09_transition_components_by_continent_wide_MTSPS_10_4km2.csv"), row.names = FALSE)

# ------------------------------------------------------------
# 12. Figure 7: Global threshold sensitivity, optional
# ------------------------------------------------------------
if (file.exists(threshold_sensitivity_csv)) {

  sens <- read.csv(threshold_sensitivity_csv, stringsAsFactors = FALSE, check.names = FALSE)

  if (all(c("threshold_id", "ssp", "period_start", "period_end", "future_suitable_10_4_km2") %in% names(sens))) {

    fig7_df <- sens |>
      filter(threshold_id %in% c("MTSPS", "T10")) |>
      mutate(
        ssp = factor(ssp, levels = ssp_levels),
        period_label = paste0(period_start, "–", period_end),
        period_label = factor(period_label, levels = period_levels),
        area_10_4 = as.numeric(future_suitable_10_4_km2),
        label_text = sprintf("%.3f", round(area_10_4, 3))
      ) |>
      arrange(ssp, period_start, threshold_id)

    fig7_df$x_base <- match(fig7_df$period_label, period_levels)
    fig7_df$x_plot <- ifelse(fig7_df$threshold_id == "MTSPS", fig7_df$x_base - 0.06, fig7_df$x_base + 0.06)

    y_span7 <- max(fig7_df$area_10_4, na.rm = TRUE) - min(fig7_df$area_10_4, na.rm = TRUE)
    if (y_span7 == 0) y_span7 <- 1

    fig7_df$label_y <- ifelse(
      fig7_df$threshold_id == "MTSPS",
      fig7_df$area_10_4 + y_span7 * 0.055,
      fig7_df$area_10_4 + y_span7 * 0.035
    )

    p7 <- ggplot() +
      geom_line(
        data = fig7_df,
        aes(x = x_plot, y = area_10_4, group = threshold_id, linetype = threshold_id),
        colour = "black", linewidth = 0.9
      ) +
      geom_point(
        data = fig7_df,
        aes(x = x_plot, y = area_10_4, shape = threshold_id, fill = threshold_id),
        colour = "black", size = 3.8, stroke = 1.0
      ) +
      geom_text(
        data = fig7_df,
        aes(x = x_plot, y = label_y, label = label_text),
        family = base_font, size = 4.8, colour = "black", vjust = 0
      ) +
      facet_wrap(~ ssp, ncol = 2) +
      scale_x_continuous(
        breaks = 1:4,
        labels = period_levels,
        expand = expansion(mult = c(0.12, 0.12))
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.08, 0.23))) +
      scale_linetype_manual(values = c("MTSPS" = "solid", "T10" = "dashed"), name = NULL) +
      scale_shape_manual(values = c("MTSPS" = 21, "T10" = 21), name = NULL) +
      scale_fill_manual(values = c("MTSPS" = "black", "T10" = "white"), name = NULL) +
      labs(
        title = "Future climatically suitable area under MTSPS and T10 thresholds",
        x = "Projection period",
        y = expression("Suitable area (" * 10^4 * " km"^2 * ")")
      ) +
      coord_cartesian(clip = "off") +
      theme_bw(base_family = base_font, base_size = 17) +
      theme(
        text = element_text(family = base_font, colour = "black"),
        plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 18)),
        axis.title.x = element_text(size = 19, face = "bold", margin = margin(t = 14)),
        axis.title.y = element_text(size = 19, face = "bold", margin = margin(r = 14)),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 15),
        strip.text = element_text(size = 19, face = "bold"),
        legend.position = "top",
        legend.text = element_text(size = 16),
        legend.key.width = unit(1.3, "cm"),
        panel.grid.major = element_line(colour = "grey80", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1.5, "lines"),
        plot.margin = margin(22, 30, 22, 22)
      )

    ggsave(file.path(dir_figures, "Figure_7_global_threshold_sensitivity_MTSPS_T10_10_4km2.png"),
           p7, width = 15.5, height = 11, dpi = 700, bg = "white")

    ggsave(file.path(dir_figures, "Figure_7_global_threshold_sensitivity_MTSPS_T10_10_4km2.pdf"),
           p7, width = 15.5, height = 11, bg = "white")
  }
}

# ------------------------------------------------------------
# 13. Figure 8: Continental future change linegraph under MTSPS
# ------------------------------------------------------------
fig8_df <- future_suitable_continent |>
  mutate(
    ssp = factor(ssp, levels = ssp_levels),
    period_label = factor(period_label, levels = period_levels),
    continent = factor(continent, levels = sort(unique(continent)))
  )

p8 <- ggplot(
  fig8_df,
  aes(x = period_label, y = change_10_4_km2, group = ssp, linetype = ssp, shape = ssp)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey55") +
  geom_line(colour = "black", linewidth = 0.85) +
  geom_point(colour = "black", fill = "white", size = 3.1, stroke = 0.9) +
  facet_wrap(~ continent, ncol = 2, scales = "free_y") +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted"), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  labs(
    title = "Future change in climatically suitable area by continent under MTSPS",
    x = "Projection period",
    y = expression("Change in suitable area (" * 10^4 * " km"^2 * ")")
  ) +
  theme_bw(base_family = base_font, base_size = 17) +
  theme(
    text = element_text(family = base_font, colour = "black"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 18)),
    axis.title.x = element_text(size = 19, face = "bold", margin = margin(t = 14)),
    axis.title.y = element_text(size = 19, face = "bold", margin = margin(r = 14)),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(size = 18, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.key.width = unit(1.3, "cm"),
    panel.grid.major = element_line(colour = "grey82", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    plot.margin = margin(22, 30, 22, 22)
  )

ggsave(file.path(dir_figures, "Figure_8_continental_future_change_MTSPS_linegraph.png"),
       p8, width = 15.5, height = 13.5, dpi = 700, bg = "white")

ggsave(file.path(dir_figures, "Figure_8_continental_future_change_MTSPS_linegraph.pdf"),
       p8, width = 15.5, height = 13.5, bg = "white")

# ------------------------------------------------------------
# 14. Figure 9: Representative spatial transition maps
# ------------------------------------------------------------
transition_manifest$period_label <- factor(transition_manifest$period_label, levels = period_levels)
transition_manifest$ssp <- factor(transition_manifest$ssp, levels = ssp_levels)

representative_manifest <- transition_manifest |>
  filter(
    (ssp == "SSP245" & period_label == "2081–2100") |
      (ssp == "SSP585" & period_label == "2081–2100")
  ) |>
  arrange(ssp, period_label)

if (nrow(representative_manifest) > 0) {

  map_dfs <- list()

  for (i in seq_len(nrow(representative_manifest))) {
    r <- terra::rast(representative_manifest$transition_tif[i])
    df_map <- prep_map_df(r, max_cells = max_map_cells_representative)
    df_map$transition_name <- transition_labels[as.character(df_map$value)]
    df_map$scenario <- paste(representative_manifest$ssp[i], representative_manifest$period_label[i])
    map_dfs[[i]] <- df_map
  }

  fig9_df <- bind_rows(map_dfs)
  fig9_df$transition_name <- factor(fig9_df$transition_name, levels = unname(transition_labels))

  continents_plot <- sf::st_as_sf(continents)

  p9 <- ggplot() +
    geom_raster(
      data = fig9_df,
      aes(x = x, y = y, fill = transition_name)
    ) +
    geom_sf(
      data = continents_plot,
      fill = NA,
      colour = "grey30",
      linewidth = 0.12,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ scenario, ncol = 1) +
    scale_fill_manual(values = transition_cols, name = NULL, drop = FALSE) +
    coord_sf(expand = FALSE) +
    labs(
      title = "Representative future distribution-change patterns under MTSPS"
    ) +
    theme_void(base_family = base_font, base_size = 16) +
    theme(
      text = element_text(family = base_font, colour = "black"),
      plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 12)),
      strip.text = element_text(size = 18, face = "bold"),
      legend.position = "bottom",
      legend.text = element_text(size = 15),
      legend.key.width = unit(1.1, "cm"),
      panel.spacing = unit(0.9, "lines"),
      plot.margin = margin(18, 18, 18, 18)
    )

  ggsave(file.path(dir_figures, "Figure_9_representative_transition_maps_MTSPS.png"),
         p9, width = 12.5, height = 12, dpi = 700, bg = "white")

  ggsave(file.path(dir_figures, "Figure_9_representative_transition_maps_MTSPS.pdf"),
         p9, width = 12.5, height = 12, bg = "white")
}

# ------------------------------------------------------------
# 15. Supplementary full 16-panel transition map
# ------------------------------------------------------------
full_map_dfs <- list()

for (i in seq_len(nrow(transition_manifest))) {
  r <- terra::rast(transition_manifest$transition_tif[i])
  df_map <- prep_map_df(r, max_cells = max_map_cells_full16)
  df_map$transition_name <- transition_labels[as.character(df_map$value)]
  df_map$ssp <- as.character(transition_manifest$ssp[i])
  df_map$period_label <- as.character(transition_manifest$period_label[i])
  full_map_dfs[[i]] <- df_map
}

figS_df <- bind_rows(full_map_dfs)
figS_df$transition_name <- factor(figS_df$transition_name, levels = unname(transition_labels))
figS_df$ssp <- factor(figS_df$ssp, levels = ssp_levels)
figS_df$period_label <- factor(figS_df$period_label, levels = period_levels)

continents_plot <- sf::st_as_sf(continents)

pS <- ggplot() +
  geom_raster(
    data = figS_df,
    aes(x = x, y = y, fill = transition_name)
  ) +
  geom_sf(
    data = continents_plot,
    fill = NA,
    colour = "grey30",
    linewidth = 0.08,
    inherit.aes = FALSE
  ) +
  facet_grid(ssp ~ period_label) +
  scale_fill_manual(values = transition_cols, name = NULL, drop = FALSE) +
  coord_sf(expand = FALSE) +
  labs(
    title = "Future distribution-change patterns under MTSPS across all SSPs and periods"
  ) +
  theme_void(base_family = base_font, base_size = 15) +
  theme(
    text = element_text(family = base_font, colour = "black"),
    plot.title = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 12)),
    strip.text = element_text(size = 15, face = "bold"),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.key.width = unit(1.0, "cm"),
    panel.spacing = unit(0.25, "lines"),
    plot.margin = margin(18, 18, 18, 18)
  )

ggsave(file.path(dir_figures, "Figure_S_full_16_transition_maps_MTSPS.png"),
       pS, width = 17, height = 11.5, dpi = 700, bg = "white")

ggsave(file.path(dir_figures, "Figure_S_full_16_transition_maps_MTSPS.pdf"),
       pS, width = 17, height = 11.5, bg = "white")

# ------------------------------------------------------------
# 16. Final report
# ------------------------------------------------------------
cat("\n============================================================\n")
cat("SECTION 3.4 OUTPUTS CREATED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Output root:\n", output_root, "\n\n", sep = "")
cat("Derived rasters:\n")
cat("  Binary rasters:     ", dir_binary, "\n", sep = "")
cat("  Transition rasters: ", dir_transition, "\n\n", sep = "")
cat("Tables:\n")
cat("  ", dir_tables, "\n\n", sep = "")
cat("Figures:\n")
cat("  ", dir_figures, "\n\n", sep = "")
cat("Recommended manuscript items:\n")
cat("  Fig. 7: Figure_7_global_threshold_sensitivity_MTSPS_T10_10_4km2.pdf\n")
cat("  Fig. 8: Figure_8_continental_future_change_MTSPS_linegraph.pdf\n")
cat("  Fig. 9: Figure_9_representative_transition_maps_MTSPS.pdf\n")
cat("  Table 4: 08_Table_4_continental_summary_future_change_MTSPS.csv\n")
cat("  Supplementary full map: Figure_S_full_16_transition_maps_MTSPS.pdf\n")
cat("============================================================\n")
