# =============================================================================
# 20_MANUSCRIPT_TABLES_AND_FIGURES_CAMELLIA_JAPONICA_V2_DEVICE_FIX.R
# =============================================================================
# PURPOSE
# -------
# Compile final manuscript-ready tables and figures for the revised
# Camellia japonica global SDM analysis.
#
# This script DOES NOT rerun MaxEnt or any spatial modelling step.
# It only reads the verified outputs from Scripts 13, 17A, 18, 18B, and 19,
# then produces clean tables and figures for manuscript Results,
# Supplementary Materials, and reviewer response.
#
# FINAL WORKFLOW ASSUMPTIONS
# --------------------------
# Occurrence records: 828
# Predictors: Bio2, Bio6, Bio14, Bio18, Elv
# MaxEnt tuning: FC = lph; RM = 1
# Primary threshold: MTSPS = fixed final current-model threshold
# Main future-confidence rule:
#     3/3 GCM suitable + no clamping
# Conservation priority:
#     restricted to native_wild_range_union
#
# OUTPUTS
# -------
# 1. Manuscript-ready summary tables
# 2. Publication-style PDF figures:
#      - future suitable-area trajectories under MTSPS and T10
#      - uncertainty and 3/3 GCM agreement trajectories
#      - clamping-aware high-confidence suitability
#      - native-range conservation priority bar chart
#      - final native-range map panel
# 3. Supplementary table inventory
# 4. Figure inventory and recommended manuscript placement
# =============================================================================

options(warn = 1)

cat("\n")
cat("======================================================================\n")
cat("STARTING SCRIPT 20 — MANUSCRIPT TABLES AND FIGURES\n")
cat("Species: Camellia japonica\n")
cat("======================================================================\n")

script_started_at <- Sys.time()

# =============================================================================
# A. USER SETTINGS
# =============================================================================

species_name <- "Camellia_japonica"

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# For figures.
figure_width <- 8.8
figure_height <- 5.8
map_figure_width <- 12.5
map_figure_height <- 5.2

# Use high-resolution PNG copies in addition to PDF.
write_png_copies <- TRUE
png_width_px <- 3600
png_height_px <- 2400
png_res <- 400

# Main manuscript threshold and supplementary threshold to display.
main_threshold_id <- "MTSPS"
supplementary_threshold_id <- "T10"

# =============================================================================
# B. INPUT PATHS
# =============================================================================

global_output_root <- file.path(
  project_dir,
  "global_reanalysis_outputs"
)

script13_root <- file.path(
  global_output_root,
  "13_MTSPS_4class_future_ensemble_area"
)

script17_root <- file.path(
  global_output_root,
  "17_threshold_sensitivity_MTSPS_MTP_T10"
)

script18_root <- file.path(
  global_output_root,
  "18_interGCM_uncertainty_agreement_extrapolation"
)

script18B_root <- file.path(
  global_output_root,
  "18B_clamping_quantification_extrapolation_confidence"
)

script19_root <- file.path(
  global_output_root,
  "19_climate_stable_refugia_conservation_priority"
)

script17_future_area_csv <- file.path(
  script17_root,
  "05_area_tables",
  "02_future_suitable_area_by_threshold_long.csv"
)

script17_threshold_summary_csv <- file.path(
  script17_root,
  "05_area_tables",
  "07_threshold_comparison_summary.csv"
)

script17_threshold_values_csv <- file.path(
  script17_root,
  "02_thresholds_and_current_scores",
  "01_threshold_values.csv"
)

script17_occurrence_prediction_summary_csv <- file.path(
  script17_root,
  "02_thresholds_and_current_scores",
  "02_current_prediction_summary_at_occurrences.csv"
)

script18_agreement_area_csv <- file.path(
  script18_root,
  "05_area_tables",
  "01_GCM_agreement_area_by_scenario_long.csv"
)

script18_sd_summary_csv <- file.path(
  script18_root,
  "05_area_tables",
  "03_inter_GCM_SD_summary_by_scenario.csv"
)

script18_high_conf_transition_csv <- file.path(
  script18_root,
  "05_area_tables",
  "04_high_confidence_transition_area_by_scenario.csv"
)

script18_scenario_confidence_csv <- file.path(
  script18_root,
  "05_area_tables",
  "06_scenario_confidence_summary.csv"
)

script18_extrapolation_inventory_csv <- file.path(
  script18_root,
  "01_setup",
  "04_existing_MESS_clamping_inventory.csv"
)

script18B_clamping_area_csv <- file.path(
  script18B_root,
  "06_area_tables",
  "02_clamping_GCM_count_and_high_confidence_area_by_scenario.csv"
)

script18B_priority_readiness_csv <- file.path(
  script18B_root,
  "06_area_tables",
  "05_clamping_priority_readiness_summary.csv"
)

script19_native_summary_csv <- file.path(
  script19_root,
  "07_area_tables",
  "01_native_range_and_current_suitability_summary.csv"
)

script19_strict_priority_area_csv <- file.path(
  script19_root,
  "07_area_tables",
  "04_strict_integrated_priority_area_native_range.csv"
)

script19_binary_priority_area_csv <- file.path(
  script19_root,
  "07_area_tables",
  "06_binary_priority_mask_area_summary.csv"
)

script19_strict_vs_relaxed_csv <- file.path(
  script19_root,
  "07_area_tables",
  "07_strict_vs_relaxed_priority_sensitivity_summary.csv"
)

script19_manuscript_summary_csv <- file.path(
  script19_root,
  "07_area_tables",
  "08_manuscript_ready_priority_summary.csv"
)

# Key rasters for final map panel.
script19_current_native_tif <- file.path(
  script19_root,
  "03_current_suitability_native",
  paste0(species_name, "_current_MTSPS_suitable_native_range.tif")
)

script19_strict_persistence_native_tif <- file.path(
  script19_root,
  "04_future_persistence_maps",
  paste0(species_name, "_strict_future_persistence_0_16_native_range.tif")
)

script19_strict_priority_native_tif <- file.path(
  script19_root,
  "05_integrated_priority_maps",
  paste0(species_name, "_strict_integrated_conservation_priority_native_range.tif")
)

# =============================================================================
# C. OUTPUT PATHS
# =============================================================================

output_root <- file.path(
  global_output_root,
  "20_manuscript_tables_and_figures"
)

status_dir <- file.path(output_root, "00_run_status")
table_dir <- file.path(output_root, "01_manuscript_tables")
figure_dir <- file.path(output_root, "02_manuscript_figures")
supplement_dir <- file.path(output_root, "03_supplementary_tables")
inventory_dir <- file.path(output_root, "04_inventory")
log_dir <- file.path(output_root, "05_logs")

started_marker <- file.path(status_dir, "20_STARTED.txt")
completed_marker <- file.path(status_dir, "20_COMPLETED.txt")
failed_marker <- file.path(status_dir, "20_FAILED.txt")
execution_log <- file.path(log_dir, "20_execution_log.txt")

table_model_setup_csv <- file.path(
  table_dir,
  "Table_1_model_setup_and_thresholds.csv"
)

table_threshold_sensitivity_csv <- file.path(
  table_dir,
  "Table_2_threshold_sensitivity_future_suitable_area.csv"
)

table_uncertainty_csv <- file.path(
  table_dir,
  "Table_3_GCM_agreement_and_uncertainty_summary.csv"
)

table_priority_csv <- file.path(
  table_dir,
  "Table_4_native_range_conservation_priority_summary.csv"
)

table_manuscript_key_numbers_csv <- file.path(
  table_dir,
  "Table_5_key_numbers_for_results_text.csv"
)

supplement_threshold_long_csv <- file.path(
  supplement_dir,
  "Table_S_threshold_sensitivity_all_scenarios.csv"
)

supplement_agreement_long_csv <- file.path(
  supplement_dir,
  "Table_S_GCM_agreement_all_scenarios.csv"
)

supplement_clamping_long_csv <- file.path(
  supplement_dir,
  "Table_S_clamping_confidence_all_scenarios.csv"
)

supplement_priority_strict_relaxed_csv <- file.path(
  supplement_dir,
  "Table_S_strict_vs_relaxed_priority_sensitivity.csv"
)

figure_future_area_pdf <- file.path(
  figure_dir,
  "Figure_1_future_suitable_area_threshold_sensitivity.pdf"
)

figure_gcm_agreement_pdf <- file.path(
  figure_dir,
  "Figure_2_GCM_agreement_uncertainty.pdf"
)

figure_clamping_pdf <- file.path(
  figure_dir,
  "Figure_3_clamping_aware_high_confidence_suitability.pdf"
)

figure_priority_bar_pdf <- file.path(
  figure_dir,
  "Figure_4_native_range_priority_area_summary.pdf"
)

figure_priority_map_pdf <- file.path(
  figure_dir,
  "Figure_5_native_range_priority_map_panel.pdf"
)

figure_inventory_csv <- file.path(
  inventory_dir,
  "figure_inventory_and_recommended_placement.csv"
)

table_inventory_csv <- file.path(
  inventory_dir,
  "table_inventory_and_recommended_placement.csv"
)

run_summary_csv <- file.path(
  inventory_dir,
  "script20_run_summary.csv"
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

read_required_csv <- function(path, label) {
  if (!file.exists(path)) {
    stop("\nMissing required ", label, ":\n", path, call. = FALSE)
  }

  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

check_required_files <- function(paths) {
  missing_paths <- paths[!file.exists(paths)]

  if (length(missing_paths) > 0L) {
    stop(
      "\nRequired file(s) missing:\n",
      paste(missing_paths, collapse = "\n"),
      call. = FALSE
    )
  }
}

format_km2 <- function(x, digits = 1) {
  format(round(x, digits), big.mark = ",", scientific = FALSE)
}

format_percent <- function(x, digits = 2) {
  paste0(format(round(x, digits), nsmall = digits), "%")
}

scenario_order_factor <- function(ssp, period_start, period_end) {
  ordered_labels <- as.vector(
    outer(
      c("SSP126", "SSP245", "SSP370", "SSP585"),
      c("2021-2040", "2041-2060", "2061-2080", "2081-2100"),
      paste
    )
  )

  labels <- paste0(ssp, " ", period_start, "-", period_end)
  factor(labels, levels = ordered_labels, ordered = TRUE)
}

ssp_to_numeric <- function(ssp) {
  match(ssp, c("SSP126", "SSP245", "SSP370", "SSP585"))
}

period_midpoint <- function(period_start, period_end) {
  (safe_number(period_start) + safe_number(period_end)) / 2
}

open_pdf_device <- function(path, width = figure_width, height = figure_height) {
  grDevices::pdf(path, width = width, height = height, useDingbats = FALSE)
}

copy_pdf_to_png <- function(pdf_path, png_path) {
  if (!write_png_copies) return(invisible(FALSE))

  # Use magick only if installed; otherwise silently skip PNG export.
  if (!requireNamespace("magick", quietly = TRUE)) {
    write_log("PNG copy skipped because package 'magick' is not installed: ", pdf_path)
    return(invisible(FALSE))
  }

  image <- magick::image_read_pdf(
    pdf_path,
    density = png_res
  )

  magick::image_write(
    image,
    path = png_path,
    format = "png"
  )

  invisible(TRUE)
}

save_base_plot_dual <- function(pdf_path, png_path, plot_function,
                                width = figure_width,
                                height = figure_height) {
  # Open and close each graphics device explicitly.
  # Do not use on.exit(dev.off()) here because the device is already closed
  # after each plot; a second dev.off() call would try to close the null device.
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  plot_function()
  grDevices::dev.off()

  if (write_png_copies) {
    grDevices::png(
      png_path,
      width = png_width_px,
      height = png_height_px,
      res = png_res
    )
    plot_function()
    grDevices::dev.off()
  }
}

# =============================================================================
# E. PLOTTING FUNCTIONS
# =============================================================================

plot_future_area <- function(future_area) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  thresholds_to_plot <- c(main_threshold_id, supplementary_threshold_id)

  data <- future_area[
    future_area$threshold_id %in% thresholds_to_plot,
    ,
    drop = FALSE
  ]

  ssps <- c("SSP126", "SSP245", "SSP370", "SSP585")
  threshold_lty <- c(MTSPS = 1, T10 = 2)

  y <- data$future_suitable_km2 / 1000000
  y_lim <- range(y, finite = TRUE)
  y_pad <- diff(y_lim) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3.0, 1.2), oma = c(0, 0, 1.5, 0))

  for (ssp in ssps) {
    subset_ssp <- data[data$ssp == ssp, , drop = FALSE]

    plot(
      NA,
      xlim = range(period_midpoint(subset_ssp$period_start, subset_ssp$period_end)),
      ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad),
      xlab = "Projection period midpoint",
      ylab = expression("Suitable area (" * 10^6 * " km"^2 * ")"),
      main = ssp,
      axes = FALSE
    )

    axis(1, at = period_midpoint(unique(subset_ssp$period_start), unique(subset_ssp$period_end)),
         labels = unique(paste0(subset_ssp$period_start, "-", subset_ssp$period_end)),
         las = 2, cex.axis = 0.70)
    axis(2, las = 1)
    box()
    grid(col = "grey85")

    for (threshold in thresholds_to_plot) {
      d <- subset_ssp[subset_ssp$threshold_id == threshold, , drop = FALSE]
      d <- d[order(d$period_start), , drop = FALSE]

      lines(
        period_midpoint(d$period_start, d$period_end),
        d$future_suitable_km2 / 1000000,
        type = "b",
        pch = ifelse(threshold == main_threshold_id, 16, 1),
        lty = threshold_lty[threshold],
        lwd = 2
      )
    }

    legend(
      "topleft",
      legend = c("MTSPS primary", "T10 sensitivity"),
      lty = c(1, 2),
      pch = c(16, 1),
      lwd = 2,
      bty = "n",
      cex = 0.80
    )
  }

  mtext(
    "Future climatically suitable area under primary and supplementary thresholds",
    side = 3,
    outer = TRUE,
    font = 2,
    cex = 1.05
  )
}

plot_gcm_agreement <- function(agreement, sd_summary) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  ssps <- c("SSP126", "SSP245", "SSP370", "SSP585")

  y1 <- agreement$ensemble_mean_suitable_km2 / 1000000
  y2 <- agreement$agreement_3of3_km2 / 1000000
  y_lim <- range(c(y1, y2), finite = TRUE)
  y_pad <- diff(y_lim) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  par(mfrow = c(2, 2), mar = c(4.2, 4.6, 3.0, 4.6), oma = c(0, 0, 1.5, 0))

  for (ssp in ssps) {
    a <- agreement[agreement$ssp == ssp, , drop = FALSE]
    s <- sd_summary[sd_summary$ssp == ssp, , drop = FALSE]
    a <- a[order(a$period_start), , drop = FALSE]
    s <- s[order(s$period_start), , drop = FALSE]

    x <- period_midpoint(a$period_start, a$period_end)

    plot(
      x,
      a$ensemble_mean_suitable_km2 / 1000000,
      type = "b",
      pch = 16,
      lty = 1,
      lwd = 2,
      ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad),
      xlab = "Projection period midpoint",
      ylab = expression("Area (" * 10^6 * " km"^2 * ")"),
      main = ssp,
      axes = FALSE
    )

    axis(1, at = x, labels = paste0(a$period_start, "-", a$period_end),
         las = 2, cex.axis = 0.70)
    axis(2, las = 1)
    box()
    grid(col = "grey85")

    lines(
      x,
      a$agreement_3of3_km2 / 1000000,
      type = "b",
      pch = 1,
      lty = 2,
      lwd = 2
    )

    par(new = TRUE)
    plot(
      x,
      s$inter_GCM_SD_mean,
      type = "b",
      pch = 17,
      lty = 3,
      lwd = 1.8,
      axes = FALSE,
      xlab = "",
      ylab = "",
      ylim = range(sd_summary$inter_GCM_SD_mean, finite = TRUE)
    )
    axis(4, las = 1, cex.axis = 0.75)
    mtext("Mean inter-GCM SD", side = 4, line = 3.1, cex = 0.75)

    legend(
      "topleft",
      legend = c("Ensemble suitable", "3/3 GCM suitable", "Mean inter-GCM SD"),
      lty = c(1, 2, 3),
      pch = c(16, 1, 17),
      lwd = c(2, 2, 1.8),
      bty = "n",
      cex = 0.70
    )
  }

  mtext(
    "Inter-GCM agreement and uncertainty in future suitability",
    side = 3,
    outer = TRUE,
    font = 2,
    cex = 1.05
  )
}

plot_clamping_confidence <- function(clamping) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  ssps <- c("SSP126", "SSP245", "SSP370", "SSP585")

  y <- clamping$strict_high_confidence_unclamped_km2 / 1000000
  y2 <- clamping$relaxed_high_confidence_clamp_at_most_1_km2 / 1000000
  y_lim <- range(c(y, y2), finite = TRUE)
  y_pad <- diff(y_lim) * 0.08
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1

  par(mfrow = c(2, 2), mar = c(4.2, 4.6, 3.0, 4.6), oma = c(0, 0, 1.5, 0))

  for (ssp in ssps) {
    d <- clamping[clamping$ssp == ssp, , drop = FALSE]
    d <- d[order(d$period_start), , drop = FALSE]

    x <- period_midpoint(d$period_start, d$period_end)

    plot(
      x,
      d$strict_high_confidence_unclamped_km2 / 1000000,
      type = "b",
      pch = 16,
      lty = 1,
      lwd = 2,
      ylim = c(y_lim[1] - y_pad, y_lim[2] + y_pad),
      xlab = "Projection period midpoint",
      ylab = expression("High-confidence area (" * 10^6 * " km"^2 * ")"),
      main = ssp,
      axes = FALSE
    )

    axis(1, at = x, labels = paste0(d$period_start, "-", d$period_end),
         las = 2, cex.axis = 0.70)
    axis(2, las = 1)
    box()
    grid(col = "grey85")

    lines(
      x,
      d$relaxed_high_confidence_clamp_at_most_1_km2 / 1000000,
      type = "b",
      pch = 1,
      lty = 2,
      lwd = 2
    )

    par(new = TRUE)
    plot(
      x,
      d$strict_percent_of_3of3_agreement_suitable,
      type = "b",
      pch = 17,
      lty = 3,
      lwd = 1.8,
      axes = FALSE,
      xlab = "",
      ylab = "",
      ylim = c(95, 100)
    )
    axis(4, las = 1, cex.axis = 0.75)
    mtext("Strict share of 3/3 habitat (%)", side = 4, line = 3.1, cex = 0.75)

    legend(
      "bottomleft",
      legend = c("Strict: 3/3 + no clamping", "Relaxed: 3/3 + clamping ≤1", "Strict share"),
      lty = c(1, 2, 3),
      pch = c(16, 1, 17),
      lwd = c(2, 2, 1.8),
      bty = "n",
      cex = 0.66
    )
  }

  mtext(
    "Clamping-aware high-confidence future suitability",
    side = 3,
    outer = TRUE,
    font = 2,
    cex = 1.05
  )
}

plot_priority_bar <- function(priority_area) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  d <- priority_area[priority_area$code %in% 1:5, , drop = FALSE]
  d <- d[order(d$code), , drop = FALSE]

  labels <- c(
    "Emerging\nfuture suitable",
    "Ex situ\ngermplasm",
    "Secondary\nmonitoring",
    "High-priority\nin situ",
    "Core\nrefugia"
  )

  par(mar = c(7.8, 5.0, 3.5, 1.2))

  heights <- d$area_km2 / 1000

  bars <- barplot(
    heights,
    names.arg = labels,
    las = 2,
    ylab = expression("Area (" * 10^3 * " km"^2 * ")"),
    main = "Native-range conservation priority areas under the strict rule",
    border = "grey30"
  )

  grid(nx = NA, ny = NULL, col = "grey85")
  barplot(
    heights,
    names.arg = labels,
    las = 2,
    ylab = expression("Area (" * 10^3 * " km"^2 * ")"),
    main = "Native-range conservation priority areas under the strict rule",
    border = "grey30",
    add = TRUE
  )

  text(
    bars,
    heights,
    labels = format_km2(d$area_km2, digits = 0),
    pos = 3,
    cex = 0.75
  )

  mtext(
    "Strict rule = current/future MTSPS suitability filtered by 3/3 GCM agreement and no clamping",
    side = 1,
    line = 6.6,
    cex = 0.75
  )
}

plot_priority_map_panel <- function(current_native, persistence_native, priority_native) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(mfrow = c(1, 3), mar = c(2.8, 2.8, 3.0, 1.0))

  plot(
    current_native,
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
    cex = 0.70
  )

  plot(
    persistence_native,
    col = hcl.colors(17, "YlOrRd", rev = FALSE),
    breaks = seq(-0.5, 16.5, by = 1),
    plg = FALSE,
    axes = TRUE,
    main = "B. Strict persistence"
  )
  legend(
    "bottomleft",
    legend = c("0", "4", "8", "12", "16"),
    title = "Scenarios",
    bty = "n",
    cex = 0.70
  )

  plot(
    priority_native,
    col = c("#F0F0F0", "#A6CEE3", "#E31A1C", "#FDBF6F", "#33A02C", "#006D2C"),
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5),
    plg = FALSE,
    axes = TRUE,
    main = "C. Priority class"
  )
  legend(
    "bottomleft",
    legend = c(
      "No priority",
      "Emerging",
      "Ex situ",
      "Secondary",
      "High in situ",
      "Core refugia"
    ),
    fill = c("#F0F0F0", "#A6CEE3", "#E31A1C", "#FDBF6F", "#33A02C", "#006D2C"),
    bty = "n",
    cex = 0.58
  )
}

# =============================================================================
# F. MAIN SCRIPT
# =============================================================================

run_script20 <- function() {

  for (folder in c(
    output_root, status_dir, table_dir, figure_dir, supplement_dir,
    inventory_dir, log_dir
  )) {
    ensure_dir(folder)
  }

  unlink(completed_marker, force = TRUE)
  unlink(failed_marker, force = TRUE)

  writeLines(
    c(
      "Script 20 started",
      paste("Started:", format(script_started_at, "%Y-%m-%d %H:%M:%S")),
      "Purpose: final manuscript tables and figures"
    ),
    started_marker
  )

  if (!requireNamespace("terra", quietly = TRUE)) {
    install.packages("terra", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(terra))

  write_log("Package terra loaded.")

  required_files <- c(
    script17_future_area_csv,
    script17_threshold_summary_csv,
    script17_threshold_values_csv,
    script17_occurrence_prediction_summary_csv,
    script18_agreement_area_csv,
    script18_sd_summary_csv,
    script18_high_conf_transition_csv,
    script18_scenario_confidence_csv,
    script18_extrapolation_inventory_csv,
    script18B_clamping_area_csv,
    script18B_priority_readiness_csv,
    script19_native_summary_csv,
    script19_strict_priority_area_csv,
    script19_binary_priority_area_csv,
    script19_strict_vs_relaxed_csv,
    script19_manuscript_summary_csv,
    script19_current_native_tif,
    script19_strict_persistence_native_tif,
    script19_strict_priority_native_tif
  )

  check_required_files(required_files)

  # -------------------------------------------------------------------
  # F1. Read verified outputs
  # -------------------------------------------------------------------

  future_area <- read_required_csv(script17_future_area_csv, "Script 17 future area table")
  threshold_summary <- read_required_csv(script17_threshold_summary_csv, "Script 17 threshold summary table")
  threshold_values <- read_required_csv(script17_threshold_values_csv, "Script 17 threshold values table")
  occurrence_prediction_summary <- read_required_csv(script17_occurrence_prediction_summary_csv, "Script 17 occurrence prediction summary")
  agreement_area <- read_required_csv(script18_agreement_area_csv, "Script 18 GCM agreement area table")
  sd_summary <- read_required_csv(script18_sd_summary_csv, "Script 18 inter-GCM SD summary")
  high_conf_transition <- read_required_csv(script18_high_conf_transition_csv, "Script 18 high-confidence transition table")
  scenario_confidence <- read_required_csv(script18_scenario_confidence_csv, "Script 18 scenario-confidence table")
  extrapolation_inventory <- read_required_csv(script18_extrapolation_inventory_csv, "Script 18 extrapolation inventory")
  clamping_area <- read_required_csv(script18B_clamping_area_csv, "Script 18B clamping area table")
  priority_readiness <- read_required_csv(script18B_priority_readiness_csv, "Script 18B priority readiness table")
  native_summary <- read_required_csv(script19_native_summary_csv, "Script 19 native summary")
  strict_priority_area <- read_required_csv(script19_strict_priority_area_csv, "Script 19 strict priority area")
  binary_priority_area <- read_required_csv(script19_binary_priority_area_csv, "Script 19 binary priority area")
  strict_vs_relaxed <- read_required_csv(script19_strict_vs_relaxed_csv, "Script 19 strict-relaxed comparison")
  manuscript_summary <- read_required_csv(script19_manuscript_summary_csv, "Script 19 manuscript summary")

  # Ensure numeric fields.
  for (data_name in c("future_area", "agreement_area", "sd_summary", "clamping_area")) {
    data_object <- get(data_name)
    for (column in names(data_object)) {
      if (grepl("km2|percent|SD|threshold|period", column, ignore.case = TRUE)) {
        candidate <- suppressWarnings(as.numeric(data_object[[column]]))
        if (sum(is.finite(candidate)) > 0) data_object[[column]] <- candidate
      }
    }
    assign(data_name, data_object)
  }

  # -------------------------------------------------------------------
  # F2. Manuscript Table 1: model setup and thresholds
  # -------------------------------------------------------------------

  mtsps_value <- threshold_values$threshold_cloglog[
    threshold_values$threshold_id == "MTSPS"
  ]
  mtp_value <- threshold_values$threshold_cloglog[
    threshold_values$threshold_id == "MTP"
  ]
  t10_value <- threshold_values$threshold_cloglog[
    threshold_values$threshold_id == "T10"
  ]

  if (length(mtsps_value) != 1L) mtsps_value <- NA_real_
  if (length(mtp_value) != 1L) mtp_value <- NA_real_
  if (length(t10_value) != 1L) t10_value <- NA_real_

  model_setup <- data.frame(
    item = c(
      "Final occurrence records",
      "Selected predictors",
      "MaxEnt feature class",
      "MaxEnt regularization multiplier",
      "Primary binary threshold",
      "MTSPS threshold",
      "Supplementary threshold: T10",
      "Supplementary threshold: MTP",
      "Future GCMs",
      "Future scenarios",
      "Future periods",
      "Primary future-confidence rule",
      "Conservation-priority spatial restriction"
    ),
    value = c(
      "828",
      "Bio2; Bio6; Bio14; Bio18; Elv",
      "lph",
      "1",
      "MTSPS",
      format(mtsps_value, scientific = TRUE, digits = 8),
      format(t10_value, scientific = TRUE, digits = 8),
      format(mtp_value, scientific = TRUE, digits = 8),
      "ACCESS-CM2; BCC-CSM2-MR; CMCC-ESM2",
      "SSP126; SSP245; SSP370; SSP585",
      "2021-2040; 2041-2060; 2061-2080; 2081-2100",
      "3/3 GCM suitable + no clamping",
      "native_wild_range_union"
    ),
    interpretation = c(
      "Final globally screened occurrence-candidate dataset",
      "Final low-collinearity climate-topographic predictor set",
      "Best calibrated MaxEnt feature combination",
      "Best calibrated regularization strength",
      "Used for main binary maps and future comparisons",
      "Primary threshold applied to current and future projections",
      "Used as supplementary threshold-sensitivity check",
      "Detected as zero and retained only as an extreme lower-bound check",
      "Three-GCM uncertainty framework",
      "Four CMIP6 SSP pathways",
      "Four future time slices",
      "Used for strict climate-stable refugia and priority analysis",
      "Main conservation-priority interpretation is native-range only"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(model_setup, table_model_setup_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # F3. Manuscript Table 2: threshold sensitivity
  # -------------------------------------------------------------------

  threshold_main <- threshold_summary[
    threshold_summary$threshold_id %in% c("MTSPS", "T10", "MTP"),
    ,
    drop = FALSE
  ]

  threshold_sensitivity <- threshold_main[
    ,
    c(
      "threshold_id",
      "threshold_cloglog",
      "current_suitable_km2",
      "current_percent_suitable",
      "minimum_future_suitable_km2",
      "maximum_future_suitable_km2",
      "minimum_net_change_km2",
      "maximum_net_change_km2",
      "scenarios_with_net_gain",
      "scenarios_with_net_loss"
    ),
    drop = FALSE
  ]

  threshold_sensitivity$note <- ifelse(
    threshold_sensitivity$threshold_id == "MTP",
    "Degenerate zero-threshold result; do not use for substantive spatial interpretation.",
    "Retained for substantive threshold-sensitivity interpretation."
  )

  write.csv(threshold_sensitivity, table_threshold_sensitivity_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # F4. Manuscript Table 3: GCM agreement and uncertainty
  # -------------------------------------------------------------------

  uncertainty_summary <- merge(
    agreement_area[
      ,
      c(
        "scenario",
        "ssp",
        "period_start",
        "period_end",
        "ensemble_mean_suitable_km2",
        "agreement_3of3_km2",
        "ensemble_suitable_but_not_3of3_km2",
        "percent_ensemble_suitable_with_3of3_agreement"
      ),
      drop = FALSE
    ],
    sd_summary[
      ,
      c(
        "scenario",
        "inter_GCM_SD_mean",
        "inter_GCM_SD_max",
        "mean_SD_agreement_3of3"
      ),
      drop = FALSE
    ],
    by = "scenario",
    all.x = TRUE,
    sort = FALSE
  )

  uncertainty_summary <- uncertainty_summary[
    order(uncertainty_summary$ssp, uncertainty_summary$period_start),
    ,
    drop = FALSE
  ]

  write.csv(uncertainty_summary, table_uncertainty_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # F5. Manuscript Table 4: native conservation priority
  # -------------------------------------------------------------------

  strict_priority_for_manuscript <- strict_priority_area[
    strict_priority_area$code %in% 1:5,
    ,
    drop = FALSE
  ]

  priority_table <- strict_priority_for_manuscript[
    ,
    c("code", "class_name", "area_km2", "area_10_4_km2", "percent_native_range"),
    drop = FALSE
  ]

  priority_table$future_rule <- "strict: 3/3 GCM suitable + no clamping"
  priority_table$native_range_restriction <- "native_wild_range_union"

  write.csv(priority_table, table_priority_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # F6. Key numbers for Results text
  # -------------------------------------------------------------------

  get_metric_area <- function(metric_name) {
    value <- manuscript_summary$area_km2[manuscript_summary$metric == metric_name]
    if (length(value) != 1L) return(NA_real_)
    safe_number(value)
  }

  native_area <- get_metric_area("Native range area within strict common domain")
  current_native <- get_metric_area("Current MTSPS suitable native range")
  strict_core <- get_metric_area("Strict core climate-stable refugia")
  strict_high <- get_metric_area("Strict high-priority in situ zones including core refugia")
  strict_exsitu <- get_metric_area("Strict ex situ germplasm collection priority")
  strict_emerging <- get_metric_area("Strict emerging future suitable native areas")

  mtsps_threshold_row <- threshold_sensitivity[
    threshold_sensitivity$threshold_id == "MTSPS",
    ,
    drop = FALSE
  ]

  t10_threshold_row <- threshold_sensitivity[
    threshold_sensitivity$threshold_id == "T10",
    ,
    drop = FALSE
  ]

  key_numbers <- data.frame(
    result_item = c(
      "MTSPS future suitable area range",
      "MTSPS net suitable-area change range",
      "T10 future suitable area range",
      "T10 net suitable-area change range",
      "Lowest future ensemble-suitable area",
      "Highest future ensemble-suitable area",
      "Lowest 3/3-GCM suitable area",
      "Highest 3/3-GCM suitable area",
      "Lowest percent of ensemble suitable area with 3/3 agreement",
      "Highest percent of ensemble suitable area with 3/3 agreement",
      "Native range area",
      "Current suitable native range",
      "Strict core climate-stable refugia",
      "Strict high-priority in situ zones including core refugia",
      "Strict ex situ germplasm collection priority",
      "Strict emerging future suitable native areas"
    ),
    value = c(
      paste0(
        format_km2(mtsps_threshold_row$minimum_future_suitable_km2),
        "–",
        format_km2(mtsps_threshold_row$maximum_future_suitable_km2),
        " km²"
      ),
      paste0(
        format_km2(mtsps_threshold_row$minimum_net_change_km2),
        "–",
        format_km2(mtsps_threshold_row$maximum_net_change_km2),
        " km²"
      ),
      paste0(
        format_km2(t10_threshold_row$minimum_future_suitable_km2),
        "–",
        format_km2(t10_threshold_row$maximum_future_suitable_km2),
        " km²"
      ),
      paste0(
        format_km2(t10_threshold_row$minimum_net_change_km2),
        "–",
        format_km2(t10_threshold_row$maximum_net_change_km2),
        " km²"
      ),
      paste0(format_km2(min(agreement_area$ensemble_mean_suitable_km2)), " km²"),
      paste0(format_km2(max(agreement_area$ensemble_mean_suitable_km2)), " km²"),
      paste0(format_km2(min(agreement_area$agreement_3of3_km2)), " km²"),
      paste0(format_km2(max(agreement_area$agreement_3of3_km2)), " km²"),
      format_percent(min(agreement_area$percent_ensemble_suitable_with_3of3_agreement)),
      format_percent(max(agreement_area$percent_ensemble_suitable_with_3of3_agreement)),
      paste0(format_km2(native_area), " km²"),
      paste0(format_km2(current_native), " km²"),
      paste0(format_km2(strict_core), " km²"),
      paste0(format_km2(strict_high), " km²"),
      paste0(format_km2(strict_exsitu), " km²"),
      paste0(format_km2(strict_emerging), " km²")
    ),
    suggested_results_sentence = c(
      "Under the primary MTSPS threshold, future ensemble suitable area ranged from the lower to upper value reported here.",
      "All MTSPS scenarios showed net suitable-area gain relative to the current baseline.",
      "The T10 sensitivity threshold produced a comparable range of future suitable area.",
      "The T10 sensitivity threshold preserved the direction of net suitable-area change.",
      "This is the minimum ensemble-suitable area across the 16 future projections.",
      "This is the maximum ensemble-suitable area across the 16 future projections.",
      "This is the minimum area supported by full 3/3 GCM agreement.",
      "This is the maximum area supported by full 3/3 GCM agreement.",
      "This lower percentage shows why uncertainty-aware interpretation is required.",
      "This upper percentage shows the maximum 3/3 agreement share of ensemble-suitable habitat.",
      "Conservation-priority interpretation was restricted to this native-wild range area.",
      "This is the currently suitable native-range area under MTSPS.",
      "These cells remained strict suitable across all 16 future projections.",
      "These cells retained strict suitability in at least 12 of 16 projections.",
      "These currently suitable cells had low strict future persistence and are suitable for ex situ collection priority.",
      "These are current-unsuitable native-range cells projected to gain strict future suitability in at least 12 of 16 projections."
    ),
    stringsAsFactors = FALSE
  )

  write.csv(key_numbers, table_manuscript_key_numbers_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # F7. Supplementary table copies
  # -------------------------------------------------------------------

  write.csv(future_area, supplement_threshold_long_csv, row.names = FALSE)
  write.csv(agreement_area, supplement_agreement_long_csv, row.names = FALSE)
  write.csv(clamping_area, supplement_clamping_long_csv, row.names = FALSE)
  write.csv(strict_vs_relaxed, supplement_priority_strict_relaxed_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # F8. Figures
  # -------------------------------------------------------------------

  write_log("Creating manuscript figures.")

  save_base_plot_dual(
    pdf_path = figure_future_area_pdf,
    png_path = sub("\\.pdf$", ".png", figure_future_area_pdf),
    plot_function = function() plot_future_area(future_area),
    width = figure_width,
    height = figure_height
  )

  save_base_plot_dual(
    pdf_path = figure_gcm_agreement_pdf,
    png_path = sub("\\.pdf$", ".png", figure_gcm_agreement_pdf),
    plot_function = function() plot_gcm_agreement(agreement_area, sd_summary),
    width = figure_width,
    height = figure_height
  )

  save_base_plot_dual(
    pdf_path = figure_clamping_pdf,
    png_path = sub("\\.pdf$", ".png", figure_clamping_pdf),
    plot_function = function() plot_clamping_confidence(clamping_area),
    width = figure_width,
    height = figure_height
  )

  save_base_plot_dual(
    pdf_path = figure_priority_bar_pdf,
    png_path = sub("\\.pdf$", ".png", figure_priority_bar_pdf),
    plot_function = function() plot_priority_bar(strict_priority_area),
    width = figure_width,
    height = figure_height
  )

  current_native_raster <- terra::rast(script19_current_native_tif)
  persistence_native_raster <- terra::rast(script19_strict_persistence_native_tif)
  priority_native_raster <- terra::rast(script19_strict_priority_native_tif)

  save_base_plot_dual(
    pdf_path = figure_priority_map_pdf,
    png_path = sub("\\.pdf$", ".png", figure_priority_map_pdf),
    plot_function = function() {
      plot_priority_map_panel(
        current_native = current_native_raster,
        persistence_native = persistence_native_raster,
        priority_native = priority_native_raster
      )
    },
    width = map_figure_width,
    height = map_figure_height
  )

  # -------------------------------------------------------------------
  # F9. Inventories
  # -------------------------------------------------------------------

  figure_inventory <- data.frame(
    figure_id = c("Figure 1", "Figure 2", "Figure 3", "Figure 4", "Figure 5"),
    figure_file_pdf = c(
      figure_future_area_pdf,
      figure_gcm_agreement_pdf,
      figure_clamping_pdf,
      figure_priority_bar_pdf,
      figure_priority_map_pdf
    ),
    recommended_placement = c(
      "Main text or Supplementary; useful for reviewer response on threshold sensitivity",
      "Main text; uncertainty analysis and GCM agreement",
      "Supplementary or Main text if reviewer emphasis on extrapolation is strong",
      "Main text; conservation priority Results",
      "Main text; final conservation priority map panel"
    ),
    main_message = c(
      "MTSPS and T10 preserve the same direction of future suitable-area change.",
      "Only part of ensemble-suitable area is supported by full 3/3 GCM agreement.",
      "Most 3/3-GCM suitable habitat remains unclamped; strict rule is defensible.",
      "Native-range strict priority areas are dominated by core refugia and high-priority in situ zones.",
      "Strict future persistence and integrated priority classes are spatially explicit within the native range."
    ),
    stringsAsFactors = FALSE
  )

  figure_inventory$figure_file_pdf <- normalizePath(
    figure_inventory$figure_file_pdf,
    winslash = "/",
    mustWork = TRUE
  )

  if (write_png_copies) {
    figure_inventory$figure_file_png <- normalizePath(
      sub("\\.pdf$", ".png", figure_inventory$figure_file_pdf),
      winslash = "/",
      mustWork = FALSE
    )
  }

  write.csv(figure_inventory, figure_inventory_csv, row.names = FALSE)

  table_inventory <- data.frame(
    table_id = c(
      "Table 1",
      "Table 2",
      "Table 3",
      "Table 4",
      "Table 5",
      "Table S1",
      "Table S2",
      "Table S3",
      "Table S4"
    ),
    table_file = c(
      table_model_setup_csv,
      table_threshold_sensitivity_csv,
      table_uncertainty_csv,
      table_priority_csv,
      table_manuscript_key_numbers_csv,
      supplement_threshold_long_csv,
      supplement_agreement_long_csv,
      supplement_clamping_long_csv,
      supplement_priority_strict_relaxed_csv
    ),
    recommended_placement = c(
      "Methods or Supplementary",
      "Results or Supplementary",
      "Results",
      "Results",
      "Internal writing aid; not necessarily for submission",
      "Supplementary",
      "Supplementary",
      "Supplementary",
      "Supplementary"
    ),
    purpose = c(
      "Documents final occurrence, predictor, tuning and threshold settings.",
      "Shows threshold sensitivity and MTP degeneracy.",
      "Summarises GCM agreement and inter-GCM SD.",
      "Reports native-range conservation priority areas.",
      "Provides exact numbers for Results text.",
      "Full scenario-level threshold-sensitivity table.",
      "Full scenario-level GCM agreement table.",
      "Full scenario-level clamping-confidence table.",
      "Strict versus relaxed conservation-priority sensitivity."
    ),
    stringsAsFactors = FALSE
  )

  table_inventory$table_file <- normalizePath(
    table_inventory$table_file,
    winslash = "/",
    mustWork = TRUE
  )

  write.csv(table_inventory, table_inventory_csv, row.names = FALSE)

  # -------------------------------------------------------------------
  # F10. Run summary
  # -------------------------------------------------------------------

  run_summary <- data.frame(
    species = species_name,
    script = "20_manuscript_tables_and_figures",
    completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    manuscript_tables_created = 5L,
    supplementary_tables_created = 4L,
    figures_created = 5L,
    output_root = normalizePath(output_root, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )

  write.csv(run_summary, run_summary_csv, row.names = FALSE)

  writeLines(
    c(
      "SCRIPT 20 COMPLETED SUCCESSFULLY",
      paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      paste("Output root:", normalizePath(output_root, winslash = "/", mustWork = TRUE)),
      paste("Main tables:", table_dir),
      paste("Figures:", figure_dir),
      paste("Supplementary tables:", supplement_dir)
    ),
    completed_marker
  )

  write_log("Script completed successfully.")

  cat("\n======================================================================\n")
  cat("SCRIPT 20 COMPLETED — MANUSCRIPT TABLES AND FIGURES\n")
  cat("======================================================================\n")
  cat("Main manuscript tables:\n")
  cat(table_model_setup_csv, "\n")
  cat(table_threshold_sensitivity_csv, "\n")
  cat(table_uncertainty_csv, "\n")
  cat(table_priority_csv, "\n")
  cat(table_manuscript_key_numbers_csv, "\n")
  cat("\nFigures:\n")
  cat(figure_future_area_pdf, "\n")
  cat(figure_gcm_agreement_pdf, "\n")
  cat(figure_clamping_pdf, "\n")
  cat(figure_priority_bar_pdf, "\n")
  cat(figure_priority_map_pdf, "\n")
  cat("\nInventories:\n")
  cat(figure_inventory_csv, "\n")
  cat(table_inventory_csv, "\n")
  cat("\nOutput folder:\n")
  cat(output_root, "\n")

  invisible(list(
    model_setup = model_setup,
    threshold_sensitivity = threshold_sensitivity,
    uncertainty_summary = uncertainty_summary,
    priority_table = priority_table,
    key_numbers = key_numbers,
    figure_inventory = figure_inventory,
    table_inventory = table_inventory
  ))
}

# =============================================================================
# G. EXECUTE WITH EXPLICIT FAILURE REPORTING
# =============================================================================

tryCatch(
  {
    run_script20()
  },
  error = function(error) {

    dir.create(status_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

    failure_message <- paste0(
      "SCRIPT 20 FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(error), "\n\n",
      "Execution log: ", execution_log, "\n"
    )

    writeLines(failure_message, failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 20 FAILED\n")
    cat("======================================================================\n")
    cat(failure_message, "\n")

    stop(error)
  }
)
