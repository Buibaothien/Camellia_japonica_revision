# =====================================================================
# REDRAW KUENM CALIBRATION FIGURE — Camellia japonica
# Global workflow: OR10, 31 FC x 8 RM
# =====================================================================
#
# This script DOES NOT rerun kuenm.
# It reads the existing kuenm_ceval output files:
#   calibration_results.csv
#   selected_models.csv
#
# It automatically:
# - uses the global Camellia japonica kuenm workspace;
# - detects the actual omission-rate column (expected OR10);
# - identifies selected models from selected_models.csv;
# - identifies non-significant partial-ROC models when p-value data exist;
# - expands the X-axis based on the actual log(AICc) range;
# - exports a high-resolution PNG and editable PDF.
# =====================================================================

# ---------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

# This must match the workspace created by:
# kuenm_Camellia_japonica_global_OR10_like_Crinum.R
workdir <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "09_kuenm_global_OR10"
)

eval_dir <- file.path(workdir, "05_calibration_results")

all_models_csv <- file.path(eval_dir, "calibration_results.csv")
selected_models_csv <- file.path(eval_dir, "selected_models.csv")

output_png <- file.path(
  eval_dir,
  "Camellia_japonica_global_OR10_kuenm_calibration_extended_x.png"
)

output_pdf <- file.path(
  eval_dir,
  "Camellia_japonica_global_OR10_kuenm_calibration_extended_x.pdf"
)

if (!file.exists(all_models_csv)) {
  stop(
    "calibration_results.csv was not found:\n",
    all_models_csv,
    "\n\nRun kuenm_ceval successfully before drawing the figure."
  )
}

if (!file.exists(selected_models_csv)) {
  stop(
    "selected_models.csv was not found:\n",
    selected_models_csv,
    "\n\nRun kuenm_ceval successfully before drawing the figure."
  )
}

# ---------------------------------------------------------------------
# 2. FIGURE SETTINGS
# ---------------------------------------------------------------------

# Additional empty space on the left and right of the actual log(AICc) range.
# Increase these values to make the X-axis longer.
x_padding_left <- 0.060
x_padding_right <- 0.060

# Leave NULL for automatic X-axis limits.
# To force a specific range, use for example:
# manual_x_limits <- c(8.80, 10.20)
manual_x_limits <- NULL

# Figure output specifications.
png_width_px <- 3200
png_height_px <- 2100
png_resolution_dpi <- 600

# ---------------------------------------------------------------------
# 3. HELPER FUNCTIONS
# ---------------------------------------------------------------------

normalise_name <- function(x) {
  tolower(gsub("[^a-z0-9]", "", x))
}

find_column <- function(data, candidates, object_label, required = TRUE) {
  original_names <- names(data)
  normalized_names <- normalise_name(original_names)
  normalized_candidates <- normalise_name(candidates)

  direct_matches <- match(normalized_candidates, normalized_names, nomatch = 0)
  direct_matches <- direct_matches[direct_matches > 0]

  if (length(direct_matches) > 0) {
    return(original_names[direct_matches[1]])
  }

  # Allow flexible matching for package-version differences.
  for (candidate in normalized_candidates) {
    partial_matches <- grep(candidate, normalized_names, fixed = TRUE)

    if (length(partial_matches) > 0) {
      return(original_names[partial_matches[1]])
    }
  }

  if (required) {
    stop(
      "Could not identify a required column in ", object_label, ".\n",
      "Available columns:\n",
      paste(original_names, collapse = ", ")
    )
  }

  NA_character_
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

# ---------------------------------------------------------------------
# 4. READ KUENM RESULTS
# ---------------------------------------------------------------------

all_models <- read.csv(
  all_models_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

selected_models <- read.csv(
  selected_models_csv,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

model_col <- find_column(
  all_models,
  c("Model", "model_name"),
  "calibration_results.csv"
)

aicc_col <- find_column(
  all_models,
  c("AICc"),
  "calibration_results.csv"
)

# The global Camellia workflow is OR10. The order below prioritizes OR10,
# but it also supports other kuenm output formats.
omission_col <- find_column(
  all_models,
  c(
    "Omission_rate_at_10%",
    "Omission_rate_at_10",
    "Omission_rate_10",
    "Omission_rate_at_5%",
    "Omission_rate_at_5",
    "Omission_rate"
  ),
  "calibration_results.csv"
)

pval_col <- find_column(
  all_models,
  c(
    "pval_pROC",
    "p_value_pROC",
    "partial_roc_p_value",
    "pROC_p_value"
  ),
  "calibration_results.csv",
  required = FALSE
)

selected_model_col <- find_column(
  selected_models,
  c("Model", "model_name"),
  "selected_models.csv",
  required = FALSE
)

# ---------------------------------------------------------------------
# 5. PREPARE PLOTTING TABLES
# ---------------------------------------------------------------------

all_models[[aicc_col]] <- safe_numeric(all_models[[aicc_col]])
all_models[[omission_col]] <- safe_numeric(all_models[[omission_col]])

all_models <- all_models[
  is.finite(all_models[[aicc_col]]) &
    is.finite(all_models[[omission_col]]) &
    all_models[[aicc_col]] > 0,
  ,
  drop = FALSE
]

if (nrow(all_models) == 0) {
  stop(
    "No valid rows remain after checking AICc and omission-rate values in calibration_results.csv."
  )
}

# Use selected_models.csv directly when it includes AICc and omission values.
# Otherwise, recover the selected rows through Model IDs.
if (
  aicc_col %in% names(selected_models) &&
  omission_col %in% names(selected_models)
) {
  selected_models[[aicc_col]] <- safe_numeric(selected_models[[aicc_col]])
  selected_models[[omission_col]] <- safe_numeric(selected_models[[omission_col]])

  selected_plot <- selected_models[
    is.finite(selected_models[[aicc_col]]) &
      is.finite(selected_models[[omission_col]]) &
      selected_models[[aicc_col]] > 0,
    ,
    drop = FALSE
  ]
} else {
  if (is.na(selected_model_col)) {
    stop(
      "selected_models.csv has neither matching AICc/omission columns nor a Model column."
    )
  }

  selected_ids <- as.character(selected_models[[selected_model_col]])

  selected_plot <- all_models[
    as.character(all_models[[model_col]]) %in% selected_ids,
    ,
    drop = FALSE
  ]
}

if (nrow(selected_plot) == 0) {
  warning(
    "No selected model rows could be matched. The figure will still be drawn, ",
    "but blue selected-model symbols will be absent."
  )
}

# For kuenm, p-value > 0.05 is treated as non-significant partial ROC.
if (!is.na(pval_col)) {
  all_models[[pval_col]] <- safe_numeric(all_models[[pval_col]])

  non_significant <- all_models[
    is.na(all_models[[pval_col]]) | all_models[[pval_col]] > 0.05,
    ,
    drop = FALSE
  ]
} else {
  non_significant <- all_models[0, , drop = FALSE]

  warning(
    "No partial-ROC p-value column was found. ",
    "The non-significant-model overlay will not be shown."
  )
}

x_all <- log(all_models[[aicc_col]])
y_all <- all_models[[omission_col]]

x_selected <- if (nrow(selected_plot) > 0) log(selected_plot[[aicc_col]]) else numeric(0)
y_selected <- if (nrow(selected_plot) > 0) selected_plot[[omission_col]] else numeric(0)

x_nonsig <- if (nrow(non_significant) > 0) log(non_significant[[aicc_col]]) else numeric(0)
y_nonsig <- if (nrow(non_significant) > 0) non_significant[[omission_col]] else numeric(0)

# ---------------------------------------------------------------------
# 6. AXIS SETTINGS
# ---------------------------------------------------------------------

if (is.null(manual_x_limits)) {
  x_limits <- c(
    min(x_all, na.rm = TRUE) - x_padding_left,
    max(x_all, na.rm = TRUE) + x_padding_right
  )
} else {
  x_limits <- manual_x_limits
}

y_range <- range(y_all, na.rm = TRUE)
y_span <- diff(y_range)
y_padding <- max(0.01, y_span * 0.08)

y_limits <- c(
  max(0, y_range[1] - y_padding),
  y_range[2] + y_padding
)

x_span <- diff(x_limits)

x_tick_step <- if (x_span <= 0.12) {
  0.02
} else if (x_span <= 0.30) {
  0.05
} else {
  0.10
}

x_ticks <- seq(
  from = ceiling(x_limits[1] / x_tick_step) * x_tick_step,
  to = floor(x_limits[2] / x_tick_step) * x_tick_step,
  by = x_tick_step
)

omission_number <- regmatches(
  omission_col,
  regexpr("[0-9]+", omission_col)
)

if (length(omission_number) == 0 || !nzchar(omission_number)) {
  omission_number <- "10"
}

y_axis_label <- paste0(
  "Omission rates at ",
  omission_number,
  "% threshold value"
)

# ---------------------------------------------------------------------
# 7. DRAW FIGURE
# ---------------------------------------------------------------------

draw_kuenm_plot <- function() {
  par(
    mar = c(5.7, 6.6, 1.0, 1.4),
    las = 1,
    cex.axis = 1.15,
    cex.lab = 1.30
  )

  plot(
    x_all,
    y_all,
    xlim = x_limits,
    ylim = y_limits,
    xlab = "Natural logarithm of AICc",
    ylab = y_axis_label,
    pch = 1,
    cex = 1.12,
    lwd = 1.40,
    col = "grey35",
    axes = FALSE
  )

  axis(
    side = 1,
    at = x_ticks,
    labels = format(x_ticks, nsmall = 2)
  )

  axis(side = 2)
  box(lwd = 1.50)

  if (nrow(non_significant) > 0) {
    points(
      x_nonsig,
      y_nonsig,
      pch = 19,
      cex = 1.03,
      col = "red2"
    )
  }

  if (nrow(selected_plot) > 0) {
    points(
      x_selected,
      y_selected,
      pch = 17,
      cex = 1.70,
      col = "dodgerblue3"
    )
  }

  legend(
    "topright",
    legend = c(
      "Selected models",
      "Non-significant models",
      "All candidate models"
    ),
    pch = c(17, 19, 1),
    pt.cex = c(1.50, 1.03, 1.12),
    col = c("dodgerblue3", "red2", "grey35"),
    bty = "n",
    cex = 0.95,
    inset = 0.018
  )
}

# ---------------------------------------------------------------------
# 8. EXPORT FIGURES
# ---------------------------------------------------------------------

png(
  filename = output_png,
  width = png_width_px,
  height = png_height_px,
  res = png_resolution_dpi
)
draw_kuenm_plot()
dev.off()

pdf(
  file = output_pdf,
  width = 11.0,
  height = 7.3,
  useDingbats = FALSE
)
draw_kuenm_plot()
dev.off()

# ---------------------------------------------------------------------
# 9. REPORT
# ---------------------------------------------------------------------

cat("\n======================================================================\n")
cat("CAMELLIA JAPONICA KUENM FIGURE CREATED\n")
cat("======================================================================\n")
cat("Input folder:\n", normalizePath(eval_dir, winslash = "/"), "\n", sep = "")
cat("PNG:\n", normalizePath(output_png, winslash = "/"), "\n", sep = "")
cat("PDF:\n", normalizePath(output_pdf, winslash = "/"), "\n", sep = "")

cat("\nX-axis limits used:\n")
print(x_limits)

cat("\nSelected model(s):\n")

columns_to_print <- unique(c(
  model_col,
  "Mean_AUC_ratio",
  pval_col,
  omission_col,
  aicc_col,
  "Delta_AICc",
  "delta_AICc",
  "W_AICc",
  "N_parameters",
  "num_parameters"
))

columns_to_print <- columns_to_print[
  !is.na(columns_to_print) &
    columns_to_print %in% names(selected_plot)
]

if (nrow(selected_plot) > 0 && length(columns_to_print) > 0) {
  print(selected_plot[, columns_to_print, drop = FALSE])
} else {
  cat("No selected rows were available for printing.\n")
}
