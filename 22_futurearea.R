# ============================================================
# FIGURE 7
# Future suitable area under MTSPS and T10 thresholds
# Camellia japonica
# Improved layout and Times New Roman styling
# ============================================================

rm(list = ls())

# -----------------------------
# 1. Load packages
# -----------------------------
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(dplyr)

# -----------------------------
# 2. Define Times New Roman font for Windows
# -----------------------------
windowsFonts(TNR = windowsFont("Times New Roman"))

# -----------------------------
# 3. Input and output paths
# -----------------------------
input_csv <- "E:/PREDICTION TREE/Camellia_japonica/global_reanalysis_outputs/20_manuscript_tables_and_figures/03_supplementary_tables/Table_S_threshold_sensitivity_all_scenarios.csv"

output_dir <- "E:/PREDICTION TREE/Camellia_japonica/global_reanalysis_outputs/20_manuscript_tables_and_figures/02_manuscript_figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

output_png <- file.path(output_dir, "Figure_7_future_suitable_area_MTSPS_T10_10_4km2_TNR.png")
output_pdf <- file.path(output_dir, "Figure_7_future_suitable_area_MTSPS_T10_10_4km2_TNR.pdf")

# -----------------------------
# 4. Read data
# -----------------------------
if (!file.exists(input_csv)) {
  stop("Input file not found:\n", input_csv)
}

df <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

# -----------------------------
# 5. Prepare data
# -----------------------------
plot_df <- df %>%
  filter(threshold_id %in% c("MTSPS", "T10")) %>%
  mutate(
    ssp = factor(ssp, levels = c("SSP126", "SSP245", "SSP370", "SSP585")),
    period_label = paste0(period_start, "\u2013", period_end),
    period_label = factor(
      period_label,
      levels = c("2021–2040", "2041–2060", "2061–2080", "2081–2100")
    ),
    area_10_4 = as.numeric(future_suitable_10_4_km2),
    label_text = sprintf("%.3f", round(area_10_4, 3))
  ) %>%
  arrange(ssp, period_start, threshold_id)

if (nrow(plot_df) == 0) {
  stop("No rows found for MTSPS and T10.")
}

# -----------------------------
# 6. Position settings
# -----------------------------
period_levels <- c("2021–2040", "2041–2060", "2061–2080", "2081–2100")
plot_df$x_base <- match(plot_df$period_label, period_levels)

# Slight horizontal separation between the two thresholds
plot_df$x_plot <- ifelse(
  plot_df$threshold_id == "MTSPS",
  plot_df$x_base - 0.06,
  plot_df$x_base + 0.06
)

# Larger vertical separation so labels are above points
y_span <- max(plot_df$area_10_4, na.rm = TRUE) - min(plot_df$area_10_4, na.rm = TRUE)
if (y_span == 0) y_span <- 1

plot_df$label_y <- ifelse(
  plot_df$threshold_id == "MTSPS",
  plot_df$area_10_4 + y_span * 0.045,
  plot_df$area_10_4 + y_span * 0.030
)

# -----------------------------
# 7. Create plot
# -----------------------------
p <- ggplot() +

  geom_line(
    data = plot_df,
    aes(
      x = x_plot,
      y = area_10_4,
      group = threshold_id,
      linetype = threshold_id
    ),
    colour = "black",
    linewidth = 0.9
  ) +

  geom_point(
    data = plot_df,
    aes(
      x = x_plot,
      y = area_10_4,
      shape = threshold_id,
      fill = threshold_id
    ),
    colour = "black",
    size = 3.8,
    stroke = 1.0
  ) +

  geom_text(
    data = plot_df,
    aes(
      x = x_plot,
      y = label_y,
      label = label_text
    ),
    family = "TNR",
    size = 4.6,
    colour = "black",
    vjust = 0
  ) +

  facet_wrap(~ ssp, ncol = 2) +

  scale_x_continuous(
    breaks = 1:4,
    labels = period_levels,
    expand = expansion(mult = c(0.10, 0.10))
  ) +

  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.20))
  ) +

  scale_linetype_manual(
    values = c("MTSPS" = "solid", "T10" = "dashed"),
    name = NULL
  ) +

  scale_shape_manual(
    values = c("MTSPS" = 21, "T10" = 21),
    name = NULL
  ) +

  scale_fill_manual(
    values = c("MTSPS" = "black", "T10" = "white"),
    name = NULL
  ) +

  labs(
    title = "Future climatically suitable area under MTSPS and T10 thresholds",
    x = "Projection period",
    y = expression("Suitable area (" * 10^4 * " km"^2 * ")")
  ) +

  coord_cartesian(clip = "off") +

  theme_bw(base_family = "TNR", base_size = 16) +
  theme(
    text = element_text(family = "TNR", colour = "black"),
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 18)
    ),
    axis.title.x = element_text(
      size = 18,
      face = "bold",
      margin = margin(t = 14)
    ),
    axis.title.y = element_text(
      size = 18,
      face = "bold",
      margin = margin(r = 14)
    ),
    axis.text.x = element_text(
      size = 14,
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.text.y = element_text(size = 14),
    strip.text = element_text(
      size = 18,
      face = "bold"
    ),
    legend.position = "top",
    legend.text = element_text(size = 15),
    legend.key.width = unit(1.2, "cm"),
    panel.grid.major = element_line(colour = "grey80", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.4, "lines"),
    plot.margin = margin(20, 28, 20, 20)
  )

# -----------------------------
# 8. Save figure
# -----------------------------
ggsave(
  filename = output_png,
  plot = p,
  width = 15,
  height = 10.5,
  dpi = 700,
  bg = "white"
)

ggsave(
  filename = output_pdf,
  plot = p,
  width = 15,
  height = 10.5,
  bg = "white"
)

# -----------------------------
# 9. Message
# -----------------------------
cat("\n============================================================\n")
cat("FIGURE 7 CREATED SUCCESSFULLY\n")
cat("============================================================\n")
cat("Input file:\n", input_csv, "\n\n", sep = "")
cat("Output PNG:\n", output_png, "\n\n", sep = "")
cat("Output PDF:\n", output_pdf, "\n\n", sep = "")
cat("============================================================\n")