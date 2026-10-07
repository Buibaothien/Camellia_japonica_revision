# =============================================================================
# 18A_CREATE_FOCAL_UNCERTAINTY_FIGURES_CAMELLIA_JAPONICA.R
# =============================================================================
# PURPOSE
# -------
# Create the two focal uncertainty figures that Script 18 skipped because of a
# terra plotting-legend compatibility issue. This script does NOT redo any
# agreement, uncertainty, or area analysis. It only reads existing Script 18
# rasters and writes figures for:
#
#   SSP245 2081-2100
#   SSP585 2081-2100
#
# Each figure contains:
#   Panel A: number of GCMs (0-3) predicting suitability at MTSPS
#   Panel B: inter-GCM standard deviation of continuous suitability
#   Panel C: high-confidence suitable habitat (3/3 GCM agreement)
# =============================================================================

options(warn = 1)

cat("\n======================================================================\n")
cat("STARTING SCRIPT 18A — CREATE FOCAL UNCERTAINTY FIGURES\n")
cat("Species: Camellia japonica\n")
cat("======================================================================\n")

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

script18_root <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "18_interGCM_uncertainty_agreement_extrapolation"
)

agreement_dir <- file.path(script18_root, "02_GCM_agreement_maps")
high_confidence_dir <- file.path(script18_root, "03_high_confidence_maps")
figure_dir <- file.path(script18_root, "06_figures")
status_dir <- file.path(script18_root, "00_run_status")
log_dir <- file.path(script18_root, "08_logs")

started_marker <- file.path(status_dir, "18A_FIGURES_STARTED.txt")
completed_marker <- file.path(status_dir, "18A_FIGURES_COMPLETED.txt")
failed_marker <- file.path(status_dir, "18A_FIGURES_FAILED.txt")
log_file <- file.path(log_dir, "18A_focal_figure_log.txt")

species_name <- "Camellia_japonica"

focal_scenarios <- data.frame(
  ssp = c("SSP245", "SSP585"),
  period_start = c(2081L, 2081L),
  period_end = c(2100L, 2100L),
  stringsAsFactors = FALSE
)

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

write_log <- function(...) {
  line <- paste0(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ",
    paste0(..., collapse = "")
  )
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

assert_exists <- function(path, label) {
  if (!file.exists(path)) {
    stop("\nMissing ", label, ":\n", path, call. = FALSE)
  }
}

create_figure <- function(agreement, sd_raster, high_confidence,
                          scenario_label, output_pdf) {

  grDevices::pdf(output_pdf, width = 14, height = 5.5, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  par(
    mfrow = c(1, 3),
    mar = c(2.7, 2.8, 3.0, 1.0),
    oma = c(0, 0, 0, 0)
  )

  # Panel A: agreement count
  plot(
    agreement,
    col = c("#F0F0F0", "#FEE08B", "#FDAE61", "#1A9850"),
    breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5),
    plg = FALSE,
    main = "A. GCM agreement",
    axes = TRUE
  )
  legend(
    "bottomleft",
    legend = c("0/3", "1/3", "2/3", "3/3"),
    fill = c("#F0F0F0", "#FEE08B", "#FDAE61", "#1A9850"),
    title = "GCMs predicting\nsuitability",
    bty = "n",
    cex = 0.78,
    inset = 0.02
  )

  # Panel B: standard deviation
  plot(
    sd_raster,
    main = "B. Inter-GCM SD",
    axes = TRUE,
    plg = FALSE
  )
  sd_range <- terra::global(sd_raster, c("min", "max"), na.rm = TRUE)
  legend(
    "bottomleft",
    legend = paste0(
      "Range: ",
      format(round(sd_range[1, 1], 4), nsmall = 4),
      "–",
      format(round(sd_range[1, 2], 4), nsmall = 4)
    ),
    bty = "n",
    cex = 0.78,
    inset = 0.02
  )

  # Panel C: high-confidence suitable habitat
  plot(
    high_confidence,
    col = c("#F0F0F0", "#2166AC"),
    breaks = c(-0.5, 0.5, 1.5),
    plg = FALSE,
    main = "C. High-confidence suitability",
    axes = TRUE
  )
  legend(
    "bottomleft",
    legend = c("No", "Yes (3/3 GCMs)"),
    fill = c("#F0F0F0", "#2166AC"),
    title = "MTSPS suitability",
    bty = "n",
    cex = 0.78,
    inset = 0.02
  )

  mtext(
    scenario_label,
    side = 3,
    outer = FALSE,
    line = 0.6,
    font = 2,
    cex = 1.15
  )
}

tryCatch(
  {
    ensure_dir(figure_dir)
    ensure_dir(status_dir)
    ensure_dir(log_dir)

    unlink(completed_marker, force = TRUE)
    unlink(failed_marker, force = TRUE)

    writeLines(
      c(
        "Script 18A started",
        paste("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
      ),
      started_marker
    )

    if (!requireNamespace("terra", quietly = TRUE)) {
      install.packages("terra", repos = "https://cloud.r-project.org")
    }
    suppressPackageStartupMessages(library(terra))

    write_log("terra loaded.")

    output_records <- list()

    for (i in seq_len(nrow(focal_scenarios))) {
      ssp <- focal_scenarios$ssp[i]
      start <- focal_scenarios$period_start[i]
      end <- focal_scenarios$period_end[i]
      scenario_label <- paste0(ssp, " ", start, "-", end)

      agreement_path <- file.path(
        agreement_dir,
        paste0(
          species_name, "_GCM_suitable_agreement_",
          ssp, "_", start, "_", end, "_MTSPS.tif"
        )
      )

      sd_path <- file.path(
        project_dir,
        "global_reanalysis_outputs",
        "12_future_ensemble_3GCM",
        "02_ensemble_SD",
        paste0(
          species_name, "_future_ensemble_",
          ssp, "_", start, "_", end, "_SD_cloglog.tif"
        )
      )

      # Discover SD raster robustly if its exact filename differs.
      if (!file.exists(sd_path)) {
        sd_candidates <- list.files(
          file.path(
            project_dir,
            "global_reanalysis_outputs",
            "12_future_ensemble_3GCM"
          ),
          pattern = paste0(ssp, ".*", start, ".*", end, ".*(sd|SD).*\\.(tif|tiff)$"),
          full.names = TRUE,
          recursive = TRUE,
          ignore.case = TRUE
        )

        if (length(sd_candidates) == 1L) {
          sd_path <- sd_candidates[1]
        } else {
          stop(
            "\nCould not uniquely find the ensemble SD raster for ",
            scenario_label, ". Found ", length(sd_candidates), " candidates.",
            call. = FALSE
          )
        }
      }

      high_confidence_path <- file.path(
        high_confidence_dir,
        paste0(
          species_name, "_high_confidence_3of3_GCM_suitable_",
          ssp, "_", start, "_", end, "_MTSPS.tif"
        )
      )

      output_pdf <- file.path(
        figure_dir,
        paste0(
          "Figure_S_inter_GCM_agreement_SD_",
          ssp, "_", start, "_", end, ".pdf"
        )
      )

      assert_exists(agreement_path, paste0("agreement raster for ", scenario_label))
      assert_exists(sd_path, paste0("ensemble SD raster for ", scenario_label))
      assert_exists(high_confidence_path, paste0("high-confidence raster for ", scenario_label))

      agreement <- terra::rast(agreement_path)
      sd_raster <- terra::rast(sd_path)
      high_confidence <- terra::rast(high_confidence_path)

      same_grid <- terra::compareGeom(
        agreement, sd_raster,
        crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
        stopOnError = FALSE
      ) && terra::compareGeom(
        agreement, high_confidence,
        crs = TRUE, ext = TRUE, rowcol = TRUE, res = TRUE,
        stopOnError = FALSE
      )

      if (!same_grid) {
        stop(
          "\nGrid mismatch among focal figure rasters for ",
          scenario_label,
          call. = FALSE
        )
      }

      create_figure(
        agreement = agreement,
        sd_raster = sd_raster,
        high_confidence = high_confidence,
        scenario_label = scenario_label,
        output_pdf = output_pdf
      )

      write_log("Created figure: ", output_pdf)

      output_records[[length(output_records) + 1]] <- data.frame(
        scenario = scenario_label,
        agreement_tif = normalizePath(agreement_path, winslash = "/", mustWork = TRUE),
        SD_tif = normalizePath(sd_path, winslash = "/", mustWork = TRUE),
        high_confidence_tif = normalizePath(high_confidence_path, winslash = "/", mustWork = TRUE),
        figure_pdf = normalizePath(output_pdf, winslash = "/", mustWork = TRUE),
        stringsAsFactors = FALSE
      )
    }

    figure_manifest <- do.call(rbind, output_records)

    manifest_path <- file.path(
      figure_dir,
      "focal_uncertainty_figure_manifest.csv"
    )

    write.csv(figure_manifest, manifest_path, row.names = FALSE)

    writeLines(
      c(
        "SCRIPT 18A COMPLETED SUCCESSFULLY",
        paste("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
        paste("Figure manifest:", manifest_path)
      ),
      completed_marker
    )

    write_log("Script completed successfully.")

    cat("\n======================================================================\n")
    cat("SCRIPT 18A COMPLETED — FOCAL UNCERTAINTY FIGURES\n")
    cat("======================================================================\n")
    cat("Figure 1:\n", figure_manifest$figure_pdf[1], "\n", sep = "")
    cat("Figure 2:\n", figure_manifest$figure_pdf[2], "\n", sep = "")
  },
  error = function(error) {
    ensure_dir(status_dir)
    ensure_dir(log_dir)

    message_text <- paste0(
      "SCRIPT 18A FAILED\n",
      "Time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n",
      conditionMessage(error), "\n\n",
      "Log: ", log_file, "\n"
    )

    writeLines(message_text, failed_marker)

    cat("\n======================================================================\n")
    cat("SCRIPT 18A FAILED\n")
    cat("======================================================================\n")
    cat(message_text, "\n")

    stop(error)
  }
)
