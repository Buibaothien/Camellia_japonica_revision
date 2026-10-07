# ============================================================================
# 02_coordinatecleaner_global_without_sea.R
# ============================================================================
# Purpose
# -------
# Apply CoordinateCleaner to the global metadata-QC candidates from Step 01.
#
# Important methodological choice
# -------------------------------
# This script intentionally DOES NOT run CoordinateCleaner's "seas" test.
# In the previous workflow, automated sea checks generated implausibly large
# numbers of flags. The next stage uses ESA WorldCover and excludes water
# contexts directly, which is more appropriate for the present workflow.
#
# Tests applied:
#   capitals     - country capital coordinates
#   centroids    - country/province centroid coordinates
#   equal        - equal latitude/longitude
#   gbif         - GBIF headquarters coordinates
#   institutions - biodiversity institution coordinates
#   zeros        - zero coordinates and vicinity of 0/0
#
# Tests NOT applied:
#   seas, outliers, urban, ranges, countries, duplicates
#
# Input
# -----
# global_reanalysis_outputs/01_occurrence_QC/01_metadata_QC_candidates.csv
#
# Outputs
# -------
# global_reanalysis_outputs/02_coordinatecleaner/
# ============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ----------------------------------------------------------------------------
# 1. Packages
# ----------------------------------------------------------------------------

required_packages <- c("CoordinateCleaner", "readr", "dplyr")

missing_packages <- setdiff(required_packages, rownames(installed.packages()))

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(CoordinateCleaner)
  library(readr)
  library(dplyr)
})

# ----------------------------------------------------------------------------
# 2. User settings
# ----------------------------------------------------------------------------

project_dir <- "E:/PREDICTION TREE/Camellia_japonica"

input_csv <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "01_occurrence_QC",
  "01_metadata_QC_candidates.csv"
)

output_dir <- file.path(
  project_dir,
  "global_reanalysis_outputs",
  "02_coordinatecleaner"
)

# CoordinateCleaner tests to apply.
# Do NOT add "seas" here.
tests_to_apply <- c(
  "capitals",
  "centroids",
  "equal",
  "gbif",
  "institutions",
  "zeros"
)

# Defaults are explicitly stated for reproducibility.
capitals_radius_m <- 10000
centroids_radius_m <- 1000
institutions_radius_m <- 100
zeros_radius_degrees <- 0.5

# ----------------------------------------------------------------------------
# 3. Helper functions
# ----------------------------------------------------------------------------

stop_with_message <- function(message) {
  stop(paste0("\nERROR:\n", message), call. = FALSE)
}

find_summary_column <- function(column_names) {
  exact_candidates <- c(".summary", "summary", ".Summary")

  exact_match <- exact_candidates[
    exact_candidates %in% column_names
  ]

  if (length(exact_match) > 0) {
    return(exact_match[1])
  }

  possible_match <- grep(
    "summary",
    column_names,
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(possible_match) == 1) {
    return(possible_match)
  }

  stop_with_message(
    paste0(
      "Could not identify CoordinateCleaner summary column.\n",
      "Available columns:\n",
      paste(column_names, collapse = ", ")
    )
  )
}

safe_logical <- function(x) {
  if (is.logical(x)) {
    return(x)
  }

  tolower(as.character(x)) %in% c("true", "t", "1", "yes")
}

# ----------------------------------------------------------------------------
# 4. Validate input and read occurrence candidates
# ----------------------------------------------------------------------------

if (!file.exists(input_csv)) {
  stop_with_message(
    paste0(
      "Step 01 output was not found:\n",
      input_csv,
      "\n\nRun 01_global_occurrence_QC_v4_fix_csv_reader.py first."
    )
  )
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

occurrences <- readr::read_csv(
  input_csv,
  show_col_types = FALSE
)

required_columns <- c(
  "record_id",
  "decimalLongitude",
  "decimalLatitude"
)

missing_columns <- setdiff(required_columns, names(occurrences))

if (length(missing_columns) > 0) {
  stop_with_message(
    paste0(
      "Input file is missing required column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

occurrences <- occurrences %>%
  mutate(
    decimalLongitude = as.numeric(decimalLongitude),
    decimalLatitude = as.numeric(decimalLatitude),
    .input_row = row_number()
  ) %>%
  filter(
    is.finite(decimalLongitude),
    is.finite(decimalLatitude)
  )

if (nrow(occurrences) == 0) {
  stop_with_message("No valid coordinate records remain in the Step 01 output.")
}

# CoordinateCleaner does not require a species column when "outliers" and
# "duplicates" are not requested. A placeholder avoids version-specific issues.
occurrences_for_cleaner <- occurrences %>%
  mutate(
    species_for_cleaner = "Camellia japonica"
  )

# ----------------------------------------------------------------------------
# 5. Run CoordinateCleaner
# ----------------------------------------------------------------------------

message("\n======================================================================")
message("GLOBAL COORDINATECLEANER: STEP 02")
message("Input records: ", nrow(occurrences_for_cleaner))
message("Tests: ", paste(tests_to_apply, collapse = ", "))
message("Sea test: NOT USED")
message("======================================================================\n")

cleaned_with_flags <- CoordinateCleaner::clean_coordinates(
  x = as.data.frame(occurrences_for_cleaner),
  lon = "decimalLongitude",
  lat = "decimalLatitude",
  species = "species_for_cleaner",
  tests = tests_to_apply,
  capitals_rad = capitals_radius_m,
  centroids_rad = centroids_radius_m,
  centroids_detail = "both",
  inst_rad = institutions_radius_m,
  zeros_rad = zeros_radius_degrees,
  value = "spatialvalid",
  verbose = TRUE
)

cleaned_with_flags <- as.data.frame(cleaned_with_flags)

summary_column <- find_summary_column(names(cleaned_with_flags))

cleaned_with_flags$coordinatecleaner_keep <- safe_logical(
  cleaned_with_flags[[summary_column]]
)

# ----------------------------------------------------------------------------
# 6. Create audit table and retained candidate table
# ----------------------------------------------------------------------------

cleaned_with_flags <- cleaned_with_flags %>%
  arrange(.input_row)

retained_candidates <- cleaned_with_flags %>%
  filter(coordinatecleaner_keep) %>%
  select(-species_for_cleaner, -all_of(summary_column)) %>%
  arrange(record_id)

# Find only the CoordinateCleaner test-flag columns for summary reporting.
possible_flag_columns <- names(cleaned_with_flags)[
  grepl("^\\.", names(cleaned_with_flags))
]

possible_flag_columns <- setdiff(
  possible_flag_columns,
  c(summary_column)
)

flag_summary <- data.frame(
  test_column = c(
    possible_flag_columns,
    "coordinatecleaner_keep"
  ),
  n_pass = c(
    vapply(
      possible_flag_columns,
      function(column_name) {
        sum(safe_logical(cleaned_with_flags[[column_name]]), na.rm = TRUE)
      },
      numeric(1)
    ),
    sum(cleaned_with_flags$coordinatecleaner_keep, na.rm = TRUE)
  ),
  n_flagged_or_removed = c(
    vapply(
      possible_flag_columns,
      function(column_name) {
        sum(!safe_logical(cleaned_with_flags[[column_name]]), na.rm = TRUE)
      },
      numeric(1)
    ),
    sum(!cleaned_with_flags$coordinatecleaner_keep, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)

overall_summary <- data.frame(
  stage = c(
    "Step_01_metadata_QC_candidates",
    "CoordinateCleaner_retained",
    "CoordinateCleaner_removed"
  ),
  n_records = c(
    nrow(cleaned_with_flags),
    sum(cleaned_with_flags$coordinatecleaner_keep, na.rm = TRUE),
    sum(!cleaned_with_flags$coordinatecleaner_keep, na.rm = TRUE)
  )
)

settings_table <- data.frame(
  parameter = c(
    "input_csv",
    "tests_applied",
    "sea_test",
    "outlier_test",
    "urban_test",
    "capital_radius_m",
    "centroid_radius_m",
    "institution_radius_m",
    "zero_radius_degrees",
    "CoordinateCleaner_version"
  ),
  value = c(
    input_csv,
    paste(tests_to_apply, collapse = "; "),
    "Not applied",
    "Not applied",
    "Not applied",
    capitals_radius_m,
    centroids_radius_m,
    institutions_radius_m,
    zeros_radius_degrees,
    as.character(utils::packageVersion("CoordinateCleaner"))
  ),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------------------------------
# 7. Save outputs
# ----------------------------------------------------------------------------

readr::write_csv(
  cleaned_with_flags,
  file.path(
    output_dir,
    "02_coordinatecleaner_all_records_with_flags.csv"
  ),
  na = ""
)

readr::write_csv(
  retained_candidates,
  file.path(
    output_dir,
    "02_coordinatecleaner_retained_candidates.csv"
  ),
  na = ""
)

readr::write_csv(
  flag_summary,
  file.path(
    output_dir,
    "02_coordinatecleaner_test_flag_summary.csv"
  ),
  na = ""
)

readr::write_csv(
  overall_summary,
  file.path(
    output_dir,
    "02_coordinatecleaner_overall_summary.csv"
  ),
  na = ""
)

readr::write_csv(
  settings_table,
  file.path(
    output_dir,
    "02_coordinatecleaner_settings.csv"
  ),
  na = ""
)

writeLines(
  c(
    "CoordinateCleaner was applied after metadata-based global QC.",
    "The sea test was intentionally not applied because prior automated sea",
    "screening produced implausibly extensive flags in this dataset.",
    "ESA WorldCover screening will subsequently exclude water-context points.",
    "The retained table is an intermediate candidate dataset only; LULC",
    "screening and 5 km spatial thinning remain to be completed."
  ),
  con = file.path(
    output_dir,
    "02_coordinatecleaner_interpretation_note.txt"
  )
)

message("\n======================================================================")
message("CoordinateCleaner completed successfully")
message("Input records: ", nrow(cleaned_with_flags))
message(
  "Retained records: ",
  sum(cleaned_with_flags$coordinatecleaner_keep, na.rm = TRUE)
)
message(
  "Flagged/removed records: ",
  sum(!cleaned_with_flags$coordinatecleaner_keep, na.rm = TRUE)
)
message("Output folder:")
message(output_dir)
message("======================================================================\n")
