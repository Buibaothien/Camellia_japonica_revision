"""
07_confirm_final_predictors_VIF_and_prepare_current_inputs.py

Purpose
-------
Confirm a parsimonious, low-collinearity primary predictor set for the global
Camellia japonica MaxEnt / maxnet model, calculate VIF from the actual
extracted values, and create a clean current-predictor folder for modelling.

Selected primary predictors
---------------------------
Bio2  : Mean diurnal temperature range
Bio6  : Minimum temperature of the coldest month
Bio14 : Precipitation of the driest month
Bio18 : Precipitation of the warmest quarter
Elv   : Elevation

Selection rationale
-------------------
- Bio6 represents winter cold limitation while avoiding redundancy with
  Bio1, Bio9 and Bio11.
- Bio2 adds thermal variability and is not strongly correlated with Bio6.
- Bio14 represents dry-season water limitation while avoiding the highly
  correlated Bio17 and Bio19.
- Bio18 represents warm-season moisture while avoiding the highly correlated
  Bio13 and Bio16.
- Elevation is retained as a topographic control and its effect will be
  interpreted cautiously because it is static in future projections.

The selection was made after reviewing both:
- correlations at final occurrence points (n = 828);
- correlations across the 20,000-cell global environmental screening sample.

This script verifies the decision using actual VIF values. It does not use
any old M-based occurrence dataset or prior predictor set.

Inputs
------
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
06_predictor_screening\06_occurrence_primary_predictor_values.csv

E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
06_predictor_screening\06_global_environmental_screening_values_complete.csv

E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
06_predictor_screening\06_primary_candidate_raster_inventory.csv

Outputs
-------
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
07_final_predictors_VIF\

Current model inputs created at:
E:\PREDICTION TREE\Camellia_japonica\model_inputs\
Primary_Global_ClimateTopo_Current\
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


# ============================================================================
# 1. USER SETTINGS
# ============================================================================

PROJECT_DIR = Path(r"E:\PREDICTION TREE\Camellia_japonica")

SCREENING_DIR = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "06_predictor_screening"
)

OCCURRENCE_VALUES_CSV = (
    SCREENING_DIR
    / "06_occurrence_primary_predictor_values.csv"
)

GLOBAL_VALUES_CSV = (
    SCREENING_DIR
    / "06_global_environmental_screening_values_complete.csv"
)

RASTER_INVENTORY_CSV = (
    SCREENING_DIR
    / "06_primary_candidate_raster_inventory.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "07_final_predictors_VIF"
)

CURRENT_MODEL_INPUT_DIR = (
    PROJECT_DIR
    / "model_inputs"
    / "Primary_Global_ClimateTopo_Current"
)

# Final primary predictor set for the global model.
SELECTED_PREDICTORS = ["Bio2", "Bio6", "Bio14", "Bio18", "Elv"]

# Interpretation threshold only. All VIF values are exported irrespective of
# this threshold; a script error is raised only if any VIF is >= 5.
VIF_THRESHOLD = 5.0


# ============================================================================
# 2. PACKAGE CHECKS
# ============================================================================

def ensure_package(import_name: str, pip_name: str | None = None) -> None:
    try:
        __import__(import_name)
    except ImportError:
        package_name = pip_name or import_name
        print(f"Installing missing package: {package_name}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", package_name]
        )


ensure_package("sklearn", "scikit-learn")
ensure_package("rasterio")

from sklearn.linear_model import LinearRegression
import rasterio


# ============================================================================
# 3. HELPER FUNCTIONS
# ============================================================================

def stop_with_message(message: str) -> None:
    raise RuntimeError(f"\nERROR:\n{message}")


def validate_input_file(path: Path, description: str) -> None:
    if not path.exists():
        stop_with_message(
            f"{description} was not found:\n{path}\n\n"
            "Run Step 06 before running this script."
        )


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [column for column in columns if column not in df.columns]

    if missing:
        stop_with_message(
            f"{label} lacks required predictor column(s): "
            + ", ".join(missing)
        )


def calculate_vif(df: pd.DataFrame, predictors: list[str], dataset_label: str):
    """
    Calculate VIF using R² from a linear regression of each variable against
    all remaining predictors. Complete rows only are used.
    """
    numeric = df[predictors].apply(pd.to_numeric, errors="coerce").dropna()

    if len(numeric) < len(predictors) + 2:
        stop_with_message(
            f"Too few complete rows for VIF calculation in {dataset_label}."
        )

    rows = []

    for predictor in predictors:
        others = [name for name in predictors if name != predictor]

        model = LinearRegression()
        model.fit(numeric[others], numeric[predictor])

        r_squared = float(model.score(numeric[others], numeric[predictor]))

        # Prevent division by zero caused by a numerical edge case.
        if r_squared >= 0.999999999:
            vif = np.inf
        else:
            vif = 1.0 / (1.0 - r_squared)

        rows.append(
            {
                "dataset": dataset_label,
                "predictor": predictor,
                "n_complete_records": len(numeric),
                "R_squared_against_remaining_predictors": r_squared,
                "VIF": vif,
                "VIF_below_5": bool(vif < VIF_THRESHOLD),
            }
        )

    return pd.DataFrame(rows)


def correlation_pairs(df: pd.DataFrame, predictors: list[str], dataset_label: str):
    numeric = df[predictors].apply(pd.to_numeric, errors="coerce")

    rows = []

    for i, predictor_1 in enumerate(predictors):
        for predictor_2 in predictors[i + 1:]:
            pair = numeric[[predictor_1, predictor_2]].dropna()

            r_value = pair[predictor_1].corr(pair[predictor_2])

            rows.append(
                {
                    "dataset": dataset_label,
                    "predictor_1": predictor_1,
                    "predictor_2": predictor_2,
                    "pearson_r": r_value,
                    "abs_pearson_r": abs(r_value),
                    "n_complete_pairs": len(pair),
                    "abs_r_below_080": bool(abs(r_value) < 0.80),
                }
            )

    return pd.DataFrame(rows)


def copy_selected_rasters(inventory: pd.DataFrame) -> pd.DataFrame:
    """
    Copy and rename selected rasters to canonical predictor-code filenames.
    """
    CURRENT_MODEL_INPUT_DIR.mkdir(parents=True, exist_ok=True)

    selected_inventory = inventory.loc[
        inventory["predictor"].isin(SELECTED_PREDICTORS)
    ].copy()

    if selected_inventory["predictor"].nunique() != len(SELECTED_PREDICTORS):
        found = sorted(selected_inventory["predictor"].astype(str).unique())
        missing = sorted(set(SELECTED_PREDICTORS) - set(found))

        stop_with_message(
            "Not all selected predictor rasters were found in the Step 06 inventory.\n"
            f"Missing: {missing}"
        )

    copy_rows = []

    for predictor in SELECTED_PREDICTORS:
        row = selected_inventory.loc[
            selected_inventory["predictor"] == predictor
        ]

        source = Path(row.iloc[0]["raster_file"])

        if not source.exists():
            stop_with_message(
                f"Raster file for {predictor} does not exist:\n{source}"
            )

        destination = CURRENT_MODEL_INPUT_DIR / f"{predictor}.tif"

        # Use copy2 to preserve timestamps and basic filesystem metadata.
        shutil.copy2(source, destination)

        with rasterio.open(destination) as src:
            copy_rows.append(
                {
                    "predictor": predictor,
                    "source_raster": str(source),
                    "destination_raster": str(destination),
                    "width": src.width,
                    "height": src.height,
                    "CRS": str(src.crs),
                    "resolution_x": src.res[0],
                    "resolution_y": src.res[1],
                    "nodata": src.nodata,
                }
            )

    return pd.DataFrame(copy_rows)


# ============================================================================
# 4. MAIN
# ============================================================================

def main() -> None:
    for path, description in [
        (OCCURRENCE_VALUES_CSV, "Occurrence predictor-value table"),
        (GLOBAL_VALUES_CSV, "Global screening predictor-value table"),
        (RASTER_INVENTORY_CSV, "Primary raster inventory"),
    ]:
        validate_input_file(path, description)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 82)
    print("FINAL PREDICTOR CONFIRMATION + VIF + CURRENT INPUT PREPARATION: STEP 07")
    print("Selected predictors: " + ", ".join(SELECTED_PREDICTORS))
    print("VIF acceptance threshold: < " + str(VIF_THRESHOLD))
    print("=" * 82)

    occurrence_values = pd.read_csv(OCCURRENCE_VALUES_CSV, low_memory=False)
    global_values = pd.read_csv(GLOBAL_VALUES_CSV, low_memory=False)
    raster_inventory = pd.read_csv(RASTER_INVENTORY_CSV, low_memory=False)

    require_columns(
        occurrence_values,
        SELECTED_PREDICTORS,
        "Occurrence predictor-value table"
    )

    require_columns(
        global_values,
        SELECTED_PREDICTORS,
        "Global predictor-value table"
    )

    require_columns(
        raster_inventory,
        ["predictor", "raster_file"],
        "Primary raster inventory"
    )

    # ------------------------------------------------------------------------
    # 4.1 VIF and pairwise-correlation checks
    # ------------------------------------------------------------------------

    occurrence_vif = calculate_vif(
        occurrence_values,
        SELECTED_PREDICTORS,
        "final_occurrence_records"
    )

    global_vif = calculate_vif(
        global_values,
        SELECTED_PREDICTORS,
        "global_environmental_screening_sample"
    )

    all_vif = pd.concat(
        [occurrence_vif, global_vif],
        ignore_index=True
    )

    occurrence_pairs = correlation_pairs(
        occurrence_values,
        SELECTED_PREDICTORS,
        "final_occurrence_records"
    )

    global_pairs = correlation_pairs(
        global_values,
        SELECTED_PREDICTORS,
        "global_environmental_screening_sample"
    )

    all_pairs = pd.concat(
        [occurrence_pairs, global_pairs],
        ignore_index=True
    )

    failed_vif = all_vif.loc[~all_vif["VIF_below_5"]]

    failed_pairs = all_pairs.loc[~all_pairs["abs_r_below_080"]]

    # ------------------------------------------------------------------------
    # 4.2 Copy canonical current model input rasters
    # ------------------------------------------------------------------------

    raster_copy_inventory = copy_selected_rasters(raster_inventory)

    # ------------------------------------------------------------------------
    # 4.3 Export transparent rationale and results
    # ------------------------------------------------------------------------

    selection_rationale = pd.DataFrame(
        {
            "predictor": SELECTED_PREDICTORS,
            "variable_name": [
                "Mean diurnal temperature range",
                "Minimum temperature of the coldest month",
                "Precipitation of the driest month",
                "Precipitation of the warmest quarter",
                "Elevation",
            ],
            "functional_dimension": [
                "Diurnal thermal variability",
                "Winter cold limitation",
                "Dry-season water limitation",
                "Warm-season moisture availability",
                "Topographic climatic context",
            ],
            "selection_rationale": [
                (
                    "Retained to represent daily thermal variability while "
                    "avoiding redundancy with Bio3 and Bio4."
                ),
                (
                    "Retained as the directly interpretable cold-extreme "
                    "variable; Bio1, Bio9 and Bio11 were excluded as "
                    "redundant thermal alternatives."
                ),
                (
                    "Retained to represent dry-season water limitation; "
                    "Bio17 and Bio19 were excluded because of strong "
                    "correlation with Bio14."
                ),
                (
                    "Retained to represent precipitation during the warmest "
                    "quarter; Bio13 and Bio16 were excluded because of "
                    "strong correlation with Bio18."
                ),
                (
                    "Retained as a static topographic control. Future "
                    "interpretation must not treat elevation as a changing "
                    "climate predictor."
                ),
            ],
            "future_projection_status": [
                "Matched CMIP6 bioclimatic projection available",
                "Matched CMIP6 bioclimatic projection available",
                "Matched CMIP6 bioclimatic projection available",
                "Matched CMIP6 bioclimatic projection available",
                "Static current elevation retained consistently",
            ],
        }
    )

    selection_rationale.to_csv(
        OUTPUT_DIR / "07_final_primary_predictor_selection_rationale.csv",
        index=False,
        encoding="utf-8-sig"
    )

    all_vif.to_csv(
        OUTPUT_DIR / "07_VIF_results_occurrences_and_global_sample.csv",
        index=False,
        encoding="utf-8-sig"
    )

    all_pairs.to_csv(
        OUTPUT_DIR / "07_selected_predictor_pairwise_correlations.csv",
        index=False,
        encoding="utf-8-sig"
    )

    raster_copy_inventory.to_csv(
        OUTPUT_DIR / "07_current_model_input_raster_inventory.csv",
        index=False,
        encoding="utf-8-sig"
    )

    summary = pd.DataFrame(
        {
            "metric": [
                "selected_predictor_count",
                "maximum_VIF_occurrence_records",
                "maximum_VIF_global_screening_sample",
                "maximum_abs_Pearson_r_occurrence_records",
                "maximum_abs_Pearson_r_global_screening_sample",
                "all_VIF_below_5",
                "all_pairwise_abs_r_below_080",
                "current_model_input_folder",
            ],
            "value": [
                len(SELECTED_PREDICTORS),
                occurrence_vif["VIF"].max(),
                global_vif["VIF"].max(),
                occurrence_pairs["abs_pearson_r"].max(),
                global_pairs["abs_pearson_r"].max(),
                bool(all_vif["VIF_below_5"].all()),
                bool(all_pairs["abs_r_below_080"].all()),
                str(CURRENT_MODEL_INPUT_DIR),
            ],
        }
    )

    summary.to_csv(
        OUTPUT_DIR / "07_predictor_selection_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    pd.DataFrame(
        {
            "parameter": [
                "occurrence_dataset",
                "global_screening_dataset",
                "selected_predictors",
                "VIF_threshold",
                "pairwise_correlation_threshold",
                "old_M_restricted_dataset_used",
                "soil_moisture_used",
                "UV_used",
                "NPP_or_runoff_used",
            ],
            "value": [
                "Strict global natural-land-cover candidates after 5 km thinning",
                "20,000 global terrestrial environmental screening cells",
                "; ".join(SELECTED_PREDICTORS),
                VIF_THRESHOLD,
                0.80,
                "No",
                "No",
                "No",
                "No",
            ],
        }
    ).to_csv(
        OUTPUT_DIR / "07_predictor_selection_settings.csv",
        index=False,
        encoding="utf-8-sig"
    )

    note = [
        "The primary global model uses Bio2, Bio6, Bio14, Bio18, and elevation.",
        "The selected set was chosen after correlation screening of occurrence and global environmental samples.",
        "VIF is calculated separately for the occurrence and global screening datasets.",
        "The model does not use static soil moisture, NPP, annual runoff, soil chemistry, or UV layers because matched future scenario layers were not included.",
        "Current rasters were copied to a clean model-input folder with canonical names.",
    ]

    (OUTPUT_DIR / "07_predictor_selection_interpretation_note.txt").write_text(
        "\n".join(note),
        encoding="utf-8"
    )

    # ------------------------------------------------------------------------
    # 4.4 Console output and safety check
    # ------------------------------------------------------------------------

    print("\n" + "=" * 82)
    print("Predictor confirmation completed")
    print("Selected predictors: " + ", ".join(SELECTED_PREDICTORS))
    print(
        "Maximum occurrence-record VIF: "
        f"{occurrence_vif['VIF'].max():.4f}"
    )
    print(
        "Maximum global-screening-sample VIF: "
        f"{global_vif['VIF'].max():.4f}"
    )
    print(
        "Maximum occurrence-record |r|: "
        f"{occurrence_pairs['abs_pearson_r'].max():.4f}"
    )
    print(
        "Maximum global-screening-sample |r|: "
        f"{global_pairs['abs_pearson_r'].max():.4f}"
    )
    print("Current input rasters copied to:")
    print(CURRENT_MODEL_INPUT_DIR)
    print("Output folder:")
    print(OUTPUT_DIR)
    print("=" * 82 + "\n")

    if not failed_vif.empty:
        stop_with_message(
            "At least one selected predictor has VIF >= 5. "
            "Do not proceed to tuning until the predictor set is revised.\n\n"
            + failed_vif.to_string(index=False)
        )

    if not failed_pairs.empty:
        stop_with_message(
            "At least one selected predictor pair has |r| >= 0.80. "
            "Do not proceed to tuning until the predictor set is revised.\n\n"
            + failed_pairs.to_string(index=False)
        )


if __name__ == "__main__":
    main()
