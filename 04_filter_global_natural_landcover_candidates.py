"""
04_filter_global_natural_landcover_candidates.py

Purpose
-------
Create the primary global occurrence dataset after ESA WorldCover screening.

This script keeps the strict natural-land-cover context candidates as the
PRIMARY occurrence set for the global SDM workflow. It also exports a broader
LULC candidate set for supplementary sensitivity analysis.

Important interpretation
------------------------
The retained records must be described in the manuscript as:

    "global natural-land-cover occurrence candidates"

They must NOT be described as confirmed wild populations, because LULC and
metadata screening cannot prove that each individual occurrence is wild.

Input
-----
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
03_worldcover\03_ESA_WorldCover_context.csv

Outputs
-------
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
04_LULC_filter\

Next step
---------
Run 05_global_spatial_thinning_5km.py using:
04_global_strict_natural_landcover_candidates.csv
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd


# ============================================================================
# 1. USER SETTINGS
# ============================================================================

PROJECT_DIR = Path(r"E:\PREDICTION TREE\Camellia_japonica")

INPUT_CSV = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "03_worldcover"
    / "03_ESA_WorldCover_context.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "04_LULC_filter"
)

# Primary choice for the main SDM workflow.
PRIMARY_RULE = "strict"

# Strict and broad thresholds are written explicitly again here, even though
# they were already calculated in Step 03. This creates an auditable record
# and independently verifies the Boolean columns from Earth Engine.
STRICT_NATURAL_MIN_TREE_SHRUB = 0.50
STRICT_NATURAL_MAX_CROPLAND = 0.20
STRICT_NATURAL_MAX_BUILT_UP = 0.05
STRICT_NATURAL_MAX_WATER = 0.10

BROAD_NATURAL_MIN_TREE_SHRUB = 0.50
BROAD_NATURAL_MAX_CROPLAND = 0.30
BROAD_NATURAL_MAX_BUILT_UP = 0.10
BROAD_NATURAL_MAX_WATER = 0.15

TREE_SHRUB_CODES = {10, 20}


# ============================================================================
# 2. HELPER FUNCTIONS
# ============================================================================

def stop_with_message(message: str) -> None:
    raise RuntimeError(f"\nERROR:\n{message}")


def as_boolean(series: pd.Series) -> pd.Series:
    """Safely convert mixed CSV Boolean representations to True/False."""
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)

    normalised = (
        series.fillna("")
        .astype(str)
        .str.strip()
        .str.lower()
    )

    return normalised.isin({"true", "t", "1", "yes", "y"})


def require_columns(df: pd.DataFrame, columns: list[str]) -> None:
    missing = [column for column in columns if column not in df.columns]

    if missing:
        stop_with_message(
            "Input file lacks required column(s):\n"
            + ", ".join(missing)
        )


# ============================================================================
# 3. MAIN
# ============================================================================

def main() -> None:
    if not INPUT_CSV.exists():
        stop_with_message(
            "Step 03 output was not found:\n"
            f"{INPUT_CSV}\n\n"
            "Run 03_sample_ESA_WorldCover_direct_VSCode_project_71a05c44.py first."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 78)
    print("GLOBAL LULC FILTER: STEP 04")
    print(f"Input: {INPUT_CSV}")
    print("=" * 78)

    df = pd.read_csv(INPUT_CSV, low_memory=False)

    required_columns = [
        "record_id",
        "decimalLongitude",
        "decimalLatitude",
        "worldcover_center_code",
        "frac_tree_cover",
        "frac_shrubland",
        "frac_cropland",
        "frac_built_up",
        "frac_permanent_water",
    ]

    require_columns(df, required_columns)

    # Standardise numeric variables.
    numeric_columns = [
        "decimalLongitude",
        "decimalLatitude",
        "worldcover_center_code",
        "frac_tree_cover",
        "frac_shrubland",
        "frac_cropland",
        "frac_built_up",
        "frac_permanent_water",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    # Ensure exactly one row per retained record.
    df = df.sort_values("record_id", kind="mergesort").drop_duplicates(
        subset="record_id",
        keep="first"
    ).reset_index(drop=True)

    # Recalculate natural-context screening rules locally.
    # This ensures that future changes in Earth Engine Boolean storage do not
    # affect reproducibility.
    df["frac_tree_shrub_natural_recalculated"] = (
        df["frac_tree_cover"].fillna(0)
        + df["frac_shrubland"].fillna(0)
    )

    df["center_is_tree_or_shrub_recalculated"] = (
        df["worldcover_center_code"]
        .round()
        .astype("Int64")
        .isin(TREE_SHRUB_CODES)
    )

    df["strict_natural_context_recalculated"] = (
        df["center_is_tree_or_shrub_recalculated"]
        & (
            df["frac_tree_shrub_natural_recalculated"]
            >= STRICT_NATURAL_MIN_TREE_SHRUB
        )
        & (df["frac_cropland"].fillna(1) <= STRICT_NATURAL_MAX_CROPLAND)
        & (df["frac_built_up"].fillna(1) <= STRICT_NATURAL_MAX_BUILT_UP)
        & (
            df["frac_permanent_water"].fillna(1)
            <= STRICT_NATURAL_MAX_WATER
        )
    )

    df["broad_natural_context_recalculated"] = (
        (
            df["frac_tree_shrub_natural_recalculated"]
            >= BROAD_NATURAL_MIN_TREE_SHRUB
        )
        & (df["frac_cropland"].fillna(1) <= BROAD_NATURAL_MAX_CROPLAND)
        & (df["frac_built_up"].fillna(1) <= BROAD_NATURAL_MAX_BUILT_UP)
        & (
            df["frac_permanent_water"].fillna(1)
            <= BROAD_NATURAL_MAX_WATER
        )
    )

    # Compare Earth Engine source flag with locally recalculated condition.
    if "strict_natural_context_candidate" in df.columns:
        df["strict_ee_flag"] = as_boolean(
            df["strict_natural_context_candidate"]
        )
    else:
        df["strict_ee_flag"] = np.nan

    if "broad_natural_context_candidate" in df.columns:
        df["broad_ee_flag"] = as_boolean(
            df["broad_natural_context_candidate"]
        )
    else:
        df["broad_ee_flag"] = np.nan

    # Use the transparent, local recalculation as final classification.
    df["strict_natural_context_final"] = (
        df["strict_natural_context_recalculated"]
    )

    df["broad_natural_context_final"] = (
        df["broad_natural_context_recalculated"]
    )

    # Classify all points for audit and later sensitivity analysis.
    df["LULC_filter_status"] = np.select(
        [
            df["strict_natural_context_final"],
            (
                ~df["strict_natural_context_final"]
                & df["broad_natural_context_final"]
            ),
        ],
        [
            "strict_primary_candidate",
            "broad_only_sensitivity_candidate",
        ],
        default="excluded_non_natural_or_uncertain_context",
    )

    df["primary_model_candidate"] = (
        df["LULC_filter_status"] == "strict_primary_candidate"
    )

    df["sensitivity_model_candidate_broad"] = (
        df["LULC_filter_status"].isin(
            [
                "strict_primary_candidate",
                "broad_only_sensitivity_candidate",
            ]
        )
    )

    # Preserve a clear reason for records that fail strict screening.
    df["strict_LULC_exclusion_reason"] = ""

    is_non_tree_shrub_center = ~df["center_is_tree_or_shrub_recalculated"]
    too_little_natural = (
        df["frac_tree_shrub_natural_recalculated"]
        < STRICT_NATURAL_MIN_TREE_SHRUB
    )
    too_much_crop = (
        df["frac_cropland"].fillna(1)
        > STRICT_NATURAL_MAX_CROPLAND
    )
    too_much_built = (
        df["frac_built_up"].fillna(1)
        > STRICT_NATURAL_MAX_BUILT_UP
    )
    too_much_water = (
        df["frac_permanent_water"].fillna(1)
        > STRICT_NATURAL_MAX_WATER
    )

    failed_strict = ~df["strict_natural_context_final"]

    df.loc[
        failed_strict & is_non_tree_shrub_center,
        "strict_LULC_exclusion_reason"
    ] += "center_not_tree_or_shrub;"

    df.loc[
        failed_strict & too_little_natural,
        "strict_LULC_exclusion_reason"
    ] += "tree_shrub_fraction_below_0.50;"

    df.loc[
        failed_strict & too_much_crop,
        "strict_LULC_exclusion_reason"
    ] += "cropland_fraction_above_0.20;"

    df.loc[
        failed_strict & too_much_built,
        "strict_LULC_exclusion_reason"
    ] += "built_up_fraction_above_0.05;"

    df.loc[
        failed_strict & too_much_water,
        "strict_LULC_exclusion_reason"
    ] += "permanent_water_fraction_above_0.10;"

    df["strict_LULC_exclusion_reason"] = (
        df["strict_LULC_exclusion_reason"]
        .str.rstrip(";")
        .replace("", np.nan)
    )

    # ------------------------------------------------------------------------
    # 4. Write primary and supplementary datasets
    # ------------------------------------------------------------------------

    primary_candidates = (
        df.loc[df["primary_model_candidate"]]
        .copy()
        .sort_values("record_id")
        .reset_index(drop=True)
    )

    broad_candidates = (
        df.loc[df["sensitivity_model_candidate_broad"]]
        .copy()
        .sort_values("record_id")
        .reset_index(drop=True)
    )

    all_audit = df.sort_values("record_id").reset_index(drop=True)

    primary_candidates.to_csv(
        OUTPUT_DIR / "04_global_strict_natural_landcover_candidates.csv",
        index=False,
        encoding="utf-8-sig"
    )

    broad_candidates.to_csv(
        OUTPUT_DIR / "04_global_broad_natural_landcover_candidates.csv",
        index=False,
        encoding="utf-8-sig"
    )

    all_audit.to_csv(
        OUTPUT_DIR / "04_global_LULC_all_candidates_audit.csv",
        index=False,
        encoding="utf-8-sig"
    )

    # ------------------------------------------------------------------------
    # 5. Summary tables
    # ------------------------------------------------------------------------

    summary = pd.DataFrame(
        {
            "stage": [
                "CoordinateCleaner_retained_input",
                "Strict_natural_landcover_primary_candidates",
                "Broad_natural_landcover_sensitivity_candidates",
                "Broad_only_candidates",
                "Excluded_non_natural_or_uncertain_context",
            ],
            "n_records": [
                len(df),
                len(primary_candidates),
                len(broad_candidates),
                int(
                    (
                        df["LULC_filter_status"]
                        == "broad_only_sensitivity_candidate"
                    ).sum()
                ),
                int(
                    (
                        df["LULC_filter_status"]
                        == "excluded_non_natural_or_uncertain_context"
                    ).sum()
                ),
            ],
        }
    )

    summary.to_csv(
        OUTPUT_DIR / "04_global_LULC_filter_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    center_class_distribution = (
        all_audit.assign(
            worldcover_center_label=all_audit.get(
                "worldcover_center_label",
                "Unknown_or_missing"
            )
        )
        .groupby(
            ["worldcover_center_code", "worldcover_center_label"],
            dropna=False
        )
        .agg(
            n_records=("record_id", "size"),
            strict_primary_candidates=(
                "primary_model_candidate",
                "sum"
            ),
            broad_sensitivity_candidates=(
                "sensitivity_model_candidate_broad",
                "sum"
            ),
        )
        .reset_index()
        .sort_values(
            ["worldcover_center_code", "worldcover_center_label"],
            na_position="last"
        )
    )

    center_class_distribution.to_csv(
        OUTPUT_DIR / "04_worldcover_center_class_distribution.csv",
        index=False,
        encoding="utf-8-sig"
    )

    agreement = pd.DataFrame(
        {
            "metric": [
                "strict_EE_flag_available",
                "strict_EE_flag_matches_local_recalculation",
                "broad_EE_flag_available",
                "broad_EE_flag_matches_local_recalculation",
            ],
            "n_records": [
                int(pd.notna(df["strict_ee_flag"]).sum()),
                int(
                    (
                        df["strict_ee_flag"]
                        == df["strict_natural_context_recalculated"]
                    ).sum()
                )
                if pd.notna(df["strict_ee_flag"]).any()
                else np.nan,
                int(pd.notna(df["broad_ee_flag"]).sum()),
                int(
                    (
                        df["broad_ee_flag"]
                        == df["broad_natural_context_recalculated"]
                    ).sum()
                )
                if pd.notna(df["broad_ee_flag"]).any()
                else np.nan,
            ],
        }
    )

    agreement.to_csv(
        OUTPUT_DIR / "04_EE_local_rule_agreement.csv",
        index=False,
        encoding="utf-8-sig"
    )

    settings = pd.DataFrame(
        {
            "parameter": [
                "primary_dataset",
                "primary_LULC_rule",
                "broad_sensitivity_rule",
                "interpretation",
                "next_step",
            ],
            "value": [
                "Strict global natural-land-cover occurrence candidates",
                (
                    "Central pixel is Tree cover or Shrubland; "
                    "Tree+Shrub fraction >= 0.50; Cropland <= 0.20; "
                    "Built-up <= 0.05; Permanent water <= 0.10"
                ),
                (
                    "Tree+Shrub fraction >= 0.50; Cropland <= 0.30; "
                    "Built-up <= 0.10; Permanent water <= 0.15"
                ),
                (
                    "Natural-land-cover context candidate only; "
                    "not confirmed wild occurrence"
                ),
                (
                    "Spatial thinning at 5 km and removal of repeated "
                    "2.5 arc-minute raster cells"
                ),
            ],
        }
    )

    settings.to_csv(
        OUTPUT_DIR / "04_LULC_filter_settings.csv",
        index=False,
        encoding="utf-8-sig"
    )

    note = [
        "The primary occurrence dataset contains strict global natural-land-cover context candidates.",
        "The LULC procedure does not prove that all retained occurrences are wild populations.",
        "The broad candidate set is retained for a supplementary sensitivity analysis.",
        "The primary dataset must next undergo 5 km spatial thinning and raster-cell deduplication.",
    ]

    (OUTPUT_DIR / "04_LULC_filter_interpretation_note.txt").write_text(
        "\n".join(note),
        encoding="utf-8"
    )

    print("\n" + "=" * 78)
    print("Global LULC filter completed successfully")
    print(f"CoordinateCleaner input candidates: {len(df)}")
    print(
        "Strict primary natural-land-cover candidates: "
        f"{len(primary_candidates)}"
    )
    print(
        "Broad sensitivity natural-land-cover candidates: "
        f"{len(broad_candidates)}"
    )
    print(
        "Excluded non-natural or uncertain-context records: "
        f"{int((~df['sensitivity_model_candidate_broad']).sum())}"
    )
    print("Output folder:")
    print(OUTPUT_DIR)
    print("=" * 78 + "\n")


if __name__ == "__main__":
    main()
