"""
03_sample_ESA_WorldCover_direct_VSCode.py

Purpose
-------
Extract ESA WorldCover v200 (2021) land-cover context directly from VS Code
through the Google Earth Engine Python API. No Earth Engine Code Editor,
Asset upload, or Google Drive export is required.

Input
-----
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
02_coordinatecleaner\02_coordinatecleaner_retained_candidates.csv

Output
------
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
03_worldcover\03_ESA_WorldCover_context.csv

Methodological interpretation
------------------------------
This script identifies "natural-land-cover context candidates"; it does not
prove that each occurrence is a wild population. Metadata, CoordinateCleaner,
LULC screening, and subsequent 5 km thinning are combined to produce the
final occurrence dataset.

Earth Engine requirements
-------------------------
1. Install:
   pip install earthengine-api pandas numpy

2. The first run will open an authentication process if Earth Engine has not
   already been authenticated for this Python environment.

3. Change EE_PROJECT below if your Earth Engine Cloud Project ID differs.
"""

from __future__ import annotations

import sys
import time
import math
import subprocess
from pathlib import Path
from typing import List, Dict, Any

import numpy as np
import pandas as pd


# ============================================================================
# 1. USER SETTINGS
# ============================================================================

PROJECT_DIR = Path(r"E:\PREDICTION TREE\Camellia_japonica")

INPUT_CSV = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "02_coordinatecleaner"
    / "02_coordinatecleaner_retained_candidates.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "03_worldcover"
)

OUTPUT_CSV = OUTPUT_DIR / "03_ESA_WorldCover_context.csv"
SUMMARY_CSV = OUTPUT_DIR / "03_ESA_WorldCover_summary.csv"
SETTINGS_CSV = OUTPUT_DIR / "03_ESA_WorldCover_settings.csv"

# Your active Google Earth Engine Cloud Project ID.
# It commonly matches the project component in an asset path such as:
# projects/ee-bbthien98/assets/...
EE_PROJECT = "project-71a05c44-154c-4b45-89d"

# Earth Engine dataset: ESA WorldCover 2021 (v200, 10 m).
WORLDCOVER_ASSET = "ESA/WorldCover/v200/2021"

# Global processing is performed in batches to avoid oversized requests.
BATCH_SIZE = 100
SCALE_M = 10
TILE_SCALE = 4

# Conservative LULC screen used later by the workflow.
STRICT_NATURAL_MIN_TREE_SHRUB = 0.50
STRICT_NATURAL_MAX_CROPLAND = 0.20
STRICT_NATURAL_MAX_BUILT_UP = 0.05
STRICT_NATURAL_MAX_WATER = 0.10

BROAD_NATURAL_MIN_TREE_SHRUB = 0.50
BROAD_NATURAL_MAX_CROPLAND = 0.30
BROAD_NATURAL_MAX_BUILT_UP = 0.10
BROAD_NATURAL_MAX_WATER = 0.15

# Retry Earth Engine calls that temporarily fail.
MAX_RETRIES = 4
RETRY_WAIT_SECONDS = 10


# ============================================================================
# 2. PACKAGE AND EARTH ENGINE SETUP
# ============================================================================

def ensure_earthengine_api() -> None:
    """Install Earth Engine API only when it is absent."""
    try:
        import ee  # noqa: F401
    except ImportError:
        print("\nInstalling earthengine-api. This may take a few minutes...")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "earthengine-api"]
        )


def initialize_earth_engine():
    """Authenticate when needed, then initialize Earth Engine."""
    ensure_earthengine_api()

    import ee

    try:
        ee.Initialize(project=EE_PROJECT)
        print(f"Earth Engine initialized with Cloud Project: {EE_PROJECT}")
    except Exception as first_error:
        print("\nEarth Engine authentication is required.")
        print("A browser-based login process will open or a link will be shown.")
        print("Use the Google account that has access to Earth Engine.\n")

        try:
            ee.Authenticate()
            ee.Initialize(project=EE_PROJECT)
            print(f"Earth Engine initialized with Cloud Project: {EE_PROJECT}")
        except Exception as second_error:
            raise RuntimeError(
                "\nCould not authenticate or initialize Google Earth Engine.\n"
                f"Cloud Project configured in this script: {EE_PROJECT}\n\n"
                "Check that this is your active Earth Engine Cloud Project. "
                "If it differs, edit EE_PROJECT near the top of the script.\n\n"
                f"Initial error: {first_error}\n"
                f"Authentication error: {second_error}"
            ) from second_error

    return ee


# ============================================================================
# 3. DATA VALIDATION AND HELPER FUNCTIONS
# ============================================================================

WORLD_COVER_CLASSES = {
    10: "Tree_cover",
    20: "Shrubland",
    30: "Grassland",
    40: "Cropland",
    50: "Built_up",
    60: "Bare_sparse_vegetation",
    70: "Snow_ice",
    80: "Permanent_water_bodies",
    90: "Herbaceous_wetland",
    95: "Mangroves",
    100: "Moss_lichen",
}


def stop_with_message(message: str) -> None:
    raise RuntimeError(f"\nERROR:\n{message}")


def retry_get_info(ee_object, description: str) -> Dict[str, Any]:
    """Retrieve a small batch of Earth Engine results with retry logic."""
    last_error = None

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return ee_object.getInfo()
        except Exception as error:
            last_error = error

            if attempt == MAX_RETRIES:
                break

            wait_seconds = RETRY_WAIT_SECONDS * attempt

            print(
                f"  Warning: Earth Engine request failed for {description} "
                f"(attempt {attempt}/{MAX_RETRIES}). "
                f"Retrying after {wait_seconds} seconds..."
            )

            time.sleep(wait_seconds)

    raise RuntimeError(
        f"Earth Engine request failed after {MAX_RETRIES} attempts: "
        f"{description}\nLast error: {last_error}"
    ) from last_error


def read_input_points() -> pd.DataFrame:
    if not INPUT_CSV.exists():
        stop_with_message(
            "Step 02 output was not found:\n"
            f"{INPUT_CSV}\n\n"
            "Run 02_coordinatecleaner_global_without_sea.R first."
        )

    df = pd.read_csv(INPUT_CSV, low_memory=False)

    required_columns = {
        "record_id",
        "decimalLongitude",
        "decimalLatitude",
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        stop_with_message(
            "Input table lacks required columns:\n"
            + ", ".join(sorted(missing_columns))
        )

    df["decimalLongitude"] = pd.to_numeric(
        df["decimalLongitude"],
        errors="coerce"
    )

    df["decimalLatitude"] = pd.to_numeric(
        df["decimalLatitude"],
        errors="coerce"
    )

    df = df[
        df["decimalLongitude"].between(-180, 180)
        & df["decimalLatitude"].between(-90, 90)
    ].copy()

    if df.empty:
        stop_with_message("No valid geographic coordinates remain in input.")

    if "lulc_buffer_m" not in df.columns:
        # The Step 01 script normally supplies this column. The fallback is
        # appropriate for records without a reported coordinate uncertainty.
        df["lulc_buffer_m"] = 250.0
    else:
        df["lulc_buffer_m"] = pd.to_numeric(
            df["lulc_buffer_m"],
            errors="coerce"
        ).fillna(250.0).clip(lower=100.0, upper=1000.0)

    df = df.drop_duplicates(subset="record_id").reset_index(drop=True)

    return df


def dataframe_to_ee_features(ee, batch: pd.DataFrame):
    """Convert one local dataframe batch into an Earth Engine FeatureCollection."""
    feature_list = []

    for row in batch.itertuples(index=False):
        properties = {
            "record_id": str(row.record_id),
            "decimalLongitude": float(row.decimalLongitude),
            "decimalLatitude": float(row.decimalLatitude),
            "lulc_buffer_m": float(row.lulc_buffer_m),
        }

        point = ee.Geometry.Point(
            [float(row.decimalLongitude), float(row.decimalLatitude)]
        )

        feature_list.append(ee.Feature(point, properties))

    return ee.FeatureCollection(feature_list)


def features_to_property_dataframe(feature_info: Dict[str, Any]) -> pd.DataFrame:
    """Convert Earth Engine FeatureCollection getInfo output into a dataframe."""
    records: List[Dict[str, Any]] = []

    for feature in feature_info.get("features", []):
        properties = feature.get("properties", {})
        records.append(properties)

    return pd.DataFrame.from_records(records)


# ============================================================================
# 4. EARTH ENGINE EXTRACTION
# ============================================================================

def extract_worldcover_context_for_batch(ee, worldcover, class_stack, batch):
    """
    Extract:
    - centre-pixel ESA WorldCover class;
    - 11 land-cover fractions in a local uncertainty-aware buffer.
    """
    point_features = dataframe_to_ee_features(ee, batch)

    # Central 10 m pixel class.
    centre_samples = worldcover.sampleRegions(
        collection=point_features,
        properties=[
            "record_id",
            "decimalLongitude",
            "decimalLatitude",
            "lulc_buffer_m",
        ],
        scale=SCALE_M,
        geometries=False,
        tileScale=TILE_SCALE,
    )

    # Local buffer based on positional uncertainty. The output geometry is
    # discarded later, so only point ID and summaries are retained.
    buffered_features = point_features.map(
        lambda feature: ee.Feature(
            feature.geometry().buffer(ee.Number(feature.get("lulc_buffer_m"))),
            feature.toDictionary()
        )
    )

    buffer_fractions = class_stack.reduceRegions(
        collection=buffered_features,
        reducer=ee.Reducer.mean(),
        scale=SCALE_M,
        crs=worldcover.projection(),
        tileScale=TILE_SCALE,
    )

    centre_info = retry_get_info(
        centre_samples,
        f"central-pixel classification for batch beginning record {batch.iloc[0]['record_id']}"
    )

    fractions_info = retry_get_info(
        buffer_fractions,
        f"buffer land-cover fractions for batch beginning record {batch.iloc[0]['record_id']}"
    )

    centre_df = features_to_property_dataframe(centre_info)
    fractions_df = features_to_property_dataframe(fractions_info)

    # The centre class is stored in the Map band.
    if "Map" in centre_df.columns:
        centre_df = centre_df.rename(
            columns={"Map": "worldcover_center_code"}
        )

    return centre_df, fractions_df


def derive_lulc_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Add labels and strict/broad natural-context screening flags locally."""
    fraction_columns = [
        "frac_tree_cover",
        "frac_shrubland",
        "frac_grassland",
        "frac_cropland",
        "frac_built_up",
        "frac_bare_sparse",
        "frac_snow_ice",
        "frac_permanent_water",
        "frac_herbaceous_wetland",
        "frac_mangroves",
        "frac_moss_lichen",
    ]

    for column in fraction_columns:
        if column not in df.columns:
            df[column] = np.nan

        df[column] = pd.to_numeric(df[column], errors="coerce")

    df["worldcover_center_code"] = pd.to_numeric(
        df.get("worldcover_center_code"),
        errors="coerce"
    )

    df["worldcover_center_label"] = df["worldcover_center_code"].map(
        WORLD_COVER_CLASSES
    ).fillna("Unknown_or_missing")

    df["frac_tree_shrub_natural"] = (
        df["frac_tree_cover"].fillna(0)
        + df["frac_shrubland"].fillna(0)
    )

    df["center_is_tree_or_shrub"] = df["worldcover_center_code"].isin(
        [10, 20]
    )

    df["strict_natural_context_candidate"] = (
        df["center_is_tree_or_shrub"]
        & (df["frac_tree_shrub_natural"] >= STRICT_NATURAL_MIN_TREE_SHRUB)
        & (df["frac_cropland"].fillna(1) <= STRICT_NATURAL_MAX_CROPLAND)
        & (df["frac_built_up"].fillna(1) <= STRICT_NATURAL_MAX_BUILT_UP)
        & (
            df["frac_permanent_water"].fillna(1)
            <= STRICT_NATURAL_MAX_WATER
        )
    )

    df["broad_natural_context_candidate"] = (
        (df["frac_tree_shrub_natural"] >= BROAD_NATURAL_MIN_TREE_SHRUB)
        & (df["frac_cropland"].fillna(1) <= BROAD_NATURAL_MAX_CROPLAND)
        & (df["frac_built_up"].fillna(1) <= BROAD_NATURAL_MAX_BUILT_UP)
        & (
            df["frac_permanent_water"].fillna(1)
            <= BROAD_NATURAL_MAX_WATER
        )
    )

    df["LULC_screen_rule_strict"] = (
        "center Tree cover/Shrubland AND tree+shrub fraction >=0.50 "
        "AND cropland <=0.20 AND built-up <=0.05 AND water <=0.10"
    )

    df["LULC_screen_rule_broad"] = (
        "tree+shrub fraction >=0.50 AND cropland <=0.30 "
        "AND built-up <=0.10 AND water <=0.15"
    )

    return df


# ============================================================================
# 5. MAIN
# ============================================================================

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print("\n" + "=" * 76)
    print("GLOBAL ESA WORLDCOVER EXTRACTION: STEP 03")
    print("Execution mode: VS Code / local Python + Earth Engine API")
    print(f"Input: {INPUT_CSV}")
    print("=" * 76)

    points = read_input_points()
    ee = initialize_earth_engine()

    worldcover = ee.Image(WORLDCOVER_ASSET).select("Map")

    # For each class, the average of a binary raster within a point buffer
    # represents the fraction of that buffer occupied by the class.
    class_stack = ee.Image.cat(
        [
            worldcover.eq(10).rename("frac_tree_cover"),
            worldcover.eq(20).rename("frac_shrubland"),
            worldcover.eq(30).rename("frac_grassland"),
            worldcover.eq(40).rename("frac_cropland"),
            worldcover.eq(50).rename("frac_built_up"),
            worldcover.eq(60).rename("frac_bare_sparse"),
            worldcover.eq(70).rename("frac_snow_ice"),
            worldcover.eq(80).rename("frac_permanent_water"),
            worldcover.eq(90).rename("frac_herbaceous_wetland"),
            worldcover.eq(95).rename("frac_mangroves"),
            worldcover.eq(100).rename("frac_moss_lichen"),
        ]
    )

    n_batches = math.ceil(len(points) / BATCH_SIZE)

    centre_parts: List[pd.DataFrame] = []
    fraction_parts: List[pd.DataFrame] = []

    print(f"\nRecords to process: {len(points)}")
    print(f"Batch size: {BATCH_SIZE}")
    print(f"Total batches: {n_batches}\n")

    for batch_index, start in enumerate(
        range(0, len(points), BATCH_SIZE),
        start=1
    ):
        end = min(start + BATCH_SIZE, len(points))
        batch = points.iloc[start:end].copy()

        print(
            f"Processing batch {batch_index}/{n_batches} "
            f"({start + 1}–{end} of {len(points)})..."
        )

        centre_df, fraction_df = extract_worldcover_context_for_batch(
            ee=ee,
            worldcover=worldcover,
            class_stack=class_stack,
            batch=batch,
        )

        centre_parts.append(centre_df)
        fraction_parts.append(fraction_df)

    centre_all = pd.concat(centre_parts, ignore_index=True)
    fractions_all = pd.concat(fraction_parts, ignore_index=True)

    # Retain only one copy of common properties when merging.
    common_fraction_columns = {
        "record_id",
        "decimalLongitude",
        "decimalLatitude",
        "lulc_buffer_m",
    }

    fraction_measurement_columns = [
        column
        for column in fractions_all.columns
        if column not in common_fraction_columns
    ]

    fractions_for_merge = fractions_all[
        ["record_id"] + fraction_measurement_columns
    ].copy()

    output = points.merge(
        centre_all[
            [
                column
                for column in [
                    "record_id",
                    "worldcover_center_code",
                ]
                if column in centre_all.columns
            ]
        ],
        on="record_id",
        how="left",
        validate="one_to_one",
    )

    output = output.merge(
        fractions_for_merge,
        on="record_id",
        how="left",
        validate="one_to_one",
    )

    output = derive_lulc_metrics(output)

    # Keep the full Step 02 metadata plus all WorldCover context fields.
    output.to_csv(OUTPUT_CSV, index=False, encoding="utf-8-sig")

    summary = pd.DataFrame(
        {
            "metric": [
                "input_coordinatecleaner_candidates",
                "records_with_center_worldcover_class",
                "strict_natural_context_candidates",
                "broad_natural_context_candidates",
                "strict_tree_cover_center",
                "strict_shrubland_center",
            ],
            "n_records": [
                len(points),
                int(output["worldcover_center_code"].notna().sum()),
                int(output["strict_natural_context_candidate"].sum()),
                int(output["broad_natural_context_candidate"].sum()),
                int(
                    (
                        output["strict_natural_context_candidate"]
                        & output["worldcover_center_code"].eq(10)
                    ).sum()
                ),
                int(
                    (
                        output["strict_natural_context_candidate"]
                        & output["worldcover_center_code"].eq(20)
                    ).sum()
                ),
            ],
        }
    )

    summary.to_csv(SUMMARY_CSV, index=False, encoding="utf-8-sig")

    settings = pd.DataFrame(
        {
            "parameter": [
                "Earth_Engine_project",
                "WorldCover_asset",
                "WorldCover_year",
                "WorldCover_resolution_m",
                "batch_size",
                "buffer_definition",
                "strict_screen",
                "broad_screen",
            ],
            "value": [
                EE_PROJECT,
                WORLDCOVER_ASSET,
                "2021",
                SCALE_M,
                BATCH_SIZE,
                (
                    "lulc_buffer_m from Step 01: reported uncertainty clipped "
                    "to 100–1000 m; 250 m where uncertainty was not reported"
                ),
                (
                    "centre=Tree/Shrub; tree+shrub>=0.50; cropland<=0.20; "
                    "built-up<=0.05; permanent water<=0.10"
                ),
                (
                    "tree+shrub>=0.50; cropland<=0.30; built-up<=0.10; "
                    "permanent water<=0.15"
                ),
            ],
        }
    )

    settings.to_csv(SETTINGS_CSV, index=False, encoding="utf-8-sig")

    print("\n" + "=" * 76)
    print("ESA WorldCover extraction completed successfully")
    print(f"Input points: {len(points)}")
    print(
        "Strict natural-land-cover candidates: "
        f"{int(output['strict_natural_context_candidate'].sum())}"
    )
    print(
        "Broad natural-land-cover candidates: "
        f"{int(output['broad_natural_context_candidate'].sum())}"
    )
    print(f"Main output:\n{OUTPUT_CSV}")
    print("=" * 76 + "\n")


if __name__ == "__main__":
    main()
