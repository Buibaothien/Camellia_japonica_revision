"""
05_global_spatial_thinning_5km_and_gridcell_dedup.py

Purpose
-------
Apply reproducible global spatial thinning to the strict global
natural-land-cover occurrence candidates retained in Step 04.

This script performs two complementary controls:

1. Global 5 km spatial thinning
   - Great-circle (haversine) distance
   - 100 deterministic random-greedy replicates
   - The replicate retaining the largest number of records is selected
   - Fixed seed = 20260702

2. Environmental raster-cell deduplication
   - Retains only one point per cell of the reference 2.5 arc-minute raster
   - This prevents duplicated environmental predictor values at the native
     modelling resolution.

Input
-----
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
04_LULC_filter\04_global_strict_natural_landcover_candidates.csv

Reference raster
----------------
The script automatically searches for Bio1.tif inside:
E:\PREDICTION TREE\Camellia_japonica\environmental_current\

Outputs
-------
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
05_spatial_thinning_5km\

The final output for subsequent modelling is:
05_global_strict_natural_5km_thinned_one_per_raster_cell.csv
"""

from __future__ import annotations

import math
import subprocess
import sys
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
    / "04_LULC_filter"
    / "04_global_strict_natural_landcover_candidates.csv"
)

ENVIRONMENTAL_DIR = PROJECT_DIR / "environmental_current"

# Leave as None for automatic search. If needed, set the exact full path:
# REFERENCE_RASTER = Path(r"E:\PREDICTION TREE\Camellia_japonica\environmental_current\Bio1.tif")
REFERENCE_RASTER = None

OUTPUT_DIR = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "05_spatial_thinning_5km"
)

MIN_DISTANCE_KM = 5.0
N_RANDOM_GREEDY_REPLICATES = 100
RANDOM_SEED = 20260702

EARTH_RADIUS_KM = 6371.0088


# ============================================================================
# 2. PACKAGE CHECKS
# ============================================================================

def ensure_package(import_name: str, pip_name: str | None = None) -> None:
    """Install a required package only when it is missing."""
    try:
        __import__(import_name)
    except ImportError:
        package_to_install = pip_name or import_name
        print(f"Installing missing package: {package_to_install}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", package_to_install]
        )


ensure_package("sklearn", "scikit-learn")
ensure_package("rasterio")


from sklearn.neighbors import BallTree
import rasterio


# ============================================================================
# 3. HELPER FUNCTIONS
# ============================================================================

def stop_with_message(message: str) -> None:
    raise RuntimeError(f"\nERROR:\n{message}")


def locate_reference_raster() -> Path:
    """Locate Bio1.tif robustly, including case variations."""
    if REFERENCE_RASTER is not None:
        candidate = Path(REFERENCE_RASTER)

        if not candidate.exists():
            stop_with_message(
                "REFERENCE_RASTER was specified but not found:\n"
                f"{candidate}"
            )

        return candidate

    if not ENVIRONMENTAL_DIR.exists():
        stop_with_message(
            "Current environmental directory was not found:\n"
            f"{ENVIRONMENTAL_DIR}"
        )

    all_tifs = [
        path
        for path in ENVIRONMENTAL_DIR.rglob("*")
        if path.is_file() and path.suffix.lower() in {".tif", ".tiff"}
    ]

    preferred_names = {
        "bio1.tif",
        "bio1.tiff",
        "wc2.1_2.5m_bio_1.tif",
        "wc2.1_2.5m_bio_1.tiff",
    }

    exact_matches = [
        path for path in all_tifs
        if path.name.lower() in preferred_names
    ]

    if len(exact_matches) == 1:
        return exact_matches[0]

    bio1_matches = [
        path
        for path in all_tifs
        if path.stem.lower() == "bio1"
    ]

    if len(bio1_matches) == 1:
        return bio1_matches[0]

    if len(exact_matches) > 1 or len(bio1_matches) > 1:
        candidates = exact_matches if len(exact_matches) > 1 else bio1_matches
        stop_with_message(
            "More than one candidate Bio1 raster was found. "
            "Set REFERENCE_RASTER manually near the top of the script.\n\n"
            + "\n".join(str(path) for path in candidates)
        )

    stop_with_message(
        "No Bio1 reference raster was found inside:\n"
        f"{ENVIRONMENTAL_DIR}\n\n"
        "Set REFERENCE_RASTER manually near the top of the script."
    )


def standardise_input(df: pd.DataFrame) -> pd.DataFrame:
    required_columns = {
        "record_id",
        "decimalLongitude",
        "decimalLatitude",
    }

    missing_columns = required_columns - set(df.columns)

    if missing_columns:
        stop_with_message(
            "Input file lacks required columns:\n"
            + ", ".join(sorted(missing_columns))
        )

    out = df.copy()

    out["decimalLongitude"] = pd.to_numeric(
        out["decimalLongitude"],
        errors="coerce"
    )

    out["decimalLatitude"] = pd.to_numeric(
        out["decimalLatitude"],
        errors="coerce"
    )

    out = out[
        out["decimalLongitude"].between(-180, 180)
        & out["decimalLatitude"].between(-90, 90)
    ].copy()

    if out.empty:
        stop_with_message("No valid geographic coordinates remain in the input.")

    # Prioritise lower reported coordinate uncertainty when ties occur.
    if "coordinateUncertaintyInMeters_numeric" in out.columns:
        out["_uncertainty_sort"] = pd.to_numeric(
            out["coordinateUncertaintyInMeters_numeric"],
            errors="coerce"
        ).fillna(np.inf)
    else:
        out["_uncertainty_sort"] = np.inf

    # Guarantee unique record IDs.
    out["record_id"] = out["record_id"].astype(str)

    if out["record_id"].duplicated().any():
        stop_with_message(
            "Duplicate record_id values were found in the Step 04 input. "
            "This should not happen. Inspect the input table before continuing."
        )

    return out.reset_index(drop=True)


def build_neighbor_lists(latitude, longitude, distance_km):
    """
    Build neighbors within the specified minimum great-circle distance.

    Coordinates passed to BallTree must be [latitude, longitude] in radians.
    """
    coords_rad = np.deg2rad(
        np.column_stack([latitude, longitude])
    )

    tree = BallTree(coords_rad, metric="haversine")

    radius_radians = distance_km / EARTH_RADIUS_KM

    neighbor_indices = tree.query_radius(
        coords_rad,
        r=radius_radians,
        return_distance=False,
        sort_results=False,
    )

    return neighbor_indices


def run_random_greedy_thinning(neighbor_lists, seed, n_replicates):
    """
    Run random-greedy thinning repeatedly and return the largest retained set.

    Every replicate is reproducible because the master random seed is fixed.
    """
    n_points = len(neighbor_lists)
    rng = np.random.default_rng(seed)

    best_selected = None
    best_count = -1
    results = []

    for replicate in range(1, n_replicates + 1):
        order = rng.permutation(n_points)
        blocked = np.zeros(n_points, dtype=bool)
        selected = []

        for index in order:
            if not blocked[index]:
                selected.append(index)
                blocked[neighbor_lists[index]] = True

        selected_count = len(selected)

        results.append(
            {
                "replicate": replicate,
                "n_retained": selected_count,
                "n_removed_due_to_5km": n_points - selected_count,
            }
        )

        if selected_count > best_count:
            best_count = selected_count
            best_selected = np.asarray(selected, dtype=int)

    replicate_summary = pd.DataFrame(results)

    selected_replicate = int(
        replicate_summary.loc[
            replicate_summary["n_retained"].idxmax(),
            "replicate"
        ]
    )

    return best_selected, selected_replicate, replicate_summary


def extract_reference_raster_cell_ids(
    df: pd.DataFrame,
    raster_path: Path
) -> pd.DataFrame:
    """
    Assign row/column/cell ID in the reference environmental raster.

    Points outside the raster extent or falling on NoData cells are marked
    invalid and removed from the final SDM occurrence dataset.
    """
    output = df.copy()

    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
        raster_transform = src.transform
        raster_height = src.height
        raster_width = src.width
        nodata_value = src.nodata

        if raster_crs is None:
            stop_with_message(
                f"Reference raster has no CRS:\n{raster_path}"
            )

        # The workflow expects WGS84 geographic input for WorldClim.
        # Coordinate transformation is added only when needed.
        if raster_crs.is_geographic:
            xs = output["decimalLongitude"].to_numpy()
            ys = output["decimalLatitude"].to_numpy()
        else:
            try:
                from pyproj import Transformer
            except ImportError:
                subprocess.check_call(
                    [sys.executable, "-m", "pip", "install", "pyproj"]
                )
                from pyproj import Transformer

            transformer = Transformer.from_crs(
                "EPSG:4326",
                raster_crs,
                always_xy=True
            )

            xs, ys = transformer.transform(
                output["decimalLongitude"].to_numpy(),
                output["decimalLatitude"].to_numpy()
            )

        rows, cols = rasterio.transform.rowcol(
            raster_transform,
            xs,
            ys,
            op=math.floor
        )

        rows = np.asarray(rows, dtype=int)
        cols = np.asarray(cols, dtype=int)

        inside_extent = (
            (rows >= 0)
            & (rows < raster_height)
            & (cols >= 0)
            & (cols < raster_width)
        )

        output["reference_raster_row"] = rows
        output["reference_raster_col"] = cols
        output["inside_reference_raster_extent"] = inside_extent

        cell_values = np.full(len(output), np.nan, dtype=float)

        valid_indices = np.where(inside_extent)[0]

        if len(valid_indices) > 0:
            sample_coords = [
                (float(xs[index]), float(ys[index]))
                for index in valid_indices
            ]

            sampled = np.asarray(
                [value[0] for value in src.sample(sample_coords)],
                dtype=float
            )

            cell_values[valid_indices] = sampled

        output["reference_raster_value"] = cell_values

        if nodata_value is not None and np.isfinite(nodata_value):
            valid_data = (
                output["inside_reference_raster_extent"]
                & np.isfinite(output["reference_raster_value"])
                & ~np.isclose(
                    output["reference_raster_value"],
                    nodata_value
                )
            )
        else:
            valid_data = (
                output["inside_reference_raster_extent"]
                & np.isfinite(output["reference_raster_value"])
            )

        output["valid_reference_raster_cell"] = valid_data

        output["reference_raster_cell_id"] = np.where(
            output["valid_reference_raster_cell"],
            (
                output["reference_raster_row"].astype(str)
                + "_"
                + output["reference_raster_col"].astype(str)
            ),
            np.nan
        )

    return output


# ============================================================================
# 4. MAIN
# ============================================================================

def main() -> None:
    if not INPUT_CSV.exists():
        stop_with_message(
            "Step 04 strict primary candidate file was not found:\n"
            f"{INPUT_CSV}\n\n"
            "Run 04_filter_global_natural_landcover_candidates.py first."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    reference_raster = locate_reference_raster()

    print("\n" + "=" * 80)
    print("GLOBAL 5 KM SPATIAL THINNING + RASTER-CELL DEDUPLICATION: STEP 05")
    print(f"Input occurrence file:\n{INPUT_CSV}")
    print(f"Reference raster:\n{reference_raster}")
    print(f"Minimum great-circle distance: {MIN_DISTANCE_KM:.1f} km")
    print(f"Random-greedy thinning replicates: {N_RANDOM_GREEDY_REPLICATES}")
    print(f"Fixed random seed: {RANDOM_SEED}")
    print("=" * 80)

    raw = pd.read_csv(INPUT_CSV, low_memory=False)
    input_df = standardise_input(raw)

    # ------------------------------------------------------------------------
    # 4.1 Global 5 km thinning
    # ------------------------------------------------------------------------

    neighbor_lists = build_neighbor_lists(
        latitude=input_df["decimalLatitude"].to_numpy(),
        longitude=input_df["decimalLongitude"].to_numpy(),
        distance_km=MIN_DISTANCE_KM,
    )

    selected_indices, selected_replicate, replicate_summary = (
        run_random_greedy_thinning(
            neighbor_lists=neighbor_lists,
            seed=RANDOM_SEED,
            n_replicates=N_RANDOM_GREEDY_REPLICATES,
        )
    )

    input_df["retained_after_5km_thinning"] = False
    input_df.loc[
        input_df.index[selected_indices],
        "retained_after_5km_thinning"
    ] = True

    input_df["selected_thinning_replicate"] = selected_replicate

    thinned_df = (
        input_df.loc[input_df["retained_after_5km_thinning"]]
        .copy()
        .sort_values("record_id")
        .reset_index(drop=True)
    )

    # ------------------------------------------------------------------------
    # 4.2 One point per environmental raster cell
    # ------------------------------------------------------------------------

    thinned_with_grid = extract_reference_raster_cell_ids(
        thinned_df,
        reference_raster
    )

    valid_grid = thinned_with_grid.loc[
        thinned_with_grid["valid_reference_raster_cell"]
    ].copy()

    # Preserve the best located record when more than one retained point falls
    # within the same 2.5 arc-minute environmental cell.
    valid_grid = valid_grid.sort_values(
        [
            "reference_raster_cell_id",
            "_uncertainty_sort",
            "record_id"
        ],
        kind="mergesort"
    )

    valid_grid["retained_after_raster_cell_dedup"] = ~valid_grid.duplicated(
        "reference_raster_cell_id",
        keep="first"
    )

    final_df = (
        valid_grid.loc[valid_grid["retained_after_raster_cell_dedup"]]
        .copy()
        .sort_values("record_id")
        .reset_index(drop=True)
    )

    # Add final retention information back to all candidates for auditing.
    all_audit = input_df.copy()

    selected_record_ids = set(thinned_df["record_id"])
    final_record_ids = set(final_df["record_id"])

    all_audit["retained_after_raster_cell_dedup"] = all_audit[
        "record_id"
    ].isin(final_record_ids)

    all_audit["final_primary_occurrence_keep"] = all_audit[
        "record_id"
    ].isin(final_record_ids)

    all_audit["thinning_status"] = np.select(
        [
            ~all_audit["retained_after_5km_thinning"],
            (
                all_audit["retained_after_5km_thinning"]
                & ~all_audit["retained_after_raster_cell_dedup"]
            ),
            all_audit["final_primary_occurrence_keep"],
        ],
        [
            "removed_by_5km_thinning",
            "removed_as_duplicate_reference_raster_cell",
            "retained_final_primary_occurrence",
        ],
        default="not_classified",
    )

    # ------------------------------------------------------------------------
    # 4.3 Export files
    # ------------------------------------------------------------------------

    final_export_columns = [
        column
        for column in final_df.columns
        if column not in {"_uncertainty_sort"}
    ]

    final_df[final_export_columns].to_csv(
        OUTPUT_DIR
        / "05_global_strict_natural_5km_thinned_one_per_raster_cell.csv",
        index=False,
        encoding="utf-8-sig"
    )

    thinned_export_columns = [
        column
        for column in thinned_with_grid.columns
        if column not in {"_uncertainty_sort"}
    ]

    thinned_with_grid[thinned_export_columns].to_csv(
        OUTPUT_DIR
        / "05_global_strict_natural_after_5km_thinning_before_grid_dedup.csv",
        index=False,
        encoding="utf-8-sig"
    )

    all_audit_export_columns = [
        column
        for column in all_audit.columns
        if column not in {"_uncertainty_sort"}
    ]

    all_audit[all_audit_export_columns].to_csv(
        OUTPUT_DIR
        / "05_global_strict_natural_all_candidates_thinning_audit.csv",
        index=False,
        encoding="utf-8-sig"
    )

    replicate_summary.to_csv(
        OUTPUT_DIR / "05_thinning_replicate_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    cell_dedup_summary = pd.DataFrame(
        {
            "metric": [
                "records_after_5km_thinning",
                "records_inside_reference_raster_with_valid_data",
                "records_outside_or_nodata_reference_raster",
                "duplicate_records_removed_within_reference_raster_cell",
                "final_one_per_raster_cell_occurrences",
                "unique_reference_raster_cells",
            ],
            "n_records": [
                len(thinned_with_grid),
                int(thinned_with_grid["valid_reference_raster_cell"].sum()),
                int((~thinned_with_grid["valid_reference_raster_cell"]).sum()),
                int(
                    valid_grid["retained_after_raster_cell_dedup"]
                    .eq(False)
                    .sum()
                ),
                len(final_df),
                final_df["reference_raster_cell_id"].nunique(),
            ],
        }
    )

    cell_dedup_summary.to_csv(
        OUTPUT_DIR / "05_raster_cell_dedup_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    overall_summary = pd.DataFrame(
        {
            "stage": [
                "strict_natural_landcover_input",
                "after_5km_global_spatial_thinning",
                "after_reference_raster_validity_check",
                "after_one_record_per_reference_raster_cell",
            ],
            "n_records": [
                len(input_df),
                len(thinned_with_grid),
                int(thinned_with_grid["valid_reference_raster_cell"].sum()),
                len(final_df),
            ],
        }
    )

    overall_summary.to_csv(
        OUTPUT_DIR / "05_spatial_thinning_overall_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    settings = pd.DataFrame(
        {
            "parameter": [
                "input_dataset",
                "minimum_distance_km",
                "distance_method",
                "random_greedy_replicates",
                "random_seed",
                "selected_replicate",
                "reference_raster",
                "reference_raster_resolution",
                "grid_cell_deduplication",
                "final_dataset_interpretation",
            ],
            "value": [
                "Strict global natural-land-cover occurrence candidates",
                MIN_DISTANCE_KM,
                "Great-circle haversine distance",
                N_RANDOM_GREEDY_REPLICATES,
                RANDOM_SEED,
                selected_replicate,
                str(reference_raster),
                "Reference raster native resolution, expected 2.5 arc-minute",
                "One retained record per valid reference raster cell",
                (
                    "Global natural-land-cover occurrence candidates after "
                    "5 km thinning and environmental-cell deduplication"
                ),
            ],
        }
    )

    settings.to_csv(
        OUTPUT_DIR / "05_spatial_thinning_settings.csv",
        index=False,
        encoding="utf-8-sig"
    )

    interpretation_note = [
        "A minimum 5 km great-circle separation was applied globally.",
        "The largest retained set across 100 random-greedy replicates was selected using a fixed seed.",
        "A second control retained only one record per valid cell of the 2.5 arc-minute reference environmental raster.",
        "The final dataset should be called global natural-land-cover occurrence candidates, not confirmed wild occurrences.",
    ]

    (OUTPUT_DIR / "05_spatial_thinning_interpretation_note.txt").write_text(
        "\n".join(interpretation_note),
        encoding="utf-8"
    )

    # ------------------------------------------------------------------------
    # 4.4 Console report
    # ------------------------------------------------------------------------

    print("\n" + "=" * 80)
    print("Global 5 km thinning and raster-cell deduplication completed")
    print(f"Input strict LULC candidates: {len(input_df)}")
    print(f"Selected thinning replicate: {selected_replicate}")
    print(f"Retained after 5 km thinning: {len(thinned_with_grid)}")
    print(
        "Retained after one-point-per-raster-cell deduplication: "
        f"{len(final_df)}"
    )
    print("Final occurrence file:")
    print(
        OUTPUT_DIR
        / "05_global_strict_natural_5km_thinned_one_per_raster_cell.csv"
    )
    print("=" * 80 + "\n")


if __name__ == "__main__":
    main()
