"""
08_create_global_background_area_weighted.py

Purpose
-------
Create a reproducible 10,000-point background sample for the GLOBAL primary
Camellia japonica MaxEnt / maxnet model.

Background domain
-----------------
The background is sampled from the global valid environmental domain, defined
as raster cells with finite, non-NoData values in all five selected primary
predictors:

    Bio2, Bio6, Bio14, Bio18, Elv

Sampling is area-weighted on the geographic grid. Random candidate cells are
drawn with probability proportional to cos(latitude), then accepted only when
all selected predictor rasters contain valid values. This avoids over-
representing high-latitude 2.5 arc-minute cells merely because the WGS84 grid
contains more cells per unit ground area at high latitudes.

This is the primary global-background design. It is appropriate because the
study's inference domain is global potential climatic suitability, not a
restricted accessible-area M.

Input
-----
E:\PREDICTION TREE\Camellia_japonica\model_inputs\
Primary_Global_ClimateTopo_Current\
    Bio2.tif
    Bio6.tif
    Bio14.tif
    Bio18.tif
    Elv.tif

Output
------
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
08_global_background\
"""

from __future__ import annotations

import math
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


# =============================================================================
# 1. USER SETTINGS
# =============================================================================

PROJECT_DIR = Path(r"E:\PREDICTION TREE\Camellia_japonica")

CURRENT_MODEL_INPUT_DIR = (
    PROJECT_DIR
    / "model_inputs"
    / "Primary_Global_ClimateTopo_Current"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "08_global_background"
)

SELECTED_PREDICTORS = ["Bio2", "Bio6", "Bio14", "Bio18", "Elv"]

N_BACKGROUND_POINTS = 10000
RANDOM_SEED = 20260704

# Candidate cells are generated in batches. Increase only if processing is
# unexpectedly slow; 25,000 is usually efficient at 2.5 arc-minute resolution.
CANDIDATE_BATCH_SIZE = 25000
MAX_BATCHES = 100


# =============================================================================
# 2. PACKAGE CHECK
# =============================================================================

def ensure_package(import_name: str, pip_name: str | None = None) -> None:
    try:
        __import__(import_name)
    except ImportError:
        name = pip_name or import_name
        print(f"Installing missing package: {name}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", name]
        )


ensure_package("rasterio")

import rasterio


# =============================================================================
# 3. HELPER FUNCTIONS
# =============================================================================

def stop_with_message(message: str) -> None:
    raise RuntimeError(f"\nERROR:\n{message}")


def locate_predictor_rasters() -> Dict[str, Path]:
    rasters: Dict[str, Path] = {}

    for predictor in SELECTED_PREDICTORS:
        candidate = CURRENT_MODEL_INPUT_DIR / f"{predictor}.tif"

        if not candidate.exists():
            stop_with_message(
                f"Required selected predictor raster is missing:\n{candidate}\n\n"
                "Run Step 07 successfully before creating background points."
            )

        rasters[predictor] = candidate

    return rasters


def compare_grids(reference_path: Path, candidate_path: Path) -> bool:
    with rasterio.open(reference_path) as ref, rasterio.open(candidate_path) as src:
        return (
            ref.crs == src.crs
            and ref.width == src.width
            and ref.height == src.height
            and ref.transform.almost_equals(src.transform, precision=1e-9)
            and all(
                abs(a - b) <= 1e-9
                for a, b in zip(ref.bounds, src.bounds)
            )
        )


def valid_sample_values(
    raster_paths: Dict[str, Path],
    coordinates: List[tuple[float, float]]
) -> Dict[str, np.ndarray]:
    """
    Extract values for every selected raster at candidate cell centres.

    The raster grids are confirmed to be aligned, so the coordinates refer to
    identical environmental cells for all predictors.
    """
    extracted: Dict[str, np.ndarray] = {}

    for predictor, path in raster_paths.items():
        with rasterio.open(path) as src:
            values = np.asarray(
                [sample[0] for sample in src.sample(coordinates)],
                dtype=float
            )

            if src.nodata is not None and np.isfinite(src.nodata):
                values[np.isclose(values, src.nodata)] = np.nan

            values[~np.isfinite(values)] = np.nan
            extracted[predictor] = values

    return extracted


# =============================================================================
# 4. MAIN
# =============================================================================

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    raster_paths = locate_predictor_rasters()
    reference_path = raster_paths["Bio2"]

    # Validate raster alignment before sampling.
    for predictor, path in raster_paths.items():
        if not compare_grids(reference_path, path):
            stop_with_message(
                "Selected current predictor rasters are not fully aligned.\n"
                f"Reference: {reference_path}\n"
                f"Misaligned raster: {path}"
            )

    print("\n" + "=" * 82)
    print("GLOBAL AREA-WEIGHTED BACKGROUND SAMPLING: STEP 08")
    print("Selected predictors: " + ", ".join(SELECTED_PREDICTORS))
    print(f"Requested background points: {N_BACKGROUND_POINTS}")
    print(f"Fixed random seed: {RANDOM_SEED}")
    print(f"Reference raster: {reference_path}")
    print("=" * 82)

    rng = np.random.default_rng(RANDOM_SEED)

    with rasterio.open(reference_path) as ref:
        if ref.crs is None or not ref.crs.is_geographic:
            stop_with_message(
                "The reference raster must use a geographic CRS for the "
                "area-weighted latitude sampling implemented here.\n"
                f"Reference CRS: {ref.crs}"
            )

        height = ref.height
        width = ref.width
        transform = ref.transform

        # Latitude at the centre of every raster row.
        rows_all = np.arange(height, dtype=np.int64)
        row_latitudes = transform.f + (rows_all + 0.5) * transform.e

        # Area of a geographic-grid cell is approximately proportional to
        # cos(latitude). The clipping prevents rare numerical edge cases near
        # poles.
        row_weights = np.cos(np.deg2rad(row_latitudes))
        row_weights = np.clip(row_weights, 0, None)

        if not np.isfinite(row_weights).all() or row_weights.sum() <= 0:
            stop_with_message(
                "Could not calculate valid latitude-based area weights."
            )

        row_probabilities = row_weights / row_weights.sum()

        # Store accepted cells keyed by row/column to ensure one point per
        # environmental raster cell.
        accepted: Dict[tuple[int, int], Dict[str, float]] = {}

        for batch_number in range(1, MAX_BATCHES + 1):
            if len(accepted) >= N_BACKGROUND_POINTS:
                break

            candidate_rows = rng.choice(
                rows_all,
                size=CANDIDATE_BATCH_SIZE,
                replace=True,
                p=row_probabilities
            )

            candidate_cols = rng.integers(
                low=0,
                high=width,
                size=CANDIDATE_BATCH_SIZE,
                endpoint=False
            )

            # Remove duplicates within the batch and cells already accepted.
            candidate_pairs = []
            seen_in_batch = set()

            for row, col in zip(candidate_rows, candidate_cols):
                key = (int(row), int(col))

                if key not in accepted and key not in seen_in_batch:
                    candidate_pairs.append(key)
                    seen_in_batch.add(key)

            if not candidate_pairs:
                continue

            rows = np.asarray(
                [pair[0] for pair in candidate_pairs],
                dtype=np.int64
            )
            cols = np.asarray(
                [pair[1] for pair in candidate_pairs],
                dtype=np.int64
            )

            xs, ys = rasterio.transform.xy(
                transform,
                rows,
                cols,
                offset="center"
            )

            coordinates = list(zip(xs, ys))

            extracted = valid_sample_values(raster_paths, coordinates)

            valid_mask = np.ones(len(candidate_pairs), dtype=bool)

            for predictor in SELECTED_PREDICTORS:
                valid_mask &= np.isfinite(extracted[predictor])

            for index, is_valid in enumerate(valid_mask):
                if not is_valid:
                    continue

                row, col = candidate_pairs[index]
                key = (row, col)

                accepted[key] = {
                    "decimalLongitude": float(xs[index]),
                    "decimalLatitude": float(ys[index]),
                    **{
                        predictor: float(extracted[predictor][index])
                        for predictor in SELECTED_PREDICTORS
                    }
                }

                if len(accepted) >= N_BACKGROUND_POINTS:
                    break

            print(
                f"Batch {batch_number:02d}: "
                f"{len(accepted):,} / {N_BACKGROUND_POINTS:,} valid unique "
                "background cells retained"
            )

    if len(accepted) < N_BACKGROUND_POINTS:
        stop_with_message(
            "The script could not obtain the requested number of valid "
            "background cells.\n"
            f"Requested: {N_BACKGROUND_POINTS}\n"
            f"Obtained: {len(accepted)}\n"
            "Increase MAX_BATCHES or inspect current predictor raster coverage."
        )

    records = []

    for index, ((row, col), values) in enumerate(
        accepted.items(),
        start=1
    ):
        records.append(
            {
                "background_id": f"bg_{index:05d}",
                "decimalLongitude": values["decimalLongitude"],
                "decimalLatitude": values["decimalLatitude"],
                "reference_raster_row": row,
                "reference_raster_col": col,
                "reference_raster_cell_id": f"{row}_{col}",
                **{
                    predictor: values[predictor]
                    for predictor in SELECTED_PREDICTORS
                },
            }
        )

    background = pd.DataFrame(records).sort_values(
        "background_id"
    ).reset_index(drop=True)

    # =============================================================================
    # 5. EXPORTS
    # =============================================================================

    main_output = (
        OUTPUT_DIR
        / "08_global_area_weighted_background_n10000.csv"
    )

    background.to_csv(
        main_output,
        index=False,
        encoding="utf-8-sig"
    )

    latitude_summary = pd.DataFrame(
        {
            "metric": [
                "minimum_latitude",
                "maximum_latitude",
                "mean_latitude",
                "median_latitude",
                "northern_hemisphere_points",
                "southern_hemisphere_points",
            ],
            "value": [
                background["decimalLatitude"].min(),
                background["decimalLatitude"].max(),
                background["decimalLatitude"].mean(),
                background["decimalLatitude"].median(),
                int((background["decimalLatitude"] >= 0).sum()),
                int((background["decimalLatitude"] < 0).sum()),
            ],
        }
    )

    latitude_summary.to_csv(
        OUTPUT_DIR / "08_background_latitude_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    settings = pd.DataFrame(
        {
            "parameter": [
                "background_domain",
                "sampling_design",
                "requested_background_points",
                "retained_background_points",
                "random_seed",
                "reference_raster",
                "selected_predictors",
                "cell_deduplication",
                "inference_scope",
            ],
            "value": [
                (
                    "Global valid environmental domain: cells with finite, "
                    "non-NoData values in Bio2, Bio6, Bio14, Bio18, and Elv"
                ),
                (
                    "Area-weighted random sampling on geographic raster grid; "
                    "row selection proportional to cos(latitude)"
                ),
                N_BACKGROUND_POINTS,
                len(background),
                RANDOM_SEED,
                str(reference_path),
                "; ".join(SELECTED_PREDICTORS),
                "One background point per reference raster cell",
                "Global potential climatic suitability",
            ],
        }
    )

    settings.to_csv(
        OUTPUT_DIR / "08_global_background_settings.csv",
        index=False,
        encoding="utf-8-sig"
    )

    note = [
        "The primary global model uses an area-weighted random background from the valid global environmental domain.",
        "The background is not restricted to an accessible-area M because the study's inference domain is global potential climatic suitability.",
        "Sampling is area-weighted to avoid a high-latitude bias caused by equal-angle WGS84 raster cells.",
        "The resulting 10,000 points will be used consistently for ENMeval tuning and final model calibration.",
    ]

    (OUTPUT_DIR / "08_global_background_interpretation_note.txt").write_text(
        "\n".join(note),
        encoding="utf-8"
    )

    print("\n" + "=" * 82)
    print("Global background sampling completed successfully")
    print(f"Background points retained: {len(background):,}")
    print("Main output:")
    print(main_output)
    print("=" * 82 + "\n")


if __name__ == "__main__":
    main()
