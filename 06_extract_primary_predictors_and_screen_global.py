"""
06_extract_primary_predictors_and_screen_global.py

Purpose
-------
Prepare the environmental screening dataset for the global Camellia japonica
reanalysis after occurrence filtering and 5 km thinning.

This script:
1. reads the final strict global natural-land-cover occurrence candidates;
2. identifies the transferable primary predictors available for both current
   and future projections:
      Bio1–Bio19 + elevation (Elv);
3. explicitly excludes static soil, NPP, runoff, and UV layers from the
   PRIMARY future-transfer model;
4. checks raster-grid consistency against Bio1;
5. extracts current raster values at all final occurrence records;
6. draws a reproducible global terrestrial environmental sample from valid
   Bio1 cells and extracts all primary predictors there;
7. calculates missingness and Pearson correlations for:
      - occurrences;
      - global environmental screening sample;
8. exports tables needed to select the final low-collinearity predictor set.

Important
---------
This script does NOT choose the final predictors automatically. Selection must
combine:
- correlations,
- missingness,
- biological interpretability,
- availability of matched future layers,
- and later VIF screening.

Inputs
------
Occurrences:
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
05_spatial_thinning_5km\
05_global_strict_natural_5km_thinned_one_per_raster_cell.csv

Current predictors:
E:\PREDICTION TREE\Camellia_japonica\environmental_current\

Outputs
-------
E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs\
06_predictor_screening\
"""

from __future__ import annotations

import math
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd


# ============================================================================
# 1. USER SETTINGS
# ============================================================================

PROJECT_DIR = Path(r"E:\PREDICTION TREE\Camellia_japonica")

OCCURRENCE_CSV = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "05_spatial_thinning_5km"
    / "05_global_strict_natural_5km_thinned_one_per_raster_cell.csv"
)

CURRENT_ENVIRONMENTAL_DIR = PROJECT_DIR / "environmental_current"

OUTPUT_DIR = (
    PROJECT_DIR
    / "global_reanalysis_outputs"
    / "06_predictor_screening"
)

# Primary predictors retained for screening.
# All must have current and future matched data.
PRIMARY_CODES = [f"Bio{i}" for i in range(1, 20)] + ["Elv"]

# Variables not allowed in the primary future-transfer model because matched
# future projections are not available in this workflow. They are documented,
# not deleted, in the output report.
EXCLUDED_STATIC_OR_UNMATCHED_CODES = {
    "Sm": "Static soil moisture excluded from the primary future-transfer model",
    "Npp": "Net primary productivity excluded because no matched future projection is used",
    "Ar": "Annual runoff excluded because no matched future projection is used",
    "Soc": "Soil organic carbon excluded because no matched future projection is used",
    "Sph": "Soil pH excluded because no matched future projection is used",
    "AmUV": "Annual mean UV-B excluded because no matched future projection is used",
    "SeaUV": "UV-B seasonality excluded because no matched future projection is used",
    "HighUV": "Maximum UV-B excluded because no matched future projection is used",
    "LowUV": "Minimum UV-B excluded because no matched future projection is used",
}

# Number of globally distributed raster cells used only for predictor screening.
# This is not yet the final MaxEnt background sample.
GLOBAL_SCREENING_SAMPLE_SIZE = 20000
RANDOM_SEED = 20260703
CORRELATION_THRESHOLD = 0.80

# Use 0.5 to compare affine transforms and bounds conservatively. All current
# rasters are expected to be perfectly aligned after preprocessing.
GRID_TOLERANCE = 1e-9


# ============================================================================
# 2. PACKAGE CHECKS
# ============================================================================

def ensure_package(import_name: str, pip_name: str | None = None) -> None:
    try:
        __import__(import_name)
    except ImportError:
        install_name = pip_name or import_name
        print(f"Installing missing package: {install_name}")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", install_name]
        )


ensure_package("rasterio")
ensure_package("pyproj")

import rasterio
from pyproj import Transformer


# ============================================================================
# 3. HELPER FUNCTIONS
# ============================================================================

def stop_with_message(message: str) -> None:
    raise RuntimeError(f"\nERROR:\n{message}")


def normalise_filename(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", text.lower())


def candidate_raster_files(directory: Path) -> List[Path]:
    if not directory.exists():
        stop_with_message(
            "Current environmental directory was not found:\n"
            f"{directory}"
        )

    files = [
        path for path in directory.rglob("*")
        if path.is_file() and path.suffix.lower() in {".tif", ".tiff"}
    ]

    if not files:
        stop_with_message(
            "No .tif or .tiff files were found in:\n"
            f"{directory}"
        )

    return sorted(files)


def detect_bio_code(path: Path) -> str | None:
    """
    Detect Bio1–Bio19 from flexible filename styles:
      Bio1.tif
      bio_01.tif
      wc2.1_2.5m_bio_1.tif
    """
    stem = path.stem.lower()

    match = re.search(
        r"(?<![a-z0-9])bio[_\-\s]?0?([1-9]|1[0-9])(?![0-9])",
        stem
    )

    if not match:
        # Fallback after separator removal.
        compact = normalise_filename(stem)
        match = re.search(r"bio0?([1-9]|1[0-9])", compact)

    if match:
        number = int(match.group(1))
        if 1 <= number <= 19:
            return f"Bio{number}"

    return None


def detect_elevation(path: Path) -> bool:
    compact = normalise_filename(path.stem)

    return (
        compact in {"elv", "elev", "elevation", "dem"}
        or compact.startswith("elevation")
        or compact.startswith("worldclimelevation")
    )


def classify_nonprimary_raster(path: Path) -> str | None:
    compact = normalise_filename(path.stem)

    known_patterns = {
        "Sm": ["sm", "soilmoisture"],
        "Npp": ["npp", "netprimaryproductivity"],
        "Ar": ["ar", "annualrunoff", "runoff"],
        "Soc": ["soc", "soilorganiccarbon"],
        "Sph": ["sph", "soilph"],
        "AmUV": ["amuv", "annualmeanuv"],
        "SeaUV": ["seauv", "seasonalityuv"],
        "HighUV": ["highuv", "highestuv", "maximumuv"],
        "LowUV": ["lowuv", "lowestuv", "minimumuv"],
    }

    for code, patterns in known_patterns.items():
        if any(pattern in compact for pattern in patterns):
            return code

    return None


def choose_best_duplicate(candidates: List[Path], code: str) -> Path:
    """
    Select an exact canonical filename where possible. Otherwise stop rather
    than silently selecting a potentially wrong raster.
    """
    canonical_names = {
        f"{code.lower()}.tif",
        f"{code.lower()}.tiff",
    }

    exact = [
        path for path in candidates
        if path.name.lower() in canonical_names
    ]

    if len(exact) == 1:
        return exact[0]

    if len(candidates) == 1:
        return candidates[0]

    stop_with_message(
        f"More than one raster looks like {code}. "
        "Resolve this manually by moving duplicates out of the input folder.\n\n"
        + "\n".join(str(path) for path in candidates)
    )


def locate_primary_rasters() -> Tuple[Dict[str, Path], pd.DataFrame]:
    """
    Identify Bio1–Bio19 + elevation, and report all other raster layers.
    """
    files = candidate_raster_files(CURRENT_ENVIRONMENTAL_DIR)

    matches: Dict[str, List[Path]] = {code: [] for code in PRIMARY_CODES}
    excluded_rows = []

    for path in files:
        bio_code = detect_bio_code(path)

        if bio_code is not None:
            matches[bio_code].append(path)
            continue

        if detect_elevation(path):
            matches["Elv"].append(path)
            continue

        nonprimary_code = classify_nonprimary_raster(path)

        if nonprimary_code is not None:
            excluded_rows.append(
                {
                    "raster_file": str(path),
                    "detected_code": nonprimary_code,
                    "status": "excluded_from_primary_model",
                    "reason": EXCLUDED_STATIC_OR_UNMATCHED_CODES[nonprimary_code],
                }
            )
        else:
            excluded_rows.append(
                {
                    "raster_file": str(path),
                    "detected_code": "Unknown",
                    "status": "not_used_in_primary_model",
                    "reason": "Not a recognised Bio1–Bio19 or elevation raster",
                }
            )

    selected = {}

    for code in PRIMARY_CODES:
        if len(matches[code]) == 0:
            stop_with_message(
                f"Required primary raster {code} was not found in:\n"
                f"{CURRENT_ENVIRONMENTAL_DIR}"
            )

        selected[code] = choose_best_duplicate(matches[code], code)

    excluded_report = pd.DataFrame(
        excluded_rows,
        columns=["raster_file", "detected_code", "status", "reason"]
    )

    return selected, excluded_report


def compare_grid(reference_path: Path, candidate_path: Path) -> Dict[str, object]:
    """Compare CRS, dimensions, transform, and bounds against the reference."""
    with rasterio.open(reference_path) as ref, rasterio.open(candidate_path) as src:
        same_crs = ref.crs == src.crs
        same_width = ref.width == src.width
        same_height = ref.height == src.height
        same_transform = ref.transform.almost_equals(
            src.transform,
            precision=GRID_TOLERANCE
        )

        same_bounds = all(
            abs(a - b) <= GRID_TOLERANCE
            for a, b in zip(ref.bounds, src.bounds)
        )

        return {
            "reference_raster": str(reference_path),
            "raster_file": str(candidate_path),
            "reference_crs": str(ref.crs),
            "raster_crs": str(src.crs),
            "reference_width": ref.width,
            "raster_width": src.width,
            "reference_height": ref.height,
            "raster_height": src.height,
            "same_crs": same_crs,
            "same_width": same_width,
            "same_height": same_height,
            "same_transform": same_transform,
            "same_bounds": same_bounds,
            "grid_aligned": all(
                [
                    same_crs,
                    same_width,
                    same_height,
                    same_transform,
                    same_bounds,
                ]
            ),
        }


def transform_points_for_raster(
    longitude: np.ndarray,
    latitude: np.ndarray,
    raster_crs
) -> Tuple[np.ndarray, np.ndarray]:
    """Transform WGS84 occurrence coordinates to each raster CRS if necessary."""
    if raster_crs is None:
        stop_with_message("A raster has no CRS.")

    if raster_crs.to_epsg() == 4326:
        return longitude, latitude

    transformer = Transformer.from_crs(
        "EPSG:4326",
        raster_crs,
        always_xy=True
    )

    x, y = transformer.transform(longitude, latitude)

    return np.asarray(x), np.asarray(y)


def sample_raster_at_lonlat(
    raster_path: Path,
    longitude: np.ndarray,
    latitude: np.ndarray
) -> np.ndarray:
    """Extract values from one raster at WGS84 longitude/latitude coordinates."""
    with rasterio.open(raster_path) as src:
        x, y = transform_points_for_raster(
            longitude,
            latitude,
            src.crs
        )

        values = np.asarray(
            [sample[0] for sample in src.sample(zip(x, y))],
            dtype=float
        )

        if src.nodata is not None and np.isfinite(src.nodata):
            values[np.isclose(values, src.nodata)] = np.nan

        values[~np.isfinite(values)] = np.nan

    return values


def generate_global_screening_coordinates(
    reference_raster_path: Path,
    sample_size: int,
    seed: int
) -> pd.DataFrame:
    """
    Draw a reproducible random sample from valid cells of Bio1.

    This sample is used only for predictor correlation screening. It is not
    yet the final MaxEnt background sample.
    """
    rng = np.random.default_rng(seed)

    with rasterio.open(reference_raster_path) as src:
        reference = src.read(1, masked=True)
        valid = ~np.ma.getmaskarray(reference)

        data = np.asarray(reference.filled(np.nan), dtype=float)
        valid &= np.isfinite(data)

        valid_rows, valid_cols = np.where(valid)

        if len(valid_rows) == 0:
            stop_with_message(
                "No valid terrestrial cells were found in the Bio1 reference raster."
            )

        n_draw = min(sample_size, len(valid_rows))

        selected = rng.choice(
            len(valid_rows),
            size=n_draw,
            replace=False
        )

        rows = valid_rows[selected]
        cols = valid_cols[selected]

        xs, ys = rasterio.transform.xy(
            src.transform,
            rows,
            cols,
            offset="center"
        )

        xs = np.asarray(xs, dtype=float)
        ys = np.asarray(ys, dtype=float)

        if src.crs is None:
            stop_with_message(
                f"Reference raster has no CRS:\n{reference_raster_path}"
            )

        if src.crs.to_epsg() == 4326:
            longitude = xs
            latitude = ys
        else:
            transformer = Transformer.from_crs(
                src.crs,
                "EPSG:4326",
                always_xy=True
            )
            longitude, latitude = transformer.transform(xs, ys)
            longitude = np.asarray(longitude, dtype=float)
            latitude = np.asarray(latitude, dtype=float)

    return pd.DataFrame(
        {
            "screening_sample_id": [
                f"global_screen_{index:05d}"
                for index in range(1, len(longitude) + 1)
            ],
            "decimalLongitude": longitude,
            "decimalLatitude": latitude,
            "reference_raster_row": rows,
            "reference_raster_col": cols,
        }
    )


def calculate_missingness(df: pd.DataFrame, predictor_codes: List[str], dataset_label: str):
    rows = []

    for code in predictor_codes:
        values = pd.to_numeric(df[code], errors="coerce")
        n_total = len(values)
        n_valid = int(values.notna().sum())
        n_missing = int(values.isna().sum())

        rows.append(
            {
                "dataset": dataset_label,
                "predictor": code,
                "n_records": n_total,
                "n_valid": n_valid,
                "n_missing": n_missing,
                "missing_percent": 100 * n_missing / n_total,
                "minimum": values.min(skipna=True),
                "maximum": values.max(skipna=True),
                "mean": values.mean(skipna=True),
                "standard_deviation": values.std(skipna=True),
            }
        )

    return pd.DataFrame(rows)


def correlation_to_long_form(
    value_df: pd.DataFrame,
    predictor_codes: List[str],
    dataset_label: str
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Return correlation matrix and unique pair table.
    """
    matrix = value_df[predictor_codes].apply(
        pd.to_numeric,
        errors="coerce"
    ).corr(method="pearson")

    matrix_out = (
        matrix.reset_index()
        .rename(columns={"index": "predictor"})
    )

    pairs = []

    for i, predictor_a in enumerate(predictor_codes):
        for predictor_b in predictor_codes[i + 1:]:
            valid_pair = value_df[
                [predictor_a, predictor_b]
            ].dropna()

            if len(valid_pair) < 3:
                r_value = np.nan
            else:
                r_value = valid_pair[predictor_a].corr(
                    valid_pair[predictor_b],
                    method="pearson"
                )

            pairs.append(
                {
                    "dataset": dataset_label,
                    "predictor_1": predictor_a,
                    "predictor_2": predictor_b,
                    "pearson_r": r_value,
                    "abs_pearson_r": abs(r_value) if pd.notna(r_value) else np.nan,
                    "n_complete_pairs": len(valid_pair),
                    "abs_r_ge_threshold": (
                        abs(r_value) >= CORRELATION_THRESHOLD
                        if pd.notna(r_value)
                        else False
                    ),
                }
            )

    return matrix_out, pd.DataFrame(pairs)


# ============================================================================
# 4. MAIN
# ============================================================================

def main() -> None:
    if not OCCURRENCE_CSV.exists():
        stop_with_message(
            "Step 05 final occurrence file was not found:\n"
            f"{OCCURRENCE_CSV}\n\n"
            "Run 05_global_spatial_thinning_5km_and_gridcell_dedup.py first."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    rasters, excluded_report = locate_primary_rasters()

    reference_raster = rasters["Bio1"]

    print("\n" + "=" * 82)
    print("GLOBAL PRIMARY PREDICTOR EXTRACTION AND SCREENING: STEP 06")
    print(f"Occurrence file:\n{OCCURRENCE_CSV}")
    print(f"Reference raster:\n{reference_raster}")
    print(f"Primary predictors found: {len(rasters)}")
    print(f"Global screening sample size: {GLOBAL_SCREENING_SAMPLE_SIZE}")
    print(f"Fixed random seed: {RANDOM_SEED}")
    print("=" * 82)

    occurrences = pd.read_csv(OCCURRENCE_CSV, low_memory=False)

    required_occurrence_columns = {
        "record_id",
        "decimalLongitude",
        "decimalLatitude",
    }

    missing_occurrence_columns = required_occurrence_columns - set(
        occurrences.columns
    )

    if missing_occurrence_columns:
        stop_with_message(
            "Occurrence file lacks required columns:\n"
            + ", ".join(sorted(missing_occurrence_columns))
        )

    occurrences["decimalLongitude"] = pd.to_numeric(
        occurrences["decimalLongitude"],
        errors="coerce"
    )
    occurrences["decimalLatitude"] = pd.to_numeric(
        occurrences["decimalLatitude"],
        errors="coerce"
    )

    occurrences = occurrences[
        occurrences["decimalLongitude"].between(-180, 180)
        & occurrences["decimalLatitude"].between(-90, 90)
    ].copy()

    if occurrences.empty:
        stop_with_message("No valid occurrence coordinates remain.")

    occurrences = occurrences.sort_values("record_id").reset_index(drop=True)

    # ------------------------------------------------------------------------
    # 4.1 Raster grid QC
    # ------------------------------------------------------------------------

    grid_rows = []

    for code, raster_path in rasters.items():
        grid_row = compare_grid(reference_raster, raster_path)
        grid_row["predictor"] = code
        grid_rows.append(grid_row)

    grid_qc = pd.DataFrame(grid_rows).sort_values("predictor")

    if not grid_qc["grid_aligned"].all():
        misaligned = grid_qc.loc[
            ~grid_qc["grid_aligned"],
            ["predictor", "raster_file"]
        ]

        stop_with_message(
            "One or more primary rasters are not aligned with Bio1. "
            "Align all current predictor rasters before proceeding.\n\n"
            + misaligned.to_string(index=False)
        )

    # ------------------------------------------------------------------------
    # 4.2 Extract values at final occurrences
    # ------------------------------------------------------------------------

    occurrence_values = occurrences.copy()

    occurrence_longitude = occurrence_values["decimalLongitude"].to_numpy()
    occurrence_latitude = occurrence_values["decimalLatitude"].to_numpy()

    print("\nExtracting predictor values at final occurrence records...")

    for code, raster_path in rasters.items():
        print(f"  {code}")
        occurrence_values[code] = sample_raster_at_lonlat(
            raster_path,
            occurrence_longitude,
            occurrence_latitude
        )

    # ------------------------------------------------------------------------
    # 4.3 Draw global terrestrial environmental sample
    # ------------------------------------------------------------------------

    print("\nDrawing global terrestrial environmental screening sample...")

    global_screen = generate_global_screening_coordinates(
        reference_raster_path=reference_raster,
        sample_size=GLOBAL_SCREENING_SAMPLE_SIZE,
        seed=RANDOM_SEED
    )

    global_longitude = global_screen["decimalLongitude"].to_numpy()
    global_latitude = global_screen["decimalLatitude"].to_numpy()

    print("Extracting predictor values at global screening cells...")

    for code, raster_path in rasters.items():
        print(f"  {code}")
        global_screen[code] = sample_raster_at_lonlat(
            raster_path,
            global_longitude,
            global_latitude
        )

    # Keep only sample locations with complete primary predictor values for
    # the global correlation analysis. The original raw sample is retained too.
    global_screen_complete = global_screen.dropna(
        subset=list(rasters.keys())
    ).copy()

    # ------------------------------------------------------------------------
    # 4.4 Missingness, correlations, and reports
    # ------------------------------------------------------------------------

    predictor_codes = list(rasters.keys())

    missingness = pd.concat(
        [
            calculate_missingness(
                occurrence_values,
                predictor_codes,
                "final_occurrence_records"
            ),
            calculate_missingness(
                global_screen,
                predictor_codes,
                "global_environmental_screening_sample"
            ),
        ],
        ignore_index=True
    )

    occurrence_matrix, occurrence_pairs = correlation_to_long_form(
        occurrence_values,
        predictor_codes,
        "final_occurrence_records"
    )

    global_matrix, global_pairs = correlation_to_long_form(
        global_screen_complete,
        predictor_codes,
        "global_environmental_screening_sample"
    )

    all_pairs = pd.concat(
        [occurrence_pairs, global_pairs],
        ignore_index=True
    )

    high_pairs = all_pairs.loc[
        all_pairs["abs_r_ge_threshold"]
    ].sort_values(
        ["dataset", "abs_pearson_r"],
        ascending=[True, False]
    )

    # ------------------------------------------------------------------------
    # 4.5 Exports
    # ------------------------------------------------------------------------

    raster_inventory = pd.DataFrame(
        [
            {
                "predictor": code,
                "raster_file": str(path),
                "status": "primary_candidate_for_screening",
                "reason": (
                    "Current and matched future projection layers expected; "
                    "included for correlation and VIF screening"
                ),
            }
            for code, path in rasters.items()
        ]
    )

    raster_inventory.to_csv(
        OUTPUT_DIR / "06_primary_candidate_raster_inventory.csv",
        index=False,
        encoding="utf-8-sig"
    )

    excluded_report.to_csv(
        OUTPUT_DIR / "06_excluded_static_or_unmatched_rasters.csv",
        index=False,
        encoding="utf-8-sig"
    )

    grid_qc.to_csv(
        OUTPUT_DIR / "06_primary_predictor_grid_alignment_QC.csv",
        index=False,
        encoding="utf-8-sig"
    )

    occurrence_values.to_csv(
        OUTPUT_DIR / "06_occurrence_primary_predictor_values.csv",
        index=False,
        encoding="utf-8-sig"
    )

    global_screen.to_csv(
        OUTPUT_DIR / "06_global_environmental_screening_values_raw.csv",
        index=False,
        encoding="utf-8-sig"
    )

    global_screen_complete.to_csv(
        OUTPUT_DIR / "06_global_environmental_screening_values_complete.csv",
        index=False,
        encoding="utf-8-sig"
    )

    missingness.to_csv(
        OUTPUT_DIR / "06_predictor_missingness_and_range_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    occurrence_matrix.to_csv(
        OUTPUT_DIR / "06_pearson_correlation_matrix_occurrences.csv",
        index=False,
        encoding="utf-8-sig"
    )

    global_matrix.to_csv(
        OUTPUT_DIR / "06_pearson_correlation_matrix_global_sample.csv",
        index=False,
        encoding="utf-8-sig"
    )

    all_pairs.to_csv(
        OUTPUT_DIR / "06_all_pearson_correlation_pairs.csv",
        index=False,
        encoding="utf-8-sig"
    )

    high_pairs.to_csv(
        OUTPUT_DIR / "06_high_correlation_pairs_abs_ge_080.csv",
        index=False,
        encoding="utf-8-sig"
    )

    settings = pd.DataFrame(
        {
            "parameter": [
                "final_occurrence_dataset",
                "n_final_occurrences",
                "primary_predictors",
                "primary_predictor_count",
                "excluded_predictor_policy",
                "global_screening_sample_size_requested",
                "global_screening_sample_size_complete",
                "global_screening_sample_seed",
                "pearson_correlation_threshold",
                "reference_raster",
            ],
            "value": [
                "Strict global natural-land-cover candidates after 5 km thinning",
                len(occurrence_values),
                "; ".join(predictor_codes),
                len(predictor_codes),
                (
                    "Soil moisture, NPP, runoff, soil carbon, soil pH, and UV "
                    "layers were excluded from the primary future-transfer "
                    "model because scenario-matched future layers are absent."
                ),
                GLOBAL_SCREENING_SAMPLE_SIZE,
                len(global_screen_complete),
                RANDOM_SEED,
                CORRELATION_THRESHOLD,
                str(reference_raster),
            ],
        }
    )

    settings.to_csv(
        OUTPUT_DIR / "06_predictor_screening_settings.csv",
        index=False,
        encoding="utf-8-sig"
    )

    note = [
        "This step screens only bioclimatic variables and elevation for the primary future-transfer model.",
        "Static soil moisture, NPP, runoff, soil properties, and UV layers were excluded because they lack matched future scenario projections in this workflow.",
        "Final predictor selection must retain ecologically interpretable variables, avoid strong pairwise correlations, and be verified with VIF in the next step.",
        "The global environmental screening sample is for correlation assessment only and is not yet the final MaxEnt background sample.",
    ]

    (OUTPUT_DIR / "06_predictor_screening_interpretation_note.txt").write_text(
        "\n".join(note),
        encoding="utf-8"
    )

    # ------------------------------------------------------------------------
    # 4.6 Console report
    # ------------------------------------------------------------------------

    print("\n" + "=" * 82)
    print("Primary predictor extraction and screening completed successfully")
    print(f"Final occurrence records: {len(occurrence_values)}")
    print(
        "Global screening cells with complete values: "
        f"{len(global_screen_complete)}"
    )
    print(f"Primary predictor candidates: {len(predictor_codes)}")
    print(
        "High-correlation pairs (|r| >= "
        f"{CORRELATION_THRESHOLD:.2f}): {len(high_pairs)}"
    )
    print("Output folder:")
    print(OUTPUT_DIR)
    print("=" * 82 + "\n")


if __name__ == "__main__":
    main()
