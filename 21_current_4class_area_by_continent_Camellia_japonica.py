# -*- coding: utf-8 -*-
"""
Calculate area of 4 current habitat-suitability classes by continent.

Input:
1) Continent boundary shapefile
2) Four-class suitability raster

Output:
- Long-format CSV
- Wide-format CSV
- Manuscript-ready CSV with areas in million km2 and percentages

Run in VS Code:
    python 21_current_4class_area_by_continent_Camellia_japonica.py

Required packages:
    pip install geopandas rasterio pyproj pandas numpy shapely
"""

import os
import warnings
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import geometry_mask
from pyproj import Geod

warnings.filterwarnings("ignore", category=UserWarning)

# ============================================================
# 1. INPUT PATHS
# ============================================================

CONTINENT_SHP = r"D:\Prediction_tree\continents\World_Continents (2).shp"

RASTER_4CLASS = (
    r"E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs"
    r"\10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
    r"\05_current_suitability_maps"
    r"\Camellia_japonica_current_final_full_data_MTSPS_4classes.tif"
)

OUT_DIR = (
    r"E:\PREDICTION TREE\Camellia_japonica\global_reanalysis_outputs"
    r"\10_final_MaxEnt_RM1_lph_CV100_MTSPS_4classes_FIXED_v2"
    r"\07_current_area_by_continent"
)

os.makedirs(OUT_DIR, exist_ok=True)

# If your shapefile has a specific continent-name field, set it here.
# Example: CONTINENT_FIELD = "CONTINENT"
# Leave as None to auto-detect.
CONTINENT_FIELD = None

# Raster class definitions
CLASS_LABELS = {
    1: "Unsuitable",
    2: "Low suitability",
    3: "Moderate suitability",
    4: "High suitability",
}

SUITABLE_CLASSES = [2, 3, 4]


# ============================================================
# 2. HELPER FUNCTIONS
# ============================================================

def detect_continent_field(gdf: gpd.GeoDataFrame) -> str:
    """Auto-detect field containing continent names."""
    candidate_fields = [
        "CONTINENT", "continent", "Continent",
        "NAME", "Name", "name",
        "REGION", "Region", "region",
        "CONT_NAME", "CONTINENT_NAME", "continent_name",
    ]

    for field in candidate_fields:
        if field in gdf.columns:
            return field

    non_geom_cols = [c for c in gdf.columns if c != gdf.geometry.name]
    raise ValueError(
        "Could not auto-detect continent field. Available fields are:\n"
        f"{non_geom_cols}\n"
        "Please set CONTINENT_FIELD manually in the script."
    )


def repair_geometries(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Repair invalid geometries using buffer(0)."""
    gdf = gdf.copy()
    gdf = gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()
    if not gdf.is_valid.all():
        gdf["geometry"] = gdf.geometry.buffer(0)
    return gdf[~gdf.geometry.is_empty & gdf.geometry.notna()].copy()


def calculate_geodesic_cell_area_by_row(transform, height, crs) -> np.ndarray:
    """
    Return cell area by raster row in km2.

    For geographic CRS, cell area is calculated geodesically with WGS84 ellipsoid.
    This avoids the reviewer issue that WGS84 cells vary in size by latitude.
    """
    if crs is None:
        raise ValueError("Raster has no CRS. Please define CRS before calculating area.")

    if abs(transform.b) > 1e-12 or abs(transform.d) > 1e-12:
        raise ValueError("Rotated rasters are not supported by this script.")

    if crs.is_geographic:
        geod = Geod(ellps="WGS84")
        cell_areas_km2 = np.zeros(height, dtype="float64")

        x_left = transform.c
        x_right = transform.c + transform.a

        for row in range(height):
            y_top = transform.f + row * transform.e
            y_bottom = transform.f + (row + 1) * transform.e

            lons = [x_left, x_right, x_right, x_left]
            lats = [y_top, y_top, y_bottom, y_bottom]

            area_m2, _ = geod.polygon_area_perimeter(lons, lats)
            cell_areas_km2[row] = abs(area_m2) / 1_000_000.0

        return cell_areas_km2

    # Projected raster: assumes metre units.
    pixel_area_m2 = abs(transform.a * transform.e)
    return np.full(height, pixel_area_m2 / 1_000_000.0, dtype="float64")


def clean_class_array(data: np.ndarray, nodata) -> np.ndarray:
    """Convert class raster to integer 1–4 and set other cells to 0."""
    arr = data.copy()

    if nodata is not None:
        arr = np.where(arr == nodata, 0, arr)

    arr = np.rint(arr).astype("int16")

    valid_mask = np.isin(arr, list(CLASS_LABELS.keys()))
    arr = np.where(valid_mask, arr, 0).astype("int16")

    return arr


# ============================================================
# 3. MAIN PROCESS
# ============================================================

def main():
    print("Reading raster...")
    with rasterio.open(RASTER_4CLASS) as src:
        raster_crs = src.crs
        transform = src.transform
        height = src.height
        width = src.width
        nodata = src.nodata
        raster_shape = (height, width)
        class_arr = clean_class_array(src.read(1), nodata)

    print(f"Raster CRS: {raster_crs}")
    print(f"Raster size: {width} columns x {height} rows")
    print(f"Raster nodata: {nodata}")
    print(f"Unique class values after cleaning: {np.unique(class_arr)}")

    print("\nCalculating geodesic cell areas by row...")
    cell_area_by_row_km2 = calculate_geodesic_cell_area_by_row(
        transform=transform,
        height=height,
        crs=raster_crs
    )

    print("Reading continent boundary...")
    continents = gpd.read_file(CONTINENT_SHP)
    continents = repair_geometries(continents)

    continent_field = CONTINENT_FIELD or detect_continent_field(continents)
    print(f"Using continent field: {continent_field}")

    continents = continents[[continent_field, "geometry"]].copy()
    continents[continent_field] = continents[continent_field].astype(str).str.strip()

    # Remove Antarctica by default because it is not meaningful for C. japonica modelling.
    continents = continents[
        (continents[continent_field].notna()) &
        (continents[continent_field] != "") &
        (~continents[continent_field].str.lower().isin(["antarctica", "antarctique"]))
    ].copy()

    if continents.crs is None:
        raise ValueError("Continent shapefile has no CRS. Please define its CRS before running.")

    if continents.crs != raster_crs:
        print("Reprojecting continent polygons to raster CRS...")
        continents = continents.to_crs(raster_crs)

    continents = continents.dissolve(by=continent_field, as_index=False)
    continents = repair_geometries(continents)

    print("\nContinents to process:")
    for cont in continents[continent_field].tolist():
        print(f" - {cont}")

    results = []
    valid_raster = class_arr > 0

    for _, row in continents.iterrows():
        continent_name = row[continent_field]
        geom = row.geometry

        print(f"\nProcessing continent: {continent_name}")

        continent_mask = geometry_mask(
            geometries=[geom],
            out_shape=raster_shape,
            transform=transform,
            invert=True,
            all_touched=False
        )

        continent_valid = continent_mask & valid_raster
        class_areas = {}

        for class_code, class_name in CLASS_LABELS.items():
            class_mask = continent_valid & (class_arr == class_code)
            pixel_count_by_row = np.count_nonzero(class_mask, axis=1)
            area_km2 = float(np.sum(pixel_count_by_row * cell_area_by_row_km2))
            class_areas[class_code] = area_km2

            results.append({
                "continent": continent_name,
                "class_code": class_code,
                "class_name": class_name,
                "area_km2": area_km2,
                "area_million_km2": area_km2 / 1_000_000.0,
            })

        total_classified_area = sum(class_areas.values())
        total_suitable_area = sum(class_areas[c] for c in SUITABLE_CLASSES)

        results.append({
            "continent": continent_name,
            "class_code": 99,
            "class_name": "Total suitable",
            "area_km2": total_suitable_area,
            "area_million_km2": total_suitable_area / 1_000_000.0,
        })

        results.append({
            "continent": continent_name,
            "class_code": 100,
            "class_name": "Total classified area",
            "area_km2": total_classified_area,
            "area_million_km2": total_classified_area / 1_000_000.0,
        })

        print(f"  Total classified area: {total_classified_area:,.2f} km2")
        print(f"  Total suitable area:   {total_suitable_area:,.2f} km2")

    long_df = pd.DataFrame(results)

    total_area_lookup = (
        long_df[long_df["class_code"] == 100]
        .set_index("continent")["area_km2"]
        .to_dict()
    )

    long_df["continent_total_area_km2"] = long_df["continent"].map(total_area_lookup)
    long_df["percent_of_continent"] = np.where(
        long_df["continent_total_area_km2"] > 0,
        long_df["area_km2"] / long_df["continent_total_area_km2"] * 100.0,
        np.nan
    )

    long_csv = os.path.join(OUT_DIR, "current_4class_area_by_continent_long.csv")
    long_df.to_csv(long_csv, index=False, encoding="utf-8-sig")

    wide_km2 = long_df.pivot_table(
        index="continent",
        columns="class_name",
        values="area_km2",
        aggfunc="sum"
    ).reset_index()

    wide_pct = long_df.pivot_table(
        index="continent",
        columns="class_name",
        values="percent_of_continent",
        aggfunc="sum"
    ).reset_index()

    ordered_cols = [
        "continent",
        "Unsuitable",
        "Low suitability",
        "Moderate suitability",
        "High suitability",
        "Total suitable",
        "Total classified area",
    ]

    wide_km2 = wide_km2[[c for c in ordered_cols if c in wide_km2.columns]]
    wide_pct = wide_pct[[c for c in ordered_cols if c in wide_pct.columns]]

    wide_km2_csv = os.path.join(OUT_DIR, "current_4class_area_by_continent_wide_km2.csv")
    wide_pct_csv = os.path.join(OUT_DIR, "current_4class_area_by_continent_wide_percent.csv")

    wide_km2.to_csv(wide_km2_csv, index=False, encoding="utf-8-sig")
    wide_pct.to_csv(wide_pct_csv, index=False, encoding="utf-8-sig")

    manuscript_rows = []

    for cont in wide_km2["continent"]:
        row_km2 = wide_km2[wide_km2["continent"] == cont].iloc[0]
        row_pct = wide_pct[wide_pct["continent"] == cont].iloc[0]

        out_row = {"Continent": cont}
        for col in ordered_cols[1:]:
            if col in wide_km2.columns and col in wide_pct.columns:
                million_km2 = row_km2[col] / 1_000_000.0
                pct = row_pct[col]
                out_row[col] = f"{million_km2:.3f} ({pct:.2f}%)"
        manuscript_rows.append(out_row)

    manuscript_df = pd.DataFrame(manuscript_rows)
    manuscript_csv = os.path.join(
        OUT_DIR,
        "Table_current_suitability_area_by_continent_manuscript_ready.csv"
    )
    manuscript_df.to_csv(manuscript_csv, index=False, encoding="utf-8-sig")

    global_summary = (
        long_df[long_df["class_code"].isin(list(CLASS_LABELS.keys()) + [99, 100])]
        .groupby(["class_code", "class_name"], as_index=False)["area_km2"]
        .sum()
    )

    global_total = float(
        global_summary.loc[global_summary["class_code"] == 100, "area_km2"].iloc[0]
    )

    global_summary["area_million_km2"] = global_summary["area_km2"] / 1_000_000.0
    global_summary["percent_of_global_classified_area"] = (
        global_summary["area_km2"] / global_total * 100.0
    )

    global_csv = os.path.join(OUT_DIR, "current_4class_global_summary_from_continents.csv")
    global_summary.to_csv(global_csv, index=False, encoding="utf-8-sig")

    print("\nFinished.")
    print("Output files:")
    print(f"  1) {long_csv}")
    print(f"  2) {wide_km2_csv}")
    print(f"  3) {wide_pct_csv}")
    print(f"  4) {manuscript_csv}")
    print(f"  5) {global_csv}")

    print("\nManuscript-ready preview:")
    print(manuscript_df.to_string(index=False))


if __name__ == "__main__":
    main()
