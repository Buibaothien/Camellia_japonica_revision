"""
01_global_occurrence_QC_v2_auto_find_input.py

Global taxonomic, coordinate, GBIF-issue, coordinate-uncertainty, and
cultivation/ex-situ metadata screening for Camellia japonica.

This script DOES NOT:
- restrict records to a native-range polygon;
- use accessible area M;
- perform LULC filtering yet;
- perform spatial thinning yet.

It creates a global candidate table that will be passed to:
1) CoordinateCleaner,
2) ESA WorldCover LULC screening,
3) 5 km spatial thinning.

Key improvement over the previous version:
The raw GBIF file is found automatically anywhere inside the project folder.
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd


# =====================================================================
# 1. USER SETTINGS
# =====================================================================

PROJECT_DIR = Path(r"E:\PREDICTION TREE\Camellia_japonica")

# Leave as None to automatically locate the occurrence file.
# If automatic search finds more than one plausible file, set the exact path here.
RAW_OCCURRENCE_FILE = Path(r"E:\PREDICTION TREE\Camellia_japonica\data\0074456-251120083545085.csv")

# The script first looks for this exact historical file name, then searches
# CSV / TSV / TXT files with GBIF-style columns across the project folder.
PREFERRED_RAW_FILENAME = "0074456-251120083545085.csv"

OUT_DIR = PROJECT_DIR / "global_reanalysis_outputs" / "01_occurrence_QC"

TARGET_SPECIES = "camellia japonica"
MAX_REPORTED_COORD_UNCERTAINTY_M = 5000.0


# =====================================================================
# 2. FILTERING RULES
# =====================================================================

# Strong textual evidence of cultivation, managed growing, or ex-situ context.
# "park" is deliberately not included by itself because it can refer to natural areas.
CULTIVATION_PATTERNS = [
    r"\bbotani[ck]al garden\b",
    r"\barboretum\b",
    r"\bnursery\b",
    r"\bgreenhouse\b",
    r"\bcultivat(?:ed|ion|ar)\b",
    r"\bplanted\b",
    r"\bornamental\b",
    r"\bgarden\b",
    r"\byard\b",
    r"\bresidential\b",
    r"\bcampus\b",
    r"\bprivate collection\b",
    r"\bex[\s-]?situ\b",
    r"\blandscape plant(?:ing)?\b",
    r"\burban park\b",
]

# GBIF issue strings considered fatal for spatial modelling.
FATAL_ISSUE_TOKENS = {
    "COORDINATE_INVALID",
    "COORDINATE_OUT_OF_RANGE",
    "ZERO_COORDINATE",
    "COUNTRY_COORDINATE_MISMATCH",
    "CONTINENT_COORDINATE_MISMATCH",
    "PRESUMED_SWAPPED_COORDINATES",
    "COORDINATE_REPROJECTION_SUSPICIOUS",
    "GEODETIC_DATUM_INVALID",
    "TAXON_MATCH_NONE",
}


# =====================================================================
# 3. HELPER FUNCTIONS
# =====================================================================

def normalise_column_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def find_column(df: pd.DataFrame, candidates, required=True):
    normalised_columns = {
        normalise_column_name(column): column
        for column in df.columns
    }

    for candidate in candidates:
        candidate_key = normalise_column_name(candidate)
        if candidate_key in normalised_columns:
            return normalised_columns[candidate_key]

    if required:
        raise KeyError(
            f"Missing required column. Tried: {candidates}\n"
            f"Available columns: {list(df.columns)}"
        )

    return None


def read_occurrence_table(path: Path) -> pd.DataFrame:
    """Read CSV/TSV occurrence data, including GBIF downloads saved as .csv."""
    return pd.read_csv(
        path,
        sep=None,
        engine="python",
        dtype=str,
        encoding_errors="replace",
    )


def file_has_occurrence_columns(path: Path) -> bool:
    """Check whether a candidate table looks like an occurrence dataset."""
    try:
        preview = pd.read_csv(
            path,
            sep=None,
            engine="python",
            dtype=str,
            encoding_errors="replace",
            nrows=5,
        )
    except Exception:
        return False

    columns = {normalise_column_name(column) for column in preview.columns}

    has_lon = any(name in columns for name in ["decimallongitude", "longitude", "lon"])
    has_lat = any(name in columns for name in ["decimallatitude", "latitude", "lat"])
    has_taxon = any(
        name in columns
        for name in ["species", "scientificname", "taxonkey", "acceptedtaxonkey"]
    )

    return has_lon and has_lat and has_taxon


def locate_raw_occurrence_file(project_dir: Path) -> Path:
    """
    Find the original occurrence file automatically.

    Search order:
    1. User-specified RAW_OCCURRENCE_FILE, if defined.
    2. Exact preferred file name in the project root.
    3. Exact preferred file name anywhere under the project folder.
    4. Any CSV/TSV/TXT under the project folder containing occurrence columns.
    """
    if RAW_OCCURRENCE_FILE is not None:
        specified_path = Path(RAW_OCCURRENCE_FILE)

        if not specified_path.is_absolute():
            specified_path = project_dir / specified_path

        if not specified_path.exists():
            raise FileNotFoundError(
                "RAW_OCCURRENCE_FILE was manually specified but was not found:\n"
                f"{specified_path}"
            )

        return specified_path

    root_candidate = project_dir / PREFERRED_RAW_FILENAME
    if root_candidate.exists():
        return root_candidate

    exact_name_candidates = sorted(
        path
        for path in project_dir.rglob(PREFERRED_RAW_FILENAME)
        if "global_reanalysis_outputs" not in path.parts
    )

    if len(exact_name_candidates) == 1:
        return exact_name_candidates[0]

    if len(exact_name_candidates) > 1:
        candidate_text = "\n".join(str(path) for path in exact_name_candidates)

        raise FileNotFoundError(
            "More than one file with the expected GBIF name was found.\n"
            "Set RAW_OCCURRENCE_FILE at the top of this script to the correct file.\n\n"
            f"Candidates:\n{candidate_text}"
        )

    allowed_extensions = {".csv", ".tsv", ".txt"}

    broad_candidates = [
        path
        for path in project_dir.rglob("*")
        if path.is_file()
        and path.suffix.lower() in allowed_extensions
        and "global_reanalysis_outputs" not in path.parts
    ]

    valid_occurrence_candidates = [
        path for path in broad_candidates if file_has_occurrence_columns(path)
    ]

    if len(valid_occurrence_candidates) == 1:
        return valid_occurrence_candidates[0]

    if len(valid_occurrence_candidates) == 0:
        raise FileNotFoundError(
            "No suitable occurrence CSV/TSV/TXT file could be found anywhere under:\n"
            f"{project_dir}\n\n"
            "Place the original GBIF occurrence file inside this project folder, "
            "or set RAW_OCCURRENCE_FILE to its exact full path at the top of the script."
        )

    candidate_text = "\n".join(str(path) for path in valid_occurrence_candidates)

    raise FileNotFoundError(
        "More than one plausible occurrence table was found.\n"
        "Set RAW_OCCURRENCE_FILE at the top of this script to the correct file.\n\n"
        f"Candidates:\n{candidate_text}"
    )


def as_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def norm_text(series: pd.Series) -> pd.Series:
    return (
        series.fillna("")
        .astype(str)
        .str.lower()
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )


def contains_fatal_issue(value) -> bool:
    if pd.isna(value):
        return False

    issue_tokens = re.split(r"[;,|]", str(value).upper())

    return any(
        token.strip() in FATAL_ISSUE_TOKENS
        for token in issue_tokens
    )


def first_present(df: pd.DataFrame, columns) -> pd.Series:
    usable_columns = [column for column in columns if column in df.columns]

    if not usable_columns:
        return pd.Series("", index=df.index, dtype="object")

    return df[usable_columns].fillna("").astype(str).agg(" | ".join, axis=1)


# =====================================================================
# 4. MAIN WORKFLOW
# =====================================================================

def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    raw_file = locate_raw_occurrence_file(PROJECT_DIR)

    print("\n" + "=" * 72)
    print("GLOBAL OCCURRENCE QC: STEP 01")
    print(f"Input file detected:\n{raw_file}")
    print("=" * 72)

    raw = read_occurrence_table(raw_file)
    raw.columns = [str(column).strip() for column in raw.columns]
    raw["_source_row"] = np.arange(1, len(raw) + 1)

    lon_col = find_column(raw, ["decimalLongitude", "longitude", "lon"])
    lat_col = find_column(raw, ["decimalLatitude", "latitude", "lat"])

    species_col = find_column(raw, ["species"], required=False)
    scientific_name_col = find_column(raw, ["scientificName"], required=False)
    taxon_key_col = find_column(raw, ["taxonKey", "acceptedTaxonKey"], required=False)
    source_id_col = find_column(raw, ["gbifID", "key", "occurrenceID"], required=False)

    uncertainty_col = find_column(
        raw,
        ["coordinateUncertaintyInMeters"],
        required=False
    )

    issues_col = find_column(raw, ["issues"], required=False)
    geospatial_issues_col = find_column(raw, ["hasGeospatialIssues"], required=False)
    establishment_col = find_column(
        raw,
        ["establishmentMeans", "degreeOfEstablishment"],
        required=False
    )
    basis_col = find_column(raw, ["basisOfRecord"], required=False)
    year_col = find_column(raw, ["year"], required=False)
    event_date_col = find_column(raw, ["eventDate"], required=False)

    df = raw.copy()

    # ----------------------------------------------------------
    # 4.1 Coordinates and taxonomic identity
    # ----------------------------------------------------------

    df["decimalLongitude"] = as_num(df[lon_col])
    df["decimalLatitude"] = as_num(df[lat_col])

    if species_col is not None:
        species_name = norm_text(df[species_col])
    elif scientific_name_col is not None:
        species_name = norm_text(df[scientific_name_col])
    else:
        raise KeyError(
            "Neither species nor scientificName is available for taxonomic filtering."
        )

    exact_species = species_name.eq(TARGET_SPECIES)

    scientific_name_match = pd.Series(False, index=df.index)

    if scientific_name_col is not None:
        scientific_name_normalised = norm_text(df[scientific_name_col])

        scientific_name_match = scientific_name_normalised.str.match(
            r"^camellia japonica(\b|$)",
            na=False
        )

        # Remove hybrid expressions such as Camellia japonica × ...
        scientific_name_match &= ~scientific_name_normalised.str.contains(
            r"\bx\b|×",
            regex=True,
            na=False
        )

    df["flag_taxon_target"] = exact_species | scientific_name_match

    df["flag_valid_coordinates"] = (
        df["decimalLongitude"].between(-180, 180)
        & df["decimalLatitude"].between(-90, 90)
        & ~(df["decimalLongitude"].eq(0) & df["decimalLatitude"].eq(0))
    )

    # ----------------------------------------------------------
    # 4.2 Coordinate uncertainty and GBIF issue fields
    # ----------------------------------------------------------

    if uncertainty_col is not None:
        df["coordinateUncertaintyInMeters_numeric"] = as_num(df[uncertainty_col])
    else:
        df["coordinateUncertaintyInMeters_numeric"] = np.nan

    # Missing coordinate uncertainty is retained but will remain transparent
    # in the audit output. Reported uncertainty >5 km is excluded.
    df["flag_reported_uncertainty_acceptable"] = (
        df["coordinateUncertaintyInMeters_numeric"].isna()
        | (
            df["coordinateUncertaintyInMeters_numeric"]
            <= MAX_REPORTED_COORD_UNCERTAINTY_M
        )
    )

    df["uncertainty_status"] = np.where(
        df["coordinateUncertaintyInMeters_numeric"].isna(),
        "not_reported",
        np.where(
            df["coordinateUncertaintyInMeters_numeric"]
            <= MAX_REPORTED_COORD_UNCERTAINTY_M,
            "reported_le_5km",
            "reported_gt_5km",
        ),
    )

    if issues_col is not None:
        df["flag_no_fatal_gbif_issue"] = ~df[issues_col].apply(
            contains_fatal_issue
        )
    else:
        df["flag_no_fatal_gbif_issue"] = True

    if geospatial_issues_col is not None:
        declared_issue = norm_text(df[geospatial_issues_col])

        df["flag_no_declared_geospatial_issue"] = ~declared_issue.isin(
            ["true", "1", "yes"]
        )
    else:
        df["flag_no_declared_geospatial_issue"] = True

    # ----------------------------------------------------------
    # 4.3 Cultivation and ex-situ metadata screening
    # ----------------------------------------------------------

    text_fields = [
        column
        for column in [
            find_column(df, ["locality"], required=False),
            find_column(df, ["verbatimLocality"], required=False),
            find_column(df, ["occurrenceRemarks"], required=False),
            find_column(df, ["habitat"], required=False),
            find_column(df, ["fieldNotes"], required=False),
            find_column(df, ["recordedBy"], required=False),
            find_column(df, ["institutionCode"], required=False),
            find_column(df, ["collectionCode"], required=False),
        ]
        if column is not None
    ]

    combined_text = norm_text(first_present(df, text_fields))
    cultivation_regex = "|".join(CULTIVATION_PATTERNS)

    df["flag_no_cultivation_text"] = ~combined_text.str.contains(
        cultivation_regex,
        regex=True,
        na=False
    )

    df["cultivation_text_flag"] = np.where(
        df["flag_no_cultivation_text"],
        "",
        "cultivation_or_exsitu_keyword",
    )

    if establishment_col is not None:
        establishment_text = norm_text(df[establishment_col])

        df["flag_no_cultivation_establishment"] = ~establishment_text.str.contains(
            r"cultivat|planted|introduced|managed|captive|ornamental",
            regex=True,
            na=False,
        )
    else:
        df["flag_no_cultivation_establishment"] = True

    if basis_col is not None:
        basis_text = norm_text(df[basis_col])

        # Living specimens and material samples are commonly ex-situ records.
        df["flag_acceptable_basis"] = ~basis_text.isin(
            ["living_specimen", "material_sample"]
        )
    else:
        df["flag_acceptable_basis"] = True

    # ----------------------------------------------------------
    # 4.4 IDs, year and duplicate coordinates
    # ----------------------------------------------------------

    if source_id_col is not None:
        df["source_record_id"] = (
            df[source_id_col].fillna("").astype(str).str.strip()
        )
    else:
        df["source_record_id"] = ""

    df["record_id"] = np.where(
        df["source_record_id"].ne(""),
        "source_" + df["source_record_id"],
        "row_" + df["_source_row"].astype(str),
    )

    # Guarantee a unique identifier even when source IDs are duplicated.
    df["record_id"] = (
        df["record_id"]
        + "_"
        + df.groupby("record_id").cumcount().astype(str)
    )

    if year_col is not None:
        df["record_year"] = as_num(df[year_col])
    elif event_date_col is not None:
        df["record_year"] = pd.to_datetime(
            df[event_date_col],
            errors="coerce"
        ).dt.year
    else:
        df["record_year"] = np.nan

    df["coord_key_6dp"] = (
        df["decimalLongitude"].round(6).astype(str)
        + "_"
        + df["decimalLatitude"].round(6).astype(str)
    )

    # At exact duplicate coordinates, retain the record with the lowest
    # reported uncertainty. Original row order breaks remaining ties.
    df["_uncertainty_sort"] = (
        df["coordinateUncertaintyInMeters_numeric"].fillna(np.inf)
    )

    df = df.sort_values(
        ["_uncertainty_sort", "_source_row"],
        kind="mergesort"
    )

    df["flag_unique_exact_coordinate"] = ~df.duplicated(
        "coord_key_6dp",
        keep="first"
    )

    # ----------------------------------------------------------
    # 4.5 Create strict metadata-QC candidate set
    # ----------------------------------------------------------

    strict_keep = (
        df["flag_taxon_target"]
        & df["flag_valid_coordinates"]
        & df["flag_reported_uncertainty_acceptable"]
        & df["flag_no_fatal_gbif_issue"]
        & df["flag_no_declared_geospatial_issue"]
        & df["flag_no_cultivation_text"]
        & df["flag_no_cultivation_establishment"]
        & df["flag_acceptable_basis"]
        & df["flag_unique_exact_coordinate"]
    )

    df["metadata_QC_keep"] = strict_keep

    # ----------------------------------------------------------
    # 4.6 Write outputs
    # ----------------------------------------------------------

    keep_columns = [
        "record_id",
        "source_record_id",
        "decimalLongitude",
        "decimalLatitude",
        "coordinateUncertaintyInMeters_numeric",
        "uncertainty_status",
        "record_year",
        "coord_key_6dp",
        "metadata_QC_keep",
        "cultivation_text_flag",
        "flag_taxon_target",
        "flag_valid_coordinates",
        "flag_reported_uncertainty_acceptable",
        "flag_no_fatal_gbif_issue",
        "flag_no_declared_geospatial_issue",
        "flag_no_cultivation_text",
        "flag_no_cultivation_establishment",
        "flag_acceptable_basis",
        "flag_unique_exact_coordinate",
    ]

    # Preserve selected raw metadata columns for traceability.
    for column in (
        [species_col, scientific_name_col, taxon_key_col, issues_col,
         establishment_col, basis_col]
        + text_fields
    ):
        if column is not None and column not in keep_columns:
            keep_columns.append(column)

    keep_columns = [
        column for column in keep_columns if column in df.columns
    ]

    candidates = (
        df.loc[df["metadata_QC_keep"], keep_columns]
        .copy()
        .sort_values("record_id")
        .reset_index(drop=True)
    )

    # Buffer radius will be used later for ESA WorldCover screening.
    # It respects reported precision but avoids sampling a single 10 m pixel.
    candidates["lulc_buffer_m"] = np.where(
        candidates["coordinateUncertaintyInMeters_numeric"].notna(),
        candidates["coordinateUncertaintyInMeters_numeric"].clip(
            lower=100,
            upper=1000
        ),
        250.0,
    )

    audit = pd.DataFrame(
        {
            "stage": [
                "raw_input",
                "target_taxon",
                "valid_coordinates",
                "reported_uncertainty_le_5km_or_missing",
                "no_fatal_GBIF_issue",
                "no_declared_geospatial_issue",
                "no_cultivation_or_exsitu_metadata",
                "acceptable_basisOfRecord",
                "unique_exact_coordinate",
                "metadata_QC_candidates",
            ],
            "n_records": [
                len(raw),
                int(df["flag_taxon_target"].sum()),
                int(
                    (
                        df["flag_taxon_target"]
                        & df["flag_valid_coordinates"]
                    ).sum()
                ),
                int(
                    (
                        df["flag_taxon_target"]
                        & df["flag_valid_coordinates"]
                        & df["flag_reported_uncertainty_acceptable"]
                    ).sum()
                ),
                int(
                    (
                        df["flag_taxon_target"]
                        & df["flag_valid_coordinates"]
                        & df["flag_reported_uncertainty_acceptable"]
                        & df["flag_no_fatal_gbif_issue"]
                    ).sum()
                ),
                int(
                    (
                        df["flag_taxon_target"]
                        & df["flag_valid_coordinates"]
                        & df["flag_reported_uncertainty_acceptable"]
                        & df["flag_no_fatal_gbif_issue"]
                        & df["flag_no_declared_geospatial_issue"]
                    ).sum()
                ),
                int(
                    (
                        df["flag_taxon_target"]
                        & df["flag_valid_coordinates"]
                        & df["flag_reported_uncertainty_acceptable"]
                        & df["flag_no_fatal_gbif_issue"]
                        & df["flag_no_declared_geospatial_issue"]
                        & df["flag_no_cultivation_text"]
                        & df["flag_no_cultivation_establishment"]
                    ).sum()
                ),
                int(
                    (
                        df["flag_taxon_target"]
                        & df["flag_valid_coordinates"]
                        & df["flag_reported_uncertainty_acceptable"]
                        & df["flag_no_fatal_gbif_issue"]
                        & df["flag_no_declared_geospatial_issue"]
                        & df["flag_no_cultivation_text"]
                        & df["flag_no_cultivation_establishment"]
                        & df["flag_acceptable_basis"]
                    ).sum()
                ),
                int(strict_keep.sum()),
                len(candidates),
            ],
        }
    )

    candidates.to_csv(
        OUT_DIR / "01_metadata_QC_candidates.csv",
        index=False,
        encoding="utf-8-sig"
    )

    df[keep_columns + ["metadata_QC_keep"]].to_csv(
        OUT_DIR / "01_metadata_QC_all_records_audit.csv",
        index=False,
        encoding="utf-8-sig"
    )

    audit.to_csv(
        OUT_DIR / "01_metadata_QC_summary.csv",
        index=False,
        encoding="utf-8-sig"
    )

    year_summary = (
        candidates["record_year"]
        .dropna()
        .astype(int)
        .value_counts()
        .sort_index()
        .rename_axis("year")
        .reset_index(name="n_records")
    )

    year_summary.to_csv(
        OUT_DIR / "01_occurrence_year_distribution.csv",
        index=False,
        encoding="utf-8-sig"
    )

    run_info = pd.DataFrame(
        {
            "parameter": [
                "input_file",
                "target_species",
                "coordinate_uncertainty_cutoff_m",
                "geographic_filter",
                "cultivation_metadata_filter",
                "LULC_filter",
                "spatial_thinning",
            ],
            "value": [
                str(raw_file),
                "Camellia japonica",
                MAX_REPORTED_COORD_UNCERTAINTY_M,
                "Global; no native-range or M restriction",
                "Applied in this script",
                "Not yet applied; next stage",
                "Not yet applied; later stage at 5 km",
            ],
        }
    )

    run_info.to_csv(
        OUT_DIR / "01_run_information.csv",
        index=False,
        encoding="utf-8-sig"
    )

    print("\n" + "=" * 72)
    print("Global occurrence metadata QC completed successfully")
    print(f"Raw records: {len(raw)}")
    print(f"Metadata-QC candidates for CoordinateCleaner/LULC: {len(candidates)}")
    print(f"Output folder:\n{OUT_DIR}")
    print("=" * 72 + "\n")


if __name__ == "__main__":
    main()
