from __future__ import annotations
import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
import csv
import gzip
import hashlib
import json
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Iterable
import openpyxl

META_VERSION = "2026-07-25_v1.0"
MISSING = "not_available"
NOT_APPLICABLE = "not_applicable"
METADATA_ROOT = Path(__file__).resolve().parents[1]
INTERIM_DIR = METADATA_ROOT / "interim" / "geo_sample_metadata"
SOURCE_DIR = METADATA_ROOT / "source"
OUTPUT_DIR = METADATA_ROOT / "frozen" / META_VERSION
BY_DATASET_DIR = OUTPUT_DIR / "by_dataset"
MASTER_FIELDS = [
    "meta_version",
    "observation_id",
    "dataset",
    "sample_id",
    "sample_title",
    "patient_id",
    "visit_id",
    "timepoint",
    "disease",
    "disease_state",
    "disease_activity_group",
    "response",
    "response_definition",
    "tissue",
    "treatment",
    "treatment_stage",
    "SLEDAI",
    "renal_activity",
    "renal_class",
    "C3",
    "C3_unit",
    "C4",
    "C4_unit",
    "anti_dsDNA",
    "anti_dsDNA_unit",
    "proteinuria",
    "proteinuria_unit",
    "serum_creatinine",
    "serum_creatinine_unit",
    "age",
    "age_unit",
    "age_at_onset",
    "sex",
    "race",
    "cohort",
    "collection_year",
    "batch",
    "platform",
    "assay_type",
    "sc_library",
    "cell_count",
    "analysis_include",
    "analysis_role",
    "source_url",
    "source_header_file",
    "source_note",
]
GEO_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={}"


def text_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "TRUE" if value else "FALSE"
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value).strip()


def missing_code(value: object, *, not_applicable: bool = False) -> str:
    raw = text_value(value)
    if not_applicable:
        return NOT_APPLICABLE
    if raw.lower() in {
        "",
        "na",
        "n/a",
        "nd",
        "--",
        "none reported",
        "no data",
        "no val",
        "not reported",
        "not applicable",
    }:
        return MISSING
    return raw


def numeric_or_missing(value: object, *, not_applicable: bool = False) -> str:
    raw = missing_code(value, not_applicable=not_applicable)
    if raw in {MISSING, NOT_APPLICABLE}:
        return raw
    try:
        number = float(raw)
    except ValueError:
        return raw
    return str(int(number)) if number.is_integer() else format(number, ".12g")


def normalize_sex(value: object) -> str:
    raw = missing_code(value)
    if raw == MISSING:
        return raw
    mapping = {"f": "Female", "female": "Female", "m": "Male", "male": "Male"}
    return mapping.get(raw.lower(), raw)


def normalize_tissue(value: object) -> str:
    raw = missing_code(value)
    if raw == MISSING:
        return raw
    low = raw.lower()
    if "pbmc" in low or "peripheral blood mononuclear" in low:
        return "PBMC"
    if "whole blood" in low:
        return "Whole blood"
    if "glomer" in low:
        return "Kidney glomeruli"
    if "tubulo" in low or "interst" in low:
        return "Kidney tubulointerstitium"
    return raw


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: text_value(row.get(field, "")) for field in fields})


def geo_rows(dataset: str) -> list[dict[str, str]]:
    paths = sorted(INTERIM_DIR.glob(f"{dataset}*.sample_metadata.tsv"))
    if not paths:
        raise FileNotFoundError(f"No extracted GEO metadata found for {dataset}")
    rows: list[dict[str, str]] = []
    for path in paths:
        rows.extend(read_tsv(path))
    return rows


def blank_record(dataset: str, raw: dict[str, str]) -> dict[str, str]:
    record = {field: MISSING for field in MASTER_FIELDS}
    record.update(
        {
            "meta_version": META_VERSION,
            "observation_id": f"{dataset}:{raw['sample_id']}",
            "dataset": dataset,
            "sample_id": raw["sample_id"],
            "sample_title": missing_code(raw.get("sample_title")),
            "age_unit": "years",
            "analysis_include": "yes",
            "source_url": GEO_URL.format(dataset),
            "source_header_file": raw.get("source_header_file", MISSING),
        }
    )
    return record


def extract_gse135779_supplement() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    path = SOURCE_DIR / "supplementary" / "GSE135779_Nature_Supplementary_Tables.xlsx"
    workbook = openpyxl.load_workbook(path, read_only=True, data_only=True)

    def extract(sheet_name: str, header_row: int) -> dict[str, dict[str, str]]:
        sheet = workbook[sheet_name]
        header_values = next(
            sheet.iter_rows(min_row=header_row, max_row=header_row, values_only=True)
        )
        headers = [text_value(value) for value in header_values]
        result: dict[str, dict[str, str]] = {}
        for values in sheet.iter_rows(min_row=header_row + 1, values_only=True):
            key = text_value(values[0] if values else None)
            if not key:
                continue
            record = {header: text_value(value) for header, value in zip(headers, values) if header}
            result[key] = record
        return result

    clinical = extract("ST1b-Clinical information ", 4)
    sequencing = extract("ST1c-Sequencing information ", 4)
    workbook.close()
    return (clinical, sequencing)


def aggregate_gse135779_author_metadata() -> dict[str, dict[str, str]]:
    path = SOURCE_DIR / "supplementary" / "GSE135779_Meta_caSLE_processed_08092021_small.csv"
    aggregate: dict[str, dict[str, set[str]]] = defaultdict(
        lambda: {"Names": set(), "SLEDAI": set()}
    )
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            library_id = text_value(row.get("IDs"))
            if not library_id:
                continue
            aggregate[library_id]["Names"].add(text_value(row.get("Names")))
            aggregate[library_id]["SLEDAI"].add(text_value(row.get("SLEDAI")))
    result: dict[str, dict[str, str]] = {}
    for library_id, fields in aggregate.items():
        names = sorted((value for value in fields["Names"] if value))
        sledai_groups = sorted((value for value in fields["SLEDAI"] if value))
        if len(names) != 1 or len(sledai_groups) != 1:
            raise ValueError(
                f"GSE135779 author metadata is not library-constant for {library_id}: Names={names}; SLEDAI={sledai_groups}"
            )
        result[library_id] = {"Names": names[0], "SLEDAI": sledai_groups[0]}
    return result


def gse135779_records() -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    clinical, sequencing = extract_gse135779_supplement()
    author = aggregate_gse135779_author_metadata()
    rows = geo_rows("GSE135779")
    records: list[dict[str, str]] = []
    geo_names: set[str] = set()
    med_fields = [
        ("MMF", "Mycophenolate mofetil"),
        ("OS", "Oral steroids"),
        ("MTX", "Methotrexate"),
        ("Plaquenil", "Hydroxychloroquine"),
    ]
    for raw in rows:
        name = re.sub("\\s*\\[.*$", "", raw.get("sample_title", "")).strip()
        geo_names.add(name)
        clinical_row = clinical.get(name)
        sequencing_row = sequencing.get(name)
        if clinical_row is None or sequencing_row is None:
            raise ValueError(f"GSE135779 supplement match failed for {name}")
        library_match = re.search("\\[([^]]+)\\]", raw.get("sample_title", ""))
        if library_match is None:
            raise ValueError(f"GSE135779 library ID missing from title: {raw['sample_title']}")
        library_id = library_match.group(1)
        author_row = author.get(library_id)
        if author_row is None:
            raise ValueError(f"GSE135779 author metadata missing library {library_id}")
        if author_row["Names"] != name:
            raise ValueError(
                f"GSE135779 author/GEO subject mismatch for {library_id}: {author_row['Names']} != {name}"
            )
        group = clinical_row["Groups"]
        is_control = "HD" in group
        active_meds = [
            label
            for field, label in med_fields
            if text_value(clinical_row.get(field)) in {"1", "1.0"}
        ]
        medication_values = [text_value(clinical_row.get(field)) for field, _ in med_fields]
        if is_control:
            treatment = NOT_APPLICABLE
        elif active_meds:
            treatment = "; ".join(active_meds)
        elif all((value in {"0", "0.0"} for value in medication_values)):
            treatment = "No listed medication"
        else:
            treatment = MISSING
        record = blank_record("GSE135779", raw)
        record.update(
            {
                "patient_id": name,
                "visit_id": "V1",
                "timepoint": "Cross-sectional",
                "disease": "Healthy control" if is_control else "SLE",
                "disease_state": "Healthy control" if is_control else "SLE",
                "disease_activity_group": (
                    NOT_APPLICABLE if is_control else missing_code(author_row["SLEDAI"])
                ),
                "response": NOT_APPLICABLE if is_control else MISSING,
                "response_definition": NOT_APPLICABLE if is_control else MISSING,
                "tissue": "PBMC",
                "treatment": treatment,
                "treatment_stage": NOT_APPLICABLE if is_control else "Cross-sectional",
                "SLEDAI": numeric_or_missing(clinical_row.get("SLEDAI"), not_applicable=is_control),
                "renal_activity": missing_code(
                    clinical_row.get("Neph_all"), not_applicable=is_control
                ),
                "renal_class": missing_code(
                    clinical_row.get("Nephrology_Class"), not_applicable=is_control
                ),
                "C3": numeric_or_missing(clinical_row.get("C3"), not_applicable=is_control),
                "C4": numeric_or_missing(clinical_row.get("C4"), not_applicable=is_control),
                "anti_dsDNA": numeric_or_missing(
                    clinical_row.get("DSDNA"), not_applicable=is_control
                ),
                "age": numeric_or_missing(clinical_row.get("Age")),
                "sex": normalize_sex(clinical_row.get("Gender")),
                "race": missing_code(clinical_row.get("Race")),
                "cohort": "Adult" if group.startswith("a") else "Child",
                "collection_year": numeric_or_missing(clinical_row.get("Collection_year")),
                "batch": missing_code(clinical_row.get("Batch")),
                "platform": raw.get("sample_platform_id", MISSING),
                "assay_type": "scRNA-seq",
                "sc_library": library_id,
                "cell_count": numeric_or_missing(sequencing_row.get("Estimated Number of Cells")),
                "analysis_role": "single_cell_discovery",
                "source_note": "GEO sample record; Nature Supplementary Tables 1b/1c; author cell metadata aggregated to library for published SLEDAI group",
            }
        )
        records.append(record)
    exclusions = []
    for name in sorted(set(clinical) - geo_names):
        exclusions.append(
            {
                "dataset": "GSE135779",
                "sample_id": MISSING,
                "patient_id": name,
                "exclusion_level": "published_individual",
                "reason": "No public GSM/count matrix in the current GEO processed archive",
                "decision": "Excluded from frozen sample-level master; retained in source supplement",
                "source": "Nature Supplementary Table 1b versus GEO GSE135779",
            }
        )
    clinical_rows = []
    for name, row in clinical.items():
        clinical_rows.append({"Names": name, **row})
    sequencing_rows = []
    for name, row in sequencing.items():
        sequencing_rows.append({"SampleID": name, **row})
    return (
        records,
        exclusions,
        [{"table": "clinical", **row} for row in clinical_rows]
        + [{"table": "sequencing", **row} for row in sequencing_rows],
    )


def gse65391_records() -> list[dict[str, str]]:
    records = []
    for raw in geo_rows("GSE65391"):
        is_control = raw.get("characteristic_disease_state", "").lower() == "healthy"
        is_technical = raw.get("characteristic_set", "").lower() == "technical_replicate"
        visit = missing_code(raw.get("characteristic_visit"))
        timepoint = missing_code(raw.get("characteristic_cumulative_time"))
        if timepoint == MISSING:
            timepoint = f"Visit {visit}"
        record = blank_record("GSE65391", raw)
        record.update(
            {
                "patient_id": missing_code(raw.get("characteristic_subject")),
                "visit_id": f"V{visit}",
                "timepoint": timepoint,
                "disease": "Healthy control" if is_control else "SLE",
                "disease_state": missing_code(raw.get("characteristic_disease_state")),
                "disease_activity_group": missing_code(
                    raw.get("characteristic_disease_activity"), not_applicable=is_control
                ),
                "response": NOT_APPLICABLE if is_control else MISSING,
                "response_definition": NOT_APPLICABLE if is_control else MISSING,
                "tissue": "Whole blood",
                "treatment": missing_code(
                    raw.get("characteristic_treatment"), not_applicable=is_control
                ),
                "treatment_stage": "Control sampling" if is_control else "Longitudinal visit",
                "SLEDAI": numeric_or_missing(
                    raw.get("characteristic_sledai"), not_applicable=is_control
                ),
                "renal_activity": missing_code(
                    raw.get("characteristic_renal"), not_applicable=is_control
                ),
                "renal_class": missing_code(
                    raw.get("characteristic_nephritis_class"), not_applicable=is_control
                ),
                "C3": numeric_or_missing(raw.get("characteristic_c3"), not_applicable=is_control),
                "C4": numeric_or_missing(raw.get("characteristic_c4"), not_applicable=is_control),
                "anti_dsDNA": numeric_or_missing(
                    raw.get("characteristic_ds_dna"), not_applicable=is_control
                ),
                "age": numeric_or_missing(raw.get("characteristic_age")),
                "sex": normalize_sex(raw.get("characteristic_gender")),
                "race": missing_code(raw.get("characteristic_race")),
                "cohort": "Child",
                "batch": missing_code(raw.get("characteristic_batch")),
                "platform": raw.get("sample_platform_id", MISSING),
                "assay_type": "Microarray",
                "sc_library": NOT_APPLICABLE,
                "cell_count": NOT_APPLICABLE,
                "analysis_include": "no" if is_technical else "yes",
                "analysis_role": (
                    "technical_replicate" if is_technical else "pediatric_longitudinal_validation"
                ),
                "source_note": "GEO sample record; clinical variables deposited as Sample_characteristics_ch1",
            }
        )
        records.append(record)
    return records


def gse121239_records() -> list[dict[str, str]]:
    records = []
    for raw in geo_rows("GSE121239"):
        is_control = raw.get("characteristic_disease_state", "").lower() == "healthy"
        patient = missing_code(raw.get("characteristic_patient_id"))
        if patient == MISSING:
            patient = re.sub("^PBMC_", "", raw.get("sample_title", ""))
        visit_date = missing_code(raw.get("characteristic_visit_date"))
        visit_id = (
            "V1" if visit_date == MISSING else re.sub("[^A-Za-z0-9]+", "_", visit_date).strip("_")
        )
        record = blank_record("GSE121239", raw)
        record.update(
            {
                "patient_id": patient,
                "visit_id": visit_id,
                "timepoint": "Baseline control" if is_control else visit_date,
                "disease": "Healthy control" if is_control else "SLE",
                "disease_state": missing_code(raw.get("characteristic_disease_state")),
                "response": NOT_APPLICABLE if is_control else MISSING,
                "response_definition": NOT_APPLICABLE if is_control else MISSING,
                "tissue": normalize_tissue(raw.get("characteristic_tissue")),
                "treatment": NOT_APPLICABLE if is_control else MISSING,
                "treatment_stage": "Control sampling" if is_control else "Longitudinal visit",
                "SLEDAI": numeric_or_missing(
                    raw.get("characteristic_sledai"), not_applicable=is_control
                ),
                "renal_activity": NOT_APPLICABLE if is_control else MISSING,
                "cohort": "Adult",
                "platform": raw.get("sample_platform_id", MISSING),
                "assay_type": "Microarray",
                "sc_library": NOT_APPLICABLE,
                "cell_count": NOT_APPLICABLE,
                "analysis_role": "adult_longitudinal_activity_validation",
                "source_note": "GEO sample record; control patient IDs reconstructed only from explicit sample titles",
            }
        )
        records.append(record)
    return records


GSE200306_RESPONSE_DEFINITION = "End of induction (~6 months): CR = proteinuria <0.5 g/day with normalized serum creatinine; PR = >50% proteinuria reduction to 0.5-<3 g/day with stable/improved serum creatinine; NR = <50% reduction or >3 g/day with stable/worsening serum creatinine (PMID 35788359)."
SRI4_DEFINITION = "SRI-4: >=4-point SLEDAI-2K improvement, no new BILAG A or >1 new BILAG B, and no clinically important worsening in physician global assessment; GEO deposited YES/NO status was mapped to R/NR (PMID 38260685)."


def gse8888x_records(dataset: str) -> list[dict[str, str]]:
    records = []
    for raw in geo_rows(dataset):
        group = raw.get("characteristic_group", "")
        is_control = group.lower() == "normal"
        timepoint = missing_code(raw.get("characteristic_time"))
        treatment_raw = missing_code(
            raw.get("characteristic_treatment_q2w_once_in_2_weeks_q4w_once_in_4_weeks")
        )
        if is_control:
            treatment = NOT_APPLICABLE
        elif dataset == "GSE88884":
            treatment = MISSING
        elif treatment_raw.startswith("LY "):
            treatment = treatment_raw.replace("LY ", "Tabalumab ", 1)
        else:
            treatment = treatment_raw
        if is_control:
            stage = NOT_APPLICABLE
        elif timepoint == "baseline":
            stage = "Pre-treatment"
        elif timepoint.startswith("week"):
            stage = f"On-treatment {timepoint}"
        else:
            stage = timepoint
        record = blank_record(dataset, raw)
        record.update(
            {
                "patient_id": missing_code(raw.get("characteristic_subject_id")),
                "visit_id": timepoint,
                "timepoint": timepoint,
                "disease": "Healthy control" if is_control else "SLE",
                "disease_state": "Healthy control" if is_control else "SLE",
                "response": NOT_APPLICABLE if is_control else MISSING,
                "response_definition": NOT_APPLICABLE if is_control else MISSING,
                "tissue": "Whole blood",
                "treatment": treatment,
                "treatment_stage": stage,
                "SLEDAI": numeric_or_missing(
                    raw.get("characteristic_sledai_at_baseline"), not_applicable=is_control
                ),
                "C3": numeric_or_missing(
                    raw.get("characteristic_c3_at_baseline_g_l"), not_applicable=is_control
                ),
                "C3_unit": NOT_APPLICABLE if is_control else "g/L",
                "C4": numeric_or_missing(
                    raw.get("characteristic_c4_at_baseline_g_l"), not_applicable=is_control
                ),
                "C4_unit": NOT_APPLICABLE if is_control else "g/L",
                "anti_dsDNA": numeric_or_missing(
                    raw.get("characteristic_antidsdna_at_baseline_iu"), not_applicable=is_control
                ),
                "anti_dsDNA_unit": NOT_APPLICABLE if is_control else "IU (as deposited)",
                "age": numeric_or_missing(raw.get("characteristic_age_at_baseline")),
                "sex": normalize_sex(raw.get("characteristic_sex")),
                "race": missing_code(raw.get("characteristic_race")),
                "cohort": "Adult",
                "batch": missing_code(raw.get("characteristic_batch")),
                "platform": raw.get("sample_platform_id", MISSING),
                "assay_type": "Microarray",
                "sc_library": NOT_APPLICABLE,
                "cell_count": NOT_APPLICABLE,
                "analysis_role": (
                    "baseline_case_control"
                    if dataset == "GSE88884"
                    else "tabalumab_pharmacodynamics"
                ),
                "source_note": "GEO sample record from ILLUMINATE-1/2; response labels are not publicly deposited",
            }
        )
        records.append(record)
    return records
