from __future__ import annotations
import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
import csv
import gzip
import hashlib
import json
import logging
import math
import platform
import re
import statistics
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("F:\\program\\Hei\\SLE")
OUT = Path(__file__).resolve().parent
PREVIOUS = Path(
    "F:\\program\\Hei\\SLE_v1\\results\\reviewer_response\\08_trial_reporting_ancova_permutation"
)
DATASETS = {"GSE88884": None, "GSE88885": "ILLUMINATE-1", "GSE88886": "ILLUMINATE-2"}
EXPECTED = {"GSE88884": 1820, "GSE88885": 908, "GSE88886": 418}
MISSING = {"", "--", "na", "n/a", "nan", "null", "not available", "unknown"}
FIELD_PREFIXES = {
    "subject_id": "subject_id",
    "time": "visit",
    "group": "condition",
    "treatment (": "treatment_raw",
    "batch": "trial_or_batch",
    "race": "race",
    "Sex": "sex",
    "age_at_baseline": "age_at_baseline",
    "region": "region",
    "sledai_at_baseline": "sledai_at_baseline",
    "adna_at_baseline (": "anti_dsdna_status_at_baseline",
    "antidsdna_at_baseline (": "anti_dsdna_iu_at_baseline",
    "c3base (": "c3_category_at_baseline",
    "c3_at_baseline (": "c3_g_l_at_baseline",
    "c4base (": "c4_category_at_baseline",
    "c4_at_baseline (": "c4_g_l_at_baseline",
    "tissue": "tissue",
}
NUMERIC = [
    "age_at_baseline",
    "sledai_at_baseline",
    "anti_dsdna_iu_at_baseline",
    "c3_g_l_at_baseline",
    "c4_g_l_at_baseline",
]
CATEGORICAL = [
    "sex",
    "race",
    "region",
    "anti_dsdna_status_at_baseline",
    "c3_category_at_baseline",
    "c4_category_at_baseline",
]
CLINICAL = NUMERIC + CATEGORICAL
STABLE = ["condition", "trial_or_batch", "treatment_raw", "tissue"] + CLINICAL
ARMS = ["Placebo", "Tabalumab_Q2W", "Tabalumab_Q4W"]
VISITS = ["baseline", "week16", "week52"]
SOURCES: dict[str, dict] = {}
CHECKS: list[dict] = []
INVALID: list[dict] = []
CONFLICTS: list[dict] = []
PARTIAL: list[dict] = []
MISMATCHES: list[dict] = []


def write_tsv(name, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys((k for row in rows for k in row)))
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="raise"
        )
        writer.writeheader()
        writer.writerows(rows)
    return rows


def read_tsv(path):
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256(path):
    dig = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            dig.update(block)
    return dig.hexdigest()


def register(path, role):
    path = path.resolve()
    if str(path) not in SOURCES:
        SOURCES[str(path)] = {
            "path": str(path),
            "role": role,
            "bytes": path.stat().st_size,
            "mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
            "sha256_before": sha256(path),
        }
    return path


def check(name, passed, details, fatal=True):
    CHECKS.append({"check": name, "status": "PASS" if passed else "FAIL", "details": str(details)})
    logging.info("%s %s %s", "PASS" if passed else "FAIL", name, details)
    if not passed and fatal:
        write_tsv("validation_checks.tsv", CHECKS)
        raise ValueError(f"{name}: {details}")


def missing(value):
    return value is None or str(value).strip().lower() in MISSING


def alias(key):
    hits = [
        v for p, v in FIELD_PREFIXES.items() if key == p or (p.endswith("(") and key.startswith(p))
    ]
    if len(hits) != 1:
        raise ValueError(f"Unmapped or ambiguous characteristic key: {key}")
    return hits[0]


def normalize(field, value, dataset, sample_id):
    if missing(value):
        return None
    value = value.strip()
    if field in NUMERIC:
        try:
            result = float(value)
            if not math.isfinite(result):
                raise ValueError("not finite")
            return result
        except ValueError:
            INVALID.append(
                {"dataset": dataset, "sample_id": sample_id, "field": field, "raw_value": value}
            )
            return None
    if field == "condition":
        if value.lower() in {"normal", "healthy"}:
            return "Healthy"
        if value.lower() in {"sle", "systemic lupus erythematosus (sle)"}:
            return "SLE"
    if field == "visit":
        return value.lower()
    return value


def arm(raw):
    if missing(raw):
        return "Unknown"
    return {
        "Placebo": "Placebo",
        "LY 120 mg Q2W": "Tabalumab_Q2W",
        "LY 120 mg Q4W": "Tabalumab_Q4W",
    }.get(raw, "UNMAPPED:" + raw)


def parse_geo(dataset):
    path = register(
        ROOT / f"{dataset}_pipeline/data/raw/downloads/{dataset}_series_matrix.txt.gz",
        "raw GEO series matrix cache",
    )
    annotation_rows = []
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as fh:
        for line_no, line in enumerate(fh, 1):
            if line.startswith("!series_matrix_table_begin"):
                break
            cells = next(csv.reader([line], delimiter="\t"))
            if cells and cells[0].startswith("!Sample_"):
                annotation_rows.append((line_no, cells[0], cells[1:]))
    sample_ids = next(
        (values for _, key, values in annotation_rows if key == "!Sample_geo_accession")
    )
    check(
        dataset + "_sample_accessions_unique",
        len(set(sample_ids)) == len(sample_ids),
        len(sample_ids),
    )
    check(dataset + "_sample_count", len(sample_ids) == EXPECTED[dataset], len(sample_ids))
    samples = [{"dataset": dataset, "sample_id": sid} for sid in sample_ids]
    inventory, long_rows = ([], [])
    for line_no, key, values in annotation_rows:
        check_key = f"{dataset}_annotation_width_line_{line_no}"
        check(check_key, len(values) == len(samples), f"{key}: {len(values)}")
        if key != "!Sample_characteristics_ch1":
            if key in {"!Sample_title", "!Sample_source_name_ch1", "!Sample_supplementary_file"}:
                suffix = {
                    "!Sample_title": "title",
                    "!Sample_source_name_ch1": "source_name",
                    "!Sample_supplementary_file": "supplementary_file",
                }[key]
                for sample, value in zip(samples, values):
                    sample["geo_" + suffix] = value
            continue
        pairs = [value.rsplit(": ", 1) for value in values]
        check(
            f"{dataset}_key_parse_line_{line_no}",
            all((len(p) == 2 for p in pairs)),
            "split on final colon-space, preserving colon within field definition",
        )
        keys = {p[0] for p in pairs}
        check(f"{dataset}_single_key_line_{line_no}", len(keys) == 1, sorted(keys))
        raw_key = next(iter(keys))
        field = alias(raw_key)
        inventory.append(
            {
                "dataset": dataset,
                "source_file": str(path),
                "decompressed_line": line_no,
                "raw_key": raw_key,
                "normalized_field": field,
                "raw_column": "geo_raw_" + field,
                "data_type": "numeric" if field in NUMERIC else "categorical/text",
                "missing_raw_tokens": " | ".join(sorted({p[1] for p in pairs if missing(p[1])})),
                "sample_count": len(samples),
                "nonmissing_count": sum((not missing(p[1]) for p in pairs)),
            }
        )
        for sample, (_, value) in zip(samples, pairs):
            sample["geo_raw_" + field] = value
            sample[field] = normalize(field, value, dataset, sample["sample_id"])
            long_rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample["sample_id"],
                    "source_file": str(path),
                    "decompressed_line": line_no,
                    "raw_key": raw_key,
                    "normalized_field": field,
                    "raw_value": value,
                }
            )
    for sample in samples:
        sample.setdefault("treatment_raw", None)
        sample["treatment"] = arm(sample["treatment_raw"])
        treatment_text = sample["treatment_raw"] or ""
        dose_match = re.search("\\b([0-9]+(?:\\.[0-9]+)?)\\s*mg\\b", treatment_text, flags=re.I)
        interval_match = re.search("\\bQ([24])W\\b", treatment_text, flags=re.I)
        sample["stated_dose_mg"] = float(dose_match.group(1)) if dose_match else None
        sample["stated_interval_weeks"] = int(interval_match.group(1)) if interval_match else None
        sample["dose_interpretation"] = (
            "directly stated in treatment annotation"
            if dose_match
            else "not a stated active-drug dose"
        )
    return (samples, inventory, long_rows, annotation_rows)


def join_audit(dataset, samples):
    path = register(
        ROOT / f"{dataset}_pipeline/results/00_input_audit/sample_audit.tsv",
        "existing pipeline sample audit",
    )
    audit = read_tsv(path)
    index = {row["sample_id"]: row for row in audit}
    check(dataset + "_audit_unique", len(index) == len(audit), len(audit))
    raw_ids = {row["sample_id"] for row in samples}
    check(
        dataset + "_exact_GSM_join",
        raw_ids == set(index),
        json.dumps(
            {"raw_only": sorted(raw_ids - set(index)), "audit_only": sorted(set(index) - raw_ids)}
        ),
    )
    for row in samples:
        old = index[row["sample_id"]]
        for key, value in old.items():
            if key != "sample_id":
                row["audit_" + key] = value
        comparisons = {
            "subject_id": row["subject_id"],
            "time": row["visit"],
            "condition": row["condition"],
            "sample_title": row["geo_title"],
            "source_name": row["geo_source_name"],
            "supplementary_file": row["geo_supplementary_file"],
            "treatment": row["treatment"],
        }
        for field, value in comparisons.items():
            if value != old[field]:
                MISMATCHES.append(
                    {
                        "dataset": dataset,
                        "sample_id": row["sample_id"],
                        "field": field,
                        "raw_derived_value": value,
                        "sample_audit_value": old[field],
                    }
                )
        check_trial = DATASETS[dataset]
        if check_trial is not None and row["trial_or_batch"] != check_trial:
            MISMATCHES.append(
                {
                    "dataset": dataset,
                    "sample_id": row["sample_id"],
                    "field": "trial",
                    "raw_derived_value": row["trial_or_batch"],
                    "sample_audit_value": check_trial,
                }
            )
    n_bad = sum((x["dataset"] == dataset for x in MISMATCHES))
    check(dataset + "_identity_visit_arm_agreement", n_bad == 0, n_bad, fatal=False)


def unique_nonmissing(rows, field):
    return sorted({row.get(field) for row in rows if row.get(field) is not None}, key=str)


def subject_rows(samples):
    groups = defaultdict(list)
    for row in samples:
        groups[row["dataset"], row["trial_or_batch"], row["subject_id"]].append(row)
    subjects = []
    provenance = []
    for (dataset, trial, sid), rows in sorted(groups.items()):
        rows.sort(
            key=lambda row: (
                VISITS.index(row["visit"]) if row["visit"] in VISITS else 99,
                row["sample_id"],
            )
        )
        times = Counter((row["visit"] for row in rows))
        check(f"{dataset}_{sid}_unique_subject_visit", max(times.values()) == 1, json.dumps(times))
        baseline = [row for row in rows if row["visit"] == "baseline"]
        sub = {
            "dataset": dataset,
            "trial_or_batch": trial,
            "subject_id": sid,
            "subject_key": "|".join((dataset, trial, sid)),
            "n_molecular_samples": len(rows),
            "sample_ids": "|".join((row["sample_id"] for row in rows)),
            "observed_visits": "|".join((v for v in VISITS if times[v])),
            "baseline_sample_present": int(bool(baseline)),
            "week16_sample_present": int(bool(times["week16"])),
            "week52_sample_present": int(bool(times["week52"])),
            "complete_baseline_week16": int(bool(baseline) and bool(times["week16"])),
            "complete_baseline_week52": int(bool(baseline) and bool(times["week52"])),
            "baseline_sample_id": baseline[0]["sample_id"] if baseline else None,
            "baseline_annotation_context": (
                "explicit baseline molecular sample"
                if baseline
                else "baseline-labelled annotations carried on follow-up sample(s)"
            ),
            "n_conflicting_fields": 0,
        }
        for field in STABLE:
            values = unique_nonmissing(rows, field)
            base_values = unique_nonmissing(baseline, field)
            available_rows = [row for row in rows if row.get(field) is not None]
            if len(values) > 1:
                sub[field] = None
                sub["n_conflicting_fields"] += 1
                source = "conflict; no value selected"
                CONFLICTS.append(
                    {
                        "dataset": dataset,
                        "trial_or_batch": trial,
                        "subject_id": sid,
                        "field": field,
                        "values": json.dumps(values, ensure_ascii=False),
                        "sample_values": json.dumps(
                            {row["sample_id"]: row.get(field) for row in rows}, ensure_ascii=False
                        ),
                    }
                )
            elif values:
                sub[field] = values[0]
                source = (
                    "explicit baseline annotation"
                    if base_values
                    else "consistent baseline-labelled annotation on follow-up sample(s)"
                )
                if len(available_rows) < len(rows):
                    PARTIAL.append(
                        {
                            "dataset": dataset,
                            "trial_or_batch": trial,
                            "subject_id": sid,
                            "field": field,
                            "n_samples": len(rows),
                            "n_nonmissing_samples": len(available_rows),
                            "selected_source": source,
                        }
                    )
            else:
                sub[field] = None
                source = "all annotations missing / field absent"
            provenance.append(
                {
                    "dataset": dataset,
                    "trial_or_batch": trial,
                    "subject_id": sid,
                    "field": field,
                    "selected_value": sub[field],
                    "selection_source": source,
                    "n_nonmissing_distinct_values": len(values),
                    "supporting_sample_ids": "|".join((row["sample_id"] for row in available_rows)),
                }
            )
        sub["treatment"] = arm(sub["treatment_raw"])
        dose = re.search(
            "\\b([0-9]+(?:\\.[0-9]+)?)\\s*mg\\b", sub["treatment_raw"] or "", flags=re.I
        )
        interval = re.search("\\bQ([24])W\\b", sub["treatment_raw"] or "", flags=re.I)
        sub["stated_dose_mg"] = float(dose.group(1)) if dose else None
        sub["stated_interval_weeks"] = int(interval.group(1)) if interval else None
        subjects.append(sub)
    return (subjects, provenance)


def quantile(values, p):
    values = sorted(values)
    at = (len(values) - 1) * p
    lo, hi = (math.floor(at), math.ceil(at))
    return values[lo] + (at - lo) * (values[hi] - values[lo])


def describe(subjects):
    numeric, categorical, missing_rows, flow, readable = ([], [], [], [], [])
    sets = {
        "all_available_molecular": lambda r: True,
        "observed_baseline_molecular": lambda r: bool(r["baseline_sample_present"]),
        "complete_baseline_week16": lambda r: bool(r["complete_baseline_week16"]),
        "complete_baseline_week52": lambda r: bool(r["complete_baseline_week52"]),
    }
    trials = sorted({s["trial_or_batch"] for s in subjects if s["dataset"] != "GSE88884"})
    for trial in trials:
        trial_subjects = [
            s for s in subjects if s["dataset"] != "GSE88884" and s["trial_or_batch"] == trial
        ]
        for treatment in ARMS:
            members = [s for s in trial_subjects if s["treatment"] == treatment]
            for visit in ["week16", "week52"]:
                flow.append(
                    {
                        "trial": trial,
                        "treatment": treatment,
                        "followup": visit,
                        "transcriptomic_cohort_subjects": len(members),
                        "baseline_available": sum((s["baseline_sample_present"] for s in members)),
                        "followup_available": sum((s[visit + "_sample_present"] for s in members)),
                        "paired_baseline_followup": sum(
                            (s["complete_baseline_" + visit] for s in members)
                        ),
                        "baseline_missing": sum(
                            (not s["baseline_sample_present"] for s in members)
                        ),
                        "followup_missing_among_baseline": sum(
                            (
                                s["baseline_sample_present"] and (not s[visit + "_sample_present"])
                                for s in members
                            )
                        ),
                    }
                )
        for analysis_set, predicate in sets.items():
            readable_cells = {}
            for treatment in ARMS:
                members = [
                    s for s in trial_subjects if s["treatment"] == treatment and predicate(s)
                ]
                common = {
                    "trial": trial,
                    "analysis_set": analysis_set,
                    "treatment": treatment,
                    "n_subjects": len(members),
                }
                readable_cells["Subjects", treatment] = str(len(members))
                for field in CLINICAL:
                    values = [s[field] for s in members if s[field] is not None]
                    miss = len(members) - len(values)
                    base = {
                        **common,
                        "field": field,
                        "n_nonmissing": len(values),
                        "n_missing": miss,
                    }
                    missing_rows.append(
                        {**base, "percent_missing": 100 * miss / len(members) if members else None}
                    )
                    if field in NUMERIC:
                        result = {
                            **base,
                            "mean": statistics.mean(values) if values else None,
                            "sd": statistics.stdev(values) if len(values) > 1 else None,
                            "median": statistics.median(values) if values else None,
                            "q1": quantile(values, 0.25) if values else None,
                            "q3": quantile(values, 0.75) if values else None,
                            "min": min(values) if values else None,
                            "max": max(values) if values else None,
                        }
                        numeric.append(result)
                        if values:
                            readable_cells[field, treatment] = (
                                f"{result['mean']:.2f} ({result['sd']:.2f}); median {result['median']:.2f} [{result['q1']:.2f}, {result['q3']:.2f}]; n={len(values)}, missing={miss}"
                            )
                        else:
                            readable_cells[field, treatment] = f"Unavailable; missing={miss}"
                    else:
                        levels = sorted({s[field] for s in trial_subjects if s[field] is not None})
                        counts = Counter(values)
                        for level in levels:
                            percent = 100 * counts[level] / len(values) if values else None
                            categorical.append(
                                {
                                    **base,
                                    "level": level,
                                    "n": counts[level],
                                    "percent_of_nonmissing": percent,
                                    "percent_of_all_subjects": (
                                        100 * counts[level] / len(members) if members else None
                                    ),
                                }
                            )
                            readable_cells[field + ": " + level, treatment] = (
                                f"{counts[level]}/{len(values)} ({percent:.1f}%); missing={miss}"
                                if percent is not None
                                else f"Unavailable; missing={miss}"
                            )
            for label in dict.fromkeys((label for label, _ in readable_cells)):
                readable.append(
                    {
                        "trial": trial,
                        "analysis_set": analysis_set,
                        "characteristic": label,
                        **{t: readable_cells.get((label, t), "") for t in ARMS},
                    }
                )
    return (numeric, categorical, missing_rows, flow, readable)


def validate_previous_flow(flow, subjects):
    path = register(
        PREVIOUS / "trial_transcriptomic_flow.tsv", "previous trial reporting flow output"
    )
    previous = read_tsv(path)
    keys = ["trial", "treatment", "followup"]
    current_map = {tuple((row[k] for k in keys)): row for row in flow}
    check("previous_flow_row_count", len(previous) == len(flow), len(flow))
    deltas = []
    for old in previous:
        key = tuple((old[k] for k in keys))
        if key not in current_map:
            deltas.append(
                {"key": str(key), "field": "key", "previous": "present", "recovered": "absent"}
            )
            continue
        for field in old:
            if field not in keys and int(old[field]) != current_map[key][field]:
                deltas.append(
                    {
                        "key": str(key),
                        "field": field,
                        "previous": old[field],
                        "recovered": current_map[key][field],
                    }
                )
    write_tsv("previous_flow_discrepancies.tsv", deltas, ["key", "field", "previous", "recovered"])
    check("previous_flow_exact_agreement", not deltas, len(deltas), fatal=False)
    pair_path = register(
        PREVIOUS / "complete_pair_coordinate_data.tsv",
        "previous ANCOVA complete-pair subject coordinates",
    )
    pairs = read_tsv(pair_path)
    old_pairs = {(r["trial"], r["subject_id"], r["treatment"], r["time"]) for r in pairs}
    new_pairs = {
        (r["trial_or_batch"], r["subject_id"], r["treatment"], visit)
        for r in subjects
        if r["dataset"] != "GSE88884"
        for visit in ("week16", "week52")
        if r["complete_baseline_" + visit]
    }
    pair_deltas = [
        {"side": side, "trial": key[0], "subject_id": key[1], "treatment": key[2], "visit": key[3]}
        for side, keys2 in (
            ("previous_only", old_pairs - new_pairs),
            ("recovered_only", new_pairs - old_pairs),
        )
        for key in sorted(keys2)
    ]
    write_tsv(
        "previous_pair_subject_discrepancies.tsv",
        pair_deltas,
        ["side", "trial", "subject_id", "treatment", "visit"],
    )
    check(
        "previous_complete_pair_subject_agreement",
        not pair_deltas and len(old_pairs) == len(pairs),
        f"unique pairs={len(new_pairs)}; discrepancies={len(pair_deltas)}",
        fatal=False,
    )


def cross_dataset_links(samples, subjects):
    baseline_index = {
        (s["trial_or_batch"], s["subject_id"]): s for s in subjects if s["dataset"] == "GSE88884"
    }
    baseline_samples = {s["sample_id"]: s for s in samples if s["dataset"] == "GSE88884"}
    raw_samples = {s["sample_id"]: s for s in samples}
    links, differences = ([], [])
    for subject in subjects:
        if subject["dataset"] == "GSE88884":
            continue
        key = (subject["trial_or_batch"], subject["subject_id"])
        other = baseline_index.get(key)
        link = {
            "dataset": subject["dataset"],
            "trial": key[0],
            "subject_id": key[1],
            "treatment": subject["treatment"],
            "longitudinal_baseline_sample_id": subject["baseline_sample_id"],
            "gse88884_subject_match": int(other is not None),
            "gse88884_baseline_sample_id": other["baseline_sample_id"] if other else None,
            "baseline_CEL_exact_match": None,
            "n_clinical_fields_equal": 0,
            "n_clinical_fields_conflicting": 0,
            "n_clinical_fields_one_side_missing": 0,
            "n_clinical_fields_both_missing": 0,
        }
        if other:
            old_sample = baseline_samples[other["baseline_sample_id"]]
            if subject["baseline_sample_id"]:
                cur_sample = raw_samples[subject["baseline_sample_id"]]
                link["baseline_CEL_exact_match"] = int(
                    old_sample["audit_expression_column"] == cur_sample["audit_expression_column"]
                )
            for field in CLINICAL:
                value, old_value = (subject[field], other[field])
                if value is None and old_value is None:
                    link["n_clinical_fields_both_missing"] += 1
                elif value is None or old_value is None:
                    link["n_clinical_fields_one_side_missing"] += 1
                elif value == old_value:
                    link["n_clinical_fields_equal"] += 1
                else:
                    link["n_clinical_fields_conflicting"] += 1
                    differences.append(
                        {
                            "trial": key[0],
                            "subject_id": key[1],
                            "field": field,
                            "longitudinal_annotation": value,
                            "gse88884_annotation": old_value,
                        }
                    )
        links.append(link)
    return (links, differences)


def search_fields(dataset, annotation_rows, inventory):
    patterns = {
        "SRI5": "\\bSRI[\\s_-]*5\\b",
        "SRI4": "\\bSRI[\\s_-]*4\\b",
        "SRI_unspecified": "\\bSRI\\b",
        "responder": "responder|response[ _-]*(?:status|label)|clinical[ _-]*response",
        "BILAG": "\\bBILAG\\b",
        "physician_global_PGA": "physician.{0,20}global|\\bPGA\\b",
        "steroid_or_prednisone": "steroid|predni[sz]o|glucocorticoid",
        "active_drug_mg_dose": "\\b\\d+(?:\\.\\d+)?\\s*mg\\b",
        "Q2W_Q4W_frequency": "\\bQ[24]W\\b",
    }
    result = []
    for label, pattern in patterns.items():
        regex = re.compile(pattern, re.I)
        matches = [
            (ln, key, vals)
            for ln, key, vals in annotation_rows
            if regex.search(key + " " + " ".join(vals))
        ]
        matching_keys = sorted({f["raw_key"] for f in inventory if regex.search(f["raw_key"])})
        values = sorted({v for _, _, vals in matches for v in vals if regex.search(v)})
        result.append(
            {
                "dataset": dataset,
                "query": label,
                "regex": pattern,
                "scope": "all local !Sample_* metadata rows before expression table; field keys and values",
                "matched_annotation_rows": len(matches),
                "decompressed_lines": "|".join((str(x[0]) for x in matches)),
                "matched_characteristic_keys": "|".join(matching_keys),
                "matching_value_examples": " | ".join(values[:10]),
                "interpretation": (
                    "found in local sample annotation"
                    if matches
                    else "not found in these locally cached sample annotations; no claim about other repositories"
                ),
            }
        )
    return result


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[logging.FileHandler(OUT / "run.log", mode="w", encoding="utf-8")],
    )
    logging.info("BEGIN %s", Path(__file__).resolve())
    versions = {
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "dependencies": "Python standard library only",
        "run_started_local": datetime.now().astimezone().isoformat(),
        "script_sha256": sha256(Path(__file__)),
        "source_policy": "read only; no network; no model fitting or manuscript edits",
    }
    (OUT / "software_and_run.json").write_text(
        json.dumps(versions, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    shared = Path(__file__)
    register(
        PREVIOUS / "analysis_set_definitions.tsv",
        "previous manuscript model analysis-set definitions",
    )
    samples, inventory, raw_long, search = ([], [], [], [])
    for dataset in DATASETS:
        sample, inv, raw, headers = parse_geo(dataset)
        join_audit(dataset, sample)
        samples.extend(sample)
        inventory.extend(inv)
        raw_long.extend(raw)
        search.extend(search_fields(dataset, headers, inv))
    write_tsv("characteristic_field_mapping.tsv", inventory)
    write_tsv("sample_characteristics_all_raw_long.tsv", raw_long)
    write_tsv("sample_level_clinical_with_audit.tsv", samples)
    write_tsv(
        "sample_audit_mismatches.tsv",
        MISMATCHES,
        ["dataset", "sample_id", "field", "raw_derived_value", "sample_audit_value"],
    )
    write_tsv(
        "invalid_numeric_annotations.tsv", INVALID, ["dataset", "sample_id", "field", "raw_value"]
    )
    check("all_audit_joins_and_identity_checks", not MISMATCHES, len(MISMATCHES))
    check("all_nonmissing_numeric_annotations_parse", not INVALID, len(INVALID))
    subjects, provenance = subject_rows(samples)
    write_tsv("subject_level_baseline_all_datasets.tsv", subjects)
    write_tsv(
        "subject_level_baseline_longitudinal_trials.tsv",
        [s for s in subjects if s["dataset"] != "GSE88884"],
    )
    write_tsv("subject_field_provenance.tsv", provenance)
    write_tsv(
        "within_subject_annotation_conflicts.tsv",
        CONFLICTS,
        ["dataset", "trial_or_batch", "subject_id", "field", "values", "sample_values"],
    )
    write_tsv(
        "within_subject_partial_missingness.tsv",
        PARTIAL,
        [
            "dataset",
            "trial_or_batch",
            "subject_id",
            "field",
            "n_samples",
            "n_nonmissing_samples",
            "selected_source",
        ],
    )
    check("within_subject_baseline_annotation_consistency", not CONFLICTS, len(CONFLICTS))
    numeric, categorical, miss, flow, readable = describe(subjects)
    write_tsv("baseline_numeric_by_trial_arm_and_analysis_set.tsv", numeric)
    write_tsv("baseline_categorical_by_trial_arm_and_analysis_set.tsv", categorical)
    write_tsv("subject_missingness_by_trial_arm_and_analysis_set.tsv", miss)
    write_tsv("molecular_sample_flow_recovered.tsv", flow)
    write_tsv("baseline_by_treatment_arm_readable.tsv", readable)
    validate_previous_flow(flow, subjects)
    links, differences = cross_dataset_links(samples, subjects)
    write_tsv("gse88884_to_longitudinal_subject_linkage.tsv", links)
    write_tsv(
        "cross_dataset_clinical_discrepancies.tsv",
        differences,
        ["trial", "subject_id", "field", "longitudinal_annotation", "gse88884_annotation"],
    )
    check(
        "matched_cross_dataset_clinical_values_agree",
        not differences,
        len(differences),
        fatal=False,
    )
    write_tsv("SRI_response_and_dose_field_search.tsv", search)
    sample_missing = []
    for dataset in DATASETS:
        for condition in sorted({s["condition"] for s in samples if s["dataset"] == dataset}):
            members = [
                s for s in samples if s["dataset"] == dataset and s["condition"] == condition
            ]
            for field in CLINICAL:
                sample_missing.append(
                    {
                        "dataset": dataset,
                        "condition": condition,
                        "field": field,
                        "n_samples": len(members),
                        "n_nonmissing": sum((s[field] is not None for s in members)),
                        "n_missing": sum((s[field] is None for s in members)),
                    }
                )
    write_tsv("sample_level_missingness_by_dataset_condition.tsv", sample_missing)
    cohort_counts = []
    for dataset in DATASETS:
        for trial in sorted({s["trial_or_batch"] for s in subjects if s["dataset"] == dataset}):
            for condition in sorted(
                {
                    s["condition"]
                    for s in subjects
                    if s["dataset"] == dataset and s["trial_or_batch"] == trial
                }
            ):
                members = [
                    s
                    for s in subjects
                    if s["dataset"] == dataset
                    and s["trial_or_batch"] == trial
                    and (s["condition"] == condition)
                ]
                cohort_counts.append(
                    {
                        "dataset": dataset,
                        "trial_or_batch": trial,
                        "condition": condition,
                        "n_unique_subjects": len(members),
                        "n_samples": sum((s["n_molecular_samples"] for s in members)),
                        "n_observed_baseline_samples": sum(
                            (s["baseline_sample_present"] for s in members)
                        ),
                    }
                )
    write_tsv("cohort_counts.tsv", cohort_counts)
    dose_counts = Counter(
        (
            (
                s["dataset"],
                s["trial_or_batch"],
                s["treatment_raw"],
                s["treatment"],
                s["stated_dose_mg"],
                s["stated_interval_weeks"],
            )
            for s in subjects
        )
    )
    write_tsv(
        "treatment_dose_annotation_summary.tsv",
        [
            {
                "dataset": k[0],
                "trial_or_batch": k[1],
                "treatment_raw": k[2],
                "treatment": k[3],
                "stated_dose_mg": k[4],
                "stated_interval_weeks": k[5],
                "n_subjects": n,
            }
            for k, n in dose_counts.items()
        ],
    )
    summary = {
        "datasets": {},
        "raw_characteristic_cells_preserved": len(raw_long),
        "sample_audit_mismatches": len(MISMATCHES),
        "numeric_parse_failures": len(INVALID),
        "within_subject_conflicts": len(CONFLICTS),
        "cross_dataset_clinical_conflicts": len(differences),
        "gse88884_subject_matches": sum((r["gse88884_subject_match"] for r in links)),
        "exact_baseline_CEL_matches": sum((r["baseline_CEL_exact_match"] == 1 for r in links)),
        "SRI5_matched_local_sample_annotation_rows": sum(
            (r["matched_annotation_rows"] for r in search if r["query"] == "SRI5")
        ),
        "limitations": [
            "This is a descriptive analysis of publicly deposited molecular subsets, not all randomized parent-trial participants.",
            "Follow-up samples may carry baseline clinical annotations; this never creates a baseline expression sample or a complete expression pair.",
            "No clinical treatment effect, arm-balance equivalence, response-prediction, SRI5 reconstruction, or inferential test was performed.",
            "An annotation absent in three local caches is not evidence of absence from all other datasets or investigator records.",
            "Reported clinical category codes and numeric values are retained as deposited; assay upper/lower censoring was not inferred.",
        ],
    }
    for dataset in DATASETS:
        sm = [s for s in samples if s["dataset"] == dataset]
        ss = [s for s in subjects if s["dataset"] == dataset]
        summary["datasets"][dataset] = {
            "n_samples": len(sm),
            "n_subjects": len(ss),
            "n_characteristic_fields": sum((r["dataset"] == dataset for r in inventory)),
            "n_observed_baseline_subjects": sum((s["baseline_sample_present"] for s in ss)),
            "n_followup_only_subjects": sum((not s["baseline_sample_present"] for s in ss)),
            "subject_missing_clinical": {f: sum((s[f] is None for s in ss)) for f in CLINICAL},
            "clinical_fields_all_nonmissing_subjects": sum(
                (all((s[f] is not None for f in CLINICAL)) for s in ss)
            ),
        }
    for source in SOURCES.values():
        source["sha256_after"] = sha256(Path(source["path"]))
        source["unchanged"] = source["sha256_before"] == source["sha256_after"]
    write_tsv("source_manifest_sha256.tsv", SOURCES.values())
    check("all_sources_unchanged", all((r["unchanged"] for r in SOURCES.values())), len(SOURCES))
    summary["validation_failed_checks"] = [c for c in CHECKS if c["status"] != "PASS"]
    write_tsv("validation_checks.tsv", CHECKS)
    (OUT / "results_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    notes = [
        "临床基线恢复与有限描述性分析（仅本地缓存）",
        "",
        "纠正：原始GSE88884/85/86样本注释含年龄、性别、种族、地区、基线SLEDAI、anti-dsDNA和C3/C4。旧sample_audit未导出这些字段，不等于原始资料不存在。",
        f"原共享pipeline {shared} 第85–123行只读取accession/title/source/supplementary并推导subject/visit/arm；未解析characteristics。",
        "全部characteristics原样保存在sample_characteristics_all_raw_long.tsv；characteristic_field_mapping.tsv保存原字段名/源行/标准化映射。",
        "以GSM精确联结既有sample_audit；按dataset + trial + subject建一人一行。各个baseline-labelled字段先核跨访视一致；有明确baseline且非缺失时优先取该值，否则仅取一致的随访样本所携带基线注释。值冲突不平均。",
        "subject_field_provenance.tsv逐字段记录选择来源。baseline_sample_present只表示该纵向表达数据集中真实baseline样本，随访注释不增加配对数。",
        "基线描述按试验、治疗臂分列四种分子分析集：所有可用纵向分子供者、实际baseline分子供者、baseline/week16完整配对、baseline/week52完整配对。",
        "数值报告非缺失n/缺失n、均值/样本SD、中位数/线性插值Q1-Q3、范围；类别报告计数/非缺失百分比及全组百分比。无P值、无SMD、无治疗效应检验。",
        "各治疗臂原始注释明确LY 120 mg Q2W/Q4W；剂量与间隔仅按该字符串记录，不补充未注明的负荷剂量或其他药物用法。Placebo剂量留空。",
        "GSE88884与纵向队列按trial + subject联结，另核对baseline CEL。跨数据集注释核对结果保留，未将GSE88884单独的baseline表达样本加入纵向模型或配对。",
        "GSE88884是更广泛的公开分子baseline样本池，仍不等同于父试验全部随机化患者；目前不能用此描述表宣称分子子集代表父试验全体。",
        "SRI/SRI4/SRI5、response、BILAG/PGA等检索范围限定三套本地原始Sample注释，详见SRI_response_and_dose_field_search.tsv。基线SLEDAI不能重建SRI5疗效标签；未进行网络检索。",
        "原pipeline使用baseline/week16/week52全部可用样本的subject随机截距混合模型，paired-change用baseline及目标访视完整配对。本次恢复的访视/治疗臂及配对身份与既有trial reporting输出精确核对。",
        "读取的源文件前后SHA256一致性记录在source_manifest_sha256.tsv；脚本先保存后执行；仅使用Python标准库。",
        "",
        json.dumps(summary, indent=2, ensure_ascii=False),
    ]
    (OUT / "READ_ME_结果与限制.txt").write_text("\n".join(notes) + "\n", encoding="utf-8")
    logging.info("COMPLETE %s", json.dumps(summary, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
