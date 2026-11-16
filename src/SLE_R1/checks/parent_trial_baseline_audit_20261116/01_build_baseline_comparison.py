import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
from pathlib import Path
import csv
import hashlib
import json
import statistics

ROOT = Path(__file__).resolve().parent
SUBJECTS = (
    ROOT.parent
    / "targeted_recheck_20261115/clinical_baseline/subject_level_baseline_longitudinal_trials.tsv"
)
REGISTRY = ROOT / "public_sources/registered_baseline_long.csv"
MAIN = Path(__file__)
OUT = ROOT / "comparison"


def read(path, delim="\t"):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(name, rows):
    with (OUT / name).open("w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


inputs = [
    SUBJECTS,
    REGISTRY,
    MAIN,
    ROOT / "public_sources/NCT01196091.json",
    ROOT / "public_sources/NCT01205438.json",
]
before = {str(p): sha(p) for p in inputs}
subjects = read(SUBJECTS)
registered = read(REGISTRY, ",")
assert len(subjects) == len({r["subject_key"] for r in subjects}) == 456
assert all((r["n_conflicting_fields"] == "0" for r in subjects))
OUT.mkdir(exist_ok=True)
trials = {"ILLUMINATE-1": ("NCT01196091", 1164, 1138), "ILLUMINATE-2": ("NCT01205438", 1124, 1124)}
sets = {
    "all_available_molecular": None,
    "observed_baseline_molecular": "baseline_sample_present",
    "complete_baseline_week16": "complete_baseline_week16",
    "complete_baseline_week52": "complete_baseline_week52",
}
expected = {"ILLUMINATE-1": [309, 302, 296, 297], "ILLUMINATE-2": [147, 139, 128, 137]}
race_map = {
    "American Indian or Alaska Native": "amer_ind_ak_native",
    "Asian": "asian",
    "Native Hawaiian or Other Pacific Islander": None,
    "Black or African American": "afr_amer",
    "White": "white",
    "More than one race": "multiple",
    "Unknown or Not Reported": None,
}
numeric = {
    "age_at_baseline": ("Age, Continuous", "Age, years", "years", "same_label_and_unit"),
    "sledai_at_baseline": (
        "Safety of Estrogens in Lupus Erythematosus National Assessment (SELENA-SLEDAI) Score",
        "GEO SLEDAI / registry SELENA-SLEDAI",
        "GEO SLEDAI score",
        "side_by_side_definition_not_fully_verified",
    ),
    "anti_dsdna_iu_at_baseline": (
        "Anti-dsDNA Antibody Level",
        "Anti-dsDNA (original source units)",
        "GEO IU; denominator not explicit",
        "side_by_side_unit_not_fully_verified",
    ),
}


def parent_row(nct, measure, category=""):
    hits = [
        r
        for r in registered
        if r["nct"] == nct
        and r["group_title"] == "Total"
        and (r["measure"] == measure)
        and (r["category"] == category)
    ]
    assert len(hits) == 1, (nct, measure, category, len(hits))
    return hits[0]


comparisons, counts, region_rows, mappings = ([], [], [], [])
for trial, (nct, randomized, baseline_n) in trials.items():
    all_trial = [r for r in subjects if r["trial_or_batch"] == trial]
    assert len(all_trial) == expected[trial][0]
    assert set((r["race"] for r in all_trial)) <= set((v for v in race_map.values() if v))
    assert set((r["sex"] for r in all_trial)) <= {"Female", "Male"}
    for i, (analysis_set, flag) in enumerate(sets.items()):
        rows = [r for r in all_trial if flag is None or r[flag] == "1"]
        assert len(rows) == expected[trial][i]
        counts.append(
            dict(
                trial=trial,
                analysis_set=analysis_set,
                molecular_n=len(rows),
                parent_randomized_n=randomized,
                registered_baseline_n=baseline_n,
                percent_of_randomized=100 * len(rows) / randomized,
                no_baseline_expression_n=sum((r["baseline_sample_present"] == "0" for r in rows)),
            )
        )
        for field, (measure, label, unit, comparability) in numeric.items():
            p = parent_row(nct, measure)
            assert p["param_type"] == "MEAN" and p["dispersion_type"] == "STANDARD_DEVIATION"
            assert int(p["default_baseline_n"]) == baseline_n
            values = [float(r[field]) for r in rows if r[field] != ""]
            avg, sd = (statistics.mean(values), statistics.stdev(values))
            comparisons.append(
                dict(
                    trial=trial,
                    analysis_set=analysis_set,
                    characteristic=label,
                    source_field=field,
                    category="",
                    kind="continuous",
                    molecular_n=len(rows),
                    molecular_nonmissing_n=len(values),
                    molecular_count="",
                    molecular_mean=avg,
                    molecular_sd=sd,
                    molecular_percent="",
                    molecular_units=unit,
                    parent_baseline_n=baseline_n,
                    parent_measure_n=p["explicit_measure_n"],
                    parent_count="",
                    parent_mean=p["value"],
                    parent_sd=p["spread"],
                    parent_percent="",
                    parent_units=p["units"],
                    difference=avg - float(p["value"]) if field == "age_at_baseline" else "",
                    difference_units="years" if field == "age_at_baseline" else "",
                    comparability=comparability,
                    parent_population=p["population_description"],
                    source_url=p["source_url"],
                    json_pointer=p["json_pointer"],
                )
            )
        for field, measure, levels in [
            ("sex", "Sex: Female, Male", {"Female": "Female", "Male": "Male"}),
            ("race", "Race (NIH/OMB)", race_map),
        ]:
            total_parent = 0
            for label, raw_code in levels.items():
                p = parent_row(nct, measure, label)
                assert p["param_type"] == "COUNT_OF_PARTICIPANTS"
                assert int(p["default_baseline_n"]) == baseline_n
                k = sum((r[field] == raw_code for r in rows)) if raw_code else 0
                kp = int(p["value"])
                total_parent += kp
                prop, parent_prop = (100 * k / len(rows), 100 * kp / baseline_n)
                comparisons.append(
                    dict(
                        trial=trial,
                        analysis_set=analysis_set,
                        characteristic=field,
                        source_field=field,
                        category=label,
                        kind="categorical",
                        molecular_n=len(rows),
                        molecular_nonmissing_n=len(rows),
                        molecular_count=k,
                        molecular_mean="",
                        molecular_sd="",
                        molecular_percent=prop,
                        molecular_units="participants / percent",
                        parent_baseline_n=baseline_n,
                        parent_measure_n=p["explicit_measure_n"],
                        parent_count=kp,
                        parent_mean="",
                        parent_sd="",
                        parent_percent=parent_prop,
                        parent_units=p["units"],
                        difference=prop - parent_prop,
                        difference_units="percentage points",
                        comparability="explicit_label_mapping",
                        parent_population=p["population_description"],
                        source_url=p["source_url"],
                        json_pointer=p["json_pointer"],
                    )
                )
                if i == 0:
                    mappings.append(
                        dict(
                            trial=trial,
                            field=field,
                            registry_category=label,
                            GEO_code=raw_code or "(none observed; all codes accounted for)",
                            molecular_nonmissing_n=len(rows),
                            registry_n=baseline_n,
                        )
                    )
            assert total_parent == baseline_n
        for code in sorted({r["region"] for r in rows}):
            k = sum((r["region"] == code for r in rows))
            region_rows.append(
                dict(
                    trial=trial,
                    analysis_set=analysis_set,
                    source="GEO",
                    geography_original=code,
                    count=k,
                    denominator=len(rows),
                    percent=100 * k / len(rows),
                    note="Original region code; no forced country-to-region crosswalk",
                    source_url="",
                    json_pointer="",
                )
            )
    for p in registered:
        if (
            p["nct"] == nct
            and p["group_title"] == "Total"
            and (p["measure"] == "Region of Enrollment")
        ):
            region_rows.append(
                dict(
                    trial=trial,
                    analysis_set="registered_baseline",
                    source="ClinicalTrials.gov",
                    geography_original=p["class_title"],
                    count=int(p["value"]),
                    denominator=baseline_n,
                    percent=100 * int(p["value"]) / baseline_n,
                    note="Original registry country; not mapped to GEO region",
                    source_url=p["source_url"],
                    json_pointer=p["json_pointer"],
                )
            )
write_tsv("parent_vs_molecular_all_analysis_sets.tsv", comparisons)
write_tsv(
    "parent_vs_molecular_primary_309_147.tsv",
    [r for r in comparisons if r["analysis_set"] == "all_available_molecular"],
)
write_tsv("analysis_set_denominators.tsv", counts)
write_tsv("sex_race_label_mapping.tsv", mappings)
write_tsv("geography_original_labels_not_harmonized.tsv", region_rows)
availability = [
    ("Age", "yes", "yes", "direct descriptive comparison; age in years"),
    ("Sex", "yes", "yes", "direct label mapping"),
    ("Race", "yes", "yes", "explicit category mapping; no missing molecular race"),
    (
        "Disease activity",
        "GEO SLEDAI",
        "SELENA-SLEDAI",
        "side by side; exact instrument label not independently confirmed",
    ),
    (
        "Anti-dsDNA continuous",
        "GEO IU",
        "IU/mL",
        "side by side; denominator and assay equivalence not confirmed",
    ),
    ("Region", "broad codes", "country", "separate original categories; no forced mapping"),
    (
        "C3/C4",
        "yes, g/L",
        "not in registered baseline table",
        "not comparable from retrieved registry baseline",
    ),
    (
        "Anti-dsDNA positivity",
        "yes",
        "not in registered baseline table",
        "not inferred from mean anti-dsDNA",
    ),
    ("Disease duration", "not recovered", "yes", "no individual molecular field to compare"),
    (
        "Background steroid use/dose",
        "not recovered",
        "not in registered baseline table",
        "tabalumab120mg is not steroid dose",
    ),
    (
        "Ethnicity / PGA / BILAG / quality of life",
        "not recovered",
        "yes",
        "not substituted using race or SLEDAI",
    ),
]
write_tsv(
    "field_availability_and_comparability.tsv",
    [dict(field=a, molecular=b, parent=c, action=d) for a, b, c, d in availability],
)
lines = [
    "分子亚组与登记母试验基线：描述性对照",
    "I1: 分子309，登记基线ITT1138（随机化1164）；I2: 分子147，登记基线1124。",
    "连续变量为mean (SD)；分类变量为n/N (%)。差值不含显著性检验。",
    "登记连续变量未单独报告测量有效N；表级N不等于已证实每项完整。",
    "",
]
for trial in trials:
    lines.append(trial)
    for r in comparisons:
        if r["trial"] != trial or r["analysis_set"] != "all_available_molecular":
            continue
        label = r["characteristic"] + (" / " + r["category"] if r["category"] else "")
        if r["kind"] == "continuous":
            mol = f"{r['molecular_mean']:.2f} ({r['molecular_sd']:.2f}); n={r['molecular_nonmissing_n']}"
            par = f"{float(r['parent_mean']):.2f} ({float(r['parent_sd']):.2f})"
        else:
            mol = f"{r['molecular_count']}/{r['molecular_n']} ({r['molecular_percent']:.2f}%)"
            par = f"{r['parent_count']}/{r['parent_baseline_n']} ({r['parent_percent']:.2f}%)"
        diff = (
            f"; difference={r['difference']:+.2f} {r['difference_units']}"
            if r["difference"] != ""
            else "; 原标签并列，未计算差值"
        )
        lines.append(f"{label}: molecular {mol}; parent {par}{diff}")
    lines.append("")
(OUT / "母试验基线对照_可读版.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
after = {str(p): sha(p) for p in inputs}
assert before == after
validation = dict(
    subjects=456,
    comparison_rows=len(comparisons),
    primary_rows=sum((r["analysis_set"] == "all_available_molecular" for r in comparisons)),
    analysis_set_rows=len(counts),
    inputs_unchanged=before == after,
    input_sha256=before,
    no_tests_or_models_fitted=True,
    formal_manuscript_unchanged=True,
)
(OUT / "build_validation.json").write_text(
    json.dumps(validation, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps({k: v for k, v in validation.items() if k != "input_sha256"}, ensure_ascii=False))
