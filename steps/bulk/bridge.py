import os
import csv
import shutil
import sys
from pathlib import Path

root = Path(os.environ["SLE_RUN_ROOT"]).resolve()
pairs = {
    "hta": [
        (
            "SLE/GSE88884_pipeline/results/00_input_audit/sample_audit.tsv",
            "SLE_v1/pipeline/3.2/data/metadata/GSE88884_sample_audit.tsv",
        ),
        (
            "SLE/GSE88884_pipeline/data/raw/downloads/GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz",
            "SLE_v1/pipeline/3.2/data/raw/GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz",
        ),
        *[
            (
                f"SLE/GSE{number}_pipeline/data/processed/analysis_results.rds",
                f"SLE_v1/pipeline/3.3/data/processed/GSE{number}_analysis_results.rds",
            )
            for number in ["88884", "88885", "88886"]
        ],
        *[
            (
                f"SLE/GSE{number}_pipeline/results/02_program_scores/sample_program_scores.tsv",
                f"SLE_v1/pipeline/3.4/data/processed/ILLUMINATE{trial}_sample_program_scores.tsv",
            )
            for number, trial in [("88885", 1), ("88886", 2)]
        ],
    ],
    "grouping": [
        (
            "SLE_v1/results/3.3/01_grouping/frozen_grouping_model.rds",
            "SLE_v1/pipeline/3.4/data/processed/frozen_grouping_model.rds",
        )
    ],
    "cross_cohort": [
        (
            "SLE_v1/results/3.2/03_preservation/cross_cohort_module_membership.tsv",
            "SLE/SLE_mainline_research/02_cross_cohort_module_preservation/results/01_preservation/cross_cohort_module_membership.tsv",
        ),
        (
            "SLE_v1/results/3.2/03_preservation/cross_cohort_module_membership.tsv",
            "SLE_v1/pipeline/3.5/data/processed/cross_cohort_module_membership.tsv",
        ),
    ],
}
for source_relative, target_relative in pairs[sys.argv[1]]:
    source = root / source_relative
    target = root / target_relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)
if sys.argv[1] == "hta":
    coverage_rows = []
    columns = None
    for number in ["88884", "88885", "88886"]:
        path = root / f"SLE/GSE{number}_pipeline/results/02_program_scores/program_coverage.tsv"
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            current_columns = reader.fieldnames + ["dataset"]
            if columns is not None and columns != current_columns:
                raise ValueError("Inconsistent platform coverage schemas")
            columns = current_columns
            rows = [
                {**row, "dataset": "GSE" + number}
                for row in reader
                if row["analysis_role"] == "primary"
            ]
            if len(rows) != 3:
                raise ValueError("Each bulk dataset must contain three frozen primary programs")
            coverage_rows.extend(rows)
    target = root / "SLE_R1/checks/r1_revision_20261116/supplement/frozen_bulk_program_coverage.tsv"
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(coverage_rows)
