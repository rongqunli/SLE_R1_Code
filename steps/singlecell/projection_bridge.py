from pathlib import Path
import csv
import os
import shutil

root = Path(os.environ["SLE_RUN_ROOT"])
source = root / "SLE/SLE_mainline_research/05_single_cell_module_projection/results/01_projection"
target = root / "SLE_v1/pipeline/3.5/data/processed"
target.mkdir(parents=True, exist_ok=True)
for name in [
    "donor_celltype_program_scores.tsv",
    "donor_celltype_proportions.tsv",
    "frozen_program_gene_matching.tsv",
]:
    shutil.copy2(source / name, target / name)
shutil.copy2(
    root / "SLE/GSE135779_pipeline/results/05_pseudobulk/pseudobulk_all_genes.rds",
    target / "pseudobulk_all_genes.rds",
)
shutil.copy2(
    root
    / "SLE/SLE_mainline_research/02_cross_cohort_module_preservation/results/01_preservation/cross_cohort_module_membership.tsv",
    target / "cross_cohort_module_membership.tsv",
)
with (target / "frozen_program_gene_matching.tsv").open(encoding="utf-8-sig", newline="") as handle:
    program_rows = list(csv.DictReader(handle, delimiter="\t"))
with (target / "cross_cohort_module_membership.tsv").open(
    encoding="utf-8-sig", newline=""
) as handle:
    membership_rows = list(csv.DictReader(handle, delimiter="\t"))
members = {}
for row in membership_rows:
    members.setdefault(row["module"], set()).add(row["gene"])
expected = {"magenta": 157, "cyan": 35, "greenyellow": 64}
if any((len(members.get(module, set())) != count for module, count in expected.items())):
    raise ValueError("The module membership does not match the frozen pediatric programs")
columns = [
    "source_module",
    "program",
    "requested_genes",
    "matched_genes",
    "match_fraction",
    "passes_scoring_gate",
    "matched_genes_in_membership",
]
audit_rows = []
for row in program_rows:
    matched = [gene for gene in row["matched_gene_list"].split(";") if gene]
    result = {key: row[key] for key in columns[:-1]}
    result["matched_genes_in_membership"] = sum(
        (gene in members[row["source_module"]] for gene in matched)
    )
    if result["matched_genes_in_membership"] != int(row["matched_genes"]):
        raise ValueError("The scoring genes do not match the frozen program membership")
    audit_rows.append(result)
audit = root / "SLE_v1/results/3.5/00_input_audit/program_definition_audit.tsv"
audit.parent.mkdir(parents=True, exist_ok=True)
with audit.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(audit_rows)
