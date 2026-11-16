from pathlib import Path
import csv
import hashlib
import json
import os
import shutil

ROOT = Path(os.environ["SLE_RUN_ROOT"]).resolve()
OUT = ROOT / "SLE_R1/supplementary_data"
MAPPING = [
    (
        "S1/trial_transcriptomic_flow.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/trial_transcriptomic_flow.tsv",
    ),
    (
        "S1/trial_flow_parent_and_molecular.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/trial_flow_parent_and_molecular.tsv",
    ),
    (
        "S1/parent_trial_randomization.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/parent_trial_randomization.tsv",
    ),
    (
        "S2/ancova_coordinate_effects_trial_specific.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_coordinate_effects_trial_specific.tsv",
    ),
    (
        "S2/ancova_coordinate_effects_fixed_effect_meta.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_coordinate_effects_fixed_effect_meta.tsv",
    ),
    (
        "S2/all_available_mixed_model_audit.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/all_available_mixed_model_audit.tsv",
    ),
    (
        "S2/baseline_coordinate_balance_by_arm.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/baseline_coordinate_balance_by_arm.tsv",
    ),
    (
        "S2/baseline_coordinate_pairwise_smd.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/baseline_coordinate_pairwise_smd.tsv",
    ),
    (
        "S2/coordinate_standardization.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/coordinate_standardization.tsv",
    ),
    (
        "S2/trial_mixed_model_coefficients.tsv",
        "SLE_v1/results/3.4/01_overall/trial_mixed_model_coefficients.tsv",
    ),
    ("S2/paired_change_summary.tsv", "SLE_v1/results/3.4/01_overall/paired_change_summary.tsv"),
    (
        "S2/ILLUMINATE-1_mixed_model_results.tsv",
        "SLE_v1/results/3.4/01_overall/ILLUMINATE-1/mixed_model_results.tsv",
    ),
    (
        "S2/ILLUMINATE-2_mixed_model_results.tsv",
        "SLE_v1/results/3.4/01_overall/ILLUMINATE-2/mixed_model_results.tsv",
    ),
    (
        "S2/cross_trial_fixed_effect_meta.tsv",
        "SLE_v1/results/3.4/03_cross_trial/cross_trial_fixed_effect_meta.tsv",
    ),
    (
        "S2/cross_trial_pharmacodynamic_meta.tsv",
        "SLE_v1/results/3.4/03_cross_trial/cross_trial_pharmacodynamic_meta.tsv",
    ),
    ("S2/comparison_families.tsv", "SLE_v1/results/3.4/03_cross_trial/comparison_families.tsv"),
    (
        "S3/vector_components_bootstrap.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/vector_components_bootstrap.tsv",
    ),
    (
        "S4/ancova_bcell_baseline_interactions_trial_specific.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_bcell_baseline_interactions_trial_specific.tsv",
    ),
    (
        "S4/ancova_bcell_baseline_interactions_fixed_effect_meta.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_bcell_baseline_interactions_fixed_effect_meta.tsv",
    ),
    (
        "S5/cohort_counts.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/clinical_baseline/cohort_counts.tsv",
    ),
    (
        "S5/molecular_sample_flow_recovered.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/clinical_baseline/molecular_sample_flow_recovered.tsv",
    ),
    (
        "S5/trial_baseline_overlap_audit.tsv",
        "SLE_v1/results/3.3/01_grouping/trial_baseline_overlap_audit.tsv",
    ),
    ("S5/trial_overlap_summary.tsv", "SLE_v1/results/3.3/02_summary/trial_overlap_summary.tsv"),
    (
        "S5/complete_case_flow.tsv",
        "SLE_v1/results/reviewer_response/02_clinical_landscape_validation/complete_case_flow.tsv",
    ),
    (
        "S5/model_definitions.tsv",
        "SLE_v1/results/reviewer_response/02_clinical_landscape_validation/model_definitions.tsv",
    ),
    (
        "S5/analysis_set_definitions.tsv",
        "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/analysis_set_definitions.tsv",
    ),
    (
        "S6/module_selection_audit.tsv",
        "SLE_v1/results/3.1/07_module_selection_audit/Table_S2_module_selection_audit.tsv",
    ),
    ("S6/module_stability.tsv", "SLE_v1/results/3.1/06_wgcna/module_stability.tsv"),
    (
        "S6/six_modules_membership_hubs.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/six_modules_membership_hubs.tsv",
    ),
    (
        "S6/six_module_GO_BP_enrichment.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/six_module_GO_BP_enrichment.tsv",
    ),
    (
        "S6/six_module_top_GO_terms.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/six_module_top_GO_terms.tsv",
    ),
    (
        "S6/plasmablast_marker_locations.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/plasmablast_marker_locations.tsv",
    ),
    (
        "S6/full_56_module_trait_FDR.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/full_56_module_trait_FDR.tsv",
    ),
    (
        "S6/plasmablast_proxy_clinical_sensitivity.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/plasmablast_proxy_clinical_sensitivity.tsv",
    ),
    (
        "S6/selected_modules_GO_KEGG_overrepresentation.tsv",
        "SLE_v1/results/3.1/07_wgcna_functional_annotation/selected_modules_GO_KEGG_overrepresentation.tsv",
    ),
    (
        "S6/selected_modules_kME_ranked_GO_KEGG_GSEA.tsv",
        "SLE_v1/results/3.1/07_wgcna_functional_annotation/selected_modules_kME_ranked_GO_KEGG_GSEA.tsv",
    ),
    (
        "S6/selected_annotation_terms_for_Figure1.tsv",
        "SLE_v1/results/3.1/07_wgcna_functional_annotation/selected_annotation_terms_for_Figure1.tsv",
    ),
    (
        "S7/frozen_bulk_program_coverage.tsv",
        "SLE_R1/checks/r1_revision_20261116/supplement/frozen_bulk_program_coverage.tsv",
    ),
    (
        "S7/pediatric_single_cell_program_coverage.tsv",
        "SLE_v1/results/3.5/00_input_audit/program_definition_audit.tsv",
    ),
    (
        "S7/frozen_program_membership.tsv",
        "SLE_R1/checks/r1_revision_20261116/supplement/frozen_program_membership.tsv",
    ),
    (
        "S10/Supplementary_Table_S10.tsv",
        "SLE_R1/Try/R2_13_14_15_single_cell/QC_10_15_20_donor_effects.tsv",
    ),
    ("S11/all_tf_results.tsv", "SLE_R1/checks/tf_recheck_20261115/all_tf_results.tsv"),
    ("S11/tf_eligibility.tsv", "SLE_R1/checks/tf_recheck_20261115/tf_eligibility.tsv"),
    ("S11/ifn_gene_audit.tsv", "SLE_R1/checks/tf_recheck_20261115/ifn_gene_audit.tsv"),
    (
        "S11/four_tf_primary_summary.tsv",
        "SLE_R1/checks/tf_recheck_20261115/four_tf_primary_summary.tsv",
    ),
    (
        "S11/four_tf_model_comparison.tsv",
        "SLE_R1/checks/tf_recheck_20261115/four_tf_model_comparison.tsv",
    ),
    ("S11/leave_one_donor_out.tsv", "SLE_R1/checks/tf_recheck_20261115/leave_one_donor_out.tsv"),
    ("S11/bootstrap_audit.tsv", "SLE_R1/checks/tf_recheck_20261115/bootstrap_audit.tsv"),
    (
        "S11/irf7_stat1_difference_summary.tsv",
        "SLE_R1/checks/tf_recheck_20261115/irf7_stat1_difference_summary.tsv",
    ),
    ("S11/model_designs.tsv", "SLE_R1/checks/tf_recheck_20261115/model_designs.tsv"),
    ("S11/candidate_regulator_resource.txt", "SLE_R1/Try/hs_hgnc_tfs.txt"),
    (
        "S12/gse121239_patient_means_corrected.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/association_results/gse121239_patient_means_corrected.tsv",
    ),
    (
        "S12/gse121239_longitudinal_corrected.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/association_results/gse121239_longitudinal_corrected.tsv",
    ),
    (
        "S12/CCL5_paired_sensitivity_differences.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/association_results/CCL5_paired_sensitivity_differences.tsv",
    ),
    (
        "S12/program_gene_overlap_audit.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/association_results/program_gene_overlap_audit.tsv",
    ),
    (
        "S12/gse135779_crosscell_all121.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/association_results/gse135779_crosscell_all121.tsv",
    ),
    (
        "S12/gse135779_systemic_IFN_audit.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/association_results/gse135779_systemic_IFN_audit.tsv",
    ),
    (
        "S12/program_membership.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/gzmk_inputs/program_membership.tsv",
    ),
    (
        "S12/program_coverage.tsv",
        "SLE_R1/checks/targeted_recheck_20261115/gzmk_inputs/program_coverage.tsv",
    ),
    (
        "S13/participant_split_summary.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/participant_split_summary.tsv",
    ),
    (
        "S13/participant_split_results.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/participant_split_results.tsv",
    ),
    (
        "S13/selected_axis_heldout_tests.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/selected_axis_heldout_tests.tsv",
    ),
    ("S14/benchmark_gene_coverage.tsv", "SLE_R1/Try/R2_05_08_adult/benchmark_gene_coverage.tsv"),
    ("S14/cv_summary.tsv", "SLE_R1/Try/R2_05_08_adult/cv_summary.tsv"),
    ("S14/repeated_cv_metrics.tsv", "SLE_R1/Try/R2_05_08_adult/repeated_cv_metrics.tsv"),
    ("S14/heldout_predictions.tsv", "SLE_R1/Try/R2_05_08_adult/heldout_predictions.tsv"),
    (
        "S14/heldout_metrics_trial_separate.tsv",
        "SLE_R1/Try/R2_05_08_adult/heldout_metrics_trial_separate.tsv",
    ),
    (
        "S14/heldout_benchmark_differences.tsv",
        "SLE_R1/Try/R2_05_08_adult/heldout_benchmark_differences.tsv",
    ),
    (
        "S15/medication_field_values.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/medication_field_values.tsv",
    ),
    (
        "S15/medication_adjusted_module_traits.tsv",
        "SLE_R1/Try/R2_02_03_10_12_modules/medication_adjusted_module_traits.tsv",
    ),
    (
        "S16/parent_vs_molecular_primary_309_147.tsv",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/comparison/parent_vs_molecular_primary_309_147.tsv",
    ),
    (
        "S16/parent_vs_molecular_all_analysis_sets.tsv",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/comparison/parent_vs_molecular_all_analysis_sets.tsv",
    ),
    (
        "S16/analysis_set_denominators.tsv",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/comparison/analysis_set_denominators.tsv",
    ),
    (
        "S16/sex_race_label_mapping.tsv",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/comparison/sex_race_label_mapping.tsv",
    ),
    (
        "S16/geography_original_labels_not_harmonized.tsv",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/comparison/geography_original_labels_not_harmonized.tsv",
    ),
    (
        "S16/field_availability_and_comparability.tsv",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/comparison/field_availability_and_comparability.tsv",
    ),
    (
        "S16/NCT01196091.json",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/public_sources/NCT01196091.json",
    ),
    (
        "S16/NCT01205438.json",
        "SLE_R1/checks/parent_trial_baseline_audit_20261116/public_sources/NCT01205438.json",
    ),
    (
        "S17/joint_fdr_results.tsv",
        "SLE_R1/checks/paperwide_fdr_sensitivity_20261115/unassessed14_recheck/joint_update/joint_fdr_results.tsv",
    ),
    (
        "S17/family_summary.tsv",
        "SLE_R1/checks/paperwide_fdr_sensitivity_20261115/unassessed14_recheck/joint_update/family_summary.tsv",
    ),
    (
        "S17/newly_assessed_14.tsv",
        "SLE_R1/checks/paperwide_fdr_sensitivity_20261115/unassessed14_recheck/joint_update/newly_assessed_14.tsv",
    ),
    (
        "S17/threshold_changes_among_previous_1643.tsv",
        "SLE_R1/checks/paperwide_fdr_sensitivity_20261115/unassessed14_recheck/joint_update/threshold_changes_among_previous_1643.tsv",
    ),
    (
        "S17/heldout_final_results.tsv",
        "SLE_R1/checks/paperwide_fdr_sensitivity_20261115/unassessed14_recheck/heldout_inference/heldout_final_results.tsv",
    ),
    ("S18/adult_soft_threshold.tsv", "SLE_R1/Try/R2_05_08_adult/adult_soft_threshold.tsv"),
    (
        "S18/independent_adult_module_membership.tsv",
        "SLE_R1/Try/R2_05_08_adult/independent_adult_module_membership.tsv",
    ),
    (
        "S18/adult_pediatric_module_overlaps.tsv",
        "SLE_R1/Try/R2_05_08_adult/adult_pediatric_module_overlaps.tsv",
    ),
    ("S18/best_module_overlaps.tsv", "SLE_R1/Try/R2_05_08_adult/best_module_overlaps.tsv"),
    (
        "S18/module_preservation_statistics.tsv",
        "SLE_v1/results/3.2/03_preservation/module_preservation_statistics.tsv",
    ),
    (
        "S18/module_membership_concordance.tsv",
        "SLE_v1/results/3.2/03_preservation/module_membership_concordance.tsv",
    ),
]
MANIFEST = []


def read(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return (reader.fieldnames, list(reader))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def record(target, sources, operation):
    destination = OUT / target
    rows = len(read(destination)[1]) if destination.suffix == ".tsv" else ""
    MANIFEST.append(
        dict(
            table=target.split("/")[0],
            file=target,
            rows=rows,
            sha256=sha(destination),
            source_files=[str(p.relative_to(ROOT)) for p in sources],
            source_sha256=[sha(p) for p in sources],
            operation=operation,
        )
    )


def tf_methods():
    tf = ROOT / "SLE_R1/checks/tf_recheck_20261115"
    paths = [
        tf / "formal_inputs/formal_nk_metadata.tsv",
        tf / "tf_eligibility.tsv",
        tf / "ifn_gene_audit.tsv",
        tf / "four_tf_primary_summary.tsv",
        tf / "irf7_stat1_difference_summary.tsv",
        tf / "model_designs.tsv",
        ROOT / "SLE_R1/Try/hs_hgnc_tfs.txt",
    ]
    metadata = read(paths[0])[1]
    eligibility = read(paths[1])[1]
    ifn = read(paths[2])[1]
    summary = read(paths[3])[1]
    differences = read(paths[4])[1]
    designs = read(paths[5])[1]
    resource_n = len(paths[6].read_text().split())
    first = next((row for row in summary if row["model"] == "all_adjusted"))
    primary = [row for row in ifn if row["primary_score_member"].lower() == "true"]
    candidates = len(eligibility)
    cells = sum((int(row["n_cells"]) for row in metadata))
    paragraphs = [
        f"Raw counts from {cells:,} NK cells were aggregated to {len(metadata)} pediatric donors. TMM normalization was estimated over genes with nonzero total counts; log2 CPM used prior count 2. The candidate resource contains {resource_n:,} human TF/cofactor genes. Eligibility requires raw CPM >=1 in at least five donors and nonconstant normalized expression.",
        f"The common primary IFN outcome contains {len(primary)} genes after excluding all eligible candidate regulators from the measurable frozen IFN module. Gene-wise centers and scales are estimated over the full pediatric donor set and retained for SLE-only models. The complete requested-candidate eligibility table has {candidates} rows.",
        "Associations are Pearson correlations of rank-scale residuals. Primary covariates are disease, batch, sex, ranked age and ranked log mean raw UMI per NK cell. SLE-only analyses replace disease with ranked continuous SLEDAI and use complete cases. Approximate two-sided partial-correlation t tests and restricted Freedman-Lane residual permutations are reported separately. Permutations use disease-by-batch strata in the full cohort and batch strata in SLE-only analyses, 99,999 draws and plus-one probabilities. BH correction includes the complete eligible candidate family within each model and scoring variant.",
        "Donor bootstrap intervals use 2,000 draws, re-ranking and refitting each draw while retaining normalization, gene-standardization weights and the candidate set. Rank-deficient or unestimable draws are rejected explicitly. Percentile confidence intervals condition on the valid draws. Leave-one-donor-out fits retain all eligible candidates.",
        "Scoring variants are exclude_all_tested_TFs (primary), full_frozen, exclude_self and exclude_four_candidates. The latter removes measured frozen-module members among IRF7, STAT1, STAT2 and IRF9. Resource version, normalization inputs, full model designs, eligibility and resampling diagnostics are saved with the analysis.",
    ]
    for name in ["all_adjusted", "sle_adjusted"]:
        rows = [row for row in summary if row["model"] == name]
        one = rows[0]
        design = next((row for row in designs if row["model"] == name))
        diff = next((row for row in differences if row["model"] == name))
        paragraphs.append(
            f"{name}: n={design['n']}, nuisance design rank={design['nuisance_rank']}, partial df={design['partial_df']}, valid bootstrap draws={one['bootstrap_valid']}, rejected draws={one['bootstrap_rejected']}. IRF7-minus-STAT1 rho difference={float(diff['observed_difference']):.6g}; conditional 95% interval [{float(diff['bootstrap_95_low']):.6g}, {float(diff['bootstrap_95_high']):.6g}]."
        )
    path = OUT / "S11/expanded_methods.txt"
    path.write_text("\n\n".join(paragraphs) + "\n", encoding="utf-8")
    record("S11/expanded_methods.txt", paths, "methods and this-run resampling summary")


def fdr_methods():
    path = (
        ROOT
        / "SLE_R1/checks/paperwide_fdr_sensitivity_20261115/unassessed14_recheck/joint_update/joint_fdr_results.tsv"
    )
    rows = read(path)[1]
    families = {row["family"] for row in rows}
    count = lambda key: sum((float(row[key]) < 0.05 for row in rows))
    changed = sum(
        (
            any((row[key].lower() == "true" for key in row if key.startswith("changed_005_")))
            for row in rows
        )
    )
    paragraphs = [
        f"This supplementary sensitivity includes {len(rows):,} assessed hypotheses in {len(families)} complete families. It preserves the original family-specific results and does not replace them in the main manuscript. Joint BH and BY values are supplied. Gene-level differential expression and annotation enrichment retain their own complete local families. Alternative scores, medication, QC and SLE-only sensitivity analyses are not counted again as new primary hypotheses.",
        "The two baseline B-score/B-cell-abundance tests use two-sided asymptotic Spearman t inference and an independent 99,999-permutation plus-one sensitivity. Twelve held-out comparisons use 99,999 null-centered studentized participant-bootstrap draws with discovery-fitted predictions fixed; pooled analyses resample within trial.",
        "Joint method sensitivities replace every eligible TF approximate-test probability together with its restricted residual-permutation probability; baseline-correlation permutation variants are separate and combined. Complete candidate sets are preserved for TF, adult-module and cross-cell searches.",
        f"The primary joint analysis has {count('q_joint_BH')} BH and {count('q_joint_BY')} BY discoveries at q<0.05. Replacing TF probabilities has {count('q_joint_TFperm_BH')} BH and {count('q_joint_TFperm_BY')} BY discoveries. Completing the fourteen formerly unassessed entries changes threshold classifications in {changed} of the other 1,643 entries across primary and TF-permutation BH/BY comparisons.",
        "The conservative former state is reconstructed from the same current inputs by setting only the fourteen assessed entries to probability 1. No historical inventory or old corrected-results table is used. All source records and runtime input hashes are retained.",
    ]
    (OUT / "S17/expanded_methods.txt").write_text("\n\n".join(paragraphs) + "\n", encoding="utf-8")
    record("S17/expanded_methods.txt", [path], "methods and this-run family summary")


def main():
    missing = [str(ROOT / source) for _, source in MAPPING if not (ROOT / source).is_file()]
    if missing:
        raise FileNotFoundError("Missing this-run supplementary inputs:\n" + "\n".join(missing))
    OUT.mkdir(parents=True, exist_ok=True)
    for target, relative in MAPPING:
        source = ROOT / relative
        destination = OUT / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        if target in [
            "S5/analysis_set_definitions.tsv",
            "S10/Supplementary_Table_S10.tsv",
            "S11/model_designs.tsv",
        ]:
            fields, rows = read(source)
            if target == "S5/analysis_set_definitions.tsv":
                rows = [row for row in rows if "selectivity" not in " ".join(row.values()).lower()]
            elif target == "S10/Supplementary_Table_S10.tsv":
                assert len(rows) == 6
                rows.sort(
                    key=lambda row: (
                        ["NK", "CD14_Monocyte"].index(row["cell_type"]),
                        int(row["cutoff"]),
                    )
                )
            else:
                fields = ["model", "n", "SLE", "HD", "nuisance_rank", "partial_df", "covariates"]
                rows = [{key: row[key] for key in fields} for row in rows]
            write(destination, fields, rows)
            operation = (
                "preserve numerical text; select retained rows/columns or presentation order"
            )
        else:
            shutil.copyfile(source, destination)
            assert sha(source) == sha(destination)
            operation = "copy this-run output or staged reference resource"
        record(target, [source], operation)
    tf_methods()
    fdr_methods()
    assert len(MANIFEST) == 91 and len({row["file"] for row in MANIFEST}) == 91
    assert all(((OUT / row["file"]).is_file() for row in MANIFEST))
    fields = ["table", "file", "rows", "sha256"]
    write(OUT / "file_index.tsv", fields, [{key: row[key] for key in fields} for row in MANIFEST])
    (OUT / "runtime_source_manifest.json").write_text(
        json.dumps(MANIFEST, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
