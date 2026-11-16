from pathlib import Path
import os
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(os.environ["SLE_RUN_ROOT"])
V1 = ROOT / "SLE_v1"
CHECKS = ROOT / "SLE_R1/checks"
FDR = CHECKS / "paperwide_fdr_sensitivity_20261115"
OUT = FDR / "unassessed14_recheck/joint_update"
ROWS = []
SOURCES = {}


def adjust(values, by=False):
    p = np.asarray(values, dtype=float)
    assert np.isfinite(p).all() and ((p > 0) & (p <= 1)).all()
    order = np.argsort(p, kind="stable")
    n = len(p)
    multiplier = np.sum(1 / np.arange(1, n + 1)) if by else 1.0
    ranked = np.minimum.accumulate((p[order] * n * multiplier / np.arange(1, n + 1))[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(1, ranked)
    return q


def collect(
    family,
    relative,
    ids,
    pcol,
    qcol,
    effect,
    expected,
    filters=None,
    tfperm=None,
    baseperm=None,
    method="",
):
    path = ROOT / relative
    SOURCES[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    data = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    data["source_file_line"] = np.arange(2, len(data) + 2)
    for field, accepted in (filters or {}).items():
        data = data[data[field].isin(accepted)]
    assert len(data) == expected, (family, len(data), expected)
    assert not data.duplicated(ids).any(), family
    for _, row in data.iterrows():
        test_id = "|".join((row[x] for x in ids))
        p = float(row[pcol])
        perm_tf = float(row[tfperm]) if tfperm else p
        perm_base = float(row[baseperm]) if baseperm else p
        ROWS.append(
            dict(
                hypothesis_id=family + "::" + test_id,
                family=family,
                test_id=test_id,
                raw_p_text=row[pcol],
                p_for_adjustment=p,
                p_TF_permutation_sensitivity=perm_tf,
                p_baseline_permutation_sensitivity=perm_base,
                p_TF_and_baseline_permutation_sensitivity=perm_tf if tfperm else perm_base,
                original_local_q=float(row[qcol]) if qcol else np.nan,
                original_q_column=qcol or "",
                effect=float(row[effect]) if effect else np.nan,
                source_file=str(path),
                source_file_line=int(row.source_file_line),
                source_p_column=pcol,
                method=method or pcol,
                valid_p=True,
            )
        )


def main():
    r1 = "SLE_R1/Try/"
    results = "SLE_v1/results/"
    rr = results + "reviewer_response/"
    collect(
        "child_module_traits_56",
        r1 + "R2_02_03_10_12_modules/full_56_module_trait_FDR.tsv",
        ["module", "term"],
        "p_value",
        "global_fdr_56",
        "estimate",
        56,
    )
    collect(
        "adult_module_traits_56",
        results + "3.2/03_preservation/adult_module_clinical_associations.tsv",
        ["module", "trait"],
        "p_value",
        "fdr",
        "beta",
        56,
    )
    collect(
        "celltype_IFN_effects_9",
        rr + "04_cd14_nk_unified/all_celltype_unified_ifn_effects.tsv",
        ["cell_type"],
        "p_value",
        "fdr_cell_effects",
        "estimate",
        9,
    )
    collect(
        "multicellular_factors_10",
        results
        + "advanced_program_landscape/02_multicellular_tensor/factor_disease_associations.tsv",
        ["factor"],
        "p_value",
        "fdr",
        "estimate",
        10,
    )
    collect(
        "trial_coordinate_effects_24",
        rr + "08_trial_reporting_ancova_permutation/ancova_coordinate_effects_trial_specific.tsv",
        ["trial", "time", "dose", "coordinate"],
        "p_value",
        "fdr",
        "estimate",
        24,
    )
    collect(
        "trial_baseline_interactions_8",
        rr
        + "08_trial_reporting_ancova_permutation/ancova_bcell_baseline_interactions_trial_specific.tsv",
        ["trial", "time", "dose"],
        "p_value",
        "fdr",
        "estimate",
        8,
    )
    collect(
        "continuous_vs_groups_4",
        rr + "02_clinical_landscape_validation/formal_cv_continuous_vs_four_group.tsv",
        ["trait"],
        "corrected_resampled_p",
        None,
        "mean_cv_r2_difference_continuous_minus_four",
        4,
    )
    collect(
        "adult_network_overlap_48",
        r1 + "R2_05_08_adult/adult_pediatric_module_overlaps.tsv",
        ["pediatric_module", "adult_module"],
        "p_value",
        "fdr_all_overlaps",
        "jaccard",
        48,
    )
    collect(
        "adult_frozen_GSEA_3",
        results + "3.2/02_deg_gsea/fgsea_all_collections.tsv",
        ["pathway"],
        "pval",
        "padj",
        "NES",
        3,
        {"collection": ["FROZEN"]},
    )
    collect(
        "clinical_topology_12",
        results
        + "exploratory_program_landscape/01_pairwise_clinical/program_high_vs_immune_low_clinical_contrasts.tsv",
        ["high_group", "trait"],
        "p_value",
        "fdr",
        "adjusted_standardized_effect",
        12,
        {"trait": ["SLEDAI", "C3", "C4", "anti_dsDNA"]},
    )
    collect(
        "trial_Bcell_fraction_8",
        rr + "07_bcell_abundance/abundance_adjusted_treatment_models.tsv",
        ["trial", "time", "dose"],
        "p_value",
        "fdr",
        "estimate",
        8,
        {"model": ["quantiseq_fraction_change"]},
    )
    checks = "SLE_R1/checks/"
    collect(
        "corrected_TF_1278",
        checks + "tf_recheck_20261115/all_tf_results.tsv",
        ["TF"],
        "p_value",
        "fdr",
        "rho",
        1278,
        {"model": ["all_adjusted"], "score": ["exclude_all_tested_TFs"]},
        tfperm="permutation_p",
    )
    gz = checks + "targeted_recheck_20261115/association_results/"
    collect(
        "GZMK_patient_means_3",
        gz + "gse121239_patient_means_corrected.tsv",
        ["relation"],
        "p_t_approx",
        "q_BH_within_family",
        "rho",
        3,
        {"score_variant": ["noCCL5"]},
    )
    collect(
        "GZMK_within_patient_3",
        gz + "gse121239_longitudinal_corrected.tsv",
        ["relation"],
        "p_wild_cluster",
        "q_BH_within_family",
        "rho",
        3,
        {"score_variant": ["noCCL5"]},
    )
    collect(
        "GZMK_crosscell_121",
        gz + "gse135779_crosscell_all121.tsv",
        ["CD14_program", "GZMK_program"],
        "p_t_approx",
        "q_BH_121",
        "rho",
        121,
    )
    assessed = checks + "paperwide_fdr_sensitivity_20261115/unassessed14_recheck/"
    collect(
        "heldout_performance_12",
        assessed + "heldout_inference/heldout_final_results.tsv",
        ["validation_set", "trait"],
        "p_primary",
        None,
        "delta_r2",
        12,
    )
    collect(
        "baseline_B_fraction_correlation_2",
        assessed + "baseline_correlations/baseline_correlation_final_results.tsv",
        ["trial"],
        "p_primary",
        None,
        "rho",
        2,
        baseperm="p_permutation_plus_one",
    )
    data = (
        pd.DataFrame(ROWS)
        .sort_values(["family", "hypothesis_id"], kind="stable")
        .reset_index(drop=True)
    )
    assert len(data) == 1657 and data.family.nunique() == 17
    assert not data.hypothesis_id.duplicated().any()
    assert not data.duplicated(["source_file", "source_file_line", "source_p_column"]).any()
    configs = [
        ("p_for_adjustment", "q_joint"),
        ("p_TF_permutation_sensitivity", "q_joint_TFperm"),
        ("p_baseline_permutation_sensitivity", "q_joint_baselineperm"),
        ("p_TF_and_baseline_permutation_sensitivity", "q_joint_TF_and_baselineperm"),
    ]
    for pcol, prefix in configs:
        for suffix, by in [("BH", False), ("BY", True)]:
            data[prefix + "_" + suffix] = adjust(data[pcol], by)
    data["q_local_BH_check"] = np.nan
    data["q_local_TFperm_BH_check"] = np.nan
    for _, index in data.groupby("family").groups.items():
        data.loc[index, "q_local_BH_check"] = adjust(data.loc[index, "p_for_adjustment"])
        data.loc[index, "q_local_TFperm_BH_check"] = adjust(
            data.loc[index, "p_TF_permutation_sensitivity"]
        )
    data["joint_family_size"] = len(data)
    data["was_unassessed"] = data.family.isin(
        ["heldout_performance_12", "baseline_B_fraction_correlation_2"]
    )
    data["family_display"] = data.family
    data["evaluation_status"] = np.where(
        data.was_unassessed, "assessed_supplementary_recheck", "assessed"
    )
    for pcol, prefix in configs[:2]:
        previous_p = data[pcol].to_numpy(float).copy()
        previous_p[data.was_unassessed] = 1.0
        for suffix, by in [("BH", False), ("BY", True)]:
            qcol = prefix + "_" + suffix
            data["previous_" + qcol] = adjust(previous_p, by)
            data["changed_005_" + qcol] = ~data.was_unassessed & (
                (data[qcol] < 0.05) != (data["previous_" + qcol] < 0.05)
            )
    OUT.mkdir(parents=True, exist_ok=True)
    data.to_csv(OUT / "joint_fdr_results.tsv", sep="\t", index=False, float_format="%.17g")
    data[data.was_unassessed].to_csv(
        OUT / "newly_assessed_14.tsv", sep="\t", index=False, float_format="%.17g"
    )
    changed = data[[x for x in data if x.startswith("changed_005_")]].any(axis=1)
    data[changed].to_csv(
        OUT / "threshold_changes_among_previous_1643.tsv",
        sep="\t",
        index=False,
        float_format="%.17g",
    )
    summaries = []
    for family, block in data.groupby("family", sort=True):
        row = dict(
            family=family,
            family_display=family,
            n_registered=len(block),
            n_assessed=len(block),
            n_unassessed=0,
            n_newly_assessed=int(block.was_unassessed.sum()),
            n_q_local_BH_below_005=int((block.q_local_BH_check < 0.05).sum()),
        )
        for _, prefix in configs:
            for suffix in ["BH", "BY"]:
                col = prefix + "_" + suffix
                row["n_" + col + "_below_005"] = int((block[col] < 0.05).sum())
        summaries.append(row)
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "family_summary.tsv", sep="\t", index=False)
    (OUT / "source_sha256.json").write_text(json.dumps(SOURCES, indent=2), encoding="utf-8")
    membership = pd.read_csv(
        ROOT / "SLE_R1/Try/R2_02_03_10_12_modules/six_modules_membership_hubs.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    selected = membership[membership.module.isin(["magenta", "greenyellow", "cyan"])].copy()
    assert selected.groupby("module").size().to_dict() == {
        "cyan": 35,
        "greenyellow": 64,
        "magenta": 157,
    }
    supplement = CHECKS / "r1_revision_20261116/supplement"
    supplement.mkdir(parents=True, exist_ok=True)
    selected.to_csv(
        supplement / "frozen_program_membership.tsv", sep="\t", index=False, float_format="%.17g"
    )


if __name__ == "__main__":
    main()
