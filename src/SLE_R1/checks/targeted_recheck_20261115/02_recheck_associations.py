import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
from pathlib import Path
import argparse
import hashlib
import json
import platform
import sys
import time
import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent
INP = ROOT / "gzmk_inputs"
OUT = ROOT / "association_results"
COLS = [
    "CD8_GZMK_STATE",
    "CYTOTOXICITY_SHARED",
    "GREENYELLOW_B_CELL_COMPLEMENT",
    "CYAN_NEUTROPHIL_IMMATURE",
    "MAGENTA_IFN_ANTIVIRAL",
    "CD8_GZMK_STATE_noCCL5",
    "CYTOTOXICITY_SHARED_noCCL5",
]
MODELS = [
    (v, label, x, y)
    for v, gx, y in [("noCCL5", 5, 6), ("full", 0, 1)]
    for label, x in [("GZMK", gx), ("B_cell", 2), ("Neutrophil", 3)]
]


def save(df, name):
    df.to_csv(OUT / name, sep="\t", index=False, float_format="%.17g")


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def bh_fixed(p):
    p = np.asarray(p, float)
    order = np.argsort(np.where(np.isfinite(p), p, 1.0), kind="stable")
    ordered = np.where(np.isfinite(p[order]), p[order], 1.0)
    q = np.minimum.accumulate((ordered * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(q, 1)
    return out


def add_families(df, pcol):
    df["q_BH_within_family"] = np.nan
    df["family_size"] = 3
    for variant, idx in df.groupby("score_variant", sort=False).groups.items():
        df.loc[idx, "q_BH_within_family"] = bh_fixed(df.loc[idx, pcol])
        df.loc[idx, "family"] = (
            df.loc[idx, "analysis"] + "_" + variant + "_three_prespecified_relations"
        )
    return df


def adjusted_ranks(values, groups=None):
    ranks = rankdata(np.asarray(values, float), axis=0, method="average")
    n = len(ranks)
    if groups is None:
        within = ranks - ranks.mean(axis=0)
        controls_rank = 1
    else:
        groups = np.asarray(groups, int)
        m = int(groups.max()) + 1
        counts = np.bincount(groups, minlength=m)
        if np.any(counts == 0):
            raise ValueError("empty cluster label")
        sums = np.zeros((m, ranks.shape[1]))
        np.add.at(sums, groups, ranks)
        within = ranks - (sums / counts[:, None])[groups]
        controls_rank = m
    z = within[:, 4]
    zz = float(z @ z)
    if zz <= 1e-12:
        raise ValueError("no residual IFN rank variation")
    residuals = within - np.outer(z, z @ within / zz)
    return (ranks, within, z, residuals, controls_rank + 1)


def corr(x, y):
    den = np.sqrt(x @ x * (y @ y))
    if den <= 1e-12:
        return np.nan
    return float(np.clip(x @ y / den, -1, 1))


def approx_t(r, df):
    if not np.isfinite(r) or df <= 0:
        return (np.nan, np.nan)
    t = r * np.sqrt(df / max(1 - r * r, np.finfo(float).tiny))
    return (float(t), float(2 * stats.t.sf(abs(t), df)))


def six_r(values, groups=None):
    _, _, _, residuals, _ = adjusted_ranks(values, groups)
    return np.array([corr(residuals[:, x], residuals[:, y]) for _, _, x, y in MODELS])


def between_models(values):
    _, _, _, res, kr = adjusted_ranks(values)
    rows = []
    for variant, label, x, y in MODELS:
        r = corr(res[:, x], res[:, y])
        t, p = approx_t(r, len(values) - kr - 1)
        rows.append(
            dict(
                analysis="patient_means",
                score_variant=variant,
                relation=label,
                n_patients=len(values),
                rho=r,
                t_reference=t,
                df_reference=len(values) - kr - 1,
                p_t_approx=p,
                controls_rank=kr,
            )
        )
    return add_families(pd.DataFrame(rows), "p_t_approx")


def cluster_model(x, y, z, groups, controls_rank, signs):
    n = len(x)
    m = int(groups.max()) + 1
    k = controls_rank + 1
    xx = float(x @ x)
    zz = float(z @ z)
    if n <= k or m < 3 or xx <= 1e-12 or (float(y @ y) <= 1e-12):
        raise ValueError("non-estimable clustered model")
    beta = float(x @ y / xx)
    error = y - beta * x
    correction = m / (m - 1) * (n - 1) / (n - k)
    cluster_scores = np.bincount(groups, weights=x * error, minlength=m)
    se = float(np.sqrt(correction * (cluster_scores @ cluster_scores)) / xx)
    t_obs = beta / se
    a = np.bincount(groups, weights=x * y, minlength=m)
    d = np.bincount(groups, weights=x * x, minlength=m)
    c = np.bincount(groups, weights=x * z, minlength=m)
    zy = np.bincount(groups, weights=z * y, minlength=m)
    beta_star = signs @ a / xx
    zcoef_star = signs @ zy / zz
    cluster_star = signs * a - zcoef_star[:, None] * c - beta_star[:, None] * d
    se_star = np.sqrt(correction * np.sum(cluster_star**2, axis=1)) / xx
    t_star = beta_star / se_star
    if not np.isfinite(t_star).all():
        raise ValueError("nonfinite wild-cluster bootstrap statistic")
    b = int(np.sum(np.abs(t_star) >= abs(t_obs)))
    return (
        dict(
            beta_rank=beta,
            se_CR1=se,
            t_CR1=t_obs,
            df_CR1=m - 1,
            p_CR1_reference=float(2 * stats.t.sf(abs(t_obs), m - 1)),
            full_design_rank=k,
            CR1_correction=correction,
            wild_exceedances=b,
            wild_B=len(signs),
            p_wild_cluster=(b + 1) / (len(signs) + 1),
        ),
        t_star,
    )


def within_models(values, groups, signs):
    _, _, z, res, kr = adjusted_ranks(values, groups)
    rows = []
    draws = {}
    for variant, label, x, y in MODELS:
        info, t_star = cluster_model(res[:, x], res[:, y], z, groups, kr, signs)
        rows.append(
            dict(
                analysis="within_patient",
                score_variant=variant,
                relation=label,
                n_patients=int(groups.max()) + 1,
                n_visits=len(values),
                rho=corr(res[:, x], res[:, y]),
                **info,
            )
        )
        draws[f"{variant}_{label}"] = t_star
    return (add_families(pd.DataFrame(rows), "p_wild_cluster"), draws)


def bootstrap_between(values, indices):
    out = np.full((len(indices), len(MODELS)), np.nan)
    rejected = []
    for b, idx in enumerate(indices):
        try:
            out[b] = six_r(values[idx])
        except ValueError as e:
            rejected.append(dict(analysis="patient_means", replicate=b, reason=str(e)))
    return (out, rejected)


def bootstrap_within(values, groups, indices):
    blocks = [values[groups == g] for g in range(int(groups.max()) + 1)]
    out = np.full((len(indices), len(MODELS)), np.nan)
    rejected = []
    for b, sampled in enumerate(indices):
        selected = [blocks[g] for g in sampled]
        vals = np.concatenate(selected, axis=0)
        boot_groups = np.repeat(np.arange(len(selected)), [len(x) for x in selected])
        try:
            out[b] = six_r(vals, boot_groups)
        except ValueError as e:
            rejected.append(dict(analysis="within_patient", replicate=b, reason=str(e)))
    return (out, rejected)


def attach_ci(df, draws):
    df["ci_lo"] = np.nanquantile(draws, 0.025, axis=0)
    df["ci_hi"] = np.nanquantile(draws, 0.975, axis=0)
    df["bootstrap_valid"] = np.isfinite(draws).sum(axis=0)
    df["bootstrap_B"] = len(draws)
    df["ci_method"] = "patient_resampling_percentile_rerank_refit"
    return df


def gene_overlap_audit():
    members = pd.read_csv(INP / "program_membership.tsv", sep="\t")
    keep = members.matched.astype(str).str.lower().eq("true")
    members = members[keep]
    rows = []
    for variant in ["full", "noCCL5"]:
        retained = (
            members
            if variant == "full"
            else members[members.retained_in_noCCL5.astype(str).str.lower().eq("true")]
        )
        gene_sets = {p: set(g.gene) for p, g in retained.groupby("program")}
        for i, p1 in enumerate(sorted(gene_sets)):
            for p2 in sorted(gene_sets)[i + 1 :]:
                common = sorted(gene_sets[p1] & gene_sets[p2])
                rows.append(
                    dict(
                        score_variant=variant,
                        program1=p1,
                        program2=p2,
                        n1=len(gene_sets[p1]),
                        n2=len(gene_sets[p2]),
                        n_shared=len(common),
                        shared_genes=";".join(common),
                    )
                )
        expected = {"CCL5"} if variant == "full" else set()
        if gene_sets["CD8_GZMK_STATE"] & gene_sets["CYTOTOXICITY_SHARED"] != expected:
            raise ValueError("unexpected shared members in GZMK/CYTO")
    save(pd.DataFrame(rows), "program_gene_overlap_audit.tsv")


def crosscell(b_boot, seed):
    source = pd.read_csv(INP / "gse135779_child_sle_existing_scores.tsv", sep="\t")
    key = ["subject_id", "cell_type_analysis", "program"]
    if source.duplicated(key).any():
        raise ValueError("duplicate donor-celltype-program score")
    pivot = source.pivot(
        index="subject_id", columns=["cell_type_analysis", "program"], values="mean_score"
    )
    ifn = pivot.loc[:, pivot.columns.get_level_values(1) == "IFN_I"]
    z = ifn.mean(axis=1)
    aud = pd.DataFrame(
        {
            "subject_id": pivot.index,
            "n_celltypes_for_systemic_IFN": ifn.notna().sum(axis=1),
            "systemic_IFN": z,
        }
    ).reset_index(drop=True)
    save(aud, "gse135779_systemic_IFN_audit.tsv")
    cyto_pair = pivot.loc[
        :, [("CD8_T_GZMK", "CYTOTOXICITY"), ("CD8_T_Cytotoxic", "CYTOTOXICITY")]
    ].dropna()
    rr, pp = stats.spearmanr(cyto_pair.iloc[:, 0], cyto_pair.iloc[:, 1])
    save(
        pd.DataFrame(
            [
                dict(
                    x="CD8_T_GZMK_CYTOTOXICITY",
                    y="CD8_T_Cytotoxic_CYTOTOXICITY",
                    n_donors=len(cyto_pair),
                    rho_unadjusted=rr,
                    p_descriptive=pp,
                    interpretation="different_endpoint_not_evidence_of_functional_independence",
                )
            ]
        ),
        "gse135779_rho009_variable_audit.tsv",
    )
    programs = pd.read_csv(INP / "gse135779_requested_programs.tsv", sep="\t")["program"].tolist()
    candidate = [p for p in programs if p != "IFN_I"]
    if len(candidate) != 11 or len(pivot) != 33:
        raise ValueError("unexpected cross-cell family or donor count")
    xall = pivot.loc[:, [("CD14_Monocyte", p) for p in candidate]].to_numpy()
    yall = pivot.loc[:, [("CD8_T_GZMK", p) for p in candidate]].to_numpy()
    all_values = np.column_stack([xall, yall, z.to_numpy()])
    save(
        pd.DataFrame(
            all_values,
            columns=["CD14_" + p for p in candidate]
            + ["GZMK_" + p for p in candidate]
            + ["systemic_IFN"],
        ).assign(subject_id=pivot.index),
        "gse135779_model_inputs.tsv",
    )
    rows = []
    for i, p1 in enumerate(candidate):
        for j, p2 in enumerate(candidate):
            a = all_values[:, [i, 11 + j, 22]]
            a = a[np.isfinite(a).all(axis=1)]
            n = len(a)
            r = p = np.nan
            kr = 0
            reason = ""
            if n < 4:
                reason = "insufficient_complete_donors"
            else:
                ranked = rankdata(a, axis=0)
                c = np.column_stack([np.ones(n), ranked[:, 2]])
                kr = np.linalg.matrix_rank(c)
                er = ranked[:, :2] - c @ np.linalg.lstsq(c, ranked[:, :2], rcond=None)[0]
                r = corr(er[:, 0], er[:, 1])
                _, p = approx_t(r, n - kr - 1)
                if not np.isfinite(r):
                    reason = "zero_residual_variation"
            rows.append(
                dict(
                    CD14_program=p1,
                    GZMK_program=p2,
                    n_donors=n,
                    controls_rank=kr,
                    df_reference=n - kr - 1,
                    rho=r,
                    p_t_approx=p,
                    reason=reason,
                )
            )
    tab = pd.DataFrame(rows)
    tab["q_BH_121"] = bh_fixed(tab["p_t_approx"])
    tab["family_size"] = 121
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(all_values), (b_boot, len(all_values)))
    draws = np.full((b_boot, 121), np.nan)
    if np.isfinite(all_values).all():
        for b, idx in enumerate(indices):
            rk = rankdata(all_values[idx], axis=0)
            c = np.column_stack([np.ones(len(idx)), rk[:, -1]])
            er = rk[:, :-1] - c @ np.linalg.lstsq(c, rk[:, :-1], rcond=None)[0]
            norm = np.sqrt(np.sum(er**2, axis=0))
            norm[norm <= 1e-12] = np.nan
            er = er / norm
            draws[b] = (er[:, :11].T @ er[:, 11:]).ravel()
    else:
        for b, idx in enumerate(indices):
            a = all_values[idx]
            for h in range(121):
                sub = a[:, [h // 11, 11 + h % 11, 22]]
                sub = sub[np.isfinite(sub).all(axis=1)]
                if len(sub) < 4:
                    continue
                rk = rankdata(sub, axis=0)
                c = np.column_stack([np.ones(len(rk)), rk[:, -1]])
                er = rk[:, :2] - c @ np.linalg.lstsq(c, rk[:, :2], rcond=None)[0]
                draws[b, h] = corr(er[:, 0], er[:, 1])
    tab["ci_lo"] = np.nanquantile(draws, 0.025, axis=0)
    tab["ci_hi"] = np.nanquantile(draws, 0.975, axis=0)
    tab["bootstrap_valid"] = np.isfinite(draws).sum(axis=0)
    tab["bootstrap_B"] = b_boot
    save(tab, "gse135779_crosscell_all121.tsv")
    np.savez_compressed(
        OUT / "gse135779_bootstrap_draws.npz",
        indices=indices,
        rho=draws,
        CD14_program=tab.CD14_program.to_numpy(dtype=str),
        GZMK_program=tab.GZMK_program.to_numpy(dtype=str),
    )
    return tab


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--wild", type=int, default=19999)
    parser.add_argument("--seed", type=int, default=20261115)
    args = parser.parse_args()
    start = time.time()
    OUT.mkdir(exist_ok=True)
    sources = [
        INP / "sle_visit_scores.tsv",
        INP / "gse135779_child_sle_existing_scores.tsv",
        INP / "gse135779_requested_programs.tsv",
        Path(__file__),
    ]
    before = {str(p): sha(p) for p in sources}
    gene_overlap_audit()
    dat = pd.read_csv(INP / "sle_visit_scores.tsv", sep="\t", dtype={"patient_id": str})
    if dat.sample_id.duplicated().any():
        raise ValueError("duplicate sample id")
    dat[COLS] = dat[COLS].apply(pd.to_numeric, errors="coerce")
    eligible = dat["patient_id"].notna() & dat[COLS].notna().all(axis=1)
    save(dat.loc[~eligible, ["sample_id", "patient_id"] + COLS], "excluded_visits.tsv")
    dat = dat.loc[eligible].sort_values(["patient_id", "sample_id"]).reset_index(drop=True)
    assert dat.disease_state.eq("Systemic Lupus Erythematosus").all()
    vals = dat[COLS].to_numpy(float)
    means = dat.groupby("patient_id", sort=True)[COLS].mean()
    save(means.reset_index(), "gse121239_patient_mean_scores.tsv")
    counts = dat.groupby("patient_id").size().rename("n_visits").reset_index()
    counts["included_longitudinal"] = counts.n_visits.ge(2)
    save(counts, "patient_inclusion.tsv")
    means_array = means.to_numpy(float)
    groups_all = pd.factorize(dat.patient_id, sort=True)[0]
    long = dat[dat.patient_id.isin(counts.loc[counts.included_longitudinal, "patient_id"])].copy()
    groups = pd.factorize(long.patient_id, sort=True)[0]
    long_vals = long[COLS].to_numpy(float)
    save(long[["sample_id", "patient_id"] + COLS], "gse121239_longitudinal_inputs.tsv")
    print(
        f"Inputs: {len(means)} patients, {len(dat)} visits; longitudinal {groups.max() + 1} patients, {len(long)} visits",
        flush=True,
    )
    rng = np.random.default_rng(args.seed)
    mean_indices = rng.integers(0, len(means), (args.bootstrap, len(means)))
    m_long = int(groups.max()) + 1
    long_indices = rng.integers(0, m_long, (args.bootstrap, m_long))
    signs = rng.choice(np.array([-1.0, 1.0]), size=(args.wild, m_long))
    mean_tab = between_models(means_array)
    within_tab, wild_draws = within_models(long_vals, groups, signs)
    print("Computing patient-resampling confidence intervals", flush=True)
    boot_mean, rejected_mean = bootstrap_between(means_array, mean_indices)
    boot_long, rejected_long = bootstrap_within(long_vals, groups, long_indices)
    mean_tab = attach_ci(mean_tab, boot_mean)
    within_tab = attach_ci(within_tab, boot_long)
    save(mean_tab, "gse121239_patient_means_corrected.tsv")
    save(within_tab, "gse121239_longitudinal_corrected.tsv")
    save(
        pd.DataFrame(rejected_mean + rejected_long, columns=["analysis", "replicate", "reason"]),
        "bootstrap_rejections.tsv",
    )
    np.savez_compressed(
        OUT / "gse121239_resampling_audit.npz",
        mean_indices=mean_indices,
        longitudinal_indices=long_indices,
        wild_signs=signs.astype(np.int8),
        mean_rho=boot_mean,
        longitudinal_rho=boot_long,
        **wild_draws,
    )
    delta = []
    for name, table, boot in [
        ("patient_means", mean_tab, boot_mean),
        ("within_patient", within_tab, boot_long),
    ]:
        for j, label in enumerate(["GZMK", "B_cell", "Neutrophil"]):
            d = boot[:, j] - boot[:, 3 + j]
            delta.append(
                dict(
                    analysis=name,
                    relation=label,
                    delta_noCCL5_minus_full=float(table.iloc[j].rho - table.iloc[3 + j].rho),
                    ci_lo=float(np.nanquantile(d, 0.025)),
                    ci_hi=float(np.nanquantile(d, 0.975)),
                    bootstrap_valid=int(np.isfinite(d).sum()),
                )
            )
    save(pd.DataFrame(delta), "CCL5_paired_sensitivity_differences.tsv")
    descriptive = []
    for variant, y in [("full", 1), ("noCCL5", 6)]:
        r, p = stats.spearmanr(means_array[:, 4], means_array[:, y])
        descriptive.append(
            dict(
                analysis="patient_means",
                score_variant=variant,
                relation="IFN_vs_CYTO_unadjusted",
                rho=r,
                p_approx_descriptive=p,
                n_patients=len(means),
                interpretation="not_test_of_independence",
            )
        )
    save(pd.DataFrame(descriptive), "IFN_unadjusted_descriptive.tsv")
    print("Rechecking all 121 prespecified cross-cell candidates", flush=True)
    cross = crosscell(args.bootstrap, args.seed + 1)
    after = {str(p): sha(p) for p in sources}
    if before != after:
        raise RuntimeError("source changed during run")
    save(
        pd.DataFrame(
            [
                {"path": k, "sha256_before": v, "sha256_after": after[k], "unchanged": True}
                for k, v in before.items()
            ]
        ),
        "input_and_script_sha256.tsv",
    )
    metadata = dict(
        arguments=vars(args),
        python=sys.version,
        numpy=np.__version__,
        pandas=pd.__version__,
        scipy=scipy.__version__,
        platform=platform.platform(),
        elapsed_seconds=time.time() - start,
        n_bulk_patients=len(means),
        n_bulk_visits=len(dat),
        n_longitudinal_patients=m_long,
        n_longitudinal_visits=len(long),
        mean_bootstrap_rejections=len(rejected_mean),
        longitudinal_bootstrap_rejections=len(rejected_long),
        crosscell_family_size=len(cross),
        crosscell_estimable=int(cross.rho.notna().sum()),
    )
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(mean_tab.to_string(index=False), flush=True)
    print(within_tab.to_string(index=False), flush=True)
    print(cross.sort_values("p_t_approx").head(8).to_string(index=False), flush=True)
    print(f"DONE in {time.time() - start:.1f}s", flush=True)


if __name__ == "__main__":
    main()
