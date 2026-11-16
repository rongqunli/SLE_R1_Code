import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
from pathlib import Path
import json
import hashlib
import platform
import sys
import importlib.metadata as im
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent
INP = ROOT / "formal_inputs"
OLD = Path("F:\\program\\Hei\\SLE_R1\\Try")
FOCUS = ["IRF7", "STAT1", "STAT2", "IRF9"]
SEED = 20261115


def sha256(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for data in iter(lambda: f.read(1048576), b""):
            h.update(data)
    return h.hexdigest()


def zcols(x):
    sd = np.std(x, axis=0, ddof=1)
    valid = sd > 1e-12
    z = np.zeros_like(x, dtype=float)
    z[:, valid] = (x[:, valid] - np.mean(x[:, valid], axis=0)) / sd[valid]
    return z


def design(meta, kind, expected_columns=None):
    blocks = [np.ones((len(meta), 1))]
    labels = ["Intercept"]
    cats = []
    cont = []
    if kind == "all_adjusted":
        cats = ["condition", "batch", "sex"]
        cont = ["age", "log_mean_umi_per_cell"]
    elif kind == "sle_adjusted":
        cats = ["batch", "sex"]
        cont = ["age", "log_mean_umi_per_cell", "SLEDAI"]
    elif kind == "all_technical":
        cats = ["condition", "batch"]
        cont = ["log_mean_umi_per_cell"]
    elif kind != "raw":
        raise ValueError(kind)
    for col in cats:
        dum = pd.get_dummies(meta[col].astype(str), prefix=col, drop_first=True, dtype=float)
        if dum.shape[1]:
            blocks.append(dum.to_numpy())
            labels.extend(dum.columns.tolist())
    for col in cont:
        v = meta[col].to_numpy(float)
        if not np.isfinite(v).all():
            raise ValueError(f"Missing {col} in design")
        if np.ptp(v) > 0:
            ranked = stats.rankdata(v, method="average")
            blocks.append(((ranked - ranked.mean()) / ranked.std())[:, None])
            labels.append("rank_" + col)
    c = np.column_stack(blocks)
    if expected_columns is not None:
        tmp = pd.DataFrame(c, columns=labels).reindex(columns=expected_columns, fill_value=0.0)
        c, labels = (tmp.to_numpy(), list(expected_columns))
    return (c, labels)


def residual(x, c):
    return x - c @ np.linalg.lstsq(c, x, rcond=None)[0]


def bh_fixed_family(p):
    good = np.isfinite(p)
    q = multipletests(np.where(good, p, 1.0), method="fdr_bh")[1]
    q[~good] = np.nan
    return q


def fit(x, y, c):
    rx = stats.rankdata(x, axis=0, method="average")
    ry = stats.rankdata(y, axis=0, method="average")
    ex = residual(rx, c)
    ey = residual(ry, c)
    ex -= ex.mean(axis=0)
    ey -= ey.mean(axis=0)
    nx = np.linalg.norm(ex, axis=0)
    ny = np.linalg.norm(ey, axis=0)
    denom = nx * ny
    r = np.divide(
        np.sum(ex * ey, axis=0), denom, out=np.full(x.shape[1], np.nan), where=denom > 1e-10
    )
    r = np.clip(r, -1, 1)
    df = len(x) - np.linalg.matrix_rank(c) - 1
    t = r * np.sqrt(df / np.maximum(1 - r * r, 1e-15))
    p = 2 * stats.t.sf(np.abs(t), df)
    q = bh_fixed_family(p)
    return (r, p, q, df, ex, ey)


def scores(expr, ifn_idx, tf_names, gene_index, eligible):
    z = zcols(expr[:, ifn_idx])
    idx_to_position = {idx: j for j, idx in enumerate(ifn_idx)}
    full = z.mean(axis=1)
    own = np.tile(full[:, None], (1, len(tf_names)))
    k = len(ifn_idx)
    for j, tf in enumerate(tf_names):
        pos = idx_to_position.get(gene_index[tf])
        if pos is not None:
            own[:, j] = (z.sum(axis=1) - z[:, pos]) / (k - 1)
    common_positions = [j for j, idx in enumerate(ifn_idx) if idx not in eligible]
    focus_positions = [
        j
        for j, idx in enumerate(ifn_idx)
        if idx not in {gene_index[g] for g in FOCUS if g in gene_index}
    ]
    return {
        "exclude_all_tested_TFs": z[:, common_positions].mean(axis=1)[:, None],
        "full_frozen": full[:, None],
        "exclude_self": own,
        "exclude_four_candidates": z[:, focus_positions].mean(axis=1)[:, None],
    }


def permute_all(x, y, c, meta, kind, nperm=99999):
    r, p, q, df, ex, ey = fit(x, y, c)
    yr = stats.rankdata(y[:, 0])
    fitted = yr - ey[:, 0]
    projector = np.eye(len(c)) - c @ np.linalg.pinv(c)
    valid = np.isfinite(r)
    xn = ex[:, valid] / np.linalg.norm(ex[:, valid], axis=0)
    if not np.isfinite(ey).all() or np.linalg.norm(ey) < 1e-10:
        raise RuntimeError("Permutation outcome is not estimable")
    strata = (
        meta["batch"].astype(str) + "|" + meta["condition"].astype(str)
        if kind == "all_adjusted"
        else meta["batch"].astype(str)
    )
    groups = [np.flatnonzero(strata.to_numpy() == v) for v in sorted(strata.unique())]
    rng = np.random.default_rng(SEED + (0 if kind == "all_adjusted" else 1))
    hits = np.zeros(int(valid.sum()), dtype=np.int64)
    for start in range(0, nperm, 2000):
        n = min(2000, nperm - start)
        ix = np.tile(np.arange(len(meta)), (n, 1))
        for g in groups:
            ix[:, g] = g[np.argsort(rng.random((n, len(g))), axis=1)]
        yp = fitted[None, :] + ey[:, 0][ix]
        ep = yp @ projector.T
        ep -= ep.mean(axis=1)[:, None]
        ep_norm = np.linalg.norm(ep, axis=1)
        if not np.isfinite(ep_norm).all() or np.any(ep_norm < 1e-10):
            raise RuntimeError("Degenerate permutation residual encountered")
        ep /= ep_norm[:, None]
        rp = ep @ xn
        hits += np.sum(np.abs(rp) >= np.abs(r[valid])[None, :] - 1e-12, axis=0)
    pp = np.full(len(r), np.nan)
    pp[valid] = (hits + 1) / (nperm + 1)
    qq = bh_fixed_family(pp)
    return (pp, qq, {str(v): int((strata == v).sum()) for v in sorted(strata.unique())})


def self_checks():
    rng = np.random.default_rng(7)
    x = rng.normal(size=(50, 3))
    y = rng.normal(size=(50, 1))
    c = np.ones((50, 1))
    r, p, q, *_ = fit(x, y, c)
    assert np.allclose(r, [stats.spearmanr(x[:, j], y[:, 0]).statistic for j in range(3)])
    z = rng.normal(size=50)
    c = np.column_stack([np.ones(50), stats.rankdata(z)])
    r, *_ = fit(x, y, c)
    for j in range(3):
        rxy = stats.spearmanr(x[:, j], y[:, 0]).statistic
        rxz = stats.spearmanr(x[:, j], z).statistic
        ryz = stats.spearmanr(y[:, 0], z).statistic
        expected = (rxy - rxz * ryz) / np.sqrt((1 - rxz**2) * (1 - ryz**2))
        assert np.isclose(r[j], expected)
    xx = np.column_stack([np.ones(50), x[:, 0]])
    mm = pd.DataFrame({"batch": ["B1"] * 50, "condition": ["SLE"] * 50})
    pp, qq, _ = permute_all(xx, y, np.ones((50, 1)), mm, "sle_adjusted", nperm=19)
    assert np.isnan(pp[0]) and np.isnan(qq[0]) and np.isfinite(pp[1])


def main():
    self_checks()
    m = pd.read_csv(INP / "formal_nk_metadata.tsv", sep="\t").set_index("donor", drop=False)
    e = pd.read_csv(INP / "formal_nk_tmm_logcpm.tsv", sep="\t", index_col=0)
    d = pd.read_csv(INP / "formal_nk_gene_detection.tsv", sep="\t").set_index("gene")
    assert e.index.is_unique and m.index.is_unique and (set(m.index) == set(e.columns))
    e = e.loc[:, m.index]
    expr = e.to_numpy().T
    genes = e.index.tolist()
    gi = {g: i for i, g in enumerate(genes)}
    frozen = list(
        dict.fromkeys(pd.read_csv(INP / "formal_frozen_ifn_members.tsv", sep="\t")["gene"].tolist())
    )
    requested_tf = list(dict.fromkeys((OLD / "hs_hgnc_tfs.txt").read_text().split()))
    screen = []
    for tf in requested_tf:
        present = tf in gi
        detected = present and int(d.loc[tf, "donors_raw_cpm_ge1"]) >= 5
        variable = present and np.std(expr[:, gi[tf]]) > 1e-12
        screen.append(
            dict(
                TF=tf,
                present=present,
                detected_cpm1_in5=detected,
                variable=variable,
                eligible=detected and variable,
            )
        )
    screen = pd.DataFrame(screen)
    screen.to_csv(ROOT / "tf_eligibility.tsv", sep="\t", index=False)
    tfs = screen.loc[screen.eligible, "TF"].tolist()
    ti = [gi[g] for g in tfs]
    ifn = [g for g in frozen if g in gi and np.std(expr[:, gi[g]]) > 1e-12]
    ii = [gi[g] for g in ifn]
    eligible = set(ti)
    focus_indices = [tfs.index(g) for g in FOCUS]
    annotation = pd.DataFrame({"gene": frozen})
    annotation["matched"] = annotation.gene.isin(gi)
    annotation["variable"] = annotation.gene.isin(ifn)
    annotation["tested_TF"] = annotation.gene.isin(tfs)
    annotation["primary_score_member"] = annotation.variable & ~annotation.tested_TF
    annotation.to_csv(ROOT / "ifn_gene_audit.tsv", sep="\t", index=False)
    print(
        "Inputs:",
        len(m),
        "donors;",
        int(m.n_cells.sum()),
        "cells;",
        len(tfs),
        "eligible TFs;",
        len(ifn),
        "variable IFN;",
        sum(annotation.primary_score_member),
        "primary score genes",
        flush=True,
    )
    for col in ["age", "SLEDAI", "log_mean_umi_per_cell"]:
        m[col] = pd.to_numeric(m[col], errors="coerce")
    req = ["condition", "batch", "sex", "age", "log_mean_umi_per_cell"]
    mask_all = m[req].notna().all(axis=1).to_numpy()
    mask_sle = (m.condition.eq("SLE") & m[req + ["SLEDAI"]].notna().all(axis=1)).to_numpy()
    model_defs = [
        ("all_adjusted", mask_all, "all_adjusted"),
        ("sle_adjusted", mask_sle, "sle_adjusted"),
        ("all_raw", mask_all, "raw"),
        ("sle_raw", mask_sle, "raw"),
        ("all_technical", mask_all, "all_technical"),
    ]
    all_rows = []
    designs = []
    score_rows = []
    cache = {}
    score_defs = scores(expr, ii, tfs, gi, eligible)
    observed = pd.DataFrame(
        {name: val[:, 0] for name, val in score_defs.items() if name != "exclude_self"},
        index=m.index,
    )
    observed.to_csv(ROOT / "donor_ifn_scores.tsv", sep="\t")
    for name, mask, kind in model_defs:
        meta = m.iloc[np.flatnonzero(mask)].copy()
        c, cl = design(meta, kind)
        x = expr[mask][:, ti]
        rank = int(np.linalg.matrix_rank(c))
        if rank != c.shape[1]:
            raise RuntimeError(f"{name} design rank deficient {rank}/{c.shape[1]}")
        designs.append(
            dict(
                model=name,
                n=len(meta),
                SLE=int(meta.condition.eq("SLE").sum()),
                HD=int(meta.condition.eq("HD").sum()),
                nuisance_rank=rank,
                partial_df=len(meta) - rank - 1,
                covariates=";".join(cl),
                donors=";".join(meta.donor),
            )
        )
        pd.DataFrame(c, columns=cl, index=meta.index).to_csv(ROOT / f"design_{name}.tsv", sep="\t")
        for score_name, ys in score_defs.items():
            y = ys[mask]
            r, p, q, df, ex, ey = fit(x, y, c)
            table = pd.DataFrame(
                dict(
                    TF=tfs,
                    model=name,
                    score=score_name,
                    n=len(meta),
                    df=df,
                    rho=r,
                    p_value=p,
                    fdr=q,
                )
            )
            table["estimable"] = np.isfinite(r)
            table["BH_family_size"] = len(tfs)
            table["positive_rank"] = table.rho.rank(ascending=False, method="min")
            table["absolute_rank"] = table.rho.abs().rank(ascending=False, method="min")
            if name in ["all_adjusted", "sle_adjusted"] and score_name == "exclude_all_tested_TFs":
                print("Permutation:", name, flush=True)
                pp, qq, strata = permute_all(x, y, c, meta, kind)
                table["permutation_p"] = pp
                table["permutation_fdr"] = qq
                cache[name] = (mask, kind, table.copy())
                (ROOT / f"permutation_{name}.json").write_text(
                    json.dumps(
                        {
                            "B": 99999,
                            "seed": SEED + (name == "sle_adjusted"),
                            "strata": strata,
                            "method": "Freedman-Lane on rank-scale outcome residuals; approximate nuisance-adjusted restricted permutation",
                        },
                        indent=2,
                    ),
                    encoding="utf-8",
                )
            all_rows.append(table)
        print("Done", name, "n", len(meta), "df", len(meta) - rank - 1, flush=True)
    all_results = pd.concat(all_rows, ignore_index=True)
    all_results.to_csv(ROOT / "all_tf_results.tsv", sep="\t", index=False)
    all_results[all_results.TF.isin(FOCUS)].to_csv(
        ROOT / "four_tf_model_comparison.tsv", sep="\t", index=False
    )
    pd.DataFrame(designs).to_csv(ROOT / "model_designs.tsv", sep="\t", index=False)
    boot_all = []
    lodo_all = []
    ci_all = []
    difference = []
    boot_audit = []
    for name, (mask, kind, table) in cache.items():
        ix = np.flatnonzero(mask)
        mm = m.iloc[ix]
        ee = expr[ix]
        fixed_y = score_defs["exclude_all_tested_TFs"][mask]
        _, base_columns = design(mm, kind)
        rng = np.random.default_rng(SEED + 100 + (name == "sle_adjusted"))
        groups = [
            np.flatnonzero(mm.condition.to_numpy() == v) for v in sorted(mm.condition.unique())
        ]
        rejected = 0
        for b in range(2000):
            bi = np.concatenate([rng.choice(g, size=len(g), replace=True) for g in groups])
            bm = mm.iloc[bi].reset_index(drop=True)
            be = ee[bi]
            cy = fixed_y[bi]
            cc, _ = design(bm, kind, base_columns)
            audit = dict(
                model=name,
                replicate=b,
                unique_donors=int(bm.donor.nunique()),
                observed_sex_levels=";".join(sorted(bm.sex.unique())),
                observed_batches=";".join(sorted(bm.batch.unique())),
                nuisance_rank=int(np.linalg.matrix_rank(cc)),
                expected_rank=cc.shape[1],
            )
            if audit["nuisance_rank"] != cc.shape[1]:
                rejected += 1
                boot_audit.append(dict(audit, status="rejected_rank_deficient"))
                continue
            rr, _, _, df, _, _ = fit(be[:, [gi[g] for g in FOCUS]], cy, cc)
            if not np.isfinite(rr).all():
                rejected += 1
                boot_audit.append(dict(audit, status="rejected_unestimable_correlation"))
                continue
            boot_audit.append(dict(audit, status="valid"))
            boot_all.extend(
                (
                    dict(model=name, replicate=b, TF=tf, rho=float(rr[j]))
                    for j, tf in enumerate(FOCUS)
                )
            )
            difference.append(
                dict(model=name, replicate=b, delta_IRF7_minus_STAT1=float(rr[0] - rr[1]))
            )
        for j in range(len(ix)):
            keep = np.arange(len(ix)) != j
            lm = mm.iloc[np.flatnonzero(keep)]
            le = ee[keep]
            ly = fixed_y[keep]
            lc, _ = design(lm, kind, base_columns)
            if np.linalg.matrix_rank(lc) != lc.shape[1]:
                raise RuntimeError(
                    f"Leave-one-out design not estimable: {name}, {mm.iloc[j].donor}"
                )
            rr, pp, qq, df, _, _ = fit(le[:, ti], ly, lc)
            ranks = pd.Series(rr).rank(ascending=False, method="min").to_numpy()
            for tf, pos in zip(FOCUS, focus_indices):
                lodo_all.append(
                    dict(
                        model=name,
                        omitted_donor=mm.iloc[j].donor,
                        TF=tf,
                        rho=rr[pos],
                        p_value=pp[pos],
                        fdr=qq[pos],
                        positive_rank=ranks[pos],
                        df=df,
                    )
                )
        bb = pd.DataFrame([r for r in boot_all if r["model"] == name])
        ll = pd.DataFrame([r for r in lodo_all if r["model"] == name])
        dd = pd.DataFrame([r for r in difference if r["model"] == name])
        for tf in FOCUS:
            base = table.set_index("TF").loc[tf]
            bvals = bb.loc[bb.TF == tf, "rho"]
            lvals = ll.loc[ll.TF == tf]
            ci_all.append(
                dict(
                    model=name,
                    TF=tf,
                    rho=base.rho,
                    fdr=base.fdr,
                    permutation_fdr=base.permutation_fdr,
                    positive_rank=base.positive_rank,
                    bootstrap_low=bvals.quantile(0.025),
                    bootstrap_high=bvals.quantile(0.975),
                    bootstrap_valid=len(bvals),
                    bootstrap_rejected=rejected,
                    lodo_min=lvals.rho.min(),
                    lodo_max=lvals.rho.max(),
                    lodo_positive_n=int((lvals.rho > 0).sum()),
                    lodo_n=len(lvals),
                    lodo_rank_min=lvals.positive_rank.min(),
                    lodo_rank_max=lvals.positive_rank.max(),
                )
            )
        print(
            name,
            "bootstrap valid",
            2000 - rejected,
            "IRF7-STAT1 diff 95%",
            dd.delta_IRF7_minus_STAT1.quantile([0.025, 0.975]).to_list(),
            flush=True,
        )
    pd.DataFrame(boot_all).to_csv(ROOT / "bootstrap_four_tfs.tsv", sep="\t", index=False)
    pd.DataFrame(boot_audit).to_csv(ROOT / "bootstrap_audit.tsv", sep="\t", index=False)
    pd.DataFrame(difference).to_csv(ROOT / "bootstrap_irf7_minus_stat1.tsv", sep="\t", index=False)
    pd.DataFrame(lodo_all).to_csv(ROOT / "leave_one_donor_out.tsv", sep="\t", index=False)
    summary = pd.DataFrame(ci_all)
    summary.to_csv(ROOT / "four_tf_primary_summary.tsv", sep="\t", index=False)
    differences = []
    for name, block in pd.DataFrame(difference).groupby("model", sort=False):
        current = summary[summary.model.eq(name)].set_index("TF")
        differences.append(
            dict(
                model=name,
                observed_difference=current.loc["IRF7", "rho"] - current.loc["STAT1", "rho"],
                bootstrap_95_low=block.delta_IRF7_minus_STAT1.quantile(0.025),
                bootstrap_95_high=block.delta_IRF7_minus_STAT1.quantile(0.975),
                bootstrap_valid=len(block),
                bootstrap_attempted=2000,
            )
        )
    pd.DataFrame(differences).to_csv(
        ROOT / "irf7_stat1_difference_summary.tsv", sep="\t", index=False
    )
    inputs = [
        INP / "formal_nk_tmm_logcpm.tsv",
        INP / "formal_nk_metadata.tsv",
        INP / "formal_nk_gene_detection.tsv",
        INP / "formal_frozen_ifn_members.tsv",
        OLD / "hs_hgnc_tfs.txt",
        Path(__file__),
    ]
    versions = {p: im.version(p) for p in ["numpy", "pandas", "scipy", "statsmodels"]}
    (ROOT / "run_manifest.json").write_text(
        json.dumps(
            {
                "python": sys.version,
                "platform": platform.platform(),
                "seed": SEED,
                "versions": versions,
                "input_sha256": {str(p): sha256(p) for p in inputs},
                "primary_tf_count": len(tfs),
                "primary_ifn_count": int(sum(annotation.primary_score_member)),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
