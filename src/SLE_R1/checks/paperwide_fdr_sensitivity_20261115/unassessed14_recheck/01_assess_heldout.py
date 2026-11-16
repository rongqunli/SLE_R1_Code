import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
from pathlib import Path
import hashlib
import json
import platform
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "heldout_inference"
INPUT = ROOT / "heldout_inputs" / "patient_predictions.tsv"
ORIGINAL = Path(
    "F:/program/Hei/SLE_v1/results/reviewer_response/02_clinical_landscape_validation/heldout_continuous_vs_four_group.tsv"
)
B = 99999
B_CI = 5000
AUDIT_N = 1000
BATCH = 512
TRAITS = ["SLEDAI", "C3", "C4", "anti_dsDNA"]
SETS = ["ILLUMINATE1_validation", "ILLUMINATE2_validation", "pooled_trial_linked_validation"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_exact(path):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def save(frame, name):
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.17g")


def cp_interval(k, n):
    return (
        0.0 if k == 0 else stats.beta.ppf(0.025, k, n - k + 1),
        1.0 if k == n else stats.beta.ppf(0.975, k + 1, n - k),
    )


def main():
    OUT.mkdir(exist_ok=True)
    (OUT / "audit_indices").mkdir(exist_ok=True)
    plan_hash = sha(Path(__file__))
    input_hash = sha(INPUT)
    assert "PASS:" in (ROOT / "heldout_inputs" / "RESTORATION_STATUS.txt").read_text()
    dat = read_exact(INPUT)
    for c in ["observed", "pred_continuous", "pred_four_group"]:
        dat[c] = dat[c].map(float)
    assert len(dat) == 1761 and dat.subject_id.nunique() == 441
    assert not dat.duplicated(["trait", "subject_id"]).any()
    assert set(dat.trial) == set(SETS[:2])
    old = read_exact(ORIGINAL)
    old["source_file_line"] = np.arange(2, len(old) + 2)
    records, diagnostics, bootstrap_frames, ci_frames, audit_manifest = ([], [], [], [], [])
    for ti, trait in enumerate(TRAITS, 1):
        for si, validation_set in enumerate(SETS, 1):
            key = f"{ti:02d}_{si:02d}_{trait}"
            sub = dat[dat.trait.eq(trait)].copy()
            if validation_set != SETS[2]:
                sub = sub[sub.trial.eq(validation_set)].copy()
            sub = sub.reset_index(drop=True)
            y = sub.observed.to_numpy(float)
            cont, four = (sub.pred_continuous.to_numpy(float), sub.pred_four_group.to_numpy(float))
            loss = (y - four) ** 2 - (y - cont) ** 2
            assert np.isfinite(loss).all() and np.var(y) > 0
            trial_names = [s for s in SETS[:2] if s in set(sub.trial)]
            stratum_idx = [np.flatnonzero(sub.trial.eq(t).to_numpy()) for t in trial_names]
            loss_h = [loss[ii] for ii in stratum_idx]
            weights = np.array([len(x) / len(loss) for x in loss_h])
            ns = np.array([len(x) for x in loss_h])
            assert (ns >= 2).all() and ns.sum() == len(y)
            mu = sum((w * x.mean() for w, x in zip(weights, loss_h)))
            variances = np.array([w * w * x.var(ddof=1) / len(x) for w, x in zip(weights, loss_h)])
            se = np.sqrt(variances.sum())
            assert se > 0 and np.isfinite(se)
            t_obs = mu / se
            df = variances.sum() ** 2 / np.sum(variances**2 / (ns - 1))
            p_t = 2 * stats.t.sf(abs(t_obs), df)
            r2_four = 1 - np.sum((y - four) ** 2) / np.sum((y - y.mean()) ** 2)
            r2_cont = 1 - np.sum((y - cont) ** 2) / np.sum((y - y.mean()) ** 2)
            delta = r2_cont - r2_four
            assert abs(delta - mu / np.var(y)) < 1e-12
            legacy = old[old.trait.eq(trait) & old.validation_set.eq(validation_set)].iloc[0]
            assert abs(delta - float(legacy.observed_delta_r2_continuous_minus_four)) < 1e-12
            line = int(legacy.source_file_line)
            hypothesis_id = (
                f"heldout_performance_12_unassessed::clinical_heldout_bootstrap_12:{line}"
            )
            sub["loss_difference"] = loss
            sub["analysis_key"] = key
            save(sub, f"audit_indices/{key}_patient_order.tsv")
            seed = 20261115 + 1000 * ti + si
            rng = np.random.default_rng(seed)
            null_t = np.empty(B)
            audit_parts = [[] for _ in loss_h]
            for start in range(0, B, BATCH):
                size = min(BATCH, B - start)
                centered_means = np.zeros(size)
                boot_var = np.zeros(size)
                for h, (w, d) in enumerate(zip(weights, loss_h)):
                    idx = rng.integers(0, len(d), size=(size, len(d)))
                    draw = d[idx] - mu
                    centered_means += w * draw.mean(axis=1)
                    boot_var += w * w * draw.var(axis=1, ddof=1) / len(d)
                    if start < AUDIT_N:
                        audit_parts[h].append(idx[: min(size, AUDIT_N - start)])
                with np.errstate(divide="ignore", invalid="ignore"):
                    null_t[start : start + size] = centered_means / np.sqrt(boot_var)
            invalid = int((~np.isfinite(null_t)).sum())
            assert invalid == 0, "Invalid draws must be investigated, not silently discarded"
            exceedances = int(np.count_nonzero(np.abs(null_t) >= abs(t_obs)))
            p = (exceedances + 1) / (B + 1)
            mc_lo, mc_hi = cp_interval(exceedances, B)
            for h, name in enumerate(trial_names):
                idx = np.concatenate(audit_parts[h])
                assert idx.shape == (AUDIT_N, len(loss_h[h]))
                filename = f"audit_indices/{key}_stratum_{h + 1}_indices.tsv.gz"
                pd.DataFrame(idx).to_csv(
                    OUT / filename, sep="\t", index=False, header=False, compression="gzip"
                )
                audit_manifest.append(
                    dict(
                        analysis_key=key,
                        trait=trait,
                        validation_set=validation_set,
                        stratum=h + 1,
                        trial=name,
                        n=len(loss_h[h]),
                        weight=weights[h],
                        draws=AUDIT_N,
                        index_base=0,
                        indices_file=filename,
                        patient_order_file=f"audit_indices/{key}_patient_order.tsv",
                    )
                )
            bootstrap_frames.append(
                pd.DataFrame(dict(analysis_key=key, replicate=np.arange(1, B + 1), t_star=null_t))
            )
            rng_ci = np.random.default_rng(seed + 1000000)
            ci_delta = np.empty(B_CI)
            for start in range(0, B_CI, BATCH):
                size = min(BATCH, B_CI - start)
                sum_y = np.zeros(size)
                sum_y2 = np.zeros(size)
                sum_loss = np.zeros(size)
                for ii in stratum_idx:
                    idx = ii[rng_ci.integers(0, len(ii), size=(size, len(ii)))]
                    sum_y += y[idx].sum(axis=1)
                    sum_y2 += (y[idx] ** 2).sum(axis=1)
                    sum_loss += loss[idx].sum(axis=1)
                denom = sum_y2 - sum_y**2 / len(y)
                assert (denom > 0).all()
                ci_delta[start : start + size] = sum_loss / denom
            ci_lo, ci_hi = np.quantile(ci_delta, [0.025, 0.975])
            ci_frames.append(
                pd.DataFrame(
                    dict(analysis_key=key, replicate=np.arange(1, B_CI + 1), delta_r2=ci_delta)
                )
            )
            records.append(
                dict(
                    hypothesis_id=hypothesis_id,
                    analysis_key=key,
                    trait=trait,
                    validation_set=validation_set,
                    n=len(y),
                    n_trial1=int(sum(sub.trial.eq(SETS[0]))),
                    n_trial2=int(sum(sub.trial.eq(SETS[1]))),
                    r2_continuous=r2_cont,
                    r2_four_group=r2_four,
                    delta_r2=delta,
                    mean_loss_difference=mu,
                    se_mean_loss_difference=se,
                    t_observed=t_obs,
                    p_primary=p,
                    p_primary_method="null_centered_studentized_patient_bootstrap_99999",
                    bootstrap_B=B,
                    seed=seed,
                    bootstrap_exceedances=exceedances,
                    invalid_draws=invalid,
                    mc_tail_probability_ci025=mc_lo,
                    mc_tail_probability_ci975=mc_hi,
                    mc_p_se=np.sqrt(p * (1 - p) / (B + 1)),
                    p_t_reference=p_t,
                    t_reference_df=df,
                    delta_r2_ci025=ci_lo,
                    delta_r2_ci975=ci_hi,
                    ci_draws=B_CI,
                    ci_seed=seed + 1000000,
                    ci_method="percentile_patient_triplets_trial_stratified_if_pooled",
                    legacy_raw_p=float(legacy.bootstrap_two_sided_p),
                    legacy_ci025=float(legacy.bootstrap_delta_q025),
                    legacy_ci975=float(legacy.bootstrap_delta_q975),
                    source_file_line=line,
                )
            )
            for name, ii in zip(trial_names + ["ALL"], stratum_idx + [np.arange(len(loss))]):
                d = loss[ii]
                centered = d - d.mean()
                diagnostics.append(
                    dict(
                        analysis_key=key,
                        trait=trait,
                        validation_set=validation_set,
                        stratum=name,
                        n=len(d),
                        mean_loss_difference=d.mean(),
                        sd_loss_difference=d.std(ddof=1),
                        skewness=stats.skew(d, bias=False),
                        max_centered_squared_share=np.max(centered**2) / np.sum(centered**2),
                        max_centered_absolute_share=np.max(np.abs(centered))
                        / np.sum(np.abs(centered)),
                        min_loss=d.min(),
                        max_loss=d.max(),
                        null_t_q001=np.quantile(null_t, 0.001),
                        null_t_q025=np.quantile(null_t, 0.025),
                        null_t_q975=np.quantile(null_t, 0.975),
                        null_t_q999=np.quantile(null_t, 0.999),
                    )
                )
            print(
                f"{trait} | {validation_set} | n={len(y)} delta={delta:.9g} p={p:.8g} MC95%=[{mc_lo:.6g},{mc_hi:.6g}]",
                flush=True,
            )
    result = pd.DataFrame(records)
    assert len(result) == 12 and (not result.hypothesis_id.duplicated().any())
    save(result, "heldout_final_results.tsv")
    save(pd.DataFrame(diagnostics), "loss_diagnostics.tsv")
    save(pd.DataFrame(audit_manifest), "audit_index_manifest.tsv")
    save(pd.concat(bootstrap_frames, ignore_index=True), "null_t_distribution.tsv.gz")
    save(pd.concat(ci_frames, ignore_index=True), "delta_ci_distribution.tsv.gz")
    assert sha(INPUT) == input_hash and sha(Path(__file__)) == plan_hash
    metadata = dict(
        input_path=str(INPUT),
        input_sha256=input_hash,
        plan_sha256=plan_hash,
        script_sha256=sha(__file__),
        B=B,
        B_CI=B_CI,
        batch_size=BATCH,
        audit_n=AUDIT_N,
        rng="numpy.random.default_rng (PCG64)",
        python=sys.version,
        numpy=np.__version__,
        pandas=pd.__version__,
        platform=platform.platform(),
        n_tests=len(result),
        invalid_draws=int(result.invalid_draws.sum()),
        estimand="Expected paired squared prediction loss difference, conditional on fixed fitted functions",
        pooled_target="Fixed observed trial weights; patient resampling within trial",
    )
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
