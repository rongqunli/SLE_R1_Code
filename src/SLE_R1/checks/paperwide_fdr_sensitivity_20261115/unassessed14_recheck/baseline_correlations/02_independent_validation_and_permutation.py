from __future__ import annotations
import os as _sle_runtime_os

if not _sle_runtime_os.environ.get("SLE_RUN_ROOT"):
    raise RuntimeError("Run this analysis through run.py with an external workspace")
import hashlib
import json
import platform
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import stats

OUT = Path(__file__).resolve().parent
SOURCE = Path("F:/program/Hei/SLE_v1/results/reviewer_response/07_bcell_abundance")
INPUTS = [
    SOURCE / "sample_level_b_abundance_and_program_scores.tsv",
    SOURCE / "b_program_abundance_correlations.tsv",
    Path("F:/program/Hei/SLE_v1/analysis/reviewer_bcell_abundance_analysis.R"),
]
PERMUTATIONS = 99999
SEED = 20261115
TRIALS = ("ILLUMINATE-1", "ILLUMINATE-2")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    inputs_before = {str(p): sha256(p) for p in INPUTS}
    samples = pd.read_csv(INPUTS[0], sep="\t", float_precision="round_trip")
    old = pd.read_csv(INPUTS[1], sep="\t", float_precision="round_trip")
    r_results = pd.read_csv(OUT / "baseline_correlation_R_recheck.tsv", sep="\t")
    pairs_r = pd.read_csv(OUT / "baseline_analysis_pairs.tsv", sep="\t")
    records = []
    child_seeds = np.random.SeedSequence(SEED).spawn(len(TRIALS))
    for trial, child_seed in zip(TRIALS, child_seeds):
        baseline = samples.loc[(samples.trial == trial) & (samples.time == "baseline")].copy()
        assert baseline.subject_id.notna().all()
        assert baseline.subject_id.is_unique and baseline.sample_id.is_unique
        keep = np.isfinite(baseline.official_b_score) & np.isfinite(baseline.quantiseq_b_fraction_z)
        data = baseline.loc[keep].copy()
        x = data.official_b_score.to_numpy()
        y = data.quantiseq_b_fraction_z.to_numpy()
        n = len(x)
        rx = stats.rankdata(x)
        ry = stats.rankdata(y)
        assert len(np.unique(rx)) == n and len(np.unique(ry)) == n
        rx -= rx.mean()
        ry -= ry.mean()
        denominator = np.linalg.norm(rx) * np.linalg.norm(ry)
        rho = np.dot(rx, ry) / denominator
        t_statistic = rho * np.sqrt((n - 2) / (1 - rho * rho))
        p = 2 * stats.t.sf(abs(t_statistic), n - 2)
        scipy_test = stats.spearmanr(x, y)
        r = r_results.loc[r_results.trial == trial].iloc[0]
        original = old.loc[
            (old.trial == trial)
            & (old.subset == "baseline")
            & (old.proxy == "quantiseq_b_fraction_z")
            & (old.method == "spearman")
        ].iloc[0]
        assert n == int(original["n"]) == int(r.n_complete)
        assert abs(rho - original.correlation) < 1e-13
        assert abs(rho - r.rho_reproduced) < 1e-13
        assert abs(np.log(p) - np.log(r.p_primary)) < 1e-10
        assert abs(np.log(p) - np.log(scipy_test.pvalue)) < 1e-10
        rp = pairs_r.loc[pairs_r.trial == trial]
        assert list(rp.sample_id) == list(data.sample_id)
        assert list(rp.subject_id) == list(data.subject_id)
        assert np.array_equal(rp.rank_b_score.to_numpy(), rx + (n + 1) / 2)
        assert np.array_equal(rp.rank_quantiseq.to_numpy(), ry + (n + 1) / 2)
        rng = np.random.Generator(np.random.PCG64(child_seed))
        permuted_rhos = np.empty(PERMUTATIONS, dtype=np.float64)
        for start in range(0, PERMUTATIONS, 512):
            stop = min(start + 512, PERMUTATIONS)
            permuted_y = rng.permuted(np.broadcast_to(ry, (stop - start, n)), axis=1)
            permuted_rhos[start:stop] = (permuted_y * rx).sum(axis=1) / denominator
        exceedances = int(np.count_nonzero(np.abs(permuted_rhos) >= abs(rho) - 1e-12))
        p_permutation = (exceedances + 1) / (PERMUTATIONS + 1)
        np.savez_compressed(OUT / f"permutation_rhos_{trial}.npz", rho=permuted_rhos)
        records.append(
            {
                "hypothesis_id": r.hypothesis_id,
                "trial": trial,
                "n_baseline_rows": len(baseline),
                "n": n,
                "n_unique_subjects": data.subject_id.nunique(),
                "n_excluded_nonfinite": int((~keep).sum()),
                "rho": rho,
                "p_original": original.p_value,
                "p_R_default_reproduced": r.p_default_reproduced,
                "p_primary": p,
                "p_primary_method": "two_sided_Spearman_asymptotic_t_exact_FALSE_df_n_minus_2",
                "p_scipy_spearmanr": scipy_test.pvalue,
                "t_statistic": t_statistic,
                "df": n - 2,
                "permutations": PERMUTATIONS,
                "permutation_exceedances": exceedances,
                "p_permutation_plus_one": p_permutation,
                "permutation_resolution": 1 / (PERMUTATIONS + 1),
                "max_abs_permuted_rho": np.abs(permuted_rhos).max(),
                "permutation_root_seed": SEED,
                "permutation_spawn_key": str(child_seed.spawn_key),
                "default_zero_cause": r.default_zero_cause,
                "sample_identity_and_ranks_match_R": True,
            }
        )
    result = pd.DataFrame(records)
    result.to_csv(
        OUT / "baseline_correlation_final_results.tsv", sep="\t", index=False, float_format="%.17g"
    )
    inputs_after = {str(p): sha256(p) for p in INPUTS}
    assert inputs_before == inputs_after, "An input source changed during the analysis."
    metadata = {
        "status": "PASS",
        "input_sha256_before": inputs_before,
        "input_sha256_after": inputs_after,
        "inputs_unchanged": inputs_before == inputs_after,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "rng": "numpy.Generator(PCG64(SeedSequence(20261115).spawn(2)[trial_index]))",
        "trial_order": list(TRIALS),
        "batch_size": 512,
        "permutation_null": "exchangeability of baseline patient score pairing within each trial",
        "interpretation": "observational association; no conditional or causal claim",
        "analysis_plan_sha256": sha256(Path(__file__)),
        "scripts_sha256": {
            str(p): sha256(p) for p in (OUT / "01_reproduce_and_diagnose.R", Path(__file__))
        },
    }
    (OUT / "validation_and_provenance.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )
    print(result.to_string(index=False))
    print("\nPASS: independent n/rho/p/sample-rank validation and all source integrity checks.")


if __name__ == "__main__":
    main()
