"""
Coupling scan around a grown network: python scripts/sweeps/scan_frozen.py RUN_DIR [RUN_DIR ...] [--out FIG]

For each run, freeze the grown connectivity M, multiply it by g over a range around 1, and
measure (frozen LIF/branching dynamics, same model as the run):
  * sigma(g)  = mean column sum of min(1, g * spike_input * M)  (expected offspring per spike)
  * MR branching ratio, mean rate
  * avalanche size / duration exponents
  * dynamic range Delta(g)
This shows where the self-organised network (g = 1) sits relative to the critical point and
to the dynamic-range optimum.
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import sys
from types import SimpleNamespace

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from growing.analysis import avalanches, probes  # noqa: E402

G_DEFAULT = [0.5, 0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.4]


def load_run(run_dir):
    cfg = json.load(open(os.path.join(run_dir, "config.json")))
    lif = SimpleNamespace(**cfg["lif"])
    M = np.load(os.path.join(run_dir, "M.npy"))
    return lif, M


def sigma_of(M, lif, g=1.0):
    if getattr(lif, "spike_mode", "threshold") == "branching":
        P = np.clip(g * lif.spike_input * M, 0, 1)
        return float(P.sum(0).mean())
    return float(np.max(np.abs(np.linalg.eigvals(g * lif.spike_input * M))))


def lambda_of(M, lif, g=1.0):
    """Largest |eigenvalue| of the effective transmission matrix at gain g."""
    W = g * lif.spike_input * np.asarray(M, float)
    if getattr(lif, "spike_mode", "threshold") == "branching":
        W = np.clip(W, 0, 1)
    return float(np.max(np.abs(np.linalg.eigvals(W))))


def gain_for_lambda(M, lif, target, g_hi=64.0, tol=1e-4):
    """Gain g with lambda_of(M, lif, g) = target (bisection; Λ is non-decreasing in g)."""
    lo, hi = 0.0, 1.0
    while lambda_of(M, lif, hi) < target and hi < g_hi:
        lo, hi = hi, hi * 2
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if lambda_of(M, lif, mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def scan(run_dir, gs=G_DEFAULT, steps=300_000, seed=0, ext_rate=None, lambda_targets=None):
    """Coupling scan. With lambda_targets, the gains are chosen so that the largest eigenvalue of
    the transmission matrix takes the given values (compare networks at matched Λ)."""
    lif, M = load_run(run_dir)
    if ext_rate is not None:          # weaker drive -> better separated avalanches (h -> 0 limit)
        lif.ext_rate = ext_rate
    if lambda_targets is not None:
        gs = [gain_for_lambda(M, lif, t) for t in lambda_targets]
    rows = []
    for g in gs:
        rng = np.random.default_rng(seed)
        S = probes.lif_run(g * M, lif, steps, rng=rng)
        a = avalanches.analyze(S, 1, 50.0, 0.02)
        sm = avalanches.summary(a)
        dr = probes.dynamic_range(g * M, lif, seed=seed)
        rows.append(dict(g=g, sigma=sigma_of(M, lif, g), lambda_max=lambda_of(M, lif, g), rate=float(S.mean() / lif.dt), m=sm["branching_ratio_m"],
                         tau=sm["size_exponent_tau"], alpha=sm["duration_exponent_alpha"], gamma=sm["gamma_fit"],
                         max_size=sm["max_size"], n_av=sm["n_avalanches"], dr=dr["delta_db"]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+")
    ap.add_argument("--out", default=None)
    ap.add_argument("--steps", type=int, default=300_000)
    ap.add_argument("--ext-rate", type=float, default=None, help="override the external event rate h")
    ap.add_argument("--g", type=float, nargs="+", default=G_DEFAULT)
    ap.add_argument("--match-lambda", type=float, nargs="+", default=None, metavar="LAMBDA",
                    help="instead of --g, scale each network so the largest eigenvalue of the "
                         "transmission matrix equals these values (geometry compared at matched gain)")
    a = ap.parse_args()
    allrows = {}
    for r in a.runs:
        rows = scan(r, gs=a.g, steps=a.steps, ext_rate=a.ext_rate, lambda_targets=a.match_lambda)
        allrows[os.path.basename(r.rstrip("/"))] = rows
        print(os.path.basename(r.rstrip("/"))[:100])
        for x in rows:
            print("  " + "  ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in x.items()))
    if a.out:
        fig, axes = plt.subplots(1, 4, figsize=(18, 4))
        for i, (name, rows) in enumerate(allrows.items()):
            sig = [x["sigma"] for x in rows]
            for ax, key, lab in zip(axes, ["dr", "tau", "alpha", "rate"],
                                    ["dynamic range Δ (dB)", "size exponent τ", "duration exponent α", "rate"]):
                ax.plot(sig, [x[key] for x in rows], "o-", ms=3, color=f"C{i % 10}", alpha=0.8)
                s1 = min(rows, key=lambda x: abs(x["g"] - 1.0))
                ax.plot(s1["sigma"], s1[key], "k*", ms=12)
                ax.set_xlabel("σ (offspring per spike)"); ax.set_title(lab)
                ax.axvline(1.0, color="0.6", lw=0.8)
        axes[1].axhline(1.5, color="r", ls="--", lw=0.8); axes[2].axhline(2.0, color="r", ls="--", lw=0.8)
        axes[3].set_yscale("log")
        fig.suptitle("Coupling scan around grown networks (★ = self-organised network, g = 1)")
        fig.tight_layout(); fig.savefig(a.out, dpi=110)
        json.dump(allrows, open(os.path.splitext(a.out)[0] + ".json", "w"), indent=1)


if __name__ == "__main__":
    main()
