"""
F1 analysis: morphology knobs -> measured geometry, and selection of settings for F2/F3.

    python scripts/analysis/fractal_calibration.py [SWEEP=f1_morph_calibration] [--targets 1.2 1.35 1.5 1.65]
                                  [--metric D_neuron_mean] [--cable-tol 0.25]

Writes reports/<sweep>_calibration.png and prints the chosen settings as JSON "runs" entries to
paste into f2a/f2b/f3 specs (replacing the FROM_F1 placeholders). Settings are chosen so the measured
D hits each target while cable length stays within ±cable-tol of the median of the candidates
(geometry varied at roughly matched wiring material).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repository root
sys.path.insert(0, HERE)
KNOBS = ["param.branch.trajectory_deviation", "param.branch.branch_angle", "param.branch.D",
         "param.branch.branch_length_scale", "param.branch.randomness"]
GEOM = ["D_neuron_mean", "D_neuron_skel_mean", "D_network_fill", "lacunarity_neuron_16", "tortuosity_mean",
        "branches_per_neuron", "cable_length_mean", "mean_length"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name", nargs="?", default="f1_morph_calibration")
    ap.add_argument("--targets", type=float, nargs="+", default=[1.2, 1.35, 1.5, 1.65])
    ap.add_argument("--metric", default="D_neuron_mean")
    ap.add_argument("--cable-tol", type=float, default=0.25)
    a = ap.parse_args()
    df = pd.read_csv(os.path.join(HERE, "sweeps", a.name, "results.csv"))
    if "error" in df:
        df = df[df["error"].isna()]
    knobs = [k for k in KNOBS if k in df and df[k].nunique() > 1]
    geom = [g for g in GEOM if g in df]

    fig, axes = plt.subplots(len(geom), len(knobs), figsize=(3.2 * len(knobs), 2.4 * len(geom)), squeeze=False)
    for i, g in enumerate(geom):
        for j, k in enumerate(knobs):
            s = df.groupby(k)[g].agg(["mean", "std"])
            axes[i, j].errorbar(s.index, s["mean"], s["std"], marker="o", capsize=2)
            if i == len(geom) - 1:
                axes[i, j].set_xlabel(k.split(".")[-1])
            if j == 0:
                axes[i, j].set_ylabel(g, fontsize=8)
    fig.suptitle(f"{a.name}: morphology knobs → measured geometry (mean ± sd over the other knobs)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "reports", f"{a.name}_calibration.png"), dpi=90)

    print("correlation of geometry with knobs:")
    print(df[knobs + geom].corr().loc[geom, knobs].round(2).to_string())
    print(f"\nrange of {a.metric}: {df[a.metric].min():.3f} – {df[a.metric].max():.3f}")

    med = df["cable_length_mean"].median()
    tol = a.cable_tol
    pool = df[(df["cable_length_mean"] - med).abs() <= tol * med]
    while len(pool) < 2 * len(a.targets) and tol < 10:      # widen until enough candidates
        tol *= 1.5
        pool = df[(df["cable_length_mean"] - med).abs() <= tol * med]
    print(f"{len(pool)}/{len(df)} settings within ±{tol:.0%} of median cable length {med:.0f} px")
    chosen = []
    for t in a.targets:
        if pool.empty:
            break
        row = pool.iloc[(pool[a.metric] - t).abs().argsort().iloc[0]]
        chosen.append({k.replace("param.", ""): float(row[k]) for k in knobs}
                      | {"_measured": f"{a.metric}={row[a.metric]:.3f}, cable={row['cable_length_mean']:.0f}, "
                                      f"tort={row.get('tortuosity_mean', np.nan):.2f}"})
    print("\nsuggested 'runs' entries for F2/F3 specs:")
    print(json.dumps(chosen, indent=1))


if __name__ == "__main__":
    main()
