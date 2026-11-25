"""
Summarise sweep results: python scripts/analysis/analyze_sweeps.py SWEEP_NAME [--x PARAM] [--out reports/]

Produces a multi-panel figure of criticality / function metrics against the swept parameter
(mean ± sd over seeds) and prints a compact table.
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

PANELS = [
    ("frozen_branching_ratio_m", "branching ratio m (frozen)", (0, 1.05)),
    ("branching_ratio_m", "branching ratio m (plastic run)", (0, 1.05)),
    ("frozen_size_exponent_tau", "size exponent τ (frozen)", None),
    ("frozen_duration_exponent_alpha", "duration exponent α (frozen)", None),
    ("frozen_gamma_fit", "γ fit vs predicted (frozen)", None),
    ("frozen_max_size", "largest avalanche (frozen)", "log"),
    ("dynamic_range_db", "dynamic range Δ (dB)", None),
    ("rate_over_tfr", "mean rate (plastic) / tfr", None),
    ("lambda_max", "λ_max(W)", None),
    ("mean_degree", "mean degree", None),
    ("mean_length", "mean arbor radius (px)", None),
    ("clustering", "clustering coefficient", None),
]


def load(name):
    df = pd.read_csv(os.path.join(HERE, "sweeps", name, "results.csv"))
    if "error" in df:
        bad = df["error"].notna()
        if bad.any():
            print(f"{bad.sum()} failed runs:", df.loc[bad, "error"].unique()[:3])
        df = df[~bad]
    tfr = df["param.lif.tfr"] if "param.lif.tfr" in df else 0.4
    df = df.assign(tfr=tfr, rate_over_tfr=df["mean_rate"] / tfr)
    return df


def summarize(df, x, group=None):
    keys = [x] + ([group] if group else [])
    num = df.select_dtypes("number")
    cols = [c for c, *_ in PANELS if c in num] + ["frozen_gamma_predicted", "param.lif.tfr"]
    cols = [c for c in dict.fromkeys(cols) if c in df]
    return df.groupby(keys)[cols].agg(["mean", "std"])


def plot(df, x, out, title, group=None, logx=True, theory_m=None):
    fig, axes = plt.subplots(3, 4, figsize=(17, 11))
    groups = [(None, df)] if group is None else list(df.groupby(group))
    for ax, (col, label, ylim) in zip(axes.ravel(), PANELS):
        for gi, (gval, d) in enumerate(groups):
            if col not in d:
                continue
            g = d.groupby(x)[col].agg(["mean", "std"]).reset_index()
            yv = g["mean"]
            lab = None if gval is None else f"{group.split('.')[-1]}={gval}"
            ax.errorbar(g[x], yv, yerr=g["std"].fillna(0), marker="o", ms=4, capsize=2, label=lab, color=f"C{gi}")
            ax.scatter(d[x], d[col], s=8, alpha=0.3, color=f"C{gi}")
            if col == "frozen_gamma_fit" and "frozen_gamma_predicted" in d:
                gp = d.groupby(x)["frozen_gamma_predicted"].mean()
                ax.plot(gp.index, gp.values, "k:", lw=1, label="predicted (α-1)/(τ-1)" if gi == 0 else None)
        if col in ("frozen_branching_ratio_m", "branching_ratio_m") and theory_m is not None:
            xs = np.sort(df[x].unique())
            ax.plot(xs, theory_m(xs), "r--", lw=1, label="1 - h/tfr")
            ax.axhline(1.0, color="0.6", lw=0.8)
        if col == "frozen_size_exponent_tau":
            ax.axhline(1.5, color="r", ls="--", lw=0.8, label="critical 3/2")
        if col == "frozen_duration_exponent_alpha":
            ax.axhline(2.0, color="r", ls="--", lw=0.8, label="critical 2")
        if col == "rate_over_tfr":
            ax.axhline(1.0, color="0.6", lw=0.8)
        if logx:
            ax.set_xscale("log")
        if ylim == "log":
            ax.set_yscale("log")
        elif ylim:
            ax.set_ylim(*ylim)
        ax.set_title(label, fontsize=10)
        ax.set_xlabel(x.split(".")[-1])
        if ax.get_legend_handles_labels()[0]:
            ax.legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out, dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--x", default=None)
    ap.add_argument("--group", default=None)
    ap.add_argument("--linear", action="store_true")
    ap.add_argument("--theory", action="store_true", help="overlay m = 1 - ext_rate/tfr")
    ap.add_argument("--out", default=os.path.join(HERE, "reports"))
    a = ap.parse_args()
    df = load(a.name)
    with open(os.path.join(HERE, "sweeps", a.name, "spec.json")) as f:
        spec = json.load(f)
    x = a.x or ("param." + list(spec.get("grid", {"seed": 0}).keys())[0])
    os.makedirs(a.out, exist_ok=True)
    theory = None
    if a.theory:
        tfr = float(df["tfr"].iloc[0])
        theory = lambda h: np.clip(1 - np.asarray(h) / tfr, 0, 1)
    plot(df, x, os.path.join(a.out, f"{a.name}.png"), a.name, group=a.group, logx=not a.linear, theory_m=theory)
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.precision", 3)
    s = summarize(df, x, a.group)
    print(s)
    s.to_csv(os.path.join(a.out, f"{a.name}_summary.csv"))


if __name__ == "__main__":
    main()
