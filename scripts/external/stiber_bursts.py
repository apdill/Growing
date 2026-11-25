"""
Avalanche-style analysis of a Stiber-lab simulated culture (Zenodo 10.5281/zenodo.17220617;
Kawasaki & Stiber 2014 model; Stiber, Gonzales & Lee 2026, arXiv:2602.22364).

The deposit contains only *burst epochs*: runs of 10 ms bins with > 50 network spikes
(allBinnedBurstInfo.csv) and, per burst, per-neuron spike counts in 10 ms frames (allFrames.mat,
MATLAB v7.3 cell array, 9160 × [bins × 10000]). Inter-burst activity is not included, so a
low-threshold avalanche analysis is impossible; a burst is exactly a thresholded avalanche with
θ = 50 spikes / 10 ms (Poil et al. 2012 style). What we can measure:
  * P(S), P(T) of bursts, fixed-window and KS power-law fits, size–duration scaling γ
  * how these evolve over development (quartiles of simulated time)
  * within-burst amplification: ratio of successive frame counts A(t+1)/A(t) in the rising phase
    (an upper bound on the branching ratio during propagation) and the peak-normalised profile
  * participation: fraction of the 10 000 neurons that spike in each burst

usage: python scripts/external/stiber_bursts.py external/stiber_2025/tR_1.0--fE_0.98 [--frames]
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
from growing.analysis import avalanches as av  # noqa: E402

BIN_S = 0.010
N_NEURONS = 10_000
DIV_TOTAL = 28.0           # simulated development spans the equivalent of 28 days in vitro


def frame_stats(mat_path, max_bursts=None):
    import h5py
    rows = []
    with h5py.File(mat_path, "r") as f:
        refs = f["allFrames"][:, 0]
        for k, ref in enumerate(refs[:max_bursts]):
            X = np.asarray(f[ref], dtype=np.int32)          # (bins, neurons)
            A = X.sum(1)                                     # network spikes per 10 ms frame
            pk = int(np.argmax(A))
            rise = A[:pk + 1]
            ratios = rise[1:] / np.maximum(rise[:-1], 1)
            rows.append(dict(burst=k, frames=len(A), frame_total=int(A.sum()), peak=int(A.max()), peak_idx=pk,
                             participation=float((X.sum(0) > 0).mean()),
                             spikes_per_participant=float(A.sum() / max((X.sum(0) > 0).sum(), 1)),
                             rise_ratio_med=float(np.median(ratios)) if len(ratios) else np.nan,
                             rise_ratio_first=float(ratios[0]) if len(ratios) else np.nan,
                             A=A.tolist()))
    return rows


def fits(S, T):
    sf, tf = av.fit_discrete_power_law(S), av.fit_discrete_power_law(T)
    sw = av.fit_truncated_power_law(S, int(np.percentile(S, 1)), int(np.percentile(S, 99)))
    tw = av.fit_truncated_power_law(T, int(max(np.percentile(T, 1), 1)), int(np.percentile(T, 99)))
    sc = av.size_vs_duration_exponent(S, T)
    return dict(n=len(S), S_mean=float(S.mean()), S_cv=float(S.std() / S.mean()), T_mean=float(T.mean()),
                T_cv=float(T.std() / T.mean()), tau_ks=sf["alpha"], tau_ks_xmin=sf["xmin"], alpha_ks=tf["alpha"],
                tau_window=sw["alpha"], alpha_window=tw["alpha"], gamma=sc["gamma"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sim_dir")
    ap.add_argument("--frames", action="store_true", help="also read allFrames.mat (slow, ~4 GB streamed)")
    ap.add_argument("--max-bursts", type=int, default=None)
    a = ap.parse_args()
    name = os.path.basename(a.sim_dir.rstrip("/"))
    out = os.path.join(HERE, "reports", f"stiber_{name}")
    df = pd.read_csv(os.path.join(a.sim_dir, "allBinnedBurstInfo.csv"))
    df.columns = [c.strip() for c in df.columns]
    S = df["totalSpikeCount"].to_numpy(int)
    T = df["width(bins)"].to_numpy(int)
    t_bin = df["startBin#"].to_numpy(float)
    ibi = np.diff(t_bin) * BIN_S
    t_end = t_bin.max()
    div = t_bin / t_end * DIV_TOTAL                     # assumes the run ends at 28 DIV
    q = np.digitize(div, [7, 14, 21])                    # weeks 1-4

    summary = dict(sim=name, all=fits(S, T),
                   by_week={f"week{w + 1}": fits(S[q == w], T[q == w]) for w in range(4) if (q == w).sum() > 30},
                   first_burst_div=float(div.min()), n_bursts=int(len(S)),
                   ibi_s_by_week={f"week{w + 1}": float(np.median(ibi[q[1:] == w])) for w in range(4) if (q[1:] == w).sum() > 5},
                   S_over_N_median=float(np.median(S) / N_NEURONS))
    fr = None
    if a.frames:
        fr = pd.DataFrame(frame_stats(os.path.join(a.sim_dir, "allFrames.mat"), a.max_bursts))
        fr["week"] = q[fr["burst"].to_numpy()] + 1
        summary["frames"] = {f"week{w}": dict(participation=float(g.participation.median()),
                                              spikes_per_participant=float(g.spikes_per_participant.median()),
                                              rise_ratio_med=float(g.rise_ratio_med.median()),
                                              rise_ratio_first=float(g.rise_ratio_first.median()))
                             for w, g in fr.groupby("week")}

    fig, ax = plt.subplots(2, 3, figsize=(17, 9))
    for w in range(4):
        m = q == w
        if m.sum() < 30:
            continue
        for axx, x in ((ax[0, 0], S[m]), (ax[0, 1], T[m])):
            c, p = av._log_binned_pdf(x, 30)
            axx.loglog(c, p, "o-", ms=3, lw=0.8, color=f"C{w}", label=f"week {w + 1} (n={m.sum()})")
    x = np.array([1e3, 1e5]); ax[0, 0].loglog(x, 1e2 * x ** -1.5, "k--", lw=1, label="−3/2")
    ax[0, 0].set_title("burst size S (spikes) — thresholded avalanches θ = 50 / 10 ms"); ax[0, 0].legend(fontsize=7)
    ax[0, 1].set_title("burst duration T (10 ms bins)"); ax[0, 1].legend(fontsize=7)
    ax[0, 2].loglog(T, S, ".", ms=1, alpha=0.3)
    ax[0, 2].set_title(f"S vs T, γ = {summary['all']['gamma']:.2f}"); ax[0, 2].set_xlabel("T"); ax[0, 2].set_ylabel("S")
    ax[1, 0].plot(div, S, ".", ms=1, alpha=0.4); ax[1, 0].set_xlabel("DIV (equivalent)"); ax[1, 0].set_ylabel("S")
    ax[1, 0].set_title("burst size over development")
    ax[1, 1].semilogy(div[1:], ibi, ".", ms=1, alpha=0.4); ax[1, 1].set_xlabel("DIV (equivalent)")
    ax[1, 1].set_ylabel("inter-burst interval (s)"); ax[1, 1].set_title("inter-burst interval")
    if fr is not None:
        for w, g in fr.groupby("week"):
            prof = [np.interp(np.linspace(0, 1, 40), np.linspace(0, 1, len(A)), np.asarray(A) / max(A)) for A in g.A[:500]]
            ax[1, 2].plot(np.linspace(0, 1, 40), np.mean(prof, 0), color=f"C{w - 1}",
                          label=f"week {w}: rise ×{g.rise_ratio_med.median():.2f}/10 ms, participation {g.participation.median():.0%}")
        ax[1, 2].set_title("mean burst profile (normalised)"); ax[1, 2].legend(fontsize=7)
    fig.suptitle(f"Stiber-lab simulated culture {name}: bursts as thresholded avalanches")
    fig.tight_layout(); fig.savefig(out + ".png", dpi=100)
    json.dump(summary, open(out + ".json", "w"), indent=1, default=float)
    print(json.dumps(summary, indent=1, default=lambda v: round(float(v), 4)))


if __name__ == "__main__":
    main()
