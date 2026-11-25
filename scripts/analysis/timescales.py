"""
Growth-timescale analysis of a sweep (python scripts/analysis/timescales.py SWEEP [--key param] [--group param]):
pooled P(S) per value of `key`, a fixed-range slope fit (S in [2, 200]), and the sigma-conditioned
(quasi-criticality) test from finite_size.py.
"""
from __future__ import annotations
import argparse, glob, json, os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repository root
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from growing.analysis import avalanches as av  # noqa: E402
from finite_size import sigma_conditioned  # noqa: E402


def get(cfg, dotted):
    sec, f = dotted.split(".")
    return cfg[sec][f]


def slope(sizes, lo=2, hi=200):
    c, p = av._log_binned_pdf(sizes)
    m = (c >= lo) & (c <= hi) & (p > 0)
    return -np.polyfit(np.log(c[m]), np.log(p[m]), 1)[0] if m.sum() > 3 else np.nan


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("name"); ap.add_argument("--key", default="lif.length_growth_rate")
    ap.add_argument("--group", default="lif.rate_tau"); ap.add_argument("--burn-in", type=float, default=0.2)
    a = ap.parse_args()
    data = {}
    for d in sorted(glob.glob(os.path.join(HERE, "sweeps", a.name, "runs", "*"))):
        if not os.path.exists(os.path.join(d, "metrics.json")):
            continue
        cfg = json.load(open(os.path.join(d, "config.json")))
        S = np.load(os.path.join(d, "spikes.npz"))["spikes"]; ct = np.load(os.path.join(d, "traces.npz"))["coupling_trace"]
        r = av.analyze(S, 1, 50.0, a.burn_in)
        late = ct[len(ct) // 2:, 1]
        data.setdefault((get(cfg, a.key), get(cfg, a.group)), []).append(dict(
            sizes=r["sizes"], durs=r["durations"], m=r["branching"]["m"], cond=sigma_conditioned(S, ct, a.burn_in),
            sd=float(np.std(late)), rate=float(S.mean() / cfg["lif"]["dt"])))
    keys = sorted({k for k, _ in data}); groups = sorted({g for _, g in data})
    fig, ax = plt.subplots(1, 3, figsize=(17, 4.6)); rows = []
    for i, k in enumerate(keys):
        for j, g in enumerate(groups):
            rs = data.get((k, g), [])
            if not rs:
                continue
            s = np.concatenate([r["sizes"] for r in rs]); T = np.concatenate([r["durs"] for r in rs])
            row = dict(key=k, group=g, slope_S=slope(s), slope_T=slope(T, 2, 60), sigma_sd=np.mean([r["sd"] for r in rs]),
                       m=np.mean([r["m"] for r in rs]), rate=np.mean([r["rate"] for r in rs]),
                       hi_lo=np.mean([r["cond"][2]["mean_S"] / r["cond"][0]["mean_S"] for r in rs]),
                       sc=float((s.astype(float) ** 2).mean() / s.mean()))
            rows.append(row)
            c, p = av._log_binned_pdf(s)
            ax[j].loglog(c, p, "o-", ms=3, lw=0.8, color=f"C{i}",
                         label=f"{a.key.split('.')[-1]}={k}: slope {row['slope_S']:.2f}, sd σ={row['sigma_sd']:.2f}")
    for j, g in enumerate(groups):
        x = np.array([1, 1000.0]); ax[j].loglog(x, 0.3 * x ** -1.5, "k--", lw=1, label="−3/2")
        ax[j].set_title(f"P(S) during growth, {a.group.split('.')[-1]}={g}"); ax[j].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "reports", f"{a.name}_timescales.png"), dpi=110)
    for r in rows:
        print("  ".join(f"{k}={v:.3g}" for k, v in r.items()))
    json.dump(rows, open(os.path.join(HERE, "reports", f"{a.name}_timescales.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
