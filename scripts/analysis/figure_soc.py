"""Overlay avalanche statistics of plastic runs: python scripts/analysis/figure_soc.py OUT.png RUN_DIR:label [RUN_DIR:label ...]"""
import sys, os, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from growing.analysis import avalanches as av

out, specs = sys.argv[1], sys.argv[2:]
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
for i, sp in enumerate(specs):
    d, lab = sp.rsplit(":", 1)
    S = np.load(os.path.join(d, "spikes.npz"))["spikes"]
    r = av.analyze(S, 1, 50.0, 0.2); sm = av.summary(r)
    for ax, key in zip(axes[:2], ["sizes", "durations"]):
        c, p = av._log_binned_pdf(r[key]); ax.loglog(c, p, "o-", ms=3, lw=0.8, color=f"C{i}",
            label=f"{lab}: {'τ' if key == 'sizes' else 'α'}={sm['size_exponent_tau' if key == 'sizes' else 'duration_exponent_alpha']:.2f}")
    sc = r["scaling"]; axes[2].loglog(sc["T"], sc["meanS"], "o", ms=3, color=f"C{i}",
            label=f"{lab}: γ={sc['gamma']:.2f} (pred {r['gamma_predicted']:.2f}), m={r['branching']['m']:.3f}")
x = np.array([1, 300.0]); axes[0].loglog(x, 0.3 * x ** -1.5, "k--", lw=1, label="slope -3/2")
x = np.array([1, 100.0]); axes[1].loglog(x, 0.3 * x ** -2.0, "k--", lw=1, label="slope -2")
for ax, t in zip(axes, ["P(S) avalanche size", "P(T) avalanche duration (steps)", "<S>(T)"]):
    ax.set_title(t); ax.legend(fontsize=7)
fig.suptitle("Avalanches during self-organised growth (plastic runs, after burn-in)")
fig.tight_layout(); fig.savefig(out, dpi=110)
