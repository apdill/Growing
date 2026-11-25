"""Cross-sweep view of the avalanche-body exponent (fixed window, reports/trunc_exponents.csv) against
coupling fluctuations, N and drive: python cross_sweep.py -> reports/cross_sweep_tau.png"""
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
t = pd.read_csv(os.path.join(HERE, "reports", "trunc_exponents.csv"))
parts = []
for sw in t.sweep.unique():
    r = pd.read_csv(os.path.join(HERE, "sweeps", sw, "results.csv"))[["run_id", "sigma_late_std", "sigma_late_mean", "mean_rate", "dynamic_range_db"]]
    parts.append(t[t.sweep == sw].merge(r, left_on="run", right_on="run_id", how="left"))
df = pd.concat(parts)
df = df[df["lif.spike_mode"] == "branching"]
fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
for i, (sw, d) in enumerate(df.groupby("sweep")):
    ax[0].plot(d.sigma_late_std, d.tau_trunc, "o", color=f"C{i}", label=sw, alpha=0.8)
    ax[1].plot(d.sigma_late_std, d.m, "o", color=f"C{i}", label=sw, alpha=0.8)
    ax[2].plot(d.tau_trunc, d.dynamic_range_db, "o", color=f"C{i}", label=sw, alpha=0.8)
ax[0].axhline(1.5, color="r", ls="--", lw=0.8); ax[0].set_xlabel("s.d. of σ(t) during growth"); ax[0].set_ylabel("τ, MLE on S ∈ [2, 200]")
ax[1].set_xlabel("s.d. of σ(t) during growth"); ax[1].set_ylabel("MR branching ratio m (plastic)")
ax[2].set_xlabel("τ, MLE on S ∈ [2, 200]"); ax[2].set_ylabel("dynamic range Δ (dB, final network)")
for a in ax: a.legend(fontsize=7)
fig.suptitle("Branching-model sweeps: avalanche-body exponent vs coupling fluctuations")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "reports", "cross_sweep_tau.png"), dpi=110)
print(df[["sigma_late_std", "tau_trunc", "m", "dynamic_range_db"]].corr().round(2))
