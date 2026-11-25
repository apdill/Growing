"""s7 (fixed h) vs s9 (h ~ 1/N): exponents in a fixed window and distance from m = 1 vs N."""
import json, os, numpy as np, pandas as pd, matplotlib, sys
matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, R); from growing.analysis import avalanches as av
t = pd.read_csv(os.path.join(R, "reports", "trunc_exponents.csv"))
fig, ax = plt.subplots(1, 4, figsize=(20, 4.4))
for i, (sw, lab) in enumerate([("s7_finite_size", "s7: h fixed (h/r* = 0.01)"), ("s9_fss_scaled_drive", "s9: N·h fixed (h/r* = 0.5/N)")]):
    g = t[t.sweep == sw].groupby("morph.num_neurons")[["tau_trunc", "alpha_trunc", "m"]].agg(["mean", "std"])
    N = g.index.values
    ax[0].errorbar(N, g["tau_trunc"]["mean"], g["tau_trunc"]["std"], marker="o", capsize=2, color=f"C{i}", label=lab)
    ax[1].errorbar(N, g["alpha_trunc"]["mean"], g["alpha_trunc"]["std"], marker="o", capsize=2, color=f"C{i}", label=lab)
    ax[2].loglog(N, 1 - g["m"]["mean"], "o-", color=f"C{i}", label=lab)
    fss = json.load(open(os.path.join(R, "reports", f"{sw}_fss.json")))
    ax[3].loglog([r["N"] for r in fss["rows"]], [r["sc"] for r in fss["rows"]], "o-", color=f"C{i}",
                 label=f"{lab}: ~N^{fss['D_sc']:.2f}")
# Galton-Watson benchmark in the same windows
rng = np.random.default_rng(0)
def gw():
    s = a = T = 1
    while a and s < 1e5:
        a = rng.poisson(a); s += a; T += 1 if a else 0
    return s, T
x = np.array([gw() for _ in range(20000)])
tb, ab = av.fit_truncated_power_law(x[:, 0], 2, 200)["alpha"], av.fit_truncated_power_law(x[:, 1], 2, 50)["alpha"]
ax[0].axhline(tb, color="r", ls="--", lw=1, label=f"critical branching, same window ({tb:.2f})")
ax[1].axhline(ab, color="r", ls="--", lw=1, label=f"critical branching, same window ({ab:.2f})")
ax[2].loglog([25, 200], [0.02, 0.0025], "k:", lw=1, label="h/r* (s9)")
for a_, ttl in zip(ax, ["size exponent τ (MLE, S ∈ [2,200])", "duration exponent α (MLE, T ∈ [2,50])",
                        "1 − m (MR estimator, plastic run)", "avalanche cutoff ⟨S²⟩/⟨S⟩"]):
    a_.set_title(ttl); a_.set_xlabel("N"); a_.set_xscale("log"); a_.legend(fontsize=7)
fig.suptitle("Approaching the SOC limit: growing networks become critical as N → ∞ only if the drive is scaled down (h ∝ 1/N)")
fig.tight_layout(); fig.savefig(os.path.join(R, "reports", "s9_vs_s7.png"), dpi=110)
