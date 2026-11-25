"""
Finite-size analysis of a sweep over num_neurons: python finite_size.py SWEEP_NAME [--out reports/]

Per N (seeds pooled), from the plastic runs after burn-in:
  * P(S), P(T) and the collapse S^tau P(S) vs S / N^D, with D from the cutoff scaling
    s_c(N) = <S^2>/<S> ~ N^D
  * fluctuations of the coupling sigma(t) during growth vs N (does the network's distance to the
    critical point shrink with system size?)
  * MR branching ratio m and exponents vs N
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import avalanches as av  # noqa: E402


def load_runs(name, burn_in):
    runs = {}
    for d in sorted(glob.glob(os.path.join(HERE, "sweeps", name, "runs", "*"))):
        if not os.path.exists(os.path.join(d, "metrics.json")):
            continue
        cfg = json.load(open(os.path.join(d, "config.json")))
        N = int(cfg["morph"]["num_neurons"])
        S = np.load(os.path.join(d, "spikes.npz"))["spikes"]
        r = av.analyze(S, 1, 50.0, burn_in)
        ct = np.load(os.path.join(d, "traces.npz"))["coupling_trace"]
        late = ct[int(len(ct) * 0.5):, 1] if len(ct) else np.array([np.nan])
        runs.setdefault(N, []).append(dict(dir=d, res=r, sigma=late, ct=ct, rate=S.mean(),
                                           cond=sigma_conditioned(S, ct, burn_in)))
    return runs


def sigma_conditioned(S, ct, burn_in):
    """Split avalanches (bin 1, threshold 0) by the coupling sigma at their onset into tertiles.
    A self-organised *quasi*-critical network sweeps through sub- and supercritical states, so the
    statistics differ strongly between tertiles; a genuinely critical one would not."""
    start = int(burn_in * len(S)); A = S[start:].sum(1)
    e = np.diff(np.r_[0, (A > 0).astype(int), 0]); st = np.flatnonzero(e == 1); en = np.flatnonzero(e == -1)
    sizes = np.array([A[a:b].sum() for a, b in zip(st, en)], float)
    sig = np.interp(st + start, ct[:, 0], ct[:, 1])
    q = np.quantile(sig, [1 / 3, 2 / 3])
    out = []
    for m in (sig < q[0], (sig >= q[0]) & (sig < q[1]), sig >= q[1]):
        out.append(dict(sigma=float(sig[m].mean()), mean_S=float(sizes[m].mean()),
                        sc=float((sizes[m] ** 2).mean() / sizes[m].mean())))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("name")
    ap.add_argument("--burn-in", type=float, default=0.2)
    ap.add_argument("--out", default=os.path.join(HERE, "reports"))
    a = ap.parse_args()
    runs = load_runs(a.name, a.burn_in)
    Ns = np.array(sorted(runs))
    rows = []
    for N in Ns:
        sizes = np.concatenate([r["res"]["sizes"] for r in runs[N]])
        durs = np.concatenate([r["res"]["durations"] for r in runs[N]])
        sf = av.fit_discrete_power_law(sizes)
        df_ = av.fit_discrete_power_law(durs)
        sig = np.concatenate([r["sigma"] for r in runs[N]])
        sig_std = np.mean([np.std(r["sigma"]) for r in runs[N]])
        rows.append(dict(N=int(N), n_av=len(sizes), sc=float((sizes.astype(float) ** 2).mean() / sizes.mean()),
                         s999=float(np.percentile(sizes, 99.9)), smax=int(sizes.max()),
                         tau=sf["alpha"], tau_xmin=sf["xmin"], alpha=df_["alpha"],
                         m=float(np.mean([r["res"]["branching"]["m"] for r in runs[N]])),
                         sigma_mean=float(sig.mean()), sigma_std=float(sig_std),
                         sigma_cv=float(sig_std / sig.mean()),
                         rate_per_step=float(np.mean([r["rate"] for r in runs[N]])),
                         meanS_high_over_low=float(np.mean([r["cond"][2]["mean_S"] / r["cond"][0]["mean_S"] for r in runs[N]])),
                         sigma_low=float(np.mean([r["cond"][0]["sigma"] for r in runs[N]])),
                         sigma_high=float(np.mean([r["cond"][2]["sigma"] for r in runs[N]])),
                         sizes=sizes, durs=durs))
    logN = np.log(Ns)
    D = np.polyfit(logN, np.log([r["sc"] for r in rows]), 1)[0]
    D999 = np.polyfit(logN, np.log([r["s999"] for r in rows]), 1)[0]
    b_sig = np.polyfit(logN, np.log([r["sigma_std"] for r in rows]), 1)[0]
    b_m = np.polyfit(logN, np.log([1 - r["m"] for r in rows]), 1)[0]
    tau_ref = 1.5

    fig, ax = plt.subplots(2, 3, figsize=(17, 9.5))
    for i, r in enumerate(rows):
        c, p = av._log_binned_pdf(r["sizes"])
        ax[0, 0].loglog(c, p, "o-", ms=3, lw=0.8, color=f"C{i}", label=f"N={r['N']}: τ={r['tau']:.2f}")
        ax[0, 1].loglog(c / r["N"] ** D, c ** tau_ref * p, "o-", ms=3, lw=0.8, color=f"C{i}", label=f"N={r['N']}")
        c, p = av._log_binned_pdf(r["durs"])
        ax[0, 2].loglog(c, p, "o-", ms=3, lw=0.8, color=f"C{i}", label=f"N={r['N']}: α={r['alpha']:.2f}")
    x = np.array([1, 1000.0]); ax[0, 0].loglog(x, 0.3 * x ** -1.5, "k--", lw=1, label="−3/2")
    x = np.array([1, 200.0]); ax[0, 2].loglog(x, 0.3 * x ** -2.0, "k--", lw=1, label="−2")
    ax[0, 0].set_title("P(S), plastic runs (seeds pooled)")
    ax[0, 1].set_title(f"collapse: S^{tau_ref} P(S) vs S / N^D,  D = {D:.2f}")
    ax[0, 2].set_title("P(T)")
    ax[1, 0].loglog(Ns, [r["sc"] for r in rows], "o-", label=f"<S²>/<S> ~ N^{D:.2f}")
    ax[1, 0].loglog(Ns, [r["s999"] for r in rows], "s-", label=f"99.9th pct ~ N^{D999:.2f}")
    ax[1, 0].set_title("avalanche cutoff vs N"); ax[1, 0].set_xlabel("N")
    ax[1, 1].loglog(Ns, [r["sigma_std"] for r in rows], "o-", label=f"s.d. σ(t) ~ N^{b_sig:.2f}")
    ax[1, 1].loglog(Ns, [1 - r["m"] for r in rows], "s-", label=f"1 − m ~ N^{b_m:.2f}")
    ax[1, 1].loglog(Ns, rows[0]["sigma_std"] * (Ns / Ns[0]) ** -0.5, "k:", lw=1, label="N^−1/2")
    ax[1, 1].loglog(Ns, [r["meanS_high_over_low"] / 100 for r in rows], "^-",
                    label="<S>(σ top third) / <S>(σ bottom third) ÷ 100")
    ax[1, 1].set_title("distance from / fluctuation around the critical point"); ax[1, 1].set_xlabel("N")
    for i, N in enumerate(Ns):
        for r in runs[N][:1]:
            ct = r["ct"]
            ax[1, 2].plot(ct[:, 0] / ct[-1, 0], ct[:, 1], lw=0.6, color=f"C{i}", label=f"N={N}")
    ax[1, 2].set_title("σ(t) during growth (seed 1)"); ax[1, 2].set_xlabel("fraction of run")
    for axx in ax.ravel():
        if axx.get_legend_handles_labels()[0]:
            axx.legend(fontsize=7)
    fig.suptitle(f"{a.name}: finite-size scaling of avalanches in self-organised growing networks")
    fig.tight_layout()
    out = os.path.join(a.out, f"{a.name}_fss.png")
    fig.savefig(out, dpi=110)
    table = [{k: v for k, v in r.items() if k not in ("sizes", "durs")} for r in rows]
    json.dump(dict(rows=table, D_sc=D, D_999=D999, sigma_std_exp=b_sig, one_minus_m_exp=b_m),
              open(os.path.splitext(out)[0] + ".json", "w"), indent=1)
    for t in table:
        print("  ".join(f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}" for k, v in t.items()))
    print(f"D(<S2>/<S>)={D:.3f}  D(99.9)={D999:.3f}  sd(sigma)~N^{b_sig:.3f}  (1-m)~N^{b_m:.3f}")


if __name__ == "__main__":
    main()
