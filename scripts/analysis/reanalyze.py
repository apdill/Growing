"""Recompute fixed-window exponents (avalanches.fit_truncated_power_law) for all runs of the given
sweeps, plastic and frozen: python scripts/analysis/reanalyze.py s6_soc_long s7_finite_size ... -> reports/trunc_exponents.csv"""
import glob, json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # repository root
sys.path.insert(0, HERE)
from growing.analysis import avalanches as av


def flat(cfg):
    return {f"{s}.{k}": v for s, d in cfg.items() if isinstance(d, dict) for k, v in d.items()}


rows = []
for name in sys.argv[1:]:
    for d in sorted(glob.glob(os.path.join(HERE, "sweeps", name, "runs", "*"))):
        if not os.path.exists(os.path.join(d, "spikes.npz")):
            continue
        cfg = flat(json.load(open(os.path.join(d, "config.json"))))
        met = json.load(open(os.path.join(d, "metrics.json")))
        S = np.load(os.path.join(d, "spikes.npz"))["spikes"]
        r = av.analyze(S, 1, 50.0, 0.2)
        row = dict(sweep=name, run=os.path.basename(d), seed=met.get("seed"),
                   tau_trunc=r["size_fit_trunc"]["alpha"], tau_trunc_se=r["size_fit_trunc"]["alpha_se"],
                   alpha_trunc=r["duration_fit_trunc"]["alpha"], tau_ks=r["size_fit"]["alpha"],
                   xmin_ks=r["size_fit"]["xmin"], m=r["branching"]["m"], n_av=len(r["sizes"]),
                   gamma_fit=r["scaling"]["gamma"])
        fz = os.path.join(d, "frozen_avalanches.npz")
        if os.path.exists(fz):
            z = np.load(fz)
            row.update(frozen_tau_trunc=av.fit_truncated_power_law(z["sizes"], 2, 200)["alpha"],
                       frozen_alpha_trunc=av.fit_truncated_power_law(z["durations"], 2, 50)["alpha"])
        for k in ("lif.ext_rate", "lif.tfr", "morph.num_neurons", "lif.length_growth_rate", "lif.rate_tau",
                  "branch.D", "branch.branch_length_scale", "lif.spike_mode", "branch.pass_through_probability",
                  "morph.contact_ref_px", "lif.refractory_steps"):
            row[k] = cfg.get(k)
        rows.append(row)
        print(name, row["run"][:40], {k: round(v, 3) if isinstance(v, float) else v for k, v in row.items() if "trunc" in k})
df = pd.DataFrame(rows)
out = os.path.join(HERE, "reports", "trunc_exponents.csv")
if os.path.exists(out):
    old = pd.read_csv(out); df = pd.concat([old[~old["sweep"].isin(df["sweep"].unique())], df])
df.to_csv(out, index=False)
