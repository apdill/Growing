"""
Headless parameter sweeps of the branch-morphology plasticity simulation.

Usage
-----
    python scripts/sweeps/sweep.py specs/SPEC.json [--workers 12] [--dry-run]

SPEC.json:
    {
      "name": "drive_vs_coupling",
      "base":  {"lif.duration": 7500, "morph.num_neurons": 50},          # overrides for every run
      "grid":  {"lif.input_noise": [0.1, 0.2], "lif.spike_input": [2, 4]},  # Cartesian product
      "runs":  [{"lif.tfr": 0.2}, {"lif.tfr": 0.4}],                       # optional explicit list
      "seeds": [1, 2]
    }

Keys are "<cfg>.<field>" with cfg in {morph, branch, lif, analysis, axon}. Each run writes
sweeps/<name>/runs/<run_id>/{config.json, metrics.json, spikes.npz, M.npy, traces.npz,
morphology.png, avalanches.png}; all metrics are collected in sweeps/<name>/results.csv.
Runs already finished (metrics.json present) are skipped, so a sweep can be resumed.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import random
import sys
import time
import traceback
from multiprocessing import get_context

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from growing.paths import SWEEPS_DIR  # noqa: E402

FROZEN_STEPS_DEFAULT = 500_000


def _jsonable(v):
    if isinstance(v, (np.floating, float)):
        return None if not np.isfinite(v) else float(v)
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, np.ndarray):
        return [_jsonable(x) for x in v.tolist()]
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in v.items()}
    return v


def expand(spec):
    base = spec.get("base", {})
    runs = []
    grid = spec.get("grid", {})
    explicit = spec.get("runs", [{}])
    keys = list(grid.keys())
    combos = list(itertools.product(*[grid[k] for k in keys])) if keys else [()]
    for combo in combos:
        for extra in explicit:
            for seed in spec.get("seeds", [1]):
                ov = dict(base)
                ov.update(dict(zip(keys, combo)))
                ov.update(extra)
                ov["seed"] = seed
                runs.append(ov)
    return runs


def run_id(ov):
    parts = []
    for k in sorted(ov):
        v = ov[k]
        short = k.split(".")[-1]
        parts.append(f"{short}={v:g}" if isinstance(v, (int, float)) and not isinstance(v, bool) else f"{short}={v}")
    rid = "__".join(parts).replace("/", "-")
    if len(rid) > 180:
        # keep it readable but unique: varying keys stay visible via the hash of the full id
        import hashlib
        rid = rid[:150] + "__" + hashlib.md5(rid.encode()).hexdigest()[:10]
    return f"seed={ov.get('seed', 1)}__" + rid


def run_one(args):
    ov, out_root, burn_in = args
    rid = run_id(ov)
    out = os.path.join(out_root, "runs", rid)
    if os.path.exists(os.path.join(out, "metrics.json")):
        with open(os.path.join(out, "metrics.json")) as f:
            return json.load(f)
    os.makedirs(out, exist_ok=True)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from growing.sim import run_sim as rs
        from growing.analysis import avalanches, probes

        from growing.sim.axon_nets import AxonConfig, render_axo_frame
        morph, branch, lif, analysis = rs.default_configs()
        axon = AxonConfig()
        cfgs = dict(morph=morph, branch=branch, lif=lif, analysis=analysis, axon=axon)
        seed = int(ov.get("seed", 1))
        morph.morphology_seed = seed
        lif.animate = False
        for k, v in ov.items():
            if isinstance(v, str) and v.startswith("FROM_"):
                raise ValueError(f"placeholder {k}={v} not filled in")
            if k in ("seed", "frozen_steps", "fractal", "fractal_dynamic") or k.startswith("_"):
                continue
            c, field = k.split(".", 1)
            if not hasattr(cfgs[c], field):
                raise KeyError(f"unknown config field {k}")
            setattr(cfgs[c], field, v)
        with open(os.path.join(out, "config.json"), "w") as f:
            json.dump(_jsonable({c: vars(o) for c, o in cfgs.items()} | {"overrides": ov}), f, indent=1)

        random.seed(seed); np.random.seed(seed)
        t0 = time.time()
        sim = rs.BranchMorphPlasticLIFSim(morph, lif, branch, axon)
        res = sim.simulate(analysis_cfg=None)
        wall = time.time() - t0

        S = np.asarray(res["spikes_list"], bool)
        R = np.asarray(res["firing_rates_list"], float)
        L = np.asarray(res["length_list"], float)
        M = res["connectivity_matrix"]
        steps = len(S)
        b0 = int(steps * burn_in)

        av = avalanches.analyze(S, bin_steps=1, threshold_percentile=50.0, burn_in_fraction=burn_in)
        avs = avalanches.summary(av)
        # Avalanches in coarser bins (robustness of exponents to the bin choice)
        av4 = avalanches.analyze(S, bin_steps=4, threshold_percentile=50.0, burn_in_fraction=burn_in)
        avs4 = {f"bin4_{k}": v for k, v in avalanches.summary(av4).items()
                if k in ("n_avalanches", "size_exponent_tau", "duration_exponent_alpha", "gamma_fit",
                         "gamma_predicted", "branching_ratio_m", "max_size", "threshold")}
        wb = probes.windowed_branching(S, n_windows=10)
        # Quenched statistics: freeze the grown network and record a long LIF-only run.
        frozen_steps = int(ov.get("frozen_steps", FROZEN_STEPS_DEFAULT))
        frz = {}
        if frozen_steps > 0:
            t1 = time.time()
            Sf = probes.lif_run(M, lif, frozen_steps, rng=np.random.default_rng(seed + 1000))
            avf = avalanches.analyze(Sf, bin_steps=1, threshold_percentile=50.0, burn_in_fraction=0.02)
            frz = {f"frozen_{k}": v for k, v in avalanches.summary(avf).items()}
            frz["frozen_spike_prob_per_step"] = float(Sf.mean())
            frz["frozen_rate"] = float(Sf.mean() / lif.dt)
            frz["frozen_wall_s"] = time.time() - t1
            fig = avalanches.plot(avf, dt=lif.dt, title="frozen network: " + rid[:100])
            fig.savefig(os.path.join(out, "avalanches_frozen.png"), dpi=80); plt.close(fig)
            np.savez_compressed(os.path.join(out, "frozen_avalanches.npz"), sizes=avf["sizes"],
                                durations=avf["durations"], threshold=avf["threshold"])
        dr = probes.dynamic_range(M, lif, seed=seed)
        cs = probes.coupling_spectrum(M, lif.spike_input)
        gs = probes.graph_stats(M)
        total_px = int(np.count_nonzero(sim.branch_net.filled_mask))

        late_rates = R[b0:].mean(0)
        ct = np.asarray(res.get("coupling_trace", np.zeros((0, 2))))
        late_ct = ct[ct[:, 0] >= b0, 1] if len(ct) else np.array([np.nan])
        sigma_stats = dict(sigma_late_mean=float(np.mean(late_ct)), sigma_late_std=float(np.std(late_ct)),
                           sigma_final=float(late_ct[-1]) if len(late_ct) else np.nan)
        metrics = dict(
            run_id=rid, seed=seed, wall_s=wall, steps=steps, N=sim.N,
            mean_rate=float(late_rates.mean()), rate_cv_across_neurons=float(late_rates.std() / max(late_rates.mean(), 1e-12)),
            frac_neurons_within_10pct=float(np.mean(np.abs(late_rates - lif.tfr) < 0.1 * lif.tfr)),
            spike_prob_per_step=float(S[b0:].mean()),
            mean_length=float(L[-1].mean()), max_length=float(L[-1].max()),
            length_drift_last_quarter=float(L[-1].mean() - L[int(0.75 * steps)].mean()),
            arbor_pixels=total_px, arbor_coverage=total_px / float(morph.width * morph.height),
            **avs, **avs4, **frz,
            m_windows=wb["m"].tolist(), m_first_window=float(wb["m"][0]) if len(wb["m"]) else np.nan,
            dynamic_range_db=dr["delta_db"], dr_F0=dr["F0"], dr_Fmax=dr["Fmax"],
            **cs, **gs, **sigma_stats, **res.get("axon_stats", {}),
        )
        if ov.get("fractal", True):
            # Morphology geometry: box-counting D (FracStack, min counts, unweighted OLS) per neuron
            # and for the whole network, lacunarity, tortuosity, connection-distance decay,
            # modularity and the largest eigenvalue of the transmission matrix.
            try:
                from growing.analysis import fractal
                t2 = time.time()
                fm, per = fractal.morphology_metrics(sim, M, lif, dynamic=bool(ov.get("fractal_dynamic", True)), seed=seed)
                fm["fractal_wall_s"] = time.time() - t2
                metrics.update(fm)
                np.savez_compressed(os.path.join(out, "morphology_metrics.npz"), **per)
            except Exception:
                metrics["fractal_error"] = traceback.format_exc()[-500:]
        metrics.update({f"param.{k}": v for k, v in ov.items()})

        np.savez_compressed(os.path.join(out, "spikes.npz"), spikes=S, dt=lif.dt)
        np.save(os.path.join(out, "M.npy"), M)
        np.savez_compressed(os.path.join(out, "traces.npz"), rates=R[::10].astype(np.float32),
                            lengths=L[::10].astype(np.float32), dr_stim=dr["stim_rates"], dr_F=dr["F"],
                            m_windows=wb["m"], window_edges=wb["window_edges"],
                            coupling_trace=res.get("coupling_trace", np.zeros((0, 2))))
        try:
            fig = avalanches.plot(av, dt=lif.dt, title=rid[:120])
            fig.savefig(os.path.join(out, "avalanches.png"), dpi=80); plt.close(fig)
            if sim.axons is not None:
                frame = render_axo_frame(sim.branch_net, sim.axons)
            else:
                frame = sim.branch_net.render_frame(mask_type="filled", colormap="summer", show_connections=False,
                                                    show_fusions=False, show_soma_connections=False,
                                                    background_color="black", figsize=(6, 6), dpi=100)
            import imageio.v2 as imageio
            imageio.imwrite(os.path.join(out, "morphology.png"), frame)
        except Exception as e:  # plotting must never lose a run
            print("plot failed:", e)
        with open(os.path.join(out, "metrics.json"), "w") as f:
            json.dump(_jsonable(metrics), f, indent=1)
        return _jsonable(metrics)
    except Exception:
        err = traceback.format_exc()
        with open(os.path.join(out, "error.txt"), "w") as f:
            f.write(err)
        return dict(run_id=rid, error=err.splitlines()[-1], **{f"param.{k}": v for k, v in ov.items()})


def write_csv(results, path):
    keys = []
    for r in results:
        for k in r:
            if k not in keys and not isinstance(r[k], (list, dict)):
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k) for k in keys})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("spec")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--burn-in", type=float, default=0.4)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    with open(a.spec) as f:
        spec = json.load(f)
    out_root = os.path.join(SWEEPS_DIR, spec["name"])
    os.makedirs(out_root, exist_ok=True)
    with open(os.path.join(out_root, "spec.json"), "w") as f:
        json.dump(spec, f, indent=1)
    runs = expand(spec)
    print(f"{spec['name']}: {len(runs)} runs, {a.workers} workers", flush=True)
    if a.dry_run:
        for r in runs:
            print(run_id(r))
        return
    results = []
    ctx = get_context("spawn")
    with ctx.Pool(a.workers, maxtasksperchild=1) as pool:
        for i, r in enumerate(pool.imap_unordered(run_one, [(ov, out_root, a.burn_in) for ov in runs])):
            results.append(r)
            fmt = lambda k: (f"{r.get(k):.3g}" if isinstance(r.get(k), (int, float)) else str(r.get(k)))
            status = r.get("error") or (f"rate={fmt('mean_rate')} m={fmt('branching_ratio_m')} | frozen: tau={fmt('frozen_size_exponent_tau')} "
                                        f"alpha={fmt('frozen_duration_exponent_alpha')} m={fmt('frozen_branching_ratio_m')} "
                                        f"nav={fmt('frozen_n_avalanches')} | DR={fmt('dynamic_range_db')} wall={fmt('wall_s')}")
            print(f"[{i+1}/{len(runs)}] {r['run_id'][:90]} :: {status}", flush=True)
            write_csv(results, os.path.join(out_root, "results.csv"))
    print("done", flush=True)


if __name__ == "__main__":
    main()
