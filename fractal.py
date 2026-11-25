"""
Fractal / morphology / graph-geometry measurements of a grown network (plan: plans/FRACTAL_SWEEPS.md, step 0).

Box counting uses FracStack (../FracStack, github.com/apdill/FracStack) with
    use_min_count=True     (minimum box count over grid offsets)
    unweighted fits        (ordinary least squares on log N vs log size)
Two estimates are reported for every mask:
    D_fixed : OLS slope over a fixed, common size range (comparable across conditions)
    D_dyn   : FracStack dynamic_boxcount (best-R² range spanning >= `min_decade_span` decades)

FracStack's package __init__ also imports its 3-D module, which needs torch; this environment's
`torch` is an empty namespace package, so the 2-D modules are loaded directly (see _fracstack()).
Set FRACSTACK_PATH to override the location of the FracStack repository.

Main entry point: morphology_metrics(sim, M, lif) -> flat dict of scalars (+ per-neuron arrays).
"""
from __future__ import annotations

import contextlib
import importlib
import io
import os
import sys
import types

import numpy as np

_FS = None
FRACSTACK_PATH = os.environ.get(
    "FRACSTACK_PATH", os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "FracStack"))

NEURON_FIXED_RANGE = (4, 96)      # px: dendritic-arbor scales (arbor spans are ~300-800 px)
NETWORK_MIN_SIZE = 4               # px; network fixed range is [4, min(H, W) / 8]
MAX_SIZE_FRACTION = 8              # largest box <= extent / 8 limits finite-object edge bias
                                   # (validation: carpet 1.844 vs 1.893, square 1.957, Koch 1.289 vs 1.262)
LACUNARITY_SIZES = (4, 8, 16, 32, 64)


# --------------------------------------------------------------------------- FracStack access
def _fracstack():
    """Return FracStack's boxcount module without running the package __init__ (torch-free)."""
    global _FS
    if _FS is not None:
        return _FS
    pkg_dir = os.path.join(FRACSTACK_PATH, "fracstack")
    if not os.path.isdir(pkg_dir):
        raise ImportError(f"FracStack not found at {FRACSTACK_PATH} (set FRACSTACK_PATH)")
    if "fracstack" not in sys.modules or not hasattr(sys.modules["fracstack"], "__path__"):
        pkg = types.ModuleType("fracstack")
        pkg.__path__ = [pkg_dir]
        sys.modules["fracstack"] = pkg
    _FS = importlib.import_module("fracstack.boxcount")
    return _FS


def _crop(mask, pad=2):
    ys, xs = np.nonzero(mask)
    if len(ys) == 0:
        return None
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, mask.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, mask.shape[1])
    return np.ascontiguousarray(mask[y0:y1, x0:x1].astype(np.uint8))


def box_dimension(mask, fixed_range=None, dynamic=True, num_sizes=40, num_offsets=20,
                  min_decade_span=1.0, seed=0) -> dict:
    """Box-counting dimension of a binary mask (min counts, unweighted OLS fits).

    fixed_range: (min_size, max_size) in px for D_fixed; default (4, min(shape)//4).
    dynamic: also run FracStack dynamic_boxcount (best-R² window of >= min_decade_span decades)."""
    fs = _fracstack()
    out = dict(D_fixed=np.nan, R2_fixed=np.nan, fixed_lo=np.nan, fixed_hi=np.nan,
               D_dyn=np.nan, R2_dyn=np.nan, dyn_lo=np.nan, dyn_hi=np.nan, dyn_decades=np.nan,
               n_pixels=int(np.count_nonzero(mask)))
    m = _crop(np.asarray(mask) > 0)
    if m is None or out["n_pixels"] < 20:
        return out
    lo, hi = fixed_range if fixed_range is not None else (NETWORK_MIN_SIZE, max(m.shape) // MAX_SIZE_FRACTION)
    hi = int(min(hi, max(max(m.shape) // MAX_SIZE_FRACTION, lo * 4)))
    if hi >= lo * 4:
        padded = fs.pad_image_for_boxcounting(m, hi, pad_factor=1.5)
        with contextlib.redirect_stdout(io.StringIO()):
            sizes, counts = fs.boxcount(padded, mode="D0", min_size=int(lo), max_size=int(hi),
                                        num_sizes=num_sizes, num_offsets=num_offsets,
                                        use_min_count=True, seed=seed)
        try:
            _, _, d, _, r2, _, _ = fs.compute_dimension(sizes, counts, mode="D0", use_weighted_fit=False,
                                                        use_bootstrap_ci=False)
            out.update(D_fixed=float(d), R2_fixed=float(r2), fixed_lo=int(lo), fixed_hi=int(hi))
        except Exception:
            pass
    if dynamic:
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                best = fs.dynamic_boxcount(m, mode="D0", global_min_size=2,
                                           global_max_size=max(max(m.shape) // MAX_SIZE_FRACTION, 16),
                                           min_decade_span=min_decade_span, min_points=8, num_sizes=num_sizes,
                                           num_offsets=num_offsets, show_progress=False, use_min_count=True,
                                           seed=seed, use_bootstrap_ci=False)
            out.update(D_dyn=float(best["D_value"]), R2_dyn=float(best["R2"]),
                       dyn_lo=float(best["optimal_range"][0]), dyn_hi=float(best["optimal_range"][1]),
                       dyn_decades=float(best["decade_span"]))
        except Exception as e:
            out["dyn_error"] = repr(e)[:200]
    return out


# --------------------------------------------------------------------------- lacunarity
def lacunarity(mask, sizes=LACUNARITY_SIZES) -> dict:
    """Gliding-box lacunarity Λ(r) = E[M²]/E[M]² over all r×r boxes (Allain & Cloitre 1991)."""
    m = (np.asarray(mask) > 0).astype(np.float64)
    c = _crop(m > 0, pad=max(sizes))
    if c is None:
        return {f"lac_{r}": np.nan for r in sizes}
    I = np.pad(c.astype(np.float64), ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    out = {}
    for r in sizes:
        if r >= min(c.shape):
            out[f"lac_{r}"] = np.nan
            continue
        M = I[r:, r:] - I[:-r, r:] - I[r:, :-r] + I[:-r, :-r]
        mu = M.mean()
        out[f"lac_{r}"] = float((M ** 2).mean() / mu ** 2) if mu > 0 else np.nan
    return out


# --------------------------------------------------------------------------- branch geometry
def branch_geometry(neuron) -> dict:
    """Tortuosity (path length / chord) per branch (all segments sharing a tip_id, ordered by path
    distance), branch count and total cable length of one neuron."""
    by_tip = {}
    total = 0.0
    for seg in getattr(neuron, "branch_data", []):
        pts = np.asarray(seg.get("points", []), float)
        if len(pts) < 2:
            continue
        L = float(np.linalg.norm(np.diff(pts, axis=0), axis=1).sum())
        total += L
        by_tip.setdefault(seg.get("tip_id", -1), []).append((seg.get("path_start_dist", 0.0), pts, L))
    torts = []
    for segs in by_tip.values():
        segs.sort(key=lambda s: s[0])
        L = sum(s[2] for s in segs)
        chord = float(np.linalg.norm(segs[-1][1][-1] - segs[0][1][0]))
        if L > 10 and chord > 1e-6:
            torts.append(L / chord)
    return dict(tortuosity=float(np.mean(torts)) if torts else np.nan, n_branches=len(by_tip),
                cable_length=total)


# --------------------------------------------------------------------------- graph geometry
def soma_positions(net):
    P = []
    for n in net.neurons:
        p = getattr(getattr(n, "soma", None), "position", None)
        P.append(np.asarray(p, float) if p is not None else np.array([np.nan, np.nan]))
    return np.array(P)


def connection_distance(M, pos, n_bins=12) -> dict:
    """P(connection) vs soma distance; exponential decay length λ_c from a log-linear fit."""
    A = (np.asarray(M) > 0)
    A = A | A.T
    N = len(A)
    iu = np.triu_indices(N, 1)
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1)[iu]
    a = A[iu]
    out = dict(conn_decay_len=np.nan, conn_p0=np.nan, mean_conn_dist=float(d[a].mean()) if a.any() else np.nan)
    if a.sum() < 5:
        return out
    edges = np.linspace(0, np.percentile(d, 95), n_bins + 1)
    idx = np.digitize(d, edges) - 1
    ctr, p = [], []
    for b in range(n_bins):
        sel = idx == b
        if sel.sum() >= 5:
            ctr.append(0.5 * (edges[b] + edges[b + 1])); p.append(a[sel].mean())
    ctr, p = np.array(ctr), np.array(p)
    ok = p > 0
    if ok.sum() >= 3:
        slope, icpt = np.polyfit(ctr[ok], np.log(p[ok]), 1)
        if slope < 0:
            out.update(conn_decay_len=float(-1 / slope), conn_p0=float(np.exp(icpt)))
    out["conn_prob_bins"] = ctr.tolist(); out["conn_prob"] = p.tolist()
    return out


def modularity(M, seed=0) -> dict:
    import networkx as nx
    W = np.asarray(M, float)
    W = W + W.T
    np.fill_diagonal(W, 0)
    G = nx.from_numpy_array(W)
    if G.number_of_edges() == 0:
        return dict(modularity=np.nan, n_modules=0)
    comms = nx.algorithms.community.louvain_communities(G, weight="weight", seed=seed)
    return dict(modularity=float(nx.algorithms.community.modularity(G, comms, weight="weight")),
                n_modules=len(comms))


def lambda_max(M, lif) -> float:
    """Largest |eigenvalue| of the effective transmission matrix (P = min(1, spike_input·M) in
    branching mode, spike_input·M otherwise) — the quantity that sets criticality (Larremore 2011)."""
    W = lif.spike_input * np.asarray(M, float)
    if getattr(lif, "spike_mode", "threshold") == "branching":
        W = np.clip(W, 0, 1)
    return float(np.max(np.abs(np.linalg.eigvals(W))))


# --------------------------------------------------------------------------- aggregate
def morphology_metrics(sim, M, lif, dynamic=True, seed=0) -> tuple[dict, dict]:
    """Scalar metrics for sweep results + per-neuron arrays (saved separately)."""
    net = sim.branch_net
    H, W = net.filled_mask.shape
    somas = net.somas_mask
    per = dict(D_fixed=[], D_dyn=[], R2_fixed=[], D_skel=[], D_skel_dyn=[], tortuosity=[], n_branches=[],
               cable_length=[], lac_16=[])
    for i, n in enumerate(net.neurons):
        # image-like mask (filled dendrites with thickness + own soma), comparable to published
        # box-counting of neuron images; and the arbor skeleton alone (branching pattern only)
        image = (net._branch_filled_counts[i] > 0) | (somas == n.id_int)
        skel = net._branch_skeleton_counts[i] > 0
        bd = box_dimension(image, fixed_range=NEURON_FIXED_RANGE, dynamic=dynamic, seed=seed)
        bs = box_dimension(skel, fixed_range=NEURON_FIXED_RANGE, dynamic=dynamic, seed=seed)
        bg = branch_geometry(n)
        per["D_fixed"].append(bd["D_fixed"]); per["D_dyn"].append(bd["D_dyn"]); per["R2_fixed"].append(bd["R2_fixed"])
        per["D_skel"].append(bs["D_fixed"]); per["D_skel_dyn"].append(bs["D_dyn"])
        per["tortuosity"].append(bg["tortuosity"]); per["n_branches"].append(bg["n_branches"])
        per["cable_length"].append(bg["cable_length"])
        per["lac_16"].append(lacunarity(net._branch_filled_counts[i] > 0, (16,))["lac_16"])
    per = {k: np.asarray(v, float) for k, v in per.items()}

    net_skel = (net.skeleton_mask > 0) | (somas > 0)
    net_fill = (net.filled_mask > 0) | (somas > 0)
    hi = min(H, W) // MAX_SIZE_FRACTION
    nd_s = box_dimension(net_skel, fixed_range=(NETWORK_MIN_SIZE, hi), dynamic=dynamic, seed=seed)
    nd_f = box_dimension(net_fill, fixed_range=(NETWORK_MIN_SIZE, hi), dynamic=dynamic, seed=seed)
    lac = lacunarity(net_fill)
    pos = soma_positions(net)
    cd = connection_distance(M, pos)
    out = dict(
        D_neuron_mean=float(np.nanmean(per["D_fixed"])), D_neuron_sd=float(np.nanstd(per["D_fixed"])),
        D_neuron_dyn_mean=float(np.nanmean(per["D_dyn"])), D_neuron_R2_mean=float(np.nanmean(per["R2_fixed"])),
        D_neuron_skel_mean=float(np.nanmean(per["D_skel"])), D_neuron_skel_dyn_mean=float(np.nanmean(per["D_skel_dyn"])),
        D_network_skel=nd_s["D_fixed"], D_network_skel_R2=nd_s["R2_fixed"], D_network_skel_dyn=nd_s["D_dyn"],
        D_network_fill=nd_f["D_fixed"], D_network_fill_R2=nd_f["R2_fixed"], D_network_fill_dyn=nd_f["D_dyn"],
        D_network_dyn_lo=nd_f["dyn_lo"], D_network_dyn_hi=nd_f["dyn_hi"],
        tortuosity_mean=float(np.nanmean(per["tortuosity"])),
        branches_per_neuron=float(np.nanmean(per["n_branches"])),
        cable_length_mean=float(np.nanmean(per["cable_length"])),
        lacunarity_neuron_16=float(np.nanmean(per["lac_16"])),
        **{f"lacunarity_net_{k[4:]}": v for k, v in lac.items()},
        conn_decay_len=cd["conn_decay_len"], conn_p0=cd["conn_p0"], mean_conn_dist=cd["mean_conn_dist"],
        lambda_max_P=lambda_max(M, lif),
        **modularity(M, seed=seed),
    )
    per["soma_pos"] = pos
    per["conn_prob_bins"] = np.asarray(cd.get("conn_prob_bins", []))
    per["conn_prob"] = np.asarray(cd.get("conn_prob", []))
    return out, per


# --------------------------------------------------------------------------- validation
def _koch(level=5, size=1024):
    import cv2
    pts = np.array([[0.0, 0.0], [1.0, 0.0]])
    for _ in range(level):
        new = []
        for a, b in zip(pts[:-1], pts[1:]):
            d = (b - a) / 3
            p1, p2 = a + d, a + 2 * d
            rot = np.array([[0.5, -np.sqrt(3) / 2], [np.sqrt(3) / 2, 0.5]])
            new += [a, p1, p1 + rot @ d, p2]
        new.append(pts[-1]); pts = np.array(new)
    img = np.zeros((size // 3 + 40, size + 40), np.uint8)
    q = (pts * (size - 1) + 20).astype(np.int32)
    q[:, 1] = img.shape[0] - 20 - q[:, 1]
    cv2.polylines(img, [q.reshape(-1, 1, 2)], False, 1, 1)
    return img


def _sierpinski_carpet(level=5):
    c = np.ones((1, 1), np.uint8)
    for _ in range(level):
        z = np.zeros_like(c)
        c = np.block([[c, c, c], [c, z, c], [c, c, c]])
    return c


def validate():
    """Known-D test images: line (1), Koch curve (1.262), Sierpinski carpet (1.893), filled square (2)."""
    tests = {"line (1.000)": np.pad(np.ones((1, 900), np.uint8), 50),
             "koch (1.262)": _koch(),
             "carpet (1.893)": _sierpinski_carpet(6),
             "square (2.000)": np.pad(np.ones((800, 800), np.uint8), 50)}
    rows = {}
    for name, img in tests.items():
        r = box_dimension(img, fixed_range=None)
        rows[name] = r
        print(f"{name:16s} D_fixed={r['D_fixed']:.3f} (R²={r['R2_fixed']:.4f}, {r['fixed_lo']}-{r['fixed_hi']} px)  "
              f"D_dyn={r['D_dyn']:.3f} (R²={r['R2_dyn']:.4f}, {r['dyn_lo']}-{r['dyn_hi']} px)")
    return rows


if __name__ == "__main__":
    validate()
