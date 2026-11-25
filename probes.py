"""
Post-hoc probes of a grown network (morphology / connectivity frozen).

* dynamic_range: drive the frozen LIF network with an extra Poisson stimulus of rate s
  (kicks of size ext_kick per neuron) over a log-spaced range of s, measure the mean firing
  rate F(s), and compute Delta = 10 log10(s_90 / s_10) (Kinouchi & Copelli, Nat Phys 2006),
  where s_x is the stimulus giving F = F_0 + x% (F_max - F_0). Criticality is expected to
  maximise Delta.
* coupling_spectrum: largest eigenvalue of the effective coupling W = spike_input * M.
* graph_stats: degree, clustering, characteristic path length, number of components.
* windowed_branching: MR branching ratio in consecutive windows of the run, to follow how
  the network's dynamical state develops while it grows.
"""
from __future__ import annotations

import numpy as np

import avalanches


def lif_run(M, lif, steps: int, stim_rate: float = 0.0, stim_kick: float | None = None,
            rng: np.random.Generator | None = None, V0=None):
    """Simulate the LIF network of run_sim with fixed coupling M (same update rule).
    Returns the (steps, N) spike raster."""
    rng = rng or np.random.default_rng()
    N = M.shape[0]
    dt = lif.dt
    V = rng.uniform(-1, 1, N) if V0 is None else np.array(V0, float)
    spikes = np.zeros(N, bool)
    out = np.zeros((steps, N), bool)
    kick = lif.ext_kick if stim_kick is None else stim_kick
    for t in range(steps):
        rec = lif.spike_input * (M @ spikes.astype(float))
        if lif.external_drive == "poisson":
            drive = lif.ext_kick * rng.poisson(lif.ext_rate * dt, N)
            dV = -lif.leak * (V - lif.rest)
        else:
            drive = 0.0
            dV = lif.input_noise * rng.random(N) - lif.leak * (V - lif.rest)
        stim = kick * rng.poisson(stim_rate * dt, N) if stim_rate > 0 else 0.0
        V = V + rec + drive + stim + dV * dt
        spikes = V >= lif.spike_threshold
        V[spikes] = lif.rest
        V[V < lif.rest] = lif.rest
        out[t] = spikes
    return out


def dynamic_range(M, lif, stim_rates=None, steps: int = 1500, warmup: int = 200, seed: int = 0):
    rng = np.random.default_rng(seed)
    if stim_rates is None:
        stim_rates = np.logspace(-4, 1, 21)
    F = []
    for s in stim_rates:
        S = lif_run(M, lif, steps + warmup, stim_rate=float(s), rng=rng)
        F.append(S[warmup:].mean() / lif.dt)          # spikes per neuron per time unit
    F = np.array(F)
    F0, Fmax = F[0], F[-1]
    out = dict(stim_rates=np.asarray(stim_rates), F=F, F0=float(F0), Fmax=float(Fmax), delta_db=np.nan)
    if Fmax - F0 <= 1e-9:
        return out
    logs = np.log10(stim_rates)

    def s_at(frac):
        target = F0 + frac * (Fmax - F0)
        Fm = np.maximum.accumulate(F)                 # monotone envelope for interpolation
        idx = np.searchsorted(Fm, target)
        if idx <= 0:
            return stim_rates[0]
        if idx >= len(Fm):
            return stim_rates[-1]
        x0, x1, y0, y1 = logs[idx - 1], logs[idx], Fm[idx - 1], Fm[idx]
        return 10 ** (x0 + (target - y0) * (x1 - x0) / max(y1 - y0, 1e-12))

    s10, s90 = s_at(0.1), s_at(0.9)
    out.update(s10=float(s10), s90=float(s90), delta_db=float(10 * np.log10(s90 / s10)))
    return out


def coupling_spectrum(M, spike_input: float):
    W = spike_input * np.asarray(M, float)
    ev = np.linalg.eigvals(W)
    return dict(lambda_max=float(np.max(np.abs(ev))), row_sum_mean=float(W.sum(1).mean()),
                row_sum_max=float(W.sum(1).max()))


def graph_stats(M, eps: float = 1e-12):
    A = (np.asarray(M) > eps).astype(int)
    np.fill_diagonal(A, 0)
    N = len(A)
    deg = A.sum(1)
    # clustering coefficient (undirected, unweighted)
    tri = np.diag(A @ A @ A) / 2.0
    poss = deg * (deg - 1) / 2.0
    clust = np.where(poss > 0, tri / np.maximum(poss, 1), 0.0)
    # BFS shortest paths
    dist = np.full((N, N), np.inf)
    for s in range(N):
        dist[s, s] = 0
        frontier = [s]; d = 0
        while frontier:
            d += 1
            nxt = []
            for u in frontier:
                for v in np.flatnonzero(A[u]):
                    if dist[s, v] == np.inf:
                        dist[s, v] = d; nxt.append(v)
            frontier = nxt
    finite = np.isfinite(dist) & ~np.eye(N, dtype=bool)
    # connected components
    seen = np.zeros(N, bool); comps = []
    for s in range(N):
        if not seen[s]:
            members = np.flatnonzero(np.isfinite(dist[s])); seen[members] = True; comps.append(len(members))
    p = deg.mean() / max(N - 1, 1)
    C_rand = p
    L_rand = np.log(N) / np.log(max(deg.mean(), 1.0001)) if deg.mean() > 1 else np.nan
    C = float(clust.mean())
    L = float(dist[finite].mean()) if finite.any() else np.nan
    sw = (C / C_rand) / (L / L_rand) if (C_rand > 0 and np.isfinite(L) and np.isfinite(L_rand) and L > 0) else np.nan
    return dict(mean_degree=float(deg.mean()), max_degree=int(deg.max()), isolated=int(np.sum(deg == 0)),
                clustering=C, path_length=L, n_components=len(comps), largest_component=int(max(comps)),
                small_world_sigma=float(sw) if np.isfinite(sw) else np.nan)


def windowed_branching(spikes, n_windows: int = 10, bin_steps: int = 1):
    S = np.asarray(spikes, bool)
    edges = np.linspace(0, len(S), n_windows + 1).astype(int)
    ms, rates = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        A = avalanches.population_activity(S[a:b], bin_steps)
        ms.append(avalanches.mr_branching_ratio(A)['m'] if len(A) > 50 else np.nan)
        rates.append(float(S[a:b].mean()))
    return dict(window_edges=edges, m=np.array(ms), spike_prob=np.array(rates))
