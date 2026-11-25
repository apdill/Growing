"""
Neuronal-avalanche analysis for continuously active spiking networks.

The LIF network here is noise-driven and never falls silent (≈10 spikes per step), so the
classic definition (avalanche = run of non-empty bins bounded by empty bins) yields a single
never-ending avalanche. We therefore use *thresholded* avalanches (Poil et al., J Neurosci
2012; Del Papa et al., PLoS ONE 2017):

    A(t)      population spike count in bins of `bin_steps` simulation steps
    theta     integer threshold = percentile of A over the analysed window
    avalanche maximal run of consecutive bins with A(t) > theta
    size S    sum over the run of (A(t) - theta)     (integer >= 1)
    duration T number of bins in the run             (integer >= 1)

Criticality diagnostics reported:
  * discrete power-law MLE for P(S) ~ S^-tau and P(T) ~ T^-alpha, with xmin chosen by
    minimising the KS distance (Clauset, Shalizi & Newman, SIAM Rev 2009)
  * <S>(T) ~ T^gamma, compared with the crackling-noise prediction
    gamma_pred = (alpha - 1) / (tau - 1)  (Sethna et al. 2001; Friedman et al. 2012)
  * branching ratio m from the multistep-regression (MR) estimator
    (Wilting & Priesemann, Nat Commun 2018); m ~ 1 is critical, m < 1 subcritical.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import zeta


# --------------------------------------------------------------------------- detection
def population_activity(spikes, bin_steps: int = 1) -> np.ndarray:
    """Spike count per bin from a (T, N) boolean spike raster."""
    S = np.asarray(spikes, dtype=bool)
    counts = S.sum(axis=1).astype(np.int64)
    bin_steps = max(1, int(bin_steps))
    n = (len(counts) // bin_steps) * bin_steps
    return counts[:n].reshape(-1, bin_steps).sum(axis=1)


def detect_avalanches(activity: np.ndarray, threshold: int):
    """Sizes (excess above threshold) and durations of supra-threshold excursions."""
    A = np.asarray(activity, dtype=np.int64)
    above = A > threshold
    if not above.any():
        return np.array([], int), np.array([], int)
    edges = np.diff(np.concatenate([[0], above.astype(np.int8), [0]]))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    # Drop runs touching the window boundaries: their true extent is unknown.
    keep = (starts > 0) & (ends < len(A))
    starts, ends = starts[keep], ends[keep]
    excess = np.concatenate([[0], np.cumsum(A - threshold)])
    sizes = excess[ends] - excess[starts]
    durations = ends - starts
    return sizes.astype(int), durations.astype(int)


# --------------------------------------------------------------------------- power-law fit
def _discrete_mle_alpha(x: np.ndarray, xmin: int) -> float:
    x = x[x >= xmin]
    n = len(x)
    slx = np.log(x).sum()
    nll = lambda a: n * np.log(zeta(a, xmin)) + a * slx
    res = minimize_scalar(nll, bounds=(1.01, 8.0), method="bounded")
    return float(res.x)


def _ks_distance(x: np.ndarray, xmin: int, alpha: float) -> float:
    x = np.sort(x[x >= xmin])
    vals, counts = np.unique(x, return_counts=True)
    emp_cdf = np.cumsum(counts) / len(x)
    model_ccdf = zeta(alpha, vals + 1) / zeta(alpha, xmin)      # P(X > v)
    return float(np.max(np.abs(emp_cdf - (1.0 - model_ccdf))))


def fit_discrete_power_law(x, min_tail: int = 50, max_xmin_candidates: int = 40) -> dict:
    """Clauset-style discrete power-law fit; returns alpha, xmin, KS distance, tail size."""
    x = np.asarray(x, dtype=np.int64)
    x = x[x >= 1]
    out = dict(alpha=np.nan, xmin=np.nan, ks=np.nan, n_tail=0, n=int(len(x)))
    if len(x) < min_tail:
        return out
    cands = np.unique(x)
    cands = cands[[np.count_nonzero(x >= c) >= min_tail for c in cands]]
    if len(cands) == 0:
        return out
    if len(cands) > max_xmin_candidates:
        cands = np.unique(np.round(np.geomspace(cands[0], cands[-1], max_xmin_candidates)).astype(int))
        cands = cands[[np.count_nonzero(x >= c) >= min_tail for c in cands]]
    best = None
    for xmin in cands:
        a = _discrete_mle_alpha(x, int(xmin))
        d = _ks_distance(x, int(xmin), a)
        if best is None or d < best[2]:
            best = (int(xmin), a, d)
    xmin, a, d = best
    tail = x[x >= xmin]
    # Standard error of the MLE (Clauset eq. 3.6, continuous approximation)
    se = (a - 1) / np.sqrt(len(tail))
    out.update(alpha=a, alpha_se=float(se), xmin=xmin, ks=d, n_tail=int(len(tail)))
    return out


# --------------------------------------------------------------------------- other diagnostics
def size_vs_duration_exponent(sizes, durations, min_count: int = 5) -> dict:
    T = np.asarray(durations); S = np.asarray(sizes, float)
    uT = np.unique(T)
    meanS = np.array([S[T == t].mean() for t in uT])
    cnt = np.array([np.count_nonzero(T == t) for t in uT])
    ok = (cnt >= min_count) & (uT >= 1) & (meanS > 0)
    out = dict(gamma=np.nan, T=uT, meanS=meanS, used=ok)
    if ok.sum() >= 3:
        gamma, logc = np.polyfit(np.log(uT[ok]), np.log(meanS[ok]), 1)
        out.update(gamma=float(gamma), log_c=float(logc))
    return out


def mr_branching_ratio(activity, k_max: int = 40) -> dict:
    """Multistep-regression estimate of the branching ratio m (Wilting & Priesemann 2018)."""
    A = np.asarray(activity, float)
    ks = np.arange(1, min(k_max, len(A) // 4) + 1)
    rk = []
    for k in ks:
        a, b = A[:-k], A[k:]
        va = a.var()
        rk.append(np.cov(a, b, bias=True)[0, 1] / va if va > 0 else np.nan)
    rk = np.array(rk)
    ok = np.isfinite(rk)
    out = dict(m=np.nan, k=ks, r_k=rk)
    if ok.sum() >= 3:
        # Least-squares fit of r_k = b * m^k on the raw coefficients (as in the MR estimator);
        # a log-linear fit is biased towards m -> 1 by the noisy near-zero tail.
        from scipy.optimize import curve_fit
        # 'exp_offset' variant (Wilting & Priesemann 2018): r_k = b*m^k + c, where the
        # offset c absorbs slow non-stationarity (here: ongoing structural plasticity).
        r1 = float(np.clip(rk[ok][0], 1e-3, 0.999))
        try:
            (b, m), _ = curve_fit(lambda k, b, m: b * m ** k, ks[ok], rk[ok], p0=(1.0, r1),
                                  bounds=([0.0, 0.0], [np.inf, 1.5]), maxfev=10000)
            out.update(m_exp=float(m))
        except Exception:
            pass
        try:
            (b, m, c), _ = curve_fit(lambda k, b, m, c: b * m ** k + c, ks[ok], rk[ok],
                                     p0=(r1, 0.5, 0.0), bounds=([0.0, 0.0, -1.0], [np.inf, 1.5, 1.0]),
                                     maxfev=20000)
            out.update(m=float(m), b=float(b), c=float(c))
        except Exception:
            out['m'] = out.get('m_exp', np.nan)
    out['r1'] = float(rk[0]) if len(rk) else np.nan
    return out


def analyze(spikes, bin_steps: int = 1, threshold_percentile: float = 50.0,
            burn_in_fraction: float = 0.3) -> dict:
    """Full avalanche analysis on the post-burn-in part of a spike raster."""
    S = np.asarray(spikes, dtype=bool)
    start = int(len(S) * float(np.clip(burn_in_fraction, 0.0, 0.95)))
    A = population_activity(S[start:], bin_steps)
    theta = int(np.floor(np.percentile(A, threshold_percentile))) if len(A) else 0
    sizes, durations = detect_avalanches(A, theta)
    res = dict(activity=A, threshold=theta, sizes=sizes, durations=durations,
               burn_in_steps=start, bin_steps=int(bin_steps),
               threshold_percentile=float(threshold_percentile),
               mean_activity=float(A.mean()) if len(A) else np.nan)
    res['size_fit'] = fit_discrete_power_law(sizes)
    res['duration_fit'] = fit_discrete_power_law(durations)
    res['scaling'] = size_vs_duration_exponent(sizes, durations)
    tau, alpha = res['size_fit']['alpha'], res['duration_fit']['alpha']
    res['gamma_predicted'] = float((alpha - 1) / (tau - 1)) if np.isfinite(tau) and np.isfinite(alpha) and tau > 1 else np.nan
    res['branching'] = mr_branching_ratio(A)
    return res


def summary(res: dict) -> dict:
    """JSON-friendly scalar summary."""
    sf, df_ = res['size_fit'], res['duration_fit']
    return dict(
        n_avalanches=int(len(res['sizes'])), threshold=int(res['threshold']),
        threshold_percentile=res['threshold_percentile'], bin_steps=res['bin_steps'],
        burn_in_steps=int(res['burn_in_steps']), mean_activity_per_bin=res['mean_activity'],
        size_exponent_tau=sf['alpha'], size_exponent_se=sf.get('alpha_se', np.nan),
        size_xmin=sf['xmin'], size_ks=sf['ks'], size_n_tail=sf['n_tail'],
        duration_exponent_alpha=df_['alpha'], duration_exponent_se=df_.get('alpha_se', np.nan),
        duration_xmin=df_['xmin'], duration_ks=df_['ks'], duration_n_tail=df_['n_tail'],
        gamma_fit=res['scaling']['gamma'], gamma_predicted=res['gamma_predicted'],
        branching_ratio_m=res['branching']['m'],
        branching_offset_c=res['branching'].get('c', np.nan),
        lag1_autocorrelation=res['branching'].get('r1', np.nan),
        max_size=int(res['sizes'].max()) if len(res['sizes']) else 0,
        max_duration=int(res['durations'].max()) if len(res['durations']) else 0,
    )


# --------------------------------------------------------------------------- plotting
def _log_binned_pdf(x, n_bins=25):
    x = np.asarray(x, float)
    x = x[x >= 1]
    if len(x) == 0 or x.max() <= 1:
        return np.array([]), np.array([])
    edges = np.unique(np.floor(np.geomspace(1, x.max() + 1, n_bins)))
    edges = np.append(edges, edges[-1] + 1) if len(edges) < 2 else edges
    hist, edges = np.histogram(x, bins=edges)
    width = np.diff(edges)
    pdf = hist / (width * len(x))
    centers = np.sqrt(edges[:-1] * np.maximum(edges[1:] - 1, edges[:-1]))
    ok = hist > 0
    return centers[ok], pdf[ok]


def plot(res: dict, dt: float = 1.0, title: str = "Avalanche analysis"):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    A, theta = res['activity'], res['threshold']

    ax = axes[0, 0]
    n_show = min(len(A), 600)
    t = (np.arange(n_show) * res['bin_steps'] + res['burn_in_steps']) * dt
    ax.plot(t, A[:n_show], lw=0.8, color='0.3')
    ax.axhline(theta, color='C3', ls='--', lw=1, label=f"threshold = {theta} ({res['threshold_percentile']:.0f}th pct)")
    ax.fill_between(t, theta, A[:n_show], where=A[:n_show] > theta, color='C3', alpha=0.25, step=None)
    ax.set_xlabel('time'); ax.set_ylabel('spikes / bin'); ax.legend(loc='upper right', fontsize=8)
    ax.set_title('Population activity (excerpt after burn-in)')

    for ax, key, fitkey, lab in [(axes[0, 1], 'sizes', 'size_fit', 'S'), (axes[1, 0], 'durations', 'duration_fit', 'T')]:
        x = res[key]
        c, p = _log_binned_pdf(x)
        ax.loglog(c, p, 'o', ms=4, color='C0', label=f'{lab} (n={len(x)})')
        f = res[fitkey]
        if np.isfinite(f['alpha']):
            xs = np.arange(int(f['xmin']), max(int(f['xmin']) + 2, int(np.max(x)) + 1))
            frac = f['n_tail'] / max(1, len(x))
            ax.loglog(xs, frac * xs ** (-f['alpha']) / zeta(f['alpha'], f['xmin']), 'r--',
                      label=f"MLE: exponent {f['alpha']:.2f} ± {f.get('alpha_se', np.nan):.2f}, "
                            f"xmin={f['xmin']}, KS={f['ks']:.3f}")
        ax.set_xlabel(lab); ax.set_ylabel(f'P({lab})'); ax.legend(fontsize=8)
        ax.set_title(f"{'Size' if lab == 'S' else 'Duration'} distribution (log-binned)")

    ax = axes[1, 1]
    sc = res['scaling']
    ax.loglog(sc['T'], sc['meanS'], 'o', ms=4, color='C2', alpha=0.6, label='<S>(T)')
    if np.isfinite(sc['gamma']):
        Tu = sc['T'][sc['used']]
        ax.loglog(Tu, np.exp(sc['log_c']) * Tu ** sc['gamma'], 'k-', label=f"fit γ = {sc['gamma']:.2f}")
        if np.isfinite(res['gamma_predicted']):
            ax.loglog(Tu, np.exp(sc['log_c']) * Tu ** res['gamma_predicted'], 'r:',
                      label=f"predicted (α-1)/(τ-1) = {res['gamma_predicted']:.2f}")
    ax.set_xlabel('T'); ax.set_ylabel('<S>'); ax.legend(fontsize=8)
    ax.set_title(f"Size vs duration   |   MR branching ratio m = {res['branching']['m']:.3f}")

    fig.suptitle(title)
    fig.tight_layout()
    return fig
