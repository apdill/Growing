"""
Shared LIF update used by the plastic simulation (run_sim) and the frozen-network probes.

One step:
    V += spike_input * M @ spikes + external drive + (-leak (V - rest) [+ uniform drive]) dt
then spiking by either
    'threshold' : V >= spike_threshold                         (legacy, hard threshold)
    'escape'    : P(spike) = 1 - exp(-dt * escape_rate0 * exp((V - spike_threshold) / escape_delta))
                  (escape noise / soft threshold; makes the population transition continuous,
                   branching-process-like, instead of the bistable quiescent/saturated switch of a
                   hard-threshold LIF with slow leak)
followed by reset to rest and an optional absolute refractory period (refractory_steps).

`rng` may be the np.random module (legacy global stream) or a np.random.Generator; the call
order of random draws for the legacy configuration ('threshold', no refractoriness) is
unchanged from the original inline code.
"""
from __future__ import annotations

import numpy as np


def _rand(rng, n):
    return rng.rand(n) if hasattr(rng, "rand") else rng.random(n)


def lif_step(V, spikes, M, cfg, rng, refractory=None, extra_input=None):
    """Advance one step in place-free style. Returns (V, spikes, refractory)."""
    N = V.shape[0]
    rec = cfg.spike_input * (M @ spikes.astype(float))
    if getattr(cfg, "external_drive", "uniform") == "poisson":
        ext = cfg.ext_kick * rng.poisson(cfg.ext_rate * cfg.dt, N)
        dV = -cfg.leak * (V - cfg.rest)
        V = V + rec + ext + dV * cfg.dt
    else:
        dV = cfg.input_noise * _rand(rng, N) - cfg.leak * (V - cfg.rest)
        V = V + rec + dV * cfg.dt
    if extra_input is not None:
        V = V + extra_input

    ref_steps = int(getattr(cfg, "refractory_steps", 0) or 0)
    if ref_steps > 0:
        if refractory is None:
            refractory = np.zeros(N, dtype=np.int32)
        blocked = refractory > 0
        V[blocked] = cfg.rest

    mode = getattr(cfg, "spike_mode", "threshold")
    if mode == "escape":
        z = np.clip((V - cfg.spike_threshold) / max(cfg.escape_delta, 1e-9), -50.0, 30.0)
        p = 1.0 - np.exp(-cfg.dt * cfg.escape_rate0 * np.exp(z))
        new_spikes = _rand(rng, N) < p
    else:
        new_spikes = V >= cfg.spike_threshold

    if ref_steps > 0:
        new_spikes &= ~blocked
        refractory = np.maximum(refractory - 1, 0)
        refractory[new_spikes] = ref_steps

    V[new_spikes] = cfg.rest
    V[V < cfg.rest] = cfg.rest
    return V, new_spikes, refractory
