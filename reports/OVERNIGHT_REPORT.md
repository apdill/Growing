# Overnight sweeps: when does the growing network become critical?

Branch-morphology homeostatic plasticity model (Kossio-style structural plasticity with realistic
dendrites). All numbers below come from `sweeps/<name>/results.csv`; figures are in `reports/`.
Each sweep is reproducible with `python sweep.py sweep_specs/<name>.json`.

## TL;DR

1. **Homeostatic growth drives the network towards criticality exactly as the driven-branching
   argument predicts.** With sparse external drive h and target rate r*, the grown network settles at
   branching ratio m ≈ 1 − h/r* (s1). The weaker the drive relative to the target, the closer to m = 1.
2. **The original neuron model cannot be critical.** The hard-threshold LIF with slow leak (legacy
   default) has a *discontinuous* (bistable) transition: frozen grown networks either fall silent or
   saturate (rate jumps 0.004 → 1.6 between 0.75× and 1× coupling). Homeostasis just keeps the
   plastic network flipping across that edge. Refractoriness turns the jump into a continuous
   transition; the cleanest model is probabilistic transmission (Haldeman–Beggs / Kinouchi–Copelli),
   in which the morphology directly sets transmission probabilities p_ij.
3. **With that model, three sandbox settings decide whether criticality is reachable:**
   * *contact granularity* (`contact_ref_px`): one dendrite crossing must change σ by a small amount,
     otherwise homeostasis can only jump between sub- and supercritical (σ s.d. 0.2–0.6 at 50 px vs
     0.02 at 200 px);
   * *pass-through probability*: with the original competition rule (tips mostly stop at the first
     contact) arbors stall at ~25 % coverage and never reach the target rate; at 0.6–0.9 they do;
   * *drive ratio* h/r*.
4. **Self-organised critical avalanches during growth (s6, long plastic runs, N = 50):** at
   h/r* = 0.01 all three seeds give τ = 1.43–1.46, α = 1.60–1.74, and the crackling-noise relation
   holds (γ_fit = 1.43–1.45 vs (α−1)/(τ−1) = 1.39–1.59), with m ≈ 0.96. At h/r* = 0.04 the network
   is subcritical (power law only up to ~100, steep cutoff); at h/r* = 0.002 it overshoots
   (system-spanning bursts, bump in P(S)).
5. **Function:** dynamic range rises towards the critical point; grown networks sit just on the
   subcritical side of the optimum (the "reverberating regime", m ≈ 0.9–0.98, also reported in vivo by
   Wilting & Priesemann 2018).

6. **Finite size (s7): is it real criticality or quasi-criticality?** The power law during growth is
   robust for N ≥ 50 (τ = 1.43–1.45, α = 1.6–1.75) and its cutoff grows with N (s_c ~ N^0.5), but
   the coupling σ(t) is *not* pinned: homeostasis makes it oscillate around the set point, and
   avalanches that start in high-σ periods are 10–50× larger than those in low-σ periods. Part of
   the pooled power law is therefore a mixture across sub- and supercritical states
   ("self-organised quasi-criticality", Bonachela & Muñoz 2009). The fluctuations shrink faster
   than N^−1/2 (s.d. σ ~ N^−0.66), and so does the mixing (high/low ⟨S⟩ ratio 50 → 11 from
   N = 25 to 200): **the larger the network, the more tightly growth holds it at one point.**
   That point stays at a fixed distance from m = 1 (1 − m ≈ 0.035 at every N), set by the drive.

7. **Function at the set point (s7 coupling scans):** in frozen copies of the grown networks
   (N = 100, 200), the peak of the dynamic range, the peak of the branching ratio and the point where
   the frozen rate equals the target all fall at the same coupling, and **the mean coupling the network
   grows to sits exactly there**. Homeostatic growth places the network at its point of maximal
   sensitivity; at N = 25–50 the curves are too noisy to say this.
8. **Timescales (s8): the clean 3/2 law is tied to *fast* plasticity.** With slow growth (σ pinned,
   s.d. 0.1) the network sits at m ≈ 0.956 and P(S) has a shallow body (slope ≈ 0.95) with a hard
   cutoff, not τ = 3/2. Faster growth (σ s.d. 0.3–0.45) gives slopes of 1.25–1.35 (KS fits τ ≈ 1.4)
   and satisfies crackling-noise scaling. The 3/2-like statistics come with strong σ sweeping
   (high/low-σ ⟨S⟩ ratio 30–80 vs ~6 when slow). Two readings, both plausible: (i) quasi-criticality
   (mixing sub/supercritical epochs), and (ii) the slow-growth network is at a stationary point that
   is near-critical but *driven* (avalanches overlap, Das & Levina 2019 "relaxed timescale
   separation"), which flattens the body. s9 (h ∝ 1/N) is set up to separate these.

<!-- S7 / S8 sections appended below -->

## What was changed in the code (all opt-in; legacy defaults are bit-identical)

| Area | Change |
|---|---|
| speed | marker clean-up only at the end (was 40 % of run time); uint16 count layers (½ memory → more parallel runs); `morphology_update_interval` k > 1 (activity every step, morphology every k steps with k·dt growth; ~2× faster, explicit separation of timescales) |
| drive | `external_drive="poisson"`: sparse kicks at rate `ext_rate` (spontaneous rate h) |
| neuron | `lif_dynamics.lif_step` shared by plastic run and probes; `refractory_steps`; `spike_mode` = `threshold` (legacy) / `escape` (soft threshold) / `branching` (probabilistic transmission) |
| analysis | `sweep.py` (parallel headless sweeps, resumable), `probes.py` (dynamic range, coupling spectrum, graph stats, windowed m, frozen-network avalanches), `scan_frozen.py` (coupling scan), `analyze_sweeps.py`, `figure_soc.py`; coupling trace σ(t) recorded during runs |

## Sweep log

### s1 — drive ratio, legacy LIF (Poisson drive, 7 levels × 3 seeds) → `reports/s1_drive_ratio.png`
Plastic m follows 1 − h/r* (slightly above it); λ_max(W) grows from 2.0 to 3.4 as h falls; dynamic
range ~16 dB → ~32–43 dB at the lowest drives. **But** frozen grown networks are bistable (frozen
rate 0.002 … 1.77), exponents steep (τ ≈ 3–8). Diagnosis by coupling scan of a frozen network:
discontinuous transition.

### s1b — legacy uniform drive
Positive mean drive makes every neuron fire on its own; recurrent propagation is negligible (m ≈ 0).
Not a route to criticality.

### s2 — neuron model (refractoriness 0/2/4, escape noise) × drive
Refractoriness/escape noise remove the bistability (frozen rate ≈ plastic rate, frozen m ≈ plastic m)
but networks settle subcritical (m 0.7–0.9): with slow leak a neuron integrates many weak inputs,
activity is sustained by mean input rather than by spike-to-spike propagation, and rate no longer
equals h/(1 − m), so homeostasis stops pinning m. Best LIF case (threshold + 2 refractory steps,
h/r* = 0.02): τ ≈ 1.6–2.9, α ≈ 2.9–3.2, avalanches up to 10⁴.

### s3 — probabilistic transmission, `contact_ref_px` = 50
Homeostasis works on average, but σ moves in large quanta (one crossing ≈ +0.2–0.5); frozen
snapshots are either dead or hot.

### s4 — contact granularity 50/200/800
200 px: σ ≈ 1 with s.d. 0.02 and frozen m ≈ plastic m — but arbors stall (~260 px, 25 % coverage)
and the rate stays far below target: tips stop at contacts (competition rule).

### s5 — pass-through 0.2/0.6/0.9 × granularity 100/200
Pass-through ≥ 0.6 lets arbors keep growing; rate reaches target and plastic m ≈ 0.94–0.99.
Coupling scan of the grown networks (`reports/s5_coupling_scan.png`): one network gives τ = 1.63,
α = 1.92 at its own coupling (and 1.53–1.69 / 1.77–2.02 over g = 0.9–1.2); dynamic range increases
with σ up to the transition and the grown networks (★) sit just below it.

### s6 — long plastic runs (60 000 time units, k = 10), h/r* = 0.04 / 0.01 / 0.002 → `reports/s6_soc_avalanches.png`

| h/r* | m | τ | α | γ fit / predicted | avalanches | max size |
|---|---|---|---|---|---|---|
| 0.04 | 0.92 | 3.1–4.2 | 2.4–4.4 | 1.35 / 0.7–1.6 | ~2 900 | ~1 000 |
| **0.01** | **0.96** | **1.43–1.46** | **1.60–1.74** | **1.43–1.45 / 1.39–1.59** | ~930 | ~2 700 |
| 0.002 | 0.98 | 1.32–1.33 | 1.41–1.43 | 1.26–1.33 / 1.28 | ~220 | ~4 500 |

Rates at target (0.043–0.049 vs 0.05) in all runs.

### s7 — finite size, N = 25/50/100/200 at constant density, h/r* = 0.01 → `reports/s7_finite_size_fss.png`, `reports/s7_coupling_scans.png`

| N | avalanches | τ (KS) | α | m | ⟨S²⟩/⟨S⟩ | s.d. σ(t) | ⟨S⟩ high-σ / low-σ tertile |
|---|---|---|---|---|---|---|---|
| 25 | 490 | 2.46 (xmin 292) | 2.67 | 0.968 | 650 | 0.52 | 50 |
| 50 | 935 | 1.45 | 1.67 | 0.963 | 990 | 0.33 | 37 |
| 100 | 1 490 | 1.43 | 1.75 | 0.964 | 1 500 | 0.21 | 16 |
| 200 | 2 710 | 1.43 | 1.60 | 0.965 | 1 980 | 0.13 | 11 |

* Cutoff grows as N^0.54 (⟨S²⟩/⟨S⟩) / N^0.51 (99.9th percentile). Caveat: h is fixed, so N·h (external
  events per step) grows with N and avalanches overlap more in large networks; part of the cutoff
  growth may be merging. s9 holds N·h fixed.
* s.d. of σ(t) ~ N^−0.66: growth-driven fluctuations of the coupling vanish with size.
* 1 − m ≈ 0.035 at every N: the distance from m = 1 is set by drive and refractoriness, not N.
* Frozen networks show steep KS exponents (τ 3–6) at every g. A single snapshot of σ is a random draw
  from the oscillation (final σ is up to 1.3 s.d. from the mean), and KS fits pick xmin in the
  cutoff; frozen exponents are not a reliable criticality test here. The scans are: rate rises smoothly
  with σ (no bistability), m peaks at ≈ 0.96 where the rate reaches target, then falls (saturation +
  refractoriness), and Δ peaks at the same place for N = 100/200 (28 dB) — where the network grew to.

### s8 — growth timescales: `length_growth_rate` 0.1/0.5/2 × `rate_tau` 50/200/800 (N = 50) → `reports/s8_timescales_timescales.png`, `reports/s8_timescales.png`

| growth rate | s.d. σ | m | slope P(S), S∈[2,200] | ⟨S⟩ high/low σ | rate / target | Δ (dB) |
|---|---|---|---|---|---|---|
| 0.1 | 0.09–0.13 | 0.955–0.962 | 0.92–0.97 | 5–7 | 0.91–0.95 | 29–33 |
| 0.5 | 0.27–0.39 | 0.956–0.974 | 1.14–1.31 | 26–64 | 0.95–0.97 | 25–31 |
| 2.0 | 0.38–0.45 | 0.961–0.976 | 1.25–1.34 | 31–82 | 0.50–0.89 | 20–32 |

* Slow growth: tight homeostasis, highest and most consistent dynamic range, but no 3/2 law.
* Fast growth: σ overshoots, power laws look "critical" (τ_KS ≈ 1.4, α ≈ 1.6, γ_fit ≈ γ_pred ≈ 1.4–1.5),
  but these are mixtures; at growth rate 2 with short rate_tau, homeostasis itself fails (rate 50 %
  of target, relaxation oscillations).
* rate_tau (activity-sensor timescale) matters less than the growth rate; a long sensor with fast
  growth gives the largest σ swings (the classic delayed-feedback oscillation).
