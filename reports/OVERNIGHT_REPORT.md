# Overnight sweeps: when does the growing network become critical?

Branch-morphology homeostatic plasticity model (Kossio-style structural plasticity with realistic
dendrites). All numbers below come from `sweeps/<name>/results.csv`; figures are in `reports/`.
Each sweep is reproducible with `python sweep.py sweep_specs/<name>.json`.

## TL;DR (updated 06:30 after s7–s10)

**Headline: homeostatic dendritic growth does self-organise to criticality, but in the proper SOC
limit — N → ∞ with the drive per neuron scaled down (h ∝ 1/N, s9). In finite networks it sits just
subcritical and wanders across the critical point (quasi-criticality). The grown network's set
point coincides with its dynamic-range optimum.**

1. **Set point.** With sparse drive h and target rate r*, the grown network settles at m ≈ 1 − h/r*
   (s1). The weaker the drive, the closer to criticality.
2. **Neuron model.** The legacy hard-threshold LIF has a discontinuous (bistable) transition and
   cannot be critical. Probabilistic transmission (Haldeman–Beggs / Kinouchi–Copelli; morphology sets
   p_ij) has a continuous one (s2–s3).
3. **Sandbox settings needed** (s3–s5): fine contact granularity (`contact_ref_px` ≥ 100–200, so one
   crossing changes σ only a little) and pass-through ≥ 0.6 (otherwise arbors stall below target).
4. **⚠ Exponent correction.** The earlier "τ = 1.43–1.46" numbers are Clauset/KS fits whose xmin
   (3–11) cuts off the body. A fixed-window MLE (S ∈ [2, 200]; `avalanches.fit_truncated_power_law`)
   shows the distributions are curved: body slope 1.0–1.3 in most finite runs. On frozen networks, KS
   picks xmin in the cutoff region (xmin 10²–10⁴, 50–200 events) → the "τ = 3–6" frozen values are
   artifacts. Window benchmark for exact critical branching: τ = 1.50, α = 1.64
   (`reports/trunc_exponents.csv`, `reanalyze.py`).
5. **Finite-size scaling — the main result (s9 vs s7, `reports/s9_vs_s7.png`).**
   * h/r* fixed (s7): m stays at 0.965 for all N; τ_window 1.01 → 1.29; cutoff ~ N^0.54.
   * N·h fixed, h/r* = 0.5/N (s9): m = 0.948 → 0.962 → 0.969 → **0.988**; τ_window 1.10 → 1.19 →
     1.32 → **1.41**; α_window 1.19 → 1.55 (critical 1.64); cutoff ~ N^1.3 (max avalanche 10⁴ at
     N = 200). **The statistics converge on critical branching as N grows.**
6. **Quasi-criticality in finite networks (s7, s8, low-drive scans).** Homeostasis makes the coupling
   σ(t) oscillate around the set point; avalanches starting in the top σ-tertile are 10–80× larger
   than in the bottom one. Frozen copies driven 10× more weakly (`reports/s8_frozen_lowdrive_pS.png`)
   put the true critical coupling about 15–20 % above the grown mean. Below it P(S) has an
   exponential cutoff; above it there is a bump of system-spanning bursts. The pooled distribution
   during growth mixes the two (Bonachela & Muñoz 2009 "SOqC"). Fluctuations shrink with N at fixed h
   (s.d. σ ~ N^−0.66) but not when h ∝ 1/N (fewer drive events → noisier rate sensor). Across all
   branching sweeps, larger σ fluctuations → steeper, more "critical-looking" body
   (`reports/cross_sweep_tau.png`).
7. **Timescales (s8).** Slow growth pins σ (s.d. 0.1): m ≈ 0.956, shallow body (τ_window ≈ 0.95),
   highest and most consistent dynamic range (29–33 dB). Fast growth: σ swings (s.d. 0.3–0.45),
   3/2-like statistics, but lower dynamic range. At growth rate 2 with a short sensor, homeostasis
   fails (rate 50 % of target, relaxation oscillations). **"Pretty" power laws and good function are
   different things here:** across all runs the dynamic range anticorrelates with the body exponent
   (r = −0.59).
8. **Function at the set point (s7 coupling scans, `reports/s7_coupling_scans.png`).** For N = 100/200,
   the frozen network's dynamic-range peak, branching-ratio peak and rate = target point coincide,
   and the mean coupling the network grows to sits right there. Growth places the network at its
   point of maximal sensitivity, on the subcritical side of the h → 0 critical point (the
   "reverberating regime", m ≈ 0.95–0.99, cf. Wilting & Priesemann 2018).
9. **Morphology barely matters (s10).** Fractal dimension D ∈ {1.3, 1.8, 2.5} × base inter-branch
   length ℓ0 ∈ {25, 60, 150} px: m 0.957–0.965, τ_window 1.17–1.24, rate on target for all nine.
   Graph structure also barely moves: degree 6–9, clustering 0.51–0.58, small-world σ ≈ 2.3–2.8.
   The critical state is set by the homeostatic rule and drive, not by arbor geometry — a
   universality-type result. Dynamic range varies 26–34 dB with no consistent trend over 2 seeds.

### Suggested next steps
* Push N further with h ∝ 1/N (N = 400 needs ~16 GB with per-neuron count layers; would need
  sparse/shared count layers first).
* A slower or better-averaged rate sensor at small h (rate_tau ∝ 1/h) to cut σ fluctuations in the
  s9 limit — tests whether quasi-criticality vanishes too.
* Topology at criticality: compare grown graphs with degree-preserving rewired controls (dynamic
  range, m) to see whether growth builds something functionally special beyond the σ value.

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

### s9 — finite-size scaling with N·h fixed (h = 0.0005·50/N) → `reports/s9_fss_scaled_drive_fss.png`, `reports/s9_vs_s7.png`

| N | h/r* | m | τ window [2,200] | α window [2,50] | τ KS | ⟨S²⟩/⟨S⟩ | max S | s.d. σ | Δ (dB) |
|---|---|---|---|---|---|---|---|---|---|
| 25 | 0.02 | 0.948 | 1.10 | 1.19 | 2.28 | 370 | 1 400 | 0.39 | 20.5 |
| 50 | 0.01 | 0.962 | 1.19 | 1.27 | 1.47 | 990 | 2 500 | 0.33 | 29.3 |
| 100 | 0.005 | 0.969 | 1.32 | 1.44 | 1.38 | 2 440 | 4 240 | 0.29 | 29.8 |
| 200 | 0.0025 | 0.988 | 1.41 | 1.55 | 1.33 | 5 320 | 10 850 | 0.30 | 28.2 |

Critical-branching benchmark in the same windows: τ = 1.50, α = 1.64. 1 − m ~ N^−0.66; cutoff ~ N^1.3
(99.9th percentile ~ N^0.98). 2 seeds per N.

### s10 — morphology scaling constants (D × ℓ0, N = 50, h/r* = 0.01) → `reports/s10_morphology.png`

| D \ ℓ0 | 25 | 60 | 150 |
|---|---|---|---|
| 1.3 | m 0.965, τ_w 1.21, Δ 27.8 | 0.964, 1.22, 30.7 | 0.961, 1.17, 34.4 |
| 1.8 | 0.961, 1.24, 26.2 | 0.961, 1.23, 28.4 | 0.959, 1.23, 32.5 |
| 2.5 | 0.964, 1.19, 32.6 | 0.962, 1.18, 30.6 | 0.957, 1.17, 25.8 |

Rates on target everywhere; arbor radius 300–375 px; coverage 26–29 %.

### Analysis tools added overnight
`finite_size.py` (FSS + σ-conditioned avalanche test), `timescales.py`, `reanalyze.py` (window
exponents for all runs → `reports/trunc_exponents.csv`), `cross_sweep.py`, `scan_frozen.py --ext-rate
--g` (low-drive coupling scans), `reports/make_*.py` (figure scripts).

### External check — Stiber-lab simulated culture (Zenodo 17220617, run tR_1.0--fE_0.98) → `reports/stiber_tR_1.0--fE_0.98.png`
Disc-growth model (van Ooyen calcium rule) + LIF with 98 % excitatory cells, short-term synaptic
depression/facilitation, endogenously active pacemakers; 10 000 neurons, 28 DIV equivalent. The
deposit only contains burst epochs (runs of 10 ms bins with > 50 network spikes = thresholded
avalanches, θ = 50), so only those were analysed (`external_analysis/stiber_bursts.py --frames`).
* **No scale-free statistics at all:** all 9 160 bursts are system-spanning and stereotyped —
  S = 42 000–56 000 spikes (CV 1–5 %), 99.5 % of the 10 000 neurons take part, ~5 spikes/neuron,
  T = 110–370 ms; P(S) occupies one or two log bins (power-law fits meaningless), S does not grow
  with T (γ ≈ 0). Quasi-periodic inter-burst interval ≈ 4.8 s throughout.
* Within a burst, activity grows by ×1.16–1.26 per 10 ms (median rising-phase ratio; ×2 at onset)
  → supercritical propagation, terminated by synaptic resource depletion.
* Development: first burst at ~4.5 DIV; bursts get larger and shorter (T 23 → 14.5 bins) and sharper,
  with a regime change at ~13–15 DIV (size distribution becomes bimodal).
* **Homeostasis is satisfied by bursting:** burst spikes alone give 0.93 Hz per neuron in the second
  half vs the outgrowth set point of 1.0 Hz. With little spontaneous drive, rate homeostasis is met by
  periodic system-wide bursts — the low-input "bursting" regime of Zierenberg et al. 2018 — not by
  criticality. Our model avoids this because probabilistic transmission without depression has a
  continuous transition, so the set point lands at m ≈ 1 − h/r* instead.
* Implication: if we add inhibition + synaptic depression for culture realism (needed to match
  Wagenaar-type data), growth may tip into this bursting regime; the drive ratio h/r* and depression
  strength then decide critical vs bursting — a concrete, testable axis for the project.
