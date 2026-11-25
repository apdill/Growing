# Plan: fractal-geometry sweeps (not yet run)

Goal: find out whether *measured* fractal geometry of the grown neurons / network influences
criticality and function in the regimes where the model is closest to critical, rather than
varying the `branch.D` knob alone (s10 did that; it only sets the per-order branch-length ratio
r = 2^(-1/D), and homeostasis resized arbors until coupling matched, hiding any effect).

Background (see `reports/Growing critical networks literature.md`, §"Fractal dimension is a weak lever"):
* Real 2-D projected dendrites: box-counting D ≈ 1.1–1.7 over ~1–1.5 decades; D > 2 is impossible
  in a plane → drop D = 2.5.
* Measured D is driven mostly by tortuosity ("weave") and branch angles, not branch-length scaling.
* No published link between dendritic D and network criticality / dynamic range → genuine gap,
  but the expected answer under homeostasis is "no effect on the set point"; effects, if any, should
  appear in *fluctuations*, *dynamic range*, *modularity / distance-dependence of connections*.

## Step 0 — code prerequisites (≈ half a day)

1. `fractal.py`
   * `box_count_D(mask, scales=None)` → D, fit range, R² (log-spaced box sizes from ~2 px to
     arbor span / 4; report the scaling range; also `sandbox_D` / mass-radius as a cross-check).
   * per-neuron D from each neuron's own mask (count layer > 0), whole-network D from the union mask,
     plus lacunarity (gliding box) — lacunarity often separates arbors with equal D.
   * validation: Koch / Vicsek / random-walk test images with known D.
2. Record in `sweep.py` metrics: `D_neuron_mean/sd`, `D_network`, `lacunarity`, `tortuosity_mean`
   (path length / end-to-end per segment), `branch_points_per_neuron`, `total_length_mean`,
   `span_mean`, connection probability vs soma distance (fit exp decay length λ_c),
   modularity (Louvain on M), Λ = largest eigenvalue of P = min(1, spike_input·M).
3. Open-loop growth arm (needed for F3): `lif.growth_rule = "open_loop"` → every neuron grows to a
   fixed `target_total_length_px` with no activity feedback; then freeze. Probes rescale the gain so
   that Λ is identical across conditions (`scan_frozen.py --match-lambda`), so geometry is compared at
   matched distance from criticality.

## F1 — does the morphology knob set move measured D? (cheap, morphology only, no dynamics)

Grow single neurons / small sparse networks (no plasticity, `lif.duration` short) over
`branch.curviness_magnitude` {0.3, 0.8, 1.2, 2.0} × `branch.branch_angle` {π/8, π/4.5, π/3} ×
`branch.D` {1.3, 1.6, 1.9} × `branch.randomness` {0.2, 0.45, 0.8}; measure box-counting D.
Output: a calibration map (knobs → measured D, lacunarity, span, total length). Pick 4–5 settings that
span measured D ≈ 1.2–1.7 at roughly matched total length. ~108 cheap runs.

## F2 — measured D in the near-critical regimes (homeostatic growth)

Use the F1 settings (D_measured ≈ 1.2 / 1.35 / 1.5 / 1.65) in two regimes:
* **F2a, large-N limit:** s9 regime (N·h fixed), N = 100 and 200, 2 seeds → 16 runs (~1 h each
  for N = 200, ~4 GB each: run ≤ 3 N=200 runs at once).
* **F2b, slow growth:** s8 regime with `lif.length_growth_rate` 0.1, `lif.rate_tau` 200, N = 50,
  4 seeds → 16 runs.
Measure: m, Λ, s.d. σ(t), σ-conditioned avalanche ratio, window exponents, dynamic range (with
coupling scan), D_neuron, D_network, λ_c, modularity, clustering.
Hypotheses: set point unchanged (homeostasis absorbs geometry); fluctuation size and dynamic range
change via contact granularity (dense, tortuous arbors → more, smaller contacts → smaller σ quanta).

## F3 — geometry at matched gain (open loop, frozen)

Open-loop growth to fixed total length {600, 1200} px with the F1 settings, N = 100, 3 seeds;
freeze; rescale to Λ ∈ {0.9, 0.97, 1.0} of the critical gain; measure dynamic range, avalanche
statistics, modularity, λ_c. Isolates "does geometry matter at the same distance from criticality".

## Analysis / figures
* measured D vs knob (F1 calibration heat-maps)
* D_neuron, D_network, lacunarity vs: s.d. σ, m, Λ, dynamic range, modularity, λ_c (F2, F3)
* partial correlations controlling for total length / span

Draft specs: `sweep_specs/drafts/f1_morph_calibration.json`, `f2a_fractal_fss.json`,
`f2b_fractal_slow.json` (F3 spec after step 0.3 exists). The `fractal.*` metrics and the open-loop
rule do not exist yet — implement step 0 first, then move specs out of `drafts/` and run e.g.
`python sweep.py sweep_specs/f2b_fractal_slow.json --workers 8 --burn-in 0.2`.
Fill the `branch.*` values for F2 from the F1 calibration (placeholders marked "FROM_F1").
