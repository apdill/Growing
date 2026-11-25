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

## Step 0 — code prerequisites ✅ done (2026-11-25)

* `fractal.py` — box counting via FracStack (`../FracStack`; loaded without its torch-dependent
  package `__init__`) with **min counts over grid offsets and unweighted OLS fits**:
  * `D_fixed`: fixed common range (neurons 4–96 px; network 4 px – canvas/8); `D_dyn`: FracStack
    `dynamic_boxcount` best-R² window ≥ 1 decade. Largest box ≤ extent/8 limits finite-object edge bias.
  * validation (`python fractal.py`): line 0.981 / 0.997 (dyn), Koch 1.289 (1.262), Sierpinski carpet
    1.844 (1.893), filled square 1.957 / 1.986 → expect a small negative bias for dense objects.
  * per neuron: `D_neuron` (image-like: filled dendrites + soma, comparable to published neuron
    box-counting) and `D_neuron_skel` (arbor skeleton only); whole network: `D_network_fill`,
    `D_network_skel`; gliding-box lacunarity (neuron at 16 px, network at 4–64 px); per-branch
    tortuosity (path / chord per tip lineage); branches and cable length per neuron.
  * graph geometry: P(connection) vs soma distance → decay length `conn_decay_len`; Louvain
    `modularity`; `lambda_max_P` = largest eigenvalue of P = min(1, spike_input·M).
* `sweep.py` records all of the above in `results.csv` (+ `morphology_metrics.npz` per run with
  per-neuron arrays). Spec keys: `"fractal": false` to skip, `"fractal_dynamic": false` to skip the
  dynamic fit; `_`-prefixed keys are ignored; unfilled `FROM_*` placeholders abort the run.
  Cost: ~1.4 s per neuron (≈ 5 min for N = 200).
* `lif.growth_rule = "open_loop"` + `lif.target_total_length_px`: no activity feedback, grows at
  `length_growth_rate` to the cap. NB: in this code "total length" is the per-path length cap
  (≈ arbor radius), not cable length — cable length is measured (`cable_length_mean`).
* `scan_frozen.py --match-lambda L1 L2 ...`: picks the gain so Λ(P) hits each value. Grown networks
  sit at Λ ≈ 2.4 (refractoriness + coalescence push the effective critical Λ to ~2–2.6, not 1), so
  for F3 first locate each network's critical gain (plain `--g` scan, peak of m / dynamic range)
  and then compare at matched Λ values bracketing it (e.g. `--match-lambda 1.8 2.1 2.4 2.7`).
* `fractal_calibration.py`: F1 → knob/geometry correlation table, calibration figure, and the 4
  settings hitting target D at roughly matched cable length, printed as spec `runs` entries.

Smoke-test findings (`sweep_specs/smoke_fractal.json`):
* `branch.trajectory_deviation` is the real tortuosity knob (with trajectory planning on,
  `randomness`·`curviness_magnitude` only add small noise). ≥ 20 stalls growth (tips curl into
  their own branches and are retracted as stubs) → F1 uses 4–16.
* In 4 smoke runs D_neuron only spanned 1.30–1.34; `branch.D` raised it, deviation lowered branch
  count and cable length. Whether 1.2–1.65 is reachable is exactly what F1 answers — if not,
  widen `branch_length_scale` / `n_primary_dendrites` / `mean_branches` in F1.

## Launch sequence

```bash
mv sweep_specs/drafts/f1_morph_calibration.json sweep_specs/
python sweep.py sweep_specs/f1_morph_calibration.json --workers 12            # ~10 min
python fractal_calibration.py f1_morph_calibration --targets 1.2 1.35 1.5 1.65
# paste the printed runs into f2a/f2b/f3 (replace FROM_F1 entries), move them out of drafts/, then
python sweep.py sweep_specs/f2b_fractal_slow.json --workers 8 --burn-in 0.2   # ~1 h
python sweep.py sweep_specs/f2a_fractal_fss.json --workers 3 --burn-in 0.2    # N=200: ~4 GB/run
python sweep.py sweep_specs/f3_matched_gain.json --workers 8
python scan_frozen.py sweeps/f3_matched_gain/runs/<run> --g 0.5 0.7 0.8 0.9 1 1.1 1.2   # locate critical gain
python scan_frozen.py sweeps/f3_matched_gain/runs/* --match-lambda <values> --out reports/f3_matched.png
```

## F1 — does the morphology knob set move measured D? (cheap, morphology only, no dynamics)

9 sparse neurons, open-loop growth to a 500 px path cap: `branch.trajectory_deviation` {4, 8, 12, 16}
× `branch.branch_angle` {π/8, π/4.5, π/3} × `branch.D` {1.3, 1.6, 1.9} × `branch.branch_length_scale`
{30, 60, 120} × `branch.randomness` {0.2, 0.8} = 216 runs (~30 s each). Output: calibration map
(knobs → measured D, lacunarity, tortuosity, cable length); pick 4 settings spanning measured
D ≈ 1.2–1.65 at roughly matched cable length (`fractal_calibration.py`).

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

Open-loop growth to fixed path-length caps {250, 400} px with the F1 settings, N = 100, 3 seeds;
freeze; rescale to Λ ∈ {0.9, 0.97, 1.0} of the critical gain; measure dynamic range, avalanche
statistics, modularity, λ_c. Isolates "does geometry matter at the same distance from criticality".

## Analysis / figures
* measured D vs knob (F1 calibration heat-maps)
* D_neuron, D_network, lacunarity vs: s.d. σ, m, Λ, dynamic range, modularity, λ_c (F2, F3)
* partial correlations controlling for total length / span

Draft specs: `sweep_specs/drafts/f1_morph_calibration.json` (ready), `f2a_fractal_fss.json`,
`f2b_fractal_slow.json`, `f3_matched_gain.json` (need the F1 settings in place of the "FROM_F1"
placeholders; f2a needs its two N entries expanded to 4 morphology settings each).
