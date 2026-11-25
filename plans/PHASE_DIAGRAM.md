# Plan: critical vs bursting — when does homeostatic growth produce which? (not yet started)

**Question.** Which combinations of drive ratio h/r* and synaptic-depression strength make
homeostatic structural growth self-organise to (near-)criticality, and which to periodic
system-wide bursting? Does neuron geometry (fractal dimension, tortuosity, lacunarity) move the
boundary, and does it act only as a feedback gain?

**Why this framing.** Rate homeostasis alone sets m ≈ 1 − h/r* in a continuous-transition model
(Zierenberg et al. 2018; Kossio et al. 2018; reproduced in s1). The Stiber-lab disc-growth culture
model (E/I LIF + short-term depression; Zenodo 17220617) instead meets its rate set point through
stereotyped, network-wide bursts (all 10 000 neurons, CV of burst size 1–5 %, 0.93 Hz vs 1.0 Hz set
point; `reports/stiber_tR_1.0--fE_0.98.png`). Nobody has mapped where between these the system lands
as a function of drive and depression, and how arbor geometry and spatial wiring shift it. See
`reports/Growing critical networks literature.md` for the novelty assessment.

## Step 1 — directed axo-dendritic connectivity (separate generation file + mode)

Current connectivity (`run_sim._compute_connectivity_from_dendrites`) is dendrite–dendrite contact,
symmetrised (M = Mᵀ), with no axons. That is the main realism gap and likely inflates reciprocal
loops (reverberation, flattened avalanche bodies, σ vs m gap).

* **New module `axon_nets.py`** — leaves `branch_nets.py` and the dendritic mode untouched.
  * One axon per neuron, grown from the soma (own seed/direction), with its own growth rules:
    long, weakly branched (few collaterals, larger inter-branch length, lower branching order),
    low tortuosity, optional terminal arborisation, thin (1 px).
  * Own per-neuron raster layers (axon skeleton/filled), reusing the incremental-count machinery
    pattern but stored **per neuron in local bounding boxes** (see step 2).
  * Axons do not stop at contacts (pass-through ≈ 1); optional self-avoidance only within an axon.
  * Growth control (config): `axon_growth = "open_loop" | "homeostatic" | "coupled"`
    — open loop to a fixed length; homeostatic with its own rate rule (MSP-style: axonal and
    dendritic elements respond to activity, Butz & van Ooyen); coupled = scales with the dendritic
    target.
* **Mode switch:** `morph.connectivity_mode = "dendro_dendritic"` (default, legacy, bit-identical)
  | `"axo_dendritic"`.
  * axo_dendritic: contact_ij = pixels of neuron j's **axon** within `contact_radius_px` of neuron i's
    **dendrites/soma** → M_ij = contact_ij / contact_ref_px (directed, no symmetrisation; synapse
    number ∝ axon–dendrite overlap, i.e. Peters' rule).
  * The dendrite fractal geometry then directly sets each neuron's capture cross-section.
* Rendering: axons drawn in a distinct colour; GIF/figure support.
* Tests:
  * the legacy mode is unchanged (hash of M on a fixed seed);
  * M is asymmetric in axo mode;
  * contacts are invariant to the order in which neurons are processed;
  * smoke sweep in both modes.

## Step 2 — memory: local per-neuron layers

Per-neuron full-canvas uint16 layers cap N at ~200–400 (N = 400 ≈ 16 GB). Store each neuron's
dendrite (and axon) count layers in a growing local bounding box (re-allocate on expansion), keep
global label masks only. Target: N = 1 000–2 000 on this machine. Needed for credible finite-size
claims and for the Das & Levina 3/2 vs 5/4 question.

## Step 3 — synaptic depression (+ optional inhibition) in the branching model

* Short-term depression (Levina, Herrmann & Geisel 2007 style) as a per-synapse or per-presynaptic
  resource x_j ∈ [0, 1]:
  * each spike of j uses a fraction U of the resource;
  * the resource recovers with τ_rec;
  * the effective transmission is p_ij = min(1, spike_input · M_ij · x_j / x_ref).
  * Config: `lif.depression_U`, `lif.depression_tau_rec` (U = 0 → off, current model).
* Optional inhibitory fraction `morph.inhibitory_fraction` (Dale's principle). In the branching
  picture, an inhibitory spike lowers its targets' firing probability for one step via a factor
  (1 − q_ij). Default 0.
* Lives in `lif_dynamics.py` (shared by plastic runs and probes); legacy behaviour unchanged when
  off.

## Step 4 — re-validate the basics in the new modes

* m vs 1 − h/r* (and Zierenberg's correction).
* Λ (largest eigenvalue of P), causal m, MR m; explain the gaps (coalescence, refractoriness).
* **Timescale-separation convergence:** k = morphology_update_interval and growth rate ×{1, 0.3,
  0.1, 0.03}. The amplitude of the σ swings and the exponents must converge before the
  quasi-criticality claims stand.
* Fixed-window exponents, σ-conditioned avalanche test, dynamic range — existing tools.

## Step 5 — phase diagram

* Axes:
  * h/r* ∈ {0.002, 0.005, 0.01, 0.02, 0.05};
  * depression (U, τ_rec) along a 1-D strength line, e.g. U ∈ {0, 0.05, 0.1, 0.2, 0.4} at fixed τ_rec;
  * optional inhibitory fraction {0, 0.2}.
* N = 200 first, then N = 800 at points near the boundary; 3 seeds; axo-dendritic mode.
* Order parameters / classifiers per run:
  * burst index: fraction of spikes in system-spanning events (S > 0.5 N), CV and bimodality of
    P(S), participation, IBI regularity (CV of intervals);
  * criticality: m (MR, causal), Λ, window exponents, crackling relation, σ-conditioned ratio;
  * homeostasis: rate / target, s.d. σ(t);
  * function: dynamic range (coupling scan at the set point).
* Outputs: phase map (critical / reverberating / quasi-critical / bursting), boundary location vs
  h/r* and U, dynamic range across the map.
* Expectations:
  * low h/r* plus depression → bursting (Zierenberg low-input regime; Stiber);
  * moderate h/r*, no depression → reverberating, m ≈ 1 − h/r*;
  * depression alone can produce SOqC (Levina) — the map separates the growth and depression
    contributions.

## Step 6 — geometry as a third axis (ties in plans/FRACTAL_SWEEPS.md)

* Run F1 (knobs → measured D, tortuosity, lacunarity) in axo-dendritic mode.
* Pick 3 geometries (measured D_neuron ≈ 1.2 / 1.4 / 1.6 at matched cable).
* Repeat the phase-diagram points **near the boundary** for each geometry.
* **Gain-collapse test:**
  * measure each geometry's coupling-per-growth gain g_geo = dσ/dL (the derivative of the coupling
    trace against the mean path length near the set point; it scales with the arbor mass in the
    growth annulus, ∝ r^(D−1));
  * rescale the growth rate by g_geo and check whether the boundaries collapse.
  * Collapse → geometry acts only as the homeostatic loop gain (universality statement).
  * No collapse → spatial structure matters, e.g. burst nucleation tied to connection-distance decay
    and lacunarity (Orlandi et al. 2013 noise focusing; Stiber burst origins).
* Also: avalanche spatial extent and burst origins vs the lacunarity map.

## Step 7 (later, optional) — realism check and data

* Minimal LIF layer (E/I, Tsodyks–Markram depression, delays, ms units) on the same grown
  axo-dendritic networks. Does the phase-diagram structure survive?
* Compare developmental trajectories with Wagenaar et al. 2006 (Internet Archive mirror) and the
  remaining three Stiber runs (`external_analysis/stiber_bursts.py`).

## Deliverable / paper framing

"How homeostatic growth chooses between criticality and bursting: drive, depression and dendritic
geometry". Core results: the phase map (steps 5–6), the gain-collapse test, timescale-separation
convergence (step 4), and the dynamic-range landscape across the map.

## Rough effort

| Step | Effort |
|---|---|
| 1 (axon module + mode) | 1–2 days |
| 2 (local layers) | 1–2 days |
| 3 (depression / inhibition) | ½ day |
| 4 | 1 overnight |
| 5 | 1–2 overnights |
| 6 | F1 (~10 min) + 1 overnight |
