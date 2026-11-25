import random
import os
import json
import time
from dataclasses import dataclass


import matplotlib
matplotlib.use('Agg')  # Force non-interactive backend before importing pyplot

import cv2
import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import find_objects

from tqdm import tqdm

import avalanches
import branch_nets
from sim import Neuron as LIFNeuron


@dataclass
class MorphologyConfig:
    width: int = 900
    height: int = 900
    num_neurons: int = 24
    edge_margin: int = 60
    max_growth_steps: int | None = None
    max_growth_steps_safety_cap: int | None = 2000
    growth_dt: float | None = None  # if set, must match LIF dt
    deterministic_growth: bool = False
    morphology_seed: int = 7
    contact_radius_px: int = 6
    connectivity_gain: float = 8.0
    # 'absolute': w_ij = gain * contact_px / contact_ref_px, so more arbor overlap -> more input
    #             (Kossio-style; growth actually changes a neuron's drive).
    # 'row':      legacy row-normalisation (every connected neuron gets the same total input).
    connectivity_normalization: str = "absolute"
    contact_ref_px: float = 50.0
    prune_forget_on_degrowth: bool = True
    # Neurons below target whose tips have all been stopped (collision / edge) sprout a new
    # primary dendrite from the soma in the least crowded direction.
    resprout_blocked: bool = True
    resprout_cooldown_steps: int = 40        # base wait between sprouts; doubles per failed sprout
    max_resprouts_per_neuron: int | None = 4  # max re-sprouted primaries alive at once per neuron
    # Any primary dendrite (original or re-sprout) blocked before reaching this length without
    # having branched is retracted -- a collapsed growth cone -- so somas don't carry stubs.
    retract_blocked_stubs: bool = True
    stub_min_length_px: float = 30.0
    stub_check_interval_steps: int = 20   # how often blocked stubs are looked for
    # Prune only once the realised arbor exceeds the target by this much. Retractions below
    # one pixel are invisible on the raster masks, and pruning every step for ~0.06 px of
    # retraction dominated run time. 0 = prune every step (legacy).
    prune_quantum_px: float = 1.0


@dataclass
class LIFConfig:
    duration: float = 250.0
    dt: float = 0.1
    spike_input: float = 3.0
    input_noise: float = 0.2
    # External drive. 'uniform' (legacy): dV += input_noise * U(0,1) * dt every step, a positive
    # mean drive that makes neurons fire regularly on their own. 'poisson': sparse external
    # kicks of size ext_kick at rate ext_rate per neuron (Kossio et al.-style), which separates
    # spontaneous activity (h ~ ext_rate) from recurrent propagation; for a driven branching
    # process homeostasis then sets m ~ 1 - h / tfr.
    external_drive: str = "uniform"
    ext_rate: float = 0.02
    ext_kick: float = 2.0
    leak: float = 0.02
    spike_threshold: float = 1.0
    rest: float = -1.0
    tfr: float = 0.03
    growth_coeff: float = 0.03
    # Homeostatic rate estimate: exponential trace with time constant rate_tau (same time
    # units as dt). Replaces the lifetime average nspikes/t, which stops responding late in a run.
    rate_tau: float = 50.0
    # 'proportional': dL/dt = length_growth_rate * (1 - r/tfr)   (Kossio et al.-style, smooth)
    # 'multiplicative': legacy bang-bang L *= (1 +/- growth_coeff) every step
    growth_rule: str = "proportional"
    length_growth_rate: float = 0.5   # px per time unit at r = 0
    # Target length may lead the physically realised length by at most this many px, so a
    # blocked/finished arbor cannot accumulate an unbounded target. None disables.
    target_lead_px: float | None = 40.0
    initial_target_length_scale: float = 1.5
    # Pin each neuron's branching schedule (segment length per branch order, taper) to its
    # intrinsic arbor size (sampled from [min_total_length, max_total_length]). If False the
    # schedule follows the moving homeostatic target, which stunts arbors when the target is
    # small (all branch orders complete within a few tens of px).
    fixed_branch_schedule: bool = True
    plasticity_min_total_length: float = 0.0
    plasticity_max_total_length: float | None = 1500.0  # hard cap on target length
    morphology_update_interval: int = 15
    connectivity_update_interval: int = 1
    gif_frame_stride: int = 1
    # Morphological baseline range used to sample each neuron's initial target length.
    min_total_length: float = 80.0
    max_total_length: float = 420.0
    animate: bool = True
    gif_path: str = "branch_plasticity.gif"
    gif_fps: float = 30.0


@dataclass
class BranchNeuronConfig:
    mean_soma_radius: float = 14.0
    std_soma_radius: float = 2.0
    depth: int = 4
    D: float = 1.8
    mean_branches: int = 2
    branch_angle: float = np.pi / 4.5
    weave_type: str = "Gauss"
    randomness: float = 0.45
    curviness: str = "Gauss"
    curviness_magnitude: float = 1.2
    n_primary_dendrites: int = 4
    thickness_soma_ratio: float = 0.28
    growth_rate: float = 2.8
    variable_growth_rate: bool = False
    variable_primary_dendrites: bool = False
    pass_through_probability: float = 0.10
    fusion_probability: float = 0.05
    fusion_distance_threshold: float = 12.0
    use_splines: bool = True
    # Stochastic branching: each tip branches as a Poisson process along its path with mean
    # inter-branch length branch_length_scale * r**order, r = mean_branches**(-1/D) (the same
    # fractal length ratio as the deterministic schedule). False = legacy depth schedule.
    stochastic_branching: bool = True
    branch_length_scale: float = 60.0
    branch_min_segment_px: float = 8.0
    trajectory_planning: bool = True
    trajectory_length: float = 40.0
    trajectory_deviation: float = 10.0
    trajectory_replan_interval: float = 15.0
    trajectory_smoothness: float = 0.2
    spline_segment_length: float = 2.0
    self_avoidance: bool = False
    self_avoidance_radius: int = 4
    self_avoidance_retries: int = 3
    self_avoidance_parent_exempt_radius: int = 4
    self_avoidance_min_path_distance: float = 12.0
    self_avoidance_cone_half_angle: float = np.pi / 3


@dataclass
class AnalysisConfig:
    show_summary_plots: bool = True
    save_summary_plots: bool = True
    summary_plots_dir: str = "summary_plots"
    # Thresholded-avalanche analysis (see avalanches.py): bins of avalanche_bin_steps steps,
    # threshold = this percentile of population activity, first fraction of run discarded.
    avalanche_bin_steps: int = 1
    avalanche_threshold_percentile: float = 50.0
    avalanche_burn_in_fraction: float = 0.3
    save_final_morphology_png: bool = True
    final_morphology_png_path: str = "final_branch_network.png"
    final_render_mask_type: str = "filled"
    final_render_colormap: str = "summer"
    final_render_show_connections: bool = True
    final_render_show_fusions: bool = True
    final_render_show_soma_connections: bool = True
    final_render_show_soma_labels: bool = False
    final_render_background_color: str = "black"
    enable_profiling: bool = False
    profiling_output_dir: str = "profiling"
    profiling_summary_plot_path: str = "profiling/profiling_summary.png"
    profiling_summary_json_path: str = "profiling/profiling_summary.json"


class BranchMorphPlasticLIFSim:
    """
    Online merge of:
    - LIF dynamics / homeostatic growth rule from sim.py
    - Incremental dendrite growth from branch_nets.py

    Dendrites grow during LIF simulation steps (same dt), and degrowth is handled by
    in-place pruning/replay from recorded segment history.
    """

    def __init__(self, morph_cfg: MorphologyConfig, lif_cfg: LIFConfig, branch_cfg: BranchNeuronConfig):
        self.morph_cfg = morph_cfg
        self.lif_cfg = lif_cfg
        self.branch_cfg = branch_cfg
        self.N = morph_cfg.num_neurons

        # By design, branch growth integration uses the same dt as the LIF simulation.
        if self.morph_cfg.growth_dt is None:
            self.morph_cfg.growth_dt = self.lif_cfg.dt
        elif abs(float(self.morph_cfg.growth_dt) - float(self.lif_cfg.dt)) > 1e-12:
            raise ValueError(
                f"MorphologyConfig.growth_dt ({self.morph_cfg.growth_dt}) must match "
                f"LIFConfig.dt ({self.lif_cfg.dt}). Set growth_dt=None to auto-sync."
            )

        if self.lif_cfg.morphology_update_interval != 1:
            raise ValueError(
                "LIFConfig.morphology_update_interval must be 1 for online branch growth. "
                "Growth is integrated at every simulation step by design."
            )

        self.seed_coordinates = self._sample_seed_coordinates()
        self.initial_target_lengths = self._sample_initial_target_lengths()
        self.total_lengths = self.initial_target_lengths.copy()  # desired lengths
        self.visible_lengths = np.zeros(self.N, dtype=float)  # realized lengths from recorded growth
        self.branch_net = None
        self.M = np.zeros((self.N, self.N), dtype=float)
        self.t = 0.0
        self.tfr = lif_cfg.tfr  # compatibility with sim.Graphics helpers
        self.rates = np.zeros(self.N, dtype=float)  # exponential firing-rate trace
        self._last_resprout_step = {}
        self._resprout_failures = {}
        if self.lif_cfg.plasticity_max_total_length is None:
            print("warning: plasticity_max_total_length=None -> target lengths are uncapped")

        self.lif_neurons = [
            LIFNeuron(
                potential=2 * np.random.rand() - 1,
                position=(self.seed_coordinates[i][0], self.seed_coordinates[i][1]),
                size=float(self.total_lengths[i]),
                reset_potential=lif_cfg.rest,
            )
            for i in range(self.N)
        ]

        self._init_online_morphology()

    def _sample_seed_coordinates(self):
        # Use branch_nets placement logic once so positions remain stable across rebuilds.
        tmp_network_params = {
            "width": self.morph_cfg.width,
            "height": self.morph_cfg.height,
            "num_neurons": self.morph_cfg.num_neurons,
            "edge_margin": self.morph_cfg.edge_margin,
            "deterministic": True,
            "manual_buffer_distance": 15,
            "jitter_factor": 0.5,
        }
        neuron_params = self._base_neuron_params(total_length=self.lif_cfg.min_total_length)
        sim_props = {"dt": 1.0, "max_growth_steps": 1}

        random.seed(self.morph_cfg.morphology_seed)
        np.random.seed(self.morph_cfg.morphology_seed)
        net = branch_nets.Network(
            width=tmp_network_params["width"],
            height=tmp_network_params["height"],
            num_neurons=tmp_network_params["num_neurons"],
            neuron_params=neuron_params,
            edge_margin=tmp_network_params["edge_margin"],
            sim_props=sim_props,
            deterministic=True,
            seed_coordinates=None,
            network_params=tmp_network_params,
        )
        coords = []
        for n in net.neurons:
            coords.append((float(n.position[0]), float(n.position[1])))
        return coords

    def _sample_initial_target_lengths(self):
        lo = float(min(self.lif_cfg.min_total_length, self.lif_cfg.max_total_length))
        hi = float(max(self.lif_cfg.min_total_length, self.lif_cfg.max_total_length))
        rng = np.random.default_rng(int(self.morph_cfg.morphology_seed) + 101)
        base = rng.uniform(lo, hi, size=self.N) if hi > lo else np.full(self.N, lo, dtype=float)
        targets = np.asarray(base, dtype=float) * float(self.lif_cfg.initial_target_length_scale)
        dyn_max = self.lif_cfg.plasticity_max_total_length
        if dyn_max is not None:
            targets = np.minimum(targets, float(dyn_max))
        targets = np.maximum(targets, float(self.lif_cfg.plasticity_min_total_length))
        return targets.astype(float)

    def _base_neuron_params(self, total_length: float):
        cfg = self.branch_cfg
        return {
            "mean_soma_radius": cfg.mean_soma_radius,
            "std_soma_radius": cfg.std_soma_radius,
            "depth": cfg.depth,
            "D": cfg.D,
            "mean_branches": cfg.mean_branches,
            "branch_angle": cfg.branch_angle,
            "weave_type": cfg.weave_type,
            "randomness": cfg.randomness,
            "curviness": cfg.curviness,
            "curviness_magnitude": cfg.curviness_magnitude,
            "n_primary_dendrites": cfg.n_primary_dendrites,
            "total_length": float(total_length),
            "thickness_soma_ratio": cfg.thickness_soma_ratio,
            "growth_rate": cfg.growth_rate,
            "variable_growth_rate": cfg.variable_growth_rate,
            "variable_primary_dendrites": cfg.variable_primary_dendrites,
            "pass_through_probability": cfg.pass_through_probability,
            "fusion_probability": cfg.fusion_probability,
            "fusion_distance_threshold": cfg.fusion_distance_threshold,
            "use_splines": cfg.use_splines,
            "stochastic_branching": cfg.stochastic_branching,
            "branch_length_scale": cfg.branch_length_scale,
            "branch_min_segment_px": cfg.branch_min_segment_px,
            "trajectory_planning": cfg.trajectory_planning,
            "trajectory_length": cfg.trajectory_length,
            "trajectory_deviation": cfg.trajectory_deviation,
            "trajectory_replan_interval": cfg.trajectory_replan_interval,
            "trajectory_smoothness": cfg.trajectory_smoothness,
            "spline_segment_length": cfg.spline_segment_length,
            "self_avoidance": cfg.self_avoidance,
            "self_avoidance_radius": cfg.self_avoidance_radius,
            "self_avoidance_retries": cfg.self_avoidance_retries,
            "self_avoidance_parent_exempt_radius": cfg.self_avoidance_parent_exempt_radius,
            "self_avoidance_min_path_distance": cfg.self_avoidance_min_path_distance,
            "self_avoidance_cone_half_angle": cfg.self_avoidance_cone_half_angle,
        }

    def _set_neuron_total_length(self, dendrite, new_total_length: float):
        dendrite.total_length = max(1.0, float(new_total_length))
        if dendrite.depth <= 0:
            dendrite.branch_lengths = np.array([], dtype=float)
            return
        safe_br = max(float(dendrite.mean_branches), 1.0)
        safe_D = dendrite.D if dendrite.D > 0 else 1.0
        r_val = safe_br ** (-1.0 / safe_D)
        terms = np.array([r_val ** i for i in range(dendrite.depth)], dtype=float)
        denom = float(np.sum(terms)) if np.sum(terms) > 0 else 1.0
        dendrite.branch_lengths = (dendrite.total_length / denom) * terms

    def _init_online_morphology(self):
        # Initialize branch network with somas + active tips only; no pregrowth.
        random.seed(self.morph_cfg.morphology_seed)
        np.random.seed(self.morph_cfg.morphology_seed)

        neuron_params = self._base_neuron_params(total_length=float(self.lif_cfg.min_total_length))
        network_params = {
            "width": self.morph_cfg.width,
            "height": self.morph_cfg.height,
            "num_neurons": self.morph_cfg.num_neurons,
            "edge_margin": self.morph_cfg.edge_margin,
            "deterministic": self.morph_cfg.deterministic_growth,
            "manual_buffer_distance": 35,
            "jitter_factor": 0.0,
        }
        sim_props = {
            "dt": self.lif_cfg.dt,
            "max_growth_steps": self.morph_cfg.max_growth_steps,
            "max_growth_steps_safety_cap": self.morph_cfg.max_growth_steps_safety_cap,
            "record_history_every_n_steps": 0,
        }

        net = branch_nets.Network(
            width=network_params["width"],
            height=network_params["height"],
            num_neurons=network_params["num_neurons"],
            neuron_params=neuron_params,
            edge_margin=network_params["edge_margin"],
            sim_props=sim_props,
            deterministic=network_params["deterministic"],
            seed_coordinates=self.seed_coordinates,
            network_params=network_params,
        )

        for i, n in enumerate(net.neurons):
            if self.lif_cfg.fixed_branch_schedule:
                n.fixed_schedule_length = float(self.initial_target_lengths[i])
                self._set_neuron_total_length(n.dendrite, n.fixed_schedule_length)
            else:
                # Set each neuron's current target length at simulation start.
                self._set_neuron_total_length(n.dendrite, self.total_lengths[i])
            n.visible_total_length = float(self.total_lengths[i])
            # Keep LIF display positions aligned to branch somas.
            self.lif_neurons[i].position = (float(n.position[0]), float(n.position[1]))
            self.lif_neurons[i].size = float(self.total_lengths[i])

        self.branch_net = net
        # Replay current (empty) branch history so masks reflect the current visible lengths
        # and finalize auxiliary masks for rendering/connectivity.
        self.branch_net.apply_visible_lengths(
            self.total_lengths,
            finalize_masks=False,
            clear_interaction_annotations=False
        )
        self.visible_lengths = np.minimum(self.total_lengths, self.branch_net.get_recorded_branch_length_limits())
        for i, neuron in enumerate(self.lif_neurons):
            neuron.size = float(self.visible_lengths[i])
        self.M = self._compute_connectivity_from_dendrites()

    def _step_morphology_online(self, current_step: int, allow_growth: bool=True, update_connectivity: bool=True,
                                degrow_mask=None, grow_mask=None):
        if self.branch_net is None:
            return {}
        prof = {
            "prune_forget_s": 0.0,
            "replay_pre_s": 0.0,
            "growth_s": 0.0,
            "growth_direction_select_s": 0.0,
            "growth_self_avoidance_s": 0.0,
            "growth_collision_scan_s": 0.0,
            "growth_draw_record_s": 0.0,
            "replay_post_s": 0.0,
            "connectivity_s": 0.0,
            "self_avoidance_attempts": 0,
            "self_avoidance_rejections": 0,
            "tips_processed": 0,
            "active_tips_after": 0,
            "segments_replayed": 0,
            "segments_clipped": 0,
            "segments_skipped": 0,
        }

        did_prune_forget = False
        incremental_masks = bool(getattr(self.branch_net, "incremental_branch_mode", False))
        if self.morph_cfg.prune_forget_on_degrowth and degrow_mask is not None and np.any(degrow_mask):
            t_pf = time.perf_counter()
            pruned_lengths = self.total_lengths.copy()
            # Only force prune-and-forget for neurons currently in degrowth; others keep history.
            for i, should_prune in enumerate(np.asarray(degrow_mask, dtype=bool)):
                if not should_prune:
                    pruned_lengths[i] = np.inf
            safe_lengths = np.minimum(pruned_lengths, self.branch_net.get_recorded_branch_length_limits())
            self.branch_net.prune_and_forget_lengths(
                safe_lengths,
                finalize_masks=False,
                clear_interaction_annotations=False,
            )
            prof["prune_forget_s"] += time.perf_counter() - t_pf
            did_prune_forget = True

        # Update each branch neuron's target schedule and visible cutoff.
        for i, n in enumerate(self.branch_net.neurons):
            n.set_visible_total_length(float(self.total_lengths[i]))
            # Per-tip growth limit: target length acts as a maximum path distance from the soma.
            n.path_length_cap = float(self.total_lengths[i])

        # Rebuild masks first so collision checks see the current pruned morphology.
        if not incremental_masks and not did_prune_forget:
            t_r1 = time.perf_counter()
            self.branch_net.rebuild_masks_from_branch_history(
                finalize_masks=False,
                clear_interaction_annotations=False
            )
            prof["replay_pre_s"] += time.perf_counter() - t_r1

        # Online growth: only advance neurons that have not yet recorded enough path length.
        if allow_growth:
            t_g = time.perf_counter()
            for i, n in enumerate(self.branch_net.neurons):
                requested_len = float(self.total_lengths[i])
                recorded_len = float(getattr(n, "max_recorded_path_length", 0.0))
                wants_growth = requested_len > recorded_len + 1e-6
                if wants_growth and not n.is_growing:
                    # Neuron was pruned back to its soma (or lost all tips): re-sprout primaries.
                    n.restore_missing_primary_tips()
                if wants_growth and not n.is_growing and self.morph_cfg.resprout_blocked:
                    # Still no tips: every branch was stopped by a collision/edge. Sprout a
                    # new primary dendrite (rate-limited, backing off after failed sprouts)
                    # so the neuron can reach its target.
                    cap = self.morph_cfg.max_resprouts_per_neuron
                    live = len(getattr(n, "_resprout_ids", ()))
                    last = self._last_resprout_step.get(i, -10**9)
                    wait = self.morph_cfg.resprout_cooldown_steps * 2 ** min(self._resprout_failures.get(i, 0), 6)
                    if current_step - last >= wait and (cap is None or live < cap):
                        if n.sprout_new_primary():
                            self._last_resprout_step[i] = current_step
                            prof["resprouts"] = prof.get("resprouts", 0) + 1
                # Grow while the target is not shrinking (rate <= tfr). Each tip extends up to
                # the cap independently, so young/regrowing branches catch up with old ones.
                in_growth_phase = True if grow_mask is None else bool(grow_mask[i])
                if n.is_growing and in_growth_phase and n.has_growable_tip(requested_len):
                    n.step(self.lif_cfg.dt, current_step=current_step)
                    st = getattr(n, "last_step_stats", {})
                    prof["self_avoidance_attempts"] += int(st.get("self_avoidance_attempts", 0))
                    prof["self_avoidance_rejections"] += int(st.get("self_avoidance_rejections", 0))
                    prof["tips_processed"] += int(st.get("tips_processed", 0))
                    prof["growth_direction_select_s"] += float(st.get("growth_direction_select_s", 0.0))
                    prof["growth_self_avoidance_s"] += float(st.get("growth_self_avoidance_s", 0.0))
                    prof["growth_collision_scan_s"] += float(st.get("growth_collision_scan_s", 0.0))
                    prof["growth_draw_record_s"] += float(st.get("growth_draw_record_s", 0.0))
            prof["growth_s"] += time.perf_counter() - t_g

        # Retract re-sprouts that were blocked while still stubs.
        did_retract = False
        check_stubs = (not allow_growth) or (current_step % max(1, self.morph_cfg.stub_check_interval_steps) == 0)
        if check_stubs and (self.morph_cfg.resprout_blocked or self.morph_cfg.retract_blocked_stubs):
            for i, n in enumerate(self.branch_net.neurons):
                if self.morph_cfg.retract_blocked_stubs:
                    k = n.retract_blocked_primaries(self.morph_cfg.stub_min_length_px)
                else:
                    k = n.retract_blocked_resprouts(self.morph_cfg.stub_min_length_px)
                if k:
                    did_retract = True
                    self._resprout_failures[i] = self._resprout_failures.get(i, 0) + k

        # Forget interaction markers whose branches were pruned/retracted away. Markers are only
        # used for rendering (not the dynamics), so this runs once at the final sync.
        if not allow_growth:
            self.branch_net.prune_stale_interactions()

        # Replay after growth to apply pruning/visibility cutoffs exactly, then finalize masks.
        if not incremental_masks:
            t_r2 = time.perf_counter()
            self.branch_net.rebuild_masks_from_branch_history(
                finalize_masks=False,
                clear_interaction_annotations=False
            )
            prof["replay_post_s"] += time.perf_counter() - t_r2

        # The actual visible length cannot exceed recorded growth.
        recorded_limits = self.branch_net.get_recorded_branch_length_limits()
        self.visible_lengths = np.minimum(self.total_lengths, recorded_limits)
        for i, neuron in enumerate(self.lif_neurons):
            neuron.size = float(self.visible_lengths[i])
        if update_connectivity:
            t_c = time.perf_counter()
            self.M = self._compute_connectivity_from_dendrites()
            prof["connectivity_s"] += time.perf_counter() - t_c
        prof["active_tips_after"] = int(sum(len(n.active_tips) for n in self.branch_net.neurons))
        rstats = getattr(self.branch_net, "last_rebuild_history_stats", {})
        prof["segments_replayed"] = int(rstats.get("segments_replayed", 0))
        prof["segments_clipped"] = int(rstats.get("segments_clipped", 0))
        prof["segments_skipped"] = int(rstats.get("segments_skipped", 0))
        return prof

    def _compute_connectivity_from_dendrites(self):
        """
        Contact-based connectivity from the labelled dendrite mask.

        contact_ij = #pixels of neuron j inside neuron i's arbor dilated by contact_radius_px.
        Computed per neuron inside its own (padded) bounding box with a single bincount over
        the labels found there, instead of N^2 full-canvas mask ANDs. Equivalent to the
        previous implementation.
        """
        labels = self.branch_net.filled_mask
        n = self.N
        radius = max(1, int(self.morph_cfg.contact_radius_px))
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * radius + 1, 2 * radius + 1))
        H, W = labels.shape

        areas = np.bincount(labels.ravel(), minlength=n + 1)[1:n + 1].astype(float)
        boxes = find_objects(labels, max_label=n)
        C = np.zeros((n + 1, n + 1), dtype=float)   # C[i, j] = contact_ij (1-based labels)
        for i, box in enumerate(boxes):
            if box is None:
                continue
            ys, xs = box
            y0, y1 = max(0, ys.start - radius), min(H, ys.stop + radius)
            x0, x1 = max(0, xs.start - radius), min(W, xs.stop + radius)
            local = labels[y0:y1, x0:x1]
            dil = cv2.dilate((local == i + 1).view(np.uint8), kernel, iterations=1).view(bool)
            C[i + 1] = np.bincount(local[dil], minlength=n + 1)[:n + 1]

        contact = 0.5 * (C[1:, 1:] + C[1:, 1:].T)
        np.fill_diagonal(contact, 0.0)
        if self.morph_cfg.connectivity_normalization == "row":
            # Legacy: normalise by smaller arbor size, then row-normalise.
            M = contact / np.maximum(1.0, np.minimum(areas[:, None], areas[None, :]))
            row_sums = M.sum(axis=1, keepdims=True)
            nonzero = row_sums.squeeze() > 0
            M[nonzero] = M[nonzero] / row_sums[nonzero]
        else:
            M = contact / max(1e-9, self.morph_cfg.contact_ref_px)
        M *= self.morph_cfg.connectivity_gain
        return M

    def _update_dendrite_lengths(self):
        if self.t <= 0:
            return
        cfg = self.lif_cfg
        if cfg.growth_rule == "multiplicative":
            grow = self.rates < cfg.tfr
            self.total_lengths[grow] *= (1.0 + cfg.growth_coeff)
            self.total_lengths[~grow] *= (1.0 - cfg.growth_coeff)
        else:
            self.total_lengths += cfg.dt * cfg.length_growth_rate * (1.0 - self.rates / cfg.tfr)

        if cfg.target_lead_px is not None:
            # Don't let the target run away from what has physically been grown.
            self.total_lengths = np.minimum(self.total_lengths, self.visible_lengths + cfg.target_lead_px)

        self.total_lengths = np.clip(
            self.total_lengths,
            float(self.lif_cfg.plasticity_min_total_length),
            float(self.lif_cfg.plasticity_max_total_length)
            if self.lif_cfg.plasticity_max_total_length is not None
            else np.inf,
        )

    def _render_activity_frame(self, V, spikes):
        labels = self.branch_net.filled_mask
        somas = self.branch_net.somas_mask

        # Per-neuron colour lookup table, then one gather over the label image
        # (instead of N full-canvas `labels == id` comparisons).
        lut = np.zeros((self.N + 1, 3), dtype=np.uint8)
        v_norm = np.clip((np.asarray(V, dtype=float) + 1.0) / 2.0, 0.0, 1.0)
        lut[1:] = (255 * np.asarray(plt.cm.winter(v_norm))[:, :3]).astype(np.uint8)
        lut[1:][np.asarray(spikes, dtype=bool)] = (0, 255, 0)
        for idx, neuron in enumerate(self.lif_neurons):
            neuron.potential = float(V[idx])
            neuron.spiking = bool(spikes[idx])
        valid = labels <= self.N
        img = lut[np.where(valid, labels, 0)]

        # Brighten soma pixels slightly so somas remain visible.
        soma_mask = somas > 0
        img[soma_mask] = np.clip(img[soma_mask].astype(int) + 40, 0, 255).astype(np.uint8)

        return img

    def simulate(self, analysis_cfg: AnalysisConfig | None = None):
        iterations = int(self.lif_cfg.duration / self.lif_cfg.dt)
        V = np.array([n.potential for n in self.lif_neurons], dtype=float)
        spikes = np.array([bool(n.spiking) for n in self.lif_neurons], dtype=bool)

        spikes_list = []
        V_list = []
        firing_rates_list = []
        length_list = []
        growth_state_list = []
        frames = []
        profiling_enabled = bool(analysis_cfg and analysis_cfg.enable_profiling)
        profiling = {
            "phase_totals_s": {},
            "phase_calls": {},
            "per_step": [],
        } if profiling_enabled else None

        def prof_add(name, dt_s):
            if not profiling_enabled:
                return
            profiling["phase_totals_s"][name] = profiling["phase_totals_s"].get(name, 0.0) + float(dt_s)
            profiling["phase_calls"][name] = profiling["phase_calls"].get(name, 0) + 1

        for step in tqdm(range(iterations), desc="LIF + branch morphology"):
            step_prof = {"step": int(step)} if profiling_enabled else None
            t_step = time.perf_counter() if profiling_enabled else None
            self.t += self.lif_cfg.dt

            t0 = time.perf_counter() if profiling_enabled else None
            spike_input_contribution = self.lif_cfg.spike_input * np.dot(self.M, spikes.astype(float))
            if self.lif_cfg.external_drive == "poisson":
                ext = self.lif_cfg.ext_kick * np.random.poisson(self.lif_cfg.ext_rate * self.lif_cfg.dt, self.N)
                dV_dt = -self.lif_cfg.leak * (V - self.lif_cfg.rest)
                V += spike_input_contribution + ext + dV_dt * self.lif_cfg.dt
            else:
                dV_dt = self.lif_cfg.input_noise * np.random.rand(self.N) - self.lif_cfg.leak * (V - self.lif_cfg.rest)
                V += spike_input_contribution + dV_dt * self.lif_cfg.dt

            spikes = V >= self.lif_cfg.spike_threshold
            V[spikes] = self.lif_cfg.rest
            V[V < self.lif_cfg.rest] = self.lif_cfg.rest

            # Exponential firing-rate trace (spikes per time unit).
            alpha = min(1.0, self.lif_cfg.dt / max(self.lif_cfg.rate_tau, 1e-9))
            self.rates += alpha * (spikes.astype(float) / self.lif_cfg.dt - self.rates)
            for i, neuron in enumerate(self.lif_neurons):
                neuron.spiking = bool(spikes[i])
                if neuron.spiking:
                    neuron.nspikes += 1
                neuron.firing_rate = float(self.rates[i])
                neuron.potential = float(V[i])
            if profiling_enabled:
                dt_lif = time.perf_counter() - t0
                prof_add("lif_update", dt_lif)
                step_prof["lif_update_s"] = dt_lif

            prev_target_lengths = self.total_lengths.copy()
            t1 = time.perf_counter() if profiling_enabled else None
            self._update_dendrite_lengths()
            if profiling_enabled:
                dt_len = time.perf_counter() - t1
                prof_add("length_update", dt_len)
                step_prof["length_update_s"] = dt_len
            growth_state = np.where(
                np.array([n.firing_rate for n in self.lif_neurons], dtype=float) < self.lif_cfg.tfr,
                1,   # growth
                -1,  # degrowth (matches >= tfr rule)
            )
            degrow_mask = self.total_lengths < (prev_target_lengths - 1e-9)
            if self.morph_cfg.prune_quantum_px > 0:
                # Batch sub-pixel retractions: prune once the overshoot reaches the quantum.
                recorded = self.branch_net.get_recorded_branch_length_limits()
                degrow_mask = (self.total_lengths < recorded - self.morph_cfg.prune_quantum_px) & \
                              (self.total_lengths < prev_target_lengths - 1e-9)

            update_connectivity = ((step + 1) % max(1, self.lif_cfg.connectivity_update_interval) == 0)
            morph_prof = self._step_morphology_online(
                current_step=step,
                update_connectivity=update_connectivity,
                degrow_mask=degrow_mask,
                grow_mask=self.total_lengths >= prev_target_lengths - 1e-12,
            )
            if profiling_enabled:
                for k, v in morph_prof.items():
                    if k.endswith("_s"):
                        prof_add(k[:-2], v)
                    step_prof[k] = v

            spikes_list.append(spikes.copy())
            V_list.append(V.copy())
            firing_rates_list.append(np.array([n.firing_rate for n in self.lif_neurons], dtype=float))
            length_list.append(self.visible_lengths.copy())
            growth_state_list.append(growth_state.copy())

            if self.lif_cfg.animate and (step % max(1, self.lif_cfg.gif_frame_stride) == 0):
                t_r = time.perf_counter() if profiling_enabled else None
                frames.append(self._render_activity_frame(V, spikes))
                if profiling_enabled:
                    dt_render = time.perf_counter() - t_r
                    prof_add("render_frame", dt_render)
                    step_prof["render_frame_s"] = dt_render
                    step_prof["frame_rendered"] = 1
            elif profiling_enabled:
                step_prof["frame_rendered"] = 0

            if profiling_enabled:
                step_prof["step_total_s"] = time.perf_counter() - t_step
                step_prof["mean_visible_length"] = float(np.mean(self.visible_lengths))
                step_prof["mean_target_length"] = float(np.mean(self.total_lengths))
                profiling["per_step"].append(step_prof)

        # Final replay so morphology matches the last online update exactly.
        t_final_sync = time.perf_counter() if profiling_enabled else None
        self._step_morphology_online(current_step=iterations, allow_growth=False, update_connectivity=True)
        if profiling_enabled:
            prof_add("final_sync", time.perf_counter() - t_final_sync)

        if self.lif_cfg.animate and frames:
            t_gif = time.perf_counter() if profiling_enabled else None
            # imageio (pillow plugin) takes per-frame duration in milliseconds.
            imageio.mimsave(self.lif_cfg.gif_path, frames,
                            duration=1000.0 / max(1e-6, self.lif_cfg.gif_fps), loop=0)
            if profiling_enabled:
                prof_add("gif_encode", time.perf_counter() - t_gif)

        mean_fr = float(np.mean([n.firing_rate for n in self.lif_neurons]))
        print(f"mean firing rate: {mean_fr:.4f}")
        print(f"final mean target dendrite total length: {float(np.mean(self.total_lengths)):.2f}")
        print(f"final mean visible dendrite total length: {float(np.mean(self.visible_lengths)):.2f}")
        print(f"saved gif: {self.lif_cfg.gif_path}" if self.lif_cfg.animate else "gif disabled")
        if profiling_enabled:
            results_profile = self._finalize_profile_summary(profiling)
        else:
            results_profile = None

        return {
            "spikes_list": spikes_list,
            "V_list": V_list,
            "firing_rates_list": firing_rates_list,
            "length_list": length_list,
            "growth_state_list": growth_state_list,
            "branch_network": self.branch_net,
            "connectivity_matrix": self.M.copy(),
            "profiling": results_profile,
        }

    def _finalize_profile_summary(self, profiling):
        per_step = profiling.get("per_step", [])
        phase_totals = profiling.get("phase_totals_s", {})
        phase_calls = profiling.get("phase_calls", {})
        phase_avg_ms = {}
        for k, total in phase_totals.items():
            calls = max(1, int(phase_calls.get(k, 1)))
            phase_avg_ms[k] = 1000.0 * float(total) / calls
        summary = {
            "phase_totals_s": phase_totals,
            "phase_calls": phase_calls,
            "phase_avg_ms_per_call": phase_avg_ms,
            "per_step": per_step,
        }
        print("profiling summary (top phases by total time):")
        for name, total in sorted(phase_totals.items(), key=lambda kv: kv[1], reverse=True):
            avg_ms = phase_avg_ms.get(name, 0.0)
            calls = phase_calls.get(name, 0)
            print(f"  {name}: {total:.3f}s total, {avg_ms:.2f} ms/call, n={calls}")
        return summary


class BranchSimGraphics:
    """Plotting helpers for the branch-based plasticity simulation."""

    @staticmethod
    def plot_neuron_activity(spikes_list, V_list):
        plt.figure(figsize=(10, 5))
        plt.subplot(1, 2, 1)
        plt.imshow(1 - np.array(spikes_list).T, aspect="auto", cmap="binary")
        plt.xlabel("Time")
        plt.ylabel("Neuron Index")
        plt.title("Neuron Firings Over Time")

        plt.subplot(1, 2, 2)
        plt.imshow(np.array(V_list).T, aspect="auto", cmap="winter")
        plt.xlabel("Time")
        plt.ylabel("Neuron Index")
        plt.title("Membrane Potential Over Time")
        plt.colorbar(label="Membrane Potential")
        plt.tight_layout()

    @staticmethod
    def plot_firing_rates(tfr, firing_rates_list):
        plt.figure()
        fr_img = plt.imshow(np.array(firing_rates_list).T, aspect="auto")
        plt.colorbar(fr_img)
        fr_img.set_clim(0, 1.5 * tfr if tfr > 0 else max(1.0, float(np.nanmax(fr_img.get_array()))))
        plt.title(f"Neuron firing rates, tfr = {tfr}")
        plt.xlabel("Iteration #")
        plt.ylabel("Neuron index")

    @staticmethod
    def plot_dendrite_lengths(length_list):
        plt.figure()
        plt.imshow(np.array(length_list).T, aspect="auto")
        plt.colorbar(label="Total dendrite length")
        plt.xlabel("Iteration")
        plt.ylabel("Neuron index")
        plt.title("Dendrite Total Length (branch-based plasticity)")

    @staticmethod
    def plot_growth_state(growth_state_list):
        if not growth_state_list:
            return
        plt.figure()
        img = plt.imshow(np.array(growth_state_list).T, aspect="auto", cmap="bwr", vmin=-1, vmax=1)
        cbar = plt.colorbar(img, ticks=[-1, 1])
        cbar.ax.set_yticklabels(["Degrowth", "Growth"])
        plt.xlabel("Iteration")
        plt.ylabel("Neuron index")
        plt.title("Growth vs Degrowth State")

    @staticmethod
    def plot_connectivity_matrix(connectivity_matrix):
        plt.figure()
        plt.imshow(connectivity_matrix)
        plt.colorbar(label="Connectivity weight")
        plt.title("Final Dendrite-Derived Connectivity")
        plt.xlabel("Neuron index")
        plt.ylabel("Neuron index")

    @staticmethod
    def plot_avalanches(spikes_list, dt, bin_steps=1, threshold_percentile=50.0, burn_in_fraction=0.3):
        """Thresholded avalanche analysis + criticality diagnostics. Returns the result dict."""
        res = avalanches.analyze(spikes_list, bin_steps=bin_steps, threshold_percentile=threshold_percentile,
                                 burn_in_fraction=burn_in_fraction)
        if len(res["sizes"]) == 0:
            print("No avalanches detected; skipping avalanche plots.")
            return res
        avalanches.plot(res, dt=dt, title="Thresholded avalanches (branch-based sim)")
        summ = avalanches.summary(res)
        print("avalanche analysis: "
              f"n={summ['n_avalanches']}, tau={summ['size_exponent_tau']:.2f}±{summ['size_exponent_se']:.2f} "
              f"(xmin={summ['size_xmin']}), alpha={summ['duration_exponent_alpha']:.2f}±{summ['duration_exponent_se']:.2f} "
              f"(xmin={summ['duration_xmin']}), gamma fit={summ['gamma_fit']:.2f} vs predicted {summ['gamma_predicted']:.2f}, "
              f"branching ratio m={summ['branching_ratio_m']:.3f}")
        return res


def _plot_summary(results, sim_obj, analysis_cfg: AnalysisConfig):
    existing_figs = set(plt.get_fignums())
    g = BranchSimGraphics()
    g.plot_neuron_activity(results["spikes_list"], results["V_list"])
    g.plot_firing_rates(sim_obj.tfr, results["firing_rates_list"])
    g.plot_dendrite_lengths(results["length_list"])
    g.plot_growth_state(results["growth_state_list"])
    g.plot_connectivity_matrix(results["connectivity_matrix"])
    av_res = g.plot_avalanches(
        results["spikes_list"],
        dt=sim_obj.lif_cfg.dt,
        bin_steps=analysis_cfg.avalanche_bin_steps,
        threshold_percentile=analysis_cfg.avalanche_threshold_percentile,
        burn_in_fraction=analysis_cfg.avalanche_burn_in_fraction,
    )
    new_figs = [plt.figure(num) for num in plt.get_fignums() if num not in existing_figs]

    if analysis_cfg.save_summary_plots and new_figs:
        os.makedirs(analysis_cfg.summary_plots_dir, exist_ok=True)
        filenames = [
            "neuron_activity.png",
            "firing_rates.png",
            "dendrite_lengths.png",
            "growth_vs_degrowth_state.png",
            "connectivity_matrix.png",
            "avalanche_analysis.png",
        ]
        for i, fig in enumerate(new_figs):
            name = filenames[i] if i < len(filenames) else f"plot_{i+1:02d}.png"
            fig.savefig(os.path.join(analysis_cfg.summary_plots_dir, name), dpi=150, bbox_inches="tight")
        print(f"saved summary plots to: {analysis_cfg.summary_plots_dir}")
    if analysis_cfg.save_summary_plots:
        # Keep the raw raster so the avalanche analysis can be redone without re-simulating:
        #   avalanches.analyze(np.load(".../spikes.npz")["spikes"], ...)
        os.makedirs(analysis_cfg.summary_plots_dir, exist_ok=True)
        np.savez_compressed(os.path.join(analysis_cfg.summary_plots_dir, "spikes.npz"),
                            spikes=np.asarray(results["spikes_list"], dtype=bool), dt=sim_obj.lif_cfg.dt)
    if analysis_cfg.save_summary_plots and len(av_res.get("sizes", [])):
        os.makedirs(analysis_cfg.summary_plots_dir, exist_ok=True)
        with open(os.path.join(analysis_cfg.summary_plots_dir, "avalanche_summary.json"), "w") as f:
            json.dump({k: (float(v) if isinstance(v, (np.floating, float)) else v)
                       for k, v in avalanches.summary(av_res).items()}, f, indent=2)

    plt.show()


def _save_profiling_outputs(results, analysis_cfg: AnalysisConfig):
    profiling = results.get("profiling")
    if not profiling:
        return

    os.makedirs(analysis_cfg.profiling_output_dir, exist_ok=True)

    with open(analysis_cfg.profiling_summary_json_path, "w", encoding="utf-8") as f:
        json.dump(profiling, f, indent=2)
    print(f"saved profiling json: {analysis_cfg.profiling_summary_json_path}")

    phase_totals = profiling.get("phase_totals_s", {})
    per_step = profiling.get("per_step", [])
    if not phase_totals:
        return

    # Summary plot: totals bar chart + per-step timing traces.
    fig, axes = plt.subplots(2, 1, figsize=(12, 9))

    names = [k for k, _ in sorted(phase_totals.items(), key=lambda kv: kv[1], reverse=True)]
    totals = [phase_totals[k] for k in names]
    axes[0].bar(range(len(names)), totals)
    axes[0].set_xticks(range(len(names)))
    axes[0].set_xticklabels(names, rotation=35, ha="right")
    axes[0].set_ylabel("Total Time (s)")
    axes[0].set_title("Profiling Summary: Total Time by Phase")

    if per_step:
        x = np.arange(len(per_step))
        line_keys = [
            ("step_total_s", "step_total"),
            ("lif_update_s", "lif"),
            ("length_update_s", "len_update"),
            ("growth_s", "growth"),
            ("growth_direction_select_s", "growth_dir"),
            ("growth_self_avoidance_s", "growth_selfavoid"),
            ("growth_collision_scan_s", "growth_collision"),
            ("growth_draw_record_s", "growth_draw"),
            ("replay_pre_s", "replay_pre"),
            ("replay_post_s", "replay_post"),
            ("connectivity_s", "connectivity"),
            ("render_frame_s", "render"),
        ]
        for key, label in line_keys:
            y = np.array([float(s.get(key, 0.0)) * 1000.0 for s in per_step], dtype=float)
            if np.any(y > 0):
                axes[1].plot(x, y, label=label, linewidth=1.5)
        axes[1].set_xlabel("Simulation Step")
        axes[1].set_ylabel("Time (ms)")
        axes[1].set_title("Per-Step Timing (ms)")
        axes[1].grid(True, alpha=0.3)
        axes[1].legend(loc="upper right", ncol=3)

    plt.tight_layout()
    fig.savefig(analysis_cfg.profiling_summary_plot_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved profiling plot: {analysis_cfg.profiling_summary_plot_path}")


def default_configs():
    """Centralized simulation config: update values here to control the full run.
    Returns (morph_cfg, branch_cfg, lif_cfg, analysis_cfg)."""
    morph_cfg = MorphologyConfig(
        width=1024,
        height=1024,
        num_neurons=50,
        edge_margin=20,
        max_growth_steps=None,  # grow until stabilization
        max_growth_steps_safety_cap=None,  # optional guard; set None to remove
        growth_dt=None,  # auto-sync to lif_cfg.dt
        deterministic_growth=False,
        morphology_seed=np.random.randint(1000000),
        contact_radius_px=1,
        connectivity_gain=1.0,
        connectivity_normalization="absolute",  # growth -> more contact -> more input
        contact_ref_px=50.0,  # contact pixels that give a weight of 1 (x gain)
    )

    branch_cfg = BranchNeuronConfig(
        mean_soma_radius=15.0,
        std_soma_radius=3.0,
        depth=5,
        D=1.2,
        mean_branches=1.75,
        branch_angle=np.pi / 3,
        weave_type="Gauss",
        randomness=1.25,
        curviness="Gauss",
        curviness_magnitude=1.5,
        n_primary_dendrites=4,
        thickness_soma_ratio=0.28,
        growth_rate=1,
        variable_growth_rate=False,
        variable_primary_dendrites=False,
        pass_through_probability=0.20,
        fusion_probability=0.1,
        fusion_distance_threshold=12.0,
        use_splines=True,
        trajectory_planning=True,
        trajectory_length=10.0,
        trajectory_deviation=10.0,
        trajectory_replan_interval=5.0,
        trajectory_smoothness=0.2,
        spline_segment_length=2.0,  # online growth: flush geometry frequently
        self_avoidance=True,
        self_avoidance_radius=20,
        self_avoidance_retries=3,
        self_avoidance_parent_exempt_radius=4,
        self_avoidance_min_path_distance=12.0,
        self_avoidance_cone_half_angle=np.pi / 3,
    )    
    
    lif_cfg = LIFConfig(
        duration=5000.0,
        dt=0.5,
        spike_input=2.0,
        input_noise=0.5,
        leak=0.04,
        spike_threshold=1.0,
        rest=-1.0,
        tfr=0.4,
        rate_tau=50.0,  # firing-rate trace time constant (time units)
        growth_rule="proportional",  # dL/dt = length_growth_rate * (1 - rate/tfr)
        length_growth_rate=0.5,
        target_lead_px=40.0,
        plasticity_max_total_length=1500.0,
        fixed_branch_schedule=True,
        morphology_update_interval=1,
        connectivity_update_interval=5,
        gif_frame_stride=10,
        min_total_length=300.0,
        max_total_length=1024.0,
        animate=True,
        gif_path="branch_plasticity.gif",
    )
    
    analysis_cfg = AnalysisConfig(
        show_summary_plots=True,
        save_summary_plots=True,
        summary_plots_dir="summary_plots",
        avalanche_bin_steps=1,
        avalanche_threshold_percentile=50.0,
        avalanche_burn_in_fraction=0.3,
        save_final_morphology_png=True,
        final_morphology_png_path="final_branch_network.png",
        final_render_mask_type="filled",
        final_render_colormap="summer",
        final_render_show_connections=True,
        final_render_show_fusions=True,
        final_render_show_soma_connections=True,
        final_render_show_soma_labels=False,
        final_render_background_color="black",
        enable_profiling=True,
        profiling_output_dir="profiling",
        profiling_summary_plot_path="profiling/profiling_summary.png",
        profiling_summary_json_path="profiling/profiling_summary.json",
    )

    return morph_cfg, branch_cfg, lif_cfg, analysis_cfg


def main():
    morph_cfg, branch_cfg, lif_cfg, analysis_cfg = default_configs()
    sim_obj = BranchMorphPlasticLIFSim(morph_cfg, lif_cfg, branch_cfg)
    results = sim_obj.simulate(analysis_cfg=analysis_cfg)

    # Final morphology snapshot (actual branch network masks).
    if analysis_cfg.save_final_morphology_png:
        frame = sim_obj.branch_net.render_frame(
            mask_type=analysis_cfg.final_render_mask_type,
            colormap=analysis_cfg.final_render_colormap,
            show_soma_labels=analysis_cfg.final_render_show_soma_labels,
            show_connections=analysis_cfg.final_render_show_connections,
            show_fusions=analysis_cfg.final_render_show_fusions,
            show_soma_connections=analysis_cfg.final_render_show_soma_connections,
            background_color=analysis_cfg.final_render_background_color,
        )
        imageio.imwrite(analysis_cfg.final_morphology_png_path, frame)
        print(f"saved final morphology image: {analysis_cfg.final_morphology_png_path}")

    if analysis_cfg.enable_profiling:
        _save_profiling_outputs(results, analysis_cfg)
    if analysis_cfg.show_summary_plots:
        _plot_summary(results, sim_obj, analysis_cfg)


if __name__ == "__main__":
    main()
