"""
Axons and directed axo-dendritic connectivity (plans/PHASE_DIAGRAM.md, step 1).

Used only when MorphologyConfig.connectivity_mode == "axo_dendritic"; the legacy
dendro-dendritic mode (branch_nets + run_sim._compute_connectivity_from_dendrites) is untouched.

Model
-----
* One axon per neuron, initiated at the soma boundary in a random direction. A tip advances in
  steps of `step_px` with a persistent heading (Gaussian turning, s.d. `persistence` rad per step),
  so axons are long and only weakly tortuous.
* Collaterals: each tip branches as a Poisson process along its path, with mean inter-branch path
  length `collateral_length_scale`. The collateral leaves at ±`collateral_angle`, up to branch order
  `max_order`.
* Extent: every tip obeys one per-neuron path-length cap L (distance from the soma along the axon),
  as the dendrites do. Growth control (`growth`):
    - 'open_loop'   : L grows at `growth_rate` up to `length_px` (no activity feedback);
    - 'homeostatic' : dL/dt = growth_rate * (1 - r / tfr), clipped to [min_length_px, max_length_px];
    - 'coupled'     : L = length_ratio * (dendritic target length), so axons follow the dendrites.
  When L shrinks, branches are truncated at L and collaterals starting beyond L are removed.
  Regrowth then takes fresh random paths ("prune and forget", as for dendrites).
* Axons pass through everything (no collision or stopping) and are 1 px wide.
* Tips that leave the canvas stop.
* Connectivity (Peters' rule): contact_ij = number of axon pixels of neuron j that lie within
  `contact_radius_px` of neuron i's dendrites or soma. Then M_ij = contact_ij / contact_ref_px,
  directed (j -> i, the run_sim convention: input to i is M @ spikes) and not symmetrised. The
  diagonal (an axon touching its own dendrites) is zeroed. Computed as a sparse product D @ Aᵀ of
  per-neuron pixel-membership matrices.

Axon rasters are kept per neuron as sorted global pixel indices, i.e. local, not full-canvas layers.
"""
from __future__ import annotations

from dataclasses import dataclass

import cv2
import numpy as np
from scipy import sparse
from scipy.ndimage import find_objects


@dataclass
class AxonConfig:
    growth: str = "open_loop"            # 'open_loop' | 'homeostatic' | 'coupled'
    length_px: float = 600.0             # open_loop: final path-length cap
    growth_rate: float = 2.0             # px per time unit (open_loop speed / homeostatic gain)
    length_ratio: float = 2.0            # coupled: L_axon = length_ratio * dendritic target
    initial_length_px: float = 0.0
    min_length_px: float = 0.0
    max_length_px: float = 3000.0
    step_px: float = 2.0                 # tip integration step
    persistence: float = 0.06            # s.d. of heading change per step (rad)
    collateral_length_scale: float = 250.0   # mean path length between collateral branch points
    collateral_angle: float = 0.9        # rad, ± from the parent heading
    max_order: int = 2                   # 0 = main axon only
    contact_radius_px: int | None = None     # None -> MorphologyConfig.contact_radius_px
    contact_ref_px: float | None = None      # None -> MorphologyConfig.contact_ref_px
    bbox_pad_px: int = 24                # padding of the dendrite bounding boxes (owner-map bbox)
    seed_offset: int = 7919


class _Axon:
    """One neuron's axon: branches as polylines with cumulative path distance, plus active tips."""

    def __init__(self, origin, angle, soma_radius, rng, cfg: AxonConfig, width, height):
        self.cfg, self.rng, self.W, self.H = cfg, rng, width, height
        self.cap = float(cfg.initial_length_px)
        start = np.asarray(origin, float) + soma_radius * np.array([np.cos(angle), np.sin(angle)])
        self.branches = []            # dicts: pts (list of xy), path (list of float), order, parent
        self.tips = []                # dicts: branch, angle, active, next_branch
        self._new_branch(start, path0=soma_radius, angle=angle, order=0, parent=-1)
        self.dirty = True
        self._pix = np.zeros(0, np.int64)

    def _exp(self):
        return float(self.rng.exponential(self.cfg.collateral_length_scale))

    def _new_branch(self, xy, path0, angle, order, parent):
        self.branches.append(dict(pts=[np.asarray(xy, float)], path=[float(path0)], order=order, parent=parent))
        self.tips.append(dict(branch=len(self.branches) - 1, angle=float(angle), active=True,
                              next_branch=float(path0) + self._exp()))

    # ------------------------------------------------------------------ growth
    def grow(self):
        """Advance all active tips until they reach the cap (or leave the canvas)."""
        cfg = self.cfg
        changed = False
        i = 0
        while i < len(self.tips):                       # list grows as collaterals spawn
            tip = self.tips[i]
            br = self.branches[tip["branch"]]
            while tip["active"] and br["path"][-1] + cfg.step_px <= self.cap:
                tip["angle"] += self.rng.normal(0.0, cfg.persistence)
                p = br["pts"][-1] + cfg.step_px * np.array([np.cos(tip["angle"]), np.sin(tip["angle"])])
                if not (0 <= p[0] < self.W and 0 <= p[1] < self.H):
                    tip["active"] = False
                    break
                br["pts"].append(p)
                br["path"].append(br["path"][-1] + cfg.step_px)
                changed = True
                if br["path"][-1] >= tip["next_branch"]:
                    tip["next_branch"] = br["path"][-1] + self._exp()
                    if br["order"] < cfg.max_order:
                        side = 1.0 if self.rng.random() < 0.5 else -1.0
                        self._new_branch(p, br["path"][-1], tip["angle"] + side * cfg.collateral_angle,
                                         br["order"] + 1, tip["branch"])
            i += 1
        self.dirty |= changed

    def retract(self):
        """Truncate everything beyond the cap; drop collaterals that start beyond it."""
        keep = []
        remap = {}
        truncated = set()
        for bi, br in enumerate(self.branches):
            if br["path"][0] > self.cap and bi != 0:
                continue
            path = np.asarray(br["path"])
            n = int(np.searchsorted(path, self.cap, side="right"))
            if n < len(path):
                br["pts"], br["path"] = br["pts"][:max(n, 1)], br["path"][:max(n, 1)]
                self.dirty = True
                truncated.add(len(keep))
            remap[bi] = len(keep)
            keep.append(br)
        if len(keep) != len(self.branches):
            self.dirty = True
        for br in keep:                                   # parents that were dropped -> -1
            br["parent"] = remap.get(br["parent"], -1)
        new_tips = []
        for tip in self.tips:
            if tip["branch"] in remap:
                tip["branch"] = remap[tip["branch"]]
                br = keep[tip["branch"]]
                if tip["branch"] in truncated:            # a cut tip is inside the canvas again
                    tip["active"] = True
                if tip["next_branch"] <= br["path"][-1]:
                    tip["next_branch"] = br["path"][-1] + self._exp()
                new_tips.append(tip)
        self.branches, self.tips = keep, new_tips

    def set_cap(self, cap):
        cap = float(cap)
        if cap < self.cap - 1e-9:
            self.cap = cap
            self.retract()
        else:
            self.cap = cap

    # ------------------------------------------------------------------ raster / stats
    def pixels(self):
        """Sorted unique global linear indices (y * W + x) of the axon's 1-px raster."""
        if not self.dirty:
            return self._pix
        pts_all = [np.asarray(b["pts"]) for b in self.branches if len(b["pts"]) >= 2]
        if not pts_all:
            self._pix = np.zeros(0, np.int64)
        else:
            cat = np.concatenate(pts_all)
            x0, y0 = np.floor(cat.min(0)).astype(int) - 1
            x1, y1 = np.ceil(cat.max(0)).astype(int) + 2
            x0, y0 = max(x0, 0), max(y0, 0)
            canvas = np.zeros((min(y1, self.H) - y0, min(x1, self.W) - x0), np.uint8)
            polys = [np.round(p - [x0, y0]).astype(np.int32).reshape(-1, 1, 2) for p in pts_all]
            cv2.polylines(canvas, polys, False, 1, 1)
            ys, xs = np.nonzero(canvas)
            self._pix = np.unique((ys + y0).astype(np.int64) * self.W + (xs + x0))
        self.dirty = False
        return self._pix

    def cable_length(self):
        return float(sum(b["path"][-1] - b["path"][0] for b in self.branches))

    def max_path(self):
        return float(max(b["path"][-1] for b in self.branches))


class AxonNetwork:
    """Axons for all neurons of a branch_nets.Network, plus directed contact connectivity."""

    def __init__(self, branch_net, cfg: AxonConfig, morph_cfg, seed: int = 0):
        self.net = branch_net
        self.cfg = cfg
        self.H, self.W = branch_net.filled_mask.shape
        self.N = len(branch_net.neurons)
        self.radius = int(cfg.contact_radius_px if cfg.contact_radius_px is not None else morph_cfg.contact_radius_px)
        self.ref = float(cfg.contact_ref_px if cfg.contact_ref_px is not None else morph_cfg.contact_ref_px)
        rng = np.random.default_rng(seed + cfg.seed_offset)
        self.axons = []
        for n in branch_net.neurons:
            soma = n.dendrite.soma
            self.axons.append(_Axon(soma.position[:2], rng.uniform(0, 2 * np.pi), soma.radius,
                                    np.random.default_rng(rng.integers(1 << 31)), cfg, self.W, self.H))
        self.targets = np.full(self.N, float(cfg.initial_length_px))

    # ------------------------------------------------------------------ growth control
    def update_targets(self, rates, tfr, dt, dendrite_targets=None):
        cfg = self.cfg
        if cfg.growth == "open_loop":
            self.targets = np.minimum(self.targets + dt * cfg.growth_rate, cfg.length_px)
        elif cfg.growth == "homeostatic":
            self.targets = self.targets + dt * cfg.growth_rate * (1.0 - np.asarray(rates) / tfr)
        elif cfg.growth == "coupled":
            if dendrite_targets is None:
                raise ValueError("coupled axon growth needs the dendritic targets")
            self.targets = cfg.length_ratio * np.asarray(dendrite_targets, float)
        else:
            raise ValueError(f"unknown axon growth rule {cfg.growth!r}")
        self.targets = np.clip(self.targets, cfg.min_length_px, cfg.max_length_px)

    def grow(self):
        for ax, L in zip(self.axons, self.targets):
            ax.set_cap(L)
            ax.grow()

    # ------------------------------------------------------------------ connectivity
    def _axon_matrix(self):
        rows, cols = [], []
        for j, ax in enumerate(self.axons):
            p = ax.pixels()
            rows.append(np.full(len(p), j, np.int64)); cols.append(p)
        r, c = np.concatenate(rows), np.concatenate(cols)
        return sparse.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(self.N, self.H * self.W))

    def _dendrite_matrix(self):
        """D[i, pixel] = 1 if the pixel lies within `radius` of neuron i's dendrites or soma."""
        net = self.net
        labels = net.filled_mask
        somas = net.somas_mask
        boxes = find_objects(labels, max_label=self.N)
        sboxes = find_objects(somas, max_label=self.N)
        r = self.radius
        pad = r + int(self.cfg.bbox_pad_px)
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
        rows, cols = [], []
        for i, n in enumerate(net.neurons):
            bx = [b for b in (boxes[i], sboxes[i]) if b is not None]
            if not bx:
                continue
            y0 = max(min(b[0].start for b in bx) - pad, 0); y1 = min(max(b[0].stop for b in bx) + pad, self.H)
            x0 = max(min(b[1].start for b in bx) - pad, 0); x1 = min(max(b[1].stop for b in bx) + pad, self.W)
            own = net._branch_filled_counts[i][y0:y1, x0:x1] > 0
            own |= somas[y0:y1, x0:x1] == n.id_int
            dil = cv2.dilate(own.view(np.uint8), kernel, iterations=1).astype(bool)
            ys, xs = np.nonzero(dil)
            rows.append(np.full(len(ys), i, np.int64))
            cols.append((ys + y0).astype(np.int64) * self.W + (xs + x0))
        r_, c_ = np.concatenate(rows), np.concatenate(cols)
        return sparse.csr_matrix((np.ones(len(r_), np.float32), (r_, c_)), shape=(self.N, self.H * self.W))

    def contacts(self):
        """contact[i, j] = # axon pixels of j within the contact radius of i's dendrites/soma."""
        C = (self._dendrite_matrix() @ self._axon_matrix().T).toarray().astype(float)
        np.fill_diagonal(C, 0.0)
        return C

    def connectivity(self):
        return self.contacts() / max(self.ref, 1e-9)

    # ------------------------------------------------------------------ outputs
    def mask(self):
        m = np.zeros(self.H * self.W, bool)
        for ax in self.axons:
            m[ax.pixels()] = True
        return m.reshape(self.H, self.W)

    def stats(self):
        return dict(axon_cable_mean=float(np.mean([a.cable_length() for a in self.axons])),
                    axon_reach_mean=float(np.mean([a.max_path() for a in self.axons])),
                    axon_branches_mean=float(np.mean([len(a.branches) for a in self.axons])),
                    axon_target_mean=float(np.mean(self.targets)))


def render_axo_frame(branch_net, axons: AxonNetwork | None, dendrite_cmap="summer",
                     axon_color=(255, 80, 200), background=(0, 0, 0)):
    """Canvas-resolution RGB image: dendrites coloured by neuron, somas, axons on top."""
    import matplotlib.pyplot as plt
    H, W = branch_net.filled_mask.shape
    img = np.zeros((H, W, 3), np.uint8); img[:] = background
    N = len(branch_net.neurons)
    lut = (255 * np.asarray(plt.get_cmap(dendrite_cmap)(np.linspace(0, 1, max(N, 1))))[:, :3]).astype(np.uint8)
    lab = branch_net.filled_mask.astype(np.int64)
    sel = (lab > 0) & (lab <= N)
    img[sel] = lut[lab[sel] - 1]
    if axons is not None:
        img[axons.mask()] = axon_color
    som = branch_net.somas_mask.astype(np.int64)
    s = (som > 0) & (som <= N)
    img[s] = (0.6 * lut[som[s] - 1] + 0.4 * 255).astype(np.uint8)
    return img
