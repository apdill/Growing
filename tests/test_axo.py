"""Checks for the axo-dendritic mode: python tests/test_axo.py"""
import os, random, sys, time
import cv2
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from growing.sim import run_sim as rs
from growing.sim.axon_nets import AxonConfig, render_axo_frame


def build(duration=400, N=12, size=700, axon=None, seed=4):
    morph, branch, lif, _ = rs.default_configs()
    morph.num_neurons, morph.width, morph.height, morph.morphology_seed = N, size, size, seed
    morph.connectivity_mode = "axo_dendritic"
    lif.duration, lif.animate = duration, False
    lif.external_drive, lif.spike_mode, lif.ext_rate, lif.morphology_update_interval = "poisson", "branching", 0.0005, 5
    random.seed(seed); np.random.seed(seed)
    sim = rs.BranchMorphPlasticLIFSim(morph, lif, branch, axon or AxonConfig(length_px=400, growth_rate=3.0))
    t0 = time.time(); res = sim.simulate(analysis_cfg=None)
    return sim, res, time.time() - t0


def brute_contacts(sim):
    """contact[i, j] by explicit full-canvas masks (slow reference)."""
    net, ax = sim.branch_net, sim.axons
    r = ax.radius
    k = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2 * r + 1, 2 * r + 1))
    N = sim.N
    amask = [np.zeros(ax.H * ax.W, bool) for _ in range(N)]
    for j, a in enumerate(ax.axons):
        amask[j][a.pixels()] = True
    C = np.zeros((N, N))
    for i, n in enumerate(net.neurons):
        own = (net._branch_filled_counts[i] > 0) | (net.somas_mask == n.id_int)
        dil = cv2.dilate(own.astype(np.uint8), k).astype(bool).ravel()
        for j in range(N):
            if i != j:
                C[i, j] = np.count_nonzero(dil & amask[j])
    return C


if __name__ == "__main__":
    sim, res, wall = build()
    M = res["connectivity_matrix"]
    print(f"wall {wall:.1f}s  axon stats {res['axon_stats']}")
    print(f"nonzero {np.count_nonzero(M)}/{M.size - len(M)}  asymmetry |M-Mt|/|M| = "
          f"{np.abs(M - M.T).sum() / max(np.abs(M).sum(), 1e-12):.3f}  diag zero {np.all(np.diag(M) == 0)}")
    C = sim.axons.contacts(); B = brute_contacts(sim)
    print("sparse == brute-force:", np.array_equal(C, B), " max |diff|", np.abs(C - B).max())
    t0 = time.time(); sim.axons.connectivity(); print(f"connectivity time {1000 * (time.time() - t0):.1f} ms (N={sim.N})")
    # retraction: shrink caps and regrow
    ax0 = sim.axons.axons[0]; L0 = ax0.max_path(); ax0.set_cap(L0 / 2); ax0.grow()
    print(f"retraction: reach {L0:.0f} -> {ax0.max_path():.0f} (cap {ax0.cap:.0f}), branches {len(ax0.branches)}")
    img = render_axo_frame(sim.branch_net, sim.axons)
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "reports", "axo_test_frame.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    import imageio.v2 as imageio; imageio.imwrite(out, img); print("frame ->", out)
