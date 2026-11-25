"""Fingerprint of a short default-mode (dendro-dendritic) run: hashes of M, spikes and lengths.
python tests/legacy_fingerprint.py  -> prints one line; must be identical before/after changes that
claim to leave the legacy mode untouched."""
import hashlib, os, random, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from growing.sim import run_sim as rs

morph, branch, lif, analysis = rs.default_configs()
morph.num_neurons, morph.width, morph.height, morph.morphology_seed = 8, 512, 512, 3
lif.duration, lif.animate = 150, False
lif.external_drive, lif.spike_mode, lif.ext_rate, lif.morphology_update_interval = "poisson", "branching", 0.0005, 5
random.seed(3); np.random.seed(3)
sim = rs.BranchMorphPlasticLIFSim(morph, lif, branch)
res = sim.simulate(analysis_cfg=None)
h = lambda a: hashlib.md5(np.ascontiguousarray(np.asarray(a)).tobytes()).hexdigest()[:12]
print("M", h(np.round(res["connectivity_matrix"], 10)), "spikes", h(np.asarray(res["spikes_list"], bool)),
      "len", h(np.round(np.asarray(res["length_list"]), 6)))
