import sys, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, '/home/apd/Projects/Growing')
import avalanches as av, probes, scan_frozen as sf
from timescales import slope
d = glob.glob('/home/apd/Projects/Growing/sweeps/s8_timescales/runs/*c7f9501526')[0]
lif, M = sf.load_run(d); lif.ext_rate = 0.00005
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for i, g in enumerate([0.8, 0.95, 1.05, 1.1, 1.15]):
    S = probes.lif_run(g * M, lif, 2_000_000, rng=np.random.default_rng(1))
    r = av.analyze(S, 1, 50.0, 0.02); s = r['sizes']; T = r['durations']
    f = r['size_fit']
    print(f"g={g} σ={sf.sigma_of(M, lif, g):.2f} n={len(s)} slope[2,200]={slope(s):.2f} slopeT={slope(T,2,60):.2f} KS τ={f['alpha']:.2f} xmin={f['xmin']} tail={f['n_tail']} m={r['branching']['m']:.3f} rate={S.mean()/lif.dt:.4f}", flush=True)
    c, p = av._log_binned_pdf(s); ax[0].loglog(c, p, 'o-', ms=3, color=f'C{i}', label=f'g={g}')
    c, p = av._log_binned_pdf(T); ax[1].loglog(c, p, 'o-', ms=3, color=f'C{i}', label=f'g={g}')
x = np.array([1, 1e3]); ax[0].loglog(x, 0.3 * x ** -1.5, 'k--', label='−3/2'); ax[1].loglog(x[:1].tolist()+[100], 0.3*np.array([1,100.])**-2, 'k--', label='−2')
for a in ax: a.legend(fontsize=7)
ax[0].set_title('P(S), frozen network, h = 5e-5'); ax[1].set_title('P(T)')
fig.tight_layout(); fig.savefig('/home/apd/Projects/Growing/reports/s8_frozen_lowdrive_pS.png', dpi=100)
