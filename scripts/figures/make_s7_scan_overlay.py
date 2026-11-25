import json, glob, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__)))))
R = __import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.dirname(__import__('os').path.abspath(__file__)))) + '/'
fss = {r['N']: r for r in json.load(open(R + 'reports/s7_finite_size_fss.json'))['rows']}
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
items = []
for f in glob.glob(R + 'reports/s7_scan_*.json'):
    name, rows = list(json.load(open(f)).items())[0]
    d = glob.glob(R + 'sweeps/s7_finite_size/runs/' + name)[0]
    N = json.load(open(d + '/config.json'))['morph']['num_neurons']
    items.append((N, rows))
for i, (N, rows) in enumerate(sorted(items)):
    sig = np.array([r['sigma'] for r in rows])
    for a, k in zip(ax, ['rate', 'm', 'dr']):
        a.plot(sig, [r[k] for r in rows], 'o-', ms=4, color=f'C{i}', label=f'N={N}')
    for a in ax:
        a.axvspan(fss[N]['sigma_low'], fss[N]['sigma_high'], color=f'C{i}', alpha=0.08)
        a.axvline(fss[N]['sigma_mean'], color=f'C{i}', ls=':', lw=1)
ax[0].axhline(0.05, color='k', ls='--', lw=0.8, label='target rate'); ax[0].set_yscale('log')
ax[1].axhline(1, color='0.6', lw=0.8)
for a, t in zip(ax, ['frozen firing rate', 'MR branching ratio m (frozen)', 'dynamic range Δ (dB)']):
    a.set_title(t); a.set_xlabel('σ = mean Σ_i p_ij (coupling scaled by g)'); a.legend(fontsize=7)
fig.suptitle('Coupling scans of the grown s7 networks; dotted = mean σ during growth, band = low/high tertile means')
fig.tight_layout(); fig.savefig(R + 'reports/s7_coupling_scans.png', dpi=110)
