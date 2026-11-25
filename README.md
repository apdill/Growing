# Growing

Homeostatic structural plasticity in 2-D networks of neurons with stochastic fractal dendrites:
neurites grow/retract to hold a target firing rate, connectivity comes from arbor contacts, and the
resulting dynamics are analysed for criticality (avalanches, branching ratio, dynamic range) and
geometry (box-counting dimension, lacunarity, wiring statistics).

## Layout

```
growing/                    importable package
  paths.py                  repository locations (specs/, sweeps/, reports/, ...)
  sim/
    branch_nets.py          dendritic morphology generation (stochastic fractal branching)
    run_sim.py              plastic network simulation (configs, growth rules, connectivity)
    lif_dynamics.py         shared neuron update (threshold / escape-noise LIF, branching model)
    lif_neuron.py           display neuron used by run_sim
  analysis/
    avalanches.py           avalanche detection, power-law fits, MR branching ratio
    probes.py               frozen-network probes (dynamic range, spectrum, graph stats)
    fractal.py              box counting via FracStack (min counts, OLS), lacunarity, geometry
scripts/
  run_single.py             one simulation with default configs
  sweeps/sweep.py           parallel, resumable parameter sweeps   (specs/*.json -> sweeps/<name>/)
  sweeps/scan_frozen.py     coupling / matched-Λ scans of grown networks
  analysis/                 sweep summaries, finite-size, timescales, re-analysis, calibration
  figures/                  figure scripts for reports
  external/                 analysis of third-party data (e.g. Stiber-lab simulations)
specs/                      sweep specifications (drafts/ = not yet runnable)
```

Outputs (`sweeps/`, `reports/`, `summary_plots/`, `external/`) and planning material are local only.

## Usage

```bash
python scripts/run_single.py
python scripts/sweeps/sweep.py specs/s9_fss_scaled_drive.json --workers 8 --burn-in 0.2
python scripts/analysis/finite_size.py s9_fss_scaled_drive
python -m growing.analysis.fractal          # box-counting validation on known fractals
```

Fractal analysis needs FracStack (github.com/apdill/FracStack) next to this repository, or
`FRACSTACK_PATH` pointing to it.
