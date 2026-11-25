# Critiques, confounds and best practice for claims of neuronal criticality (state of the debate 2020–2026)

Verification legend used throughout:
- **[V]**: verified this session. Title, venue and DOI were confirmed by search or the Semantic Scholar API. For abstracts, the claim matches the retrieved abstract.
- **[M]**: from prior knowledge and not re-checked this session. The citation is very likely correct, but check the DOI and details before quoting.
- **[U]**: uncertain. Details are partly or wholly unconfirmed.

Project pitfalls cross-referenced below: (P1) a KS-chosen xmin gave misleading exponents; (P2) pooled power laws were mixtures of sub- and supercritical epochs; (P3) avalanches overlap under sparse drive; (P4) σ (mean column sum of transmission probabilities) differs from the measured m.

---

## Q1. Power laws without criticality: the main critiques

### Takeaway
Power-law avalanche distributions, and even the crackling-noise scaling relation, are now widely accepted as **necessary but not sufficient** evidence of criticality. Several well-understood non-critical mechanisms reproduce them:
- thresholded noise;
- a slowly varying common drive;
- neutral (causal) dynamics;
- latent-variable coupling;
- subsampling and overlapping electrodes;
- non-critical E/I networks.

A credible claim must exclude these mechanisms with matched null models and independent evidence of proximity to a phase transition.

### Cited findings
**Touboul & Destexhe 2010** [M], *PLoS One* 5(2):e8982, "Can power-law scaling and neuronal avalanches arise from stochastic dynamics?", doi:10.1371/journal.pone.0008982 — https://doi.org/10.1371/journal.pone.0008982
- They show that thresholded stochastic processes and simple non-critical network models generate apparent power-law avalanche distributions that pass standard statistical tests. They also argue that the original LFP-peak analysis produces spurious power laws largely independent of the underlying spiking.
- A later summary states that "the method of analyzing the amplitude of negative LFP peaks was shown to produce spurious power laws scalings regardless of the spike activity of cells" — [arXiv 1503.08033 summary](https://arxiv.org/html/1503.08033v5) [V, snippet].

**Touboul & Destexhe 2017** [V title; M venue], "Power-law statistics and universal scaling in the absence of criticality", *Phys. Rev. E* 95, 012413, doi:10.1103/PhysRevE.95.012413 — [arXiv 1503.08033](https://arxiv.org/abs/1503.08033)
- The paper argues that power-law avalanche statistics and universal-looking scaling "can arise naturally in neuronal systems that are not at criticality", notably in the **synchronous-irregular** regime of balanced networks. It concludes that "power-law distributions of avalanches do not necessarily reveal that the network is operating at criticality" ([arXiv](https://arxiv.org/html/1503.08033v5), [V via search snippet]).
- [M] It also reports that a stochastic surrogate (thresholded Ornstein–Uhlenbeck-type noise) satisfies the crackling-noise scaling relation.

**Destexhe & Touboul 2021** [V], "Is there sufficient evidence for criticality in cortical systems?", *eNeuro* 8(2): ENEURO.0551-20.2021, doi:10.1523/ENEURO.0551-20.2021 — https://doi.org/10.1523/ENEURO.0551-20.2021 ; [PMC8059881](https://pmc.ncbi.nlm.nih.gov/articles/PMC8059881); [arXiv 2002.08813](https://arxiv.org/abs/2002.08813)
- **Correction to the brief:** the paper appeared in *eNeuro*, not *Biology* [V].
- It tests two non-critical systems: a **Brunel (sparse E/I) network** and a **stochastic surrogate**. Both pass the criticality tests (power laws plus the scaling relation between exponents) used in a prior study. That study is [M] Fontenele et al. 2019 PRL 122:208101, doi:10.1103/PhysRevLett.122.208101. The authors conclude that the methodology "may not be enough to establish criticality" and that the question "remains open" [V, abstract-level via search].

**Priesemann & Shriki 2018** [V], "Can a time varying external drive give rise to apparent criticality in neural systems?", *PLoS Comput Biol* 14(5):e1006081, doi:10.1371/journal.pcbi.1006081 — https://doi.org/10.1371/journal.pcbi.1006081
- Homogeneous Poisson activity cannot produce power laws. **Slowly varying inhomogeneous Poisson activity** with no interactions *can* produce "approximate power laws with cutoffs that are comparable to those typically observed in experiments". They argue that telling true from apparent criticality ultimately "requires causal interventions, which allow separating internal system properties from externally imposed ones" [V, abstract].
- This bears directly on P2: a pooled mixture of epochs with different effective rates or m is the same mechanism.

**Martinello, Hidalgo, Maritan, di Santo, Plenz & Muñoz 2017** [M], "Neutral theory and scale-free neural dynamics", *Phys. Rev. X* 7, 041071, doi:10.1103/PhysRevX.7.041071 — https://doi.org/10.1103/PhysRevX.7.041071
- [M] When avalanches are defined *causally* (tracking which activation caused which), neutral dynamics among many coexisting, causally independent avalanches yields scale-free avalanche distributions. This holds even when the overall network is **not** at its critical point, for example in an active phase.
- The paper warns that causally defined avalanche statistics and time-binned statistics can disagree. Scale-free causal avalanches are therefore not by themselves evidence of criticality.
- Relevance: P3, and any simulation that uses ground-truth lineage to define avalanches.

**Nonnenmacher, Behrens, Berens, Bethge & Macke 2017** [M], "Signatures of criticality arise from random subsampling in simple population models", *PLoS Comput Biol* 13(10):e1005718, doi:10.1371/journal.pcbi.1005718 — https://doi.org/10.1371/journal.pcbi.1005718
- [M] This concerns the "thermodynamic" signature of criticality (a diverging specific heat of maximum-entropy models fitted to binarized population activity). Divergence appears generically from random subsampling of simple models that have shared input/population-rate fluctuations, so it is a poor criticality test.
- It applies to statistical-mechanics claims rather than avalanche claims, but it is often cited alongside them.

**Morrell, Sederberg & Nemenman 2021** [V], "Latent dynamical variables produce signatures of spatiotemporal criticality in large biological systems", *Phys. Rev. Lett.* 126, 118302, doi:10.1103/PhysRevLett.126.118302 — [accepted ms](https://link.aps.org/accepted/10.1103/PhysRevLett.126.118302)
- Non-interacting neurons coupled to multiple slow latent fields reproduce the scaling seen under **phenomenological renormalization-group coarse-graining** (Meshulam et al. style), within experimental error [V, search summary].

**Morrell, Nemenman & Sederberg 2024** [V], "Neural criticality from effective latent variables", *eLife* 12:RP89337, doi:10.7554/eLife.89337 — https://elifesciences.org/articles/89337 ; [arXiv 2301.00759](https://arxiv.org/abs/2301.00759)
- The model is extended to **avalanche criticality**. Populations coupled to multiple latent variables show avalanche power laws and scaling "without fine-tuning of model parameters", over a broader parameter range than a single quasi-static latent variable [V].
- Author order differs between the two papers. Confirm it on the publisher pages [V].

**Neto, Spitzner & Priesemann 2022** [V], "Sampling effects and measurement overlap can bias the inference of neuronal avalanches", *PLoS Comput Biol* 18(11):e1010678, doi:10.1371/journal.pcbi.1010678 — [PMC9733887](https://pmc.ncbi.nlm.nih.gov/articles/PMC9733887)
- Overlapping electrode fields (coarse LFP/MEA sampling) correlate the signals and can make **non-critical** branching networks look critical (power-law avalanches). Spike-based sampling avoids this specific bias. This may explain conflicting coarse versus spike-level results [V].

**Faqeeh, Osat, Radicchi & Gleeson 2019** [M], "Emergence of power laws in noncritical neuronal systems", *Phys. Rev. E* 100, 010401(R), doi:10.1103/PhysRevE.100.010401
- [M] Network-structural mechanisms produce power-law avalanches away from the critical point.

**Self-organized quasi-criticality** [M]: Bonachela, de Franciscis, Torres & Muñoz 2010, "Self-organization without conservation: are neuronal avalanches generically critical?", *J. Stat. Mech.* P02015, doi:10.1088/1742-5468/2010/02/P02015
- [M] Non-conservative SOC with homeostatic or synaptic-resource recovery does **not** sit exactly at the critical point. It hovers around it, with excursions to both sides. Pooled avalanche statistics then show apparent power laws with sub- and supercritical contributions. This is the textbook precedent for P2.

**Sipling, Zhang & Di Ventra 2026** [V], "A Critical Assessment of the Brain Criticality Hypothesis", arXiv:2604.21071 (v1 22 Apr 2026, v3 29 May 2026) — https://arxiv.org/abs/2604.21071
- Opinion piece. It argues that the nature of cortical criticality "is still unclear". Coupling of neurons to slowly varying energetic resources (memory) can alone create a "memory-induced long-range order" phase with scale-invariant correlations **without** criticality. It also notes that not all measured avalanche exponents are 3/2 [V, abstract and snippet].

**Kanders, Lorimer & Stoop 2017** [M], *Chaos* 27:047408, "Avalanche and edge-of-chaos criticality do not necessarily co-occur in neural networks", doi:10.1063/1.4978998
- [M] "Critical" in the avalanche sense and "critical" in the edge-of-chaos sense are distinct.

### Inferences
- For this project, every Q1 mechanism maps onto a concrete null that can be simulated cheaply:
  - inhomogeneous Poisson with the matched rate trajectory (Priesemann & Shriki);
  - latent-field model with matched single-neuron rates and pairwise correlations (Morrell);
  - an E/I or Brunel-type network known to be non-critical (Destexhe & Touboul);
  - homeostatic hovering (Bonachela et al.).
  
  Reviewers familiar with Destexhe & Touboul 2021 will expect at least one such "passes your tests but is not critical" control.
- P2 is a known artefact class. Reporting conditional distributions (avalanches binned by instantaneous m) is a direct and convincing answer.
- The latent-variable critique is strongest for coarse-grained/RG signatures and for avalanches under shared slow input. A homeostatic model with a slowly varying homeostatic variable is structurally similar to "latent variable" coupling. Referees may ask whether the homeostatic variable itself acts as the latent field.

### Gaps
- I did not re-read the full texts of Touboul & Destexhe 2010/2017, Martinello 2017 or Nonnenmacher 2017 this session. Their mechanism summaries are from memory [M].
- I found no formal published rebuttal to Morrell et al. 2024, and no 2025–2026 response paper to Destexhe & Touboul 2021 beyond Hengen & Shew 2025 (Q3).

---

## Q2. Best-practice tests: fitting, scaling relations, thresholding, subsampling, binning

### Takeaway
The current standard is a **battery** of tests, none sufficient alone:
- principled truncated-power-law MLE fitting with model comparison and goodness of fit;
- three independently estimated exponents that satisfy the crackling-noise relation, plus shape collapse with a matching exponent;
- robustness to bin size and threshold;
- explicit handling of subsampling (MR estimator, subsampling scaling);
- in models, **finite-size scaling** and an independently located phase transition.

### Cited findings

#### Crackling noise and shape collapse
- **Sethna, Dahmen & Myers 2001** [M], "Crackling noise", *Nature* 410:242–250, doi:10.1038/35065675. Derives the exponent relation ⟨S⟩(T) ∝ T^{1/σνz} and shows that mean avalanche shapes at different durations collapse onto a universal scaling function.
- **Friedman, Ito, Brinkman, Shimono, DeVille, Dahmen, Beggs & Butler 2012** [M], "Universal critical dynamics in high resolution neuronal avalanche data", *Phys. Rev. Lett.* 108, 208102, doi:10.1103/PhysRevLett.108.208102.
  - Introduces three tests for neuronal avalanches: P(S) ∝ S^−τ, P(T) ∝ T^−α, and ⟨S⟩(T) ∝ T^γ, with the **scaling relation γ = (α−1)/(τ−1)**, plus **shape collapse** with the same γ.
  - [M] Cultures matched the relation and collapse, but with γ ≈ 1.3, not the mean-field directed-percolation (MF-DP) value of 2.
  - Mean-field DP / critical branching predicts τ = 3/2, α = 2, γ = 2 [M; standard].
- **Necessary, not sufficient:** Destexhe & Touboul 2021 [V] show non-critical systems satisfying the exponent relation (Q1).
- **Ma, Turrigiano, Wessel & Hengen 2019** [M], *Neuron* 104:655, doi:10.1016/j.neuron.2019.08.031. Defined the **deviation from criticality coefficient (DCC)** = |γ_predicted − γ_fitted|. The Hengen/Shew group still uses DCC and a newer whole-time-series temporal scale-invariance metric "D2". They report the two usually agree but diverge in states with regular silences such as NREM sleep — [Brain Inspired BI 216 transcript, 2025](https://www.thetransmitter.org/wp-content/uploads/2025/07/BI-216-transcript-hengen-shew-1.pdf) [V transcript].
- **Capek et al. 2023** [U], *Nat. Commun.* 14:2555 (?), "Parabolic avalanche scaling in the synchronization of cortical cell assemblies", doi:10.1038/s41467-023-37976-x [U]. [M] Reports γ ≈ 2 (parabolic shapes) in properly resolved cell-assembly data and argues lower γ values reflect methodological issues. Verify before citing.

#### Fitting standards
- **Clauset, Shalizi & Newman 2009** [M], "Power-law distributions in empirical data", *SIAM Review* 51:661–703, doi:10.1137/070710111. Covers:
  - MLE for the exponent;
  - x_min chosen by minimizing the KS distance;
  - a semi-parametric bootstrap goodness-of-fit p-value (p > 0.1 for plausibility);
  - likelihood-ratio (Vuong) tests against lognormal, exponential, stretched exponential and truncated power law.
  
  The KS-chosen x_min can select a short, unrepresentative tail range when data are curved (P1). The method assumes no upper cutoff, which neuronal avalanche data always have.
- **Deluca & Corral 2013** [M], "Fitting and goodness-of-fit test of non-truncated and truncated power-law distributions", *Acta Geophysica* 61:1351–1394, doi:10.2478/s11600-013-0154-9. Covers MLE for **doubly truncated** power laws on [a, b], a KS-based goodness-of-fit with Monte Carlo p-values, and selecting the range by **maximizing the fitted range (log b/a)** subject to p > p_c. This is preferred over pure KS-minimization because it avoids tiny fit windows.
- **Marshall, Timme, Bennett, Ripp, Lautzenhiser & Beggs 2016** [M], "Analysis of power laws, shape collapses, and neural complexity: new techniques and MATLAB support via the NCC Toolbox", *Front. Physiol.* 7:250, doi:10.3389/fphys.2016.00250. Implements Deluca–Corral truncated MLE, the scaling relation and shape collapse. It recommends reporting exponents with their fit ranges and checking consistency among τ, α and γ.
- **Alternative toolkit** [M]: Alstott, Bullmore & Plenz 2014, *PLoS One* 9:e85777, "powerlaw: a Python package", doi:10.1371/journal.pone.0085777. Also implements the Clauset method and LLR comparisons, but defaults to KS x_min.
- **Field-level note** [V snippet]: "no standardized protocol of analysis has been established so far" — [Front. Neural Circuits 2022 mini-review](https://www.frontiersin.org/articles/10.3389/fncir.2022.911245) (authors not confirmed). Tian et al., "Theoretical foundations of studying criticality in the brain", *Network Neuroscience* 2022 ([arXiv 2306.05635](https://arxiv.org/abs/2306.05635)) [V title; M venue] reviews estimators and their biases at length.

#### Thresholding
- **Villegas, di Santo, Burioni & Muñoz 2019** [V], "Time-series thresholding and the definition of avalanche size", *Phys. Rev. E* 100, 012133, doi:10.1103/PhysRevE.100.012133 — https://doi.org/10.1103/PhysRevE.100.012133
  - Thresholding continuous signals hides a "peril". Shown analytically and numerically for an Ornstein–Uhlenbeck process: with **relatively large thresholds**, and with avalanche size defined as **total area above zero** (rather than above threshold), "true asymptotic scaling behavior is not seen". Observations are dominated by transients, which "leads to misinterpretations of the resulting scaling regimes" [V, abstract].
  - Recommendation [V]: define size consistently with the threshold (area above threshold) and verify that exponents are stable as the threshold varies.
- A 2013 study [V snippet, arXiv 1309.3535, "Identification of criticality in neuronal avalanches II"] notes that thresholded stochastic processes yield apparent power laws that only stringent testing rejects, while stringent testing can also wrongly reject true power laws — https://arxiv.org/pdf/1309.3535

#### Subsampling
- **Priesemann, Munk & Wibral 2009** [M], "Subsampling effects in neuronal avalanche distributions recorded in vivo", *BMC Neurosci.* 10:40, doi:10.1186/1471-2202-10-40. Subsampling a critical system distorts P(S) away from a power law, so a subsampled critical system can look subcritical or supercritical.
- **Priesemann et al. 2014** [V via PMC listing; M details], "Spike avalanches in vivo suggest a driven, slightly subcritical brain state", *Front. Syst. Neurosci.* 8:108, doi:10.3389/fnsys.2014.00108 — [PMC4068003](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4068003/). [M] Compares in vivo spike avalanches with subsampled models under varying drive and **bin sizes**. The data are best matched by a driven, slightly subcritical regime, not a critical one.
- **Levina & Priesemann 2017** [M], "Subsampling scaling", *Nat. Commun.* 8:15140, doi:10.1038/ncomms15140. For critical systems, avalanche distributions from different sampling fractions *p* collapse when rescaled, P_sub(s) ∝ p^{…}F(s/p). This allows inferring full-system criticality from subsampled data, and failure to collapse indicates non-criticality.
- **Wilting & Priesemann 2018** [M], "Inferring collective dynamical states from widely unobserved systems", *Nat. Commun.* 9:2325, doi:10.1038/s41467-018-04725-4. The **multistep regression (MR) estimator** fits lagged regression slopes r_k = b·m^k. The prefactor b absorbs subsampling bias, so m is estimated consistently. Conventional one-step regression (or ratio of descendants/ancestors) is biased towards 0 under subsampling.
  - Applied in vivo, it gives m ≈ 0.98, an autocorrelation time of hundreds of ms ("reverberating").
  - Toolbox: Spitzner et al. 2021, "MR. Estimator, a toolbox to determine intrinsic timescales from subsampled spiking activity", *PLoS One* 16:e0249447, doi:10.1371/journal.pone.0249447 [M].
- **Neto et al. 2022** [V]: measurement overlap (Q1).
- **Carvalho et al. 2021** [M], *Front. Neural Circuits* 14:576727, "Subsampled directed-percolation models explain scaling relations experimentally observed in the brain", doi:10.3389/fncir.2020.576727. [M] Subsampling changes the individual exponents τ and α but approximately preserves the γ scaling relation. This helps explain non-MF experimental exponents.

#### Binning
- **Beggs & Plenz 2003** [M], "Neuronal avalanches in neocortical circuits", *J. Neurosci.* 23:11167–11177, doi:10.1523/JNEUROSCI.23-35-11167.2003. The size exponent depended on bin width Δt (α(Δt) varied systematically). They proposed Δt = the mean inter-event interval (IEI) of the pooled population as the natural bin, which gave τ ≈ 3/2 and σ ≈ 1.
- **Priesemann 2014** [M]: in vivo spike avalanches change shape strongly with Δt. Under a driven, non-separated-timescale regime, large bins merge independent cascades, so the Δt-dependence is itself diagnostic.
- **Das & Levina 2019** [V], "Critical neuronal models with relaxed timescale separation", *Phys. Rev. X* 9, 021062 — [arXiv 1808.04196](https://arxiv.org/abs/1808.04196). An analytic and numerical treatment of what happens when external input arrives during an ongoing avalanche (no strict separation of timescales). They show how bin size and drive rate interact with the apparent exponent, recovering 1.5 for perfect separation [V, summary]. This bears directly on P3.

### Inferences
- P1 fix: KS-minimized x_min is no longer considered sufficient for neural data. Reviewers will expect:
  - doubly truncated MLE (Deluca–Corral / NCC);
  - an exponent-stability plot versus x_min and x_max;
  - likelihood-ratio comparisons with lognormal, exponential and truncated power law;
  - a stated fit range spanning at least 1.5–2 decades.

  Curved distributions under fixed-window MLE should be reported as evidence *against* (or for distance from) criticality, not fitted away.
- P3 fix: in simulation the ground-truth causal tree is available. Report both causal avalanches (each external seed's lineage) and binned avalanches. Quantify the overlap fraction (the probability that a new external input lands while another cascade is active, ≈ 1 − exp(−h·N·⟨T⟩)). Show results in a regime where overlap is small or corrected (Das & Levina 2019). Note that Martinello et al. 2017 warn that causal avalanches can be scale-free off-criticality, so they do not settle the question alone.
- The scaling relation should be tested with **three independently fitted** exponents (τ, α, and γ from ⟨S⟩(T)), plus γ from shape collapse. Agreement with a specific universality class (MF-DP for random sparse networks) is more persuasive than agreement with the relation alone.
- In models, the gold standard that experiments cannot provide is **finite-size scaling**: cutoffs scale with N (S_c ∝ N^{…}), exponents are N-independent, and shapes collapse across N.

### Gaps
- The exact γ ranges reported across 2020–2026 meta-analyses were not retrieved. Whether Hengen & Shew 2025 tabulate γ across 140 datasets was not confirmed.
- Detailed parameter recommendations (minimum decades, sample sizes) from the NCC toolbox paper were not re-read [M].

---

## Q3. Reviews of the field's status: what is currently accepted (2018–2026)

### Takeaway
The field has largely moved from "the brain is critical" towards these positions:
- the cortex operates **near** a critical point, usually slightly subcritical or "reverberating", with a state-dependent distance to criticality;
- that distance (not a binary verdict) is the useful observable;
- power laws alone are not accepted as evidence.

Proponents (Hengen & Shew 2025, Plenz 2021, O'Byrne & Jerbi 2022) argue that criticality is a homeostatically regulated setpoint. Critics (Destexhe & Touboul 2021, Morrell et al. 2021/2024, Sipling et al. 2026) maintain that current tests cannot rule out non-critical explanations.

### Cited findings
- **Muñoz 2018** [M], "Colloquium: Criticality and dynamical scaling in living systems", *Rev. Mod. Phys.* 90, 031001, doi:10.1103/RevModPhys.90.031001.
  - [M] A broad physics review. It stresses that the type of phase transition matters (e.g., absorbing-state/DP versus synchronization transitions; di Santo et al. 2018 PNAS argue avalanches emerge at the edge of synchronization). It separates "dynamical" from "statistical" criticality and warns about the pitfalls covered in Q1/Q2.
- **Wilting & Priesemann 2019** [V], "25 years of criticality in neuroscience — established results, open controversies, novel concepts", *Curr. Opin. Neurobiol.* 58:105–111, doi:10.1016/j.conb.2019.08.002.
  - "Experimental results are still contradictory".
  - They propose a "parsimonious solution" (the reverberating, slightly subcritical regime) that "(i) unifies the contradictory experimental results, (ii) avoids disadvantages of a critical state, and (iii) enables rapid, adaptive tuning of network properties to task requirements" [V, abstract].
- **Plenz, Ribeiro, Miller, Kells, Vakili & Capek 2021** [V], "Self-organized criticality in the brain", *Front. Phys.* 9:639389, doi:10.3389/fphy.2021.639389.
  - A proponent review of isolated cortex preparations (organotypic cultures, slices). Avalanches show "robust power law scaling in avalanche sizes with a slope of −3/2 and a critical branching parameter of 1".
  - Dynamics are "homeostatically regulated, exhibit separation of time scales". Avalanches are embedded in up-states and nested θ/γ oscillations [V, abstract].
- **O'Byrne & Jerbi 2022** [V], "How critical is brain criticality?", *Trends Neurosci.* 45(11):820–837, doi:10.1016/j.tins.2022.08.007.
  - Reviews edge-of-chaos and **near-criticality** data. It makes "a case for the **distance to criticality** as a useful metric for probing cognitive states and mental illness" [V, abstract].
- **Beggs 2022** [M/U], *The Cortex and the Critical Point: Understanding the Power of Emergence*, MIT Press, ISBN 9780262544030 [U ISBN] — https://mitpress.mit.edu/ (exact URL not retrieved).
  - [M] A book-length proponent account. Includes the quasicriticality framing (Q4) and discusses critiques.
- **Hengen & Shew 2025** [V metadata], "Is criticality a unified setpoint of brain function?", *Neuron* 113(16):2582–2598.e2, doi:10.1016/j.neuron.2025.05.020 — https://www.cell.com/neuron/fulltext/S0896-6273(25)00391-5 ; preprint [bioRxiv 10.1101/2024.09.02.610815](https://www.biorxiv.org/content/10.1101/2024.09.02.610815).
  - Argues criticality is a homeostatic setpoint of healthy brain function [V].
  - Abstract (via search) reports a **meta-analysis of 140 datasets (2003–2024)**. It concludes "a long-standing controversy is the product of a methodological choice with no bearing on underlying dynamics" [V phrase, via search snippet].
  - **Which methodological choice is not confirmed [U].** The podcast transcript suggests it concerns how avalanches are defined (spatiotemporal versus purely temporal analysis, with the population-sum first step that "eliminates spatial correlations entirely") and bin-width requirements [V transcript content; U as the paper's specific claim].
  - The podcast host says the review surveys ">300" studies; the abstract says 140 datasets [V transcript; numbers differ].
  - Uses **DCC** and **D2** as distance-to-criticality metrics [V transcript].
- **Fosque, Shew, Ching & Hengen 2025** [V], "Two views of the brain are reconciled by a unifying principle of maximal information processing", bioRxiv doi:10.1101/2025.11.25.690580 ([PMC12697366](https://pmc.ncbi.nlm.nih.gov/articles/PMC12697366)).
  - Defines criticality for **low-dimensional** systems through **temporal scale invariance**. Finds two mechanisms at bifurcations (deterministic and noise-driven). Only some canonical bifurcations are critical, and information-processing capacity peaks at critical bifurcations [V, abstract].
- **Morales, di Santo & Muñoz 2023** [V], "Quasiuniversal scaling in mouse-brain neuronal activity stems from edge-of-instability critical dynamics", *PNAS* 120 (doi:10.1073/pnas.2208998120).
  - Combines phenomenological RG with dynamical-state inference on thousands of neurons. Finds "quasi-universal" scale invariance across regions, with all areas "near the edge of instability". The framework separates quasi-universal background activity from non-universal input-related activity [V, abstract].
- **Ma et al. 2019** [M]: rat V1 firing returns to near-criticality (DCC) after monocular deprivation, suggesting homeostatic tuning to criticality in vivo. This is the empirical basis for "criticality as homeostatic setpoint".
- **Critics 2021–2026**: Destexhe & Touboul 2021 [V]; Morrell et al. 2021/2024 [V]; Sipling, Zhang & Di Ventra 2026 [V] (all in Q1). The podcast transcript names Destexhe/Touboul and the Nemenman/Sederberg/Morrell group as the main current critics and describes the exchange as "cordial" [V transcript].

### Inferences
- The 2025–2026 consensus language that reviewers are comfortable with is "near-critical", "distance to criticality", "reverberating", "quasicritical" and "edge of instability", **not** "at criticality".
- The proponents' current emphasis on **homeostatic tuning toward criticality** (Ma 2019; Hengen & Shew 2025) makes a homeostasis-driven growth model topical. But it must engage the reverberating/subcritical alternative (Wilting & Priesemann; Zierenberg 2018, Q4) and the latent-variable critique.

### Gaps
- The full text of Hengen & Shew 2025 could not be retrieved (403 and 429 errors). Their specific methodological argument, recommended tests and responses to Touboul/Destexhe and Morrell are **unverified**.
- Beggs 2022 book contents were not checked.
- An arXiv paper titled "A minimal mechanism to generate long timescales without fine tuning" (arXiv:2608.21586, Aug 2026) appeared in search results but was not examined. Its relevance is [U].

---

## Q4. "Slightly subcritical / reverberating" versus critical, and quasicriticality under input: framings for a growth model

### Takeaway
There are two well-developed, referee-accepted framings for a homeostatic network that ends up at m slightly below 1:
1. **The reverberating / homeostasis-plus-input framing** (Zierenberg, Wilting & Priesemann 2018). Under rate homeostasis, the external input strength *sets* the distance to criticality.
2. **Quasicriticality** (Williams-García et al. 2014; Fosque et al. 2021). Input pushes the system off the critical point along a Widom line, exponents shrink, but the scaling relation approximately holds.

Both predict that the operating point depends quantitatively on drive. That gives a growth model testable, falsifiable predictions beyond "power laws appear".

### Cited findings
- **Zierenberg, Wilting & Priesemann 2018** [V], "Homeostatic plasticity and external input shape neural network dynamics", *Phys. Rev. X* 8, 031018, doi:10.1103/PhysRevX.8.031018.
  - With spike-rate homeostasis, increasing input produces a sequence of states: "bursting, to close-to-critical, reverberating and irregular states".
  - In vitro bursting is consistent with "very low network input (< 0.1%)". In vivo suggests "on the order of 1% recorded spikes are input-driven, resulting in reverberating dynamics" [V, abstract].
  - This is the most directly relevant prior model for a homeostatic growth project.
- **Wilting et al. 2018** [M], "Operating in a reverberating regime enables rapid tuning of network states to task requirements", *Front. Syst. Neurosci.* 12:55, doi:10.3389/fnsys.2018.00055. Being slightly subcritical (m ≈ 0.98–0.99) allows small changes in m to produce large changes in timescale and sensitivity, without the costs of criticality (runaway activity, slow responses).
- **Wilting & Priesemann 2019** [V]: the reverberating regime is the "parsimonious solution" (Q3).
- **Fosque, Williams-García, Beggs & Ortiz 2021** [V], "Evidence for quasicritical brain dynamics", *Phys. Rev. Lett.* 126, 098101, doi:10.1103/PhysRevLett.126.098101 — [arXiv 2010.02938](https://arxiv.org/abs/2010.02938).
  - Measured exponents vary across species, individuals and time yet approximately satisfy the scaling relation.
  - Quasicriticality explains this: stronger external input moves cortex off criticality along a **Widom line**, exponent magnitudes decrease, and the scaling relation is approximately preserved. Tested with simulations and experimental data. Editors' Suggestion [V].
- **Williams-García, Moore, Beggs & Ortiz 2014** [M], "Quasicritical brain dynamics on a nonequilibrium Widom line", *Phys. Rev. E* 90, 062714, doi:10.1103/PhysRevE.90.062714. The original quasicriticality theory: under nonzero drive there is no true critical point, only a line of maximal susceptibility.
- **Fosque et al. 2022** [M/U], *Front. Comput. Neurosci.* 16:1037550, "Quasicriticality explains variability of human neural dynamics across life span" [U DOI].
- **Branching process with immigration** [M; textbook result]: if A_{t+1} has mean m·A_t + h, the stationary activity is A* = h/(1 − m), so **m = 1 − h/A***. With a homeostatic target rate fixed and external drive h > 0, the network settles at m strictly below 1, and the distance to criticality 1 − m equals the fraction of input-driven activity. This is the core of Zierenberg et al. 2018 [M for exact form].
- **σ versus m (P4)**: Larremore, Shew & Restrepo 2011 [V], "Predicting criticality and dynamic range in complex networks: effects of topology", *Phys. Rev. Lett.* 106, 058101, doi:10.1103/PhysRevLett.106.058101.
  - "The **largest eigenvalue** of the weighted network adjacency matrix governs the network dynamic range". When it "is exactly one, the system is in a critical state and its dynamic range is maximized". More homogeneous degree distributions give higher dynamic range [V, abstract].
- **Dynamic range** [M]: Kinouchi & Copelli 2006, "Optimal dynamical range of excitable networks at criticality", *Nat. Phys.* 2:348–351, doi:10.1038/nphys289. Shew, Yang, Petermann, Roy & Plenz 2009, *J. Neurosci.* 29:15595, doi:10.1523/JNEUROSCI.3864-09.2009 (dynamic range maximized at m ≈ 1 in cultures).

### Inferences
- **Explaining P4 (σ ≠ measured m).** The mean column sum σ is a mean-field proxy. The relevant critical quantity is the **largest eigenvalue Λ** of the transmission matrix (Larremore 2011). For heterogeneous or degree-correlated grown networks, Λ can differ substantially from the mean column sum. Measured m additionally deviates from Λ because of:
  - **saturation/collisions**: a target receiving several simultaneous transmissions fires only once, i.e., 1 − ∏(1 − p) < Σp, which lowers the effective m;
  - **refractoriness**, which lowers it further;
  - **external drive contaminating descendant counts**: naive A_{t+1}/A_t includes input, so it is biased upward unless input is subtracted;
  - **bin width** relative to the transmission delay;
  - **subsampling** (biased downward unless MR is used).

  Recommended practice: report σ, Λ (and the left/right eigenvector-weighted m), causal m (true descendants per ancestor in simulation), naive regression m, and MR-estimated m. Explain the differences quantitatively.
- **Suggested framing for a growth model.** "Homeostatic growth drives networks to a **reverberating / quasicritical** regime whose distance to criticality is set by input strength (and possibly by size and topology)". This is more defensible than "self-organizes to criticality", is consistent with P2 and P4, and turns the observed deviations into the result.
- If the model shows m → 1 only as h → 0 (with finite-size scaling), that is a clean, falsifiable statement that matches Zierenberg 2018 and Fosque 2021, and it is a natural point of comparison.

### Gaps
- I did not verify the exact analytic form used by Zierenberg et al. 2018 (whether it is m = 1 − h/R or an equivalent with refractoriness corrections).
- I did not verify whether Fosque et al. 2021 give an explicit prediction for how much exponents shrink per unit of input.

---

## Q5. What makes a simulation-based criticality paper credible or novel to referees in 2025–2026

### Takeaway
Referees will require:
- the phase transition located independently of power laws (control-parameter sweep, Λ, susceptibility, dynamic range, finite-size scaling);
- matched non-critical null models that are run through the identical pipeline;
- robustness to fitting, binning, threshold and avalanche definition;
- a demonstration that conclusions survive realistic subsampling.

Novelty in 2025–2026 lies less in "growth yields power laws", which is established prior art, and more in quantitative, falsifiable predictions linking the homeostatic or growth mechanism to distance to criticality, input dependence, and developmental trajectories that can be compared with data.

### Cited findings
- **Prior art on growing networks** (all [M]; must be cited and differentiated):
  - Abbott & Rohrkemper 2007, "A simple growth model constructs critical avalanche networks", *Prog. Brain Res.* 165:13–19, doi:10.1016/S0079-6123(06)65002-4.
  - Tetzlaff, Okujeni, Egert, Wörgötter & Butz 2010, "Self-organized criticality in developing neuronal networks", *PLoS Comput. Biol.* 6(12):e1001013, doi:10.1371/journal.pcbi.1001013.
  - Kossio, Goedeke, van den Akker, Ibarz & Memmesheimer 2018, "Growing critical: self-organized criticality in a developing neural system", *Phys. Rev. Lett.* 121, 058301, doi:10.1103/PhysRevLett.121.058301.
  - Zierenberg et al. 2018 PRX [V].
- **Referee-level demand for causal and null controls**: Priesemann & Shriki 2018 [V] ("requires causal interventions"); Destexhe & Touboul 2021 [V] (non-critical controls pass tests); Morrell 2024 [V] (latent variables yield avalanche criticality without fine-tuning).
- **Distance to criticality as the observable**: O'Byrne & Jerbi 2022 [V]; Hengen & Shew 2025 [V]; DCC and D2 [V transcript].
- **Input-dependent operating point**: Zierenberg 2018 [V]; Fosque 2021 [V]; Morales et al. 2023 [V] (separating background from input-driven activity).
- **Comparison with experiment requires the experimental pipeline**: MR estimator and subsampling scaling [M]; electrode-overlap bias [V].
- **Developmental data a growth model could target**: the Hengen/Shew lab describes ferret visual cortex development, with scale-free bursts before eye opening and the later emergence of temporal scale invariance [V transcript; unpublished or in-progress, U]. In vitro maturation of cultures (Tetzlaff 2010 [M]; Pasquale et al. 2008 [M]) shows a developmental progression through sub-, super- and near-critical stages.

### Inferences
Practices most likely to make the paper credible:
1. Treat the model as a statistical-physics object. Sweep the control parameter (global gain, or Λ directly) on frozen grown networks to locate the transition by:
   - divergence of the autocorrelation time τ = −Δt/ln m;
   - peak susceptibility or variance of activity;
   - peak dynamic range;
   - finite-size scaling.

   Then show where homeostatic growth places the network relative to that transition.
2. Report distance to criticality (1 − m, 1 − Λ, DCC) as continuous quantities with confidence intervals across seeds, not binary verdicts.
3. Give a mechanistic explanation of P2 (epoch mixtures) as homeostatic hovering, citing Bonachela et al. 2010. Show conditional statistics.
4. Present the σ versus Λ versus m differences (P4) as a result.

Possible novelty angles:
- input-dependence (h) scaling laws;
- network-size dependence of the homeostatic fixed point;
- topology/heterogeneity effects on Λ versus σ;
- predictions for developmental trajectories (sub- to supercritical overshoot to near-critical) that can be compared with culture or ferret data;
- explicit tests against the latent-variable and time-varying-drive nulls.

### Gaps
- I found no published "reviewer checklist" or community standard document for criticality claims. The checklist below is synthesized from the cited methodological papers.
- I did not verify 2024–2026 simulation-based growth-model papers beyond those listed [M]. A dedicated search for 2023–2026 "growing critical network" models is recommended to confirm novelty.

---

## Q6. Checklist of recommended analyses for the project

### Takeaway
Run every item on the model and on at least two matched non-critical null models, and report all of them, including failures.

### Cited findings (basis for each item)

**A. Locating the transition (model-only gold standard)**
- [ ] Compute **Λ (largest eigenvalue)** of the transmission-probability matrix alongside σ (mean column sum) through growth — [Larremore et al. 2011](https://doi.org/10.1103/PhysRevLett.106.058101) [V].
- [ ] On frozen snapshots of grown networks, rescale weights by a global gain g and sweep g to find g_c. Measure:
  - autocorrelation time;
  - variance/susceptibility;
  - **dynamic range** Δ = 10·log10(h_0.9/h_0.1) — [Kinouchi & Copelli 2006](https://doi.org/10.1038/nphys289) [M];
  - survival probability of single-seed cascades.

  Report where the homeostatic network sits (g/g_c, 1 − Λ).
- [ ] **Finite-size scaling**: repeat for N spanning at least one decade. Exponents should be N-independent, cutoffs should scale as powers of N, and P(S) should collapse under N-rescaling.

**B. Branching ratio and m (P4)**
- [ ] Report the full set of m estimates and explain the gaps:
  - causal m (mean true descendants per active unit, from simulation lineage);
  - naive regression m, with and without external-input subtraction;
  - MR-estimated m on the full population and on random subsamples (1%, 10%) — [Wilting & Priesemann 2018](https://doi.org/10.1038/s41467-018-04725-4) [M]; [MR. Estimator toolbox](https://doi.org/10.1371/journal.pone.0249447) [M].

  Candidate explanations for the gaps: saturation, refractoriness, input, bin width.
- [ ] Check the immigration relation **m ≈ 1 − h/A*** (drive per bin / mean activity per bin) against measured m across input strengths — [Zierenberg et al. 2018](https://doi.org/10.1103/PhysRevX.8.031018) [V paper; M exact form].
- [ ] Report the intrinsic timescale τ = −Δt/ln m and its dependence on h and N.

**C. Avalanche definition and drive (P3)**
- [ ] Compute **causal avalanches** (each external seed's lineage) and **time-binned avalanches**, and compare. Report the overlap fraction (the probability that a new seed arrives while another cascade is active). Caveat: causal avalanches can be scale-free off-criticality — [Martinello et al. 2017](https://doi.org/10.1103/PhysRevX.7.041071) [M].
- [ ] Sweep drive rate h over several decades and report exponents versus h (quasicriticality predicts shrinking exponents with preserved scaling relation) — [Fosque et al. 2021](https://doi.org/10.1103/PhysRevLett.126.098101) [V]; [Das & Levina 2019](https://arxiv.org/abs/1808.04196) [V].
- [ ] **Bin-size sweep**: Δt ∈ {0.25, 0.5, 1, 2, 4}× mean IEI (and × the transmission delay). Plot τ(Δt), α(Δt) and γ(Δt) — [Beggs & Plenz 2003](https://doi.org/10.1523/JNEUROSCI.23-35-11167.2003) [M]; [Priesemann et al. 2014](https://doi.org/10.3389/fnsys.2014.00108) [M].
- [ ] If any continuous signal is thresholded: sweep the threshold, define size as area *above threshold*, and report exponent stability — [Villegas et al. 2019](https://doi.org/10.1103/PhysRevE.100.012133) [V].

**D. Fitting (P1)**
- [ ] Use **doubly truncated power-law MLE** with Deluca–Corral range selection (maximize log(x_max/x_min) subject to a goodness-of-fit p-value > 0.1) instead of KS-minimized x_min — [Deluca & Corral 2013](https://doi.org/10.2478/s11600-013-0154-9) [M]; [NCC toolbox, Marshall et al. 2016](https://doi.org/10.3389/fphys.2016.00250) [M].
- [ ] Show exponent-stability plots versus x_min and x_max, and report the number of decades fitted and the number of avalanches in the range.
- [ ] Run likelihood-ratio tests (Vuong p) against lognormal, exponential, stretched exponential and power law with exponential cutoff — [Clauset et al. 2009](https://doi.org/10.1137/070710111) [M].
- [ ] Report curvature honestly. A fixed-window MLE that reveals curvature is evidence of distance from criticality and should be quantified, e.g., by local slope d log P / d log S versus S.

**E. Scaling relations and shape collapse**
- [ ] Fit τ, α and γ_fit (from ⟨S⟩(T)) independently. Report γ_pred = (α−1)/(τ−1), **DCC** = |γ_pred − γ_fit|, and shape-collapse γ_sc, all with bootstrap confidence intervals — [Friedman et al. 2012](https://doi.org/10.1103/PhysRevLett.108.208102) [M]; [Sethna et al. 2001](https://doi.org/10.1038/35065675) [M]; [Ma et al. 2019](https://doi.org/10.1016/j.neuron.2019.08.031) [M].
- [ ] Compare with the expected universality class (MF-DP: τ = 3/2, α = 2, γ = 2 for sparse random networks). Explain deviations with subsampling, input or quasicriticality — [Carvalho et al. 2021](https://doi.org/10.3389/fncir.2020.576727) [M].

**F. Non-stationarity and mixtures (P2)**
- [ ] Use sliding-window m (MR or causal) to obtain the time series and distribution of m(t). Report the fraction of time with m < 0.95, 0.95–1.0, and > 1.0.
- [ ] Compute **conditional avalanche distributions** binned by m(t), and show how pooled power laws arise from the mixture — [Priesemann & Shriki 2018](https://doi.org/10.1371/journal.pcbi.1006081) [V]; [Bonachela et al. 2010](https://doi.org/10.1088/1742-5468/2010/02/P02015) [M].
- [ ] Analyse only post-convergence epochs, with an explicit stationarity test, separately from growth or transient epochs.

**G. Null models run through the identical pipeline**
- [ ] **Inhomogeneous Poisson** units driven by the model's own population-rate trajectory (no recurrence) — [Priesemann & Shriki 2018](https://doi.org/10.1371/journal.pcbi.1006081) [V].
- [ ] **Latent-variable model**: independent units coupled to 1–5 slow latent fields matched for rate and pairwise correlations — [Morrell et al. 2024](https://doi.org/10.7554/eLife.89337) [V].
- [ ] A **known non-critical network**: a Brunel E/I network in the synchronous-irregular regime, or the grown network with gain rescaled to Λ = 0.9 and Λ = 1.1 — [Destexhe & Touboul 2021](https://doi.org/10.1523/ENEURO.0551-20.2021) [V].
- [ ] **Shuffled controls**: degree-preserving rewiring (tests the role of topology) and spike-time jitter/shuffling (tests temporal structure).

**H. Realistic observation (for comparison with experiment)**
- [ ] **Subsample** at random (p = 0.01–0.5). Apply **subsampling scaling** collapse and the MR estimator — [Levina & Priesemann 2017](https://doi.org/10.1038/ncomms15140) [M]; [Priesemann et al. 2009](https://doi.org/10.1186/1471-2202-10-40) [M].
- [ ] If coarse "electrode" signals are simulated, include and vary the spatial overlap — [Neto et al. 2022](https://doi.org/10.1371/journal.pcbi.1010678) [V].

**I. Function and framing**
- [ ] Report dynamic range and an intrinsic timescale versus distance to criticality across growth, input strength and N.
- [ ] Frame results as reverberating / quasicritical / near-critical with a measured distance, unless finite-size scaling plus Λ → 1 plus all of A–G support criticality — [Wilting & Priesemann 2019](https://doi.org/10.1016/j.conb.2019.08.002) [V]; [O'Byrne & Jerbi 2022](https://doi.org/10.1016/j.tins.2022.08.007) [V]; [Hengen & Shew 2025](https://doi.org/10.1016/j.neuron.2025.05.020) [V].
- [ ] Cite and differentiate the growth-model prior art (Abbott & Rohrkemper 2007; Tetzlaff et al. 2010; Kossio et al. 2018; Zierenberg et al. 2018) [M except Zierenberg V].
- [ ] Robustness: at least 10 seeds per condition, bootstrap confidence intervals on all exponents and m, and a report of parameter ranges where conclusions fail.

### Inferences
- Items A, C (causal avalanches) and G are where simulation papers can exceed experimental ones. Ground-truth lineage, exact Λ and finite-size scaling are unavailable in experiments, and referees in 2025–2026 increasingly expect them to be used.
- Items F and B(ii) convert the project's existing pitfalls (P2, P4) into results. Hovering homeostasis and an input-set distance to criticality are publishable findings consistent with the current consensus.

### Gaps
- No single authoritative community checklist exists. This list is a synthesis.
- Items marked [M] in citations should be DOI-checked before inclusion in a manuscript.
