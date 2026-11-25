# Prior work: activity-dependent / homeostatic neurite growth and structural plasticity driving networks to (self-organized) criticality

Tags: **[V]** = verified by reading the full text or abstract (in this session). **[V-rev]** = verified through a review written by the same group (van Ooyen & Butz-Ostendorf 2019), not the original paper. **[U]** = uncertain: from memory or a search snippet only, not checked against the paper.

## Q1. What exactly did the core papers show?

### Takeaway
The "neuritic-field" lineage (van Ooyen & van Pelt 1994/96 -> Abbott & Rohrkemper 2007 -> Tetzlaff et al. 2010 -> Kalle Kossio et al. 2018) already established the following. Homeostatic growth or retraction of circular neurite fields, with coupling proportional to overlap area, drives an excitatory network to a near-critical branching state with s^-3/2 sizes and T^-2 durations. Kossio 2018 derived **sigma = 1 - f0/f_sat** analytically. Separately, Zierenberg et al. 2018 derived **r* = h/(1-m)**, i.e. m = 1 - h/r*, for synaptic-scaling homeostasis in a branching network. So the project's relation m ≈ 1 - h/r* is a known result in both lineages.

### Cited Findings

**van Ooyen & van Pelt 1994 / 1996 (original neuritic-field model)**
- Citations: van Ooyen A, van Pelt J (1994) "Activity-dependent outgrowth of neurons and overshoot phenomena in developing neural networks", J Theor Biol 167:27–43. van Ooyen A, van Pelt J (1996) "Complex periodic behaviour in a neural network model with activity-dependent neurite outgrowth", J Theor Biol 179:229–242. **[V]** (citation details from the reference list of Abbott & Rohrkemper 2007) — [Abbott & Rohrkemper PDF](https://zpkilpat.github.io/appm4370/papers/abbott07.pdf)
- Model, as described in the 2019 review. Each neuron has one circular neuritic field that represents both axon and dendrite. Coupling is proportional to the overlap area. The field expands when activity is below the set-point and shrinks when it is above. Activity follows a firing-rate "shunting" model. There is **no external input**: cultures self-assemble from activity generated inside the network. **[V-rev]** — [van Ooyen & Butz-Ostendorf 2019 (PDF)](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- Overshoot: networks pass "via a transient phase of high connectivity (overshoot) into a stable network with network-wide homeostasis". Neurons acquire small fields where cell density is high and large fields where it is low. **[V-rev]** — [van Ooyen & Butz-Ostendorf 2019](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- Oscillation: depending on where the activity–connectivity curves intersect, "connectivity and activity will oscillate on the time scale of growth" (ref. 71 = the 1996 paper). Otherwise growth converges without overshoot. This is a deterministic, growth-timescale oscillation in a rate model. It is not criticality-related. **[V-rev]** — [van Ooyen & Butz-Ostendorf 2019](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- A 1995 paper (van Ooyen, van Pelt & Corner, J Theor Biol 172:63–82, "Implications of activity dependent neurite outgrowth for neuronal morphology and network development") **[U]**: from memory only, not checked.

**Abbott & Rohrkemper 2007**
- Citation: Abbott LF, Rohrkemper R (2007) "A simple growth model constructs critical avalanche networks", Prog Brain Res 165:13–19 (in Cisek, Drew & Kalaska, eds., *Computational Neuroscience: Theoretical Insights into Brain Function*). DOI [10.1016/S0079-6123(06)65002-4](https://doi.org/10.1016/S0079-6123(06)65002-4). **[V]** (full text read) — [PDF](https://zpkilpat.github.io/appm4370/papers/abbott07.pdf)
- Model:
  - N = 100 neurons placed at random in a unit square. Each has one circle of radius a_i.
  - Coupling: r_i -> r_i + g·A_ij at each spike of j, where A_ij is the overlap area and g = 500 Hz. All connections are excitatory.
  - Poisson spiking: the rate relaxes to r0 = 0.1 Hz with tau_r = 5 ms. Refractory period is 20 ms. Time step dt = 1 ms.
  - Calcium C_i is incremented at each spike and decays with tau_C = 100 ms.
  - Growth rule: da_i/dt = k(C_target - C_i), with k = 0.002 s^-1 near equilibrium, which is "slow on the time scale of neuronal activity".
  - The authors say their only change to van Ooyen–van Pelt was adding Poisson spikes.

  **[V]** — [PDF](https://zpkilpat.github.io/appm4370/papers/abbott07.pdf)
- Result: at equilibrium (growth switched off while collecting statistics, which "did not make any noticeable difference"), spike-count and duration histograms follow -3/2 and -2 power laws. They note "deviations for large, rare events due to finite-size effects". There was no systematic finite-size scaling. **[V]** — [PDF](https://zpkilpat.github.io/appm4370/papers/abbott07.pdf)
- Criticality is **not generic** in their model. The size exponent depends on C_target (Fig. 4; range 0.04–1.2), and -1.5 is hit only at C_target ≈ 0.08. Higher C_target gives flat (supercritical-like) distributions; lower gives a shortage of large events. Their explanation: the spontaneous rate of 0.1 Hz cannot reach the target alone, so neurons must "borrow activity" from others. **[V]** — [PDF](https://zpkilpat.github.io/appm4370/papers/abbott07.pdf)
- They explicitly suggested future work with separate dendritic and axonal circles and distal axonal projections. They do not discuss realistic morphology. **[V]** — [PDF](https://zpkilpat.github.io/appm4370/papers/abbott07.pdf)

**Tetzlaff, Okujeni, Egert, Wörgötter & Butz 2010**
- Citation: Tetzlaff C, Okujeni S, Egert U, Wörgötter F, Butz M (2010) "Self-organized criticality in developing neuronal networks", PLoS Comput Biol 6(12): e1001013. DOI [10.1371/journal.pcbi.1001013](https://doi.org/10.1371/journal.pcbi.1001013). **[V]** (citation via the Zeraati 2021 review) — [Zeraati et al. 2021](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2021.619661/full)
- Model:
  - Two circular fields per neuron. The dendritic field follows dR_den/dt = rho_den(C_target - C). The axonal field follows the **opposite** rule, dR_axo/dt = rho_axo(C - C_target), because axons need activity to grow.
  - Stochastic threshold neurons, with a membrane potential bounded at 1.
  - About 20% of neurons are switched to inhibitory mid-development (the GABA switch).

  **[V-rev]** — [van Ooyen & Butz-Ostendorf 2019](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf); [Zeraati et al. 2021](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2021.619661/full)
- Three developmental phases:
  - Phase I: slow rise; avalanches go from Poisson-like to power law.
  - Phase II: abrupt connectivity and activity peak; the distribution is **supercritical**. It stays supercritical even with inhibitory strength 100× the excitatory strength.
  - Phase III: pruning to homeostasis; the distribution is **critical**, with a size exponent "close to -1.5" (when |w_inh| = |w_exc|).

  **[V-rev]** — [van Ooyen & Butz-Ostendorf 2019](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- Review summary of Tetzlaff: "criticality is best reached with 20% inhibitory cells". Abbott and Kossio show inhibition is not required. **[V-rev]** — [same](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- The paper also compared against developing cortical cultures on MEAs (Okujeni/Egert). The PMC abstract snippet says avalanche distributions are sensitive to inhibition and "very unstable during development". **[U]** (snippet only) — [PMC2996321](https://pmc.ncbi.nlm.nih.gov/articles/PMC2996321/)

**Butz & van Ooyen: Model of Structural Plasticity (MSP)**
- Review: Butz M, Wörgötter F, van Ooyen A (2009) "Activity-dependent structural plasticity", Brain Res Rev 60(2):287–305 **[U]** (citation from memory; DOI likely 10.1016/j.brainresrev.2008.12.023, not checked).
- MSP: Butz M, van Ooyen A (2013) "A simple rule for dendritic spine and axonal bouton formation can account for cortical reorganization after focal retinal lesions", PLoS Comput Biol 9(10): e1003259 **[U]** (citation from memory).
- MSP is not a geometric overlap model. Each neuron has counts of axonal synaptic elements (boutons) and dendritic elements (spines), which "can combine to form synapses". Element numbers grow or retract with calcium (Gaussian growth curves in later versions **[U]**). MSP was implemented in NEST. It has been used for lesion reorganization, small-world emergence and adult neurogenesis, but **"has not yet been used to study avalanche dynamics"** as of 2019. **[V-rev]** — [van Ooyen & Butz-Ostendorf 2019](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- Diaz-Pier S, Naveau M, Butz-Ostendorf M, Morrison A (2016) "Automatic generation of connectivity for large-scale neuronal network models through structural plasticity", Front Neuroanat 10:57 (MSP in NEST, used to reach target rates). **[U]** (from memory; not fetched)

**Kalle Kossio, Goedeke, van den Akker, Ibarz & Memmesheimer 2018**
- Citation: Kalle Kossio FY, Goedeke S, van den Akker B, Ibarz B, Memmesheimer R-M (2018) "Growing critical: self-organized criticality in a developing neural system", Phys Rev Lett 121, 058301. DOI [10.1103/PhysRevLett.121.058301](https://doi.org/10.1103/PhysRevLett.121.058301); arXiv:1811.02861. **[V]** (arXiv full text read; journal ref inferred from the supplement DOI in the paper) — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- Model:
  - Neurons are a linear multivariate Hawkes process. Each spike raises postsynaptic rates by an exponentially decaying kernel (tau).
  - Somata are placed uniformly at random in a plane, with disk neurite radius R_i.
  - Coupling = g·A_ij (overlap area), excitatory only.
  - Spontaneous rate f0 per neuron.
  - Growth: dR_i/dt = K(1 - (1/f_sat) Σ δ(t - t_i^spike)). The radius grows linearly and shrinks by K/f_sat per emitted spike, so the mean rate is driven to f_sat.

  **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- Key analytical result: the branching parameter at the grown steady state is **sigma = tau·g·Ā = 1 - f0/f_sat**. The model is near-critical when f0 ≪ f_sat and subcritical when f_sat is small. **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- Exponents:
  - Avalanche sizes (defined as Hawkes clusters, i.e. genealogical trees, not time-binned) follow a Borel distribution, approximately s^-3/2 with cutoff s_c = (sigma - ln sigma - 1)^-1.
  - Durations have a tail exponent of 2, from an ODE for the generating function plus closed-form approximations.

  **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- Simulations used N = 100, with F_target = 0.04 Hz (subcritical) vs 2 Hz (near-critical) and F0 = 0.01 Hz. **[V-rev]** — [van Ooyen & Butz-Ostendorf 2019, Fig. 4](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- Robustness:
  - Fluctuations of Σ_j A_ij and of rate deviations from f_sat were small (<1%).
  - Freezing growth (K = 0) after convergence leaves the statistics unchanged ("growth carries the system close to a critical point, but is not required later on").
  - Shuffling coupling columns barely matters.
  - Refractoriness with tau_ref ≈ tau is fine; longer refractoriness distorts the statistics.

  **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)

**Zierenberg, Wilting & Priesemann 2018**
- Citation: Zierenberg J, Wilting J, Priesemann V (2018) "Homeostatic plasticity and external input shape neural network dynamics", Phys Rev X 8, 031018. DOI [10.1103/PhysRevX.8.031018](https://doi.org/10.1103/PhysRevX.8.031018); arXiv:1807.01479. **[V]** — [PRX](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.031018)
- Model:
  - Branching network of N excitatory binary stochastic neurons in discrete time (Δt = 1 ms), with external Poisson input h per neuron.
  - Homeostasis is a linear negative feedback on the postsynaptic scaling factor (activation probability) α_j: Δα_j = (Δt·r* - s_j,t)(Δt/τ_hp). This is **synaptic scaling, not structural growth**.
  - Tested on several topologies (e.g. ER), with r* = 1 Hz and τ_hp ~ 10^3 s.

  **[V]** (full PDF text) — [arXiv PDF](https://arxiv.org/pdf/1807.01479)
- Key relation: the stationary rate is **r = h/(1 - m)** (their Eq. 4). Homeostasis therefore sets m = 1 - h/r*, the same relation the project found. As the input rate falls, the state moves from irregular (m ≈ 0) to fluctuating/reverberating (m → 1) to bursting (m > 1 transiently). **[V]** — [arXiv PDF](https://arxiv.org/pdf/1807.01479)
- Breakdown at low drive: the mean-field (stationary branching) description holds only while τ ≪ τ0 = τ_hp/k, where k is the number of outgoing connections. The fluctuating-to-bursting transition occurs at **h/r* ≈ Δt/τ0**. Below that, the homeostatic feedback can no longer track activity, and the dynamics become bursty with peaked avalanche-size distributions. This is a drive-dependent limit on homeostatic criticality, analogous to the project's quasi-critical oscillation. **[V]** — [arXiv PDF](https://arxiv.org/pdf/1807.01479)
- The Zeraati review describes this regime as the network "hover[ing] between sub- and supercritical dynamics" at very weak input. **[V]** — [Zeraati et al. 2021](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2021.619661/full)
- Matched to data: in vitro bursting corresponds to input below 0.1% of spikes, and in vivo to about 1% input-driven spikes. **[V]** — [PRX abstract](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.031018)

### Inferences
- The project's **m ≈ 1 - h/r*** is mathematically the same as Kossio's sigma = 1 - f0/f_sat and Zierenberg's r* = h/(1 - m). It is a consequence of rate homeostasis in any linear branching/Hawkes process and does not depend on the growth rule or morphology. It should be presented as a reproduction or consistency check, not a finding.
- Scaling h ∝ 1/N at fixed r* gives m -> 1 as N -> ∞, again a direct corollary of the same formula. The non-trivial content would be how the finite-size cutoff, the timescale separation (total drive N·h vs. avalanche duration) and the coupling fluctuations scale with N. Kossio claims their size distribution is independent of N because their avalanches are defined genealogically. With time-binned avalanches, overlap from the total drive N·h matters.
- "Coupling oscillates / quasi-critical" has precedents in three distinct mechanisms. (i) van Ooyen & van Pelt's deterministic growth-timescale oscillations. (ii) Zierenberg's low-h bursting/hovering when τ_hp/k is too short relative to the avalanche time. (iii) The self-organized quasi-criticality (SOqC) literature (Q2). The project should identify which mechanism it is seeing.

### Gaps
- I did not read the Tetzlaff 2010 or van Ooyen & van Pelt original full texts. Their details come from the 2019 review by van Ooyen & Butz-Ostendorf, i.e. the same group.
- I could not open the MSP primary papers (Butz & van Ooyen 2013; Butz, Wörgötter & van Ooyen 2009) or Diaz-Pier 2016. Their citations come from memory **[U]**.
- I did not confirm whether Kossio et al. analysed the dependence on growth rate K quantitatively. The full text states the condition 1/K ≫ tau and shows that freezing growth does not change the statistics.

## Q2. Which 2019–2026 papers build on these?

### Takeaway
I found few direct follow-ups on growth-driven criticality after 2018. The main ones are a 2019 review chapter (van Ooyen & Butz-Ostendorf), the 2021 Zeraati/Priesemann/Levina review, and Landmann/Baumgarten/Bornholdt 2021 (topological rewiring, not geometric growth). Gallinaro & Rotter's work uses MSP-like homeostatic structural plasticity for memory, not criticality. A parallel SOqC literature (Kinouchi, Girardi-Schappo, Copelli, Muñoz) gives the theory for hovering or oscillation around criticality under homeostasis in finite networks. I found no 2022–2026 paper combining neurite growth with criticality.

### Cited Findings
- **van Ooyen A, Butz-Ostendorf M (2019)**, "Homeostatic structural plasticity can build critical networks". In: Tomen N, Herrmann JM, Ernst U (eds) *The Functional Role of Critical Dynamics in Neural Systems*, Springer Series on Bio- and Neurosystems 11, pp. 117–137. DOI [10.1007/978-3-030-20965-0_7](https://doi.org/10.1007/978-3-030-20965-0_7). **[V]** — [author PDF](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
  - It reviews van Ooyen–van Pelt, Abbott & Rohrkemper, Kossio and Tetzlaff.
  - It states: "The use of circular neuritic fields in all models is a simple yet powerful way to abstract away from detailed neuronal morphology", with the drawback that the strongest connections form between neighbours.
  - Open problems listed: MSP not yet used for avalanches; analytical role of inhibition; heterogeneous set-points.
- **Zeraati R, Priesemann V, Levina A (2021)**, "Self-organization toward criticality by synaptic plasticity", Front Phys 9:619661. DOI [10.3389/fphy.2021.619661](https://doi.org/10.3389/fphy.2021.619661); arXiv:2010.07888. **[V]** — [Frontiers](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2021.619661/full)
  - Its structural-plasticity section summarises growth models: fields expand below target and shrink above it; coupling is proportional to overlap.
  - It cites Kossio as showing analytically that a target rate well above the spontaneous rate gives criticality, while a small target gives subcritical dynamics.
  - It cites Tetzlaff for the separate axon/dendrite rules and the GABA switch.
  - It calls Hebbian-inspired rewiring rules toward criticality "less biologically plausible".
  - It stresses that external input was "overlooked for a while" as a control parameter, and that power laws alone do not prove criticality.
- **Landmann S, Baumgarten L, Bornholdt S (2021)**, "Self-organized criticality in neural networks from activity-based rewiring", Phys Rev E 103, 032304. DOI [10.1103/PhysRevE.103.032304](https://doi.org/10.1103/PhysRevE.103.032304); arXiv:2009.11781. **[V]** (abstract via PubMed/ADS) — [PubMed](https://pubmed.ncbi.nlm.nih.gov/33862737/)
  - Rewiring depends on each node's average activity, and the network self-organizes to a critical state with power-law avalanches. The authors suggest it as a developmental mechanism.
  - Uncertain details: connectivity is abstract (link addition/removal, descended from Bornholdt & Rohlf 2000 adaptive Boolean networks), not spatial or morphological **[U]**.
- **Gallinaro JV, Gašparović N, Rotter S (2022)**, "Homeostatic control of synaptic rewiring in recurrent networks induces the formation of stable memory engrams", PLoS Comput Biol 18(2): e1009836. bioRxiv 10.1101/2020.03.27.011171. **[V]** (citation and abstract-level claim) — [bioRxiv](https://www.biorxiv.org/content/10.1101/2020.03.27.011171v2.full); [AGRIS](https://agris.fao.org/search/en/providers/122436/records/675ab3bb0ce2cede71cdf9b3)
  - Homeostatic structural plasticity (MSP-style synaptic elements) alone forms "silent" memory engrams.
  - It is **not about criticality**.
  - Precursor: Gallinaro & Rotter (2018) "Associative properties of structural plasticity based on firing rate homeostasis in recurrent neuronal networks", Sci Rep 8:3754 **[U]** (from memory).
  - A CNS*2024 abstract notes the 2022 model allowed only single-engram formation **[V]** (snippet) — [CNS*2024](https://link.springer.com/article/10.1007/s10827-024-00889-9)
- **Girardi-Schappo M, Brochini L, Costa AA, Carvalho TTA, Kinouchi O (2020)**, "Synaptic balance due to homeostatically self-organized quasicritical dynamics", Phys Rev Research 2, 012042. DOI [10.1103/PhysRevResearch.2.012042](https://doi.org/10.1103/PhysRevResearch.2.012042). **[V]** (snippet/abstract; author list partly from memory **[U]**) — [arXiv](https://arxiv.org/abs/2002.09117)
  - Homeostatic adaptation of synapses and thresholds tunes the network imperfectly, so it **hovers** around the critical point. The oscillation amplitude shrinks as N grows.
- **Kinouchi O, Brochini L, Costa AA, Campos JGF, Copelli M (2019)**, "Stochastic oscillations and dragon king avalanches in self-organized quasi-critical systems", Sci Rep 9:3874. DOI [10.1038/s41598-019-40473-1](https://doi.org/10.1038/s41598-019-40473-1). **[V]** (abstract)
  - In finite SOqC networks with dynamic gains or depressing synapses, demographic noise excites an almost-unstable focus. The result is stochastic oscillations of the control parameter around criticality, plus large "dragon king" avalanches.
  - This is the closest theory to the project's "coupling oscillates" observation. — [Sci Rep](https://www.nature.com/articles/s41598-019-40473-1)
- **Kinouchi O, Pazzini R, Copelli M (2020)**, "Mechanisms of self-organized quasicriticality in neuronal network models", Front Phys 8:583213. DOI 10.3389/fphy.2020.583213; arXiv:2010.10622. **[U]** (title and DOI from search results; content not read) — [arXiv](https://arxiv.org/pdf/2010.10622)
- A 2022 "Homeostatic criticality in neuronal networks" paper (Chaos, Solitons & Fractals, S0960077922000881; Girardi-Schappo/Kinouchi group) extends homeostatic quasi-criticality to the presence of external fields. **[U]** (snippet; authors and volume not checked) — [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S0960077922000881)

### Inferences
- The Memmesheimer group's later work (e.g. Kalle Kossio et al. on drifting assemblies, from memory **[U]**) moved away from growth-to-criticality. I found no direct sequel extending Kossio 2018 to non-disk morphologies, finite-size scaling or drive scaling.
- The space between geometric growth (circles) and abstract rewiring (Bornholdt, MSP element counts) is still largely open. So are finite-size and drive scaling of growth-driven criticality.

### Gaps
- Google Scholar "cited by" lists for Kossio 2018 and Abbott & Rohrkemper 2007 could not be accessed (no tool access). The 2022–2026 coverage is therefore based on keyword searches and may miss papers, especially conference abstracts, theses and non-arXiv journals.
- I found no 2023–2026 paper on growth plus criticality in these searches. This is absence of evidence, not proof.
- One search hit "Adaptive self-organized criticality in deep neural networks" (arXiv:2608.28431) concerns artificial networks. Its ID and date were not checked, and it is off-topic.

## Q3. Has anyone combined realistic or fractal dendritic morphology with homeostatic growth and criticality?

### Takeaway
I found no paper that does this. Every growth-to-criticality model uses circular or disk neuritic fields (van Ooyen, Abbott, Kossio; Tetzlaff uses separate axon and dendrite disks) or non-spatial synaptic-element counts (MSP). The 2019 review by the field's own authors names circular fields as an abstraction of morphology. Fractal or stochastic branching morphology, dendrite–dendrite/soma contact connectivity, and homeostatic growth with retraction toward criticality appear novel as a combination. This rests on the searches reported here, not an exhaustive citation crawl.

### Cited Findings
- "The use of circular neuritic fields in all models is a simple yet powerful way to abstract away from detailed neuronal morphology. A disadvantage is that it puts some constraints on the type of network topologies that can arise, as the strongest connections are usually formed between neighboring cells." **[V]** — [van Ooyen & Butz-Ostendorf 2019](https://anvooyen.home.xs4all.nl/papers/vanooyen_butz_cdns_2019.pdf)
- Abbott & Rohrkemper suggested separate dendritic and axonal circles and distal axonal circles as extensions, still circles. **[V]** — [PDF](https://zpkilpat.github.io/appm4370/papers/abbott07.pdf)
- Tangential prior work with dendrites and criticality: Gollo LL, Kinouchi O, Copelli M (2013), "Single-neuron criticality optimizes analog dendritic computation", Sci Rep 3:3222. The dendritic arbor is an excitable Cayley tree, and criticality is within one neuron. There is no growth or network homeostasis. **[V]** (title and abstract snippet) — [Sci Rep](https://www.nature.com/articles/srep03222)
- Morphologically detailed growth simulators exist without a criticality analysis, e.g. NETMORPH, CX3D/BioDynaMo, and Okujeni & Egert 2019 (eLife) on activity-dependent migration and outgrowth producing modular cultures. **[U]**: from memory, not checked in this session.

### Inferences
- What the project could claim as novel:
  1. Fractal stochastic dendritic morphology with contact-based connectivity under rate homeostasis.
  2. Showing that the 1 - h/r* attractor survives non-disk morphology and dendrite–dendrite coupling. This is expected from the general argument, so it is weak on its own.
  3. Systematic finite-size and drive scaling (h ∝ 1/N) of a growth-driven network, which none of the growth papers did; Abbott and Kossio used N = 100.
  4. Characterising the quasi-critical coupling oscillation and linking it to SOqC theory (Kinouchi 2019; Girardi-Schappo 2020) and to the Zierenberg low-h breakdown.
- Reviewers will likely point out that m = 1 - h/r* and "criticality when h/r* -> 0" are already in Kossio 2018 and Zierenberg 2018.

### Gaps
- I could not do an exhaustive search of morphology-detailed network growth models (NETMORPH, CX3D, BioDynaMo, MSP-with-geometry) for any avalanche analysis.

## Q4. What did Kossio 2018 conclude about the limits needed for criticality?

### Takeaway
Kossio et al. require (i) **low spontaneous drive relative to the homeostatic target, f0 ≪ f_sat**, because sigma = 1 - f0/f_sat; (ii) **slow growth**, 1/K ≫ tau; and (iii) approximately linear spike interactions. They state that the avalanche-size distribution is **independent of N**, so large N is not their limit. Subcriticality is attributed to finite f0/f_sat, not to finite size.

### Cited Findings
- sigma = tau·g·Ā = 1 - f0/f_sat. Criticality needs f0 ≪ f_sat, and small f_sat gives subcritical dynamics. "Our model predicts that neural networks may develop criticality already due to their growth." **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- Growth must be slow relative to activity decay (1/K ≫ tau). Once grown, growth can be frozen without changing the statistics. **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- "The size distribution is independent of neuron number." Simulations used N = 100. **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- Deviations come from strong nonlinearities and long refractory periods. tau_ref ≈ tau leaves small and medium avalanches nearly unchanged. **[V]** — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)
- Zeraati et al. paraphrase: if the target rate is "significantly larger than the spontaneous firing rate", growth brings the network to a critical point. **[V]** — [Zeraati et al. 2021](https://www.frontiersin.org/journals/physics/articles/10.3389/fphy.2021.619661/full)
- Contrast with Abbott & Rohrkemper: there, the exponent depended on C_target and -3/2 was reached only near one value (0.08). Kossio's formula explains this as sigma approaching 1 only as f0/f_sat -> 0. Abbott's flat distributions at high C_target likely reflect refractoriness and nonlinear saturation. That last point is my inference.

### Inferences
- Kossio's N-independence holds because avalanches are defined genealogically (Hawkes clusters), which separates overlapping cascades. With time-binned avalanches, as in the project and in Abbott, the total drive N·f0 sets an extra N-dependent limit: avalanches merge unless N·h·(avalanche duration) ≪ 1. Scaling h ∝ 1/N is exactly what keeps this separation, so the project's large-N result is the time-binned analogue of Kossio's f0/f_sat -> 0 limit. This inference is not stated in Kossio.

### Gaps
- I did not extract whether Kossio quantified the coupling fluctuations' dependence on K or N beyond the "<1%" statement.
