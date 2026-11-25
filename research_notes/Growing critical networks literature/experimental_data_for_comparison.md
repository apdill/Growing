# Experimental data for comparison: how avalanche statistics develop as networks grow, and homeostatic regulation of criticality

Verification legend: **[V]** = checked this session against an abstract, full text, or publisher/indexing record. **[V-abs]** = checked against the abstract only. **[U]** = from background knowledge or secondary summaries, not checked this session. Treat [U] items as leads to confirm, not as facts.

---

## Q1. Developing dissociated cultures on MEAs: how do avalanche statistics change over days in vitro (DIV)?

### Takeaway
Every longitudinal culture study found that avalanche statistics are not critical from the start. They pass through non-critical stages before reaching power-law-like (near-critical) statistics. The path differs by study. Tetzlaff et al. 2010 report initial, then supercritical (~20–22 DIV), then subcritical (~36 DIV), then critical (~58 DIV) on 59-electrode MEAs. Yada et al. 2017 report exponential, then bimodal ("supercritical" with a large-scale avalanche), then power-law-like within 4–30 DIV on HD-CMOS MEAs. Organotypic slices (Stewart & Plenz 2008) reach a stable exponent near −1.5 much sooner. The best quantitative target for a growth simulation is the qualitative sequence: sub/exponential → supercritical bump or bimodal (bursting) → power law. Timings depend on preparation and density, so they make a weak target.

### Cited Findings

**Tetzlaff, Okujeni, Egert, Wörgötter & Butz 2010** — "Self-organized criticality in developing neuronal networks", *PLoS Comput Biol* 6(12): e1001013. doi:10.1371/journal.pcbi.1001013 — [Article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001013) **[V]**
- Neonatal Wistar rat cortical cultures, ~500,000 cells per culture, on PEI-coated MEAs with 59 TiN electrodes (200 or 500 µm pitch, Multi Channel Systems, MEA1060-BC amplifier). Recordings lasted 1 h each, with a 200 Hz high-pass filter and a 3 SD spike threshold. Cultures needed ≥50 active electrodes to be included. — [PLoS article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001013) [V]
- n = 20 cultures recorded between 13 and 95 DIV. Few were tracked longitudinally: 14% had 5 time points, 7% had 4, 29% had 3, 14% had 2, and 36% had only 1. The developmental sequence is therefore partly cross-sectional. — [PLoS article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001013) [V]
- Four states with average DIV:
  - **initial**, ~19 DIV: low activity, distributions approach a power law but are unstable
  - **supercritical**, ~20 DIV (abstract) or ~22 DIV (results): a pronounced "bump" of large avalanches and significantly higher mean activity
  - **subcritical**, ~36 DIV: large avalanches depressed
  - **critical**, ~58 DIV: close power-law fit
  - State differences were significant by ANOVA. — [PLoS article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001013) [V]
- The authors paired these data with a homeostatic structural-plasticity model (Butz-type synapse growth and retraction driven by a calcium or activity set-point). It reproduced the same sequence: a connectivity overshoot followed by pruning. Their extra tests of scale-freeness were a scaling function, a Fano-factor power law, and scale-free inter-avalanche intervals. — [PLoS article](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001013) [V]
- Experimental exponents and time bins for each state are in their Table 4. The math markup did not render in my text extraction, so I could not read the numbers (see Gaps). The text says the critical-state exponent in cultures and model "matches previous results" of about −1.5. — [PMC2996321](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996321/) [V, value not extracted]
- PKC-inhibitor-treated and untreated cultures did not differ significantly in exponent, number of active electrodes, or activity (their Table 5). This is relevant because PKC manipulation changes neurite aggregation (Okujeni et al., below). — [PMC2996321](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996321/) [V]

**Yada, Mita, Sanada, Yano, Kanzaki, Bakkum, Hierlemann & Takahashi 2017** — "Development of neural population activity toward self-organized criticality", *Neuroscience* 343: 55–65 (online 30 Nov 2016). doi:10.1016/j.neuroscience.2016.11.031 **[V-abs]** (note: the author list is Yada, Mita, Sanada, … Takahashi, not "Nakajima")
- Nine dissociated rat cortical cultures were recorded from 4 to 30 DIV at single-cell resolution on high-density CMOS MEAs. — [Europe PMC record](https://europepmc.org/search?query=DOI%3A%2210.1016%2Fj.neuroscience.2016.11.031%22) [V-abs]
- Avalanche size distribution sequence:
  1. **exponential**, while cells are immature and activity is random-like
  2. **bimodal**, with one homogeneous large-scale avalanche (network burst)
  3. **power-law-like**, with diverse spatiotemporal patterns
  - The authors call this an "integration then fragmentation" two-step process and link it to axonal elongation and synaptogenesis followed by pruning and GABAergic maturation. — [Europe PMC record](https://europepmc.org/search?query=DOI%3A%2210.1016%2Fj.neuroscience.2016.11.031%22) [V-abs]
- Exact DIV boundaries between stages and the exponents were not in the abstract. [Gap]

**Pasquale, Massobrio, Bologna, Chiappalone & Martinoia 2008** — "Self-organization and neuronal avalanches in networks of dissociated cortical neurons", *Neuroscience* 153(4): 1354–1369. doi:10.1016/j.neuroscience.2008.03.050, PMID 18448256 **[V-abs]**
- Embryonic rat cortical cultures on MEAs. Activity goes from isolated spikes in the first days to highly synchronized bursts after 3–4 weeks. Sub-critical, critical, and super-critical size and duration distributions were all observed, "depending on both the age and the development of cultures". — [Europe PMC / PubMed 18448256](https://www.pubmed.gov/18448256/) [V-abs]
- Criticality coincided with **medium synchronization among bursts and a high ratio of bursting to random spiking**. Acetylcholine and bicuculline shifted the distributions consistently with this. The result depends on choosing "an appropriate time scale of observation", meaning the bin width. — [PubMed 18448256](https://www.pubmed.gov/18448256/) [V-abs]
- Specific exponents and DIV ranges were not extracted. [U] From memory, they used 60-electrode MCS MEAs, and the size exponent at criticality was near −1.5 when the bin equals the mean inter-event interval (Beggs & Plenz convention). Unverified.

**Wagenaar, Pine & Potter 2006** — "An extremely rich repertoire of bursting patterns during the development of cortical cultures", *BMC Neurosci* 7: 11. doi:10.1186/1471-2202-7-11 **[V-abs]**
- 58 rat cortical cultures at densities from 3,000 to 50,000 neurons on 30–75 mm², followed on MEAs over the first five weeks. 963 half-hour recordings. — [Europe PMC abstract](https://europepmc.org/search?query=DOI%3A%2210.1186%2F1471-2202-7-11%22) [V-abs]
- Dense cultures started bursting earlier. Stimulus responses suggest faster axonal outgrowth in dense cultures. After 2 weeks, population bursts dominated most cultures. Inter-burst intervals ranged from 1 to 300 s, and burst patterns kept changing through week 5. All but the sparsest cultures showed globally synchronized bursts. — [abstract](https://europepmc.org/search?query=DOI%3A%2210.1186%2F1471-2202-7-11%22) [V-abs]
- This paper has no avalanche analysis. It is the **raw-data source** for building avalanche-versus-DIV trajectories at several densities (see Q6).

**Stewart & Plenz 2008** — "Homeostasis of neuronal avalanches during postnatal cortex development in vitro", *J Neurosci Methods* 169(2): 405–416. doi:10.1016/j.jneumeth.2007.10.021 **[U for citation details]**
- Tetzlaff et al. summarize it as follows. In organotypic cortex cultures, LFP amplitudes grow until about 25 DIV. Avalanche distributions generally have exponent about −1.5, consistent with branching parameter 1. Transitory sub- or supercritical stages "are also observable", for example in their Fig. 3D. — [Tetzlaff et al. 2010 discussion](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1001013) [V as secondary report]

**Thivierge & Tauskela 2014** — "Development of avalanches and efficient communication in neuronal networks", *BMC Neurosci* 15(Suppl 1): P31 (CNS meeting abstract). doi:10.1186/1471-2202-15-S1-P31 **[V-abs, conference abstract]**
- 643,039 avalanches from 15 cultures grown to 35 DIV. Power-law scaling of size and duration **emerged gradually** over development. — [PMC4126459](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4126459/) [V via search summary; not peer-reviewed full paper]

**Massobrio, de Arcangelis, Pasquale, Jensen & Plenz 2015** — "Criticality as a signature of healthy neural systems", *Front Syst Neurosci* 9: 22. doi:10.3389/fnsys.2015.00022 — [Frontiers](https://doi.org/10.3389/fnsys.2015.00022) **[V citation; content U]**
- A review of in vitro (dissociated, organotypic) and in vivo evidence. It argues that deviations from criticality mark pathological states. Its abstract was not indexed, so I did not extract content-level claims.

**2019–2026 developmental work (dissociated cultures and iPSC networks)**
- **Zhang, Yaron, Akita, Shiramatsu, Chao & Takahashi 2025** — "Emergence of deviance detection in cortical cultures through maturation, criticality, and early experience", arXiv:2510.00764 (preprint, 1 Oct 2025). Dissociated cortical cultures on HD-CMOS MEAs from 10 to 35 DIV. Avalanche analysis showed "a gradual transition from subcritical to near-critical dynamics". Robust late evoked responses appeared from DIV 15–20. — [arXiv](https://arxiv.org/abs/2510.00764) [V-abs; preprint]
- **Ikeda, Akita & Takahashi 2025** — "Emergent functions of noise-driven spontaneous activity: Homeostatic maintenance of criticality and memory consolidation", arXiv:2502.10946. SNN simulations plus dissociated cultures. Repetitive stimulation transiently disrupts criticality and E/I balance, and spontaneous activity restores them. This is an in vitro homeostasis-of-criticality analogue. — [arXiv](https://arxiv.org/abs/2502.10946) [V-abs; preprint]
- **"Self-organization of high-dimensional geometry of neural activity in culture"**, bioRxiv 2025.01.14.630383. Rat cortical cultures on HD-MEAs, recorded daily. Power-law avalanches and power-law PCA spectra emerge together during development. The paper cites the integration-then-fragmentation account. — [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.01.14.630383v1) [V via search snippet only; authors not confirmed]
- **iPSC-derived dopaminergic networks** (Heiney et al., arXiv:1907.02351 and 1907.13118; *IJUC* 2021). DIV 21–56: the network passed through a period of relatively stable criticality but did not settle into it. In a second study, no network met the criticality criterion (p > 0.1) during DIV 14–51. This shows that iPSC networks need not converge. — [arXiv 1907.02351](https://arxiv.org/abs/1907.02351); [arXiv 1907.13118](https://arxiv.org/abs/1907.13118); [IJUC PDF](https://www.oldcitypublishing.com/wp-content/uploads/2021/05/IJUCv16n4p303-326Heiney.pdf) [V via search summaries]
- **E/I ratio in mouse cortical cultures, *Stem Cell Reports* 2025** ("A proper excitatory/inhibitory ratio is required to develop synchronized network activity in mouse cortical cultures"). A longitudinal study on a 4,096-channel HD-MEA. Correlated or synchronized network activity develops only at an appropriate E/I ratio. — [Cell](https://www.cell.com/stem-cell-reports/fulltext/S2213-6711(25)00250-4); [bioRxiv](https://www.biorxiv.org/content/10.1101/2025.02.28.640720v1.full) [V via search summary; whether it includes avalanche analysis is unverified]
- **Habibollahi, Kagan, Burkitt & French 2023** — "Critical dynamics arise during structured information presentation within embodied in vitro neuronal networks", *Nat Commun* **14: 5287** (not 5409). doi:10.1038/s41467-023-41020-3 [V]
  - Human iPSC-derived or E15 mouse cortical cells on HD-MEAs (DishBrain). 14 cultures, 308 sessions (192 Gameplay, 116 Rest), 50 ms bins. — [PMC10469171](https://pmc.ncbi.nlm.nih.gov/articles/PMC10469171) [V]
  - Metrics:
    - DCC = |β_pred − β_fit| with β_pred = (α−1)/(τ−1), from a truncated power-law MLE fit
    - BR = active units at t+1 / active units at t
    - shape-collapse error
  - Networks were near-critical in Gameplay, not at Rest. Example value: Gameplay DCC 0.296 ± 0.015 versus time-shuffled 0.627 ± 0.087 (p < 0.0005). — [PMC10469171](https://pmc.ncbi.nlm.nih.gov/articles/PMC10469171) [V]
  - The authors call the results "preliminary evidence at best" because the criticality metrics did not all agree. — [PMC10469171](https://pmc.ncbi.nlm.nih.gov/articles/PMC10469171) [V]
  - This paper is about stimulation state, not development.
- **Organoids**
  - **Trujillo et al. 2019**, *Cell Stem Cell* 25(4): 558–569.e7, doi:10.1016/j.stem.2019.08.002. Human cortical organoids tracked over months on MEAs. Activity rises steadily, and periodic oscillatory events become spatiotemporally irregular, resembling preterm EEG. The abstract has no avalanche or criticality analysis. — [Europe PMC / PMC6778040](https://pmc.ncbi.nlm.nih.gov/articles/PMC6778040) [V-abs]
  - **Sharf et al. 2022**, *Nat Commun* **13: 4403** (not 4888), doi:10.1038/s41467-022-32115-4. HD-CMOS MEA plus shank recordings of organoids. Functional connectivity shows many weak edges on a skeleton of fewer strong edges, and LFP phase-locks to spikes. The abstract has no avalanche or criticality analysis. — [PMC9338020](https://pmc.ncbi.nlm.nih.gov/articles/PMC9338020) [V-abs]
  - **Alam El Din, …, Habibollahi, …, Kagan, Hartung & Smirnova 2025**, "Human neural organoid microphysiological systems show the building blocks necessary for basic learning and memory", *Commun Biol* 8: 1237, doi:10.1038/s42003-025-08632-5. HD-MEA organoids aged 6–9 weeks versus 10–13 weeks, followed for about 3 weeks (2 experiments, 5–6 wells per group). — [PMC12357958](https://pmc.ncbi.nlm.nih.gov/articles/PMC12357958) [V]
    - The older group had consistently lower, tighter DCC and BR closer to 1.
    - BR in the younger group rose significantly over 3 weeks ("maturation and pursuit of criticality").
    - Theta-burst stimulation raised BR significantly.
    - This is the only organoid developmental-criticality trajectory I found. — [PMC12357958](https://pmc.ncbi.nlm.nih.gov/articles/PMC12357958) [V]
- **Kapucu**: I found no Kapucu paper on developmental *criticality*. [U] Kapucu et al. published burst-analysis methods for developing hPSC-derived networks (e.g., *Front Comput Neurosci* 2012). Not verified this session.

### Inferences
- The supercritical or bimodal "bursting" phase appears early (about 2–3 weeks in rat cortical cultures). It is the most reproducible feature across Tetzlaff, Yada, Pasquale and Wagenaar ("after two weeks network activity dominated by population bursts"). A growth model with homeostasis should produce a transient overshoot in connectivity or bursting before settling. If it reaches criticality monotonically from below, it disagrees with these data.
- The order of the "subcritical" and "supercritical" phases is not universal. Tetzlaff has supercritical then subcritical then critical. Yada has exponential then bimodal then power law. Zhang 2025 has a monotone sub → near-critical. Compare against the shape of the trajectory and the presence of a burst-dominated phase, not against exact DIV values.
- Time to criticality ranges from about 4 weeks (Yada, HD-CMOS, single-cell) to about 8 weeks (Tetzlaff, 59-electrode). Preparation, density, electrode sampling and bin choice all matter. A simulation should map model time to DIV through a growth observable such as synapses per neuron (see Q5), not assume a fixed mapping.
- Pasquale's correlate of criticality (medium burst synchronization, high burst/spike ratio) gives a second, cheaper observable to match.

### Gaps
- I could not extract Tetzlaff Table 4's numerical exponents and bin widths for each state. They are in the PLoS article (Table 4).
- Yada 2017: the exact DIV boundaries of each stage, the fitted exponents, and whether branching ratio was computed all need the full text (Elsevier, paywalled).
- Pasquale 2008: the exponents and DIV ranges need the full text.
- I found no developmental study reporting all of size exponent τ, duration exponent α, the crackling-noise relation (γ vs (α−1)/(τ−1)) and branching ratio as a function of DIV in one series. Habibollahi-style DCC analyses cover this only for stimulation state and organoid age.

---

## Q2. In vivo homeostasis of criticality: timescales and effect sizes

### Takeaway
In freely behaving rat V1, network dynamics sit near criticality and are actively restored after perturbation. After monocular deprivation, criticality is disrupted and returns within about 48 h, more than 30 h before firing rates recover (Ma et al. 2019). During normal waking, criticality degrades progressively over hours, and sleep restores it (Xu et al. 2024). Adaptation to strong sensory drive retunes visual cortex to criticality on a much faster, sub-second to seconds scale (Shew et al. 2015). Hengen & Shew 2025 frame criticality as a homeostatic setpoint and support this with a meta-analysis of 140 datasets.

### Cited Findings
- **Ma, Turrigiano, Wessel & Hengen 2019**, "Cortical circuit dynamics are homeostatically tuned to criticality in vivo", *Neuron* 104(4): 655–664.e4, doi:10.1016/j.neuron.2019.08.031. — [Europe PMC abstract](https://pmc.ncbi.nlm.nih.gov/articles/PMC6934140) [V-abs]
  - Continuous recording of V1 in freely behaving rats for **9 days**. — [V-abs]
  - Under control conditions, networks were robustly organized around criticality in both light and dark. — [V-abs]
  - Visual (monocular) deprivation severely disrupted criticality, which then **returned over 48 h**. — [V-abs]
  - Recovery of excitatory dynamics **preceded firing-rate homeostasis by >30 h**. Criticality was disrupted before firing rates fell. — [V-abs]; [Neuroscience News summary](https://neurosciencenews.com/information-processing-tuning-15035/)
  - Modeling suggests criticality is set by **inhibitory plasticity and architecture**, not by excitatory firing-rate homeostasis. — [V-abs]
  - Effect sizes (DCC or exponent magnitudes, deprivation day and age) were not extracted. [Gap]
- **Xu, Schneider, Wessel & Hengen 2024**, "Sleep restores an optimal computational regime in cortical networks", *Nat Neurosci* 27: 328–338 (online 5 Jan 2024), doi:10.1038/s41593-023-01536-9. — [Nature](https://www.nature.com/articles/s41593-023-01536-9); [PMC11272063](https://pmc.ncbi.nlm.nih.gov/articles/PMC11272063) [V]
  - Rat V1 recorded for **10–14 days** in freely behaving animals. — [V]
  - Waking experience **progressively** disrupts criticality, and sleep restores it. — [V]
  - Deviation from criticality predicts future sleep/wake behavior better than prior behavior or slow-wave activity does. — [V]
  - Source data and analysis code are public: [github.com/hengenlab/Sleep_restores_criticality](https://github.com/hengenlab/Sleep_restores_criticality). — [V via search summary]
  - The per-hour drift rate and the magnitude of change in the criticality metric were not extracted. [Gap]
- **Hengen & Shew 2025**, "Is criticality a unified setpoint of brain function?", *Neuron* 113(16): 2582–2598.e2, doi:10.1016/j.neuron.2025.05.020. Preprint: bioRxiv 10.1101/2024.09.02.610815. — [Europe PMC / PMC12374783](https://pmc.ncbi.nlm.nih.gov/articles/PMC12374783/); [bioRxiv](https://www.biorxiv.org/content/10.1101/2024.09.02.610815) [V: existence, venue, abstract]
  - A review plus a **meta-analysis of 140 datasets published 2003–2024**. — [V-abs]
  - The paper claims that "a long-standing controversy is the product of a methodological choice with no bearing on underlying dynamics". The specific choice was not extracted: the PMC and bioRxiv fetches were blocked. [Gap]
  - It frames criticality as a setpoint reached through development, plasticity and homeostasis. — [V-abs]
  - Press coverage adds a disease framing, describing neurodegeneration as erosion of criticality. This is the authors' hypothesis. — [WashU Source](https://source.washu.edu/2025/06/a-unified-theory-of-the-mind/)
- **Shew, Clawson, Pobst, Karimipanah, Wright & Wessel 2015**, "Adaptation to sensory input tunes visual cortex to criticality", *Nat Phys* 11: 659–663, doi:10.1038/nphys3370. — [doi.org](https://doi.org/10.1038/nphys3370) [V for citation and core claim; author list U]
  - Turtle visual cortex ex vivo (eye-attached), MEA LFP, visual stimulation of the retina.
  - Stimulus onset produces non-critical dynamics: very large avalanches and intense activity. Adaptation then rapidly retunes the network toward criticality, with smaller and more diverse avalanches. — [summary via search; PLoS CB follow-up PMC5469508](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5469508/) [V core claim]
  - Follow-up: Clawson, Wright, Wessel & Shew 2017, *PLoS Comput Biol* 13(5): e1005574. Adaptation toward scale-free dynamics improves stimulus discrimination but reduces detection. — [PMC5469508](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5469508/) [V]
  - [U] Adaptation timescale is on the order of ~1 s after stimulus onset. Not confirmed this session.

### Inferences
- There are three separable timescales to compare against:
  1. sub-second to seconds: sensory adaptation (Shew 2015)
  2. hours: wake-related drift and sleep-related restoration (Xu 2024)
  3. ~2 days: recovery from deprivation (Ma 2019)
  - Firing-rate homeostasis is slower than criticality restoration (>30 h lag in Ma 2019). A model in which criticality emerges only as a by-product of a firing-rate set-point would predict the opposite order, so Ma 2019 is a direct test for a firing-rate-homeostasis growth model.
- Ma 2019 emphasizes inhibition. A purely excitatory growth model with firing-rate homeostasis may fail to reproduce "criticality restored before rates".

### Gaps
- Quantitative effect sizes were not extracted:
  - Ma 2019: DCC or κ during deprivation versus baseline
  - Xu 2024: DCC drift per hour of wake
- The specific methodological choice that Hengen & Shew identify as resolving the controversy is not recorded here.
- Shew 2015's exact adaptation time constant is not recorded here.

---

## Q3. Dynamic range maximized at criticality: experimental evidence

### Takeaway
Two experimental studies show dynamic range peaks near criticality. One uses organotypic cortex cultures with pharmacologically shifted E/I balance (Shew et al. 2009). The other uses anesthetized rat somatosensory cortex in vivo (Gautam et al. 2015). Both give a tuning curve of dynamic range against distance from criticality that a simulation can compare to qualitatively.

### Cited Findings
- **Shew, Yang, Petermann, Roy & Plenz 2009**, "Neuronal avalanches imply maximum dynamic range in cortical networks at criticality", *J Neurosci* 29(49): 15595–15600, doi:10.1523/JNEUROSCI.3864-09.2009. arXiv:0906.0527. — [Europe PMC abstract](https://doi.org/10.1523/JNEUROSCI.3864-09.2009); [arXiv](https://arxiv.org/abs/0906.0527) [V-abs]
  - Cortex slice (organotypic) cultures on planar MEAs. E/I ratio changed pharmacologically, giving a "network tuning curve for stimulus processing as a function of distance from criticality" that agrees with simulations. Dynamic range is maximized where avalanches occur. — [V-abs]
  - [U] Distance from criticality was measured by κ, which compares the measured size distribution to a reference −3/2 power law. Dynamic range peaked at κ ≈ 1. AP5/DNQX (reduced excitation) gave κ < 1, and PTX (reduced inhibition) gave κ > 1. This matches my recollection of the paper but was not verified this session (PDF text extraction unavailable).
- **Gautam, Hoang, McClanahan, Grady & Shew 2015**, "Maximizing sensory dynamic range by tuning the cortical state to criticality", *PLoS Comput Biol* 11(12): e1004576, doi:10.1371/journal.pcbi.1004576. — [PLoS](https://doi.org/10.1371/journal.pcbi.1004576) [V-abs]
  - Urethane-anesthetized rat barrel cortex, MEA, whisker stimulation. Cortical state was varied naturally and by manipulating inhibition. Dynamic range was maximized at criticality, "near the tipping point midway between" weakly and strongly correlated extremes. — [V-abs]
  - The model needed **activity-dependent depression** to reproduce the data. — [V-abs]

### Inferences
- Both studies define distance from criticality with a size-distribution metric (κ) rather than branching ratio. A simulation that reports both κ and BR could be placed on their tuning curves.
- The Gautam result means a simulation without synaptic depression may get the location of the dynamic-range peak wrong.

### Gaps
- I did not extract numerical dynamic-range values in dB, or peak widths, from either paper.

---

## Q4. Micropatterned and modular cultures: does geometry shape avalanche or dynamical statistics?

### Takeaway
Engineered modular cultures (Yamamoto et al. 2018, 2023) and self-organized clustered cultures (Okujeni et al. 2017, 2019) show that mesoscale architecture changes synchrony, burst initiation and dynamical richness. Richness peaks at intermediate coupling or clustering. Yamamoto 2023 has public data and code. Okujeni & Egert 2019 is the closest experimental analogue of growth plus homeostasis: activity-dependent outgrowth and migration produce different modularity, yet Ca²⁺ influx converges to a common set-point. None of these papers directly reports avalanche exponents against geometry.

### Cited Findings
- **Yamamoto, Moriya, Ide, Hayakawa, Akima, Sato, Kubota, Tanii, Niwano, Teller, Soriano & Hirano-Iwata 2018**, "Impact of modular organization on dynamical richness in cortical networks", *Sci Adv* 4(11): eaau4914, doi:10.1126/sciadv.aau4914. — [Europe PMC abstract](https://doi.org/10.1126/sciadv.aau4914) [V-abs]
  - Four interconnected micropatterned modules of rat cortical neurons.
  - Highest dynamical richness occurs at a **critical connectivity at the verge of physical disconnection**. Stronger coupling gives persistent coherence, and weaker coupling localizes activity within modules.
  - In the model, coherence arises from a trade-off between connectivity and "subquorum firing". — [V-abs]
- **Yamamoto, Spitzner, Takemuro, Buendía, Murota, Morante, Konno, Sato, Hirano-Iwata, Levina, Priesemann, Muñoz, Zierenberg & Soriano 2023**, "Modular architecture facilitates noise-driven control of synchrony in neuronal networks", *Sci Adv* 9(34): eade1755, doi:10.1126/sciadv.ade1755. — [PMC10456864](https://pmc.ncbi.nlm.nih.gov/articles/PMC10456864) [V]
  - Rat cortical neurons on micropatterned substrates with three modular topologies (single-bond, triple-bond, merged). GCaMP6s calcium imaging at **DIV 10–11**, with optogenetic (ChrimsonR) asynchronous stimulation. — [V]
  - Modularity increases sensitivity to noise-like stimulation: stimulation decorrelates and desynchronizes modular networks but not merged ones. The mechanism is persistent depletion of synaptic resources in stimulated neurons. Reported measures include event sizes, correlation coefficients and "functional complexity". — [V]
  - **Data are public on G-Node GIN, doi:10.12751/g-node.t77b3p. Code is on Zenodo, doi:10.5281/zenodo.7962283, and on GitHub.** — [PMC10456864 Data availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC10456864) [V]
- **Okujeni, Kandler & Egert 2017**, "Mesoscale architecture shapes initiation and richness of spontaneous network activity", *J Neurosci* 37(14): 3972–3987, doi:10.1523/JNEUROSCI.2552-16.2017. — [Europe PMC abstract](https://doi.org/10.1523/JNEUROSCI.2552-16.2017) [V-abs]
  - Rat cortical neurons self-arrange into clusters with high neurite density and fasciculated long-range axons. PKC inhibition gave uniform networks, and PKC activation gave clustered ones with bundled axons.
  - Clustered networks were more active, had more diverse patterns, stronger synaptic input and larger membrane fluctuations, and initiated more bursts, but recruitment was incomplete. **Moderately clustered networks were optimal.** — [V-abs]
- **Okujeni & Egert 2019**, "Self-organization of modular network architecture by activity-dependent neuronal migration and outgrowth", *eLife* 8: e47996, doi:10.7554/eLife.47996. — [eLife](https://doi.org/10.7554/eLife.47996) [V-abs]
  - Simple rules for neurite outgrowth and cell migration explain variation in modularity. Activity and burst rate increased with clustering, while peak in-burst rates were higher in homogeneous networks.
  - **Long-term Ca²⁺ influx converged during maturation despite very different architectures, neurite extent, connectivity, modularity and activity**, which indicates homeostasis toward a common Ca²⁺ set-point. — [V-abs]
  - This is directly relevant to a contact-based growth model with homeostasis.

### Inferences
- The modular-culture studies support a non-monotonic dependence on geometry: an optimum at intermediate coupling (Yamamoto 2018) or intermediate clustering (Okujeni 2017). A 2D growth simulation that varies spatial clustering could test whether its avalanche statistics show a similar optimum.
- Okujeni & Egert 2019 is a strong qualitative benchmark. The homeostatic variable converges across architectures while structure differs. A firing-rate-homeostasis growth model should show the same convergence of rate or Ca proxy across initial geometries, even when avalanche statistics do not converge.

### Gaps
- I found no study that reports avalanche size or duration exponents as a function of micropattern geometry or measured neurite morphology. Yamamoto 2023 data (GIN) could be reanalysed for this, but calcium imaging limits temporal resolution.
- Data availability for the Okujeni papers was not checked.

---

## Q5. Neurite growth and synapse counts during culture development (for a "contact granularity" parameter)

### Takeaway
In rat cortical cultures, synapses per neuron rise from under 5 at 7 DIV to about 150–400 by 21 DIV, a 30–80-fold increase, and depend inversely on plating density (Cullen et al. 2010). Synaptophysin levels rise steadily through 28 DIV (Harrill et al. 2015). Synapse growth over weeks 1–3 therefore spans roughly two orders of magnitude, which coincides with the onset of bursting and the supercritical or bimodal phase.

### Cited Findings
- **Cullen, Gilroy, Irons & LaPlaca 2010**, "Synapse-to-neuron ratio is inversely related to neuronal density in mature neuronal cultures", *Brain Res* 1359: 44–55, doi:10.1016/j.brainres.2010.08.058. — [PMC2969840](https://pmc.ncbi.nlm.nih.gov/articles/PMC2969840) [V-abs; volume and pages U]
  - Rat cortical neurons at 10–5000 neurons/mm², immunocytochemistry over 21 DIV.
  - **7 DIV: under 5 synapses per neuron**, independent of density.
  - **21 DIV: 30–80-fold increase.** About 400 synapses per neuron at under 500 neurons/mm², and about 150 at 500–4500 neurons/mm² (p < 0.05).
  - High-density networks were up to 18 µm thick at 21 DIV, so planar cultures are partly 3D. — [V-abs]
  - The BioNumbers entry lists about 388 and 157 synapses per neuron for the same groups. — [BioNumbers 106651](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=106651)
- **Harrill, Chen, Streifel, Yang, Mundy & Lein 2015**, "Ontogeny of biochemical, morphological and functional parameters of synaptogenesis in primary cultures of rat hippocampal and cortical neurons", *Mol Brain* 8: 10, doi:10.1186/s13041-015-0099-9. — [PMC4339650](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4339650/) [V-abs]
  - Over 28 DIV, synaptophysin rose steadily in both hippocampal and cortical cultures.
  - High-content imaging showed more synapses in hippocampal than cortical neurons and differing inhibitory synapse ontogeny.
  - Cortical cultures had higher firing and bursting rates on MEAs, possibly from weaker inhibitory tone. — [V-abs]
- Synaptic release is detectable by **3–4 DIV** in high-density cultures. — [PubMed 9741482](https://pubmed.ncbi.nlm.nih.gov/9741482/) [V via search snippet; citation details U]
- [U] Ichikawa et al. 1993 (*Neurosci Res* 16: 95–103) reported EM synapse density in rat cortical cultures rising through about 3 weeks and then plateauing. van Pelt et al. 2004 (*IEEE TBME*) linked the plateau of burst firing to synapse maturation. Neither was verified this session.
- Wagenaar et al. 2006: responses to electrical stimulation suggest axonal outgrowth proceeds faster in denser cultures. — [abstract](https://europepmc.org/search?query=DOI%3A%2210.1186%2F1471-2202-7-11%22) [V-abs]
- Okujeni & Egert 2019 relate neurite extent and connectivity to clustering, which gives an experimental link between morphology and connectivity. — [eLife](https://doi.org/10.7554/eLife.47996) [V-abs]

### Inferences
- Model time could be calibrated by matching synapses per neuron:
  - under 5 at the stage equivalent to "7 DIV"
  - about 150–400 at "21 DIV"
- The inverse dependence on density (more synapses per neuron when sparse) is a testable prediction for a contact-based model. With fixed contact-to-synapse rules, does synapse count per neuron fall with density as neurites overlap?
- If the model's "contact granularity" sets synapses per contact, the 30–80× rise from 7 to 21 DIV bounds how fast effective coupling should grow over the window in which bursting appears (around 2 weeks).

### Gaps
- I found no single study that pairs synapse counts with avalanche statistics in the same cultures over DIV.
- I found no per-neuron dendritic length or branching time course for cortical cultures this session.

---

## Q6. Public datasets

### Takeaway
The best public datasets for developmental trajectories are:
- **Wagenaar/Potter developing cortical cultures**: 58 cultures, about 5 weeks, several densities, daily spontaneous recordings, about 45 GB
- **Hengen lab** sleep and criticality data and code (in vivo, 10–14 days)
- **Yamamoto 2023** modular cultures (GIN)
- **Habibollahi 2023** DishBrain (OSF)

I did not confirm any DANDI or CRCNS dataset that gives a dissociated-culture developmental avalanche series.

### Cited Findings
- **Wagenaar/Potter "Network activity of developing cortical cultures in vitro"**:
  - Organized by density (dense, small, sparse, small-sparse "smsp", ultra-sparse "ulsp"), plating batch (1–7), culture number (1–6) and DIV.
  - Daily and overnight recordings, both spontaneous and stimulus-evoked, as Matlab (v7+) or compact text.
  - wget list files are provided for bulk download. The full archive is about 45 GB.
  - The site notes some files may be lost or damaged. Use requires citing Wagenaar et al. 2006.
  - Links: [help page](https://potterlab.bme.gatech.edu/development-data/html/help.html); [dense daily stim](https://potterlab.bme.gatech.edu/development-data/html/daily.stim.dense.full.html); [overnight spont dense](https://potterlab.bme.gatech.edu/development-data/html/overnight.spont.dense.full.html); [mirror](https://neurodatasharing.bme.gatech.edu/development-data/html/daily.spont.small.text.html) [V that the pages exist; download links not tested]
- **Hengen lab, Xu et al. 2024**: source data and analysis code at [github.com/hengenlab/Sleep_restores_criticality](https://github.com/hengenlab/Sleep_restores_criticality). [V via search summary; whether raw spike data are included is U]
- **Yamamoto et al. 2023**: data and processed results at [G-Node GIN doi:10.12751/g-node.t77b3p](https://doi.org/10.12751/g-node.t77b3p); code at [Zenodo doi:10.5281/zenodo.7962283](https://doi.org/10.5281/zenodo.7962283). [V]
- **Habibollahi et al. 2023**: all data on OSF at [osf.io/ncvpq](https://osf.io/ncvpq/?view_only=8fc5fc5aad254fce92a79390ae84b81c), plus Source Data. [V]
- **Alam El Din et al. 2025** (organoid criticality versus age): only figure source data (Supplementary Data 1–2) and RNA-seq ([GEO GSE301725](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301725)). Raw HD-MEA data are available "on request". — [PMC12357958](https://pmc.ncbi.nlm.nih.gov/articles/PMC12357958) [V]
- [U] **CRCNS**:
  - Beggs-lab organotypic culture datasets of spontaneous spiking on a 512-electrode array (Ito, Timme et al.; CRCNS ssc-series) suit mature-network avalanche comparison but are not developmental.
  - Hengen lab in vivo datasets may also be on CRCNS.
  - Neither was verified this session.
- [U] **Charlesworth, Cotterill, Morton, Grant & Eglen 2015**, *Neural Dev* 10: 1, doi:10.1186/s13064-014-0028-0. Developmental profiles of cortical and hippocampal culture activity, with data released publicly (CARMEN or GitHub). **Cotterill et al. 2016**, *J Biomol Screen*, characterized early network development in multiwell MEAs. Both are likely sources of DIV-resolved spike trains. Not verified.

### Inferences
- The Wagenaar dataset is the only verified, open, multi-density, developmental spike dataset found. Its density range (3,000–50,000 neurons) is useful because density changes both bursting onset and synapse counts per neuron (Cullen 2010). Avalanche statistics against DIV and density can be computed from it directly. With 59-channel sampling, subsampling effects on exponents must be handled.
- HD-MEA developmental data (Yada 2017; Zhang 2025; 2025 bioRxiv) would be better for single-cell avalanche statistics, but I found no public deposits.

### Gaps
- I found no DANDI or Zenodo deposit of a dissociated-culture developmental HD-MEA dataset with avalanche analysis. The DANDI search was not run directly; it is still worth a query for "MaxOne" or "HD-MEA development".
- Data availability for Yada 2017, Tetzlaff 2010, Pasquale 2008 and the Takahashi-group preprints was not determined.
