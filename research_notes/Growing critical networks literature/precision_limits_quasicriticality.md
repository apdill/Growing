# Precision limits of homeostatic / structural-plasticity criticality: quasi-criticality, fluctuations, timescale separation, drive, finite size

Scope note: about 20 tool calls. Full texts were read (via ar5iv HTML) for Zierenberg et al. 2018 PRX, Das & Levina 2019 PRX, Kinouchi et al. 2019 Sci Rep, Kossio et al. 2018 PRL, and Buendía et al. 2020 Front Phys. The other papers were checked only through abstracts or search snippets, and each is marked accordingly. "VERIFIED" means I read the claim in the source during this session. "UNCERTAIN" means it comes from a snippet, from secondary citation, or from my prior knowledge, and needs checking before it is cited.

## Q1. Self-organized quasi-criticality (SOqC): what is established about fluctuations around the critical point and about finite-size scaling?

### Takeaway
SOqC theory, mainly Bonachela & Muñoz 2009 and Bonachela et al. 2010, establishes qualitatively that non-conserving or homeostatic feedback does not converge to the critical point. The system hovers with finite-amplitude excursions to both sides. Bonachela et al. claim the excursion amplitude does not vanish with N. Kinouchi et al. 2019 give the most quantitative treatment: the mean-field fixed-point offset is O(1/τ) in the adaptation timescale τ, and σ fluctuations plateau below a τ-dependent size N(τ). None of these papers derives a general law for fluctuation amplitude as a function of the size of the discrete plasticity step.

### Cited Findings
- **Bonachela, de Franciscis, Torres & Muñoz 2010**, "Self-organization without conservation: are neuronal avalanches generically critical?", J. Stat. Mech. P02015, DOI 10.1088/1742-5468/2010/02/P02015, arXiv:1001.3256. The abstract says the Levina–Herrmann–Geisel (LHG) model is analogous to forest-fire and earthquake models, which lack true criticality. Without fine-tuning, such networks are "either subcritical or supercritical", although the pseudo-critical region is relatively broad. VERIFIED (abstract) — [Strathclyde record](https://pureportal.strath.ac.uk/en/publications/self-organization-without-conservation-are-neuronal-avalanches-ge/); [arXiv](https://arxiv.org/abs/1001.3256)
- Kinouchi et al. 2019 summarise Bonachela et al. 2010 as claiming that in the LHG model the oscillation amplitude "basically does not change with N" and stays non-zero in the thermodynamic limit, which Bonachela et al. treat as a core ingredient of SOqC. VERIFIED as a secondary citation; I did not read Bonachela 2010 itself — [Kinouchi et al. 2019, ar5iv](https://ar5iv.labs.arxiv.org/html/1803.05537)
- **Bonachela & Muñoz 2009**, "Self-organization without conservation: true or just apparent scale-invariance?", J. Stat. Mech. P09009, DOI 10.1088/1742-5468/2009/09/P09009. This is the founding SOqC paper. As reviewed by Buendía et al. 2020: with a non-vanishing dissipation term, "the average energy does not self-tune to a critical value", but "keeps hovering around it, with excursions of finite amplitude to both sides" of the transition. The energy distribution is "peaked around E_c" but much broader than in conserving SOC. Sweeping through the transition "might suffice to induce approximate or 'dirty' scaling behavior, but not strict 'bona fide' scale invariance". For large enough systems with offline driving, dynamical percolation (DP) dominates. VERIFIED via the Buendía review; I did not read the 2009 original — [Buendía et al. 2020, ar5iv](https://ar5iv.arxiv.org/html/2006.03020)
- Buendía et al. 2020 state that the ratio of driving to dissipation constants, h/ε (and also E_max), determines the state: dissipation-dominated is subcritical, strongly driven is supercritical, and intermediate hovers near criticality (SOqC). They also note that "the choice of parameters (and not only system size) can determine the 'cleanliness' of the observed scaling behavior." The review gives no closed-form expression for distance or fluctuation amplitude as a function of h, ε or N. VERIFIED — [Buendía et al. 2020](https://ar5iv.arxiv.org/html/2006.03020)
- **Kinouchi, Brochini, Costa, Campos & Copelli 2019**, "Stochastic oscillations and dragon king avalanches in self-organized quasi-critical systems", Sci. Rep. 9:3874, DOI 10.1038/s41598-019-40473-1, arXiv:1803.05537. They do a linear stability analysis of mean-field fixed points with adaptive gain or synapses. Main results (VERIFIED from full text):
  - Simplified gain model: Γ* = Γ_c/(1 − 2/τ) (so Γ* ≈ 1 + 2/τ), with eigenvalue modulus |λ±| = 1 − O(1/τ) (0.990, 0.998 and 0.999 at τ = 100, 500, 1000) and oscillation frequency ω ≃ τ^(−1/2).
  - LHG-type gain model: Γ* = (2A + τu)/(2 + τu) ≈ Γ_c + 2(A−1)/(τu).
  - Cellular automaton (CA) with LHG synapses: σ* = 1 + (A−1)/(1 + uτ), and on the complete graph |λ±| = 1 − (A − 1 + u/2)/(uτ).
  - The fixed point is therefore a weakly stable focus, close to a Neimark–Sacker-like bifurcation. Demographic noise excites "stochastic oscillations" that hit the absorbing state, producing both power-law avalanches and "dragon king" events.
  - [ar5iv](https://ar5iv.labs.arxiv.org/html/1803.05537); [Sci Rep / PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6405991/)
- Kinouchi et al. 2019, finite size (VERIFIED):
  - In the CA model ⟨σ⟩ → σ* as N grows. The standard deviation Δσ stays roughly constant for small N and starts to decrease once N exceeds a τ-dependent threshold N(τ). Larger τ extends this plateau to larger N.
  - For any finite τ, the stochastic oscillations (SO) "do not survive up to the thermodynamic limit".
  - LHG used τ = τ0·N, and under that scaling both fluctuations and damping decay with N.
  - External noise "never vanishes in the thermodynamic limit", so the authors argue SO and dragon kings would persist in real systems.
  - No explicit N^(−a) law for distance or fluctuation is given.
  - [ar5iv](https://ar5iv.labs.arxiv.org/html/1803.05537)
- **Levina, Herrmann & Geisel 2007**, "Dynamical synapses causing self-organized criticality in neural networks", Nat. Phys. 3:857–860, DOI 10.1038/nphys758. This is the prototypical adaptive SOC (aSOC) model, with depressing synapses. Its synaptic recovery time is scaled as τ = τ0·N so that the fixed point converges to criticality as N → ∞. VERIFIED via Kinouchi 2019's description; the DOI is from prior knowledge and is UNCERTAIN — [Kinouchi et al. 2019](https://ar5iv.labs.arxiv.org/html/1803.05537)
- **Brochini, Costa, Abadi, Roque, Stolfi & Kinouchi 2016**, "Phase transitions and self-organized criticality in networks of stochastic spiking neurons", Sci. Rep. 6:35831, DOI 10.1038/srep35831. This paper introduces the stochastic neuron model with gain adaptation. Mean-field exponents are 3/2 (size) and 2 (duration) on the static critical line, and adaptive gain gives SOqC. VERIFIED as a citation in Kinouchi 2019; the DOI is constructed from the volume and article number and is UNCERTAIN — [Kinouchi et al. 2019](https://ar5iv.labs.arxiv.org/html/1803.05537)

### Inferences
- The project's two observations are both qualitatively predicted by SOqC/aSOC theory: σ oscillating around a set point, and avalanches started at high σ being far larger. Kinouchi 2019 explicitly links supercritical excursions of an oscillating fixed point to dragon-king avalanches.
- In Kinouchi 2019 the mean offset is O(1/τ) and the fluctuation Δσ shrinks with N only once N exceeds N(τ). This is the closest existing quantitative "precision" result, but it is tied to continuous adaptation variables plus demographic noise. It does not involve a quantized plasticity step.

### Gaps
- I did not read the full text of Bonachela & Muñoz 2009 or Bonachela et al. 2010, so I cannot confirm whether they give an explicit scaling exponent for excursion amplitude versus N or versus the dissipation rate.

## Q2. Zierenberg, Wilting & Priesemann 2018 (PRX): is m = 1 − h/r* explicit? Do they analyse fluctuations, finite size, or the homeostatic timescale?

### Takeaway
Yes. m = 1 − h/r* is their explicit mean-field result (Eq. 7), so distance from criticality is linear in the drive. They analyse τ_hp only as a stability condition (bursting versus a stable fixed point) and fix N = 10⁴. There is no finite-size scaling and no theory of fluctuation amplitude in m.

### Cited Findings
- **Zierenberg, Wilting & Priesemann 2018**, "Homeostatic plasticity and external input shape neural network dynamics", Phys. Rev. X 8:031018, DOI 10.1103/PhysRevX.8.031018, arXiv:1807.01479. All points below are VERIFIED from the ar5iv full text — [ar5iv](https://ar5iv.labs.arxiv.org/html/1807.01479); [arXiv abstract](https://arxiv.org/abs/1807.01479)
  - Eq. 4 gives r = h/(1−m). Eq. 7 gives the mean-field solution m = 1 − h/r* and τ = −Δt/ln(1 − h/r*). Equivalently h/r = 1 − e^(−Δt/τ).
  - The homeostatic rule (Eq. 6) is local and deterministic. The mean-field analysis assumes m is approximately constant.
  - In the bursting regime, "the (rapid) homeostatic feedback causes notable changes of the network branching parameter". Fig. 7 reports m fluctuating around the prediction, with "a distinct saw-tooth pattern for small input rates".
  - Stability condition (Eqs. 8–10): with effective timescale τ′ = τ_hp/k̄, a stable fixed point requires τ ≪ τ′·|Δt r* − A_t/N|⁻¹. The transition to bursting occurs at h/r* ≈ Δt/τ′, where "the dynamics can be considered close-to-critical".
  - In the bursting regime, the maximal autocorrelation time is ≃ τ′ and "the homeostatic timescale τ_hp influences the resulting dynamics".
  - Finite size: N = 10⁴ only. A finite network shows a peak in the avalanche-size distribution at large sizes. In the bursting regime, "system size and topology impact the network dynamics more strongly". There is no N-scan and no h ∝ 1/N analysis.
- **Kossio, Goedeke, van den Akker, Ibarz & Memmesheimer 2018**, "Growing critical: self-organized criticality in a developing neural system", Phys. Rev. Lett. 121:058301, arXiv:1811.02861. This is the closest prior model to the project: neurite disks grow, overlap area sets coupling, and growth stops when the rate reaches f_sat. VERIFIED from the ar5iv full text; the PRL number is inferred from the supplement link — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861); [arXiv](https://arxiv.org/abs/1811.02861)
  - The stationary branching parameter is σ = 1 − f0/f_sat (Eq. 3), the same structure as Zierenberg et al.'s m = 1 − h/r*.
  - Growth is assumed slow (1/K ≫ τ). Simulations use N = 100.
  - Fluctuations of Σ_j A_ij and deviations of f̄_i from f_sat are reported as small (<1%). There is no formal analysis of σ fluctuations, of the effect of growth rate K, or of finite-size scaling.
  - The avalanche size distribution is P(s) ∝ s^(−3/2) e^(−(σ − ln σ − 1)s) (Eq. 5), with cutoff s_c = (σ − ln σ − 1)⁻¹. The subcritical distribution is "independent of neuron number".

### Inferences
- The project's finding that m ≈ 0.965 at fixed h is constant in N is what m = 1 − h/r* predicts. Moving m toward 1 by setting h ∝ 1/N is, in this framework, simply lowering h/r*.
- A useful check: does 1 − m track h/r* exactly across N = 25–200? Going from 0.035 to 0.012 is roughly a 3× change for an 8× change in N. Any departure from h/r* could reflect coalescence (Q7) or estimator bias rather than a new effect.
- Zierenberg 2018's "saw-tooth" m(t) is the same phenomenon as the project's σ oscillation. They do not quantify its amplitude as a function of the plasticity step.
- Kossio 2018 is the most important prior for the "growing" framing. Its coupling is a continuous overlap area, so contact granularity is absent by construction.

### Gaps
- I did not check whether later Priesemann-group papers (e.g. Zierenberg, Wilting, Priesemann & Levina, "Tailored ensembles of neural networks optimize sensitivity to stimulus statistics", Phys. Rev. Research 2:013115, 2020; UNCERTAIN, not fetched) analyse dynamic range versus m in a way that bears on the project's dynamic-range finding.

## Q3. Is there an analytic treatment of how plasticity rate or timescale ratio sets the distance from criticality or the fluctuation amplitude?

### Takeaway
Partially. Kinouchi 2019 gives explicit 1/τ formulas for the fixed-point offset and the damping of the focus. Zierenberg 2018 gives a τ_hp stability threshold. Buendía 2020 gives qualitative h/ε phase logic. Menesse et al. 2022 show that timescale separation between coupling and threshold adaptation allows SOqC under input. Kuehn 2012 analysed timescale and noise optimality in adaptive SOC networks. I found no paper deriving fluctuation amplitude as a function of the discrete plasticity increment.

### Cited Findings
- **Buendía, di Santo, Bonachela & Muñoz 2020**, "Feedback mechanisms for self-organization to the edge of a phase transition", Front. Phys. 8:333, DOI 10.3389/fphy.2020.00333, arXiv:2006.03020. They state that in the limit ε, h → 0 (keeping timescale separation) "the system approaches true self-organization". The text "does not analytically relate the recovery/driving rate to distance from criticality" beyond the h/ε phase logic. VERIFIED (full text read; the DOI is constructed and UNCERTAIN) — [ar5iv](https://ar5iv.arxiv.org/html/2006.03020); [DOAJ record](https://doaj.org/article/a44e5390c46d455481610a35c1ba7d69)
- Related work by Buendía, di Santo, Villegas, Burioni & Muñoz, "Self-organized bistability and its possible relevance for brain dynamics" (Phys. Rev. Research 2020, arXiv:1911.05382). Per a search snippet, in the limit of infinitely slow synaptic dynamics SOC and SOB become identical, but the timescales needed "do not seem biologically plausible". UNCERTAIN (snippet only) — [arXiv](https://arxiv.org/abs/1911.05382v2)
- **Kinouchi et al. 2019**: explicit O(1/τ) offsets and damping (see Q1). The distance from criticality scales as 1/(uτ), where u is the depletion fraction per spike and τ the recovery time. VERIFIED — [ar5iv](https://ar5iv.labs.arxiv.org/html/1803.05537)
- **Menesse, Marin, Girardi-Schappo & Kinouchi 2022**, "Homeostatic criticality in neuronal networks", Chaos, Solitons & Fractals 156:111877 (published 2022), arXiv:2109.02718. With adequate timescale separation between coupling-strength and firing-threshold dynamics, near-criticality (SOqC) can be reached and sustained even under external input. Thresholds adapt to cancel the input ("h decreases towards zero"), similar to perfect adaptation. VERIFIED (abstract/snippets); the volume and article number are inferred from the RePEc URL and are UNCERTAIN — [arXiv](https://arxiv.org/abs/2109.02718v2); [ar5iv](https://ar5iv.labs.arxiv.org/html/2109.02718); [RePEc](https://ideas.repec.org/a/eee/chsofr/v156y2022ics0960077922000881.html)
- **Kuehn 2012**, "Time-scale and noise optimality in self-organized critical adaptive networks", arXiv:1110.5433 (Phys. Rev. E 85:026103 per my prior knowledge, UNCERTAIN). The paper argues that self-organization must act on a slower timescale than node dynamics, and analyses an optimal timescale ratio and noise level for adaptive-network SOC. Only the snippet was VERIFIED; I did not verify the details — [arXiv](https://arxiv.org/pdf/1110.5433)
- **Hesse & Gross 2014**, "Self-organized criticality as a fundamental property of neural systems", Front. Syst. Neurosci. 8:166, DOI 10.3389/fnsys.2014.00166. This is a review linking neural criticality to network SOC (e.g. Bornholdt & Rohlf 2000). I did not read the full text, so I could not confirm whether it treats rate or fluctuation precision. Citation VERIFIED — [Frontiers PDF](https://www.frontiersin.org/articles/10.3389/fnsys.2014.00166/pdf); [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4171833)
- **Droste, Do & Gross 2013**, "Analytical investigation of self-organized criticality in neural networks", J. R. Soc. Interface (arXiv:1203.4942). This is a related analytic adaptive-network SOC paper. Only the search hit was seen, so its content is UNCERTAIN — [arXiv](https://arxiv.org/pdf/1203.4942)
- **Gross 2021**, "Not one, but many critical states: a dynamical systems perspective", Front. Neural Circuits 15:614268. This is from prior knowledge only and was not fetched (UNCERTAIN).
- **Muñoz 2018**, "Colloquium: Criticality and dynamical scaling in living systems", Rev. Mod. Phys. 90:031001, DOI 10.1103/RevModPhys.90.031001. From prior knowledge, not fetched (UNCERTAIN).

### Inferences
- The existing analytic offsets scale with the plasticity rate (1/τ, or u/τ per spike), not with a quantum of coupling change. Treating "one contact changes σ by Δσ_c" as the controlling small parameter appears not to have been done. Such a treatment would be a random walk or limit cycle on a lattice of σ values, with feedback, compared against τ·rate.
- The project found that slow growth pins σ but gives curved distributions, while fast growth gives 3/2-like distributions. This connects to Buendía 2020's remark that "parameters (not only size)" set cleanliness of scaling, and to Kinouchi 2019's SO mechanism, where sweeping can mimic power laws. I found no paper linking these regimes to dynamic range.

### Gaps
- I could not verify the full contents of Gross 2021, Hesse & Gross 2014, Kuehn 2012 or Droste et al. 2013 within the tool budget.

## Q4. Has anyone studied the effect of synapse discreteness (integer synapse counts, quantized connectivity changes) on attaining criticality?

### Takeaway
I found no source that analyses synapse-number quantization as a limit on the precision of self-organized criticality. Related work uses synapse-number-based structural plasticity (e.g. Butz & van Ooyen MSP models; an eLife 2023 reviewed preprint) or continuous overlap coupling (Kossio 2018), but does not frame a granularity-driven precision limit. This gap is my assessment from a limited search; it is not a confirmed absence.

### Cited Findings
- A structural-plasticity review (2015, PMC4923203) frames homeostatic structural plasticity as regulation of synapse number. It contains no criticality-precision analysis — [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC4923203)
- An eLife reviewed preprint (88376) proposes a biphasic, synapse-number-based homeostatic structural plasticity rule. It is not about criticality precision; snippet only, UNCERTAIN — [eLife](https://elifesciences.org/reviewed-preprints/88376)
- Tetzlaff et al. 2010, "Self-organized criticality in developing neuronal networks", PLoS Comput. Biol., report that cultures move toward criticality during development. This is empirical and modelling work, with no granularity analysis (snippet only) — [PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2996321/)
- Kossio et al. 2018 use continuous overlap areas and report coupling fluctuations below 1%. They do not address discreteness — [ar5iv](https://ar5iv.arxiv.org/html/1811.02861)

### Inferences
- Within the sources checked, a "contact-granularity sets the σ fluctuation floor" result appears novel. The analogous continuous-variable result (O(1/τ) offset plus noise-driven oscillation) exists in Kinouchi 2019.

### Gaps
- I did not search the structural-plasticity criticality literature specifically for granularity, e.g. Butz, van Ooyen, Gallinaro & Rotter, or Tetzlaff's MSP-plus-criticality papers. A targeted search on "Gallinaro Rotter structural plasticity criticality" and "Butz van Ooyen avalanches" is recommended before claiming novelty.

## Q5. Finite-size scaling in the SOC limit (h → 0, N → ∞), the h ∝ 1/N scaling, and the Widom-line / quasicriticality picture

### Takeaway
The Ortiz–Beggs quasicriticality framework treats nonzero drive as moving the system along a nonequilibrium Widom line, where susceptibility peaks and exponents shrink while scaling relations approximately hold. Das & Levina 2019 is the explicit treatment of input scaling with N. LHG-type models scale the recovery time τ ∝ N to recover criticality. I found no paper that does a full finite-size analysis of homeostatic m with h ∝ 1/N.

### Cited Findings
- **Fosque, Williams-García, Beggs & Ortiz 2021**, "Evidence for quasicritical brain dynamics", Phys. Rev. Lett. 126:098101, DOI 10.1103/PhysRevLett.126.098101, arXiv:2010.02938. As stimulation increases, the cortex departs from criticality along a Widom line, defined as the line of peak dynamical susceptibility. Avalanche exponents and the branching ratio decrease while the dynamical scaling relation approximately holds. VERIFIED (abstract/snippets) — [arXiv](https://arxiv.org/abs/2010.02938v1); [IU ScholarWorks](https://scholarworks.indianapolis.iu.edu/items/f88b365e-b600-4ca0-b643-f0a87a0cfcd1)
- **Williams-García, Moore, Beggs & Ortiz 2014**, "Quasicritical brain dynamics on a nonequilibrium Widom line", Phys. Rev. E 90:062714. This is the original quasicriticality paper, from prior knowledge and not fetched (UNCERTAIN). It introduced the cortical branching model with drive, where peak susceptibility sits at branching ratio < 1 under nonzero drive.
- Kinouchi 2019 states that LHG recovers criticality in the thermodynamic limit by scaling τ = τ0·N. This is a "slow-rate ∝ N" scaling analogous to the project's h ∝ 1/N. VERIFIED — [ar5iv](https://ar5iv.labs.arxiv.org/html/1803.05537)
- Das & Levina 2019 define input regimes by their scaling with N (see Q6). VERIFIED — [ar5iv](https://ar5iv.labs.arxiv.org/html/1808.04196)

### Inferences
- In the quasicritical picture, the project's fixed-h result (m ≈ 0.965 for all N) is a point on the Widom line set by h. The N·h-fixed protocol moves the system along that line toward the critical point as N grows, which is the standard way to approach the SOC limit. This specific protocol for a homeostatic or growing network, with measured m(N) and exponents, was not found in the literature checked. Its conceptual basis, however (m = 1 − h/r*, τ ∝ N in LHG), is established.

### Gaps
- I did not verify whether Williams-García 2014 includes finite-size scaling of the Widom line.

## Q6. Das & Levina 2019 (PRX, relaxed timescale separation): what exactly do they show about exponents under nonzero drive?

### Takeaway
External input arriving during avalanches changes the asymptotic avalanche-size exponent from 3/2 to 5/4 in the "moderate" input regime. The 3/2 law survives below a crossover size of about φ⁻², and the distribution becomes Gaussian under strong input. Criticality, in the sense of a power law, does not require tuning φ.

### Cited Findings
- **Das & Levina 2019**, "Critical neuronal models with relaxed timescale separation", Phys. Rev. X 9:021062, DOI 10.1103/PhysRevX.9.021062, arXiv:1808.04196. All points below are VERIFIED from the full text — [ar5iv](https://ar5iv.labs.arxiv.org/html/1808.04196); [arXiv](https://arxiv.org/abs/1808.04196v1)
  - Levels model at critical M = N: with no input P(k) ~ k^(−1.5); with small input φ = o(1/N) still k^(−1.5); with moderate input φ = Θ(1), P(k) ~ k^(−1.25) (Theorem 6.7, bounded ratio); with large input φ = Θ(N), the distribution is Gaussian.
  - Below a crossover of about φ⁻² the 1.5 law persists; above it 1.25 holds. The authors call this multifractal-like.
  - The authors state "ϕ does not need to be externally tuned to achieve criticality".
  - Branching-network simulations (N = 10⁵) confirm the 1.5 → 1.25 change. Duration exponents also change.
  - They cite an experiment in which the exponent moved from 1.5 to about 1.3 during a task.
  - Internal inconsistency in how φ scales: the abstract says the "moderate input regime" is input scaling as 1/N, while the body uses φ = Θ(1). The likely reason is that φ is a network-level total, so the per-neuron rate scales as 1/N. That reading is my inference and is UNCERTAIN.
  - The paper does not use m or h and does not analyse homeostasis.

### Inferences
- If moderate input really means a per-neuron rate ∝ 1/N, then the project's N·h-fixed protocol sits in Das & Levina's moderate regime. They predict a 5/4 tail for large avalanches when avalanches are defined without separating merged activity. The project instead sees exponents approach critical-branching values. This could come from a different avalanche definition (e.g. single-seed tracking versus binned activity) or from small N (≤ 200), where the crossover at about φ⁻² may not be reached. This should be checked against the project's avalanche definition.

### Gaps
- Their companion PRX 2021 (Neto, Das, Levina et al., title not confirmed) argues that the branching ratio and susceptibility fail to establish criticality without timescale separation. I did not fetch it; only the snippet was seen.

## Q7. Zierenberg, Wilting, Priesemann & Levina 2020 (PRE, coalescence): microscopic versus macroscopic branching

### Takeaway
In finite networks, coalescence (one node activated by several sources at once) makes the macroscopic branching-parameter estimate differ from the microscopic model parameter. The authors derive a universal nonlinear scaling function and a corrected estimator valid for all system sizes.

### Cited Findings
- **Zierenberg, Wilting, Priesemann & Levina 2020**, "Description of spreading dynamics by microscopic network models and macroscopic branching processes can differ due to coalescence", Phys. Rev. E 101:022301, DOI 10.1103/PhysRevE.101.022301, arXiv:1905.10402. Coalescence underlies the discrepancy. A universal nonlinear scaling function for the conditional expectation of successive activity quantifies the gap, and a nonlinear estimator recovers the microscopic parameter "for all system sizes". VERIFIED (abstract) — [arXiv](https://arxiv.org/pdf/1905.10402); [MPI record](https://is.mpg.de/publications/item_3214532)

### Inferences
- With N = 25–200, coalescence is likely large. The project's measured m (0.965 → 0.988) may be biased relative to the microscopic σ. Part of the N-dependence under N·h fixed could reflect reduced coalescence and not only the reduced h/r*. The project should report both the σ from contacts and the m from regression, and should compare against the Zierenberg 2020 estimator.

### Gaps
- I did not read the exact functional form of the scaling function, so I cannot say how large the bias is at N = 25.

---

## Overall assessment: open versus answered (for the "precision limit" analysis)

| Project question | Status | Basis |
|---|---|---|
| Homeostasis sets the mean distance from criticality via drive, 1 − m = h/r* | **Answered** | Zierenberg 2018 Eq. 7 (VERIFIED); Kossio 2018 σ = 1 − f0/f_sat for growing networks (VERIFIED) |
| The homeostatic variable oscillates and hovers, never sits exactly at the critical point; super/subcritical epochs and dragon kings | **Answered qualitatively** | Bonachela & Muñoz 2009 / Bonachela 2010 (SOqC); Kinouchi 2019 (stochastic oscillations, dragon kings); Zierenberg 2018 saw-tooth m(t) |
| Fixed-point offset and damping versus plasticity timescale | **Partly answered** (O(1/τ) for continuous adaptive gain or synapses) | Kinouchi 2019 (VERIFIED); Zierenberg 2018 τ_hp stability bound |
| Fluctuation amplitude versus N | **Partly answered and contested.** Bonachela 2010 claims amplitude is roughly independent of N (LHG); Kinouchi 2019 shows a plateau until N(τ), then a decrease, and that it vanishes for τ ∝ N | Kinouchi 2019 (VERIFIED); Bonachela 2010 (secondary) |
| Effect of nonzero drive on exponents | **Answered** for drive during avalanches (3/2 → 5/4, crossover about φ⁻²) and for Widom-line quasicriticality (exponents shrink) | Das & Levina 2019 (VERIFIED); Fosque 2021 (abstract VERIFIED) |
| Estimator bias of m at small N | **Answered in principle** (coalescence) | Zierenberg 2020 (abstract VERIFIED) |
| σ fluctuation floor set by **contact granularity** (quantized coupling changes) | **Appears open.** No source found | Limited search; see Q4 gaps |
| Growth-rate trade-off: slow growth gives pinned σ, curved distributions and best dynamic range; fast growth gives 3/2-like distributions and worse dynamic range | **Appears open.** Related qualitative points exist ("dirty" scaling from sweeping; parameters set cleanliness) but there is no analysis linking the plasticity rate or step to both avalanche shape and dynamic range | Buendía 2020; Kinouchi 2019 |
| Systematic finite-size scaling of a homeostatic or growing network under h ∝ 1/N (m(N), exponents(N)) | **Appears largely open as a concrete study.** The conceptual pieces exist (m = 1 − h/r*, LHG τ ∝ N, Das & Levina's input-scaling regimes), but none does this for structural or homeostatic plasticity, and Das & Levina's 5/4 prediction for moderate input is a possible tension to address | Zierenberg 2018; Das & Levina 2019; Kinouchi 2019 |

Novelty verdict (inference): a "precision-limit" paper would be novel if it is framed around (i) quantized coupling increments as the source of the σ fluctuation floor, (ii) a joint scaling law combining granularity Δσ_c, plasticity rate and N, and (iii) the trade-off between avalanche-shape cleanliness and dynamic range. It must explicitly position itself against Kinouchi 2019 (O(1/τ) offsets, N(τ) plateau), Zierenberg 2018 (m = 1 − h/r*, τ_hp bound), Kossio 2018 (growing networks reach σ = 1 − f0/f_sat), Bonachela 2009/2010 (SOqC), Das & Levina 2019 (drive changes exponents) and Zierenberg 2020 (estimator bias). The mean-offset result alone (1 − m ∝ h) is not novel.
