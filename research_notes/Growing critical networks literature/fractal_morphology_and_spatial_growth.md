# Fractal Morphology of Neurons, Spatial Network-Growth Models, and Links to Network Dynamics

Verification tags:
- **[V]**: confirmed this session from the publisher, PubMed, arXiv or a repository page (the existence of the paper and the stated claim).
- **[V-exists]**: the paper's existence and bibliographic data were confirmed this session, but the specific numeric or content claim comes from prior knowledge and was not re-read.
- **[PK]**: from prior knowledge only. The DOI follows the usual publisher pattern but was not re-checked this session. Check these before quoting.

Search budget: about 16 tool calls. PubMed Central blocked automated fetching, so several full texts could not be read.

---

## Q1. Measured fractal dimensions (D) of real dendrites and neurons: typical values, methods and developmental change

### Takeaway
Box-counting D of 2D projected neuronal arbors usually falls between about 1.1 and 1.7. Sparse or simple cells sit near the bottom (about 1.1–1.3). Cortical and hippocampal pyramidal arbors fall around 1.3–1.5. Dense retinal ganglion cells and cultured neurons reach about 1.6–1.8, in the range predicted by diffusion-limited aggregation (DLA). These values hold over only about 1–1.5 decades of scale, depend heavily on method (box size range, skeleton versus filled image, projection), and are disputed: some authors argue that retinal arbors are space-filling (D ≈ 2 down to a cutoff) rather than fractal. Individual dendrite branches are barely fractal (D ≈ 1.03). The "fractality" belongs to the arbor as a whole.

### Cited Findings
**Foundational methods papers**
- Smith, Marks, Lange, Sheriff & Neale (1989) introduced box-counting and related fractal analysis for 2D images of cultured neurons and other cells. "A fractal analysis of cell images", *J Neurosci Methods* 27:173–180. **[PK]** https://doi.org/10.1016/0165-0270(89)90099-0
- Smith, Lange & Marks (1996) reviewed fractal methods (dimension, lacunarity, multifractals) for cell morphology. They stress that biological fractals hold only over a limited range and that results depend on method. "Fractal methods and results in cellular morphology — dimensions, lacunarity and multifractals", *J Neurosci Methods* 69:123–136. **[PK]** https://doi.org/10.1016/S0165-0270(96)00080-5
- Neale, Bowers & Smith (1993) quantified early dendrite development in spinal-cord cultures, including fractal measures. "Early dendrite development in spinal cord cell cultures: a quantitative study", *J Neurosci Res* 34:54–66. **[V-exists]** https://scilit.com/publications/cd9510afba02a6c4fc8ba492d9914913. My recollection, not verified: D rises over the first days in vitro as arbors elaborate.

**Retinal neurons (Caserta, Eldred, Fernández, Stanley)**
- Caserta, Stanley, Eldred, Daccord, Hausman & Nittmann (1990) studied chick-embryo retinal neurons cultured in media of different viscosity. They reported D consistent with DLA, which the Werner review summarises as "approximately 1.7". "Physical mechanisms underlying neurite outgrowth: a quantitative analysis of neuronal shape", *Phys Rev Lett* 64:95–98. **[V-exists]** https://doi.org/10.1103/PhysRevLett.64.95. Value summarised in Werner 2010, *Front Physiol* 1:15 **[V]**: https://arxiv.org/pdf/0910.2741. The often-quoted exact figure "1.68" could **not** be verified this session.
- Kleinfeld et al. published a critical Comment on Caserta et al. (1990) in *PRL*, disputing the interpretation. **[V-exists]** https://ouci.dntb.gov.ua/works/4E3Vz629
- Caserta, Eldred, Fernández, Hausman, Stanford, Buldyrev, Schwarzer & Stanley (1995) measured D of physiologically characterised cat retinal ganglion cells in 2D and 3D, comparing several methods. "Determination of fractal dimension of physiologically characterized neurons in two and three dimensions", *J Neurosci Methods* 56:133–144. **[V-exists]** https://acnpsearch.tweb-dev.unibo.it/singlejournalindex/27250595 ; **[PK]** DOI https://doi.org/10.1016/0165-0270(94)00115-W
- Fernández, Eldred, Ammermüller, Block, von Bloh & Kolb (1994) found that D differs systematically between turtle-retina cell classes (amacrine, ganglion, horizontal, bipolar). "Complexity and scaling properties of amacrine, ganglion, horizontal, and bipolar cells in the turtle retina", *J Comp Neurol* 347:397–408. **[V-exists]** https://pik-potsdam.de/~bloh/abstracts/fernandez94.html
- A search snippet reported a retinal ganglion cell Hausdorff dimension of H = 1.73 ± 0.12 (n = 18, range 1.54–1.94). The exact source paper could not be identified, so treat this as **uncertain**.
- Montague & Friedlander (1989) studied isolated mammalian retinal neurons in culture and argued for an "intrinsic growth strategy". "Expression of an intrinsic growth strategy by mammalian retinal neurons", *PNAS* 86:7223–7227. **[V-exists]** https://beam.vt.edu/content/dam/fbri_vtc_vt_edu/publications/montague-publications/MontagueFriedlander1989PNAS.pdf. My recollection, unverified: D stayed roughly constant as cells grew, meaning scale-invariant growth.

**Skeptical counterpoint**
- Panico & Sterling (1995) compared retinal neuron and vessel arbors with synthetic fractal and non-fractal patterns. They concluded the biological patterns "are not fractal": they show D → 2 (space-filling down to a cutoff) and are better described as "quasi-regular lattices". Non-fractal growth algorithms produced patterns that looked as fractal as the real ones. "Retinal neurons and vessels are not fractal but space-filling", *J Comp Neurol* 361:479–490. **[V]** https://doi.org/10.1002/cne.903610311 (abstract: https://scilit.com/publications/e7e486a91c6926cff6b90cf4d5788162)

**Reviews and methodological caveats**
- Jelinek & Fernández (1998) reviewed how reliable and useful D is for neurons. Key cautions: D depends on image resolution, the box-size range, projection of 3D cells to 2D, thresholding and skeletonisation, and the choice of method (box-counting, mass-radius, dilation). D is a useful descriptor for classifying cells, but not proof of a fractal generating process. "Neurons and fractals: how reliable and useful are calculations of fractal dimensions?", *J Neurosci Methods* 81:9–18. **[PK]** https://doi.org/10.1016/S0165-0270(98)00021-1
- Fernández & Jelinek (2001) cover the same points with more detail. "Use of fractal theory in neuroscience: methods, advantages, and potential problems", *Methods* 24:309–321. **[PK]** https://doi.org/10.1006/meth.2001.1201
- Werner (2010) states that different methods give noticeably different D for the same trees. *Front Physiol* 1:15. **[V]** https://arxiv.org/pdf/0910.2741

**Milošević, Ristanović et al. (Belgrade group)**
- In a series of papers from about 2005 to 2012, this group applied box-counting (including "modified" box-counting) to Golgi-impregnated spinal dorsal-horn and other neurons. They used D to separate cell classes and laminae, reporting 2D values typically around 1.2–1.5. Example: Milošević, Ristanović & Stanković (2005), "Fractal analysis of the laminar organization of spinal cord neurons", *J Neurosci Methods* 146:198–204. **[PK]** https://doi.org/10.1016/j.jneumeth.2005.02.009. Treat the per-type values as **uncertain**; they were not re-read.

**Karperien & Jelinek (FracLac)**
- FracLac is the ImageJ box-counting plugin most often used for 2D cell D and lacunarity. It was applied extensively to microglia, where D runs from about 1.2 for amoeboid cells to about 1.5+ for ramified cells. Karperien, Ahammer & Jelinek (2013), "Quantitating the subtleties of microglial morphology with fractal analysis", *Front Cell Neurosci* 7:3. **[PK]** https://doi.org/10.3389/fncel.2013.00003
- Karperien & Jelinek (2016), "Box-counting fractal analysis: a primer for the clinician", in Di Ieva (ed.), *The Fractal Geometry of the Brain*, Springer. Covers grid-position sensitivity, the need for multiple grid placements, and choosing the scaling range. **[PK]** https://doi.org/10.1007/978-1-4939-3995-4_2

**Taylor group (see Q2)**
- Individual CA1 basal dendrites have a mean coastline D of about 1.03, close to a straight line. Whole arbors have D_A ≈ 1.42–1.46. Scaling holds over only about one order of magnitude (4–40 µm), about 1.5 at most. Rowland et al. (2023), *Front Netw Physiol*. **[V]** https://doi.org/10.3389/fnetp.2023.1072815

### Inferences
- A realistic 2D target is roughly **D ≈ 1.3–1.7**. D = 2.5 cannot be reached in a 2D embedding: box-counting D of a 2D image is bounded by 2. Values near 1.3 match sparse cells; 1.6–1.7 match dense retinal or cultured cells.
- Measured D depends on method as much as on biology. Any claim that the simulation "has D = x" should be checked by box-counting the generated images over the same 1–1.5 decade range used experimentally, with several grid offsets (FracLac-style), not inferred from the generative parameter.
- Developmental trends are poorly established. Cultured neurons appear to increase D over the first days in vitro (Smith/Neale), while Montague & Friedlander suggest scale-invariant growth with roughly constant D. Both claims need checking against the full texts.

### Gaps
- Exact per-cell-type D tables (Smith 1989; Caserta 1995; Milošević/Ristanović) could not be extracted because full texts were not accessible this session.
- The often-quoted Caserta et al. (1990) value of 1.68 is unverified.
- No systematic modern (post-2015) dataset applying box-counting D across NeuroMorpho.org cell types in 2D projection was found in this session's searches. One may exist.

---

## Q2. Richard Taylor's fractal-neuroscience work: neuron D, connectivity versus cost, an "optimal" D, and fractal electrodes

### Takeaway
Taylor's group (University of Oregon, with Dalrymple-Alford and Harland in Christchurch, NZ) analysed 3D reconstructions of rat CA1 pyramidal neurons and reported arbor D ≈ 1.4–1.5. They distorted the arbors (branch length, forking angles and "weaving", i.e. tortuosity) and argued that natural neurons sit near a peak of a geometric connectivity-to-cost ratio. "Connectivity" in this work is a purely geometric proximity measure on single neurons. It is **not** network dynamics, and the "optimum" is a ratio of their own design-specific metrics.

### Cited Findings
- Smith, Rowland, Harland, Moslehi, Montgomery, Schobert, Watterson, Dalrymple-Alford & Taylor (2021), "How neurons exploit fractal geometry to optimize their network connectivity", *Sci Rep* 11:2332. **[V]** https://doi.org/10.1038/s41598-021-81421-2
  - Built 3D models of rat neurons from confocal stacks and used D to summarise length, forking and weaving together. They tested about 1,600 distorted neuron models and concluded that neuron D reflects a balance in which connectivity outweighs cost for neurons with high D. **[V]**
  - Taylor's framing: unlike trees, neurons create their fractal character not only through branch-size scaling but through how branches "weave" through space. **[V]** (press coverage: https://www.sciencedaily.com/releases/2021/03/210302150116.htm)
  - The exact D ranges and the definition of the connectivity metric could not be read; PMC blocked the fetch. From prior knowledge **[PK]**, the metric counts potential contacts (proximity of arbor points to a surrounding volume, or "reach") per unit material.
- Rowland, Harland, Smith, Moslehi, Dalrymple-Alford & Taylor (2022), "Investigating fractal analysis as a diagnostic tool that probes the connectivity of hippocampal neurons", *Front Physiol* 13:932598. **[V]** https://doi.org/10.3389/fphys.2022.932598
  - Anterior thalamic lesion and enrichment study. Enriched housing plus spatial training raised arbor D and connectivity in sham rats. Lesion and non-lesion groups both sat at what the authors call an "optimal balance of connectivity against material cost". **[V]**
- Rowland, Smith, Moslehi, Harland, Dalrymple-Alford & Taylor (2023), "Neuron arbor geometry is sensitive to the limited-range fractal properties of their dendrites", *Front Netw Physiol* 3:1072815. **[V]** https://doi.org/10.3389/fnetp.2023.1072815
  - About 105 basal arbors from 33 rats. Individual dendrites ⟨D_B⟩ ≈ 1.031 by the coastline method, with a tortuosity method agreeing. Whole-arbor D_A ≈ 1.42–1.46. **[V]**
  - D_A is very sensitive to D_B: ⟨D_BC⟩ falling from 1.02 to 1.01 raised D_A from 1.42 to 1.46, roughly 4× amplification. **[V]**
  - Connectivity/cost ratios peak near the natural D_A and fall off on either side. The authors predict a generic U-shaped D_B–D_A relation for self-avoiding fractal neurons. **[V]**
  - Caveats the authors acknowledge: scaling spans only about 1 decade (4–40 µm), so D is an "effective" dimension. Claiming fractality over such a narrow range is controversial. Results cover one cell type, sex and strain. Clinical and implant applications are speculative. **[V]**
- Fractal electrodes and retinal implants:
  - Watterson, Montgomery & Taylor (2017), "Fractal electrodes as a generic interface for stimulating neurons", *Sci Rep* 7:6717. **[PK]** https://doi.org/10.1038/s41598-017-06762-3
  - Moslehi et al. (2022), "Controlled assembly of retinal cells on fractal and Euclidean electrodes", *PLoS One* 17:e0265685. **[PK]** https://doi.org/10.1371/journal.pone.0265685
  - The Smith 2021 press material explicitly links the neuron-D work to fractal-electrode design for retinal implants. **[V]** (ScienceDaily link above)
- The commonly repeated claim that "neurons have D ≈ 1.4–1.6" matches Taylor's CA1 values (about 1.42–1.46) **[V]** and the broad 2D literature range in Q1. The specific upper bound of 1.6 was not found stated as a Taylor-group result.

### Inferences
- Taylor's "optimal D for connectivity" concerns **single-neuron geometric reach per unit wiring**. It does not show that network dynamics (criticality, dynamic range) are optimised at a particular D. Using it to justify a D for network-level studies is an extrapolation.
- Taylor's result that weaving or tortuosity, not branch-length scaling, contributes much of D_A matters for the simulation. A model that sets D only through an inter-branch length ratio r = b^(−1/D) leaves out a major contributor to measured D. Tortuosity and branching angles could matter more than r.

### Gaps
- The exact definition and formula of the connectivity and cost metrics (R_PA, R_PV) and the numeric "optimal D" in Smith et al. (2021) were not verified, because the full text was inaccessible.
- No independent replication of the "optimal D" claim by another group was found.
- No Taylor-group paper connecting neuron D to network-level *dynamics* was found.

---

## Q3. Wiring-optimisation and morphology-generation theory: what it says about space-filling and scaling

### Takeaway
The dominant theory treats dendrites as minimum-spanning-tree-like structures that connect a set of targets (synapse locations) under a trade-off between total cable length and path length to the soma. It yields **power-law scaling between total length, branch-point number and arbor volume or area (Cuntz)**, constant branch density or "space-filling to a cutoff" (Teeter & Stevens, Snider et al.), and size-invariant excitability (dendritic constancy). These are *space-filling plus scaling-law* accounts, not fractal-dimension accounts. They imply that arbor **extent, target density and the balancing factor** set morphology, and that "D" is a derived descriptor.

### Cited Findings
- Cuntz, Forstner, Borst & Häusser (2010), "One rule to grow them all: a general theory of neuronal branching and its practical application", *PLoS Comput Biol* 6:e1000877. **[V-exists]** https://doi.org/10.1371/journal.pcbi.1000877 (PubMed: https://pubmed.ncbi.nlm.nih.gov/20700495/)
  - **[PK]** Arbors are generated by a greedy minimum-spanning-tree algorithm on carrier points. The cost function is total length + bf × path length to the root, where the balancing factor bf ≈ 0.2–0.7 depends on cell type. Point density within the arbor's spanning field determines shape. This reproduces many cell types.
- Cuntz, Forstner, Borst & Häusser (2011), "The TREES toolbox—probing the basis of axonal and dendritic branching", *Neuroinformatics* 9:91–96. A Matlab toolbox (with a later Python port) for generating, analysing and editing trees. **[PK]** https://doi.org/10.1007/s12021-010-9093-7
- Cuntz, Mathy & Häusser (2012), "A scaling law derived from optimal dendritic wiring", *PNAS* 109:11014–11018. **[PK]** https://doi.org/10.1073/pnas.1200430109
  - **[PK]** For optimally wired trees, total length L ∝ V^{1/3} · n^{2/3} in 3D, where n is the number of branch points or targets and V the spanning volume. They report this holds across many cell types.
  - **Uncertain:** the 2D analogue, L ∝ (A·n)^{1/2}, is my inference from the same geometry and was not checked.
- Cuntz, Bird, Mittag, Beining, Schneider, Mediavilla, Hoffmann, Deller & Jedlicka (2021), "A general principle of dendritic constancy: a neuron's size- and shape-invariant excitability", *Neuron* 109:3647–3662. **[PK]** https://doi.org/10.1016/j.neuron.2021.08.028
  - **[PK]** Firing response to distributed synaptic input is invariant to dendrite size and shape, because input conductance scales with total length and synapse number together. Morphology changes are therefore partly compensated at the level of single-neuron excitability.
- Chklovskii (2004), "Synaptic connectivity and neuronal morphology: two sides of the same coin", *Neuron* 43:609–617. **[PK]** https://doi.org/10.1016/j.neuron.2004.08.012
  - **[PK]** Argues that arbor shapes reflect the requirement to make specific connections at minimal wiring cost. Introduces the "potential synapse" or filling-fraction framework: connectivity is set by arbor density overlaps within a spine-reach distance.
- Wen, Stepanyants, Elston, Grosberg & Chklovskii (2009), "Maximization of the connectivity repertoire as a statistical principle governing the shapes of dendritic arbors", *PNAS* 106:12536–12541. **[V-exists]** Analysed basal arbors of 2,171 pyramidal neurons and reported universal scaling and self-similarity. https://doi.org/10.1073/pnas.0901530106 (PDF: https://neurophysics.ucsd.edu/courses/physics_171/Maximization%20of%20the%20connectivity%20repertoire%20as%20a%20statistical%20principle%20governing%20the%20shapes%20of%20dendritic%20arbors.pdf)
- Snider, Pillai & Stevens (2010), "A universal property of axonal and dendritic arbors", *Neuron* 66:45–56. **[V-exists]** https://doi.org/10.1016/j.neuron.2010.03.003 (ModelDB: https://modeldb.science/citations/148914)
  - Analysed axonal and dendritic arbors across cell types and species and found shared structural features suggesting universal growth principles. **[V]** (secondary description only)
  - **Correction to the brief:** this paper appeared in *Neuron*, not *PNAS*. The *PNAS* "universal" paper of that era is Wen et al. (2009) above. No *PNAS* "Universal statistics of branched networks" by Snider et al. was found.
- Teeter & Stevens (2011), "A general principle of neural arbor branch density", *Curr Biol* 21:2105–2108. **[PK]** https://doi.org/10.1016/j.cub.2011.11.013
  - **[PK]** Branch density falls off with distance from the soma, following a universal form after rescaling. Arbors are close to uniformly space-filling within their territory.
- van Pelt & Uylings branching models (the BES model). Branching probability depends on time (B, E parameters) and on the number of terminal segments (S, order dependence). The model reproduces segment-number distributions and topological asymmetry.
  - van Pelt & Uylings (2002), *Network* 13:261–281. **[PK]** https://doi.org/10.1088/0954-898X_13_3_302
  - van Pelt, Dityatev & Uylings (1997), *J Comp Neurol* 387:325–340. **[PK]** https://doi.org/10.1002/(SICI)1096-9861(19971027)387:3<325::AID-CNE1>3.0.CO;2-2

### Inferences
- Wiring theory says that arbor *territory* (span) and *target density* control total length and contact probability. Branch-length scaling exponents are secondary. In contact-based connectivity, two arbors connect in proportion to the overlap of their length densities (Peters'-rule / filling-fraction logic). This depends on density profiles integrated over the overlap region, not directly on D.
- The constancy principle suggests that even when morphology changes, neurons may be close to response-invariant. A model with fixed per-synapse weights would *not* show this compensation, which is a modelling choice worth making explicit.
- A more theory-grounded way to make "D" matter: generate arbors with a TREES-like MST + balancing-factor rule, with target-point density as the control knob, rather than a direct length-ratio rule.

### Gaps
- Numeric bf ranges per cell type and the exact Cuntz 2012 equation constants were not re-verified.
- No paper was found that derives box-counting D analytically from the MST/bf model. The link between bf and D is apparently empirical or unexplored.

---

## Q4. Spatial network-growth simulators and culture growth models (classic and 2019–2026)

### Takeaway
Several mature growth simulators exist: NETMORPH, CX3D, BioDynaMo, DeNSE and NeuroDevSim. Culture models (Orlandi et al. 2013; Okujeni & Egert 2019; Gritsun et al. 2012; Houben et al. 2025) show that **spatial embedding plus axon length and anisotropy** shape burst nucleation, propagation and modularity. **Homeostatic growth models (Abbott & Rohrkemper 2007; Tetzlaff et al. 2010; Kossio et al. 2018) explicitly "grow" networks into criticality.** In those models criticality is set by the activity-dependent growth rule, not by arbor geometry.

### Cited Findings
**General simulators**
- NETMORPH (Koene, Tijms, van Hees, Postma, de Ridder, Ramakers, van Pelt & van Ooyen 2009, *Neuroinformatics* 7:195–210). Stochastic, phenomenological growth of large 2D/3D networks with realistic morphology (based on the van Pelt BES model). Synapses are formed where axonal and dendritic segments come within a distance criterion. **[V-exists]** https://doi.org/10.1007/s12021-009-9052-3 (PubMed: https://pubmed.ncbi.nlm.nih.gov/19672726/)
- CX3D (Zubler & Douglas 2009, "A framework for modeling the growth and development of neurons and networks", *Front Comput Neurosci* 3:25). Mechanistic 3D physics with diffusing chemical cues. **[PK]** https://doi.org/10.3389/neuro.10.025.2009
- Aćimović, Mäki-Marttunen, Havela, Teppola & Linne (2011), "Modeling of neuronal growth in vitro: comparison of simulation tools NETMORPH and CX3D", *EURASIP J Bioinform Syst Biol*. **[V-exists]** https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3171303/. **[PK]** They found the resulting network statistics depend strongly on growth parameters and differ between tools.
- BioDynaMo (Breitwieser et al. 2022, "BioDynaMo: a modular platform for high-performance agent-based simulation", *Bioinformatics* 38:453–460). Agent-based simulator that includes neurite growth modules descended from Cx3D. **[PK]** https://doi.org/10.1093/bioinformatics/btab649
- DeNSE (Fardet, Quaresima & Bottani). "A generic platform to investigate growth and network formations in 2D ensembles of neurons". C++ core with a Python interface. Integrates with NEST via NNGT. **[V]** https://github.com/SENeC-Initiative/DeNSE ; docs: https://dense.readthedocs.io/en/latest/user/tutorial.html ; abstract: https://bmcneurosci.biomedcentral.com/articles/10.1186/s12868-019-0538-0/figures/19 (venue details **uncertain**: it appears as a 2019 BMC Neurosci conference supplement plus a preprint). This is the closest existing tool to a "2D growth → contact connectivity → NEST dynamics" pipeline.
- NeuroDevSim. "Efficient simulation of neural development using shared memory parallelization", bioRxiv 2022.10.17.512465. **[V-exists]** https://www.biorxiv.org/content/10.1101/2022.10.17.512465.full.pdf. Author (E. De Schutter) and journal version (*Front Neuroinform* 2023) are **[PK]**, not verified.

**Culture models linking growth or space to collective dynamics**
- Orlandi, Soriano, Alvarez-Lacalle, Teller & Casademunt (2013), "Noise focusing and the emergence of coherent activity in neuronal cultures", *Nat Phys* 9:582–590. **[PK]** https://doi.org/10.1038/nphys2686
  - **[PK]** Grew 2D model cultures with axons as persistent random walks of finite length. Bursts nucleate at a few sites where spontaneous noise is amplified ("noise focusing") by the spatial connectivity structure. Nucleation sites depend on the metric (spatial) organisation, not only on degree.
- Okujeni & Egert (2019), "Self-organization of modular network architecture by activity-dependent neuronal migration and outgrowth", *eLife* 8:e47996. **[PK]** https://doi.org/10.7554/eLife.47996
  - **[PK]** Activity-dependent cell migration and neurite outgrowth produce clustered, modular networks. Modularity changes burst initiation and propagation.
- Gritsun, le Feber & Rutten (2012), "Growth dynamics explain the development of spatiotemporal burst activity of young cultured neuronal networks in detail", *PLoS One* 7:e43352. **[V-exists]** https://doi.org/10.1371/journal.pone.0043352
- Houben, Garcia-Ojalvo & Soriano (2025), "Role of connectivity anisotropies in the dynamics of cultured neuronal networks", arXiv:2501.04427. Uses an Izhikevich network on grown connectivity with imprinted anisotropies, using mean axon length as a developmental proxy. Activity patterns depend strongly on connectivity details and noise, and anisotropy improves reconstruction of structure from activity. **[V]** https://arxiv.org/abs/2501.04427. Snippet-level claim **[uncertain]**: rich dynamics appear only at intermediate axon length.
- Percolation-type culture model ("Simple model of complex dynamics of activity patterns in developing networks of neuronal cultures", *PLoS One* 2019): connectivity acts as a bifurcation parameter (silence → bursts → whole-network activation). **[V-exists]** https://doi.org/10.1371/journal.pone.0218304 ; https://arxiv.org/abs/1812.09611
- 2026 items found but **unverified**:
  - "Spatiotemporal bursting in simulated cultures of cortical neurons" (Univ. Washington, *BioSystems* 2026?). The snippet claims burst development needs only a simple growth model without fine-tuning. https://depts.washington.edu/biocomp/pubs/biosystems-2026.pdf
  - "Regime-dependent reconstruction of structural connectivity: an in silico study of spatially constrained modular neuronal networks", arXiv:2611.25830. Uses Rayleigh-distributed axon lengths of 0.3–3 mm. https://arxiv.org/pdf/2611.25830
  - Neither full text could be read.

**Growth models that self-organise to criticality**
- Abbott & Rohrkemper (2007), "A simple growth model constructs critical avalanche networks", *Prog Brain Res* 165:13–19. Neurons grow or shrink circular neurite fields according to activity (homeostasis). The network settles at a point where avalanches are power-law distributed. **[PK]** https://doi.org/10.1016/S0079-6123(06)65002-4
- Tetzlaff, Okujeni, Egert, Wörgötter & Butz (2010), "Self-organized criticality in developing neuronal networks", *PLoS Comput Biol* 6:e1001013. Activity-dependent neurite growth driving firing-rate homeostasis yields critical avalanches during development, compared against culture recordings. **[V]** https://doi.org/10.1371/journal.pcbi.1001013
- Kalle Kossio, Goedeke, van den Akker, Ibarz & Memmesheimer (2018), "Growing critical: self-organized criticality in a developing neural system", *Phys Rev Lett* 121:058301. Spiking network with spatially growing neurite fields. Growth uses only local information and drives the network to criticality. Avalanches correspond to Hawkes-process clusters, with size and duration distributions derived analytically that match experiment. **[V]** https://arxiv.org/abs/1811.02861 ; journal DOI **[PK]** https://doi.org/10.1103/PhysRevLett.121.058301

### Inferences
- In all three "growing critical" models, the critical point is reached because **homeostatic growth tunes total effective coupling** (overlap × weight) until the branching ratio is about 1. Geometry sets only *how* that coupling is distributed. If the simulation includes any homeostatic, normalisation or weight-tuning step, or tunes weights to criticality after growth, D will be largely factored out. This could explain the null result.
- The culture literature points to **connection range (axon length), anisotropy, cell clustering or modularity and noise** as the spatial features that change dynamics: burst nucleation, propagation and pattern repertoire. Dendritic D is not among them.
- DeNSE (2D, NEST-integrated) and the Orlandi-style axon random-walk models are the closest published analogues to the user's setup and are natural baselines.

### Gaps
- No simulator paper was found that sweeps target fractal dimension and reports the effect on network dynamics.
- Most culture models represent dendrites as simple disks or short trees and give the axon the active role. Contact connectivity between explicitly branched dendritic *and* axonal arbors with controllable D seems largely unexplored.

---

## Q5. Has any study linked arbor fractal dimension or morphology to network criticality, avalanches, dynamic range, synchronisation or information capacity?

### Takeaway
**No study was found that directly links dendritic or axonal box-counting D to network-level criticality or avalanche statistics.** The adjacent literatures are:
- **Single-neuron dendritic criticality and dynamic range**: active dendritic trees as excitable branching media. Larger trees give larger dynamic range.
- **Network topology and criticality**: hierarchical-modular structure produces Griffiths phases, stretched critical regions and limited sustained activity.
- **Spatial embedding and avalanche universality**: short-range 2D coupling gives non-mean-field exponents.
- **Morphology-derived connectomes**: arbor shape determines wiring statistics.

The gap is real and may be a publishable niche. It also suggests D is a weak control parameter unless it changes degree heterogeneity, connection range or modularity.

### Cited Findings
**Single-neuron (dendritic) criticality and dynamic range**
- Gollo, Kinouchi & Copelli (2009), "Active dendrites enhance neuronal dynamic range", *PLoS Comput Biol* 5:e1000402. **[V]** https://doi.org/10.1371/journal.pcbi.1000402
  - Model active dendritic trees (excitable branching compartments) give strongly nonlinear input-output curves with dynamic range above 50 dB.
  - Predicts that larger dendritic trees give larger dynamic range, and that blocking active conductances reduces it.
  - Reproduces double-sigmoid responses seen in retinal ganglion cells.
- Gollo, Kinouchi & Copelli (2013), "Single-neuron criticality optimizes analog dendritic computation", *Sci Rep* 3:3222. **[V-exists]** https://arxiv.org/abs/1304.4676 ; DOI **[PK]** https://doi.org/10.1038/srep03222. **[PK]** Dynamic range of the dendritic tree peaks at a critical coupling between branches. This is criticality *within* the arbor; tree topology and size, not D as such, set the response.
- Kinouchi & Copelli (2006), "Optimal dynamical range of excitable networks at criticality", *Nat Phys* 2:348–351. **[PK]** https://doi.org/10.1038/nphys289
- Shew, Yang, Petermann, Roy & Plenz (2009), "Neuronal avalanches imply maximum dynamic range in cortical networks at criticality", *J Neurosci* 29:15595–15600. **[V-exists]** https://arxiv.org/abs/0906.0527 ; DOI **[PK]** https://doi.org/10.1523/JNEUROSCI.3864-09.2009
- Gautam, Hoang, McClanahan, Grady & Shew (2015), "Maximizing sensory dynamic range by tuning the cortical state to criticality", *PLoS Comput Biol* 11:e1004576. **[PK]** https://doi.org/10.1371/journal.pcbi.1004576

**Network topology, hierarchy and spatial embedding**
- Kaiser & Hilgetag (2010), "Optimal hierarchical modular topologies for producing limited sustained activation of neural networks", *Front Neuroinform* 4:8. Hierarchical modular networks sustain activity in a limited, non-explosive regime over a broader parameter range than random or small-world networks. **[PK]** https://doi.org/10.3389/fninf.2010.00008
- Kaiser & Hilgetag (2004), "Spatial growth of real-world networks", *Phys Rev E* 69:036103. Distance-dependent connection probability during growth produces small-world, clustered spatial networks. **[PK]** https://doi.org/10.1103/PhysRevE.69.036103
- Moretti & Muñoz (2013), "Griffiths phases and the stretching of criticality in brain networks", *Nat Commun* 4:2521. Hierarchical-modular, finite-dimensional topologies stretch the critical point into an extended Griffiths phase with slow dynamics and generic power laws. **[PK]** https://doi.org/10.1038/ncomms3521
- Wang & Zhou (2012), "Hierarchical modular structure enhances the robustness of self-organized criticality in neural networks", *New J Phys* 14:023005. **[PK]** https://doi.org/10.1088/1367-2630/14/2/023005
- Zhou, Zemanová, Zamora, Hilgetag & Kurths (2006), "Hierarchical organization unveiled by functional connectivity in complex brain networks", *Phys Rev Lett* 97:238103. **[PK]** https://doi.org/10.1103/PhysRevLett.97.238103
- Ódor, Dickman & Ódor (2015), "Griffiths phases and localization in hierarchical modular networks", *Sci Rep* 5:14451. **[PK]** https://doi.org/10.1038/srep14451
- Gastner & Ódor (2016), "The topology of large Open Connectome networks for the human brain", *Sci Rep* 6:27249. Reports that connectomes have finite topological dimension, which matters for Griffiths phases. **[PK]** https://doi.org/10.1038/srep27249
- Zúñiga-Galindo (2026 preprint), "Hierarchical connectivity sets the upper critical dimension of neural avalanches", arXiv:2610.36404, draft dated 26 Oct 2026. **[V]** https://arxiv.org/html/2610.36404
  - Argues that for short-range coupling in a 2D or 3D Euclidean sheet, directed percolation (DP) is not mean-field (upper critical dimension 4). So τ_S ≈ 3/2 requires hierarchical or ultrametric coupling decaying with exponent σ ≤ N/2.
  - Simulations on a binary tree: τ_S falls from about 1.5 to about 1.2 as σ increases.
  - Single-author, unreviewed preprint. Treat as **uncertain**.

**Morphology → connectome**
- Udvary, Harth, Macke, Hege, de Kock, Sakmann & Oberlaender (2022), "The impact of neuronal structure on cortical network architecture", *Cell Rep* 39:110677. Dense reconstructions of rat barrel cortex morphologies show that axo-dendritic overlap alone predicts many non-random wiring features (e.g. motif overrepresentation). **[PK]** https://doi.org/10.1016/j.celrep.2022.110677

**Fractal-avalanche overlap**
- In the avalanche literature, "fractal" usually refers to the *spatial footprint of avalanches* (fractal dimension of spread between 1 and 2), not neuron morphology. Example: Friedman et al. / Plenz-group analyses, "Scaling properties of neuronal avalanches are consistent with critical dynamics", arXiv:0912.5369. **[V]** https://arxiv.org/abs/0912.5369
- A broad extended search for "dendrite fractal dimension + avalanches + morphology simulation" returned no paper combining them. **[V]** (negative search result)

### Inferences
- Theory suggests morphology affects criticality and dynamic range mainly through:
  1. **Connection range relative to system size**: short-range 2D coupling gives 2D-DP exponents (τ_S ≈ 1.27 for 2D DP, **[PK]**) rather than mean-field 3/2.
  2. **Heterogeneity of in- and out-degree and of weights**, which sets the largest eigenvalue and the critical point.
  3. **Clustering and modularity**, which produce Griffiths-phase stretching.
  4. **Within-arbor excitability**, if dendrites are active (Gollo).

  Dendritic D affects these only indirectly, through how much contact density lands at long versus short distance from the soma.
- An informative experiment would vary D while (a) holding total cable length fixed and (b) holding arbor span fixed, then measure: distance dependence of connection probability p(d), degree variance, clustering and modularity, largest eigenvalue, and avalanche exponents across system sizes. If p(d) and the degree distribution do not change, dynamics should not change either.

### Gaps
- No direct empirical or modelling study links neuron D to avalanche exponents, dynamic range or synchrony at network level. This is a genuine gap in the literature found this session.
- No study was found testing whether *axonal* fractal or space-filling properties alter culture-level avalanches, despite culture models making the axon the main driver of connectivity.

---

## Q6. Why might varying D (1.3, 1.8, 2.5) and base inter-branch length barely change criticality or graph statistics in this simulation? (synthesis question from the project context)

### Takeaway
The literature predicts this null result for several reasons:
- In a 2D embedding, D > 2 is unachievable and box-counting D saturates near the space-filling limit.
- A per-order length ratio r = b^(−1/D) over a few branch orders changes total length and span only modestly. The effective scaling range is about one decade.
- Contact connectivity depends on density overlap and arbor span (Chklovskii; Peters' rule), not on D.
- If the dynamics are tuned or homeostatically normalised to criticality (as in Tetzlaff and Kossio), geometry is factored out.

### Cited Findings
- Real dendrite branches have D ≈ 1.03. Arbor D (about 1.42–1.46) arises mostly from weaving or tortuosity and forking angles, and is about 4× more sensitive to dendrite tortuosity than one might expect. Rowland et al. (2023). **[V]** https://doi.org/10.3389/fnetp.2023.1072815
- Neuronal fractal scaling holds only over about 1–1.5 decades. Rowland et al. (2023) **[V]** (link above); Smith et al. (1996) **[PK]** https://doi.org/10.1016/S0165-0270(96)00080-5
- Retinal arbors can be space-filling (D → 2 to a cutoff), and non-fractal growth rules mimic "fractal" measurements. Panico & Sterling (1995). **[V]** https://doi.org/10.1002/cne.903610311
- Wiring theory: total length and arbor span, set by target density and the balancing factor, control morphology. Cuntz et al. (2010) **[V-exists]** https://doi.org/10.1371/journal.pcbi.1000877 ; Cuntz et al. (2012) **[PK]** https://doi.org/10.1073/pnas.1200430109
- Growth-to-criticality models reach criticality through homeostatic control of coupling, not geometry. Tetzlaff et al. (2010) **[V]** https://doi.org/10.1371/journal.pcbi.1001013 ; Kossio et al. (2018) **[V]** https://arxiv.org/abs/1811.02861

### Inferences
**Mathematical check on the length rule.** With b mean branches per bifurcation and ratio r = b^(−1/D), each order contributes b^k · r^k · L0 = b^{k(1−1/D)} · L0 of cable.
- With b = 2, the factor b^(1−1/D) is 1.17 at D = 1.3, 1.36 at D = 1.8 and 1.52 at D = 2.5. Over 4–6 orders this does change total length by a factor of about 2–3, so D matters through total length.
- The *span* (Σ r^k · L0) changes far less: from about 1.7 L0 to about 2.3 L0 asymptotically. This is my calculation, not a sourced finding.
- If the connection rule saturates (any contact means an edge), or weights or thresholds are re-tuned per condition, these differences vanish.

**Suggested validation:**
1. Box-count the generated images (FracLac-style, several grid offsets) to check that achieved D actually differs between conditions.
2. Report total length, span, p(d) and the degree distribution per condition.
3. Add tortuosity or branching-angle control (the Taylor finding) and a TREES-style target-density control.
4. Use "D" as a realistic-range constraint (about 1.3–1.7 in 2D) rather than the main independent variable. Use connection range, clustering or modularity (Okujeni; Orlandi) and degree heterogeneity as the levers most likely to change criticality, Griffiths-phase width and dynamic range.

### Gaps
- No literature directly predicts how the specific r = b^(−1/D) rule maps onto measured box-counting D. Simulations would need to calibrate this empirically.
- Whether the simulation includes homeostatic or normalised coupling is unknown to this researcher. It determines whether a null result is expected.
