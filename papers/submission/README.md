# Submission Package (arXiv + Zenodo)

This folder describes what to include for submission to arXiv and Zenodo for all three papers:
1. **coherence_gravity_coupling.tex** - Main coherence-gravity coupling paper
2. **null_results.tex** - Null results and exclusion limits
3. **kappaR_to_BSM/curvature_em_to_bsm.tex** - BSM parameter space mapping

---

## Paper 1: Coherence-Gravity Coupling (Main Paper)

### arXiv Classification
**Primary:** gr-qc (General Relativity and Quantum Cosmology)
**Cross-list:** hep-th, physics.ins-det, quant-ph

### Files to Include
- `../coherence_gravity_coupling.tex`
- Figures from `../../results/figures/`
- `../references.bib` (if used)

### Compilation
```bash
cd papers
pdflatex coherence_gravity_coupling.tex
bibtex coherence_gravity_coupling  # if using BibTeX
pdflatex coherence_gravity_coupling.tex
pdflatex coherence_gravity_coupling.tex
```

---

## Paper 2: Null Results and Exclusion Limits

### arXiv Classification
**Primary:** gr-qc (General Relativity and Quantum Cosmology)

### Cross-List Categories
- **hep-th** (High Energy Physics - Theory): For scalar-tensor theories and modified gravity frameworks
- **physics.ins-det** (Instrumentation and Detectors): For torsion balance experimental design and precision measurement techniques
- **quant-ph** (Quantum Physics): For coherence field modeling (BEC, superconductors)

### Files to Include
- `../null_results.tex` (standalone LaTeX manuscript)
- Figures referenced from `../../results/analysis/` (PNG or PDF as needed)
  - `xi_sweep_*.png` (if referenced)
  - `material_comparison_*.png` (if referenced)
  - `curvature_limits_*.png` (if referenced)
- `../references.bib` (if switching from embedded \bibitem to BibTeX)
- Optional: `../../results/reports/analysis_tables.tex` (if tables are \input separately)

### Compilation Instructions (include in comments or README.txt if uploading ancillary files)
```bash
cd papers
pdflatex null_results.tex
pdflatex null_results.tex  # repeat for references and labels
```

---

## Paper 3: BSM Parameter Space (κ_R → Dark Photon/Axion)

### arXiv Classification
**Primary:** hep-ph (High Energy Physics - Phenomenology)

### Cross-List Categories
- **gr-qc** (General Relativity and Quantum Cosmology): For curvature–EM coupling framework
- **hep-ex** (High Energy Physics - Experiment): For experimental bounds and detection prospects
- **astro-ph.CO** (Cosmology and Nongalactic Astrophysics): For dark photon/axion dark matter connections

### Files to Include
- `../kappaR_to_BSM/curvature_em_to_bsm.tex` (main manuscript)
- `../kappaR_to_BSM/table_epsilon.tex` (dark photon table)
- `../kappaR_to_BSM/table_axion.tex` (axion coupling table)
- `../kappaR_to_BSM/figures/epsilon_vs_curvature.pdf`
- `../kappaR_to_BSM/figures/axion_vs_curvature.pdf`
- `../kappaR_to_BSM/figures/curvature_amplification.pdf`
- `../references.bib` (if used)

### Compilation Instructions
```bash
cd papers/kappaR_to_BSM
pdflatex curvature_em_to_bsm.tex
bibtex curvature_em_to_bsm  # if using BibTeX
pdflatex curvature_em_to_bsm.tex
pdflatex curvature_em_to_bsm.tex
```

**Note:** Generate figures first using:
```bash
cd ../..  # back to repository root
python scripts/generate_bsm_plots.py
```

### Suggested arXiv Abstract
Use manuscript abstract emphasizing:
- Mapping κ_R curvature–EM coupling to BSM parameter space (dark photon mixing ε, axion coupling g_aγγ)
- Curvature amplification effect across environments (lab → magnetar)
- Experimental detection prospects via existing BSM searches (APEX, CAST, etc.)
- Laboratory constraints: κ_R < 5×10^17 m^2 (B=10T, R_earth=10^-26 m^-2)

---

## Common: Compilation Instructions

### Suggested arXiv Abstract (if different from manuscript)
Use the manuscript abstract as-is, or condense to ~150 words emphasizing:
- Null result strategy and parameter space coverage
- Exclusion limits on $\kappa_R$ with representative values
- Open-source framework and data availability

---

## Common Zenodo Archive (Data DOI)
Archive the following for reproducibility:
- `results/analysis/*.json` (raw sweep outputs)
- `results/analysis/*.{png,pdf}` (plots)
- `results/reports/analysis_report.md` and `analysis_tables.tex`
- `results/reports/csv/*.csv`

### Recommended Zenodo Metadata

**Title:**  
"Null Results and Exclusion Limits for Coherence–Gravity and Curvature Couplings: Data and Analysis Code"

**Creators:**  
Coherence-Gravity Research Group (or individual contributors as appropriate)

**Description (LaTeX-compatible block):**
```latex
We present numerical null results from systematic parameter sweeps investigating:
(i) non-minimal coupling between a coherence field and spacetime curvature ($\xi R\Phi^2$), and
(ii) novel curvature–electromagnetism interactions ($\kappa_R R F_{\mu\nu}F^{\mu\nu}$).

Configurations tested include coupling strengths $\xi \in \{50, 100\}$, coherent materials 
(Rb-87 BEC, Nb superconducting cavity, YBCO cuprate with $\Phi_0$ spanning $3.65\times10^6$ 
to $6.67\times10^8$ m$^{-1}$), and electromagnetic field strengths $B \in [0.5, 10]$ T.

All measurements yield null signals consistent with numerical baselines ($|\Delta\tau| \approx 5\times10^{-13}$ N·m),
enabling derivation of exclusion limits on the curvature–EM coupling parameter $\kappa_R$.
For $B=10$ T, terrestrial Ricci curvature $R=10^{-26}$ m$^{-2}$, and experimental precision $\delta=10^{-6}$,
we constrain $\kappa_R < 5\times10^{17}$ m$^2$.

This archive includes:
- Raw sweep outputs (JSON timestamped files)
- Analysis plots (PNG/PDF)
- Consolidated tables (CSV, Markdown, LaTeX)
- Complete numerical framework (Python + pytest suite)
- Manuscript source (LaTeX)

Code repository: https://github.com/DawsonInstitute/coherence-gravity-coupling  
Manuscript preprint: arXiv:XXXX.XXXXX (to be assigned)
```

**Keywords / Subjects:**
- Non-minimal coupling
- Curvature–electromagnetism coupling
- Exclusion limits
- Tabletop precision gravimetry
- Null results
- Modified gravity
- Coherence field
- Torsion balance experiments
- Beyond-GR physics
- Numerical gravitation
- Open-source scientific software

**License:**  
MIT License (code), CC BY 4.0 (documentation/data)

**Related Identifiers:**
- GitHub Repository: https://github.com/DawsonInstitute/coherence-gravity-coupling
- Manuscript (arXiv): arXiv:XXXX.XXXXX (add after submission)
- Companion Paper DOI: (if applicable, link validation study)

**Version:**  
v1.0.0 (November 2025 release)

**Upload Type:**  
Dataset + Software

**Notes:**  
After DOI assignment, add the Zenodo DOI badge to the repository README and cite in manuscript Data Availability section.
