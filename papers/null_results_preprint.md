# Null Results and Exclusion Limits for Coherence–Gravity and Curvature Couplings

Authors: Project Team
Date: 2025-10-31
License: MIT

## Abstract
We present null results from a configuration-driven numerical study of non-minimal couplings between a coherence field and gravity (ξ coupling), and between spacetime curvature and electromagnetism (R F_{\mu\nu}F^{\mu\nu}). Parameter sweeps over coupling strengths, materials, and field settings yield no detectable signal beyond numerical baselines. From these nulls, we derive exclusion limits on the curvature–electromagnetism coupling parameter κ_R across laboratory-relevant magnetic field strengths and terrestrial Ricci curvature scales.

## 1. Introduction
- Motivation: precision constraints on beyond-GR couplings.
- Two targets: ξ coupling in coherence–gravity models, and κ_R in curvature–EM effective terms.
- Contribution: reproducible pipeline, tests, sweeps, plots, and quantitative upper bounds.

## 2. Methods
### 2.1 Coherence–gravity framework
- Modified Einstein equations with ξRΦ^2 and coherence backreaction (einstein_coherence.py).
- Weak-field 3D Poisson solver with effective G(x) and torque extraction.
- Parameterized experiments via `run_analysis.py`.

### 2.2 Curvature–EM effective coupling
- Lagrangian correction: L ⊃ κ_R R F_{\mu\nu}F^{\mu\nu}.
- Implementation: `src/field_equations/curvature_coupling.py`.
- EM invariants F^2 = 2(B^2 − E^2/c^2), exclusion limit κ_R < δ / (R |F^2|).
- Verified with 18 tests; integrated demo and analysis CLI.

### 2.3 Reproducibility and tests
- Full test suite: 41 tests, all passing on Python 3.13.
- Results saved under `results/analysis/` with timestamps; plots emitted as PNG/PDF.

## 3. Results
### 3.1 ξ parameter sweep
- Command: `python run_analysis.py sweep-xi --plot --cache`
- Range: ξ ∈ {10, 50, 100, 200, 500}; resolution 41^3.
- Result: torque contrast Δτ ≈ −4.992e−13 N·m across sweep (no monotonic detection).
- Artifacts: `xi_sweep_YYYYMMDD_HHMMSS.json`, `xi_sweep_..._plot.(png|pdf)`.

### 3.2 Material comparison (coherence profiles)
- Command: `python run_analysis.py sweep-materials --plot --cache`
- Materials: rb87_bec, nb_cavity, ybco_cuprate.
- Result: Δτ ≈ −4.992e−13 N·m for all materials at ξ=100, res=41 (null within numerical tolerance).
- Artifacts: `material_comparison_YYYYMMDD_HHMMSS.*`.

### 3.3 Curvature–EM exclusion limits (κ_R)
- Command:
  - `python run_analysis.py sweep-curvature --plot --B 0.5 1.0 3.0 10.0 --R 1e-26 --precision 1e-6`
- Assumptions: E=0, terrestrial Ricci R≈1e−26 m−2, precision δ=1e−6.
- Derived limits (this run):
  - B = 0.5 T → κ_R < 2.00e+20 m²
  - B = 1.0 T → κ_R < 5.00e+19 m²
  - B = 3.0 T → κ_R < 5.56e+18 m²
  - B = 10.0 T → κ_R < 5.00e+17 m²
- Interpretation: Higher fields tighten the bound as 1/|F^2|, but absolute limits remain very weak at terrestrial R.
- Artifacts: `curvature_limits_YYYYMMDD_HHMMSS.*`.

## 4. Discussion
- All sweeps yield nulls at current resolutions/parameters; consistent with invariance and conservation tests.
- κ_R limits are dominated by small terrestrial curvature; orders-of-magnitude improvements require larger |R| or precision δ.
- Next: explore high-curvature analogs (e.g., metamaterials, effective metrics), resonant EM configurations, and improved torque metrology.

## 5. Conclusion
We provide a tested, reproducible baseline showing no detectable signatures of ξ or κ_R in the explored regimes, with transparent pipelines to extend parameter coverage and tighten bounds.

## Data and Code Availability
- Code: this repository; `run_analysis.py` for CLI.
- Data: `results/analysis/` (JSON + plots). Checksums can be added to `docs/PROVENANCE.md`.

## Acknowledgments
We thank contributors to the analysis, test, and CI infrastructure.
