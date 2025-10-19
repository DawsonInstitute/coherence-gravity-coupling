# Reproducibility Guide

This document provides explicit commands to reproduce all results reported in the manuscript `coherence_gravity_coupling.tex`.

## Prerequisites

- Python 3.11 (pinned)
- NumPy 1.26.4, SciPy 1.14.1, Matplotlib 3.9.4, PyAMG ≥ 4.2.0
- 32 GB RAM recommended for 81³+ resolutions
- 4+ CPU cores for parallelization (--jobs 4)

## Environment Setup

```bash
git clone https://github.com/arcticoder/coherence-gravity-coupling.git
cd coherence-gravity-coupling

# Option A: conda (recommended)
conda env create -f environment.yml
conda activate cohgrav

# Option B: venv + pip
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

pytest -q  # Verify installation (23 tests, ~94s)
```

## Key Manuscript Results

### 1. Validated 61³ Signal (Section 3.1)

**Manuscript Claim**: τ_coh = 1.4 ± 0.2 × 10⁻¹² N·m at position (0.0012, 0.0182, 0.0659) m

**Data File**: `results/validation/validate_61_Rb87_ξ100_powell.json`

**Reproduction Command**:
```bash
python optimize_geometry.py --xi 100 --Phi0 3.65e6 --resolution 61 --method Powell --material Rb87
```

**Expected Output**:
- Optimal position: r_coh ≈ (0.001, 0.018, 0.066) m (±0.002 m tolerance)
- Torque: τ_coh ≈ 1.4 × 10⁻¹² N·m (±15% variation due to optimizer stochasticity)
- Runtime: ~15-30 minutes depending on CPU

### 2. Convergence Study (Section 3.2)

**Manuscript Claim**: Richardson extrapolation yields Δτ ≈ 2.6 × 10⁻¹² N·m with p ≈ 2.1

**Data File**: `CONVERGENCE_ANALYSIS.md`

**Reproduction Commands**:
```bash
# 61³ resolution
python run_geometric_cavendish.py --xi 100 --Phi0 3.65e6 --resolution 61 --use-volume-average

# 81³ resolution (requires ~30 min)
python run_geometric_cavendish.py --xi 100 --Phi0 3.65e6 --resolution 81 --use-volume-average

# 101³ resolution (requires ~2 hours, 32+ GB RAM)
python run_geometric_cavendish.py --xi 100 --Phi0 3.65e6 --resolution 101 --use-volume-average

# Richardson extrapolation (manual analysis)
# Use data from above runs to fit Δτ(N) = Δτ_∞ + A/N^p
# Expected: p ≈ 2.1, Δτ_∞ ≈ 2.6 × 10⁻¹² N·m
```

**Validation**: Results should show +13% increase from 61³→81³, +17% from 81³→101³

### 3. Production Grid Study (Section 3.3)

**Manuscript Claim**: 5×5×5 grid search finds optimal at (0, 0, -0.05) m for all materials, |Δτ| ≈ 1.099 × 10⁻¹² N·m

**Data File**: `results/production_study/production_study_20251018_204142.json` (included) or regenerate via command below

**Reproduction Command**:
```bash
python production_study.py --materials all --resolution 61 --grid-size 5 --jobs 4 --quick
```

**Expected Output**:
- Runtime: ~30-40 minutes (4 workers, with caching)
- YBCO: 125 points in ~990 s (7.9 s/point)
- Rb87: 125 points in ~930 s (7.4 s/point)
- Nb: 125 points in ~3 s (fully cached)
- Optimal position: (0, 0, -0.05) m for all materials
- Signal magnitude: |Δτ| ≈ 1.1 × 10⁻¹² N·m (within 10% of manuscript value)

### 4. Artifact Correction Analysis (Section 4.2)

**Manuscript Claim**: 41³ DE optimization showed spurious 523× enhancement; validated 61³ shows modest ~1.4× improvement

**Data File**: `VALIDATION_REPORT.md`

**Reproduction Commands**:
```bash
# 41³ grid search (baseline)
python optimize_geometry.py --xi 100 --Phi0 6.67e8 --resolution 41 --grid-search --grid-size 5

# 41³ DE refinement (artifact-prone)
python optimize_geometry.py --xi 100 --Phi0 6.67e8 --resolution 41 --method DE

# 61³ Powell validation (corrected)
python optimize_geometry.py --xi 100 --Phi0 6.67e8 --resolution 61 --method Powell
```

**Expected Behavior**:
- 41³ DE may show anomalously large enhancements (100-500×) due to grid aliasing
- 61³ Powell should show modest enhancement (~1-2×) over grid optimum
- This demonstrates resolution dependence and validates artifact correction

## Data File Manifest

All manuscript results are traceable to the following data files:

### Validation Data
- `results/validation/validate_61_Rb87_ξ100_powell.json` - Primary signal validation (Fig. 1 reference)
- `results/validation/validate_61_Nb_ξ100_powell.json` - Material universality confirmation

### Convergence Data
- `CONVERGENCE_ANALYSIS.md` - Richardson extrapolation analysis
- `results/convergence_test_20251018_193745.json` - 61³ convergence run
- `results/convergence_test_101_20251018_193916.json` - 101³ convergence run

### Production Data
- `results/production_study/production_study_20251018_204142.json` - Complete 5³ grid × 3 materials at 61³
- `LATEST_PRODUCTION_SUMMARY.md` - Human-readable summary with integration time estimates

### Figures
- `papers/figures/convergence_analysis.pdf|png` - Figure 1 (generated)
- `papers/figures/material_comparison.pdf|png` - Figure 2 (copied from results)
- `papers/figures/landscape_YBCO_z_slice.pdf|png` - Figure 3 (copied from results)

## Numerical Consistency Checks

### Conservation Laws
```bash
pytest tests/test_coherence_invariance.py::test_conservation_laws -v
```
Expected: Energy-momentum conservation verified to <1% relative error

### Volume Averaging Convergence
```bash
pytest tests/test_volume_average.py::test_volume_avg_convergence -v
```
Expected: Volume-averaged force converges to point-sample as radius → 0

### Material Universality
```bash
python run_analysis.py sweep-materials --xi 100 --resolution 61 --cache --plot
```
Expected: Identical torque profiles for YBCO/Rb87/Nb at fixed ξ (within 5% numerical noise)

## Performance Benchmarks

Expected runtimes on Intel i7-10700K (8 cores, 32 GB RAM):

| Resolution | Single Solve | Grid Search (125 pts) | Convergence Study (3 res) |
|------------|-------------|-----------------------|---------------------------|
| 41³        | ~3-5 s      | ~7 min (4 workers)    | ~15 min                   |
| 61³        | ~5-8 s      | ~32 min (4 workers)   | ~1 hour                   |
| 81³        | ~20-30 s    | ~2.5 hours            | ~6 hours                  |
| 101³       | ~1-2 min    | ~12 hours             | ~36 hours                 |

**Note**: Caching provides ~250× speedup on repeated geometries. First run at each resolution will be slower.

## Independent Verification

To verify manuscript claims independently:

1. **Clone repository** (as above)
2. **Run test suite**: `pytest tests/ -v` (23/23 should pass)
3. **Execute validation runs** (commands in §1-§4 above)
4. **Compare outputs** to data files in `results/` and claims in manuscript
5. **Report discrepancies** via GitHub issues

### Expected Variations

- **Position coordinates**: ±0.005 m tolerance (optimizer convergence stochasticity)
- **Torque magnitude**: ±15% variation (grid resolution, numerical precision)
- **Runtime**: ±30% variation (CPU architecture, thermal throttling)
- **Convergence order p**: ±0.2 variation (Richardson extrapolation fitting)

All variations should be within expected numerical precision. If discrepancies exceed these bounds, please open a GitHub issue with:
- Command executed
- Observed output
- Expected output (from manuscript/data files)
- System specifications (CPU, RAM, Python version)

## Contact

For reproducibility issues or questions:
- GitHub Issues: https://github.com/arcticoder/coherence-gravity-coupling/issues
- Email: Contact details are provided in the PDF uploaded to Zenodo (author list and corresponding email). The repository includes `papers/author_config.tex` used during LaTeX compilation.

---

*Last updated: October 18, 2025*  
*Manuscript version: coherence_gravity_coupling.tex (LaTeX source)*  
*Data snapshot: production_study_20251018_204142.json (61³ complete grid)*
