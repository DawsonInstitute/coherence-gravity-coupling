# From Curvature–EM Coupling to BSM Parameter Space

**Framework linking κ_R to dark photon and axion benchmarks**

## Overview

This paper presents the first systematic framework to connect laboratory bounds on the curvature–electromagnetism coupling κ_R to beyond-Standard-Model (BSM) parameter space, specifically dark photon kinetic mixing (ε) and axion-photon coupling (g_aγγ).

## Key Results

### Dark Photon Mapping (Robust)
- **Formula**: ε_eff ≈ C_ε (κ_R · R)
- **Status**: Dimensionally consistent, conservative
- **Physical meaning**: Quantifies Maxwell equation modifications in curved backgrounds

For κ_R ~ 10^-11 m²:
- Lab (flat space): ε_eff ~ 10^-41 (far below current limits)
- Earth surface: ε_eff ~ 10^-37 (approaching future sensitivity)
- Magnetar: ε_eff ~ 10^-17 (constrains matching coefficients)

### Axion Mapping (Model-Dependent)
- **Formula**: g_equiv ≈ C_a (κ_R · R) / Λ
- **Status**: Parametric only; requires CP-odd portal
- **Critical caveat**: κ_R F² is CP-even; direct equivalence to axion (CP-odd) requires additional physics

## Files

- `main.tex` - Main paper draft
- `table_epsilon.tex` - Dark photon results table
- `table_axion.tex` - Axion benchmark table

## Compilation

```bash
cd papers/kappaR_to_BSM
pdflatex main.tex
pdflatex main.tex  # Second pass for references
```

Or with latexmk:
```bash
latexmk -pdf main.tex
```

## Generating Tables

Tables are auto-generated from computed data:

```bash
# From repository root
python scripts/run_bsm_bounds.py          # Generate CSV data
python scripts/generate_bsm_tables.py     # Generate LaTeX tables
```

This updates `table_epsilon.tex` and `table_axion.tex` which are included in the main document.

## Experimental Context

### Dark Photon Limits
- **APEX** (2011): ε < 10^-3 for m_A' ~ 100 MeV
- **BaBar** (2014): ε ~ 10^-4 to 10^-3 for m_A' up to GeV scale

### Axion Limits
- **CAST** (2017): g_aγγ < 10^-10 GeV^-1 for m_a ~ 0.02 eV
- **ADMX** (2021): g_aγγ ~ 10^-15 to 10^-14 GeV^-1 for μeV masses

## Key Innovation

**Curvature amplification**: Laboratory κ_R bounds (seemingly irrelevant for BSM) become potent probes of BSM physics in strong-curvature astrophysical environments. This opens two complementary directions:

1. Use astrophysical observations to constrain κ_R
2. Search for enhanced BSM signatures in curved spacetime

## Code Modules

- `src/analysis/bsm_bounds_from_kappa.py` - Core mapping functions
- `src/analysis/kappa_k3_mapping.py` - Integration with torsion-EM framework
- `scripts/run_bsm_bounds.py` - Generate CSV results
- `scripts/generate_bsm_tables.py` - Create LaTeX tables

## Data Products

- `results/bsm_bounds/epsilon_equiv.csv` - Dark photon ε values
- `results/bsm_bounds/g_axion_equiv.csv` - Axion g_equiv values

## Citation

```bibtex
@article{DawsonBSM2025,
  title={From Curvature--EM Coupling to BSM Parameter Space: A Framework Linking $\kappa_R$ to Dark Photon and Axion Benchmarks},
  author={Dawson Institute Collaboration},
  year={2025},
  note={In preparation}
}
```

## Status

**Draft ready for internal review** - All core calculations complete, experimental context provided, tables generated from data.

## Next Steps

1. Optional: Add comparison plots (ε vs m_A', g vs m_a)
2. Refine curvature environments with specific geometries (Schwarzschild, FLRW)
3. Discuss CP-violating portal mechanisms in more detail
4. Add error analysis for matching coefficients
