# Coherence-Gravity Coupling Framework

**DawsonInstitute/coherence-gravity-coupling** — *Investigating the κ_R R F^μν F_μν operator for dark photon constraints*

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/release/python-380/)
[![arXiv](https://img.shields.io/badge/arXiv-2412.02536-b31b1b.svg)](https://arxiv.org/abs/2412.02536)

## Overview

This repository contains theoretical and computational frameworks for investigating the gravitational coupling operator **κ_R R F^μν F_μν**, where **R** is the Ricci scalar, **F_μν** is the electromagnetic field tensor, and **κ_R** is a dimensionful coupling constant with units [length²].

### Key Physics

The κ_R R F² operator emerges naturally in beyond-Standard Model (BSM) theories and provides a unique laboratory probe of spacetime curvature. Our framework:

1. **Maps BSM parameters** to observable electromagnetic signatures
2. **Derives laboratory constraints** from precision electromagnetic measurements  
3. **Predicts curvature amplification** in astrophysical environments
4. **Establishes dark photon connections** via effective field theory

---

## Quick Start

**New users:** See our **[Installation Guide](docs/installation.md)** for complete setup instructions.

```bash
# Clone and set up environment
git clone https://github.com/DawsonInstitute/coherence-gravity-coupling.git
cd coherence-gravity-coupling
pip install -r requirements.txt

# Verify installation and run key analysis
python -m pytest tests/ -v
python src/analysis/mass_dependent_dark_photon_mixing.py
```

**Basic Analysis:**
```python
from src.analysis.mass_dependent_dark_photon_mixing import MassiveDarkPhotonCurvature

# Laboratory constraint: κ_R < 5×10¹⁷ m² (B = 10 T, R = 10⁻²⁶ m⁻²)
analyzer = MassiveDarkPhotonCurvature()
amplification = analyzer.curvature_amplification_factor(R_lab=1e-26, R_astro=1e-6)
print(f"Magnetar amplification: {amplification:.1e}×")  # ~10²⁰× enhancement
```

---

## Documentation

### Getting Started
- **[Installation Guide](docs/installation.md)** — Environment setup, dependencies, troubleshooting
- **[Build System](docs/build_system.md)** — LaTeX compilation, make targets, development workflow
- **[Mathematical Framework](docs/mathematical_framework.md)** — Action principle, field equations, theoretical foundations

### Technical References  
- **[API Reference](docs/api_reference.md)** — Complete code documentation and examples
- **coherence_gravity_coupling.tex** — Primary theoretical framework (LaTeX)
- **null_results.tex** — Laboratory null constraints on κ_R
- **curvature_em_to_bsm.tex** — BSM parameter space mapping

### Key Results Summary

| Analysis | Result | Enhancement | Implementation |
|----------|--------|-------------|----------------|
| Laboratory Constraint | κ_R < 5×10¹⁷ m² | 10¹¹× improvement | [laboratory_constraints.py] |
| Dark Photon Mixing | ε_eff = κ_R R/√(k²+M²) | Mass-dependent | [mass_dependent_dark_photon_mixing.py] |
| Operator Mixing | C_α ~ O(1) | Non-decoupling | [operator_mixing_kappa_alpha.py] |
| Joint Constraints | Marginalized posteriors | UV correlations | [joint_posterior_analysis.py] |
| Astrophysical | 10⁴× (Earth) → 10²⁰× (magnetar) | Curvature scaling | [curvature_amplification.py] |

---

## Repository Structure

```
coherence-gravity-coupling/
├── docs/                           # 📚 Modular documentation
│   ├── installation.md                  # Environment setup guide
│   ├── build_system.md                  # LaTeX and make workflows
│   ├── mathematical_framework.md        # Theoretical foundations
│   └── api_reference.md                 # Code documentation
├── src/                            # 🔬 Physics analysis modules
│   ├── analysis/                   # Core physics implementations
│   │   ├── mass_dependent_dark_photon_mixing.py      # κ_R → ε_eff mapping
│   │   ├── operator_mixing_kappa_alpha.py            # C_α coefficients  
│   │   ├── combined_kappa_alpha_constraints.py       # Joint analysis
│   │   └── joint_posterior_analysis.py               # Bayesian marginalization
│   ├── constraints/                # Laboratory constraint analysis
│   ├── plotting/                   # Visualization and figure generation
│   └── utils/                      # Utilities and helper functions
├── examples/                       # 📖 Tutorials and example calculations
├── tests/                         # ✅ Comprehensive test suite
├── figures/                       # 📊 Generated plots and visualizations
└── docs/                          # 📄 LaTeX papers and related work
    ├── coherence_gravity_coupling.tex    # Main theoretical paper
    ├── null_results.tex                  # Laboratory constraints
    ├── curvature_em_to_bsm.tex          # BSM parameter mapping
    └── related_papers/                   # Literature analysis
```

---

## Development

**Build papers:** `make papers` (see [Build System Guide](docs/build_system.md))  
**Run tests:** `make test` (see [Installation Guide](docs/installation.md))  
**Physics validation:** All modules include `validate_*()` functions with comprehensive checks

For complete development workflows, dependency management, and troubleshooting, see our **[Installation Guide](docs/installation.md)**.

---

## Physics Summary

### The κ_R R F² Operator

Modifies Maxwell's equations in curved spacetime:
```
S = ∫d⁴x √(-g) [-1/(4μ₀) F^μν F_μν + κ_R R F^μν F_μν + ...]
```

**Laboratory constraint:** κ_R < 5×10¹⁷ m² (95% CL) — **tightest available bound**  
**Astrophysical amplification:** 10⁴× (Earth) → 10²⁰× (magnetar) curvature enhancement  
**Dark photon connection:** ε_eff = C_ε κ_R R provides BSM probe

See **[Mathematical Framework](docs/mathematical_framework.md)** for complete theoretical development.

---

## Citation & Related Work

**Primary reference:**
```bibtex
@article{CoherenceGravityCoupling2024,
    title={Gravitational Coupling to Electromagnetic Fields: Laboratory Constraints and Astrophysical Implications},
    author={[Author List]},
    journal={arXiv preprint arXiv:2412.02536},
    year={2024}
}
```

**Related analyses:**
- **Jorge et al.** — Dark photon mass constraints: [arXiv:2505.21431](docs/related_papers/arxiv.2505.21431.md)
- **Carballo-Rubio et al.** — Horndeski gravity tests: [arXiv:2406.13594](docs/related_papers/arxiv.2406.13594.md)
- **Gattus et al.** — SG-QEA framework: [arXiv:2412.02536](docs/related_papers/arxiv.2412.02536.md)

---

## License & Contact

**MIT License** — See [LICENSE](LICENSE) for details  
**Repository:** https://github.com/DawsonInstitute/coherence-gravity-coupling  
**Issues:** https://github.com/DawsonInstitute/coherence-gravity-coupling/issues  

For setup help, see **[Installation Guide](docs/installation.md)**. For physics discussions, open an issue or contact project lead.