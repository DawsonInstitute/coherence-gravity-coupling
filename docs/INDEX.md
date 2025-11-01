# Coherence-Gravity Coupling Documentation Index

**Last Updated**: October 17, 2025  
**Version**: Phase D Complete  
**Status**: ✅ Analysis Complete - Ready for Experimental Collaboration

---

## Quick Navigation

### 📖 For New Users
- **Start here**: [`README.md`](README.md) - Complete theoretical background and implementation overview
- **Quick start**: [`QUICKREF.md`](QUICKREF.md) - Commands, API examples, and troubleshooting
- **Try it**: Run `pytest tests/ -v` to verify installation (20 tests, ~94s)

### 🔬 For Researchers
- **Theory**: [`README.md`](README.md) sections on:
  - Modified field equations with non-minimal coupling $\xi R \Phi^2$
  - Weak-field limit and Poisson equation
  - Observational constraints (Solar System, binary pulsars, cosmology)
  - Dimensional analysis and physical interpretation
- **Results**: [`README.md#Key-Result-Summary`](README.md#-key-result-summary)
  - Corrected torque scales: τ_N ~ 2×10⁻¹³ N·m, Δτ ~ 1.6×10⁻¹² N·m
  - Feasibility: 9/18 configs achievable with cryo operation
  - Best case: 0.7 hr integration time (YBCO, ξ=100, cryo_moderate)

### 💻 For Developers
- **Implementation**: [`PROGRESS_SUMMARY.md`](PROGRESS_SUMMARY.md) - Detailed technical summary
- **API Reference**: [`QUICKREF.md#Python-API`](QUICKREF.md#python-api)
- **Test Suite**: `tests/` directory (6 files, 20 tests)
- **Example Scripts**: `examples/` directory
  - `refined_feasibility.py` - Noise analysis and optimization
  - `geometric_cavendish.py` - 3D Poisson solver and torque computation

### 📊 For Experimentalists
- **Feasibility**: [`README.md#Key-Result-Summary`](README.md#-key-result-summary)
- **Experimental Setup**: [`README.md#Try-It-Yourself`](README.md#try-it-yourself)
- **Noise Budget**: `examples/refined_feasibility.py --sweep`
- **Optimization**: `examples/refined_feasibility.py --optimize`
- **Requirements**:
  - Temperature: 4K (liquid He) or 77K (liquid N₂)
  - Seismic isolation: 10-100× suppression
  - Torsion constant: κ ~ 10⁻⁸ N·m/rad
  - Readout: <1 nrad/√Hz

### 📝 For Project Management
- **Session Summary**: [`SESSION_SUMMARY_Oct17_2025.md`](SESSION_SUMMARY_Oct17_2025.md)
- **Progress Tracking**: [`PROGRESS_SUMMARY.md`](PROGRESS_SUMMARY.md)
- **Todo List**: See inline task list (5 items remaining)

---

## File Structure

```
coherence-gravity-coupling/
├── README.md                          # Main documentation (33 KB)
├── QUICKREF.md                        # Quick reference guide (6 KB)
├── PROGRESS_SUMMARY.md                # Technical implementation summary (11 KB)
├── SESSION_SUMMARY_Oct17_2025.md     # Today's session summary (8 KB)
├── INDEX.md                           # This file
│
├── src/                               # Core library
│   ├── solvers/
│   │   └── poisson_3d.py             # 3D Poisson solver with G_eff(Φ)
│   └── analysis/
│       └── phi_calibration.py        # Φ₀ calibration for BEC/SC systems
│
├── examples/                          # Example scripts
│   ├── refined_feasibility.py        # ⭐ Noise analysis + optimization
│   ├── geometric_cavendish.py        # 3D geometric Cavendish simulation
│   └── figures/                      # Generated figures
│       ├── feasibility_integration_times.png
│       ├── noise_profile_sweep.png
│       └── optimized_vs_baseline.png
│
├── tests/                             # Test suite (20 tests)
│   ├── test_coherence_invariance.py  # Coherence-specific tests (5)
│   ├── test_volume_average.py        # Volume averaging tests (3)
│   ├── test_newtonian_torque_scale.py # Scale validation
│   ├── test_conservation.py          # Energy-momentum conservation
│   ├── test_field_equations.py       # Modified Einstein equations
│   └── test_interface_matching.py    # Boundary condition tests
│
└── results/                           # Simulation results
    └── geometric_cavendish_sweep.json # 18 baseline configurations
```

---

## Key Equations

### Modified Einstein Equation
$$G_{\mu\nu} + \xi \left[2(\nabla_\mu\nabla_\nu - g_{\mu\nu}\square)\Phi^2 + 2\Phi^2 G_{\mu\nu} - 4\nabla_\mu\Phi\nabla_\nu\Phi + 2g_{\mu\nu}(\nabla\Phi)^2\right] = 8\pi G T_{\mu\nu}$$

### Effective Coupling
$$G_{\text{eff}}(\Phi) = \frac{G}{1 + 8\pi G \xi \Phi^2}$$

### Weak-Field Poisson Equation (Corrected)
$$\nabla \cdot \left(\frac{G_{\text{eff}}(x)}{G} \nabla \phi\right) = 4\pi G \rho(x)$$

---

## Usage Examples

### Command Line
```bash
# Basic feasibility
python examples/refined_feasibility.py --profile cryo_moderate

# Compare all noise profiles
python examples/refined_feasibility.py --sweep

# Optimize geometry
python examples/refined_feasibility.py --profile cryo_moderate --optimize

# Run tests
pytest tests/ -v
```

### Python API
```python
from examples.geometric_cavendish import run_geometric_cavendish

# Run with volume averaging (recommended for convergence)
result = run_geometric_cavendish(
    xi=100.0,
    Phi0=6.67e8,  # YBCO
    coherent_position=(0.0, 0.0, -0.08),
    grid_resolution=61,
    use_volume_average=True,
    verbose=True
)

print(f"Δτ = {result['delta_tau']:.3e} N·m")
print(f"ΔG/G = {result['delta_G_over_G']:.2f}")
```

---

## Coherence Systems Catalog

| System | Φ₀ (m⁻¹) | Temperature | ΔG/G Sign | Notes |
|--------|----------|-------------|-----------|-------|
| **Rb87 BEC** | 3.65×10⁶ | ~100 nK | Negative (offset) | Lowest Φ₀; challenging signal |
| **Nb Cavity** | 2.63×10⁷ | ~4 K | Negative (offset) | Mid-range Φ₀ |
| **YBCO Cuprate** | 6.67×10⁸ | ~90 K | Positive (offset) | **Best signal**; highest Φ₀ |

---

## Noise Profiles

| Profile | T (K) | Seismic | Tilt | Readout | Mass | Feasible Configs |
|---------|-------|---------|------|---------|------|------------------|
| room_temp_baseline | 300 | 1× | 1× | 1× | 1× | 0/18 |
| **cryo_moderate** | 4 | 10× | 10× | 10× | 1× | **9/18** ⭐ |
| cryo_advanced | 4 | 30× | 10× | 10× | 1× | 9/18 |
| optimized | 4 | 100× | 10× | 10× | 10× | 9/18 |

**Recommendation**: Use `cryo_moderate` as baseline; realistic and achievable

---

## Current Limitations

1. **Grid convergence**: 41³→61³ shows ~220% Δτ change; use ≥81³ with volume averaging
2. **Experimental**: Cryogenic operation required (0/18 room-temp configs feasible)
3. **Integration times**: 0.7-24 hours (not <1 second as initially estimated)
4. **Coherence achievability**: Φ₀ ~ 10⁸ m⁻¹ extrapolated from condensed matter; needs validation
5. **Decoherence**: Environmental coupling not modeled; could reduce effective Φ₀

See [`README.md#Limitations`](README.md#limitations-and-numerical-considerations) for full discussion.

---

## Roadmap

### Completed ✅
1. Poisson solver normalization fix
2. Noise parameterization (4 profiles)
3. Geometry optimization sweeps
4. Trilinear interpolation
5. Volume-averaged force
6. CLI optimization integration
7. Comprehensive test suite (20 tests)
8. Documentation (README + PROGRESS + QUICKREF)

### Next Steps
1. **Refactor geometry parameters** - Remove scaling approximations
2. **Automate convergence studies** - Add `--convergence` CLI flag
3. **Solver performance** - Better preconditioning for ≥81³
4. **Domain/BC study** - Systematic sensitivity analysis
5. **Result caching** - Avoid redundant solves

See [`PROGRESS_SUMMARY.md#Pending-Work`](PROGRESS_SUMMARY.md#pending-work-prioritized) for details.

---

## Citation

If you use this framework, please cite:

```bibtex
@software{coherence_gravity_coupling_2025,
  title = {Coherence-Gravity Coupling Framework},
  author = {GitHub Copilot (Claude Sonnet 4.5)},
  year = {2025},
  month = {10},
  note = {Loop Quantum Gravity modification to Cavendish torsion balance},
  url = {https://github.com/DawsonInstitute/coherence-gravity-coupling}
}
```

---

## Support

- **Documentation**: Read [`README.md`](README.md), [`QUICKREF.md`](QUICKREF.md)
- **Issues**: Check existing tests and docs first; open GitHub issue with config details
- **Questions**: See "Open Research Questions" in [`README.md#Limitations`](README.md#limitations-and-numerical-considerations)
- **Collaboration**: Contact via GitHub for experimental implementation discussions

---

## License

MIT License - See repository root for details

---

## Changelog

### October 17, 2025 - Phase D Complete
- Added volume-averaged force computation
- Integrated `--optimize` CLI flag
- Created comprehensive test suite (20 tests)
- Documentation refresh: PROGRESS_SUMMARY, QUICKREF, limitations section
- Updated GitHub topics (16 tags)
- Session summary and index created

### [Earlier] - Phase D Implementation
- 3D Poisson solver with G_eff(Φ)
- Geometric Cavendish simulation
- Noise budget analysis
- Φ₀ calibration for BEC/SC systems
- Trilinear interpolation
- Geometry optimization sweeps

---

**End of Index**

**Quick Start**: `cat QUICKREF.md` → `pytest tests/ -v` → `python examples/refined_feasibility.py --help`
