# Coherence-Gravity Coupling Quick Reference

## Installation & Setup
```bash
cd /home/echo_/Code/asciimath/coherence-gravity-coupling
# Dependencies: numpy, scipy, matplotlib, pytest
python -m pytest tests/ -v  # Verify installation
```

## CLI Commands

### Basic Feasibility Analysis
```bash
# Room temperature baseline (no isolation)
python examples/refined_feasibility.py

# Cryogenic with moderate isolation
python examples/refined_feasibility.py --profile cryo_moderate

# Advanced cryogenic setup
python examples/refined_feasibility.py --profile cryo_advanced

# Optimized (4K + 100× seismic + 10× mass)
python examples/refined_feasibility.py --profile optimized
```

### Noise Profile Comparison
```bash
# Compare all 4 profiles across 18 configs
python examples/refined_feasibility.py --sweep

# Output: examples/figures/noise_profile_sweep.png
```

### Geometry Optimization
```bash
# Run optimization for representative configs
python examples/refined_feasibility.py --profile cryo_moderate --optimize

# Output: examples/figures/optimized_vs_baseline.png
# Shows integration time improvements from geometry optimization
```

## Python API

### Run Single Configuration
```python
from examples.geometric_cavendish import run_geometric_cavendish

result = run_geometric_cavendish(
    xi=100.0,                           # Coupling strength
    Phi0=6.67e8,                        # YBCO coherence field [m⁻¹]
    coherent_position=(0.0, 0.0, -0.08), # BEC/SC position [m]
    grid_resolution=41,                 # Grid points per dimension
    use_interpolation=True,             # Trilinear interpolation
    use_volume_average=False,           # Volume-averaged force (recommended for convergence)
    verbose=True
)

print(f"Δτ = {result['delta_tau']:.3e} N·m")
print(f"ΔG/G = {result['delta_G_over_G']:.2f}")
```

### Geometry Optimization
```python
from examples.geometric_cavendish import optimize_geometry

opt = optimize_geometry(
    xi=100.0,
    Phi0=6.67e8,
    verbose=True
)

print(f"Optimal position: {opt['final_config']['coherent_position']}")
print(f"Optimal Δτ: {opt['final_config']['delta_tau']:.3e} N·m")
```

### Noise Budget Calculation
```python
from examples.refined_feasibility import (
    NoiseProfile, 
    total_noise_budget,
    integration_time_for_snr
)

# Define custom noise profile
profile = NoiseProfile(
    name='My Setup',
    T=4.0,                        # Kelvin
    seismic_suppression=50.0,     # Isolation factor
    tilt_suppression=10.0,        # Tilt reduction
    readout_improvement=5.0,      # Detector sensitivity
    m_test_factor=2.0             # Test mass scaling
)

noise = total_noise_budget(profile)
print(f"Total noise: {noise['total']:.3e} N·m/√Hz")

# Integration time for signal
signal = 1.6e-12  # N·m (YBCO typical)
t_int_hr = integration_time_for_snr(signal, noise['total'], target_snr=5.0) / 3600
print(f"Integration time: {t_int_hr:.1f} hr")
```

## Noise Profiles

| Profile | T (K) | Seismic | Tilt | Readout | Mass | Use Case |
|---------|-------|---------|------|---------|------|----------|
| `room_temp_baseline` | 300 | 1× | 1× | 1× | 1× | Initial estimate |
| `cryo_moderate` | 4 | 10× | 10× | 10× | 1× | Realistic cryo |
| `cryo_advanced` | 4 | 30× | 10× | 10× | 1× | Advanced isolation |
| `optimized` | 4 | 100× | 10× | 10× | 10× | Best-case scenario |

## Coherence Systems

| System | Φ₀ (m⁻¹) | ξ Range | ΔG/G Sign | Notes |
|--------|----------|---------|-----------|-------|
| Rb87 BEC | 3.65e6 | 1-100 | Negative (offset) | Lowest Φ₀ |
| Nb cavity | 2.63e7 | 1-100 | Negative (offset) | Mid Φ₀ |
| YBCO cuprate | 6.67e8 | 1-100 | Positive (offset) | Highest Φ₀, largest Δτ |

## Typical Results

### Torque Scales (Corrected)
- Newtonian: τ_N ≈ 2×10^-13 N·m
- With coherence: Δτ ≈ 1.6×10^-12 N·m (YBCO ξ=100, offset)
- ΔG/G range: [-5, +8.3]

### Feasibility (Cryo Moderate Profile)
- Feasible configs: 9/18 (integration time <24 hr)
- Best case: YBCO ξ=100 offset → **0.7 hr**
- Marginal: Nb ξ=100 offset → ~8 hr
- Challenging: Rb87 requires advanced isolation

## Computational Parameters

### Grid Resolution
- **41³**: Fast (~3 s/solve), rough convergence
- **61³**: Moderate (~10 s/solve), better but not converged (220% Δτ change from 41³)
- **81³**: Slower (~30 s/solve), recommended for quantitative claims
- **101³+**: Production-grade (use volume averaging)

### Solver Options
```python
solution = solver.solve(
    rho_func,
    phi_func,
    method='cg',      # 'cg', 'bicgstab', 'jacobi', 'gauss-seidel'
    tol=1e-8,         # Residual tolerance
    max_iter=10000    # Max iterations
)
```

### Volume Averaging (Recommended for Convergence)
```python
result = run_geometric_cavendish(
    ...,
    use_volume_average=True  # Reduces grid aliasing
)
```

## Testing

### Run All Tests
```bash
pytest tests/ -v
# 8 tests, ~18 s runtime
```

### Individual Test Suites
```bash
# Newtonian scale
pytest tests/test_newtonian_torque_scale.py -v

# Coherence invariance (5 tests)
pytest tests/test_coherence_invariance.py -v

# Volume averaging (3 tests)
pytest tests/test_volume_average.py -v
```

## Output Files

### Figures (examples/figures/)
- `feasibility_integration_times.png` - Integration time vs (ξ, system, position)
- `noise_profile_sweep.png` - Feasible configs per profile
- `optimized_vs_baseline.png` - Geometry optimization improvements

### Data (results/)
- `geometric_cavendish_sweep.json` - 18 baseline configurations
- `convergence_*.json` - (Future) Multi-grid convergence data

## Common Issues

### "No feasible configs"
→ Use cryo_moderate or better profile

### "Large convergence errors"
→ Increase grid_resolution to ≥61 and/or enable use_volume_average=True

### "Slow solves"
→ Reduce grid_resolution or use faster solver method (e.g., 'bicgstab')

### "Missing geometric_cavendish_sweep.json"
→ Run: `python examples/geometric_cavendish.py` first

## Citation

If you use this code, please cite:
```
Coherence-Gravity Coupling Framework (2025)
Loop Quantum Gravity modification to Cavendish experiment
https://github.com/arcticoder/coherence-gravity-coupling
```

## Support

- Documentation: `README.md`, `PROGRESS_SUMMARY.md`
- Tests: `tests/test_*.py`
- Issues: Open a GitHub issue with config details and error messages

---

**Quick Start:**
```bash
python examples/refined_feasibility.py --profile cryo_moderate --optimize
pytest tests/ -v
```
