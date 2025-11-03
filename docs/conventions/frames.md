# Frame Conventions: Laboratory, TT, and Harmonic Gauges in EFQS

## Overview

The Extreme Field QED Simulator (EFQS) computes gravitational wave strain from electromagnetic stress-energy sources. This requires careful treatment of coordinate systems and gauge choices to ensure physical observables are correctly extracted.

## Three Key Frames

### 1. Laboratory Frame (Source-Centered)
**Coordinates**: Cartesian $(x, y, z, t)$ with origin at EM source  
**Metric**: $g_{\mu\nu} = \eta_{\mu\nu} + h_{\mu\nu}$ (weak-field expansion)  
**Usage**: 
- EM field configuration $E(x,y,z,t)$, $B(x,y,z,t)$
- Stress-energy computation $T_{\mu\nu}(x,y,z,t)$
- Quadrupole moment integration $Q_{ij}(t) = \int x^i x^j T^{00} d^3x$

**Key Point**: Laboratory frame is where sources live. No gauge choice yet - this is just the physical setup.

### 2. TT Frame (Transverse-Traceless)
**Transformation**: Apply projection operator to remove unphysical gauge modes  
**Metric**: $h_{\mu\nu}^{TT}$ with constraints:
- Transverse: $\partial_i h_{ij}^{TT} = 0$
- Traceless: $h_{ii}^{TT} = 0$
- Temporal gauge: $h_{0\mu}^{TT} = 0$

**Implementation** (`gravitational_coupling.py::strain_far_field`):
```python
if apply_TT_projection:
    n = line_of_sight / np.linalg.norm(line_of_sight)
    P_ij = delta_ij - n_i * n_j  # Transverse projector
    h_TT = P_ij Q_ij P_jk - (1/2) P_ij Q_kk delta_ij
```

**Usage**:
- Far-field GW strain observable by detectors (LIGO, LISA)
- Removes coordinate artifacts: only 2 physical polarizations (+, ×)
- Valid for $r \gg \lambda_{\text{GW}}$ (far zone)

**Key Point**: TT gauge is what detectors measure. It's the physically observable piece of $h_{\mu\nu}$.

### 3. Harmonic Gauge (Lorenz/de Donder)
**Gauge Condition**: $\partial_\mu h^{\mu\nu} = \frac{1}{2} \partial^\nu h$ (reduces Einstein equations to wave equation)  
**Wave Equation**: $\Box h_{\mu\nu} = -\frac{16\pi G}{c^4} T_{\mu\nu}$

**Usage in EFQS**:
- Not explicitly enforced in current version
- Would be needed for numerical relativity (constraint evolution)
- TT projection post-processes harmonic gauge solution

**Key Point**: Harmonic gauge simplifies field equations. TT gauge is a further specialization for radiation.

## Gauge Transformations

### Laboratory → TT
**When**: After computing $Q_{ij}(t)$, before evaluating far-field strain  
**Why**: Extract physical radiation degrees of freedom  
**How**: Apply transverse-traceless projector in direction of observer

Example (EFQS pipeline):
```python
# Laboratory frame
Q_ij = quadrupole_moment(positions, T00 * dV)

# TT projection
h_TT = strain_far_field(Q_ij, dt, R, use_tt=True, line_of_sight=np.array([0,0,1]))
```

### Harmonic → TT
**Relation**: 
- Harmonic gauge: General solution to Einstein equations
- TT gauge: Restriction to radiation modes in far zone

**Decomposition**:
$$h_{\mu\nu}^{\text{harm}} = h_{\mu\nu}^{TT} + \partial_{(\mu}\xi_{\nu)} + \text{trace terms}$$

where $\partial_{(\mu}\xi_{\nu)}$ are pure gauge (coordinate artifacts).

## Common Pitfalls

### Mistake 1: Mixing Frames in Integration
**Wrong**:
```python
Q_ij = integrate(x_i * x_j * T00_in_TT_gauge)  # ❌
```
**Right**:
```python
Q_ij = integrate(x_i * x_j * T00_in_lab_frame)  # ✅
h_TT = project_to_TT(Q_ij)
```

**Reason**: Sources $T_{\mu\nu}$ live in laboratory frame. TT projection is applied to *output* $h_{\mu\nu}$, not input.

### Mistake 2: Forgetting Line-of-Sight Dependence
**Issue**: TT projection depends on observer direction $\hat{n}$

**Example**:
```python
# Different observers see different polarizations
h_TT_z = strain_far_field(Q, dt, R, line_of_sight=[0,0,1])  # z-axis observer
h_TT_x = strain_far_field(Q, dt, R, line_of_sight=[1,0,0])  # x-axis observer
# h_TT_z ≠ h_TT_x in general!
```

### Mistake 3: Using TT Gauge in Near Zone
**Problem**: TT projection valid only for $r \gg \lambda_{\text{GW}}$

**Criterion**: 
- Far zone: $r > 10 \lambda_{\text{GW}}$ ✅ Use TT
- Near zone: $r < \lambda_{\text{GW}}$ ❌ Need full $h_{\mu\nu}$ or numerical relativity

## EFQS Configuration

### Default Settings (run_experiments.py)
```yaml
gravitational:
  enabled: true
  use_spectral_derivatives: true  # Smooth numerical derivatives
  apply_TT_projection: true        # Project to transverse-traceless
  observer_distances: [1.0, 10.0, 100.0]  # meters
```

### When to Disable TT Projection
- **Near-zone studies**: Source region ($r < \lambda$)
- **Numerical relativity**: Solving full Einstein equations
- **Gauge artifact analysis**: Comparing different gauge choices

### When to Keep TT Projection (Default)
- **Far-field GW predictions**: Detector sensitivity studies
- **Astrophysical recasts**: Mapping lab → compact objects
- **Frequency-domain analysis**: QNM extraction, spectral diagnostics

## Coherence-Gravity Coupling Extensions

### Additional Gauge Considerations

**Coherence Field $\Phi$**:
- Lives in laboratory frame (like $T_{\mu\nu}$)
- Couples to curvature $R$ via $\xi R \Phi^2$
- Gauge-invariant observable: $\Delta \tau = \int (G_{\text{eff}} - G) T^{00} / c^4 d^3x$

**Torsion Proxy** (from `torsion_dof.py`):
- Antisymmetric part of $\nabla\Phi \otimes \nabla\Phi$ mimics torsion
- Gauge behavior: Transforms like stress-energy (laboratory frame quantity)
- Observable: Duality-breaking integral $\int E \cdot (\nabla \times A_{\text{torsion}}) d^3x$

**Robin BC** (from `robin_bc_poisson.py`):
- Boundary condition: $\alpha \Phi + \beta \partial_n \Phi = 0$
- Gauge-independent (coordinate-invariant at boundary)
- Observable: $\xi$-dependent vacuum energy $E_\xi \propto (1/4 - \xi)$

## Quick Reference Table

| Quantity | Frame | Gauge | Observable? |
|----------|-------|-------|-------------|
| $E(x,t), B(x,t)$ | Lab | N/A | Yes (local) |
| $T_{\mu\nu}$ | Lab | N/A | Yes (local) |
| $Q_{ij}(t)$ | Lab | N/A | No (integrated) |
| $h_{\mu\nu}$ (general) | Any | Harmonic/etc | No (gauge-dependent) |
| $h_{\mu\nu}^{TT}$ | Lab | TT | Yes (far-field) |
| $\Phi$ | Lab | N/A | Yes (coherence) |
| $\Delta \tau$ | Lab | N/A | Yes (torque) |

## Further Reading

- Misner, Thorne, Wheeler (1973): *Gravitation*, Ch. 18 (TT gauge)
- Maggiore (2008): *Gravitational Waves*, Vol. 1, Ch. 1.4 (gauge choices)
- EFQS docs: `extreme-field-qed-simulator/docs/gravitational_coupling_theory.md`
- Coherence-gravity: `coherence-gravity-coupling/docs/EFQS_NEXT_STEPS.md`

## Contact

Questions? See `coherence-gravity-coupling/README.md` or open an issue.
