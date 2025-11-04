# API Reference

Complete documentation for the coherence-gravity-coupling codebase.

## Core Analysis Modules

### mass_dependent_dark_photon_mixing.py

**Primary Class:** `MassiveDarkPhotonCurvature`

Maps laboratory κ_R constraints to dark photon effective mixing parameters, including mass-dependent suppression and curvature amplification.

```python
from src.analysis.mass_dependent_dark_photon_mixing import MassiveDarkPhotonCurvature

analyzer = MassiveDarkPhotonCurvature()

# Key methods:
effective_mixing_massive(kappa_R, R, mass_U, k_photon)  # ε_eff calculation
curvature_amplification_factor(R_lab, R_astro)          # Environment scaling
mass_scan(mass_range, kappa_R, R)                       # Mass dependence
curvature_scan(R_range, kappa_R, mass_U)               # Curvature dependence
```

**Physics Implementation:**
- Modified dispersion relation: `ω² = k² + M_U²`
- Effective mixing: `ε_eff = C_ε κ_R R / √(k² + M_U²)`
- Validation against Jorge et al. constraints

**Key Results:**
- Laboratory: ε_eff ~ 10⁻¹⁴ (for κ_R ~ 10¹⁷ m², R ~ 10⁻²⁶ m⁻²)
- Magnetar: ~10²⁰× amplification at R ~ 10⁻⁶ m⁻²
- Mass threshold: M_U < keV for significant laboratory sensitivity

---

### operator_mixing_kappa_alpha.py

**Primary Class:** `OperatorMixing`

Computes mixing coefficients between κ_R R F² and α L F² operators via field redefinitions in curved spacetime.

```python
from src.analysis.operator_mixing_kappa_alpha import OperatorMixing

mixer = OperatorMixing()

# Key methods:
einstein_frame_transformation(kappa_R, alpha)     # Frame transformation
field_redefinition_mixing()                       # C_α coefficient
weak_field_expansion(h_metric, order=2)          # Perturbative analysis
decoupling_check(kappa_R, alpha, mass_scale)     # Decoupling test
```

**Physics Implementation:**
- Einstein frame action: `S = ∫ √(-g) [R/(16πG) + L_matter + κ_R R F² + α L F²]`
- Field redefinition: `g_μν → g_μν + δg_μν[F]`
- Mixing coefficient: `C_α ~ O(1)` (no suppression)

**Key Results:**
- Operators do NOT decouple at linear order
- C_α ~ 1 in weak field limit
- Backreaction becomes important at κ_R ~ M_Pl⁻²

---

### combined_kappa_alpha_constraints.py

**Primary Class:** `CombinedConstraint`

Joint analysis of laboratory κ_R and astrophysical α constraints with correlation models.

```python
from src.analysis.combined_kappa_alpha_constraints import CombinedConstraint

constraint = CombinedConstraint()

# Key methods:
chi_squared_combined(kappa_R, alpha, correlation_model)  # Joint χ²
constraint_ellipse(confidence_level)                     # 2D boundaries
cancellation_scenario(kappa_R, alpha)                   # Signal cancellation
UV_correlation(kappa_R, model='linear')                 # α prediction
```

**Correlation Models:**
- `'independent'`: No correlation (conservative)
- `'linear'`: α = C × κ_R + noise  
- `'string_theory'`: α ~ κ_R / M_string²

**Key Results:**
- Laboratory dominates: σ_κ / σ_α ~ 10⁻¹¹
- No cancellation for independent observables
- Strong correlations enable discrimination

---

### joint_posterior_analysis.py

**Primary Class:** `JointPosteriorAnalysis`

Bayesian marginalization of (κ_R, α) parameter space with UV completion priors.

```python
from src.analysis.joint_posterior_analysis import JointPosteriorAnalysis

posterior = JointPosteriorAnalysis()

# Key methods:
log_likelihood(kappa_R, alpha)                    # Gaussian likelihoods
log_prior(kappa_R, alpha, uv_model)              # UV completion priors
sample_posterior(uv_model, n_samples)            # MCMC/importance sampling
marginalized_kappa_R_constraint(uv_model)        # 1D marginalization
```

**UV Models:**
- `'uniform'`: Flat prior (model-independent)
- `'linear'`: Correlated prior with scatter
- `'string_theory'`: Theoretical scaling prediction

**Key Results:**
- Marginalized constraints depend strongly on UV correlations
- String theory models predict large α/κ_R ratios
- Bayesian evidence favors specific correlation structures

---

## Constraint Analysis Modules

### laboratory_constraints.py

Laboratory electromagnetic measurements constraining κ_R.

```python
from src.constraints.laboratory_constraints import LaboratoryConstraints

lab = LaboratoryConstraints()

# Key methods:
constraint_kappa_R(B_field, curvature, precision)  # Primary bound
systematic_uncertainties()                          # Error budget
future_sensitivity(next_generation=True)           # Projections
```

**Current Bound:** κ_R < 5×10¹⁷ m² (95% CL)  
**Systematic Errors:** < 10% for B ~ 10 T experiments  
**Future Reach:** ~10¹⁵ m² with next-generation magnetometry

---

### astrophysical_constraints.py

Curvature amplification in compact object environments.

```python
from src.constraints.astrophysical_constraints import AstrophysicalConstraints

astro = AstrophysicalConstraints()

# Key methods:
curvature_profiles(object_type)              # R(r) for compact objects
amplification_factor(object_type)            # Enhancement over lab
observational_prospects(instrument)          # Detection feasibility
```

**Objects:** White dwarfs (10⁸×), neutron stars (10¹⁶×), black holes (10²⁰×)  
**Instruments:** X-ray polarimetry, gravitational wave detectors  
**Status:** Conceptual frameworks; observations in development

---

## Utility Modules

### physics_constants.py

```python
from src.utils.physics_constants import PhysicsConstants

const = PhysicsConstants()

# Available constants:
const.c              # Speed of light (m/s)
const.hbar           # Reduced Planck constant (J⋅s)  
const.G_Newton       # Gravitational constant (m³/kg⋅s²)
const.epsilon_0      # Vacuum permittivity (F/m)
const.mu_0           # Vacuum permeability (H/m)
const.M_planck       # Planck mass (kg)
const.l_planck       # Planck length (m)
```

### numerical_methods.py

```python
from src.utils.numerical_methods import NumericalMethods

num = NumericalMethods()

# Integration and optimization:
num.adaptive_integrate(func, bounds, tolerance)     # Adaptive quadrature
num.chi_squared_minimize(model, data, weights)     # Parameter fitting
num.bootstrap_uncertainty(data, statistic, n_boot) # Error estimation
num.monte_carlo_sample(distribution, n_samples)    # Random sampling
```

### validation_tools.py

```python
from src.utils.validation_tools import ValidationTools

val = ValidationTools()

# Physics checks:
val.dimensional_analysis(expression, expected_units)  # Unit checking
val.symmetry_check(transformation, invariant)        # Symmetry validation
val.convergence_test(numerical_result, analytical)   # Convergence analysis
val.physical_bounds_check(parameter, bounds)         # Range validation
```

---

## Plotting Modules

### figure_generators.py

```python
from src.plotting.figure_generators import FigureGenerators

fig_gen = FigureGenerators()

# Standard plots:
fig_gen.constraint_plot(parameter_space, bounds)         # Parameter constraints
fig_gen.amplification_plot(curvatures, factors)         # Curvature scaling
fig_gen.mass_dependence_plot(masses, mixing_values)     # Mass scans
fig_gen.correlation_plot(x_data, y_data, model)         # Parameter correlations
```

**Output Formats:** PDF (publication), PNG (presentations), SVG (web)  
**Styling:** Publication-ready with LaTeX fonts and physics notation  
**Customization:** Configurable colors, markers, and annotations

### visualization_tools.py

```python
from src.plotting.visualization_tools import VisualizationTools

viz = VisualizationTools()

# Advanced visualizations:
viz.corner_plot(samples, parameters)              # Parameter posteriors
viz.contour_plot(x, y, z, levels)                # 2D parameter space
viz.animation_sequence(time_series, parameter)   # Time evolution
viz.interactive_plot(data, widgets)              # Jupyter interaction
```

---

## Testing Framework

### Physics Validation

All major modules include comprehensive `validate_*()` functions:

```python
# Example validation pattern
def validate_mass_dependent_mixing():
    """Comprehensive physics validation."""
    
    # 1. Dimensional analysis
    assert check_dimensions(ε_eff, "dimensionless")
    
    # 2. Limiting behavior
    assert abs(ε_eff(M_U=0) - ε_eff_massless) < 1e-10  # Massless limit
    assert ε_eff(M_U→∞) < 1e-15                         # Decoupling limit
    
    # 3. Literature comparison
    assert constraint_matches_Jorge_et_al(tolerance=0.1)
    
    # 4. Numerical consistency
    assert integration_converged(relative_error=1e-8)
    
    print("✅ All validation checks passed")
```

### Test Coverage

- **Unit tests:** Individual function validation
- **Integration tests:** Module interaction testing  
- **Physics tests:** Theoretical consistency checks
- **Regression tests:** Results stability monitoring
- **Performance tests:** Computational efficiency

**Run Tests:**
```bash
pytest tests/ -v                    # Full test suite
pytest tests/test_physics.py -x     # Physics validation only
pytest tests/ --cov=src --cov-report=html  # Coverage report
```

---

## Configuration

### Default Parameters

```python
# Analysis defaults (src/config/defaults.py)
DEFAULT_KAPPA_R_BOUND = 5e17      # m² (laboratory constraint)
DEFAULT_CURVATURE_LAB = 1e-26     # m⁻² (Earth's curvature)
DEFAULT_B_FIELD = 10.0            # T (laboratory magnetic field)
DEFAULT_MASS_RANGE = [1e-6, 1e-3] # eV (dark photon mass scan)

# Numerical defaults
INTEGRATION_TOLERANCE = 1e-10     # Relative error target
MONTE_CARLO_SAMPLES = 100000      # Default sample size
CONFIDENCE_LEVEL = 0.95           # Statistical confidence
```

### Custom Configuration

```python
# Override defaults
from src.config.configuration import Configuration

config = Configuration()
config.set_parameter('kappa_R_bound', 1e16)     # Tighter constraint
config.set_parameter('curvature_lab', 5e-27)    # Different environment
config.save_configuration('custom_analysis.yaml')
```

---

## Examples & Tutorials

### Basic Analysis Example

```python
"""
Example: Laboratory constraint to astrophysical prediction
"""
from src.analysis.mass_dependent_dark_photon_mixing import MassiveDarkPhotonCurvature

# Initialize analyzer
analyzer = MassiveDarkPhotonCurvature()

# Laboratory constraint
kappa_R_lab = 5e17  # m² (95% CL bound)
R_lab = 1e-26       # m⁻² (Earth curvature)
mass_U = 1e-3       # eV (dark photon mass)

# Effective mixing in laboratory
epsilon_lab = analyzer.effective_mixing_massive(
    kappa_R=kappa_R_lab, 
    R=R_lab, 
    mass_U=mass_U,
    k_photon=1e7  # Optical wavelength
)
print(f"Laboratory mixing: ε_eff = {epsilon_lab:.2e}")

# Predict magnetar enhancement
R_magnetar = 1e-6   # m⁻² (magnetar surface)
amplification = analyzer.curvature_amplification_factor(R_lab, R_magnetar)
epsilon_magnetar = epsilon_lab * amplification

print(f"Magnetar enhancement: {amplification:.1e}×")
print(f"Magnetar mixing: ε_eff = {epsilon_magnetar:.2e}")
```

### Advanced Correlation Analysis

```python
"""
Example: Joint posterior with UV correlations
"""
from src.analysis.joint_posterior_analysis import JointPosteriorAnalysis

# Initialize Bayesian analysis
posterior = JointPosteriorAnalysis()

# Sample joint posterior with string theory prior
kappa_samples, alpha_samples = posterior.sample_posterior(
    uv_model='string_theory',
    correlation_strength=0.8,
    n_samples=50000
)

# Marginalized κ_R constraint
kappa_lower, kappa_upper = posterior.marginalized_kappa_R_constraint(
    uv_model='string_theory',
    confidence_level=0.95
)

print(f"Marginalized constraint: κ_R ∈ [{kappa_lower:.1e}, {kappa_upper:.1e}] m²")
print(f"Improvement over individual: {2*5e17/(kappa_upper-kappa_lower):.1f}×")
```

### Custom Analysis Pipeline

```python
"""
Example: Custom analysis with user-defined parameters
"""
from src.analysis.mass_dependent_dark_photon_mixing import MassiveDarkPhotonCurvature
from src.constraints.laboratory_constraints import LaboratoryConstraints
from src.plotting.figure_generators import FigureGenerators

# Custom analysis pipeline
def custom_analysis(B_field_range, mass_range, output_dir):
    """Custom dark photon analysis pipeline."""
    
    # Initialize modules
    analyzer = MassiveDarkPhotonCurvature()
    lab_constraints = LaboratoryConstraints()
    plotter = FigureGenerators()
    
    results = {}
    
    # Scan magnetic field strengths
    for B_field in B_field_range:
        # Update laboratory constraint
        kappa_bound = lab_constraints.constraint_kappa_R(
            B_field=B_field,
            curvature=1e-26,
            precision='next_generation'
        )
        
        # Mass scan for this B field
        epsilon_values = []
        for mass in mass_range:
            epsilon = analyzer.effective_mixing_massive(
                kappa_R=kappa_bound,
                R=1e-26,
                mass_U=mass,
                k_photon=1e7
            )
            epsilon_values.append(epsilon)
        
        results[B_field] = {
            'kappa_bound': kappa_bound,
            'masses': mass_range,
            'epsilon_values': epsilon_values
        }
        
        # Generate plot for this B field
        fig = plotter.mass_dependence_plot(
            masses=mass_range,
            mixing_values=epsilon_values,
            title=f"B = {B_field:.1f} T"
        )
        fig.savefig(f"{output_dir}/mass_scan_B{B_field:.1f}T.pdf")
    
    return results

# Run custom analysis
B_range = [5.0, 10.0, 20.0]  # Tesla
mass_range = np.logspace(-6, -2, 100)  # eV
results = custom_analysis(B_range, mass_range, "output/")
```

---

## Development Guidelines

### Code Style

- **PEP 8 compliance** with 88-character line limit
- **Type hints** for all public functions
- **Docstring format:** NumPy style with physics context
- **Physics notation:** Use standard symbols (κ_R, ε_eff, etc.)

### Adding New Modules

1. **Create module:** `src/analysis/new_module.py`
2. **Add validation:** `validate_new_module()` function
3. **Write tests:** `tests/test_new_module.py`
4. **Update documentation:** Add to this API reference
5. **Physics review:** Validate against literature

### Performance Optimization

- **Vectorization:** Use NumPy arrays for parameter scans
- **Caching:** Cache expensive computations with `@lru_cache`
- **Profiling:** Use `cProfile` for bottleneck identification
- **Parallel processing:** Use `multiprocessing` for independent calculations

---

## Troubleshooting

### Common Issues

**Import errors:**
```python
# Ensure proper Python path
import sys
sys.path.append('/path/to/coherence-gravity-coupling/src')
```

**Numerical convergence:**
```python
# Increase integration tolerance
from src.utils.numerical_methods import NumericalMethods
num = NumericalMethods(tolerance=1e-12)
```

**Memory issues with large scans:**
```python
# Use chunked processing
def chunked_mass_scan(mass_range, chunk_size=1000):
    for i in range(0, len(mass_range), chunk_size):
        chunk = mass_range[i:i+chunk_size]
        yield process_mass_chunk(chunk)
```

### Performance Optimization

**Slow parameter scans:**
- Use vectorized operations instead of loops
- Cache intermediate results with `functools.lru_cache`
- Consider parallel processing with `multiprocessing.Pool`

**Memory usage:**
- Process data in chunks for large parameter spaces
- Use `numpy.memmap` for large arrays
- Clear intermediate variables with `del variable`

### Physics Validation Failures

**Dimensional analysis errors:**
- Check unit consistency in all expressions
- Verify conversion factors between unit systems
- Use `astropy.units` for automated checking

**Literature comparison mismatches:**
- Verify parameter definitions and conventions
- Check for different approximation schemes
- Account for experimental vs theoretical uncertainties

---

## Contact & Support

**Code issues:** Open GitHub issue with minimal reproducible example  
**Physics questions:** Contact project lead or open discussion issue  
**Feature requests:** Describe use case and physics motivation  
**Performance problems:** Include profiling output and system specs

**Response times:**
- Bug reports: 24-48 hours
- Physics discussions: 1-2 weeks
- Feature requests: Review at monthly team meetings

---

*This API reference covers version 1.0 of the coherence-gravity-coupling framework. For the latest updates, see the [GitHub repository](https://github.com/DawsonInstitute/coherence-gravity-coupling).*