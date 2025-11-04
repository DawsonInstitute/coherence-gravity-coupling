# Build System and Commands

## Make Targets

The project uses a Makefile for common development tasks:

```bash
make help          # Show all available targets
make test          # Run full test suite (23 tests, ~90s)
make quick-bench   # Quick benchmark at 41³ (~30s)
make figures       # Generate figures for main paper (coherence_gravity_coupling.tex)
make paper         # Build coherence_gravity_coupling.pdf (main paper)
make null-results  # Build null_results.pdf (curvature-EM constraints paper)
make bsm-paper     # Build curvature_em_to_bsm.pdf (BSM parameter space paper)
make analysis      # Analysis CLI for null_results.tex (curvature-EM bounds)
make optimize      # Run geometry optimization (for main paper experiments)
make cache-info    # Show cache statistics
make cache-clean   # Clear all cached results
```

## Paper Building

The project generates three main papers:

### 1. Main Paper: Coherence-Gravity Coupling
```bash
make paper
# or manually:
cd papers
pdflatex coherence_gravity_coupling.tex && bibtex coherence_gravity_coupling \
   && pdflatex coherence_gravity_coupling.tex && pdflatex coherence_gravity_coupling.tex
```

### 2. Null Results Paper: Curvature-EM Constraints  
```bash
make null-results
# or manually:
cd papers
pdflatex null_results.tex && bibtex null_results \
   && pdflatex null_results.tex && pdflatex null_results.tex
```

### 3. BSM Parameter Space Paper: κ_R to Dark Photon/Axion
```bash
make bsm-paper
# or manually:
cd papers/kappaR_to_BSM
pdflatex curvature_em_to_bsm.tex && pdflatex curvature_em_to_bsm.tex
```

**Why multiple pdflatex runs?** LaTeX requires 2-3 passes to resolve cross-references, citations, and bibliography. The sequence (pdflatex → bibtex → pdflatex → pdflatex) ensures all references are correct.

## Script Commands

### Figure Generation
```bash
# Main paper figures → papers/figures/*.pdf,*.png
python scripts/generate_figures.py

# BSM parameter space plots
python scripts/generate_bsm_plots.py

# Analysis reports (CSV/Markdown/LaTeX)
python scripts/generate_report.py --all
```

### Analysis Tools
```bash
# Curvature-EM coupling analysis
python scripts/analysis_cli.py --kappa-R 5e17 --B-field 10 --delta 1e-6

# Convergence study
python scripts/convergence_study.py --grids 41,61,81 --materials YBCO,Al

# Geometry optimization
python scripts/optimize_geometry.py --target-sensitivity 1e-12
```

### Testing
```bash
# Full test suite
pytest

# Quick smoke tests
pytest -q

# Specific test categories
pytest tests/test_solver.py      # Solver validation
pytest tests/test_physics.py     # Physics consistency
pytest tests/test_integration.py # End-to-end tests

# Coverage report
pytest --cov=src --cov-report=html
```

## Development Commands

### Code Quality
```bash
# Format code
black src/ tests/ scripts/

# Lint
flake8 src/ tests/ scripts/

# Type checking
mypy src/
```

### Performance Profiling
```bash
# Profile solver
python scripts/profile_solver.py --grid 61

# Memory usage
python scripts/memory_profile.py

# Benchmark suite
python scripts/benchmark.py --full
```

## Cache Management

The solver uses intelligent caching to speed up repeated calculations:

```bash
# Show cache statistics
make cache-info

# Clear all cached results
make cache-clean

# Clear specific cache types
python scripts/cache_manager.py --clear-solver
python scripts/cache_manager.py --clear-figures
```

## Advanced Build Options

### Parallel Compilation
```bash
# Use multiple cores for LaTeX
export LATEX_JOBS=4
make paper

# Parallel testing
pytest -n auto  # Requires pytest-xdist
```

### Debug Builds
```bash
# Enable debug mode
export DEBUG=1
python scripts/generate_figures.py

# Verbose solver output
export SOLVER_DEBUG=1
python examples/basic_solver.py
```

### Custom Configurations
```bash
# Use custom config file
python scripts/analysis_cli.py --config my_config.yaml

# Override default parameters
python scripts/generate_figures.py --xi 100 --grid 61
```

## Output Files

### Generated Figures
- `papers/figures/convergence_analysis.pdf` - Grid convergence study
- `papers/figures/material_comparison.pdf` - Material properties
- `papers/figures/landscape_YBCO_z_slice.pdf` - Solution landscape
- `papers/figures/epsilon_vs_curvature.pdf` - BSM parameter space
- `papers/figures/curvature_amplification.pdf` - Amplification factors

### Analysis Reports
- `results/reports/convergence_report.md` - Convergence analysis
- `results/reports/material_analysis.csv` - Material comparison data
- `results/reports/exclusion_limits.tex` - LaTeX tables for papers

### Compiled Papers
- `papers/coherence_gravity_coupling.pdf` - Main paper (5 pages)
- `papers/null_results.pdf` - Null results paper
- `papers/kappaR_to_BSM/curvature_em_to_bsm.pdf` - BSM paper

### Data Files
- `results/analysis/*.h5` - Numerical solution data
- `results/cache/*.pkl` - Cached computations
- `results/benchmarks/*.json` - Performance metrics