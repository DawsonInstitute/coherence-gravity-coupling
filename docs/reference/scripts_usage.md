# Scripts Usage Guide

Complete documentation of all scripts in the `scripts/` directory.

---

## Analysis Scripts

### `run_analysis.py`
**Purpose**: Main analysis CLI for curvature-EM coupling constraints  
**Paper**: null_results.tex (curvature–EM coupling null-results paper)

**Usage**:
```bash
python scripts/run_analysis.py [subcommand] [options]
```

**Subcommands**:
- `sweep` - Parameter space sweep (ξ, B, materials)
- `limit` - Compute exclusion limits on κ_R
- `plot` - Generate constraint plots
- `table` - Generate analysis tables

**Example**:
```bash
python scripts/run_analysis.py sweep --xi 100 --materials rb87_bec nb_cavity --b-field 1.0 5.0 10.0
python scripts/run_analysis.py limit --b-field 10.0 --delta 1e-6
python scripts/run_analysis.py plot --output results/analysis/
```

**Make target**: `make analysis`

---

### `optimize_geometry.py`
**Purpose**: Optimize test mass geometry for maximum signal  
**Paper**: coherence_gravity_coupling.tex (main paper)

**Usage**:
```bash
python scripts/optimize_geometry.py --resolution RESOLUTION --method METHOD
```

**Options**:
- `--resolution` - Grid resolution (default: 41)
- `--method` - Optimization method: Nelder-Mead, Powell, BFGS (default: Nelder-Mead)
- `--max-iter` - Maximum iterations (default: 100)
- `--tol` - Convergence tolerance (default: 1e-6)

**Example**:
```bash
python scripts/optimize_geometry.py --resolution 41 --method Nelder-Mead
```

**Make target**: `make optimize`

---

## Figure Generation Scripts

### `generate_figures.py`
**Purpose**: Generate all figures for main coherence-gravity coupling paper  
**Paper**: coherence_gravity_coupling.tex (main paper)

**Usage**:
```bash
python scripts/generate_figures.py [--output DIR]
```

**Generates**:
- `papers/figures/convergence_analysis.pdf` - Figure 1: Convergence study
- `papers/figures/material_comparison.pdf` - Figure 2: Material comparison
- `papers/figures/landscape_YBCO_z_slice.pdf` - Figure 3: 3D potential landscape
- PNG versions of all figures

**Example**:
```bash
python scripts/generate_figures.py
python scripts/generate_figures.py --output custom/output/dir/
```

**Make target**: `make figures`

---

### `generate_bsm_plots.py`
**Purpose**: Generate BSM parameter space plots  
**Paper**: kappaR_to_BSM/curvature_em_to_bsm.tex (BSM parameter space paper)

**Usage**:
```bash
python scripts/generate_bsm_plots.py
```

**Generates**:
- `papers/kappaR_to_BSM/figures/epsilon_vs_curvature.pdf` - Dark photon mixing ε vs curvature R
- `papers/kappaR_to_BSM/figures/axion_vs_curvature.pdf` - Axion coupling g_aγγ vs curvature R
- `papers/kappaR_to_BSM/figures/curvature_amplification.pdf` - Amplification factors across environments
- PNG versions of all figures

**Example**:
```bash
python scripts/generate_bsm_plots.py
```

**Make target**: `make bsm-figures`

---

## Report Generation Scripts

### `generate_report.py`
**Purpose**: Generate consolidated analysis reports and tables  
**Paper**: null_results.tex (curvature–EM coupling null-results paper)

**Usage**:
```bash
python scripts/generate_report.py [--all] [--csv] [--markdown] [--latex]
```

**Options**:
- `--all` - Generate all output formats (CSV, Markdown, LaTeX)
- `--csv` - Generate CSV tables only
- `--markdown` - Generate Markdown tables only
- `--latex` - Generate LaTeX tables only
- `--output DIR` - Output directory (default: results/reports/)

**Generates**:
- `results/reports/csv/sweep_results.csv`
- `results/reports/csv/kappa_r_limits.csv`
- `results/reports/analysis_report.md`
- `results/reports/analysis_tables.tex`

**Example**:
```bash
python scripts/generate_report.py --all
python scripts/generate_report.py --csv --output custom/dir/
```

**Make target**: `make report`

---

### `generate_manifest.py`
**Purpose**: Generate data manifest for reproducibility  
**Paper**: All papers (data availability section)

**Usage**:
```bash
python scripts/generate_manifest.py --output FILENAME --roots ROOT1 ROOT2 ...
```

**Options**:
- `--output` - Output CSV filename (default: data_manifest.csv)
- `--roots` - Root directories to scan (default: results, papers/figures)

**Generates**:
CSV file with columns: file_path, file_type, size_bytes, sha256_hash, timestamp

**Example**:
```bash
python scripts/generate_manifest.py --output data_manifest.csv --roots results papers/figures
```

**Make target**: `make manifest`

---

## Benchmarking Scripts

### `benchmark_solver.py`
**Purpose**: Benchmark 3D Poisson solver performance  
**Paper**: coherence_gravity_coupling.tex (performance section in methods)

**Usage**:
```bash
python scripts/benchmark_solver.py [options]
```

**Options**:
- `--resolution RESOLUTION` - Grid resolution (default: 41)
- `--runs N` - Number of benchmark runs (default: 3)
- `--materials MATERIALS` - Materials to test (default: rb87_bec)
- `--xi XI` - Coupling strength (default: 100)
- `--output FILE` - Output JSON file (default: benchmark_results.json)

**Example**:
```bash
python scripts/benchmark_solver.py --resolution 61 --runs 3 --materials rb87_bec nb_cavity --xi 100
```

**Make targets**:
- `make quick-bench` - Quick benchmark at 41³ (~30s)
- `make bench` - Full benchmark at 61³ (~2min)

---

## Utility Scripts

### `verify_release.sh`
**Purpose**: Verify repository state before release  
**Paper**: All papers (release workflow)

**Usage**:
```bash
bash scripts/verify_release.sh
```

**Checks**:
- All tests passing
- All papers compile
- All figures present
- Data manifest up to date
- No uncommitted changes
- Git tags match version

**Example**:
```bash
bash scripts/verify_release.sh
```

**Make target**: `make release-verify`

---

## Summary by Paper

### coherence_gravity_coupling.tex (Main Paper)
- `generate_figures.py` - All manuscript figures
- `optimize_geometry.py` - Geometry optimization
- `benchmark_solver.py` - Performance benchmarks

### null_results.tex (Null Results Paper)
- `run_analysis.py` - Analysis CLI (sweep, limit, plot, table)
- `generate_report.py` - Consolidated tables and reports

### kappaR_to_BSM/curvature_em_to_bsm.tex (BSM Paper)
- `generate_bsm_plots.py` - All BSM parameter space plots

### Common (All Papers)
- `generate_manifest.py` - Data availability manifest
- `verify_release.sh` - Release verification

---

## Quick Reference

**Generate all figures for all papers**:
```bash
python scripts/generate_figures.py        # Main paper
python scripts/generate_bsm_plots.py      # BSM paper
python scripts/run_analysis.py plot       # Null results plots
```

**Build all papers**:
```bash
make all-papers
# or individually:
make paper          # coherence_gravity_coupling.pdf
make null-results   # null_results.pdf
make bsm-paper      # curvature_em_to_bsm.pdf
```

**Generate all reports**:
```bash
python scripts/generate_report.py --all
python scripts/generate_manifest.py
```

**Run all benchmarks**:
```bash
make quick-bench    # Quick 41³ benchmark
make bench          # Full 61³ benchmark
```

---

## See Also

- [papers/README.md](../../papers/README.md) - Paper build instructions
- [papers/submission/README.md](../../papers/submission/README.md) - Submission guidelines
- [Makefile](../../Makefile) - All make targets
- [docs/reference/QUICKREF.md](quickref/QUICKREF.md) - Quick reference guide
