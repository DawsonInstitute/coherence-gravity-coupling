# Session Summary: Local Dev Workflow & Domain Study

**Date**: October 17, 2025  
**Tasks**: Domain/BC sensitivity study + local dev workflow setup

## Completed Work

### 1. Fixed Test Robustness ✅
**Files**: `tests/test_coherence_invariance.py`

Fixed two fragile tests:
- `test_xi_zero_invariance`: Added epsilon guard to avoid divide-by-zero when τ≈0
- `test_monotonicity_with_xi`: Allow near-zero ties when effect is below numerical precision

**Result**: All 23 tests now pass with 5 warnings (benign return-not-none style issues)

### 2. Domain Padding Sensitivity Study ✅
**Files**: `examples/domain_bc_sweep.py`

Created automated sweep tool to measure Δτ sensitivity vs domain padding:
- Computes minimum enclosing domain from geometry
- Tests multiple padding factors (default: 1.2×, 1.5×, 2.0×, 2.5×)
- Reports mean, std, and max variation in Δτ
- Recommends padding that keeps variation < 5%

**Initial findings** (41³, YBCO, ξ=100):
- Min domain: 0.260 m (tight bounding box)
- 1.5× padding: Δτ = -2.911e-13 N·m
- 2.0× padding: Δτ = -2.279e-13 N·m
- **Variation: 12.2%** (exceeds 5% target)
- Conclusion: Need padding ≥ 2.5× or 3.0× for stability

**Usage**:
```bash
python examples/domain_bc_sweep.py --resolution 41 --padding 1.5 2.0 2.5 3.0
make domain-sweep  # Convenience target
```

### 3. Local Development Workflow ✅
**Files**: `tox.ini`, `.pre-commit-config.yaml`, `Makefile`, `pyproject.toml`, `DEV_WORKFLOW.md`

Implemented complete local-only tooling (no GitHub Actions):

**Makefile targets**:
- `make test` - Full test suite (23 tests, ~90s)
- `make quick-bench` - Quick benchmark at 41³ (~30s)
- `make bench` - Full benchmark at 61³ (~2min)
- `make domain-sweep` - Domain padding study
- `make lint` - Flake8 linter
- `make format` - Auto-format with black + isort
- `make format-check` - Check formatting (no changes)
- `make install-dev` - Install dev dependencies + pre-commit hooks
- `make clean` - Remove generated files

**Tox environments**:
- `tox -e py313` - Run tests in isolated env
- `tox -e lint` - Flake8
- `tox -e format-check` - Black + isort check
- `tox -e format` - Auto-format
- `tox -e quick-bench` - Quick benchmark
- `tox` - Run all environments

**Pre-commit hooks** (auto-run on `git commit`):
- `black` - Code formatting (line-length=100)
- `isort` - Import sorting (black-compatible)
- `flake8` - Linting (max-line-length=100)
- `trailing-whitespace` - Fix trailing spaces
- `end-of-file-fixer` - Ensure final newline
- `check-yaml` - YAML syntax
- `check-merge-conflict` - Detect merge markers
- `check-added-large-files` - Prevent large commits

**Installation**:
```bash
make install-dev  # Installs pytest, black, isort, flake8, tox, pre-commit + hooks
```

## Test Results

**Full suite**: 23 passed, 5 warnings (100% success rate)

Warnings are benign (pytest prefers `assert` over `return` in tests) and can be cleaned later.

## Domain Study Findings

Current default domain = 0.6 m gives adequate padding for standard geometry but shows sensitivity:
- Tight padding (1.5×): Stronger signal but boundary effects
- Loose padding (2.0×): Weaker signal but cleaner
- Variation at 41³: ~12% between 1.5× and 2.0×

**Recommendation**: Run full sweep [1.5, 2.0, 2.5, 3.0, 4.0] at 61³ to find stable default.

## Files Created/Modified

**New files**:
1. `examples/domain_bc_sweep.py` - Domain padding sensitivity tool
2. `tox.ini` - Tox configuration for isolated testing
3. `.pre-commit-config.yaml` - Pre-commit hooks config
4. `pyproject.toml` - Modern Python packaging + tool config
5. `DEV_WORKFLOW.md` - Developer guide for local workflow

**Modified files**:
1. `tests/test_coherence_invariance.py` - Robustness fixes
2. `Makefile` - Updated targets for new workflow

## Next Steps

### Immediate (Complete domain study)

1. **Run full domain sweep at 61³** to find stable padding:
   ```bash
   python examples/domain_bc_sweep.py --resolution 61 --padding 1.5 2.0 2.5 3.0 4.0 --xi 100
   ```
   - Expect tighter convergence at higher resolution
   - Target: Find smallest padding with <5% variation
   - Update default domain_size in `run_geometric_cavendish()`

2. **Implement Neumann BC option** (optional):
   - Add `bc_type` parameter to `Poisson3DSolver`
   - Options: 'dirichlet' (φ=0), 'neumann' (∂φ/∂n=0)
   - Neumann may reduce boundary effects for isolated systems
   - Test both BC types in sweep

3. **Document recommended defaults** in README:
   - Optimal padding factor
   - BC type choice
   - Resolution guidelines

### Medium Priority (Result caching)

4. **Implement result cache**:
   - Hash: (grid_res, domain_size, xi, Phi0, geom_params)
   - Storage: Compressed NPZ for fields + JSON metadata
   - Location: `results/cache/`
   - Wire into `run_geometric_cavendish(cache=True)`
   - Benefits: Skip repeated solves in sweeps

5. **Add cache management**:
   - `make cache-clean` - Clear cache
   - `make cache-info` - Show cache size/stats
   - Cache expiration policy (e.g., 30 days)

### Optional Enhancements

6. **Physics validation tests**:
   - Verify Δτ ∝ Φ₀ in weak-coupling limit
   - Test symmetric geometry gives τ ≈ 0
   - Energy conservation check

7. **Benchmark plotting**:
   - Add `benchmark_solver.py --plot` flag
   - Generate performance vs resolution plot
   - Compare solver methods visually

8. **Quick vs accurate modes**:
   - Quick: 41³, tol=1e-6, cg+diagonal
   - Accurate: 61³ or 81³, tol=1e-8, cg+amg
   - CLI flag: `--mode quick|accurate`

## Developer Experience

**Local workflow is now complete**:
- ✅ Run tests: `make test`
- ✅ Format code: `make format` (or let pre-commit do it)
- ✅ Lint: `make lint`
- ✅ Benchmark: `make quick-bench`
- ✅ Domain study: `make domain-sweep`
- ✅ Pre-commit hooks: Auto-run on commit

**No cloud dependency** - everything runs locally.

## Command Reference

```bash
# Setup
make install-dev          # One-time setup

# Daily workflow
make test                 # Run tests
make format               # Format code
git commit                # Pre-commit hooks auto-format

# Benchmarking
make quick-bench          # Fast benchmark (41³)
make bench                # Full benchmark (61³)
make domain-sweep         # Domain study

# Tox (isolated)
tox                       # Run all checks
tox -e py313              # Just tests
tox -e lint               # Just linting

# Pre-commit (manual)
pre-commit run --all-files
```

## Status

**Completed**:
- ✅ Solver performance improvements
- ✅ Test robustness fixes
- ✅ Domain sweep tool
- ✅ Local dev workflow (tox, pre-commit, Makefile)
- ✅ Developer documentation

**In Progress**:
- 🔄 Domain study (need full sweep at 61³)

**Pending**:
- ⏳ Result caching and provenance
- ⏳ Neumann BC implementation
- ⏳ Optional enhancements (plotting, validation tests)

All tools are ready to use. The primary remaining task is running the comprehensive domain sweep at 61³ resolution to finalize the recommended padding factor.
