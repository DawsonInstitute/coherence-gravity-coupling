# Zenodo Deposit Checklist

## Pre-Upload Verification ✅

- [x] LICENSE file present (MIT)
- [x] README.md with Quick Start
- [x] REPRODUCIBILITY.md with exact environment
- [x] data_manifest.csv with SHA256 checksums (3940 entries)
- [x] release_manifest.json with commit hash and env
- [x] All tests passing (23/23)
- [x] Manuscript PDF compiled (5 pages, 266KB)
- [x] All figures present (3 PDFs + PNGs)
- [x] Release tagged (v1.0.0)

## Files to Include in Zenodo Upload

### Core Code and Documentation
- README.md
- LICENSE
- CITATION.cff
- REPRODUCIBILITY.md
- environment.yml
- requirements.txt
- data_manifest.csv
- release_manifest.json

### Source Code
- src/ (entire directory)
- examples/
- tests/
- scripts/
- *.py (root-level scripts)
- Makefile
- pyproject.toml

### Results and Artifacts
- results/ (select essential files):
  - convergence_test_*.json (2 files, ~3KB total)
  - production_study/*.json (6 files, ~150KB total)
  - production_study/*.pdf (3 material comparison + 3 landscapes)
  - analysis/*.json (select representative)

### Manuscript
- papers/coherence_gravity_coupling.pdf ⭐
- papers/coherence_gravity_coupling.tex
- papers/coherence_gravity_coupling.bib
- papers/figures/*.pdf (3 figures)
- papers/author_config.tex

### Documentation
- docs/manuscript/*.md (6 sections)

## Zenodo Metadata

### Basic Info
- **Title**: Coherence-Modulated Gravity Coupling Framework
- **Upload type**: Software + Dataset
- **Publication date**: 2025-10-19
- **Version**: v1.0.0
- **License**: MIT

### Authors (update with real info)
```
- Family name: [Your Last Name]
  Given names: [Your First Name]
  Affiliation: [Your Institution]
  ORCID: [Your ORCID]
```

### Description
```
Complete numerical framework and experimental feasibility study for detecting 
coherence-modulated gravitational coupling via non-minimal scalar-tensor theory. 

Includes:
- 3D Poisson solver with spatially-varying G_eff
- Convergence-validated results (61³→81³→101³ resolution)
- Publication-ready manuscript with LaTeX source and figures
- Full reproducibility documentation with pinned environment
- Comprehensive test suite and data manifest with checksums

Key finding: Tabletop torsion balance can detect coherence-induced torque signals 
(τ_coh ~ 1.4 × 10⁻¹² N·m) with 0.7-24 hour integration at cryogenic temperatures.
```

### Keywords
- quantum-gravity
- scalar-tensor-theory
- torsion-balance
- experimental-feasibility
- numerical-PDE
- reproducibility
- loop-quantum-gravity
- coherence-effects

### Related Identifiers
- GitHub repository: https://github.com/DawsonInstitute/coherence-gravity-coupling
- Related preprint: [arXiv ID if/when submitted]

### Grants (if applicable)
- [Funding agency and grant number]

## Commands to Create Archive

```bash
# Create clean archive (excluding cache, pycache, git)
cd /home/echo_/Code/asciimath/
tar -czf coherence-gravity-coupling-v1.0.0.tar.gz \
  --exclude='coherence-gravity-coupling/.git' \
  --exclude='coherence-gravity-coupling/results/cache' \
  --exclude='coherence-gravity-coupling/**/__pycache__' \
  --exclude='coherence-gravity-coupling/**/*.pyc' \
  --exclude='coherence-gravity-coupling/.pytest_cache' \
  coherence-gravity-coupling/

# Verify archive size (should be ~5-10 MB)
ls -lh coherence-gravity-coupling-v1.0.0.tar.gz

# Optional: Create minimal archive with just essentials
tar -czf coherence-gravity-coupling-v1.0.0-minimal.tar.gz \
  coherence-gravity-coupling/README.md \
  coherence-gravity-coupling/LICENSE \
  coherence-gravity-coupling/CITATION.cff \
  coherence-gravity-coupling/REPRODUCIBILITY.md \
  coherence-gravity-coupling/environment.yml \
  coherence-gravity-coupling/requirements.txt \
  coherence-gravity-coupling/data_manifest.csv \
  coherence-gravity-coupling/release_manifest.json \
  coherence-gravity-coupling/papers/coherence_gravity_coupling.pdf \
  coherence-gravity-coupling/papers/figures/
```

## Post-Upload Steps

1. **Verify DOI**: Update CITATION.cff with actual Zenodo DOI
2. **Update README**: Add Zenodo badge to README.md
3. **GitHub Release**: Create GitHub release pointing to Zenodo
4. **arXiv**: Upload manuscript to arXiv if desired
5. **Journal Submission**: Submit to target journal (PRL, Nat Comm, etc.)

## Zenodo Badge (after upload)

Add to README.md:
```markdown
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)
```

## Notes

- **File size limit**: Zenodo allows up to 50 GB per deposit
- **Private vs Public**: You can create a restricted-access deposit initially
- **Communities**: Consider adding to "Quantum Gravity" or "Experimental Physics" communities
- **Versioning**: Future updates create new versions with separate DOIs

---

**Status**: Ready for upload ✅  
**Created**: 2025-10-19  
**Tag**: v1.0.0
