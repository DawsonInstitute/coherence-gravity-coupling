# Installation Guide

## Quick Installation

### Prerequisites
- Python 3.8+ 
- Git
- LaTeX (for building papers)
- C++ compiler (for optional accelerated solvers)

### Environment Setup

Choose one of the following methods:

#### Option 1: Conda (Recommended)
```bash
git clone https://github.com/DawsonInstitute/coherence-gravity-coupling.git
cd coherence-gravity-coupling
conda env create -f environment.yml
conda activate cohgrav
```

#### Option 2: Virtual Environment
```bash
git clone https://github.com/DawsonInstitute/coherence-gravity-coupling.git
cd coherence-gravity-coupling
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

#### Option 3: Development Install
```bash
git clone https://github.com/DawsonInstitute/coherence-gravity-coupling.git
cd coherence-gravity-coupling
pip install -e .  # Editable install
```

## Dependencies

### Core Dependencies
- `numpy>=1.21.0` - Numerical computations
- `scipy>=1.7.0` - Scientific computing, optimization
- `matplotlib>=3.5.0` - Plotting and visualization
- `pytest>=6.0` - Testing framework
- `sympy>=1.8` - Symbolic mathematics
- `h5py>=3.1` - HDF5 data storage

### Optional Dependencies
- `numba>=0.54` - JIT compilation for performance
- `mpi4py>=3.1` - MPI parallelization
- `vtk>=9.0` - 3D visualization
- `jupyter>=1.0` - Interactive notebooks

### LaTeX Dependencies (for papers)
```bash
# Ubuntu/Debian
sudo apt-get install texlive-latex-extra texlive-fonts-recommended texlive-bibtex-extra

# macOS (with Homebrew)
brew install --cask mactex

# Windows
# Download and install MiKTeX from https://miktex.org/
```

## Verification

Test your installation:

```bash
# Quick smoke tests (~90 seconds)
pytest -q

# Generate sample figures
python scripts/generate_figures.py

# Build main paper
cd papers
make paper
```

Expected outputs:
- All tests pass
- Figures generated in `papers/figures/`
- `coherence_gravity_coupling.pdf` compiled successfully

## Performance Notes

Runtime estimates (Intel i7-10700K, 32GB RAM):
- 41³ grid: ~3-5s per solve
- 61³ grid: ~5-8s per solve  
- 81³ grid: ~20-30s per solve
- 101³ grid: ~1-2min per solve

For production runs, 61³ or higher is recommended for converged results.

## Troubleshooting

### Common Issues

**Import errors**: Ensure virtual environment is activated and dependencies installed
```bash
source .venv/bin/activate  # or conda activate cohgrav
pip list  # Verify packages installed
```

**LaTeX compilation fails**: Install missing LaTeX packages
```bash
# Find missing package
grep "File.*not found" paper.log
# Install via package manager or texlive
```

**Performance issues**: Enable numerical acceleration
```bash
pip install numba  # JIT compilation
pip install mpi4py  # MPI parallelization (if available)
```

**Test failures**: Check Python version and dependencies
```bash
python --version  # Should be 3.8+
pytest -v  # Verbose test output for debugging
```

### Getting Help

1. Check [Issues](https://github.com/DawsonInstitute/coherence-gravity-coupling/issues) for known problems
2. Search documentation in `docs/`
3. Run `make help` for available commands
4. Contact: rsherrington@dawsoninstitute.org