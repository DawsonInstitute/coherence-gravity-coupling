.PHONY: help test quick-bench bench lint format format-check clean install-dev domain-sweep cache-info cache-clean analysis
.PHONY: env figures paper null-results bsm-paper bsm-figures all-papers manifest release-verify

help:
	@echo "coherence-gravity-coupling development targets"
	@echo ""
	@echo "Testing:"
	@echo "  make test          - Run full test suite (23 tests, ~90s)"
	@echo "  make quick-bench   - Quick benchmark at 41³ (2 runs, ~30s)"
	@echo "  make bench         - Full benchmark at 61³ (3 runs, ~2min)"
	@echo "  make domain-sweep  - Domain padding sensitivity study"
	@echo ""
	@echo "Analysis:"
	@echo "  make analysis      - Run analysis CLI (prints usage; see --help for subcommands)"
	@echo "  make optimize      - Run geometry optimization"
	@echo ""
	@echo "Code quality:"
	@echo "  make lint          - Run flake8 linter"
	@echo "  make format        - Auto-format code (black + isort)"
	@echo "  make format-check  - Check formatting without changes"
	@echo ""
	@echo "Cache management:"
	@echo "  make cache-info    - Show cache statistics"
	@echo "  make cache-clean   - Clear all cached results"
	@echo ""
	@echo "Setup:"
	@echo "  make install-dev   - Install development dependencies"
	@echo "  make clean         - Remove generated files"
	@echo "  make env           - Create/activate conda env (cohgrav)"
	@echo ""
	@echo "Papers:"
	@echo "  make figures       - Generate figures for main paper (coherence_gravity_coupling.tex)"
	@echo "  make paper         - Build main paper PDF (coherence_gravity_coupling.tex)"
	@echo "  make null-results  - Build null results paper PDF (null_results.tex)"
	@echo "  make bsm-figures   - Generate BSM parameter space plots (curvature_em_to_bsm.tex)"
	@echo "  make bsm-paper     - Build BSM paper PDF (kappaR_to_BSM/curvature_em_to_bsm.tex)"
	@echo "  make all-papers    - Build all three papers"
	@echo ""
	@echo "Release:"
	@echo "  make manifest      - Generate data_manifest.csv for results and figures"
	@echo "  make release-verify- Run release verification script"

test:
	pytest -v

quick-bench:
	python scripts/benchmark_solver.py --resolution 41 --runs 2 --materials rb87_bec --xi 100

bench:
	python scripts/benchmark_solver.py --resolution 61 --runs 3 --materials rb87_bec nb_cavity --xi 100

domain-sweep:
	python examples/domain_bc_sweep.py --resolution 41 --padding 1.2 1.5 2.0 2.5

lint:
	flake8 src/ examples/ tests/ --max-line-length=100 --ignore=E501,W503,D100,D101,D102,D103,D104,D105,D107

format:
	black src/ examples/ tests/
	isort src/ examples/ tests/

format-check:
	black --check --diff src/ examples/ tests/
	isort --check-only --diff src/ examples/ tests/

clean:
	find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null || true
	find . -type f -name "*.pyc" -delete
	find . -type f -name "*.pyo" -delete
	find . -type d -name ".pytest_cache" -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name "*.egg-info" -exec rm -rf {} + 2>/dev/null || true
	rm -f .coverage coverage.xml
	rm -f benchmark_results.json domain_padding_sweep.json

cache-info:
	@python -c "from src.utils.result_cache import get_cache; get_cache().info()"

cache-clean:
	@python -c "from src.utils.result_cache import get_cache; get_cache().clear(); print('✅ Cache cleared')"

analysis:
	@python scripts/run_analysis.py

optimize:
	@python scripts/optimize_geometry.py --resolution 41 --method Nelder-Mead

install-dev:
	pip install -e .
	pip install pytest pytest-asyncio black isort flake8 tox pre-commit scipy
	pre-commit install

env:
	conda env create -f environment.yml || true
	@echo "Activate with: conda activate cohgrav"

# Generate figures for main coherence-gravity coupling paper
figures:
	python scripts/generate_figures.py

# Build main coherence-gravity coupling paper PDF
paper:
	cd papers && pdflatex coherence_gravity_coupling.tex && bibtex coherence_gravity_coupling \
	  && pdflatex coherence_gravity_coupling.tex && pdflatex coherence_gravity_coupling.tex

# Build null_results PDF (curvature–EM coupling null-results paper)
null-results:
	cd papers && pdflatex null_results.tex && bibtex null_results \
	  || true; pdflatex null_results.tex && pdflatex null_results.tex

# Generate BSM parameter space plots (dark photon/axion)
bsm-figures:
	python scripts/generate_bsm_plots.py

# Build BSM parameter space paper PDF (κ_R → dark photon/axion mapping)
bsm-paper: bsm-figures
	cd papers/kappaR_to_BSM && pdflatex curvature_em_to_bsm.tex \
	  && (bibtex curvature_em_to_bsm || true) \
	  && pdflatex curvature_em_to_bsm.tex && pdflatex curvature_em_to_bsm.tex

# Build all three papers
all-papers: paper null-results bsm-paper
	@echo "✅ All papers built successfully"

# Generate consolidated analysis tables (CSV/Markdown/LaTeX)
report:
	python scripts/generate_report.py --all

manifest:
	python scripts/generate_manifest.py --output data_manifest.csv --roots results papers/figures

release-verify:
	bash scripts/verify_release.sh
