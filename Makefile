.PHONY: help test quick-bench bench lint format format-check clean install-dev domain-sweep

help:
	@echo "coherence-gravity-coupling development targets"
	@echo ""
	@echo "Testing:"
	@echo "  make test          - Run full test suite (23 tests, ~90s)"
	@echo "  make quick-bench   - Quick benchmark at 41³ (2 runs, ~30s)"
	@echo "  make bench         - Full benchmark at 61³ (3 runs, ~2min)"
	@echo "  make domain-sweep  - Domain padding sensitivity study"
	@echo ""
	@echo "Code quality:"
	@echo "  make lint          - Run flake8 linter"
	@echo "  make format        - Auto-format code (black + isort)"
	@echo "  make format-check  - Check formatting without changes"
	@echo ""
	@echo "Setup:"
	@echo "  make install-dev   - Install development dependencies"
	@echo "  make clean         - Remove generated files"

test:
	pytest -v

quick-bench:
	python benchmark_solver.py --resolution 41 --runs 2 --materials rb87_bec --xi 100

bench:
	python benchmark_solver.py --resolution 61 --runs 3 --materials rb87_bec nb_cavity --xi 100

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

install-dev:
	pip install -e .
	pip install pytest pytest-asyncio black isort flake8 tox pre-commit
	pre-commit install
