# Makefile for Coherence-Gravity Coupling Framework

.PHONY: help install test analysis geometric-cavendish benchmark cli-test clean

help:
	@echo "Coherence-Gravity Coupling Framework"
	@echo ""
	@echo "Available targets:"
	@echo "  install              Install dependencies and package"
	@echo "  test                 Run all tests"
	@echo "  analysis             Run parameter space analysis"
	@echo "  geometric-cavendish  Run geometric Cavendish simulation"
	@echo "  benchmark            Run solver acceleration benchmark"
	@echo "  cli-test             Test CLI entry point"
	@echo "  clean                Remove generated files"

install:
	pip install -r requirements.txt
	pip install -e .

test:
	pytest tests/ -v

analysis:
	python run_analysis.py

geometric-cavendish:
	python examples/geometric_cavendish.py

benchmark:
	python benchmarks/solver_acceleration.py

cli-test:
	python -m cgc --help

clean:
	rm -rf results/*.png results/*.json __pycache__ .pytest_cache
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -type d -exec rm -rf {} + 2>/dev/null || true
