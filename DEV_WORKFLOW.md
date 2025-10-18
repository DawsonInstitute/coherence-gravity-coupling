# Local Development Workflow Guide

This repository uses **local development tools** instead of hosted CI (GitHub Actions). All checks run on your machine via `make`, `tox`, or `pre-commit`.

## Quick Start

### Install development dependencies
```bash
make install-dev
```

This installs:
- `pytest` for testing
- `black` and `isort` for formatting
- `flake8` for linting
- `tox` for automated testing
- `pre-commit` hooks

### Run tests
```bash
make test            # Full test suite (23 tests, ~90s)
make quick-bench     # Quick benchmark at 41³ (~30s)
make bench           # Full benchmark at 61³ (~2min)
make domain-sweep    # Domain padding sensitivity study
```

### Check code quality
```bash
make lint            # Run flake8 linter
make format-check    # Check formatting (no changes)
make format          # Auto-format code
```

### Clean up
```bash
make clean           # Remove __pycache__, .pyc, etc.
```

## Pre-commit Hooks

After running `make install-dev`, pre-commit hooks are installed automatically. They run on every `git commit` to:
- Format code with `black` and `isort`
- Lint with `flake8`
- Fix trailing whitespace and end-of-file issues
- Check for large files and merge conflicts

### Manual pre-commit run
```bash
pre-commit run --all-files
```

### Skip hooks (emergency only)
```bash
git commit --no-verify
```

## Tox Environments

Run multiple checks in isolated environments:

```bash
tox -e py313          # Run tests in Python 3.13
tox -e lint           # Run flake8
tox -e format-check   # Check formatting
tox -e format         # Auto-format
tox -e quick-bench    # Quick benchmark
tox                   # Run all environments
```

## Workflow Checklist

Before committing:
1. **Run tests**: `make test` or `pytest -v`
2. **Check formatting**: `make format-check` (or just commit, hooks will fix)
3. **Lint**: `make lint` (optional, hooks will check)
4. **Clean up**: `make clean` (removes temp files)

Pre-commit hooks will auto-format on commit. If hooks fail:
1. Review the changes
2. Stage the fixes: `git add -u`
3. Commit again: `git commit`

## Tools Configuration

### Black
- Line length: 100
- Target: Python 3.13
- Config: `pyproject.toml`

### isort
- Profile: black-compatible
- Line length: 100
- Config: `pyproject.toml`

### flake8
- Max line length: 100
- Ignore: E501 (line too long), W503 (line break before binary operator), D10x (docstring warnings)
- Config: `tox.ini`

### pytest
- Min version: 8.0
- Test paths: `tests/`
- Options: `-v --tb=short`
- Config: `pyproject.toml`

## Why No GitHub Actions?

This project uses **local-only tooling** to:
- Keep all checks on your machine (no cloud dependency)
- Run tests/benchmarks that may take minutes (not suitable for CI quotas)
- Maintain full control over test environments
- Support offline development

All quality gates run via `make` targets and `pre-commit` hooks.

## Troubleshooting

### Pre-commit hooks fail on first run
```bash
pre-commit run --all-files  # May need to run twice
git add -u                   # Stage auto-fixes
git commit                   # Try again
```

### Tox can't find Python 3.13
```bash
tox -e py313 --skip-missing-interpreters
```

### Want to skip formatting
```bash
git commit --no-verify       # Skip hooks (not recommended)
```

## Summary

- **Test**: `make test`
- **Format**: `make format` (or let pre-commit do it)
- **Lint**: `make lint`
- **Benchmark**: `make quick-bench` or `make bench`
- **All checks**: `tox`

Pre-commit hooks ensure code quality on every commit without manual intervention.
