#!/usr/bin/env bash
set -euo pipefail

echo "[env] Using Python: $(python --version)"
echo "[env] PIP packages:"; pip freeze | sed 's/^/  /' | head -n 20 || true

echo "[tests] Running smoke tests..."
pytest -q || { echo "Tests failed"; exit 1; }

echo "[results] Running minimal production study (exists-ok)"
python production_study.py --materials YBCO Rb87 Nb --grid-size 5 --resolution 61 --jobs 4 --quick || true

echo "[results] Running convergence tests (exists-ok)"
python run_analysis.py convergence --xi 100 --Phi0 3.65e6 --resolutions 61 81 101 || true

echo "[figures] Generating manuscript figures..."
python generate_figures.py

echo "[manifest] Building data manifest..."
python scripts/generate_manifest.py --output data_manifest.csv --roots results papers/figures

echo "[done] Artifacts regenerated and manifest written."
