#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR"

echo "[check] LICENSE present"; test -f LICENSE
echo "[check] CITATION.cff present"; test -f CITATION.cff

echo "[test] Running smoke tests"
pytest -q

echo "[figures] Verifying expected figure files"
for f in papers/figures/convergence_analysis.pdf papers/figures/material_comparison.pdf papers/figures/landscape_YBCO_z_slice.pdf; do
  if [ ! -f "$f" ]; then
    echo "[warn] Missing $f, generating figures...";
    python generate_figures.py
    break
  fi
done

echo "[manifest] Generating data_manifest.csv"
python scripts/generate_manifest.py --output data_manifest.csv --roots results papers/figures

COMMIT_SHA=$(git rev-parse HEAD || echo "unknown")
DATE_UTC=$(date -u +%Y-%m-%dT%H:%M:%SZ)
PY=$(python -c 'import sys; print(sys.version.split()[0])')

echo "[release] Writing release_manifest.json"
python - <<'PY'
import json, os
manifest = {
  "commit": os.environ.get("COMMIT_SHA", "unknown"),
  "date_utc": os.environ.get("DATE_UTC", "unknown"),
  "python": os.environ.get("PY", "unknown"),
  "env": {
    "pip_freeze": os.popen('pip freeze').read().strip().splitlines(),
  },
  "artifacts": [
    "data_manifest.csv",
    "papers/coherence_gravity_coupling.pdf",
    "papers/figures/convergence_analysis.pdf",
    "papers/figures/material_comparison.pdf",
    "papers/figures/landscape_YBCO_z_slice.pdf",
  ],
}
with open('release_manifest.json','w') as f:
  json.dump(manifest, f, indent=2)
print("release_manifest.json written")
PY

echo "[ok] verify_release completed"
