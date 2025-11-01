# Latest Production Study Summary (auto)

This file summarizes the most recent entries under `results/production_study/`.

- Latest files detected:
  - production_study_20251018_194030.json (empty results; run interrupted)
  - production_study_20251018_175309.json (41³, all three materials, full grid + DE recorded)

Highlights (41³ quick study):
- YBCO (ξ=100, Φ₀=6.67e8):
  - Grid optimum: (0.0000, 0.0000, -0.0500) m, Δτ ≈ -1.18e-12 N·m
  - DE refinement: (-6.55e-03, 2.99e-03, -5.01e-02) m, Δτ ≈ -1.70e-12 N·m (×1.44)
- Rb-87 (ξ=100, Φ₀=3.65e6):
  - Similar geometry trends; |Δτ| at 41³ not reliable for publication; use 61³+ with volume averaging.
- Nb (ξ=100, Φ₀=3.65e6):
  - Similar geometry trends; |Δτ| requires 61³+.

Notes:
- A 61³ parallel run was started and then halted; partial progress suggests ~10–20 s per solve with diagonal preconditioner. Please re-run interactively to completion to produce a non-empty JSON and refreshed plots.

Recommended next command (interactive, visible):

```bash
python production_study.py --materials all --resolution 61 --grid-size 5 --jobs 4 --quick
```

This will generate a timestamped JSON and figures in `results/production_study/`. After it completes, we’ll update this summary and fold values into the manuscript.
