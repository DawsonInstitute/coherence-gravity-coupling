#!/usr/bin/env python3
"""
Broader sweep driver with optional Gaussian Process surrogate for adaptive sampling.

- Explores nontrivial parameter space: higher B, pulsed geometries (E != 0), curvature R, and precision δ.
- Trains a surrogate model (if scikit-learn available) to propose high-value points
  (e.g., regions of strongest expected constraint or highest uncertainty).
- Saves results as timestamped JSON and optional CSV under results/sweeps/.

Usage:
    python sweeps/run_broader_sweeps.py --samples 64 --output results/sweeps/broader.json \
        --B 0.5 10.0 --R 1e-30 1e-22 --precision 1e-10 1e-4 --E 0 0

Notes:
- If scikit-learn is not installed, falls back to random Latin hypercube sampling.
- GPU: This script is CPU-bound; the underlying calculations are light. Hooks can be added later.
"""
from __future__ import annotations

import argparse
import json
import math
import random
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple

try:
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import Matern, WhiteKernel, ConstantKernel as C
    SKLEARN_OK = True
except Exception:
    SKLEARN_OK = False

# Local import
from src.field_equations.curvature_coupling import (
    electromagnetic_invariants,
)


@dataclass
class SweepParams:
    B_min: float = 0.5
    B_max: float = 10.0
    R_min: float = 1e-30
    R_max: float = 1e-22
    E_min: float = 0.0
    E_max: float = 0.0
    precision_min: float = 1e-10
    precision_max: float = 1e-4


def kappa_limit(R: float, B: float, E: float, precision: float) -> Tuple[float, float]:
    """Compute κ_R exclusion limit given (R, B, E, δ) and return (kappa, F2)."""
    F2, _ = electromagnetic_invariants(E=E, B=B)
    F2 = max(F2, 1e-30)  # avoid division by zero if E ~ B/c
    kappa = precision / (abs(R) * abs(F2))
    return kappa, F2


def latin_hypercube(n: int, dims: int) -> List[List[float]]:
    """Simple Latin hypercube in [0,1]^dims."""
    grid = [list((i + random.random()) / n for i in range(n)) for _ in range(dims)]
    for g in grid:
        random.shuffle(g)
    return [list(pt) for pt in zip(*grid)]


def denorm(pt01: List[float], bounds: List[Tuple[float, float]]) -> List[float]:
    return [lo + p * (hi - lo) for p, (lo, hi) in zip(pt01, bounds)]


def main() -> None:
    ap = argparse.ArgumentParser(description="Broader sweeps with optional GP surrogate")
    ap.add_argument("--samples", type=int, default=64, help="Initial random samples")
    ap.add_argument("--iters", type=int, default=0, help="Adaptive iterations (requires scikit-learn)")
    ap.add_argument("--output", type=Path, default=Path("results/sweeps/broader.json"))
    ap.add_argument("--B", nargs=2, type=float, metavar=("B_MIN", "B_MAX"), default=[0.5, 10.0])
    ap.add_argument("--R", nargs=2, type=float, metavar=("R_MIN", "R_MAX"), default=[1e-30, 1e-22])
    ap.add_argument("--E", nargs=2, type=float, metavar=("E_MIN", "E_MAX"), default=[0.0, 0.0])
    ap.add_argument("--precision", nargs=2, type=float, metavar=("D_MIN", "D_MAX"), default=[1e-10, 1e-4])
    args = ap.parse_args()

    bounds = [
        (args.B[0], args.B[1]),
        (args.R[0], args.R[1]),
        (args.E[0], args.E[1]),
        (args.precision[0], args.precision[1]),
    ]

    # Initial design via Latin hypercube
    pts01 = latin_hypercube(args.samples, dims=4)
    pts = [denorm(p, bounds) for p in pts01]

    data = []
    for (B, R, E, D) in pts:
        kappa, F2 = kappa_limit(R=R, B=B, E=E, precision=D)
        data.append({"B": B, "R": R, "E": E, "precision": D, "kappa_limit": kappa, "F_squared": F2})

    # Optional adaptive refinement using GP (maximize predicted improvement = highest kappa_limit)
    if args.iters > 0 and SKLEARN_OK:
        import numpy as np
        X = np.array([[d["B"], d["R"], d["E"], d["precision"]] for d in data])
        y = np.array([d["kappa_limit"] for d in data])
        kernel = C(1.0, (1e-3, 1e3)) * Matern(length_scale=[1.0, 1.0, 1.0, 1.0], nu=2.5) + WhiteKernel(noise_level=1e-12)
        gp = GaussianProcessRegressor(kernel=kernel, normalize_y=True, n_restarts_optimizer=3)

        for _ in range(args.iters):
            gp.fit(X, y)
            # Propose candidates and pick the one with the largest predicted kappa
            cand01 = latin_hypercube(64, 4)
            cand = np.array([denorm(p, bounds) for p in cand01])
            mu, sigma = gp.predict(cand, return_std=True)
            idx = int(np.argmax(mu))  # greedily explore where constraint is weakest (largest κ)
            B, R, E, D = cand[idx]
            kappa, F2 = kappa_limit(R=R, B=B, E=E, precision=D)
            data.append({"B": float(B), "R": float(R), "E": float(E), "precision": float(D), "kappa_limit": float(kappa), "F_squared": float(F2)})
            # Update training set incrementally
            X = np.vstack([X, [B, R, E, D]])
            y = np.append(y, kappa)
    elif args.iters > 0 and not SKLEARN_OK:
        print("[warn] scikit-learn not available; skipping adaptive surrogate stage")

    # Save results
    out = {
        "timestamp": datetime.now().strftime("%Y%m%d_%H%M%S"),
        "bounds": {"B": args.B, "R": args.R, "E": args.E, "precision": args.precision},
        "count": len(data),
        "data": data,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Saved sweep: {args.output} ({len(data)} points)")


if __name__ == "__main__":
    main()
