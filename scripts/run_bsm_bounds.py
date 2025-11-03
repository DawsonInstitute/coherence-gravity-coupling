#!/usr/bin/env python3
"""
Run BSM bounds sweeps from curvature–EM coupling κ_R.
Writes small CSV tables to results/bsm_bounds/.
"""
from __future__ import annotations
import os
import csv

import importlib.util
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
MOD_PATH = os.path.join(ROOT, 'src', 'analysis', 'bsm_bounds_from_kappa.py')
spec = importlib.util.spec_from_file_location('bsm_bounds_from_kappa', MOD_PATH)
mod = importlib.util.module_from_spec(spec)  # type: ignore[assignment]
spec.loader.exec_module(mod)  # type: ignore[union-attr]

DEFAULT_ENVIRONMENTS = mod.DEFAULT_ENVIRONMENTS
epsilon_sweep = mod.epsilon_sweep
axion_equiv_parametric = mod.axion_equiv_parametric

RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'bsm_bounds')
RESULTS_DIR = os.path.abspath(RESULTS_DIR)


def ensure_dir(p: str) -> None:
    os.makedirs(p, exist_ok=True)


def write_eps_csv(kappa_values):
    ensure_dir(RESULTS_DIR)
    out_path = os.path.join(RESULTS_DIR, 'epsilon_equiv.csv')
    tables = epsilon_sweep(kappa_values)
    # Flat CSV dump: env, C_eps, kappa_m2, epsilon
    rows = [("env", "C_eps", "kappa_m2", "epsilon")]
    for env_name, sub in tables.items():
        for C_key, row in sub.items():
            for kappa_str, eps in row.items():
                rows.append((env_name, C_key, kappa_str, f"{eps:.3e}"))
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerows(rows)
    print(f"Wrote {out_path}")


def write_axion_csv(kappa_values, Lambda_GeV=1e4):
    ensure_dir(RESULTS_DIR)
    out_path = os.path.join(RESULTS_DIR, 'g_axion_equiv.csv')
    rows = [("env", "C_a", "Lambda_GeV", "kappa_m2", "g_equiv_GeVinv")]
    for env_name, env in DEFAULT_ENVIRONMENTS.items():
        for C_a in (1.0, 1e-2):
            for k in kappa_values:
                g = axion_equiv_parametric(k, env.R_m2, C_a=C_a, Lambda_GeV=Lambda_GeV)
                rows.append((env_name, f"{C_a:.0e}", f"{Lambda_GeV:.2e}", f"{k:.3e}", f"{g:.3e}"))
    with open(out_path, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerows(rows)
    print(f"Wrote {out_path}")


def main():
    kappa_values = [1e-19, 1e-11, 1e0, 1e3]
    write_eps_csv(kappa_values)
    write_axion_csv(kappa_values)


if __name__ == "__main__":
    main()
