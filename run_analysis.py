#!/usr/bin/env python3
"""
Master Analysis Script for Coherence-Gravity Coupling Framework

Centralized script for running different types of analysis:
- Parameter sweeps (xi, Phi0, geometry)
- Domain sensitivity studies
- Convergence analysis
- Material comparisons

All results are cached and saved with timestamps for reproducibility.

Usage:
    python run_analysis.py --help
    python run_analysis.py sweep-xi --cache
    python run_analysis.py sweep-materials --resolution 61
    
Author: GitHub Copilot (Claude Sonnet 4.5)
Date: October 2025
License: MIT
"""

import argparse
import numpy as np
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from examples.geometric_cavendish import run_geometric_cavendish
from src.utils.result_cache import ResultCache
from src.visualization.plot_utils import (
    plot_parameter_sweep, plot_material_comparison, save_figure, plot_exclusion_limits
)
from src.field_equations.curvature_coupling import (
    CurvatureCouplingCalculator, compute_exclusion_limits
)


# Output directory
RESULTS_DIR = Path("results/analysis")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def save_results(data: Dict, name: str, description: str = "") -> Path:
    """Save analysis results with timestamp and metadata."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{name}_{timestamp}.json"
    filepath = RESULTS_DIR / filename
    
    output = {
        "timestamp": timestamp,
        "description": description,
        "data": data
    }
    
    with open(filepath, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n💾 Results saved: {filepath}")
    return filepath


def sweep_xi(
    xi_values: List[float],
    Phi0: float = 1e8,
    resolution: int = 41,
    cache: bool = True,
    verbose: bool = False
) -> Dict:
    """Sweep non-minimal coupling strength xi."""
    print(f"\n{'='*70}")
    print("XI PARAMETER SWEEP")
    print(f"{'='*70}")
    print(f"xi values: {xi_values}")
    print(f"Phi0: {Phi0:.2e} m⁻¹")
    print(f"Resolution: {resolution}³")
    print(f"Caching: {'enabled' if cache else 'disabled'}\n")
    
    results = {}
    
    for i, xi in enumerate(xi_values):
        print(f"[{i+1}/{len(xi_values)}] Running xi = {xi}...")
        
        t0 = time.time()
        result = run_geometric_cavendish(
            xi=xi,
            Phi0=Phi0,
            grid_resolution=resolution,
            cache=cache,
            verbose=verbose
        )
        t_elapsed = time.time() - t0
        
        results[f"xi_{xi}"] = {
            "xi": xi,
            "delta_tau": result["delta_tau"],
            "delta_G_over_G": result["delta_G_over_G"],
            "tau_coherent": result["tau_coherent"],
            "tau_newtonian": result["tau_newtonian"],
            "solve_time": result["solve_time_coherent"] + result["solve_time_newtonian"],
            "elapsed_time": t_elapsed
        }
        
        print(f"   Δτ = {result['delta_tau']:.3e} N·m")
        print(f"   Time: {t_elapsed:.2f} s\n")
    
    return results


def sweep_materials(
    materials: List[Dict],
    xi: float = 100.0,
    resolution: int = 41,
    cache: bool = True,
    verbose: bool = False
) -> Dict:
    """Compare different coherent materials."""
    print(f"\n{'='*70}")
    print("MATERIAL COMPARISON")
    print(f"{'='*70}")
    print(f"Materials: {[m['name'] for m in materials]}")
    print(f"xi: {xi}")
    print(f"Resolution: {resolution}³\n")
    
    results = {}
    
    for i, material in enumerate(materials):
        name = material['name']
        Phi0 = material['Phi0']
        geom_params = material.get('geom_params', {})
        
        print(f"[{i+1}/{len(materials)}] Running {name} (Phi0 = {Phi0:.2e} m⁻¹)...")
        
        t0 = time.time()
        result = run_geometric_cavendish(
            xi=xi,
            Phi0=Phi0,
            geom_params=geom_params,
            grid_resolution=resolution,
            cache=cache,
            verbose=verbose
        )
        t_elapsed = time.time() - t0
        
        results[name] = {
            "name": name,
            "Phi0": Phi0,
            "delta_tau": result["delta_tau"],
            "delta_G_over_G": result["delta_G_over_G"],
            "elapsed_time": t_elapsed
        }
        
        print(f"   Δτ = {result['delta_tau']:.3e} N·m")
        print(f"   Time: {t_elapsed:.2f} s\n")
    
    return results


def sweep_curvature_limits(
    B_values: List[float],
    R_value: float,
    precision: float = 1e-6,
    E_value: float = 0.0,
    verbose: bool = False
) -> Dict:
    """Sweep exclusion limits vs magnetic field strength for fixed Ricci scale.

    Args:
        B_values: List of magnetic field magnitudes [T]
        R_value: Ricci scalar [m^-2]
        precision: Experimental relative precision (e.g., 1e-6)
        E_value: Electric field magnitude [V/m]
        verbose: Print details

    Returns:
        Dict keyed by B value labels with kappa limits and inputs
    """
    print(f"\n{'='*70}")
    print("CURVATURE COUPLING EXCLUSION LIMITS (κ_R)")
    print(f"{'='*70}")
    print(f"B values [T]: {B_values}")
    print(f"Ricci scalar R: {R_value:.2e} m⁻²")
    print(f"Precision: {precision:.2e}\n")

    results: Dict[str, Dict] = {}
    calc = CurvatureCouplingCalculator(
        # Parameters not used for limits directly; using calculator for invariants
        params=None  # type: ignore
    )

    for i, B in enumerate(B_values):
        print(f"[{i+1}/{len(B_values)}] Evaluating B = {B} T...")
        E_vec = np.array([E_value, 0.0, 0.0])
        B_vec = np.array([0.0, B, 0.0])

        invariants = calc.electromagnetic_invariants(E_vec, B_vec)
        F_sq = invariants['F_squared']

        limits = compute_exclusion_limits(
            experimental_precision=precision,
            ricci_scale=R_value,
            field_strength=F_sq
        )

        results[f"B_{B}"] = {
            'B': B,
            'E': E_value,
            'R': R_value,
            'precision': precision,
            'kappa_limit': limits['kappa_ricci_em_limit'],
            'F_squared': F_sq
        }

        print(f"   κ_R limit < {limits['kappa_ricci_em_limit']:.2e} m²\n")

    return results


# Standard material configurations
STANDARD_MATERIALS = [
    {"name": "rb87_bec", "Phi0": 3.65e6, "geom_params": {"coherent_position": [0.0, 0.0, -0.08]}},
    {"name": "nb_cavity", "Phi0": 3.65e6, "geom_params": {"coherent_position": [0.0, 0.0, -0.08]}},
    {"name": "ybco_cuprate", "Phi0": 6.67e8, "geom_params": {"coherent_position": [0.0, 0.0, -0.08]}}
]


def main():
    """Main entry point for analysis script."""
    parser = argparse.ArgumentParser(
        description="Master analysis script for coherence-gravity coupling"
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Analysis type')
    
    # Xi sweep
    parser_xi = subparsers.add_parser('sweep-xi', help='Sweep xi parameter')
    parser_xi.add_argument('--xi', nargs='+', type=float, 
                          default=[10, 50, 100, 200, 500])
    parser_xi.add_argument('--Phi0', type=float, default=1e8)
    parser_xi.add_argument('--resolution', type=int, default=41)
    parser_xi.add_argument('--cache', action='store_true')
    parser_xi.add_argument('--plot', action='store_true', help='Generate plots')
    parser_xi.add_argument('--verbose', action='store_true')
    
    # Materials sweep
    parser_mat = subparsers.add_parser('sweep-materials', help='Compare materials')
    # Curvature coupling exclusion limits
    parser_curv = subparsers.add_parser('sweep-curvature', help='Exclusion limits vs magnetic field')
    parser_curv.add_argument('--B', nargs='+', type=float, default=[0.5, 1.0, 3.0, 10.0], help='Magnetic field strengths [T]')
    parser_curv.add_argument('--R', type=float, default=1e-26, help='Ricci scalar [m^-2]')
    parser_curv.add_argument('--precision', type=float, default=1e-6, help='Experimental relative precision')
    parser_curv.add_argument('--E', type=float, default=0.0, help='Electric field magnitude [V/m]')
    parser_curv.add_argument('--plot', action='store_true', help='Generate plots')
    parser_curv.add_argument('--verbose', action='store_true')
    parser_mat.add_argument('--xi', type=float, default=100.0)
    parser_mat.add_argument('--resolution', type=int, default=41)
    parser_mat.add_argument('--cache', action='store_true')
    parser_mat.add_argument('--plot', action='store_true', help='Generate plots')
    parser_mat.add_argument('--verbose', action='store_true')
    
    args = parser.parse_args()
    
    if args.command is None:
        parser.print_help()
        return
    
    # Run appropriate analysis
    if args.command == 'sweep-xi':
        results = sweep_xi(
            xi_values=args.xi,
            Phi0=args.Phi0,
            resolution=args.resolution,
            cache=args.cache,
            verbose=args.verbose
        )
        filepath = save_results(results, "xi_sweep", 
                    f"Xi sweep: res={args.resolution}")
        
        if args.plot:
            print("\n📊 Generating plots...")
            # Prepare data for plotting
            plot_data = []
            for key, val in results.items():
                if key.startswith('xi_'):
                    plot_data.append({
                        'param_value': val['xi'],
                        'delta_tau': val['delta_tau'],
                        'compute_time': val['elapsed_time'],
                        'cache_hit': val.get('cache_hit', False)
                    })
            output_path = filepath.parent / f"{filepath.stem}_plot"
            plot_parameter_sweep(plot_data, 'xi', output_path, title='Coherence Coupling Sweep (ξ)')
    
    elif args.command == 'sweep-materials':
        results = sweep_materials(
            materials=STANDARD_MATERIALS,
            xi=args.xi,
            resolution=args.resolution,
            cache=args.cache,
            verbose=args.verbose
        )
        filepath = save_results(results, "material_comparison",
                    f"Material comparison: xi={args.xi}, res={args.resolution}")
        
        if args.plot:
            print("\n📊 Generating plots...")
            # Prepare data for plotting
            plot_data = []
            for key, val in results.items():
                plot_data.append({
                    'material': val['name'],
                    'delta_tau': val['delta_tau'],
                    'config': val
                })
            output_path = filepath.parent / f"{filepath.stem}_plot"
            plot_material_comparison(plot_data, output_path)

    elif args.command == 'sweep-curvature':
        results = sweep_curvature_limits(
            B_values=args.B,
            R_value=args.R,
            precision=args.precision,
            E_value=args.E,
            verbose=args.verbose
        )
        filepath = save_results(results, "curvature_limits",
                    f"Curvature limits: R={args.R:.2e}, precision={args.precision:.1e}")

        if args.plot:
            print("\n📊 Generating plots...")
            plot_data = []
            for key, val in results.items():
                plot_data.append({
                    'param_value': val['B'],
                    'kappa_limit': val['kappa_limit']
                })
            output_path = filepath.parent / f"{filepath.stem}_plot"
            plot_exclusion_limits(plot_data, output_path)
    
    print(f"\n✅ Analysis complete! Results: {RESULTS_DIR}")


if __name__ == '__main__':
    main()
