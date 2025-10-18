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
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent))

from examples.geometric_cavendish import run_geometric_cavendish
from examples.domain_bc_sweep import sweep_domain_padding


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
    parser_xi.add_argument('--verbose', action='store_true')
    
    # Materials sweep
    parser_mat = subparsers.add_parser('sweep-materials', help='Compare materials')
    parser_mat.add_argument('--xi', type=float, default=100.0)
    parser_mat.add_argument('--resolution', type=int, default=41)
    parser_mat.add_argument('--cache', action='store_true')
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
        save_results(results, "xi_sweep", 
                    f"Xi sweep: res={args.resolution}")
    
    elif args.command == 'sweep-materials':
        results = sweep_materials(
            materials=STANDARD_MATERIALS,
            xi=args.xi,
            resolution=args.resolution,
            cache=args.cache,
            verbose=args.verbose
        )
        save_results(results, "material_comparison",
                    f"Material comparison: xi={args.xi}, res={args.resolution}")
    
    print(f"\n✅ Analysis complete! Results: {RESULTS_DIR}")


if __name__ == '__main__':
    main()
