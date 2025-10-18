#!/usr/bin/env python3
"""
Domain size and boundary condition sensitivity study.

Sweeps domain padding factor and BC type to determine optimal defaults
that keep Δτ variations < 5%.

Author: GitHub Copilot (Claude Sonnet 4.5)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from examples.geometric_cavendish import run_geometric_cavendish, CavendishGeometry
from src.analysis.phi_calibration import get_all_calibrations
import numpy as np
import json
from typing import Dict, List
import time


def compute_min_domain_size(geom: CavendishGeometry, padding: float = 1.5) -> float:
    """
    Compute minimum domain size to encompass all geometry.
    
    Args:
        geom: Cavendish geometry
        padding: Safety factor (default 1.5× bounding box)
    
    Returns:
        Domain size [m]
    """
    # Find bounding box of all masses
    positions = [
        geom.test1_pos,
        geom.test2_pos,
        geom.source1_pos,
        geom.source2_pos,
        np.array(geom.coherent_position)
    ]
    
    radii = [
        geom.R_test,
        geom.R_test,
        geom.R_source,
        geom.R_source,
        max(geom.coherent_volume) / 2
    ]
    
    max_extent = 0.0
    for pos, r in zip(positions, radii):
        extent = np.max(np.abs(pos)) + r
        max_extent = max(max_extent, extent)
    
    return 2 * max_extent * padding


def sweep_domain_padding(
    padding_factors: List[float] = [1.2, 1.5, 2.0, 2.5],
    test_config: Dict = None,
    resolution: int = 41,
    verbose: bool = True
) -> Dict:
    """
    Sweep domain padding to measure Δτ sensitivity.
    
    Args:
        padding_factors: List of padding multipliers
        test_config: Test configuration (xi, Phi0, geom_params)
        resolution: Grid resolution
        verbose: Print progress
    
    Returns:
        Dict with results and analysis
    """
    if test_config is None:
        # Default: YBCO at ξ=100 with offset position
        calibrations = get_all_calibrations()
        ybco = calibrations['ybco_cuprate']
        test_config = {
            'xi': 100.0,
            'Phi0': ybco.Phi0,
            'geom_params': {'coherent_position': (0.0, 0.0, -0.08)}
        }
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"DOMAIN PADDING SENSITIVITY STUDY")
        print(f"{'='*70}")
        print(f"Configuration: xi={test_config['xi']}, Phi0={test_config['Phi0']:.2e}")
        print(f"Resolution: {resolution}³")
    
    # Compute base domain size
    geom = CavendishGeometry(**test_config.get('geom_params', {}))
    base_domain = compute_min_domain_size(geom, padding=1.0)
    
    if verbose:
        print(f"Minimum enclosing size: {base_domain:.3f} m")
        print(f"\nSweeping padding factors: {padding_factors}")
    
    results = []
    
    for padding in padding_factors:
        domain_size = base_domain * padding
        
        if verbose:
            print(f"\n--- Padding {padding:.1f}× (domain = {domain_size:.3f} m) ---")
        
        t0 = time.time()
        result = run_geometric_cavendish(
            xi=test_config['xi'],
            Phi0=test_config['Phi0'],
            geom_params=test_config.get('geom_params'),
            grid_resolution=resolution,
            domain_size=domain_size,
            verbose=verbose
        )
        t_elapsed = time.time() - t0
        
        results.append({
            'padding': padding,
            'domain_size': domain_size,
            'tau_coherent': result['tau_coherent'],
            'tau_newtonian': result['tau_newtonian'],
            'delta_tau': result['delta_tau'],
            'delta_G_over_G': result['delta_G_over_G'],
            'solve_time': t_elapsed
        })
        
        if verbose:
            print(f"   Δτ = {result['delta_tau']:.6e} N·m")
            print(f"   ΔG/G = {result['delta_G_over_G']:.6e}")
    
    # Analysis
    delta_taus = np.array([r['delta_tau'] for r in results])
    mean_delta_tau = np.mean(delta_taus)
    std_delta_tau = np.std(delta_taus)
    max_variation = np.max(np.abs(delta_taus - mean_delta_tau)) / np.abs(mean_delta_tau)
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"ANALYSIS")
        print(f"{'='*70}")
        print(f"Mean Δτ: {mean_delta_tau:.6e} N·m")
        print(f"Std dev: {std_delta_tau:.6e} N·m")
        print(f"Max variation: {max_variation*100:.2f}%")
        
        if max_variation < 0.05:
            print(f"✅ Variation < 5% - all padding factors acceptable")
            print(f"   Recommended: padding = {min(padding_factors):.1f}× (smallest stable)")
        else:
            print(f"⚠️  Variation ≥ 5% - domain boundary effects present")
            stable = [r['padding'] for r in results 
                     if abs(r['delta_tau'] - mean_delta_tau) / abs(mean_delta_tau) < 0.05]
            if stable:
                print(f"   Recommended: padding ≥ {min(stable):.1f}×")
            else:
                print(f"   Recommended: padding ≥ {max(padding_factors):.1f}× (need larger)")
    
    return {
        'test_config': test_config,
        'resolution': resolution,
        'base_domain': base_domain,
        'results': results,
        'mean_delta_tau': float(mean_delta_tau),
        'std_delta_tau': float(std_delta_tau),
        'max_variation': float(max_variation),
        'recommended_padding': float(min(padding_factors)) if max_variation < 0.05 else None
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Domain padding sensitivity study')
    parser.add_argument('--resolution', type=int, default=41,
                        help='Grid resolution (default: 41)')
    parser.add_argument('--padding', nargs='+', type=float,
                        default=[1.2, 1.5, 2.0, 2.5],
                        help='Padding factors to test (default: 1.2 1.5 2.0 2.5)')
    parser.add_argument('--xi', type=float, default=100.0,
                        help='Non-minimal coupling (default: 100.0)')
    parser.add_argument('--material', type=str, default='ybco_cuprate',
                        help='Material from phi_calibration (default: ybco_cuprate)')
    parser.add_argument('--output', type=str, default='domain_padding_sweep.json',
                        help='Output JSON file (default: domain_padding_sweep.json)')
    
    args = parser.parse_args()
    
    # Load material
    calibrations = get_all_calibrations()
    if args.material not in calibrations:
        print(f"Error: Material '{args.material}' not found")
        print(f"Available: {list(calibrations.keys())}")
        return 1
    
    cal = calibrations[args.material]
    
    test_config = {
        'xi': args.xi,
        'Phi0': cal.Phi0,
        'geom_params': {'coherent_position': (0.0, 0.0, -0.08)}  # Offset config
    }
    
    # Run sweep
    results = sweep_domain_padding(
        padding_factors=args.padding,
        test_config=test_config,
        resolution=args.resolution,
        verbose=True
    )
    
    # Save results
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to {args.output}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
