#!/usr/bin/env python3
"""
Benchmark solver performance improvements.

Tests different solver methods and preconditioners to measure speedup
relative to baseline (CG with no preconditioner).

Author: GitHub Copilot (Claude Sonnet 4.5)
"""

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from examples.geometric_cavendish import run_geometric_cavendish
from src.analysis.phi_calibration import get_all_calibrations
import time
import json
from typing import Dict, List
import numpy as np


def benchmark_configuration(
    xi: float,
    Phi0: float,
    resolution: int,
    solver_method: str,
    preconditioner: str,
    n_runs: int = 3
) -> Dict:
    """
    Benchmark a single configuration with multiple runs.
    
    Returns:
        dict with 'mean_time', 'std_time', 'residuals', 'torque'
    """
    times = []
    residuals = []
    torques = []
    
    print(f"\n   Testing: {solver_method} + {preconditioner} (n={n_runs})")
    
    for i in range(n_runs):
        result = run_geometric_cavendish(
            xi=xi,
            Phi0=Phi0,
            grid_resolution=resolution,
            verbose=False,
            solver_method=solver_method,
            preconditioner=preconditioner
        )
        
        sys.path.insert(0, str(Path(__file__).parent.parent))
        residuals.append(result.get('residual_coherent', 0.0))
        torques.append(result['tau_coherent'])
        
        print(f"      Run {i+1}: {times[-1]:.2f} s, residual: {residuals[-1]:.2e}")
    
    return {
        'mean_time': float(np.mean(times)),
        'std_time': float(np.std(times)),
        'min_time': float(np.min(times)),
        'max_time': float(np.max(times)),
        'mean_residual': float(np.mean(residuals)),
        'mean_torque': float(np.mean(torques)),
        'solver_method': solver_method,
        'preconditioner': preconditioner
    }


def run_benchmark(
    resolution: int = 81,
    n_runs: int = 3,
    test_materials: List[str] = ['rb87_bec'],
    xi: float = 100.0
) -> Dict:
    """
    Run comprehensive benchmark comparing solver configurations.
    
    Args:
        resolution: Grid resolution (e.g., 81 for 81³)
        n_runs: Number of runs per configuration
        test_materials: Materials to test from phi_calibration
        xi: Non-minimal coupling strength (default: 100.0)
    
    Returns:
        Dictionary with benchmark results and analysis
    """
    print(f"\n{'='*70}")
    print(f"SOLVER PERFORMANCE BENCHMARK")
    print(f"{'='*70}")
    print(f"Resolution: {resolution}³ = {resolution**3:,} points")
    print(f"Runs per config: {n_runs}")
    print(f"ξ: {xi}")
    print(f"Materials: {', '.join(test_materials)}")
    
    calibrations = get_all_calibrations()
    
    # Test configurations
    configs = [
        ('cg', 'none'),         # Baseline
        ('cg', 'diagonal'),     # Fast preconditioner
        ('bicgstab', 'none'),   # Alternative solver
        ('bicgstab', 'diagonal'),
    ]
    
    # Try AMG if available
    try:
        import pyamg
        configs.append(('cg', 'amg'))
        print("   PyAMG available, will test AMG preconditioner")
    except ImportError:
        print("   PyAMG not available, skipping AMG tests")
    
    results = {}
    
    for material in test_materials:
        if material not in calibrations:
            print(f"\nWarning: Material '{material}' not found in calibrations, skipping")
            continue
        
        cal = calibrations[material]
        Phi0 = cal.Phi0
        
        print(f"\n{'-'*70}")
        print(f"Material: {material}")
        print(f"   ξ = {xi}")
        print(f"   Φ₀ = {Phi0:.2e} m⁻¹")
        print(f"{'-'*70}")
        
        material_results = {}
        
        for solver_method, preconditioner in configs:
            config_name = f"{solver_method}+{preconditioner}"
            
            try:
                bench_result = benchmark_configuration(
                    xi=xi,
                    Phi0=Phi0,
                    resolution=resolution,
                    solver_method=solver_method,
                    preconditioner=preconditioner,
                    n_runs=n_runs
                )
                material_results[config_name] = bench_result
                
            except Exception as e:
                print(f"      ERROR: {e}")
                material_results[config_name] = {
                    'error': str(e),
                    'solver_method': solver_method,
                    'preconditioner': preconditioner
                }
        
        results[material] = material_results
    
    # Compute speedups relative to baseline
    print(f"\n{'='*70}")
    print(f"PERFORMANCE SUMMARY")
    print(f"{'='*70}")
    
    for material, mat_results in results.items():
        print(f"\n{material}:")
        
        baseline_time = mat_results.get('cg+none', {}).get('mean_time', None)
        
        if baseline_time is None:
            print("   No baseline result available")
            continue
        
        print(f"   Baseline (cg+none): {baseline_time:.2f} s")
        print(f"\n   {'Configuration':<25} {'Time (s)':<12} {'Speedup':<10} {'Residual':<12}")
        print(f"   {'-'*70}")
        
        for config_name, result in mat_results.items():
            if 'error' in result:
                print(f"   {config_name:<25} ERROR: {result['error']}")
                continue
            
            mean_time = result['mean_time']
            speedup = baseline_time / mean_time
            residual = result['mean_residual']
            
            print(f"   {config_name:<25} {mean_time:>10.2f}  {speedup:>8.2f}x  {residual:>10.2e}")
            
            # Store speedup in result
            result['speedup'] = float(speedup)
    
    # Overall statistics
    print(f"\n{'='*70}")
    print(f"OPTIMIZATION RECOMMENDATIONS")
    print(f"{'='*70}")
    
    best_configs = {}
    for material, mat_results in results.items():
        if 'cg+none' not in mat_results:
            continue
        
        valid_results = [(k, v) for k, v in mat_results.items() if 'error' not in v]
        if not valid_results:
            continue
        
        # Find fastest configuration
        best = min(valid_results, key=lambda x: x[1]['mean_time'])
        best_configs[material] = {
            'config': best[0],
            'time': best[1]['mean_time'],
            'speedup': best[1].get('speedup', 1.0)
        }
        
        print(f"\n{material}:")
        print(f"   Best: {best[0]} ({best[1]['mean_time']:.2f} s, {best[1].get('speedup', 1.0):.2f}x)")
    
    # Check if target met
    target_speedup = 2.0
    print(f"\n{'='*70}")
    if any(v['speedup'] >= target_speedup for v in best_configs.values()):
        print(f"✅ TARGET MET: ≥{target_speedup}× speedup achieved")
    else:
        max_speedup = max(v['speedup'] for v in best_configs.values()) if best_configs else 0
        print(f"⚠️  TARGET NOT MET: Best speedup {max_speedup:.2f}× < {target_speedup}×")
    print(f"{'='*70}")
    
    return {
        'resolution': resolution,
        'n_runs': n_runs,
        'materials': test_materials,
        'results': results,
        'best_configs': best_configs,
        'target_speedup': target_speedup
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Benchmark solver performance')
    parser.add_argument('--resolution', type=int, default=81,
                        help='Grid resolution (default: 81 for 81³)')
    parser.add_argument('--runs', type=int, default=3,
                        help='Number of runs per configuration (default: 3)')
    parser.add_argument('--materials', nargs='+', default=['rb87_bec'],
                        help='Materials to test (default: rb87_bec)')
    parser.add_argument('--xi', type=float, default=100.0,
                        help='Non-minimal coupling (default: 100.0)')
    parser.add_argument('--output', type=str, default='benchmark_results.json',
                        help='Output JSON file (default: benchmark_results.json)')
    
    args = parser.parse_args()
    
    # Run benchmark
    results = run_benchmark(
        resolution=args.resolution,
        n_runs=args.runs,
        test_materials=args.materials,
        xi=args.xi
    )
    
    # Save results
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Results saved to {args.output}")


if __name__ == '__main__':
    main()
