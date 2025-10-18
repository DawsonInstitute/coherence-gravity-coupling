"""
Solver Acceleration Benchmark

Tests preconditioner performance for 3D Poisson solver with high-contrast G_eff.

Compares:
- No preconditioner (baseline)
- Incomplete LU (ILU)
- Algebraic Multigrid (AMG) via PyAMG

Author: GitHub Copilot (Claude Sonnet 4.5)
License: MIT
"""

import numpy as np
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.solvers.poisson_3d import Poisson3DSolver, Grid3D, spherical_mass_coherent_shell

# Try to import PyAMG
try:
    import pyamg
    HAS_PYAMG = True
except ImportError:
    HAS_PYAMG = False
    print("PyAMG not available - install with: pip install pyamg")


def benchmark_preconditioners(
    grid_sizes: list = [21, 31, 41, 51],
    xi: float = 100.0,
    Phi0: float = 1e7
):
    """
    Benchmark solver performance with different preconditioners.
    """
    print("\n" + "="*70)
    print("SOLVER ACCELERATION BENCHMARK")
    print("="*70)
    print(f"\nConfiguration:")
    print(f"   ξ = {xi}")
    print(f"   Φ₀ = {Phi0:.2e} m⁻¹")
    print(f"   Test case: Spherical mass with coherent shell")
    
    results = []
    
    for N in grid_sizes:
        print(f"\n{'='*70}")
        print(f"Grid Size: {N}³ = {N**3:,} DOF")
        print(f"{'='*70}")
        
        grid = Grid3D(nx=N, ny=N, nz=N, Lx=2.0, Ly=2.0, Lz=2.0)
        
        rho_func, Phi_func, meta = spherical_mass_coherent_shell(
            M=1.0,
            R_mass=0.05,
            R_shell_inner=0.15,
            R_shell_outer=0.40,
            Phi0=Phi0
        )
        
        solver = Poisson3DSolver(grid, xi=xi)
        
        # Test each preconditioner
        for precond in ['none', 'ilu', 'amg']:
            if precond == 'amg' and not HAS_PYAMG:
                print(f"\n--- {precond.upper()}: SKIPPED (not available) ---")
                continue
            
            print(f"\n--- Preconditioner: {precond.upper()} ---")
            
            try:
                t_start = time.time()
                solution = solver.solve(
                    rho_func, Phi_func,
                    method='cg',
                    tol=1e-8,
                    preconditioner=precond
                )
                t_total = time.time() - t_start
                
                info = solution.solver_info
                
                result = {
                    'grid_size': N,
                    'dof': N**3,
                    'preconditioner': precond,
                    'converged': info['converged'],
                    'solve_time': info['solve_time'],
                    'precond_time': info.get('precond_time', 0.0),
                    'total_time': t_total,
                    'residual': info['residual'],
                    'iterations': info.get('iterations', 'unknown')
                }
                
                results.append(result)
                
                print(f"   Total time: {t_total:.2f} s")
                print(f"   Residual: {info['residual']:.3e}")
                
            except Exception as e:
                print(f"   ❌ FAILED: {e}")
                results.append({
                    'grid_size': N,
                    'dof': N**3,
                    'preconditioner': precond,
                    'converged': False,
                    'error': str(e)
                })
    
    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}\n")
    
    print(f"{'Grid':<10} {'DOF':<10} {'Precond':<10} {'Time [s]':<12} {'Speedup':<10}")
    print("-"*70)
    
    # Group by grid size
    for N in grid_sizes:
        grid_results = [r for r in results if r['grid_size'] == N and 'total_time' in r]
        if not grid_results:
            continue
        
        # Find baseline (none)
        baseline = next((r for r in grid_results if r['preconditioner'] == 'none'), None)
        baseline_time = baseline['total_time'] if baseline else 1.0
        
        for r in grid_results:
            speedup = baseline_time / r['total_time']
            print(f"{N}³{'':<6} {r['dof']:<10,} {r['preconditioner']:<10} "
                  f"{r['total_time']:<12.2f} {speedup:<10.2f}×")
    
    print(f"\n{'='*70}\n")
    
    # Save results
    import json
    output_dir = Path('results')
    output_dir.mkdir(exist_ok=True)
    
    with open(output_dir / 'solver_benchmark.json', 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"✅ Results saved to results/solver_benchmark.json")
    
    return results


def plot_benchmark_results(results: list):
    """Plot benchmark results."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib not available, skipping plot")
        return
    
    # Extract data
    grid_sizes = sorted(set(r['grid_size'] for r in results if 'total_time' in r))
    preconditioners = ['none', 'ilu', 'amg']
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Solve time vs grid size
    for precond in preconditioners:
        sizes = []
        times = []
        for N in grid_sizes:
            r = next((x for x in results if x['grid_size'] == N and 
                     x['preconditioner'] == precond and 'total_time' in x), None)
            if r:
                sizes.append(N**3)
                times.append(r['total_time'])
        
        if sizes:
            ax1.loglog(sizes, times, 'o-', label=precond.upper(), linewidth=2, markersize=8)
    
    ax1.set_xlabel('Grid Size (DOF)', fontsize=12)
    ax1.set_ylabel('Total Time [s]', fontsize=12)
    ax1.set_title('Solver Performance vs Grid Size', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Speedup vs grid size
    for precond in preconditioners:
        if precond == 'none':
            continue  # Skip baseline
        
        sizes = []
        speedups = []
        for N in grid_sizes:
            r_precond = next((x for x in results if x['grid_size'] == N and 
                             x['preconditioner'] == precond and 'total_time' in x), None)
            r_none = next((x for x in results if x['grid_size'] == N and 
                          x['preconditioner'] == 'none' and 'total_time' in x), None)
            
            if r_precond and r_none:
                sizes.append(N**3)
                speedup = r_none['total_time'] / r_precond['total_time']
                speedups.append(speedup)
        
        if sizes:
            ax2.semilogx(sizes, speedups, 'o-', label=precond.upper(), linewidth=2, markersize=8)
    
    ax2.axhline(1.0, color='k', linestyle='--', alpha=0.3, label='Baseline')
    ax2.axhline(5.0, color='r', linestyle=':', alpha=0.5, label='5× target')
    ax2.set_xlabel('Grid Size (DOF)', fontsize=12)
    ax2.set_ylabel('Speedup vs No Preconditioner', fontsize=12)
    ax2.set_title('Preconditioner Speedup', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    output_dir = Path('results')
    plt.savefig(output_dir / 'solver_benchmark.png', dpi=200)
    print(f"   Saved: results/solver_benchmark.png")
    plt.close()


if __name__ == '__main__':
    # Run benchmark
    results = benchmark_preconditioners(
        grid_sizes=[21, 31, 41],  # Start with smaller sizes
        xi=100.0,
        Phi0=1e7
    )
    
    # Plot
    plot_benchmark_results(results)
    
    print("\n✅ Benchmark complete!")
