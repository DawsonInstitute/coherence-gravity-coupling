#!/usr/bin/env python3
"""
Geometry Optimization for Maximum Torque Signal

Uses scipy.optimize to find optimal geometry parameters that maximize
the coherence-modulated gravitational torque signal.

Optimization variables:
- Coherent system position (x, y, z)
- Test mass offset (optional)
- Source mass dimensions (optional)

Leverages result caching for fast iterations (~250× speedup on repeated configs).

Usage:
    python optimize_geometry.py --xi 100 --Phi0 1e8 --resolution 41
    python optimize_geometry.py --optimize-all --resolution 61
    
Author: GitHub Copilot (Claude Sonnet 4.5)
Date: October 2025
License: MIT
"""

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
from scipy.optimize import minimize, differential_evolution
import multiprocessing as mp
import platform
import sys
import subprocess
from importlib import metadata as _metadata

from examples.geometric_cavendish import run_geometric_cavendish

# Output directory
RESULTS_DIR = Path("results/optimization")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


class GeometryOptimizer:
    """Optimize geometry parameters to maximize torque signal."""
    
    def __init__(
        self,
        xi: float = 100.0,
        Phi0: float = 1e8,
        resolution: int = 41,
        cache: bool = True,
        verbose: bool = False
    ):
        """
        Initialize optimizer.
        
        Args:
            xi: Non-minimal coupling strength
            Phi0: Coherence field amplitude [m⁻¹]
            resolution: Grid resolution
            cache: Enable result caching
            verbose: Print detailed progress
        """
        self.xi = xi
        self.Phi0 = Phi0
        self.resolution = resolution
        self.cache = cache
        self.verbose = verbose
        
        # Optimization history
        self.history = []
        self.n_evals = 0
        
    def objective_position(self, params: np.ndarray) -> float:
        """
        Objective function for position optimization.
        
        Minimize -|Δτ| to maximize signal magnitude.
        
        Args:
            params: [x, y, z] position of coherent system [m]
        
        Returns:
            Negative absolute torque (for minimization)
        """
        x, y, z = params
        
        geom_params = {
            'coherent_position': [float(x), float(y), float(z)]
        }
        
        try:
            result = run_geometric_cavendish(
                xi=self.xi,
                Phi0=self.Phi0,
                geom_params=geom_params,
                grid_resolution=self.resolution,
                cache=self.cache,
                verbose=self.verbose
            )
            
            delta_tau = result['delta_tau']
            obj_value = -abs(delta_tau)
            
            self.n_evals += 1
            
            # Track best so far
            if len(self.history) == 0:
                best_so_far = obj_value
            else:
                best_so_far = min(self.history[-1]['best_so_far'], obj_value)
            
            self.history.append({
                'iteration': self.n_evals,
                'params': params.tolist(),
                'delta_tau': delta_tau,
                'objective': obj_value,
                'best_so_far': best_so_far
            })
            
            if self.n_evals % 5 == 0 or self.verbose:
                print(f"  [{self.n_evals}] pos=({x:.3f}, {y:.3f}, {z:.3f}) m, "
                      f"Δτ={delta_tau:.3e} N·m")
            
            return obj_value
            
        except Exception as e:
            print(f"  ⚠️  Evaluation failed: {e}")
            return 1e10  # Large penalty
    
    def save_trace(self, output_path: Path, method: str = "") -> None:
        """
        Save optimization trace to CSV and JSON.
        
        Args:
            output_path: Base path (without extension)
            method: Optimization method name for metadata
        """
        import csv
        
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Save CSV
        csv_path = output_path.with_suffix('.csv')
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=[
                'iteration', 'x', 'y', 'z', 'delta_tau', 'objective', 'best_so_far'
            ])
            writer.writeheader()
            for record in self.history:
                writer.writerow({
                    'iteration': record['iteration'],
                    'x': record['params'][0],
                    'y': record['params'][1],
                    'z': record['params'][2],
                    'delta_tau': record['delta_tau'],
                    'objective': record['objective'],
                    'best_so_far': record['best_so_far']
                })
        print(f"   Saved trace CSV: {csv_path}")
        
        # Save JSON
        json_path = output_path.with_suffix('.json')
        trace_data = {
            'method': method,
            'xi': self.xi,
            'Phi0': self.Phi0,
            'resolution': self.resolution,
            'n_evaluations': self.n_evals,
            'history': self.history
        }
        with open(json_path, 'w') as f:
            json.dump(trace_data, f, indent=2)
        print(f"   Saved trace JSON: {json_path}")
    
    def plot_trace(self, output_path: Optional[Path] = None, title: Optional[str] = None) -> None:
        """
        Plot optimization convergence trace.
        
        Args:
            output_path: Path to save plot (optional)
            title: Plot title (optional)
        """
        if len(self.history) == 0:
            print("   No history to plot")
            return
        
        try:
            from src.visualization.plot_utils import plot_optimization_trace
            
            fig = plot_optimization_trace(self.history, output_path, title)
            
            if output_path is None:
                import matplotlib.pyplot as plt
                plt.show()
        except ImportError as e:
            print(f"   ⚠️  Plotting unavailable: {e}")
    
    def optimize_position(
        self,
        x0: Optional[np.ndarray] = None,
        bounds: Optional[List[Tuple[float, float]]] = None,
        method: str = 'Nelder-Mead'
    ) -> Dict:
        """
        Optimize coherent system position.
        
        Args:
            x0: Initial guess [x, y, z] in meters
            bounds: [(x_min, x_max), (y_min, y_max), (z_min, z_max)]
            method: Optimization method ('Nelder-Mead', 'Powell', 'DE')
        
        Returns:
            Optimization results dictionary
        """
        if x0 is None:
            x0 = np.array([0.0, 0.0, -0.08])  # Default offset position
        
        if bounds is None:
            bounds = [(-0.15, 0.15), (-0.15, 0.15), (-0.20, 0.0)]
        
        print(f"\n{'='*70}")
        print("GEOMETRY OPTIMIZATION: Coherent System Position")
        print(f"{'='*70}")
        print(f"Method: {method}")
        print(f"Initial position: ({x0[0]:.3f}, {x0[1]:.3f}, {x0[2]:.3f}) m")
        print(f"Bounds: x∈[{bounds[0][0]:.2f}, {bounds[0][1]:.2f}], "
              f"y∈[{bounds[1][0]:.2f}, {bounds[1][1]:.2f}], "
              f"z∈[{bounds[2][0]:.2f}, {bounds[2][1]:.2f}] m")
        print(f"Physics: xi={self.xi}, Φ₀={self.Phi0:.2e} m⁻¹, res={self.resolution}³\n")
        
        # Reset history
        self.history = []
        self.n_evals = 0
        
        t0 = time.time()
        
        if method.upper() == 'DE':
            # Differential Evolution (global optimizer)
            result = differential_evolution(
                self.objective_position,
                bounds=bounds,
                maxiter=50,
                popsize=5,
                seed=42,
                disp=True
            )
        else:
            # Local optimizer
            result = minimize(
                self.objective_position,
                x0=x0,
                method=method,
                bounds=bounds if method in ['L-BFGS-B', 'TNC', 'SLSQP'] else None,
                options={'maxiter': 100, 'disp': True}
            )
        
        t_elapsed = time.time() - t0
        
        x_opt, y_opt, z_opt = result.x
        
        # Evaluate at optimal point to get full result
        geom_params_opt = {
            'coherent_position': [float(x_opt), float(y_opt), float(z_opt)]
        }
        
        result_opt = run_geometric_cavendish(
            xi=self.xi,
            Phi0=self.Phi0,
            geom_params=geom_params_opt,
            grid_resolution=self.resolution,
            cache=self.cache,
            verbose=False
        )
        
        # Evaluate at initial point for comparison
        geom_params_init = {
            'coherent_position': [float(x0[0]), float(x0[1]), float(x0[2])]
        }
        
        result_init = run_geometric_cavendish(
            xi=self.xi,
            Phi0=self.Phi0,
            geom_params=geom_params_init,
            grid_resolution=self.resolution,
            cache=self.cache,
            verbose=False
        )
        
        improvement = abs(result_opt['delta_tau']) / abs(result_init['delta_tau'])
        
        print(f"\n{'='*70}")
        print("OPTIMIZATION RESULTS")
        print(f"{'='*70}")
        print(f"Initial position: ({x0[0]:.4f}, {x0[1]:.4f}, {x0[2]:.4f}) m")
        print(f"  Δτ_init = {result_init['delta_tau']:.6e} N·m")
        print(f"\nOptimal position: ({x_opt:.4f}, {y_opt:.4f}, {z_opt:.4f}) m")
        print(f"  Δτ_opt  = {result_opt['delta_tau']:.6e} N·m")
        print(f"\n✨ Improvement: {improvement:.2f}× signal magnitude")
        print(f"Evaluations: {self.n_evals}")
        print(f"Time: {t_elapsed:.1f} s")
        print(f"Convergence: {'✅ SUCCESS' if result.success else '⚠️  PARTIAL'}")
        print(f"{'='*70}\n")
        
        return {
            'method': method,
            'xi': self.xi,
            'Phi0': self.Phi0,
            'resolution': self.resolution,
            'initial_position': x0.tolist(),
            'initial_delta_tau': result_init['delta_tau'],
            'optimal_position': result.x.tolist(),
            'optimal_delta_tau': result_opt['delta_tau'],
            'improvement_factor': improvement,
            'n_evaluations': self.n_evals,
            'elapsed_time': t_elapsed,
            'success': result.success,
            'history': self.history,
            'full_result': result_opt
        }
    
    def grid_search(
        self,
        x_range: Tuple[float, float, int],
        y_range: Tuple[float, float, int],
        z_range: Tuple[float, float, int],
        jobs: int = 1,
        show_progress: bool = True
    ) -> Dict:
        """
        Exhaustive grid search over position space.
        
        Args:
            x_range: (min, max, n_points)
            y_range: (min, max, n_points)
            z_range: (min, max, n_points)
        
        Returns:
            Grid search results with optimal position
        """
        x_vals = np.linspace(x_range[0], x_range[1], x_range[2])
        y_vals = np.linspace(y_range[0], y_range[1], y_range[2])
        z_vals = np.linspace(z_range[0], z_range[1], z_range[2])
        
        total_points = len(x_vals) * len(y_vals) * len(z_vals)
        
        print(f"\n{'='*70}")
        print("GRID SEARCH: Coherent System Position")
        print(f"{'='*70}")
        print(f"Grid: {len(x_vals)}×{len(y_vals)}×{len(z_vals)} = {total_points} points")
        print(f"x ∈ [{x_range[0]:.2f}, {x_range[1]:.2f}] m")
        print(f"y ∈ [{y_range[0]:.2f}, {y_range[1]:.2f}] m")
        print(f"z ∈ [{z_range[0]:.2f}, {z_range[1]:.2f}] m\n")
        
        results = []
        best_tau = 0
        best_pos = None
        
        t0 = time.time()

        # Prepare evaluation queue
        tasks: List[Tuple[int, float, float, float]] = []
        for i, x in enumerate(x_vals):
            for j, y in enumerate(y_vals):
                for k, z in enumerate(z_vals):
                    idx = i * len(y_vals) * len(z_vals) + j * len(z_vals) + k
                    tasks.append((idx, float(x), float(y), float(z)))

        def _evaluate_point(args: Tuple[int, float, float, float]) -> Dict:
            idx, x, y, z = args
            geom_params = {'coherent_position': [x, y, z]}
            try:
                r = run_geometric_cavendish(
                    xi=self.xi,
                    Phi0=self.Phi0,
                    geom_params=geom_params,
                    grid_resolution=self.resolution,
                    cache=self.cache,
                    verbose=False
                )
                return {
                    'idx': idx,
                    'position': [x, y, z],
                    'delta_tau': r['delta_tau']
                }
            except Exception as e:
                return {
                    'idx': idx,
                    'position': [x, y, z],
                    'delta_tau': 0.0,
                    'error': str(e)
                }

        # Optional progress bar
        try:
            from tqdm.auto import tqdm  # type: ignore
        except Exception:
            tqdm = None  # type: ignore

        if jobs and jobs > 1:
            # Use default context (fork on Linux) to allow nested function pickling
            with mp.Pool(processes=jobs) as pool:
                iterator = pool.imap_unordered(_evaluate_point, tasks, chunksize=1)
                if show_progress and tqdm is not None:
                    iterator = tqdm(iterator, total=total_points, desc="Grid evals", unit="pt")
                for item in iterator:
                    results.append(item)
                    dt = item['delta_tau']
                    if abs(dt) > abs(best_tau):
                        best_tau = dt
                        best_pos = item['position']
        else:
            iterator = tasks
            if show_progress and tqdm is not None:
                iterator = tqdm(iterator, total=total_points, desc="Grid evals", unit="pt")
            for item in iterator:
                out = _evaluate_point(item)
                results.append(out)
                dt = out['delta_tau']
                if abs(dt) > abs(best_tau):
                    best_tau = dt
                    best_pos = out['position']

        # Sort results by idx for deterministic ordering
        results_sorted = sorted(results, key=lambda r: r['idx'])
        # Compact result objects
        compact_results = [{'position': r['position'], 'delta_tau': r['delta_tau']} for r in results_sorted]
        
        t_elapsed = time.time() - t0
        
        print(f"\n{'='*70}")
        print("GRID SEARCH RESULTS")
        print(f"{'='*70}")
        print(f"Optimal position: ({best_pos[0]:.4f}, {best_pos[1]:.4f}, {best_pos[2]:.4f}) m")
        print(f"Δτ_max = {best_tau:.6e} N·m")
        print(f"Points evaluated: {total_points}")
        print(f"Time: {t_elapsed:.1f} s ({t_elapsed/total_points:.2f} s/point)")
        print(f"{'='*70}\n")
        
        return {
            'method': 'grid_search',
            'grid_shape': [len(x_vals), len(y_vals), len(z_vals)],
            'optimal_position': best_pos,
            'optimal_delta_tau': best_tau,
            'n_evaluations': total_points,
            'elapsed_time': t_elapsed,
            'all_results': compact_results
        }


def _get_env_metadata() -> Dict:
    """Collect minimal environment metadata for reproducibility."""
    def _git(cmd: List[str]) -> Optional[str]:
        try:
            out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode().strip()
            return out or None
        except Exception:
            return None

    git_sha = _git(["git", "rev-parse", "HEAD"])
    git_branch = _git(["git", "rev-parse", "--abbrev-ref", "HEAD"])

    try:
        sp_ver = _metadata.version("scipy")
    except Exception:
        sp_ver = None

    return {
        'python': sys.version.split(" ")[0],
        'numpy': np.__version__,
        'scipy': sp_ver,
        'platform': {
            'system': platform.system(),
            'release': platform.release(),
            'machine': platform.machine(),
        },
        'git': {'sha': git_sha, 'branch': git_branch},
    }


def save_optimization_results(results: Dict, name: str = "optimization") -> Path:
    """Save optimization results with timestamp."""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = f"{name}_{timestamp}.json"
    filepath = RESULTS_DIR / filename
    
    payload = dict(results)
    payload.setdefault('environment', _get_env_metadata())
    with open(filepath, 'w') as f:
        json.dump(payload, f, indent=2)
    
    print(f"💾 Results saved: {filepath}")
    return filepath


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Optimize geometry for maximum coherence-gravity signal",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--xi', type=float, default=100.0,
                       help='Non-minimal coupling strength')
    parser.add_argument('--Phi0', type=float, default=1e8,
                       help='Coherence field amplitude [m⁻¹]')
    parser.add_argument('--resolution', type=int, default=41,
                       help='Grid resolution')
    parser.add_argument('--method', type=str, default='Nelder-Mead',
                       choices=['Nelder-Mead', 'Powell', 'L-BFGS-B', 'DE'],
                       help='Optimization method')
    parser.add_argument('--initial-pos', nargs=3, type=float,
                       default=[0.0, 0.0, -0.08],
                       help='Initial position [x y z] in meters')
    parser.add_argument('--grid-search', action='store_true',
                       help='Use grid search instead of optimization')
    parser.add_argument('--grid-size', type=int, default=5,
                       help='Grid points per dimension for grid search')
    parser.add_argument('--jobs', type=int, default=1,
                       help='Parallel workers for grid search (1 = serial)')
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable result caching')
    parser.add_argument('--no-progress', action='store_true',
                       help='Disable progress bars')
    parser.add_argument('--plot', action='store_true',
                       help='Generate convergence plots')
    parser.add_argument('--save-trace', action='store_true',
                       help='Save convergence trace (CSV + JSON)')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    optimizer = GeometryOptimizer(
        xi=args.xi,
        Phi0=args.Phi0,
        resolution=args.resolution,
        cache=not args.no_cache,
        verbose=args.verbose
    )
    
    if args.grid_search:
        # Grid search
        results = optimizer.grid_search(
            x_range=(-0.10, 0.10, args.grid_size),
            y_range=(-0.10, 0.10, args.grid_size),
            z_range=(-0.15, -0.05, args.grid_size),
            jobs=args.jobs,
            show_progress=not args.no_progress
        )
        filepath = save_optimization_results(results, "grid_search")
        
        # For grid search, can optionally save landscape data
        if args.save_trace:
            print("\n📊 Grid results saved in JSON (use plot_utils for landscape viz)")
    else:
        # Gradient-based or derivative-free optimization
        x0 = np.array(args.initial_pos)
        results = optimizer.optimize_position(
            x0=x0,
            method=args.method
        )
        filepath = save_optimization_results(results, f"optimize_{args.method.lower()}")
        
        # Save convergence trace if requested
        if args.save_trace:
            print("\n📊 Saving convergence trace...")
            trace_path = filepath.parent / f"{filepath.stem}_trace"
            optimizer.save_trace(trace_path, method=args.method)
        
        # Generate plots if requested
        if args.plot:
            print("\n📊 Generating convergence plot...")
            plot_path = filepath.parent / f"{filepath.stem}_convergence"
            title = f"Optimization Convergence ({args.method}, ξ={args.xi})"
            optimizer.plot_trace(plot_path, title)
    
    print("✅ Optimization complete!")


if __name__ == '__main__':
    main()
