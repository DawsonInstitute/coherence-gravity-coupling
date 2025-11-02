#!/usr/bin/env python3
"""
Production Optimization Study

Systematic comparison of materials and geometries to map the signal landscape
and identify optimal configurations for experimental realization.

Study components:
1. Grid search for each material (YBCO, Rb-87 BEC, Nb cavity)
2. Refinement with global optimizer (DE) from grid optimum
3. Final local polish (Powell/L-BFGS-B)
4. Landscape visualization and comparison

Usage:
    python production_study.py --materials all --resolution 41
    python production_study.py --materials YBCO --quick
    
Author: GitHub Copilot (Claude Sonnet 4.5)
Date: October 2025
License: MIT
"""

import argparse
import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional
import platform
import sys
import subprocess
import numpy as np
from importlib import metadata as _metadata

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from optimize_geometry import GeometryOptimizer, save_optimization_results
from src.visualization.plot_utils import (
    plot_landscape_slice, plot_landscape_3d, plot_material_comparison, save_figure
)

# Output directory
RESULTS_DIR = Path("results/production_study")
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# Material configurations
MATERIALS = {
    'YBCO': {
        'name': 'YBCO Cuprate (90K)',
        'Phi0': 6.67e8,  # m⁻¹
        'xi': 100.0,
        'description': 'High-temperature superconductor'
    },
    'Rb87': {
        'name': 'Rb-87 BEC (100nK)',
        'Phi0': 3.65e6,  # m⁻¹
        'xi': 100.0,
        'description': 'Bose-Einstein condensate'
    },
    'Nb': {
        'name': 'Nb Cavity (9K)',
        'Phi0': 3.65e6,  # m⁻¹
        'xi': 100.0,
        'description': 'Superconducting RF cavity'
    }
}


class ProductionStudy:
    """Orchestrate comprehensive optimization study."""
    
    def __init__(self, resolution: int = 41, cache: bool = True, verbose: bool = False, jobs: int = 1):
        self.resolution = resolution
        self.cache = cache
        self.verbose = verbose
        self.jobs = jobs
        self.results = {}
    
    def run_material_study(
        self,
        material_key: str,
        grid_size: int = 5,
        run_refinement: bool = True
    ) -> Dict:
        """
        Complete optimization study for one material.
        
        Args:
            material_key: Material identifier (e.g., 'YBCO')
            grid_size: Grid points per dimension
            run_refinement: If True, refine grid optimum with DE + local
        
        Returns:
            Study results dictionary
        """
        material = MATERIALS[material_key]
        
        print(f"\n{'='*70}")
        print(f"MATERIAL STUDY: {material['name']}")
        print(f"{'='*70}")
        print(f"Φ₀ = {material['Phi0']:.2e} m⁻¹")
        print(f"ξ = {material['xi']}")
        print(f"Resolution: {self.resolution}³")
        print(f"Grid size: {grid_size}³ = {grid_size**3} points\n")
        
        optimizer = GeometryOptimizer(
            xi=material['xi'],
            Phi0=material['Phi0'],
            resolution=self.resolution,
            cache=self.cache,
            verbose=self.verbose
        )
        
        # Phase 1: Grid search
        print("Phase 1: Grid Search")
        print("-" * 70)
        
        grid_results = optimizer.grid_search(
            x_range=(-0.10, 0.10, grid_size),
            y_range=(-0.10, 0.10, grid_size),
            z_range=(-0.15, -0.05, grid_size),
            jobs=self.jobs,
            show_progress=True
        )
        
        grid_optimal_pos = grid_results['optimal_position']
        grid_optimal_tau = grid_results['optimal_delta_tau']
        
        study_results = {
            'material': material_key,
            'material_name': material['name'],
            'Phi0': material['Phi0'],
            'xi': material['xi'],
            'resolution': self.resolution,
            'grid_search': grid_results
        }
        
        if not run_refinement:
            return study_results
        
        # Phase 2: Global refinement (Differential Evolution)
        print("\nPhase 2: Global Refinement (Differential Evolution)")
        print("-" * 70)
        
        # Use grid optimum as center, search ±0.05 m around it
        de_bounds = [
            (grid_optimal_pos[0] - 0.05, grid_optimal_pos[0] + 0.05),
            (grid_optimal_pos[1] - 0.05, grid_optimal_pos[1] + 0.05),
            (grid_optimal_pos[2] - 0.05, grid_optimal_pos[2] + 0.05)
        ]
        
        de_results = optimizer.optimize_position(
            x0=np.array(grid_optimal_pos),
            bounds=de_bounds,
            method='DE'
        )
        
        study_results['de_refinement'] = de_results
        
        # Phase 3: Local polish (Powell)
        print("\nPhase 3: Local Polish (Powell)")
        print("-" * 70)
        
        powell_results = optimizer.optimize_position(
            x0=np.array(de_results['optimal_position']),
            method='Powell'
        )
        
        study_results['powell_polish'] = powell_results
        
        # Summary
        print(f"\n{'='*70}")
        print(f"STUDY SUMMARY: {material['name']}")
        print(f"{'='*70}")
        print(f"Grid search optimum:  Δτ = {grid_optimal_tau:.3e} N·m")
        print(f"  Position: ({grid_optimal_pos[0]:.4f}, {grid_optimal_pos[1]:.4f}, {grid_optimal_pos[2]:.4f}) m")
        print(f"\nDE refinement:        Δτ = {de_results['optimal_delta_tau']:.3e} N·m")
        print(f"  Position: ({de_results['optimal_position'][0]:.4f}, {de_results['optimal_position'][1]:.4f}, {de_results['optimal_position'][2]:.4f}) m")
        print(f"  Improvement: {abs(de_results['optimal_delta_tau'])/abs(grid_optimal_tau):.2f}×")
        print(f"\nPowell polish:        Δτ = {powell_results['optimal_delta_tau']:.3e} N·m")
        print(f"  Position: ({powell_results['optimal_position'][0]:.4f}, {powell_results['optimal_position'][1]:.4f}, {powell_results['optimal_position'][2]:.4f}) m")
        print(f"  Improvement: {abs(powell_results['optimal_delta_tau'])/abs(de_results['optimal_delta_tau']):.2f}×")
        print(f"\n🎯 Overall improvement: {abs(powell_results['optimal_delta_tau'])/abs(grid_optimal_tau):.2f}×")
        print(f"{'='*70}\n")
        
        return study_results
    
    @staticmethod
    def _get_env_metadata() -> Dict:
        """Collect lightweight reproducibility metadata for results files."""
        # Git information
        def _git(cmd: list[str]) -> Optional[str]:
            try:
                out = subprocess.check_output(cmd, stderr=subprocess.DEVNULL).decode().strip()
                return out or None
            except Exception:
                return None

        git_sha = _git(["git", "rev-parse", "HEAD"])  # current commit
        git_branch = _git(["git", "rev-parse", "--abbrev-ref", "HEAD"])  # branch name

        # Python and package versions
        py_ver = sys.version.split(" ")[0]
        np_ver = np.__version__
        try:
            sp_ver = _metadata.version("scipy")
        except Exception:
            sp_ver = None

        return {
            "python": py_ver,
            "numpy": np_ver,
            "scipy": sp_ver,
            "platform": {
                "system": platform.system(),
                "release": platform.release(),
                "machine": platform.machine(),
            },
            "git": {
                "sha": git_sha,
                "branch": git_branch,
            },
        }

    def generate_comparison_plots(self, all_results: List[Dict]) -> None:
        """Generate comparison plots across materials."""
        
        print("\n📊 Generating comparison plots...")
        
        # Material comparison bar chart
        comparison_data = []
        for result in all_results:
            comparison_data.append({
                'material': result['material_name'],
                'delta_tau': result.get('powell_polish', result.get('grid_search'))['optimal_delta_tau'],
                'config': result
            })
        
        output_path = RESULTS_DIR / "material_comparison"
        plot_material_comparison(comparison_data, output_path)
        
        # Landscape slices for each material
        for result in all_results:
            if 'grid_search' not in result or 'all_results' not in result['grid_search']:
                continue
            
            material_key = result['material']
            
            # Prepare grid data for plotting
            grid_data = {
                'grid_points': [],
                'objectives': []
            }
            
            for point_result in result['grid_search']['all_results']:
                grid_data['grid_points'].append(point_result['position'])
                grid_data['objectives'].append(point_result['delta_tau'])
            
            # Z-slice at z = -0.10 m
            output_path = RESULTS_DIR / f"landscape_{material_key}_z_slice"
            try:
                plot_landscape_slice(grid_data, 'z', -0.10, output_path)
            except Exception as e:
                print(f"   ⚠️  Could not plot landscape for {material_key}: {e}")
        
        print("   ✅ Comparison plots complete")
    
    def save_study_results(self, all_results: List[Dict], study_config: Optional[Dict] = None, environment: Optional[Dict] = None) -> Path:
        """Save complete study results with configuration and environment metadata."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        filename = f"production_study_{timestamp}.json"
        filepath = RESULTS_DIR / filename
        
        output = {
            'timestamp': timestamp,
            'resolution': self.resolution,
            'materials_studied': [r['material'] for r in all_results],
            'study_config': study_config or {},
            'environment': environment or self._get_env_metadata(),
            'results': all_results
        }
        
        with open(filepath, 'w') as f:
            json.dump(output, f, indent=2)
        
        print(f"\n💾 Study results saved: {filepath}")
        return filepath


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Production optimization study across materials",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument('--materials', nargs='+', default=['all'],
                       choices=['all', 'YBCO', 'Rb87', 'Nb'],
                       help='Materials to study')
    parser.add_argument('--resolution', type=int, default=41,
                       help='Grid resolution')
    parser.add_argument('--grid-size', type=int, default=5,
                       help='Grid points per dimension (5³ = 125 evaluations)')
    parser.add_argument('--jobs', type=int, default=1,
                       help='Parallel workers for grid search (1 = serial)')
    parser.add_argument('--quick', action='store_true',
                       help='Quick mode: grid search only, no refinement')
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable result caching')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose output')
    
    args = parser.parse_args()
    
    # Determine materials to study
    if 'all' in args.materials:
        materials_to_study = list(MATERIALS.keys())
    else:
        materials_to_study = args.materials
    
    print(f"\n{'='*70}")
    print("PRODUCTION OPTIMIZATION STUDY")
    print(f"{'='*70}")
    print(f"Materials: {', '.join(materials_to_study)}")
    print(f"Resolution: {args.resolution}³")
    print(f"Grid size: {args.grid_size}³ = {args.grid_size**3} points")
    print(f"Refinement: {'disabled (quick mode)' if args.quick else 'enabled (DE + Powell)'}")
    print(f"Caching: {'disabled' if args.no_cache else 'enabled'}")
    print(f"{'='*70}\n")
    
    study = ProductionStudy(
        resolution=args.resolution,
        cache=not args.no_cache,
        verbose=args.verbose,
        jobs=args.jobs
    )
    
    all_results = []
    
    for material_key in materials_to_study:
        try:
            result = study.run_material_study(
                material_key,
                grid_size=args.grid_size,
                run_refinement=not args.quick
            )
            all_results.append(result)
        except Exception as e:
            print(f"\n❌ Error studying {material_key}: {e}")
            import traceback
            traceback.print_exc()
    
    # Save results with run configuration and environment metadata
    study_config = {
        'materials': materials_to_study,
        'resolution': args.resolution,
        'grid_size': args.grid_size,
        'quick': args.quick,
        'cache': not args.no_cache,
        'jobs': args.jobs,
    }
    env_meta = study._get_env_metadata()
    study.save_study_results(all_results, study_config=study_config, environment=env_meta)
    
    # Generate comparison plots
    if len(all_results) > 1:
        study.generate_comparison_plots(all_results)
    
    print("\n✅ Production study complete!")
    print(f"📁 Results directory: {RESULTS_DIR}")


if __name__ == '__main__':
    main()
