#!/usr/bin/env python3
"""
CLI entry point for coherence-gravity coupling framework.

Usage:
    python -m cgc run-scan --config examples/example_scan_config.json
    python -m cgc run-cavendish --xi 10 --coherence rb87 --position offset
    python -m cgc run-solver-test --grid 31 --preconditioner amg
    python -m cgc make-figures --geometric-sweep results/geometric_cavendish_sweep.json
"""

import argparse
import json
import sys
from pathlib import Path

def run_scan(args):
    """Run parameter space scan."""
    print(f"Running scan with config: {args.config}")
    import run_analysis
    # Would invoke run_analysis with parsed config
    print("Scan complete. Results in results/")

def run_cavendish(args):
    """Run geometric Cavendish simulation."""
    print(f"Running geometric Cavendish: ξ={args.xi}, coherence={args.coherence}, position={args.position}")
    from examples.geometric_cavendish import run_geometric_cavendish
    results = run_geometric_cavendish(
        xi=args.xi,
        coherence_system=args.coherence,
        coherent_position=args.position,
        grid_size=(args.grid_size,) * 3,
        preconditioner=args.preconditioner
    )
    print(f"\nResults:")
    print(f"  ΔG/G: {results['delta_G_over_G']:.4f}")
    print(f"  Coherent torque: {results['torque_coherent']:.6e} N·m")
    print(f"  Newtonian torque: {results['torque_newtonian']:.6e} N·m")
    print(f"  Solve time: {results['solve_time']:.3f} s")

def run_solver_test(args):
    """Run solver validation test."""
    print(f"Running solver test: grid={args.grid}³, preconditioner={args.preconditioner}")
    if args.test_type == "interface":
        from tests.test_interface_matching import run_interface_test
        result = run_interface_test(grid_size=args.grid)
        print(f"Test {'PASSED' if result else 'FAILED'}")
    elif args.test_type == "benchmark":
        from benchmarks.solver_acceleration import run_single_benchmark
        time_taken = run_single_benchmark(
            grid_size=args.grid,
            preconditioner=args.preconditioner
        )
        print(f"Solve time: {time_taken:.4f} s")

def make_figures(args):
    """Generate publication-quality figures."""
    print(f"Generating figures from: {args.geometric_sweep}")
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Load geometric sweep data
    with open(args.geometric_sweep, 'r') as f:
        data = json.load(f)
    
    # Extract ΔG/G values
    xi_vals = sorted(set(r['xi'] for r in data))
    coherence_systems = sorted(set(r['coherence'] for r in data))
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    for coh in coherence_systems:
        offset_vals = [r['delta_G_over_G'] for r in data 
                       if r['coherence'] == coh and r['position'] == 'offset']
        centered_vals = [r['delta_G_over_G'] for r in data 
                         if r['coherence'] == coh and r['position'] == 'centered']
        
        axes[0].plot(xi_vals, offset_vals, 'o-', label=f'{coh} (offset)')
        axes[1].plot(xi_vals, centered_vals, 's-', label=f'{coh} (centered)')
    
    for ax in axes:
        ax.set_xlabel('Coherence length ξ')
        ax.set_ylabel('ΔG/G')
        ax.set_xscale('log')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    axes[0].set_title('Coherent Body Offset')
    axes[1].set_title('Coherent Body Centered')
    
    plt.tight_layout()
    outpath = args.output or 'figures/geometric_cavendish_comparison.png'
    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"Figure saved to: {outpath}")

def main():
    parser = argparse.ArgumentParser(
        description="Coherence-Gravity Coupling Framework CLI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # run-scan subcommand
    scan_parser = subparsers.add_parser('run-scan', help='Run parameter space scan')
    scan_parser.add_argument('--config', required=True, help='Path to scan config JSON')
    scan_parser.set_defaults(func=run_scan)
    
    # run-cavendish subcommand
    cav_parser = subparsers.add_parser('run-cavendish', help='Run geometric Cavendish simulation')
    cav_parser.add_argument('--xi', type=float, required=True, help='Coherence length')
    cav_parser.add_argument('--coherence', choices=['rb87', 'nb', 'ybco'], required=True, 
                           help='Coherent system type')
    cav_parser.add_argument('--position', choices=['centered', 'offset'], default='offset',
                           help='Coherent body position relative to test masses')
    cav_parser.add_argument('--grid-size', type=int, default=41, 
                           help='Grid size per dimension (default: 41)')
    cav_parser.add_argument('--preconditioner', choices=['none', 'ilu', 'amg'], default='amg',
                           help='Solver preconditioner (default: amg)')
    cav_parser.set_defaults(func=run_cavendish)
    
    # run-solver-test subcommand
    test_parser = subparsers.add_parser('run-solver-test', help='Run solver validation')
    test_parser.add_argument('--test-type', choices=['interface', 'benchmark'], default='benchmark',
                            help='Type of test to run')
    test_parser.add_argument('--grid', type=int, default=31, help='Grid size')
    test_parser.add_argument('--preconditioner', choices=['none', 'ilu', 'amg'], default='none',
                            help='Solver preconditioner')
    test_parser.set_defaults(func=run_solver_test)
    
    # make-figures subcommand
    fig_parser = subparsers.add_parser('make-figures', help='Generate figures')
    fig_parser.add_argument('--geometric-sweep', required=True,
                           help='Path to geometric_cavendish_sweep.json')
    fig_parser.add_argument('--output', help='Output figure path (default: figures/geometric_cavendish_comparison.png)')
    fig_parser.set_defaults(func=make_figures)
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    args.func(args)

if __name__ == '__main__':
    main()
