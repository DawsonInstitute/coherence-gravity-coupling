"""
Parameter space scanner for coherence-gravity coupling.

Maps (ξ, Φ₀) parameter space to effective coupling ratio G_eff/G
and energy cost reduction factors.

This is the key analysis: for what coherence amplitudes and coupling strengths
can we achieve significant curvature energy cost reduction?
"""

import numpy as np
import json
from pathlib import Path
from typing import Dict, List, Tuple
from datetime import datetime
import sys

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from field_equations.action import CoherenceGravityParams
from field_equations.weak_field import compute_energy_cost_reduction

# Try to import matplotlib, but don't fail if not available
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available, plots will not be generated")


class ParameterSpaceScanner:
    """
    Systematically explore (ξ, Φ₀) parameter space.
    
    Goal: Find regions where G_eff << G, making curvature "cheap".
    """
    
    def __init__(self, params: CoherenceGravityParams):
        """
        Initialize scanner.
        
        Args:
            params: Base parameters (G, c, etc.)
        """
        self.params = params
        
    def scan_xi_phi0_grid(self, xi_range: Tuple[float, float],
                          phi0_range: Tuple[float, float],
                          n_xi: int = 50,
                          n_phi0: int = 50) -> Dict:
        """
        Scan grid of (ξ, Φ₀) values.
        
        Args:
            xi_range: (xi_min, xi_max) for coupling strength
            phi0_range: (Phi0_min, Phi0_max) for coherence amplitude
            n_xi: Number of xi grid points
            n_phi0: Number of Phi0 grid points
            
        Returns:
            Dictionary with grid data and results
        """
        xi_vals = np.logspace(np.log10(xi_range[0]), np.log10(xi_range[1]), n_xi)
        phi0_vals = np.logspace(np.log10(phi0_range[0]), np.log10(phi0_range[1]), n_phi0)
        
        # Result arrays
        G_eff_ratio = np.zeros((n_xi, n_phi0))
        energy_reduction = np.zeros((n_xi, n_phi0))
        suppression_factor = np.zeros((n_xi, n_phi0))
        
        for i, xi in enumerate(xi_vals):
            for j, phi0 in enumerate(phi0_vals):
                result = compute_energy_cost_reduction(phi0, xi, self.params.G)
                
                G_eff_ratio[i, j] = result['G_eff_ratio']
                energy_reduction[i, j] = result['energy_cost_reduction']
                suppression_factor[i, j] = result['suppression_factor']
        
        return {
            'xi_vals': xi_vals.tolist(),
            'phi0_vals': phi0_vals.tolist(),
            'G_eff_ratio': G_eff_ratio.tolist(),
            'energy_reduction': energy_reduction.tolist(),
            'suppression_factor': suppression_factor.tolist(),
            'n_xi': n_xi,
            'n_phi0': n_phi0
        }
    
    def find_threshold_curves(self, target_reductions: List[float],
                              xi_range: Tuple[float, float],
                              n_points: int = 1000) -> Dict:
        """
        Find curves in (ξ, Φ₀) space for target G_eff/G ratios.
        
        For each target reduction factor, compute the required Φ₀(ξ).
        
        Args:
            target_reductions: List of desired G_eff/G ratios (e.g., [0.5, 0.1, 0.01])
            xi_range: Range of ξ values to scan
            n_points: Number of points per curve
            
        Returns:
            Dictionary with curves for each target
        """
        xi_vals = np.logspace(np.log10(xi_range[0]), np.log10(xi_range[1]), n_points)
        
        curves = {}
        for target in target_reductions:
            # From G_eff/G = 1/(1 + 8πGξΦ₀²) = target
            # → 8πGξΦ₀² = 1/target - 1
            # → Φ₀ = sqrt[(1/target - 1) / (8πGξ)]
            
            phi0_vals = np.sqrt((1/target - 1) / (8*np.pi*self.params.G*xi_vals))
            
            curves[f'target_{target:.1e}'] = {
                'xi': xi_vals.tolist(),
                'phi0': phi0_vals.tolist(),
                'target_ratio': target
            }
        
        return curves
    
    def assess_realizability(self, xi: float, phi0: float) -> Dict:
        """
        Assess whether (ξ, Φ₀) parameters are physically realizable.
        
        Compares required Φ₀ to known coherent systems:
        - BEC: |Ψ|² ~ 10¹⁴ cm⁻³
        - Superconductor: n_Cooper ~ 10²² cm⁻³
        - Plasma: n_e ~ 10²⁰ cm⁻³
        
        Args:
            xi: Coupling strength
            phi0: Coherence amplitude
            
        Returns:
            Dictionary with realizability assessment
        """
        # Estimate coherence amplitude in known systems
        # (These are rough order-of-magnitude estimates)
        
        # BEC condensate density
        n_BEC = 1e14 * 1e6  # Convert cm⁻³ to m⁻³
        Phi_BEC = np.sqrt(n_BEC)  # Simple estimate: Φ ~ √n
        
        # Superconductor Cooper pair density
        n_Cooper = 1e22 * 1e6
        Phi_SC = np.sqrt(n_Cooper)
        
        # High-density plasma
        n_plasma = 1e20 * 1e6
        Phi_plasma = np.sqrt(n_plasma)
        
        # Compare to required Phi0
        gap_BEC = phi0 / Phi_BEC
        gap_SC = phi0 / Phi_SC
        gap_plasma = phi0 / Phi_plasma
        
        # Find closest system
        gaps = {'BEC': gap_BEC, 'Superconductor': gap_SC, 'Plasma': gap_plasma}
        closest_system = min(gaps.items(), key=lambda x: x[1])[0]
        min_gap = gaps[closest_system]
        
        return {
            'required_phi0': phi0,
            'xi': xi,
            'BEC_phi0': Phi_BEC,
            'SC_phi0': Phi_SC,
            'plasma_phi0': Phi_plasma,
            'gap_BEC': gap_BEC,
            'gap_SC': gap_SC,
            'gap_plasma': gap_plasma,
            'closest_system': closest_system,
            'min_gap_factor': min_gap,
            'realizable': min_gap <= 10.0  # Within 10× is "potentially realizable"
        }
    
    def benchmark_configurations(self) -> List[Dict]:
        """
        Test a set of benchmark (ξ, Φ₀) configurations.
        
        Returns:
            List of configuration results
        """
        benchmarks = [
            # (ξ, Φ₀, description)
            (1.0, 1e10, "Weak coupling, low coherence"),
            (1.0, 1e15, "Weak coupling, BEC-scale"),
            (1.0, 1e20, "Weak coupling, superconductor-scale"),
            (10.0, 1e15, "Strong coupling, BEC-scale"),
            (100.0, 1e15, "Very strong coupling, BEC-scale"),
            (1.0, 1e25, "Weak coupling, extreme coherence"),
        ]
        
        results = []
        for xi, phi0, desc in benchmarks:
            reduction = compute_energy_cost_reduction(phi0, xi, self.params.G)
            realizability = self.assess_realizability(xi, phi0)
            
            results.append({
                'description': desc,
                'xi': float(xi),
                'phi0': float(phi0),
                'G_eff_ratio': float(reduction['G_eff_ratio']),
                'energy_reduction': float(reduction['energy_cost_reduction']),
                'gap_factor': float(realizability['min_gap_factor']),
                'closest_system': realizability['closest_system'],
                'realizable': bool(realizability['realizable'])
            })
        
        return results


def run_parameter_scan(output_dir: str = "results"):
    """
    Run comprehensive parameter space scan and save results.
    
    Args:
        output_dir: Directory to save results
    """
    print("=" * 70)
    print("COHERENCE-GRAVITY COUPLING PARAMETER SPACE SCAN")
    print("=" * 70)
    
    params = CoherenceGravityParams()
    scanner = ParameterSpaceScanner(params)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # 1. Grid scan
    print("\n1. Scanning (ξ, Φ₀) grid...")
    grid_results = scanner.scan_xi_phi0_grid(
        xi_range=(0.01, 1000.0),
        phi0_range=(1e10, 1e30),
        n_xi=100,
        n_phi0=100
    )
    
    # 2. Threshold curves
    print("2. Computing threshold curves...")
    curves = scanner.find_threshold_curves(
        target_reductions=[0.5, 0.1, 0.01, 1e-6, 1e-12, 1e-24],
        xi_range=(0.01, 1000.0)
    )
    
    # 3. Benchmark configurations
    print("3. Testing benchmark configurations...")
    benchmarks = scanner.benchmark_configurations()
    
    print("\nBenchmark Results:")
    print("-" * 70)
    for i, bench in enumerate(benchmarks):
        print(f"{i+1}. {bench['description']}")
        print(f"   ξ = {bench['xi']:.1e}, Φ₀ = {bench['phi0']:.1e}")
        print(f"   G_eff/G = {bench['G_eff_ratio']:.3e} (energy ↓ {1/bench['G_eff_ratio']:.1e}×)")
        print(f"   Gap to {bench['closest_system']}: {bench['gap_factor']:.1e}×")
        print(f"   Realizable: {'✅ YES' if bench['realizable'] else '❌ NO'}")
        print()
    
    # 4. Find optimal realizable configuration
    print("4. Finding optimal realizable configuration...")
    realizable = [b for b in benchmarks if b['realizable']]
    if realizable:
        best = min(realizable, key=lambda x: x['G_eff_ratio'])
        print(f"\n✅ BEST REALIZABLE CONFIGURATION:")
        print(f"   {best['description']}")
        print(f"   G_eff/G = {best['G_eff_ratio']:.3e}")
        print(f"   Energy cost reduction: {1/best['G_eff_ratio']:.1e}×")
    else:
        print("\n❌ No realizable configurations found in benchmark set")
        print("   Need to explore higher coherence amplitudes or coupling strengths")
    
    # 5. Generate plots if matplotlib available
    if HAS_MATPLOTLIB:
        print("\n5. Generating plots...")
        plot_parameter_space(grid_results, curves, output_path)
        plot_benchmark_comparison(benchmarks, output_path)
        print(f"   ✅ Plots saved to {output_path}/")
    else:
        print("\n5. Skipping plots (matplotlib not available)")
    
    # 6. Save results with metadata
    timestamp = datetime.now().isoformat()
    metadata = {
        'timestamp': timestamp,
        'G': params.G,
        'c': params.c,
        'hbar': params.hbar,
        'scan_parameters': {
            'xi_range': [0.01, 1000.0],
            'phi0_range': [1e10, 1e30],
            'n_xi': 100,
            'n_phi0': 100
        }
    }
    
    with open(output_path / "parameter_scan.json", "w") as f:
        json.dump({
            'metadata': metadata,
            'grid_scan': grid_results,
            'threshold_curves': curves,
            'benchmarks': benchmarks
        }, f, indent=2)
    
    print(f"\n✅ Results saved to {output_path / 'parameter_scan.json'}")
    print("=" * 70)
    
    return grid_results, curves, benchmarks


def plot_parameter_space(grid_results: Dict, curves: Dict, output_path: Path):
    """
    Generate 2D parameter space plots.
    
    Args:
        grid_results: Grid scan data
        curves: Threshold curves
        output_path: Directory to save plots
    """
    if not HAS_MATPLOTLIB:
        return
    
    xi_vals = np.array(grid_results['xi_vals'])
    phi0_vals = np.array(grid_results['phi0_vals'])
    G_eff_ratio = np.array(grid_results['G_eff_ratio'])
    
    # Create meshgrid
    Xi, Phi0 = np.meshgrid(xi_vals, phi0_vals, indexing='ij')
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Plot 1: G_eff/G ratio (log scale)
    im1 = ax1.contourf(Xi, Phi0, G_eff_ratio, levels=50, 
                       norm=mcolors.LogNorm(vmin=1e-30, vmax=1.0),
                       cmap='RdYlGn_r')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Coupling ξ', fontsize=12)
    ax1.set_ylabel('Coherence Φ₀ [m⁻¹]', fontsize=12)
    ax1.set_title('Effective Coupling Ratio G_eff/G', fontsize=14, fontweight='bold')
    
    # Add threshold curves
    for target_key, curve_data in curves.items():
        ax1.plot(curve_data['xi'], curve_data['phi0'], 'k--', linewidth=2,
                label=f"G_eff/G = {curve_data['target_ratio']:.0e}")
    
    ax1.legend(fontsize=8, loc='upper left')
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label('G_eff / G', fontsize=11)
    
    # Plot 2: Energy reduction factor (inverted)
    energy_reduction = 1.0 / G_eff_ratio
    im2 = ax2.contourf(Xi, Phi0, energy_reduction, levels=50,
                       norm=mcolors.LogNorm(vmin=1.0, vmax=1e30),
                       cmap='plasma')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Coupling ξ', fontsize=12)
    ax2.set_ylabel('Coherence Φ₀ [m⁻¹]', fontsize=12)
    ax2.set_title('Energy Cost Reduction Factor', fontsize=14, fontweight='bold')
    cbar2 = plt.colorbar(im2, ax=ax2)
    cbar2.set_label('E_standard / E_coherent', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_path / 'parameter_space.png', dpi=200)
    print(f"   Saved: parameter_space.png")
    plt.close()


def plot_benchmark_comparison(benchmarks: List[Dict], output_path: Path):
    """
    Plot benchmark configurations comparison.
    
    Args:
        benchmarks: List of benchmark results
        output_path: Directory to save plots
    """
    if not HAS_MATPLOTLIB:
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Extract data
    labels = [b['description'] for b in benchmarks]
    G_ratios = [b['G_eff_ratio'] for b in benchmarks]
    gaps = [b['gap_factor'] for b in benchmarks]
    realizable = [b['realizable'] for b in benchmarks]
    
    colors = ['green' if r else 'red' for r in realizable]
    
    # Plot 1: G_eff/G ratios
    x = np.arange(len(benchmarks))
    bars1 = ax1.bar(x, G_ratios, color=colors, alpha=0.7, edgecolor='black')
    ax1.set_yscale('log')
    ax1.set_ylabel('G_eff / G', fontsize=12)
    ax1.set_title('Effective Coupling Suppression', fontsize=13, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(range(1, len(benchmarks)+1))
    ax1.set_xlabel('Configuration #', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    ax1.axhline(1.0, color='k', linestyle='--', linewidth=1, alpha=0.5)
    
    # Plot 2: Gap factors
    bars2 = ax2.bar(x, gaps, color=colors, alpha=0.7, edgecolor='black')
    ax2.set_yscale('log')
    ax2.set_ylabel('Gap to Nearest System', fontsize=12)
    ax2.set_title('Realizability Gap', fontsize=13, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(range(1, len(benchmarks)+1))
    ax2.set_xlabel('Configuration #', fontsize=12)
    ax2.axhline(10.0, color='orange', linestyle=':', linewidth=2, 
                label='Realizability threshold (10×)')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_path / 'benchmark_comparison.png', dpi=150)
    print(f"   Saved: benchmark_comparison.png")
    plt.close()
    
    # Create detailed text summary
    with open(output_path / 'benchmark_summary.txt', 'w') as f:
        f.write("="*70 + "\n")
        f.write("COHERENCE-GRAVITY COUPLING: BENCHMARK CONFIGURATIONS\n")
        f.write("="*70 + "\n\n")
        
        for i, b in enumerate(benchmarks):
            f.write(f"Configuration {i+1}: {b['description']}\n")
            f.write(f"  ξ = {b['xi']:.2e}\n")
            f.write(f"  Φ₀ = {b['phi0']:.2e} m⁻¹\n")
            f.write(f"  G_eff/G = {b['G_eff_ratio']:.6e}\n")
            f.write(f"  Energy reduction: {1/b['G_eff_ratio']:.3e}×\n")
            f.write(f"  Closest system: {b['closest_system']}\n")
            f.write(f"  Gap factor: {b['gap_factor']:.3e}×\n")
            f.write(f"  Realizable: {'YES ✅' if b['realizable'] else 'NO ❌'}\n")
            f.write("\n" + "-"*70 + "\n\n")
    
    print(f"   Saved: benchmark_summary.txt")


if __name__ == "__main__":
    run_parameter_scan()
