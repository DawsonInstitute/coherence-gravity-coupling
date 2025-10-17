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
import matplotlib.pyplot as plt
import sys

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from field_equations.action import CoherenceGravityParams
from field_equations.einstein_coherence import compute_energy_cost_reduction


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
        target_reductions=[0.5, 0.1, 0.01, 1e-6],
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
    
    # Save results
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    with open(output_path / "parameter_scan.json", "w") as f:
        json.dump({
            'grid_scan': grid_results,
            'threshold_curves': curves,
            'benchmarks': benchmarks
        }, f, indent=2)
    
    print(f"\n✅ Results saved to {output_path / 'parameter_scan.json'}")
    print("=" * 70)
    
    return grid_results, curves, benchmarks


if __name__ == "__main__":
    run_parameter_scan()
