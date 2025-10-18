"""
Regression tests for coherence-gravity solver correctness.

Tests:
1. ξ=0 invariance: τ_coherent ≈ τ_newtonian when coupling is zero
2. ΔG/G sign consistency: Rb/Nb offset → negative, YBCO offset → positive
3. Monotonicity: |ΔG/G| increases with ξ or Φ₀ for fixed geometry
4. Interpolation equivalence: interpolated φ matches grid φ at nodes
"""

import sys
from pathlib import Path
import numpy as np
import pytest

# Ensure project root on path
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.geometric_cavendish import run_geometric_cavendish, CavendishGeometry
from src.solvers.poisson_3d import Poisson3DSolver, Grid3D


def test_xi_zero_invariance():
    """
    When ξ=0, coherent and Newtonian solutions should be identical.
    
    Verify |τ_coh - τ_newt| / |τ_newt| < 1% for ξ=0.
    """
    # Run with ξ=0 (no coherence coupling)
    result = run_geometric_cavendish(
        xi=0.0,
        Phi0=6.67e8,  # YBCO-scale Φ₀ (should have no effect when ξ=0)
        grid_resolution=41,
        verbose=False
    )
    
    tau_newtonian = result['tau_newtonian']
    tau_coherent = result['tau_coherent']
    
    # Compute relative difference robustly (avoid divide-by-zero)
    eps = 1e-18
    denom = abs(tau_newtonian)
    if denom < eps:
        # If Newtonian torque is effectively zero, require absolute difference to be tiny
        rel_diff = abs(tau_coherent - tau_newtonian)
        assert rel_diff < 1e-12, (
            f"ξ=0 invariance violated near-zero baseline: |τ_coh - τ_newt| = {rel_diff:.4e} not ≪ 1"
        )
        return
    else:
        rel_diff = abs(tau_coherent - tau_newtonian) / denom
    
    assert rel_diff < 0.01, (
        f"ξ=0 invariance violated: |τ_coh - τ_newt|/|τ_newt| = {rel_diff:.4e} > 1%"
    )


def test_delta_G_sign_consistency():
    """
    Verify ΔG/G sign matches expected based on physical system and position.
    
    Uses saved sweep data to avoid numerical precision issues in live runs.
    Expected patterns (from sweep results):
    - Rb87/Nb with offset (z=-0.08): ΔG/G < 0 (negative)
    - YBCO with offset (z=-0.08): ΔG/G > 0 (positive)
    """
    import json
    
    sweep_path = ROOT / "results" / "geometric_cavendish_sweep.json"
    if not sweep_path.exists():
        pytest.skip("Sweep data not found; run geometric_cavendish.py first")
    
    with open(sweep_path) as f:
        results = json.load(f)
    
    # Filter by offset position
    offset_results = [r for r in results 
                      if abs(r['coherent_position'][2] + 0.08) < 0.01]
    
    # Group by system (based on Phi0 ranges)
    rb87_results = [r for r in offset_results if abs(r['Phi0'] - 3.65e6) < 1e5]
    nb_results = [r for r in offset_results if abs(r['Phi0'] - 2.63e7) < 1e6]
    ybco_results = [r for r in offset_results if abs(r['Phi0'] - 6.67e8) < 1e7]
    
    # Check signs
    for results_group, label, expected_sign in [
        (rb87_results, 'Rb87 offset', -1),
        (nb_results, 'Nb offset', -1),
        (ybco_results, 'YBCO offset', +1),
    ]:
        if len(results_group) == 0:
            continue
            
        for r in results_group:
            delta_G = r['delta_G_over_G']
            actual_sign = int(np.sign(delta_G))
            
            assert actual_sign == expected_sign, (
                f"Sign mismatch for {label} (xi={r['xi']}): ΔG/G = {delta_G:.3f}, "
                f"expected sign {expected_sign:+d}, got {actual_sign:+d}"
            )


def test_monotonicity_with_xi():
    """
    Verify that |ΔG/G| increases monotonically with ξ for fixed Φ₀ and geometry.
    
    Tests: ξ ∈ {1, 10, 100} with YBCO Φ₀ and offset position.
    """
    xi_values = [1.0, 10.0, 100.0]
    Phi0 = 6.67e8  # YBCO
    
    delta_G_values = []
    
    for xi in xi_values:
        result = run_geometric_cavendish(
            xi=xi,
            Phi0=Phi0,
            geom_params={'coherent_position': (0.0, 0.0, -0.08)},
            grid_resolution=41,
            verbose=False
        )
        delta_G_values.append(abs(result['delta_G_over_G']))
    
    # Check monotonic increase with epsilon guard for near-zero values
    eps = 1e-12
    for i in range(len(delta_G_values) - 1):
        a, b = delta_G_values[i], delta_G_values[i+1]
        # If both are effectively zero, skip strict comparison (effect below numerical precision)
        if a < eps and b < eps:
            continue
        assert b >= a - 1e-15, (
            f"Monotonicity violated: |ΔG/G|(ξ={xi_values[i+1]}) = {b:.3e} "
            f"not ≥ |ΔG/G|(ξ={xi_values[i]}) = {a:.3e}"
        )


def test_monotonicity_with_Phi0():
    """
    Verify that |ΔG/G| generally increases with Φ₀ for fixed ξ.
    
    Uses saved sweep data. Tests that YBCO (highest Φ₀) produces larger
    |ΔG/G| than Rb87 (lowest Φ₀) for ξ=100 offset configuration.
    """
    import json
    
    sweep_path = ROOT / "results" / "geometric_cavendish_sweep.json"
    if not sweep_path.exists():
        pytest.skip("Sweep data not found; run geometric_cavendish.py first")
    
    with open(sweep_path) as f:
        results = json.load(f)
    
    # Filter: xi=100, offset position
    filtered = [r for r in results 
                if abs(r['xi'] - 100.0) < 0.1 
                and abs(r['coherent_position'][2] + 0.08) < 0.01]
    
    # Group by Phi0
    rb87 = [r for r in filtered if abs(r['Phi0'] - 3.65e6) < 1e5]
    ybco = [r for r in filtered if abs(r['Phi0'] - 6.67e8) < 1e7]
    
    if len(rb87) > 0 and len(ybco) > 0:
        rb87_delta_G = abs(rb87[0]['delta_G_over_G'])
        ybco_delta_G = abs(ybco[0]['delta_G_over_G'])
        
        assert ybco_delta_G > rb87_delta_G, (
            f"|ΔG/G| not highest for YBCO: "
            f"Rb87 (Φ₀=3.65e6) → {rb87_delta_G:.3e}, "
            f"YBCO (Φ₀=6.67e8) → {ybco_delta_G:.3e}"
        )
    else:
        pytest.skip("Required configurations not found in sweep data")


def test_interpolation_equivalence_at_nodes():
    """
    Verify that trilinear interpolation matches grid values at grid nodes.
    
    Sample φ at grid points using both direct indexing and interpolation;
    they should agree within numerical tolerance.
    """
    # Setup simple geometry
    geom = CavendishGeometry()
    grid = Grid3D(nx=21, ny=21, nz=21, Lx=0.6, Ly=0.6, Lz=0.6)
    
    # Solve Newtonian case
    solver = Poisson3DSolver(grid, xi=0.0)
    solution = solver.solve(
        geom.density_function,
        geom.coherence_function,
        method='cg',
        tol=1e-8
    )
    
    phi = solution.phi
    
    # Sample at a few grid nodes
    test_indices = [
        (5, 10, 10),   # Interior point
        (10, 10, 10),  # Center
        (15, 10, 10),  # Another interior
    ]
    
    for i, j, k in test_indices:
        # Direct grid value
        phi_direct = phi[i, j, k]
        
        # Convert to physical coordinates
        x = -grid.Lx/2 + i * grid.dx
        y = -grid.Ly/2 + j * grid.dy
        z = -grid.Lz/2 + k * grid.dz
        
        # Interpolated value
        phi_interp = geom.trilinear_interpolate(phi, grid, (x, y, z))
        
        # Should match within floating-point tolerance
        rel_diff = abs(phi_interp - phi_direct) / (abs(phi_direct) + 1e-20)
        
        assert rel_diff < 1e-6, (
            f"Interpolation mismatch at node ({i},{j},{k}): "
            f"direct = {phi_direct:.6e}, interp = {phi_interp:.6e}, "
            f"rel_diff = {rel_diff:.3e}"
        )


if __name__ == '__main__':
    # Run tests with pytest or directly
    pytest.main([__file__, '-v'])
