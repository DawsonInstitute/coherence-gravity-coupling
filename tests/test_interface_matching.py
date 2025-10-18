"""
Interface Matching Validation for Discontinuous G_eff

Tests numerical accuracy at sharp G_eff interfaces by comparing to
analytic 1D slab solutions.

Physical matching conditions at interface z=z₀:
1. Continuity of potential: φ(z₀⁻) = φ(z₀⁺)
2. Continuity of normal flux: G_eff(z₀⁻) ∂φ/∂z|_{z₀⁻} = G_eff(z₀⁺) ∂φ/∂z|_{z₀⁺}

Author: GitHub Copilot (Claude Sonnet 4.5)
License: MIT
"""

import numpy as np
import sys
from pathlib import Path
from typing import Tuple, Callable

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.solvers.poisson_3d import Poisson3DSolver, Grid3D

# Constants
G_SI = 6.674e-11  # m³/(kg·s²)

try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def analytic_slab_solution(
    z: np.ndarray,
    rho_top: float,
    rho_slab: float,
    rho_bot: float,
    G_eff_top: float,
    G_eff_slab: float,
    G_eff_bot: float,
    z_interface_1: float,
    z_interface_2: float
) -> np.ndarray:
    """
    Analytic 1D Poisson solution for layered slab:
    
    Region 1 (z < z₁): G_eff = G₁, ρ = ρ₁
    Region 2 (z₁ < z < z₂): G_eff = G₂, ρ = ρ₂
    Region 3 (z > z₂): G_eff = G₃, ρ = ρ₃
    
    In each region: d²φ/dz² = 4πG ρ / G_eff
    
    For uniform ρ in each layer:
    φ(z) = A_i z² + B_i z + C_i  in region i
    
    Matching conditions at interfaces give the constants.
    """
    phi = np.zeros_like(z)
    
    # Simplification: assume ρ=0 except in slab
    # Slab: ρ = rho_slab, others: ρ = 0
    
    # Region 3 (z > z₂): φ₃ = 0 (reference)
    # Region 2 (slab): d²φ/dz² = 4πG ρ_slab / G_eff_slab
    #   φ₂ = (2πG ρ_slab / G_eff_slab) z² + B₂ z + C₂
    # Region 1 (z < z₁): d²φ/dz² = 0 → φ₁ = B₁ z + C₁
    
    # Source term in slab
    k_slab = 2.0 * np.pi * G_SI * rho_slab / G_eff_slab
    
    # Boundary condition at z → +∞: φ = 0
    # At z = z₂:
    #   φ₂(z₂) = 0 (continuity with region 3)
    #   G_eff_slab * dφ₂/dz|_{z₂} = G_eff_top * dφ₃/dz|_{z₂} = 0
    
    # Working backward from z₂:
    # φ₂(z₂) = 0
    # dφ₂/dz|_{z₂} = 0 (no flux above)
    # dφ₂/dz = 2 k_slab z + B₂ → B₂ = -2 k_slab z₂
    # φ₂(z) = k_slab (z² - 2z₂ z + z₂²) = k_slab (z - z₂)²
    
    # At z = z₁:
    # φ₁(z₁) = φ₂(z₁) = k_slab (z₁ - z₂)²
    # G_eff_bot * dφ₁/dz = G_eff_slab * dφ₂/dz|_{z₁}
    # dφ₂/dz|_{z₁} = 2 k_slab (z₁ - z₂)
    # B₁ = (G_eff_slab / G_eff_bot) * 2 k_slab (z₁ - z₂)
    # C₁ = φ₂(z₁) - B₁ z₁
    
    B1 = (G_eff_slab / G_eff_bot) * 2.0 * k_slab * (z_interface_1 - z_interface_2)
    phi_at_z1 = k_slab * (z_interface_1 - z_interface_2)**2
    C1 = phi_at_z1 - B1 * z_interface_1
    
    # Evaluate
    for i, zi in enumerate(z):
        if zi < z_interface_1:
            # Region 1 (below slab)
            phi[i] = B1 * zi + C1
        elif zi < z_interface_2:
            # Region 2 (inside slab)
            phi[i] = k_slab * (zi - z_interface_2)**2
        else:
            # Region 3 (above slab)
            phi[i] = 0.0
    
    return phi


def test_1d_slab_interface():
    """
    Test 3D solver on 1D slab problem with sharp G_eff interface.
    """
    print("\n" + "="*70)
    print("INTERFACE MATCHING VALIDATION")
    print("1D Slab with Discontinuous G_eff")
    print("="*70)
    
    # Parameters
    z_interface_1 = -0.10  # Bottom of slab [m]
    z_interface_2 = 0.10   # Top of slab [m]
    
    rho_slab = 1000.0  # kg/m³ (water density)
    Phi_slab = 1e7     # m⁻¹ (BEC scale)
    xi = 100.0
    
    G_eff_slab = G_SI / (1.0 + 8.0 * np.pi * G_SI * xi * Phi_slab**2)
    G_eff_out = G_SI
    
    print(f"\nConfiguration:")
    print(f"   Slab: z ∈ [{z_interface_1}, {z_interface_2}] m")
    print(f"   ρ_slab = {rho_slab} kg/m³")
    print(f"   Φ_slab = {Phi_slab:.2e} m⁻¹")
    print(f"   ξ = {xi}")
    print(f"   G_eff inside slab: {G_eff_slab/G_SI:.3e} × G")
    print(f"   G_eff outside slab: {G_eff_out/G_SI:.3e} × G")
    print(f"   Contrast ratio: {G_eff_out/G_eff_slab:.1f}×")
    
    # 3D solver (1D problem via symmetry)
    grid = Grid3D(
        nx=5, ny=5, nz=81,  # Fine in z-direction only
        Lx=0.1, Ly=0.1, Lz=0.8
    )
    
    def rho_func(x, y, z):
        if z_interface_1 < z < z_interface_2:
            return rho_slab
        return 0.0
    
    def Phi_func(x, y, z):
        if z_interface_1 < z < z_interface_2:
            return Phi_slab
        return 0.0
    
    print(f"\n3D Solver:")
    print(f"   Grid: {grid.nx}×{grid.ny}×{grid.nz}")
    solver = Poisson3DSolver(grid, xi=xi)
    solution = solver.solve(rho_func, Phi_func, method='cg', tol=1e-10)
    
    # Extract z-profile along center
    i_c = grid.nx // 2
    j_c = grid.ny // 2
    
    z_vals = []
    phi_numeric = []
    G_eff_vals = []
    
    for k in range(grid.nz):
        x, y, z = grid.coord(i_c, j_c, k)
        z_vals.append(z)
        phi_numeric.append(solution.phi[i_c, j_c, k])
        G_eff_vals.append(solution.G_eff[i_c, j_c, k])
    
    z_vals = np.array(z_vals)
    phi_numeric = np.array(phi_numeric)
    G_eff_vals = np.array(G_eff_vals)
    
    # Analytic solution
    phi_analytic = analytic_slab_solution(
        z_vals,
        rho_top=0.0,
        rho_slab=rho_slab,
        rho_bot=0.0,
        G_eff_top=G_eff_out,
        G_eff_slab=G_eff_slab,
        G_eff_bot=G_eff_out,
        z_interface_1=z_interface_1,
        z_interface_2=z_interface_2
    )
    
    # Compute error
    error = np.abs(phi_numeric - phi_analytic)
    relative_error = error / (np.abs(phi_analytic).max() + 1e-20)
    
    max_error = error.max()
    max_rel_error = relative_error.max()
    rms_error = np.sqrt(np.mean(error**2))
    
    print(f"\nError Analysis:")
    print(f"   Max absolute error: {max_error:.3e} m²/s²")
    print(f"   Max relative error: {max_rel_error:.3e}")
    print(f"   RMS error: {rms_error:.3e} m²/s²")
    
    # Check interface matching
    k1 = np.argmin(np.abs(z_vals - z_interface_1))
    k2 = np.argmin(np.abs(z_vals - z_interface_2))
    
    # Potential continuity
    phi_below_1 = phi_numeric[max(0, k1-1)]
    phi_above_1 = phi_numeric[min(grid.nz-1, k1+1)]
    delta_phi_1 = abs(phi_above_1 - phi_below_1) / grid.dz
    
    phi_below_2 = phi_numeric[max(0, k2-1)]
    phi_above_2 = phi_numeric[min(grid.nz-1, k2+1)]
    delta_phi_2 = abs(phi_above_2 - phi_below_2) / grid.dz
    
    print(f"\nInterface Matching:")
    print(f"   Interface 1 (z={z_interface_1}):")
    print(f"      Δφ/Δz: {delta_phi_1:.3e} m/s²")
    print(f"   Interface 2 (z={z_interface_2}):")
    print(f"      Δφ/Δz: {delta_phi_2:.3e} m/s²")
    
    # Flux continuity (G_eff * dφ/dz should be continuous)
    dz = grid.dz
    flux_below_1 = G_eff_vals[k1-1] * (phi_numeric[k1] - phi_numeric[k1-2]) / (2*dz)
    flux_above_1 = G_eff_vals[k1+1] * (phi_numeric[k1+2] - phi_numeric[k1]) / (2*dz)
    flux_mismatch_1 = abs(flux_above_1 - flux_below_1) / (abs(flux_below_1) + 1e-20)
    
    flux_below_2 = G_eff_vals[k2-1] * (phi_numeric[k2] - phi_numeric[k2-2]) / (2*dz)
    flux_above_2 = G_eff_vals[k2+1] * (phi_numeric[k2+2] - phi_numeric[k2]) / (2*dz)
    flux_mismatch_2 = abs(flux_above_2 - flux_below_2) / (abs(flux_below_2) + 1e-20)
    
    print(f"   Flux mismatch at interface 1: {flux_mismatch_1:.3e}")
    print(f"   Flux mismatch at interface 2: {flux_mismatch_2:.3e}")
    
    # Pass/fail criteria
    tolerance_abs = 1e-6  # m²/s²
    tolerance_rel = 1e-2  # 1%
    tolerance_flux = 0.1  # 10%
    
    pass_abs = max_error < tolerance_abs
    pass_rel = max_rel_error < tolerance_rel
    pass_flux = max(flux_mismatch_1, flux_mismatch_2) < tolerance_flux
    
    print(f"\nTest Results:")
    print(f"   Absolute error < {tolerance_abs}: {'✅ PASS' if pass_abs else '❌ FAIL'}")
    print(f"   Relative error < {tolerance_rel}: {'✅ PASS' if pass_rel else '❌ FAIL'}")
    print(f"   Flux matching < {tolerance_flux}: {'✅ PASS' if pass_flux else '❌ FAIL'}")
    
    overall_pass = pass_abs or pass_rel and pass_flux
    
    if overall_pass:
        print(f"\n✅ OVERALL: PASS")
    else:
        print(f"\n❌ OVERALL: FAIL")
    
    assert overall_pass, f"Interface matching failed: max_error={max_error:.3e}, max_rel_error={max_rel_error:.3e}, flux_mismatch={max(flux_mismatch_1, flux_mismatch_2):.3e}"
    
    # Plot
    if HAS_MATPLOTLIB:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Potential
        axes[0,0].plot(z_vals, phi_analytic, 'k-', linewidth=2, label='Analytic')
        axes[0,0].plot(z_vals, phi_numeric, 'r--', linewidth=1.5, label='Numeric')
        axes[0,0].axvline(z_interface_1, color='gray', linestyle=':', alpha=0.5)
        axes[0,0].axvline(z_interface_2, color='gray', linestyle=':', alpha=0.5)
        axes[0,0].set_xlabel('z [m]')
        axes[0,0].set_ylabel('φ [m²/s²]')
        axes[0,0].set_title('Gravitational Potential')
        axes[0,0].legend()
        axes[0,0].grid(True, alpha=0.3)
        
        # Error
        axes[0,1].plot(z_vals, error, 'b-', linewidth=1.5)
        axes[0,1].axvline(z_interface_1, color='gray', linestyle=':', alpha=0.5)
        axes[0,1].axvline(z_interface_2, color='gray', linestyle=':', alpha=0.5)
        axes[0,1].set_xlabel('z [m]')
        axes[0,1].set_ylabel('|φ_numeric - φ_analytic| [m²/s²]')
        axes[0,1].set_title('Absolute Error')
        axes[0,1].set_yscale('log')
        axes[0,1].grid(True, alpha=0.3)
        
        # G_eff profile
        axes[1,0].plot(z_vals, G_eff_vals / G_SI, 'g-', linewidth=2)
        axes[1,0].axvline(z_interface_1, color='gray', linestyle=':', alpha=0.5)
        axes[1,0].axvline(z_interface_2, color='gray', linestyle=':', alpha=0.5)
        axes[1,0].set_xlabel('z [m]')
        axes[1,0].set_ylabel('G_eff / G')
        axes[1,0].set_title('Effective Coupling Profile')
        axes[1,0].set_yscale('log')
        axes[1,0].grid(True, alpha=0.3)
        
        # Relative error
        axes[1,1].plot(z_vals, relative_error, 'm-', linewidth=1.5)
        axes[1,1].axvline(z_interface_1, color='gray', linestyle=':', alpha=0.5)
        axes[1,1].axvline(z_interface_2, color='gray', linestyle=':', alpha=0.5)
        axes[1,1].axhline(tolerance_rel, color='r', linestyle='--', alpha=0.5, label=f'Tolerance ({tolerance_rel})')
        axes[1,1].set_xlabel('z [m]')
        axes[1,1].set_ylabel('Relative Error')
        axes[1,1].set_title('Relative Error')
        axes[1,1].set_yscale('log')
        axes[1,1].legend()
        axes[1,1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        output_dir = Path('results')
        output_dir.mkdir(exist_ok=True)
        plt.savefig(output_dir / 'interface_validation.png', dpi=200)
        print(f"\n   Saved: results/interface_validation.png")
        plt.close()
    
    print("="*70 + "\n")


if __name__ == '__main__':
    success = test_1d_slab_interface()
    sys.exit(0 if success else 1)
