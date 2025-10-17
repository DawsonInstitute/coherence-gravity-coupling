"""
Example: Point mass with coherent shell.

Demonstrates how a coherent shell around a mass can create different
effective gravitational coupling inside vs outside.

This is a proof-of-concept for spatially-varying G_eff.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from field_equations.action import CoherenceGravityParams
from solvers.static_spherical import StaticSphericalSolver


def plot_coherent_shell_potential():
    """
    Plot gravitational potential for mass with coherent shell.
    
    Shows how coherence transition creates modified potential profile.
    """
    params = CoherenceGravityParams(xi=10.0)
    solver = StaticSphericalSolver(params)
    
    # Setup
    M_core = 1e3  # kg
    R_shell = 1.0  # m (shell radius)
    
    # Coherence values
    Phi0_inside = 1e15  # Strong coherence inside shell
    Phi0_outside = 0.0  # No coherence outside
    
    r = np.linspace(0.1, 5.0, 200)  # m
    
    # Solve
    result = solver.coherent_shell_example(
        r, M_core, R_shell, Phi0_inside, Phi0_outside
    )
    
    # Standard GR for comparison
    Phi_standard = -params.G * M_core / r
    
    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # Potential
    ax1.plot(r, result['potential'], 'b-', linewidth=2, label='With coherent shell')
    ax1.plot(r, Phi_standard, 'k--', linewidth=2, label='Standard GR')
    ax1.axvline(R_shell, color='r', linestyle=':', label=f'Shell at r={R_shell}m')
    ax1.set_xlabel('Radius r [m]', fontsize=12)
    ax1.set_ylabel('Gravitational Potential Φ [J/kg]', fontsize=12)
    ax1.set_title('Point Mass with Coherent Shell', fontsize=14, fontweight='bold')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Potential ratio
    ratio = result['potential'] / Phi_standard
    ax2.plot(r, ratio, 'g-', linewidth=2)
    ax2.axvline(R_shell, color='r', linestyle=':', label=f'Shell at r={R_shell}m')
    ax2.axhline(1.0, color='k', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Radius r [m]', fontsize=12)
    ax2.set_ylabel('Φ_coherent / Φ_standard', fontsize=12)
    ax2.set_title('Potential Suppression Factor', fontsize=14, fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('coherent_shell_potential.png', dpi=150)
    print("✅ Saved plot to coherent_shell_potential.png")
    
    # Print summary
    print("\n" + "="*70)
    print("COHERENT SHELL ANALYSIS")
    print("="*70)
    print(f"Core mass: {M_core:.1e} kg")
    print(f"Shell radius: {R_shell:.2f} m")
    print(f"Coherence inside: Φ₀ = {Phi0_inside:.1e} m⁻¹")
    print(f"Coherence outside: Φ₀ = {Phi0_outside:.1e} m⁻¹")
    print()
    print(f"G_eff inside: {result['G_int']:.3e} m³/(kg·s²)")
    print(f"G_eff outside: {result['G_ext']:.3e} m³/(kg·s²)")
    print(f"Suppression ratio: G_in/G_out = {result['suppression_ratio']:.3e}")
    print()
    print(f"At shell boundary (r={R_shell}m):")
    idx = np.argmin(np.abs(r - R_shell))
    print(f"  Φ_coherent = {result['potential'][idx]:.3e} J/kg")
    print(f"  Φ_standard = {Phi_standard[idx]:.3e} J/kg")
    print(f"  Suppression = {ratio[idx]:.3e}×")
    print("="*70)


def compare_coherence_levels():
    """
    Compare potentials for different coherence amplitudes.
    """
    params = CoherenceGravityParams(xi=1.0)
    solver = StaticSphericalSolver(params)
    
    M = 1e3  # kg
    r = np.linspace(0.1, 10.0, 200)
    
    # Different coherence levels
    Phi0_values = [0.0, 1e10, 1e12, 1e14, 1e15]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for Phi0 in Phi0_values:
        result = solver.point_mass_analytic(r, M, Phi0)
        suppression = result['suppression']
        
        label = f'Φ₀ = {Phi0:.0e}, G_eff/G = {suppression:.2e}'
        ax.plot(r, result['potential'], linewidth=2, label=label)
    
    ax.set_xlabel('Radius r [m]', fontsize=12)
    ax.set_ylabel('Gravitational Potential [J/kg]', fontsize=12)
    ax.set_title('Point Mass: Effect of Coherence Amplitude', fontsize=14, fontweight='bold')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('coherence_levels_comparison.png', dpi=150)
    print("✅ Saved plot to coherence_levels_comparison.png")


def test_force_profile():
    """
    Test gravitational force profile with coherence.
    """
    params = CoherenceGravityParams(xi=10.0)
    solver = StaticSphericalSolver(params)
    
    M = 1e3
    R = 1.0
    Phi0 = 1e15
    
    r = np.linspace(0.1, 5.0, 200)
    
    result = solver.uniform_sphere_analytic(r, M, R, Phi0)
    
    # Standard GR force
    F_std = np.where(r < R, 
                     params.G * M * r / R**3,
                     params.G * M / r**2)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.plot(r, result['force'], 'b-', linewidth=2, label='With coherence')
    ax.plot(r, F_std, 'k--', linewidth=2, label='Standard GR')
    ax.axvline(R, color='r', linestyle=':', linewidth=2, label=f'Surface r={R}m')
    
    ax.set_xlabel('Radius r [m]', fontsize=12)
    ax.set_ylabel('Gravitational Force [N/kg]', fontsize=12)
    ax.set_title('Force Profile: Uniform Sphere with Coherence', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('force_profile.png', dpi=150)
    print("✅ Saved plot to force_profile.png")
    
    print(f"\nG_eff/G = {result['G_eff']/params.G:.3e}")
    print(f"Force suppression: {result['G_eff']/params.G:.3e}×")


if __name__ == "__main__":
    print("="*70)
    print("COHERENCE-MODULATED GRAVITY: SPHERICAL EXAMPLES")
    print("="*70)
    print()
    
    print("1. Plotting coherent shell potential...")
    plot_coherent_shell_potential()
    print()
    
    print("2. Comparing coherence levels...")
    compare_coherence_levels()
    print()
    
    print("3. Testing force profile...")
    test_force_profile()
    print()
    
    print("="*70)
    print("✅ ALL EXAMPLES COMPLETE")
    print("="*70)
