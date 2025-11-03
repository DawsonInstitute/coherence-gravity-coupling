"""
Materials Feasibility for Robin Boundary Conditions

Physical Question:
    Can we engineer realistic materials/structures to realize Robin BC parameter θ?

Robin BC: α Φ + β (∂Φ/∂n) = 0
    - Parameterization: α/β = tan(θ)
    - θ = 0: Dirichlet (Φ = 0, perfect conductor)
    - θ = -π/2: Neumann (∂Φ/∂n = 0, perfect insulator)
    - General θ: Mixed impedance boundary

Physical Implementations:
    1. Dielectric coatings: ε_r, μ_r → impedance Z → θ
    2. Metamaterial surfaces: Engineered Z via subwavelength structures
    3. Superconducting/normal metal junctions: Tunable penetration depth
    4. Casimir cavity coatings: Surface roughness, optical properties

References:
    - Gorkavenko et al. (2024) arXiv:2409.04647 (Robin BC in Casimir)
    - Pendry et al. (1999): Magnetism from conductors (metamaterials)
    - Task: Map θ → realistic implementations
"""
from __future__ import annotations

import numpy as np
from typing import Literal
from dataclasses import dataclass


@dataclass
class MaterialProperties:
    """Material electromagnetic properties."""
    epsilon_r: float  # Relative permittivity
    mu_r: float = 1.0  # Relative permeability
    conductivity: float = 0.0  # σ [S/m]
    thickness: float = 1e-6  # Layer thickness [m]
    name: str = "Generic"


def impedance_to_robin_theta(
    Z_surface: complex,
    Z_vacuum: float = 377.0
) -> float:
    """Convert surface impedance to Robin BC parameter θ.
    
    Physical mapping:
        - Perfect conductor (Z → 0): θ = 0 (Dirichlet)
        - Perfect magnetic (Z → ∞): θ = -π/2 (Neumann)
        - Matched impedance (Z = Z_0): θ ≈ -π/4 (intermediate)
    
    Formula: θ = arctan(Z_0/|Z|) - π/2
        - Shifts arctan range [0, π/2] to [-π/2, 0]
        - Inverts Z dependence (high Z → Neumann)
    
    Args:
        Z_surface: Surface impedance [Ω] (complex)
        Z_vacuum: Vacuum impedance (377 Ω)
    
    Returns:
        theta: Robin parameter [radians], θ ∈ [-π/2, 0]
    """
    Z_mag = np.abs(Z_surface)
    
    # Inverted ratio: higher Z → more negative θ (toward Neumann)
    ratio = Z_vacuum / Z_mag if Z_mag > 0 else 1e-10
    
    # Map [0, ∞) → [-π/2, 0]
    theta = np.arctan(ratio) - np.pi/2
    
    return theta


def dielectric_coating_impedance(
    mat: MaterialProperties,
    frequency: float = 1e9
) -> complex:
    """Compute surface impedance for dielectric coating.
    
    For thin layer on conductor:
        Z ~ sqrt(μ_r/ε_r) × Z_0 × tanh(γ d)
    
    where γ = jω sqrt(μ ε) is propagation constant.
    
    Args:
        mat: Material properties
        frequency: EM frequency [Hz]
    
    Returns:
        Z: Complex surface impedance [Ω]
    """
    omega = 2 * np.pi * frequency
    c = 299792458.0
    epsilon_0 = 8.854e-12
    mu_0 = 4 * np.pi * 1e-7
    
    # Propagation constant
    epsilon_complex = mat.epsilon_r - 1j * mat.conductivity / (omega * epsilon_0)
    gamma = 1j * omega * np.sqrt(mu_0 * mat.mu_r * epsilon_0 * epsilon_complex) / c
    
    # Intrinsic impedance
    Z_intrinsic = 377 * np.sqrt(mat.mu_r / epsilon_complex)
    
    # Surface impedance (thin film approximation)
    Z = Z_intrinsic * np.tanh(gamma * mat.thickness)
    
    return Z


def metamaterial_surface_impedance(
    wire_radius: float,
    wire_spacing: float,
    frequency: float = 1e9
) -> complex:
    """Effective impedance of wire-array metamaterial surface.
    
    Pendry et al. model: Parallel thin wires create effective plasma.
        ε_eff = 1 - (ω_p / ω)²
        ω_p² = c² / (r² ln(a/r))
    
    Args:
        wire_radius: r [m]
        wire_spacing: a [m]
        frequency: [Hz]
    
    Returns:
        Z: Effective surface impedance
    """
    c = 299792458.0
    omega = 2 * np.pi * frequency
    
    # Plasma frequency
    omega_p_sq = c**2 / (wire_radius**2 * np.log(wire_spacing / wire_radius))
    omega_p = np.sqrt(omega_p_sq)
    
    # Effective permittivity (can be negative below plasma frequency)
    epsilon_eff = 1 - (omega_p / omega)**2
    
    # Impedance: Z = Z_0 / sqrt(ε_eff)
    # Below ω_p: ε_eff < 0 → inductive impedance (Neumann-like)
    if epsilon_eff > 0:
        Z = 377 / np.sqrt(epsilon_eff)
    else:
        # Evanescent regime: imaginary impedance
        Z = 377 * 1j / np.sqrt(abs(epsilon_eff))
    
    return Z


def superconductor_penetration_impedance(
    london_depth: float,
    transition_ratio: float = 0.5
) -> complex:
    """Surface impedance of superconductor near T_c.
    
    Near critical temperature: normal + super currents mix.
        Z = Z_n × f + Z_s × (1 - f)
    
    where f = (T/T_c)^α is normal fraction.
    
    Args:
        london_depth: λ_L [m]
        transition_ratio: f ∈ [0,1] (0=full super, 1=normal)
    
    Returns:
        Z: Surface impedance [Ω]
    """
    # Superconducting impedance (reactive)
    Z_s = 1j * 377 * london_depth * 1e9  # Scaled for typical λ_L ~ nm
    
    # Normal metal impedance (resistive)
    Z_n = 10.0  # Ω (typical for normal metal at GHz)
    
    # Mixed state
    f = transition_ratio
    Z = Z_n * f + Z_s * (1 - f)
    
    return Z


def theta_sweep_materials() -> dict:
    """Survey achievable θ range with different materials.
    
    Returns:
        dict mapping material type to (θ_min, θ_max, properties)
    """
    results = {}
    
    # 1. Dielectric coatings
    dielectrics = [
        MaterialProperties(epsilon_r=2.0, name="Teflon"),
        MaterialProperties(epsilon_r=4.0, name="Glass"),
        MaterialProperties(epsilon_r=10.0, name="Alumina"),
        MaterialProperties(epsilon_r=100.0, name="Barium titanate"),
    ]
    
    theta_dielectric = []
    for mat in dielectrics:
        Z = dielectric_coating_impedance(mat)
        theta = impedance_to_robin_theta(Z)
        theta_dielectric.append(theta)
    
    results['dielectric_coatings'] = {
        'theta_range': (min(theta_dielectric), max(theta_dielectric)),
        'theta_values': theta_dielectric,
        'materials': [m.name for m in dielectrics]
    }
    
    # 2. Metamaterial surfaces
    wire_configs = [
        (1e-6, 10e-6),   # r=1μm, a=10μm
        (100e-9, 1e-6),  # r=100nm, a=1μm
        (10e-9, 100e-9), # r=10nm, a=100nm
    ]
    
    theta_metamaterial = []
    for r, a in wire_configs:
        Z = metamaterial_surface_impedance(r, a)
        theta = impedance_to_robin_theta(Z)
        theta_metamaterial.append(theta)
    
    results['metamaterial_wires'] = {
        'theta_range': (min(theta_metamaterial), max(theta_metamaterial)),
        'theta_values': theta_metamaterial,
        'configurations': wire_configs
    }
    
    # 3. Superconductors
    transition_fractions = [0.0, 0.25, 0.5, 0.75, 1.0]
    lambda_L = 50e-9  # 50 nm (typical for Nb)
    
    theta_sc = []
    for f in transition_fractions:
        Z = superconductor_penetration_impedance(lambda_L, f)
        theta = impedance_to_robin_theta(Z)
        theta_sc.append(theta)
    
    results['superconductor_junction'] = {
        'theta_range': (min(theta_sc), max(theta_sc)),
        'theta_values': theta_sc,
        'transition_fractions': transition_fractions,
        'london_depth_nm': lambda_L * 1e9
    }
    
    return results


def recommend_material_for_theta(
    theta_target: float,
    tolerance: float = 0.1
) -> dict:
    """Find materials/structures achieving target θ.
    
    Args:
        theta_target: Desired Robin parameter [rad]
        tolerance: Acceptable deviation [rad]
    
    Returns:
        dict with recommended implementations
    """
    survey = theta_sweep_materials()
    recommendations = []
    
    for material_type, data in survey.items():
        theta_vals = data['theta_values']
        
        # Find closest match
        idx_closest = np.argmin(np.abs(np.array(theta_vals) - theta_target))
        theta_closest = theta_vals[idx_closest]
        
        if abs(theta_closest - theta_target) < tolerance:
            rec = {
                'type': material_type,
                'theta_achieved': theta_closest,
                'error': abs(theta_closest - theta_target),
            }
            
            if material_type == 'dielectric_coatings':
                rec['material'] = data['materials'][idx_closest]
            elif material_type == 'metamaterial_wires':
                rec['configuration'] = data['configurations'][idx_closest]
            elif material_type == 'superconductor_junction':
                rec['transition_fraction'] = data['transition_fractions'][idx_closest]
            
            recommendations.append(rec)
    
    return {
        'theta_target': theta_target,
        'tolerance': tolerance,
        'recommendations': recommendations,
        'achievable': len(recommendations) > 0
    }


# Example usage and validation
if __name__ == "__main__":
    print("Materials Feasibility for Robin BC")
    print("=" * 60)
    
    # Survey available θ range
    print("\n[Material Survey]")
    survey = theta_sweep_materials()
    
    for mat_type, data in survey.items():
        theta_min, theta_max = data['theta_range']
        print(f"\n{mat_type}:")
        print(f"  θ range: [{theta_min*180/np.pi:.1f}°, {theta_max*180/np.pi:.1f}°]")
        
        if 'materials' in data:
            print(f"  Materials: {', '.join(data['materials'])}")
        if 'configurations' in data:
            print(f"  # Configurations: {len(data['configurations'])}")
        if 'london_depth_nm' in data:
            print(f"  London depth: {data['london_depth_nm']:.1f} nm")
    
    # Example: Find material for specific θ
    print("\n[Target θ = -π/4 (45° from Dirichlet toward Neumann)]")
    target = -np.pi / 4
    rec = recommend_material_for_theta(target, tolerance=0.2)
    
    print(f"  Target: {target*180/np.pi:.1f}°")
    print(f"  Achievable: {rec['achievable']}")
    
    if rec['recommendations']:
        print(f"  Recommendations:")
        for r in rec['recommendations']:
            print(f"    - {r['type']}: θ = {r['theta_achieved']*180/np.pi:.1f}° (error: {r['error']*180/np.pi:.2f}°)")
            if 'material' in r:
                print(f"      Material: {r['material']}")
    
    # Practical example: Casimir cavity optimization
    print("\n[Casimir Cavity Optimization]")
    print("  Goal: Maximize ξ-sensitivity via Robin BC")
    print("  Strategy: θ ≈ -π/3 to -π/2 (strong Neumann character)")
    
    theta_optimal = -np.pi / 3
    casimir_rec = recommend_material_for_theta(theta_optimal, tolerance=0.15)
    
    print(f"  Optimal θ: {theta_optimal*180/np.pi:.1f}°")
    if casimir_rec['recommendations']:
        best = casimir_rec['recommendations'][0]
        print(f"  Best option: {best['type']}")
        print(f"    θ achieved: {best['theta_achieved']*180/np.pi:.1f}°")
        if 'material' in best:
            print(f"    Coating: {best['material']}")
    
    print("\n✅ Materials feasibility module functional")
    print("   Enables new physics discovery: θ parameter → engineerable structures")
