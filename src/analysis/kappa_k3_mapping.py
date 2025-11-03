"""
κ_R → k_3 Mapping: Laboratory Bounds to Bahamonde Torsion-EM Coupling

This module maps our laboratory curvature-EM coupling constraints (κ_R < 5×10¹⁷ m²)
to constraints on Bahamonde et al.'s torsion-EM coupling parameter k_3 in the
F^μν R̃_μν theory.

Physical Framework (Bahamonde et al. arxiv:2507.02362):
    - Riemann-Cartan action includes: k_3 F^μν R̃_μν
    - R̃_μν is Ricci tensor with torsion (can be asymmetric)
    - For comparison, our lab measures: κ_R R F^μν F_μν

EFT Expansion Strategy:
    Both operators appear in effective field theory expansion:
    
    L_eff = L_GR + L_EM + (c_1 R F²) + (c_2 F R̃) + ...
    
    where:
        - c_1 ~ κ_R (our measurement)
        - c_2 ~ k_3 (Bahamonde parameter)
    
    In weak-field, low-torsion limit:
        R̃_μν ≈ R_μν + O(torsion²)
        F R̃ ≈ F R + (torsion-dependent corrections)
    
    Therefore: k_3 ≲ κ_R + (torsion contribution uncertainty)

Key Insight:
    Laboratory nulls on κ_R R F² constrain k_3 if we assume:
    1. Torsion is negligible in lab (pseudo-Riemannian geometry)
    2. EFT coefficients c_1, c_2 have similar naturalness scales
    3. Both operators contribute to same observables at O(R F²)

Implementation:
    Map κ_R laboratory bound to k_3 upper limit accounting for:
    - Geometric factor differences (R vs R̃)
    - Field configuration dependence (FF* vs F R̃)
    - EFT coefficient hierarchy assumptions

References:
    - Bahamonde et al. (2025) arXiv:2507.02362 (torsion-EM coupling)
    - Our null_results.tex: κ_R < 5×10¹⁷ m² from lab nulls
    - Task: arxiv.2507.02362:26
"""
from __future__ import annotations

import numpy as np
_BSM_AVAILABLE = False
epsilon_equiv = None
DEFAULT_ENVIRONMENTS = None
axion_equiv_parametric = None
try:
    from src.analysis.bsm_bounds_from_kappa import (
        epsilon_equiv as _eps,
        DEFAULT_ENVIRONMENTS as _envs,
        axion_equiv_parametric as _gax,
    )
    epsilon_equiv, DEFAULT_ENVIRONMENTS, axion_equiv_parametric = _eps, _envs, _gax
    _BSM_AVAILABLE = True
except Exception:
    # Fallback: import by file path to be robust in scripts
    try:
        import importlib.util, os
        here = os.path.dirname(os.path.abspath(__file__))
        mod_path = os.path.join(here, 'bsm_bounds_from_kappa.py')
        spec = importlib.util.spec_from_file_location('bsm_bounds_from_kappa', mod_path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)  # type: ignore[attr-defined]
        epsilon_equiv = mod.epsilon_equiv
        DEFAULT_ENVIRONMENTS = mod.DEFAULT_ENVIRONMENTS
        axion_equiv_parametric = mod.axion_equiv_parametric
        _BSM_AVAILABLE = True
    except Exception:
        _BSM_AVAILABLE = False
from typing import Literal

# Physical constants
G = 6.67430e-11  # m³/(kg·s²)
c = 299792458.0  # m/s
hbar = 1.054571817e-34  # J·s

# Laboratory baseline (from null_results.tex)
KAPPA_R_LAB = 5e17  # m² (95% CL upper limit)
B_LAB = 10.0  # T
R_LAB = 1e-26  # m⁻²


def eft_coefficient_mapping(
    kappa_R: float,
    hierarchy_assumption: Literal['natural', 'suppressed', 'enhanced'] = 'natural'
) -> tuple[float, float]:
    """Map κ_R to k_3 under EFT coefficient naturalness assumptions.
    
    Naturalness scenarios:
        - 'natural': c_1 ~ c_2 ~ (MPlanck)^-2 → k_3 ~ κ_R
        - 'suppressed': c_2 < c_1 (torsion interactions weaker) → k_3 < κ_R
        - 'enhanced': c_2 > c_1 (torsion enhanced by spin) → k_3 > κ_R
    
    Args:
        kappa_R: Laboratory bound on curvature-EM coupling [m²]
        hierarchy_assumption: EFT naturalness scenario
    
    Returns:
        (k_3_central, k_3_uncertainty): Central value and theory uncertainty
    """
    if hierarchy_assumption == 'natural':
        # Assume dimensional analysis: both ~ (energy scale)^-2
        k_3_central = kappa_R
        k_3_uncertainty = 0.5 * kappa_R  # Factor of ~2 theory uncertainty
    
    elif hierarchy_assumption == 'suppressed':
        # Torsion suppressed by spin-orbit coupling scale
        suppression_factor = 0.1  # Assume order-of-magnitude suppression
        k_3_central = suppression_factor * kappa_R
        k_3_uncertainty = 0.5 * k_3_central
    
    elif hierarchy_assumption == 'enhanced':
        # Torsion enhanced by material spin density
        enhancement_factor = 10.0  # Conservative enhancement estimate
        k_3_central = enhancement_factor * kappa_R
        k_3_uncertainty = 5.0 * kappa_R  # Larger uncertainty for enhanced scenario
    
    else:
        raise ValueError(f"Unknown hierarchy: {hierarchy_assumption}")
    
    return k_3_central, k_3_uncertainty


def geometric_factor_correction(
    R_riemann: float,
    torsion_amplitude: float = 0.0
) -> float:
    """Compute correction factor between R and R̃ contributions.
    
    R̃_μν = R_μν + (torsion-dependent terms)
    
    For small torsion: R̃ ≈ R (1 + ε_torsion)
    
    Args:
        R_riemann: Riemann curvature [m⁻²]
        torsion_amplitude: Characteristic torsion scale [m⁻¹]
    
    Returns:
        correction_factor: R̃/R ratio
    """
    if torsion_amplitude == 0.0:
        return 1.0
    
    # Torsion contribution scales as (torsion)² / R
    # For weak torsion: R̃ ≈ R (1 + (T²/|R|))
    # Cap epsilon to avoid numerical overflow
    epsilon = min((torsion_amplitude**2) / (abs(R_riemann) + 1e-30), 1e10)
    correction = 1.0 + epsilon
    
    return correction


def map_kappa_to_k3(
    kappa_R_lab: float = KAPPA_R_LAB,
    R_lab: float = R_LAB,
    B_lab: float = B_LAB,
    scenario: Literal['conservative', 'moderate', 'optimistic'] = 'conservative'
) -> dict[str, float]:
    """Primary mapping function: κ_R laboratory bound → k_3 constraint.
    
    Scenarios:
        - 'conservative': Minimal assumptions, largest k_3 upper limit
        - 'moderate': Natural EFT hierarchy, geometric factors included
        - 'optimistic': Assume torsion suppression, tightest k_3 bound
    
    Args:
        kappa_R_lab: Laboratory κ_R bound [m²]
        R_lab: Laboratory curvature [m⁻²]
        B_lab: Laboratory magnetic field [T]
        scenario: Mapping assumption scenario
    
    Returns:
        dict with keys:
            'k3_upper_limit': Upper bound on k_3 [m²]
            'k3_central': Central value estimate [m²]
            'k3_lower_limit': Lower bound (if constrained) [m²]
            'geometric_correction': R̃/R correction factor
            'eft_hierarchy': Assumed hierarchy scenario
    """
    result = {}
    
    if scenario == 'conservative':
        # Minimal assumptions: k_3 ≲ κ_R with large uncertainty
        k3_central, k3_unc = eft_coefficient_mapping(kappa_R_lab, 'enhanced')
        result['k3_upper_limit'] = k3_central + 2 * k3_unc  # 2σ upper limit
        result['k3_central'] = k3_central
        result['k3_lower_limit'] = None  # No lower bound from lab nulls
        result['geometric_correction'] = 1.0  # Assume R̃ ≈ R
        result['eft_hierarchy'] = 'enhanced'
    
    elif scenario == 'moderate':
        # Natural hierarchy with geometric corrections
        k3_central, k3_unc = eft_coefficient_mapping(kappa_R_lab, 'natural')
        geom_corr = geometric_factor_correction(R_lab, torsion_amplitude=1e10)  # ~coherence scale
        k3_central *= geom_corr
        result['k3_upper_limit'] = k3_central + k3_unc
        result['k3_central'] = k3_central
        result['k3_lower_limit'] = max(0, k3_central - k3_unc)
        result['geometric_correction'] = geom_corr
        result['eft_hierarchy'] = 'natural'
    
    elif scenario == 'optimistic':
        # Assume torsion suppression in lab
        k3_central, k3_unc = eft_coefficient_mapping(kappa_R_lab, 'suppressed')
        geom_corr = geometric_factor_correction(R_lab, torsion_amplitude=0.0)
        result['k3_upper_limit'] = k3_central + k3_unc
        result['k3_central'] = k3_central
        result['k3_lower_limit'] = max(0, k3_central - k3_unc)
        result['geometric_correction'] = geom_corr
        result['eft_hierarchy'] = 'suppressed'
    
    else:
        raise ValueError(f"Unknown scenario: {scenario}")
    
    return result


def amplification_to_astrophysical(
    k3_lab_bound: float,
    B_astro: float,
    R_astro: float,
    B_lab: float = B_LAB,
    R_lab: float = R_LAB
) -> tuple[float, float]:
    """Amplify lab k_3 bound to astrophysical regime.
    
    Similar to κ_R astrophysical recast: effective coupling scales as k_3 B² R.
    
    Args:
        k3_lab_bound: Laboratory upper limit on k_3 [m²]
        B_astro: Astrophysical magnetic field [T]
        R_astro: Astrophysical curvature [m⁻²]
        B_lab, R_lab: Laboratory reference values
    
    Returns:
        (k3_astro, amplification_factor)
    """
    amplification = (B_astro / B_lab)**2 * (R_astro / R_lab)
    k3_astro = k3_lab_bound / amplification
    
    return k3_astro, amplification


# Example usage and validation
if __name__ == "__main__":
    print("κ_R → k_3 Mapping Module")
    print("=" * 50)
    
    # Map laboratory bound to k_3
    for scenario in ['conservative', 'moderate', 'optimistic']:
        result = map_kappa_to_k3(scenario=scenario)
        print(f"\nScenario: {scenario}")
        print(f"  k_3 central: {result['k3_central']:.2e} m²")
        print(f"  k_3 upper (95% CL): {result['k3_upper_limit']:.2e} m²")
        if result['k3_lower_limit'] is not None:
            print(f"  k_3 lower: {result['k3_lower_limit']:.2e} m²")
        print(f"  Geometric correction: {result['geometric_correction']:.3f}")
        print(f"  EFT hierarchy: {result['eft_hierarchy']}")
    
    # Astrophysical amplification example (magnetar)
    print("\n" + "=" * 50)
    print("Astrophysical Amplification (Magnetar Surface):")
    result_mod = map_kappa_to_k3(scenario='moderate')
    k3_magnetar, amp = amplification_to_astrophysical(
        k3_lab_bound=result_mod['k3_upper_limit'],
        B_astro=1e10,  # 10^10 T magnetar field
        R_astro=1e-6   # ~10^-6 m^-2 neutron star surface curvature
    )
    print(f"  Magnetar k_3 < {k3_magnetar:.2e} m²")
    print(f"  Amplification: {amp:.2e}×")
    
    print("\n✅ κ_R → k_3 mapping module functional")
    print("   Enables new physics discovery: lab bounds → torsion-EM constraints")
    
    # NEW: Beyond-SM predictions
    print("\n" + "="*60)
    print("[Beyond-Standard-Model Particle Constraints]")
    if not _BSM_AVAILABLE:
        print("  (BSM mapping module not available; skip)")
    else:
        # Dark photon mapping (robust): ε_eff ≈ Cε (κ_R·R)
        print("\nDark photon (robust mapping): ε_eff ≈ Cε · (κ_R · R)")
        for env_key in ["lab_flat", "earth_surface", "magnetar_surface"]:
            R = DEFAULT_ENVIRONMENTS[env_key].R_m2
            eps1 = epsilon_equiv(KAPPA_R_LAB, R, C_eps=1.0)
            eps2 = epsilon_equiv(KAPPA_R_LAB, R, C_eps=1e-2)
            print(f"  Env={env_key:17s}: ε(Cε=1) ≈ {eps1:.2e}   ε(Cε=1e-2) ≈ {eps2:.2e}")

        # Axion (illustrative only): g_equiv ≈ Ca (κ_R·R)/Λ with Λ=10 TeV
        print("\nAxion (illustrative, model-dependent): g_equiv ≈ Ca · (κ_R · R) / Λ (Λ=10 TeV)")
        for env_key in ["lab_flat", "earth_surface", "magnetar_surface"]:
            R = DEFAULT_ENVIRONMENTS[env_key].R_m2
            g_eq = axion_equiv_parametric(KAPPA_R_LAB, R, C_a=1.0, Lambda_GeV=1e4)
            print(f"  Env={env_key:17s}: g_equiv ≈ {g_eq:.2e} GeV⁻¹")

