"""
Astrophysical Recast Helper: Map Laboratory κ_R Bounds to Compact-Object Regimes

This module provides utilities to propagate our laboratory curvature–EM coupling bounds
(κ_R < 5×10¹⁷ m² from null_results.tex) to astrophysical environments with strong fields
and large curvature (magnetars, black hole horizons).

Usage:
    from utils.astrophysical_recast import recast_to_magnetar, recast_to_black_hole
    
    # Magnetar surface
    kappa_astro, amplification = recast_to_magnetar(B_surface=1e10)
    
    # Black hole horizon
    kappa_bh, amp_bh = recast_to_black_hole(M_solar_masses=10, B_tesla=1e8)

Physical Scaling:
    The effective coupling strength S = κ_R × B² × R appears in observables.
    If lab measures S_lab as upper limit, then astrophysical environments with
    (B_astro, R_astro) constrain:
        κ_R < S_lab / (B_astro² × R_astro)
    
    Amplification factor: A = (B_astro/B_lab)² × (R_astro/R_lab)
    Improved bound: κ_R_astro = κ_R_lab / A
"""
from __future__ import annotations
import numpy as np

# Physical constants (SI)
G = 6.67430e-11  # m³/(kg·s²)
c = 299792458.0  # m/s
M_sun = 1.98847e30  # kg

# Laboratory baseline from null_results.tex
KAPPA_R_LAB = 5e17  # m² (95% CL upper limit)
B_LAB = 10.0  # T (representative field strength)
R_LAB = 1e-26  # m⁻² (terrestrial Ricci curvature)

S_LAB = KAPPA_R_LAB * B_LAB**2 * R_LAB  # Dimensionless effective coupling


def horizon_curvature(M_kg: float) -> float:
    """Compute Ricci curvature at Schwarzschild horizon.
    
    For Schwarzschild metric, R ~ r_s^-2 where r_s = 2GM/c².
    
    Args:
        M_kg: Black hole mass in kilograms
    
    Returns:
        R: Ricci curvature scalar at horizon [m⁻²]
    """
    r_s = 2 * G * M_kg / c**2
    return 1 / r_s**2


def recast_to_magnetar(
    B_surface: float = 1e10,
    R_surface: float | None = None,
    M_ns: float = 1.4 * M_sun,
    R_ns: float = 12e3
) -> tuple[float, float]:
    """Propagate lab κ_R bound to magnetar environment.
    
    Magnetar parameters:
        - Surface field: B ~ 10⁸–10¹⁰ T (confirmed via cyclotron lines)
        - Interior field: B ~ 10¹⁰–10¹² T (theoretical)
        - Radius: R_ns ~ 10–12 km
        - Mass: M_ns ~ 1.4–2.0 M_sun
    
    Args:
        B_surface: Magnetic field at surface [T]
        R_surface: Ricci curvature at surface [m⁻²]. If None, estimate from R_ns.
        M_ns: Neutron star mass [kg]
        R_ns: Neutron star radius [m]
    
    Returns:
        (kappa_R_constrained, amplification_factor)
    """
    if R_surface is None:
        # Rough estimate: R ~ GM/(R_ns² c²) for surface curvature
        R_surface = G * M_ns / (R_ns**2 * c**2)
    
    amplification = (B_surface / B_LAB)**2 * (R_surface / R_LAB)
    kappa_R_constrained = KAPPA_R_LAB / amplification
    
    return kappa_R_constrained, amplification


def recast_to_black_hole(
    M_solar_masses: float = 10.0,
    B_tesla: float = 1e8
) -> tuple[float, float]:
    """Propagate lab κ_R bound to black hole environment.
    
    Typical scenarios:
        - Stellar-mass BH: M ~ 3–20 M_sun
        - Supermassive BH: M ~ 10⁶–10¹⁰ M_sun
        - Magnetized BH candidates: B ~ 10⁴–10¹⁰ T (theoretical/indirect)
    
    Args:
        M_solar_masses: Black hole mass in solar masses
        B_tesla: Magnetic field strength near horizon [T]
    
    Returns:
        (kappa_R_constrained, amplification_factor)
    """
    M_kg = M_solar_masses * M_sun
    R_horizon = horizon_curvature(M_kg)
    
    amplification = (B_tesla / B_LAB)**2 * (R_horizon / R_LAB)
    kappa_R_constrained = KAPPA_R_LAB / amplification
    
    return kappa_R_constrained, amplification


def qnm_frequency_shift(
    kappa_R: float,
    B: float,
    R: float,
    baseline_omega: float = 1e3
) -> float:
    """Estimate fractional QNM frequency shift from κ_R coupling.
    
    Following Karimabadi et al. (arXiv:2508.13820), assume linear scaling:
        Δω/ω ∝ κ_R × B² × R
    
    Args:
        kappa_R: Coupling parameter [m²]
        B: Magnetic field [T]
        R: Ricci curvature [m⁻²]
        baseline_omega: Characteristic QNM frequency [Hz]
    
    Returns:
        |Δω/ω|: Fractional frequency shift (dimensionless)
    """
    # Dimensionless coupling strength
    S = kappa_R * B**2 * R
    # Assume O(1) proportionality constant for order-of-magnitude estimate
    return abs(S)


def print_scenario_summary(
    name: str,
    M_solar: float | None = None,
    B: float = 1e10,
    R: float | None = None,
    is_magnetar: bool = False
):
    """Print formatted summary for astrophysical scenario."""
    print(f"\n{'='*70}")
    print(f"SCENARIO: {name}")
    print(f"{'='*70}")
    
    if is_magnetar:
        kappa, amp = recast_to_magnetar(B_surface=B, R_surface=R)
        if R is None:
            R = G * 1.4 * M_sun / (12e3**2 * c**2)  # default NS curvature
    else:
        if M_solar is None:
            raise ValueError("Must specify M_solar for BH scenario")
        kappa, amp = recast_to_black_hole(M_solar_masses=M_solar, B_tesla=B)
        R = horizon_curvature(M_solar * M_sun)
    
    shift = qnm_frequency_shift(kappa, B, R)
    
    print(f"  Field: B = {B:.2e} T")
    if M_solar is not None:
        print(f"  Mass: M = {M_solar:.1f} M_sun")
    print(f"  Curvature: R = {R:.2e} m⁻²")
    print(f"  Amplification vs. lab: {amp:.2e}")
    print(f"  Constrained κ_R: {kappa:.2e} m²")
    print(f"  Improvement over lab: {KAPPA_R_LAB / kappa:.2e}×")
    print(f"  Predicted |Δω/ω|: {shift:.2e}")
    
    # Detectability assessment
    if shift > 1e-2:
        status = "DETECTABLE (conservative LIGO)"
    elif shift > 1e-3:
        status = "Marginally detectable (optimistic LIGO)"
    else:
        status = "Below current GW sensitivity"
    print(f"  Status: {status}")
    print(f"{'='*70}")


if __name__ == "__main__":
    print("=" * 70)
    print("Astrophysical Recast: Laboratory κ_R Bounds")
    print("=" * 70)
    print(f"\nLaboratory baseline:")
    print(f"  κ_R < {KAPPA_R_LAB:.2e} m² (95% CL)")
    print(f"  B_lab = {B_LAB} T")
    print(f"  R_lab = {R_LAB:.2e} m⁻²")
    print(f"  S_lab = {S_LAB:.2e} (effective coupling)")
    
    # Scenario A: Weak-field stellar BH
    print_scenario_summary(
        "Stellar BH + Weak Field",
        M_solar=10,
        B=1e8,
        is_magnetar=False
    )
    
    # Scenario B: Magnetar-strength low-mass BH
    print_scenario_summary(
        "Low-Mass BH + Magnetar Field",
        M_solar=3,
        B=1e10,
        is_magnetar=False
    )
    
    # Scenario C: Magnetar surface
    print_scenario_summary(
        "Magnetar Surface (canonical)",
        B=1e10,
        is_magnetar=True
    )
    
    # Scenario D: Hypermagnetic compact object
    print_scenario_summary(
        "Magnetar Interior (extreme)",
        B=1e12,
        R=1e5,  # compressed NS core
        is_magnetar=True
    )
    
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("Laboratory nulls remain relevant for weak-field environments.")
    print("Magnetar-BH systems (B ~ 10¹⁰ T, M ~ 1–10 M_sun) offer")
    print("amplification factors ~10¹⁵–10²⁵, potentially detectable in")
    print("next-generation GW observatories (Einstein Telescope, Cosmic Explorer).")
    print("="*70)
