"""
Physical calibration of coherence field Φ to experimental systems.

Maps observable quantities (densities, correlation lengths, healing lengths)
to the coherence field amplitude Φ [m⁻¹] used in the theory.

This is critical for grounding the (ξ, Φ₀) parameter space in reality.
"""

import numpy as np
from typing import Dict, Optional
from dataclasses import dataclass


# Physical constants
HBAR = 1.054571817e-34  # J·s
M_E = 9.1093837015e-31  # kg (electron mass)
K_B = 1.380649e-23      # J/K
C = 299792458.0         # m/s


@dataclass
class SystemCalibration:
    """
    Calibration result mapping physical system to Φ₀.
    
    Attributes:
        system_name: Human-readable name
        Phi0: Coherence field amplitude [m⁻¹]
        uncertainty_factor: Multiplicative uncertainty (e.g., 2.0 means 2× error bars)
        basis: Physical quantity used for calibration
        notes: Additional context
    """
    system_name: str
    Phi0: float
    uncertainty_factor: float
    basis: str
    notes: str = ""


class BECCalibration:
    """
    Calibrate coherence field from Bose-Einstein Condensate parameters.
    
    Physical basis:
    - Healing length ξ_h = 1/√(8πna_s) where n is density, a_s is scattering length
    - Coherence length ℓ_c ~ ξ_h in homogeneous BEC
    - Order parameter |ψ|² = n
    
    Proposed mapping: Φ ≈ 1/ℓ_c (inverse correlation length)
    Alternative: Φ ~ √n (number density scaling)
    """
    
    @staticmethod
    def from_density_scattering(n: float, a_s: float, 
                                 method: str = "healing") -> SystemCalibration:
        """
        Calibrate Φ from BEC density and scattering length.
        
        Args:
            n: Number density [m⁻³]
            a_s: s-wave scattering length [m]
            method: "healing" (Φ = 1/ξ_h) or "sqrt_density" (Φ ~ √n)
            
        Returns:
            SystemCalibration with Φ₀ estimate
        """
        if method == "healing":
            # Healing length ξ_h = 1/√(8πna_s)
            xi_h = 1.0 / np.sqrt(8 * np.pi * n * a_s)
            Phi0 = 1.0 / xi_h
            basis = f"ξ_h = {xi_h*1e9:.2f} nm from n={n:.2e} m⁻³, a_s={a_s*1e9:.2f} nm"
            uncertainty = 2.0  # Factor of 2 uncertainty in mapping
            
        elif method == "sqrt_density":
            # Dimensional analysis: Φ ~ √n
            Phi0 = np.sqrt(n)
            basis = f"√n scaling from n={n:.2e} m⁻³"
            uncertainty = 3.0  # Higher uncertainty for this ad-hoc scaling
            
        else:
            raise ValueError(f"Unknown method: {method}")
        
        return SystemCalibration(
            system_name=f"BEC (generic, {method})",
            Phi0=Phi0,
            uncertainty_factor=uncertainty,
            basis=basis,
            notes=f"Typical lab BEC, method={method}"
        )
    
    @staticmethod
    def rubidium_87_typical() -> SystemCalibration:
        """
        Typical ⁸⁷Rb BEC parameters.
        
        Standard parameters:
        - n ~ 10¹⁴ cm⁻³ = 10²⁰ m⁻³
        - a_s = 5.3 nm (5.3e-9 m)
        - T ~ 100 nK
        """
        n = 1e20  # m⁻³ (10¹⁴ cm⁻³)
        a_s = 5.3e-9  # m
        
        result = BECCalibration.from_density_scattering(n, a_s, "healing")
        result.system_name = "⁸⁷Rb BEC (typical)"
        result.notes = "Standard dilute ultracold atomic gas; n=10²⁰ m⁻³, a_s=5.3 nm"
        
        return result
    
    @staticmethod
    def sodium_23_typical() -> SystemCalibration:
        """
        Typical ²³Na BEC parameters.
        
        - n ~ 10¹⁴ cm⁻³
        - a_s = 2.8 nm
        """
        n = 1e20  # m⁻³
        a_s = 2.8e-9  # m
        
        result = BECCalibration.from_density_scattering(n, a_s, "healing")
        result.system_name = "²³Na BEC (typical)"
        result.notes = "Dilute sodium BEC; n=10²⁰ m⁻³, a_s=2.8 nm"
        
        return result
    
    @staticmethod
    def high_density_bec() -> SystemCalibration:
        """
        High-density BEC (optimistic case).
        
        - n ~ 10²² m⁻³ (pushing limits of dilute-gas approximation)
        - a_s ~ 5 nm
        """
        n = 1e22  # m⁻³
        a_s = 5e-9  # m
        
        result = BECCalibration.from_density_scattering(n, a_s, "healing")
        result.system_name = "BEC (high density)"
        result.notes = "Near limit of dilute-gas regime; n=10²² m⁻³"
        result.uncertainty_factor = 3.0  # Higher uncertainty at boundary
        
        return result


class SuperconductorCalibration:
    """
    Calibrate coherence field from superconductor parameters.
    
    Physical basis:
    - Cooper pair density n_s
    - Coherence length ξ_SC (Pippard/BCS coherence length)
    - London penetration depth λ_L
    
    Proposed mapping: Φ ≈ 1/ξ_SC
    """
    
    @staticmethod
    def from_coherence_length(xi_SC: float, 
                              material: str = "generic") -> SystemCalibration:
        """
        Calibrate from superconductor coherence length.
        
        Args:
            xi_SC: Coherence length [m]
            material: Material name for labeling
            
        Returns:
            SystemCalibration
        """
        Phi0 = 1.0 / xi_SC
        
        return SystemCalibration(
            system_name=f"Superconductor ({material})",
            Phi0=Phi0,
            uncertainty_factor=2.0,
            basis=f"ξ_SC = {xi_SC*1e9:.1f} nm",
            notes=f"Coherence length for {material}"
        )
    
    @staticmethod
    def aluminum_thin_film() -> SystemCalibration:
        """
        Aluminum superconducting thin film.
        
        - T_c = 1.2 K
        - ξ_SC ~ 1600 nm (at T << T_c)
        """
        xi_SC = 1600e-9  # m
        result = SuperconductorCalibration.from_coherence_length(xi_SC, "Al")
        result.notes = "Aluminum thin film, T_c=1.2K, ξ~1.6 μm"
        return result
    
    @staticmethod
    def niobium_cavity() -> SystemCalibration:
        """
        Niobium superconducting RF cavity.
        
        - T_c = 9.2 K
        - ξ_SC ~ 38 nm (type-II, clean limit)
        """
        xi_SC = 38e-9  # m
        result = SuperconductorCalibration.from_coherence_length(xi_SC, "Nb")
        result.notes = "Niobium RF cavity, T_c=9.2K, ξ~38 nm (type-II)"
        return result
    
    @staticmethod
    def high_tc_ybco() -> SystemCalibration:
        """
        YBCO high-T_c cuprate (optimistic).
        
        - T_c ~ 90 K
        - ξ_ab ~ 1-2 nm (ab-plane, very short)
        """
        xi_SC = 1.5e-9  # m (in-plane)
        result = SuperconductorCalibration.from_coherence_length(xi_SC, "YBCO")
        result.notes = "YBCO high-T_c cuprate, ξ_ab~1.5 nm (very short)"
        result.uncertainty_factor = 3.0  # High uncertainty for cuprates
        return result


class PlasmaCalibration:
    """
    Calibrate coherence field from plasma correlation length.
    
    Physical basis:
    - Debye length λ_D = √(ε₀k_BT / ne²)
    - Correlation length in strongly coupled plasmas
    """
    
    @staticmethod
    def from_debye_length(lambda_D: float) -> SystemCalibration:
        """
        Use Debye length as correlation scale.
        
        Args:
            lambda_D: Debye length [m]
            
        Returns:
            SystemCalibration
        """
        Phi0 = 1.0 / lambda_D
        
        return SystemCalibration(
            system_name="Plasma (Debye screening)",
            Phi0=Phi0,
            uncertainty_factor=5.0,  # Very uncertain mapping
            basis=f"λ_D = {lambda_D*1e6:.2f} μm",
            notes="Highly speculative; plasma coherence != condensate coherence"
        )
    
    @staticmethod
    def typical_lab_plasma() -> SystemCalibration:
        """
        Typical laboratory plasma.
        
        - n_e ~ 10¹⁶ m⁻³
        - T_e ~ 10 eV
        - λ_D ~ 70 μm
        """
        n_e = 1e16  # m⁻³
        T_e = 10 * 1.602e-19  # J (10 eV)
        
        # Debye length
        lambda_D = np.sqrt(8.854e-12 * T_e / (n_e * (1.602e-19)**2))
        
        result = PlasmaCalibration.from_debye_length(lambda_D)
        result.system_name = "Lab plasma (typical)"
        result.notes = f"n_e=10¹⁶ m⁻³, T_e=10 eV, λ_D={lambda_D*1e6:.1f} μm"
        
        return result


def get_all_calibrations() -> Dict[str, SystemCalibration]:
    """
    Generate all standard calibration presets.
    
    Returns:
        Dictionary mapping preset name to SystemCalibration
    """
    calibrations = {
        # BECs
        "rb87_bec": BECCalibration.rubidium_87_typical(),
        "na23_bec": BECCalibration.sodium_23_typical(),
        "bec_high_density": BECCalibration.high_density_bec(),
        
        # Superconductors
        "al_film": SuperconductorCalibration.aluminum_thin_film(),
        "nb_cavity": SuperconductorCalibration.niobium_cavity(),
        "ybco_cuprate": SuperconductorCalibration.high_tc_ybco(),
        
        # Plasma (speculative)
        "lab_plasma": PlasmaCalibration.typical_lab_plasma(),
    }
    
    return calibrations


def print_calibration_table():
    """
    Print a formatted table of all calibrations.
    """
    calibrations = get_all_calibrations()
    
    print("="*90)
    print("COHERENCE FIELD CALIBRATION TO EXPERIMENTAL SYSTEMS")
    print("="*90)
    print()
    print(f"{'System':<30} {'Φ₀ [m⁻¹]':<15} {'Uncertainty':<12} {'Basis':<30}")
    print("-"*90)
    
    for key, cal in calibrations.items():
        phi_str = f"{cal.Phi0:.2e}"
        unc_str = f"±{cal.uncertainty_factor:.1f}×"
        print(f"{cal.system_name:<30} {phi_str:<15} {unc_str:<12} {cal.basis[:30]:<30}")
    
    print("-"*90)
    print()
    print("PHYSICAL INTERPRETATION:")
    print("- BECs: Φ ≈ 1/ξ_h where ξ_h is healing length")
    print("- Superconductors: Φ ≈ 1/ξ_SC where ξ_SC is coherence length")
    print("- Plasmas: Φ ≈ 1/λ_D (highly speculative)")
    print()
    print("TYPICAL RANGES:")
    print("- Lab BECs: Φ₀ ~ 10⁶ - 10⁸ m⁻¹")
    print("- Superconductors: Φ₀ ~ 10⁷ - 10⁹ m⁻¹ (Nb cavity up to 10⁹)")
    print("- High-T_c cuprates: Φ₀ ~ 10⁹ m⁻¹ (pushing limits)")
    print()
    print("PREVIOUS CLAIM:")
    print("- \"BEC-scale = 10¹⁵ m⁻¹\" is NOT JUSTIFIED")
    print("- Real BECs: Φ₀ ~ 10⁷ m⁻¹ (factor of 10⁸ lower!)")
    print("="*90)


if __name__ == "__main__":
    print_calibration_table()
    
    # Show impact on G_eff
    print("\nIMPACT ON G_eff/G RATIOS:")
    print("-"*90)
    
    G = 6.674e-11  # m³/(kg·s²)
    xi_values = [1.0, 10.0, 100.0]
    
    calibrations = get_all_calibrations()
    
    for xi in xi_values:
        print(f"\nξ = {xi}")
        print(f"{'System':<30} {'Φ₀ [m⁻¹]':<15} {'G_eff/G':<15} {'Energy reduction':<15}")
        print("-"*80)
        
        for key, cal in calibrations.items():
            Phi0 = cal.Phi0
            denominator = 1.0 + 8*np.pi*G*xi*Phi0**2
            ratio = 1.0 / denominator
            energy_reduction = 1.0 / ratio
            
            print(f"{cal.system_name:<30} {Phi0:.2e}    {ratio:.2e}    {energy_reduction:.2e}×")
    
    print("="*90)
