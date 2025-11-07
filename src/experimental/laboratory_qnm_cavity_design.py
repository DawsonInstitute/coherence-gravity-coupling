#!/usr/bin/env python3
"""
Laboratory quasi-normal mode (QNM) cavity design for detecting κ_R ≠ 0.

This module implements detailed engineering calculations for a superconducting
microwave cavity with tunable effective curvature, designed to detect frequency
shifts proportional to κ_R R_eff.

Physical principle:
    In curved spacetime, cavity resonance frequencies shift:
    
    ω_QNM = ω_0 + δω_curv + O(κ_R R)
    
    where:
    - ω_0: Flat-space resonance (10 GHz for X-band microwave)
    - δω_curv: GR curvature correction (~nHz for Earth's gravity)
    - κ_R R term: New physics contribution (target: mHz resolution)

NEW PHYSICS signature:
    If Δf ∝ R linearly with slope proportional to κ_R, this constitutes
    smoking-gun evidence for curvature-EM coupling.

Engineering constraints:
    - Quality factor Q > 10^9 (superconducting niobium cavity)
    - Frequency stability δf/f < 10^-15 (hydrogen maser reference)
    - Tunable R via rotating mass or centrifuge: 10^-12 to 10^-8 m^-2
    - Operating temperature: T < 4.2 K (liquid helium)
    - Vibration isolation: <1 nm rms displacement

Author: DawsonInstitute
Date: November 6, 2025
"""

import numpy as np
from typing import Tuple, Dict, List, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar
from scipy.constants import c, epsilon_0, mu_0, h, k as k_B


@dataclass
class CavityParameters:
    """Physical parameters for superconducting microwave cavity."""
    
    # Geometry
    radius: float = 0.05  # Cavity radius [m] (5 cm typical for X-band)
    length: float = 0.10  # Cavity length [m] (10 cm)
    wall_thickness: float = 5e-3  # Wall thickness [m] (5 mm niobium)
    
    # Material properties (Niobium Nb, superconducting)
    critical_temp: float = 9.2  # K (Nb Tc)
    london_penetration: float = 40e-9  # London penetration depth [m]
    coherence_length: float = 40e-9  # Superconducting coherence length [m]
    surface_resistance: float = 1e-9  # Ohm (at T << Tc, 10 GHz)
    
    # Operating conditions
    operating_temp: float = 2.0  # K (pumped liquid helium)
    mode_number: Tuple[int, int, int] = (0, 1, 0)  # TM_010 mode (azimuthally symmetric)
    
    # Performance targets
    target_Q: float = 1e9  # Quality factor
    target_stability: float = 1e-15  # Fractional frequency stability δf/f
    
    def compute_unloaded_Q(self) -> float:
        """
        Compute theoretical unloaded quality factor.
        
        For superconducting cavity:
        Q_0 = (μ_0 ω λ_L) / R_s
        
        where λ_L is London penetration depth, R_s surface resistance.
        """
        omega = self.compute_resonance_frequency()
        Q = (mu_0 * omega * self.london_penetration) / self.surface_resistance
        return Q
    
    def compute_resonance_frequency(self) -> float:
        """
        Compute fundamental resonance frequency for TM_010 mode.
        
        For cylindrical cavity:
        f_010 = (x_01 * c) / (2π * radius)
        
        where x_01 = 2.405 is first zero of J_0 Bessel function.
        """
        x_01 = 2.405  # First zero of J_0
        f = (x_01 * c) / (2 * np.pi * self.radius)
        return f


@dataclass
class CurvatureTuningMechanism:
    """Mechanisms for generating tunable effective curvature in laboratory."""
    
    # Rotating platform parameters
    rotation_radius: float = 1.0  # Radius of rotation [m]
    max_angular_velocity: float = 10.0  # Max rotation rate [rad/s]
    mass_distribution: str = "point"  # "point", "ring", or "disk"
    
    # Mass parameters
    rotating_mass: float = 100.0  # kg (test mass on rotating platform)
    center_mass: float = 1000.0  # kg (central heavy mass for static curvature)
    
    def compute_effective_curvature(self, omega: float) -> float:
        """
        Compute effective Ricci scalar from rotating frame.
        
        In rotating frame with angular velocity Ω:
        R_eff ≈ -12 Ω² / c² (for r << c/Ω)
        
        Plus static contribution from central mass:
        R_static ≈ -2 G M / (c² r³)
        
        Args:
            omega: Angular velocity [rad/s]
            
        Returns:
            Effective Ricci scalar [m^-2]
        """
        # Rotational contribution (dominant)
        R_rotation = -12 * omega**2 / c**2
        
        # Static gravitational contribution
        G = 6.674e-11  # m^3 kg^-1 s^-2
        R_static = -2 * G * self.center_mass / (c**2 * self.rotation_radius**3)
        
        # Total effective curvature
        R_eff = R_rotation + R_static
        
        return R_eff
    
    def rotation_rate_for_target_curvature(self, R_target: float) -> float:
        """
        Compute required rotation rate to achieve target curvature.
        
        Args:
            R_target: Desired curvature [m^-2]
            
        Returns:
            Required angular velocity [rad/s]
        """
        # Solve: R_rotation = -12 Ω² / c² for Ω
        omega_squared = -R_target * c**2 / 12
        
        if omega_squared < 0:
            raise ValueError(f"Cannot achieve R_target = {R_target} with rotation (need R < 0)")
        
        omega = np.sqrt(omega_squared)
        
        return omega
    
    def centrifugal_acceleration(self, omega: float) -> float:
        """Compute centrifugal acceleration at rotation radius."""
        a = omega**2 * self.rotation_radius
        return a
    
    def is_safe_operation(self, omega: float, 
                         max_g_force: float = 10.0) -> bool:
        """
        Check if rotation rate is safe for equipment.
        
        Args:
            omega: Proposed angular velocity [rad/s]
            max_g_force: Maximum allowable acceleration [in units of g]
            
        Returns:
            True if safe, False otherwise
        """
        a = self.centrifugal_acceleration(omega)
        g = 9.81  # m/s^2
        
        g_force = a / g
        
        return g_force <= max_g_force


class QNMFrequencyShift:
    """
    Calculate QNM frequency shift due to κ_R coupling.
    
    Key formula:
        Δf / f_0 = (κ_R R_eff / 4π) × (|E|² + |B|²) / (ε_0 E_vac²)
    
    where E_vac = m_e c² / (e ℓ_Compton) is vacuum field scale.
    """
    
    def __init__(self, cavity: CavityParameters):
        self.cavity = cavity
        
        # Fundamental constants
        self.e = 1.602e-19  # C (elementary charge)
        self.m_e = 9.109e-31  # kg (electron mass)
        self.l_compton = 2.426e-12  # m (Compton wavelength)
        
        # Vacuum field scale
        self.E_vac = (self.m_e * c**2) / (self.e * self.l_compton)  # V/m
        
    def compute_cavity_field_energy(self, 
                                    Q_factor: float,
                                    input_power: float = 1e-3) -> float:
        """
        Compute stored EM field energy in cavity.
        
        Args:
            Q_factor: Quality factor
            input_power: Input power [W]
            
        Returns:
            Stored energy |E|² + |B|² [J/m³]
        """
        # Resonance frequency
        f_0 = self.cavity.compute_resonance_frequency()
        omega_0 = 2 * np.pi * f_0
        
        # Cavity volume
        V = np.pi * self.cavity.radius**2 * self.cavity.length
        
        # Stored energy at resonance
        # U = Q P_in / ω_0
        U = Q_factor * input_power / omega_0
        
        # Energy density
        u = U / V
        
        # Field energy density: u = (ε_0/2)|E|² + (1/2μ_0)|B|²
        # For TM mode: |E|² ≈ |B|²c², so u ≈ ε_0 |E|²
        E_squared = u / epsilon_0
        
        return E_squared
    
    def predict_frequency_shift(self,
                               kappa_R: float,
                               R_eff: float,
                               Q_factor: float,
                               input_power: float = 1e-3) -> float:
        """
        Predict absolute frequency shift Δf due to κ_R coupling.
        
        Args:
            kappa_R: Curvature-EM coupling [m²]
            R_eff: Effective curvature scalar [m^-2]
            Q_factor: Cavity quality factor
            input_power: Input microwave power [W]
            
        Returns:
            Frequency shift Δf [Hz]
        """
        # Base frequency
        f_0 = self.cavity.compute_resonance_frequency()
        
        # Cavity field energy density
        E_squared = self.compute_cavity_field_energy(Q_factor, input_power)
        
        # Fractional frequency shift
        # δf/f ≈ (κ_R R / 4π) × (E² / E_vac²)
        delta_f_over_f = (kappa_R * R_eff / (4 * np.pi)) * (E_squared / self.E_vac**2)
        
        # Absolute shift
        delta_f = delta_f_over_f * f_0
        
        return delta_f
    
    def minimum_detectable_kappa_R(self,
                                   R_eff: float,
                                   Q_factor: float,
                                   frequency_resolution: float = 1e-3,
                                   input_power: float = 1e-3) -> float:
        """
        Compute minimum detectable κ_R given experimental parameters.
        
        Args:
            R_eff: Effective curvature [m^-2]
            Q_factor: Cavity quality factor
            frequency_resolution: Minimum resolvable Δf [Hz]
            input_power: Input power [W]
            
        Returns:
            Minimum κ_R [m²] detectable at 3σ
        """
        # For Δf = frequency_resolution, solve for κ_R
        f_0 = self.cavity.compute_resonance_frequency()
        E_squared = self.compute_cavity_field_energy(Q_factor, input_power)
        
        # κ_R,min = (4π Δf / f_0) × (E_vac² / E²) / R_eff
        kappa_R_min = (4 * np.pi * frequency_resolution / f_0) * \
                      (self.E_vac**2 / E_squared) / abs(R_eff)
        
        # 3σ detection threshold
        kappa_R_min *= 3
        
        return kappa_R_min


def design_optimization_study(kappa_R_target: float = 1e17) -> Dict:
    """
    Optimize cavity design parameters to maximize sensitivity to κ_R.
    
    Strategy:
    1. Maximize Q (superconducting cavity: Q ~ 10^9 - 10^11)
    2. Maximize R_eff (rotation + heavy central mass)
    3. Maximize input power (limited by cavity heating)
    4. Optimize cavity radius for highest field concentration
    
    Args:
        kappa_R_target: Target κ_R to detect [m²]
        
    Returns:
        Dictionary of optimized parameters
    """
    
    print("\n" + "="*60)
    print("LABORATORY QNM CAVITY DESIGN OPTIMIZATION")
    print("="*60)
    print(f"Target: Detect κ_R = {kappa_R_target:.2e} m² at 3σ\n")
    
    # Scan over cavity radius
    radii = np.linspace(0.01, 0.20, 50)  # 1 cm to 20 cm
    sensitivities = []
    
    for radius in radii:
        cavity = CavityParameters(radius=radius)
        qnm = QNMFrequencyShift(cavity)
        
        # Assume optimistic parameters
        Q = 5e9  # Achievable with Nb cavity at 2K
        R_eff = -1e-9  # From rotation at ~1 rad/s
        freq_res = 1e-3  # Hz (1 mHz, achievable with long integration)
        
        kappa_min = qnm.minimum_detectable_kappa_R(
            R_eff=R_eff,
            Q_factor=Q,
            frequency_resolution=freq_res,
            input_power=1e-3
        )
        
        sensitivities.append(kappa_min)
    
    # Find optimal radius
    optimal_idx = np.argmin(sensitivities)
    optimal_radius = radii[optimal_idx]
    optimal_sensitivity = sensitivities[optimal_idx]
    
    # Generate final design
    cavity_design = CavityParameters(radius=optimal_radius)
    tuning = CurvatureTuningMechanism()
    qnm_final = QNMFrequencyShift(cavity_design)
    
    # Performance metrics
    f_0 = cavity_design.compute_resonance_frequency()
    Q_unloaded = cavity_design.compute_unloaded_Q()
    
    # Curvature scan
    R_min = -1e-11  # m^-2 (static gravity only)
    R_max = -1e-7   # m^-2 (rotation at max safe rate)
    
    omega_max = tuning.rotation_rate_for_target_curvature(R_max)
    g_force_max = tuning.centrifugal_acceleration(omega_max) / 9.81
    
    # Predicted signal
    delta_f_max = qnm_final.predict_frequency_shift(
        kappa_R=kappa_R_target,
        R_eff=R_max,
        Q_factor=Q_unloaded,
        input_power=1e-3
    )
    
    design_summary = {
        "cavity_radius": optimal_radius,
        "cavity_length": cavity_design.length,
        "resonance_frequency": f_0,
        "quality_factor": Q_unloaded,
        "operating_temp": cavity_design.operating_temp,
        "curvature_range": (R_min, R_max),
        "rotation_rate_max": omega_max,
        "g_force_max": g_force_max,
        "predicted_signal": delta_f_max,
        "min_detectable_kappa_R": optimal_sensitivity,
        "is_safe": tuning.is_safe_operation(omega_max, max_g_force=10.0),
        "cost_estimate_USD": estimate_cost(cavity_design, tuning)
    }
    
    return design_summary


def estimate_cost(cavity: CavityParameters, 
                 tuning: CurvatureTuningMechanism) -> float:
    """
    Rough cost estimate for laboratory QNM experiment.
    
    Major components:
    - Superconducting cavity: $50K (niobium fabrication)
    - Cryostat (liquid helium): $100K
    - Rotating platform with slip rings: $75K
    - Frequency stabilization (H-maser reference): $150K
    - Vibration isolation table: $30K
    - Data acquisition / control: $20K
    - Installation / commissioning: $75K
    
    Total: ~$500K
    """
    
    costs = {
        "superconducting_cavity": 50000,
        "cryostat_system": 100000,
        "rotating_platform": 75000,
        "frequency_reference": 150000,
        "vibration_isolation": 30000,
        "daq_control": 20000,
        "installation": 75000
    }
    
    total = sum(costs.values())
    
    return total


def plot_sensitivity_analysis(output_file: str = "figures/qnm_cavity_sensitivity.pdf"):
    """Generate comprehensive sensitivity plots."""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Sensitivity vs cavity radius
    ax = axes[0, 0]
    radii = np.linspace(0.01, 0.20, 100)
    kappa_min_list = []
    
    for r in radii:
        cavity = CavityParameters(radius=r)
        qnm = QNMFrequencyShift(cavity)
        
        kappa_min = qnm.minimum_detectable_kappa_R(
            R_eff=-1e-9,
            Q_factor=5e9,
            frequency_resolution=1e-3
        )
        kappa_min_list.append(kappa_min)
    
    ax.semilogy(radii * 100, kappa_min_list, 'b-', linewidth=2)
    ax.axhline(y=5e17, color='r', linestyle='--', label='Lab limit (95% CL)')
    ax.set_xlabel("Cavity radius [cm]")
    ax.set_ylabel(r"Min detectable $\kappa_R$ [m²]")
    ax.set_title("Sensitivity vs Cavity Size")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Frequency shift vs curvature
    ax = axes[0, 1]
    R_values = np.logspace(-12, -7, 100)
    delta_f_list = []
    
    cavity = CavityParameters(radius=0.05)
    qnm = QNMFrequencyShift(cavity)
    
    for R in -R_values:  # Negative for Schwarzschild-type
        delta_f = qnm.predict_frequency_shift(
            kappa_R=1e17,
            R_eff=R,
            Q_factor=5e9
        )
        delta_f_list.append(abs(delta_f))
    
    ax.loglog(R_values, delta_f_list, 'g-', linewidth=2)
    ax.axhline(y=1e-3, color='orange', linestyle='--', label='Freq. resolution (1 mHz)')
    ax.set_xlabel(r"Effective curvature $|R_{\rm eff}|$ [m$^{-2}$]")
    ax.set_ylabel(r"Frequency shift $|\Delta f|$ [Hz]")
    ax.set_title(r"Signal vs Curvature ($\kappa_R = 10^{17}$ m²)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 3: Q-factor requirements
    ax = axes[1, 0]
    Q_values = np.logspace(7, 11, 100)
    kappa_min_Q = []
    
    for Q in Q_values:
        kappa_min = qnm.minimum_detectable_kappa_R(
            R_eff=-1e-9,
            Q_factor=Q,
            frequency_resolution=1e-3
        )
        kappa_min_Q.append(kappa_min)
    
    ax.loglog(Q_values, kappa_min_Q, 'm-', linewidth=2)
    ax.axhline(y=5e17, color='r', linestyle='--', label='Lab limit')
    ax.axvline(x=1e9, color='b', linestyle=':', label='Nb cavity (achievable)')
    ax.set_xlabel("Quality factor Q")
    ax.set_ylabel(r"Min detectable $\kappa_R$ [m²]")
    ax.set_title("Sensitivity vs Quality Factor")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 4: SNR vs integration time
    ax = axes[1, 1]
    integration_times = np.logspace(0, 6, 100)  # 1 s to ~10 days
    snr_list = []
    
    # Assume thermal noise limited
    # SNR ∝ √(t_int)
    
    delta_f_signal = qnm.predict_frequency_shift(
        kappa_R=1e17,
        R_eff=-1e-9,
        Q_factor=5e9
    )
    
    # Frequency noise: δf ~ 1 Hz / √(t_int)
    for t in integration_times:
        freq_noise = 1.0 / np.sqrt(t)
        snr = abs(delta_f_signal) / freq_noise
        snr_list.append(snr)
    
    ax.loglog(integration_times / 3600, snr_list, 'c-', linewidth=2)
    ax.axhline(y=3, color='k', linestyle='--', label='3σ threshold')
    ax.set_xlabel("Integration time [hours]")
    ax.set_ylabel("Signal-to-noise ratio")
    ax.set_title(r"SNR vs Measurement Duration ($\kappa_R = 10^{17}$ m²)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved sensitivity analysis to {output_file}")


def main():
    """Run full QNM cavity design study."""
    
    # Target: Lab limit κ_R < 5×10^17 m²
    design = design_optimization_study(kappa_R_target=1e17)
    
    print("\n" + "="*60)
    print("OPTIMIZED DESIGN PARAMETERS")
    print("="*60)
    
    print(f"\nCavity specifications:")
    print(f"  Radius: {design['cavity_radius']*100:.2f} cm")
    print(f"  Length: {design['cavity_length']*100:.2f} cm")
    print(f"  Resonance frequency: {design['resonance_frequency']/1e9:.3f} GHz")
    print(f"  Quality factor: {design['quality_factor']:.2e}")
    print(f"  Operating temperature: {design['operating_temp']:.1f} K")
    
    print(f"\nCurvature tuning:")
    print(f"  Range: {design['curvature_range'][0]:.2e} to {design['curvature_range'][1]:.2e} m^-2")
    print(f"  Max rotation rate: {design['rotation_rate_max']:.2f} rad/s")
    print(f"  Max g-force: {design['g_force_max']:.1f} g")
    print(f"  Safe operation: {design['is_safe']}")
    
    print(f"\nSensitivity:")
    print(f"  Predicted signal (κ_R = 10^17 m²): {design['predicted_signal']*1e3:.3f} mHz")
    print(f"  Minimum detectable κ_R: {design['min_detectable_kappa_R']:.2e} m²")
    
    comparison = design['min_detectable_kappa_R'] / 5e17
    print(f"  Compared to lab limit: {comparison:.2%} of current bound")
    
    print(f"\nCost estimate: ${design['cost_estimate_USD']:,.0f} USD")
    
    print("\n" + "="*60)
    print("EXPERIMENTAL PROTOCOL")
    print("="*60)
    
    print("""
1. Cool cavity to T = 2.0 K using pumped liquid helium
2. Lock frequency reference to hydrogen maser (stability < 10^-15)
3. Scan rotation rate from 0 to max safe value (~10 rad/s)
4. Measure cavity resonance frequency at each R_eff
5. Fit Δf vs R: If slope ≠ 0 with >3σ → κ_R detected
6. Cross-check: Reverse rotation direction, signal should flip sign

Timeline:
  - Design phase: 6 months
  - Fabrication: 12 months
  - Commissioning: 6 months
  - Data collection: 3 months
  - Analysis & publication: 3 months
  
Total: ~2.5 years from funding to first results

Collaborations:
  - NIST (superconducting cavity expertise)
  - PTB (frequency metrology)
  - University labs (rotating platform engineering)
""")
    
    # Generate plots
    plot_sensitivity_analysis()
    
    print("\n" + "="*60)
    print("NEW PHYSICS DISCOVERY POTENTIAL")
    print("="*60)
    
    if design['min_detectable_kappa_R'] < 5e17:
        print(f"\n✓ This design CAN IMPROVE on current lab limits!")
        print(f"  Expected sensitivity: {design['min_detectable_kappa_R']:.2e} m²")
        print(f"  Current limit: 5.00e+17 m²")
        improvement = 5e17 / design['min_detectable_kappa_R']
        print(f"  Improvement factor: {improvement:.1f}×")
    else:
        print(f"\n✗ This design does NOT reach current lab limit")
        print(f"  Would need Q × R_eff product increased by:")
        shortfall = design['min_detectable_kappa_R'] / 5e17
        print(f"  Factor of {shortfall:.1f}×")
        print("\n  Options to improve:")
        print("  - Use sapphire cavity (Q ~ 10^10)")
        print("  - Increase rotation rate (stronger centrifuge)")
        print("  - Add heavy central mass (static R contribution)")


if __name__ == "__main__":
    import os
    os.makedirs("figures", exist_ok=True)
    main()
