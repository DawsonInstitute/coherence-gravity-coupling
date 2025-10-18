#!/usr/bin/env python3
"""
Refined feasibility estimate for geometric Cavendish experiment.

Incorporates:
1. Geometric torque predictions from examples/geometric_cavendish.py
2. Realistic noise sources (seismic, thermal, tilt, Casimir)
3. SNR calculation with integration time scaling
4. Parameter space scan for optimal (xi, Phi_0, geometry)
"""

import numpy as np
import json
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.constants import k as k_B, hbar, c, G as G_newton

# Physical constants
m_test = 10e-3  # Test mass 10 mg (typical for micro-torsion balance)
L_arm = 0.10    # Torsion bar arm length 10 cm
kappa = 1e-8    # Torsion constant ~10 nN·m/rad (ultra-soft fiber)
omega_0 = np.sqrt(kappa / (m_test * L_arm**2))  # Natural frequency
T = 300  # Temperature in K

# Noise sources (all in units of torque N·m for 1 Hz bandwidth)

def seismic_noise(freq_Hz=1e-3):
    """
    Seismic acceleration noise: ~10^-9 g/sqrt(Hz) at 1 mHz.
    Couples to test mass as F_seismic = m * a_seismic.
    """
    a_seismic = 1e-9 * 9.81  # m/s^2/sqrt(Hz)
    F_seismic = m_test * a_seismic
    tau_seismic = F_seismic * L_arm
    return tau_seismic

def thermal_noise():
    """
    Brownian motion of torsion fiber: tau_th = sqrt(4 k_B T kappa / omega_0).
    This is the fundamental thermal limit.
    """
    # Fluctuation-dissipation theorem
    # For Q ~ 10^4 typical of good torsion fibers
    Q = 1e4
    tau_th = np.sqrt(4 * k_B * T * kappa / (omega_0 * Q))
    return tau_th

def tilt_coupling_noise():
    """
    Earth tides + building tilt: ~1 μrad/sqrt(hr) at mHz.
    Couples to gravitational gradient.
    """
    tilt_rad = 1e-6 / np.sqrt(3600)  # rad/sqrt(Hz)
    # Tilt creates spurious torque via gradient of local g
    # tau_tilt ~ m * g * L * theta
    g_local = 9.81
    tau_tilt = m_test * g_local * L_arm * tilt_rad
    return tau_tilt

def casimir_patch_potential_noise():
    """
    Electrostatic patch potentials on test mass and chamber walls.
    Typical: ~100 mV over 1 cm² patches → fluctuating force.
    """
    V_patch = 0.1  # Volts
    d_gap = 1e-3   # 1 mm gap to chamber
    epsilon_0 = 8.854e-12  # F/m
    A_patch = 1e-4  # 1 cm²
    
    # Fluctuating capacitive force: F ~ epsilon_0 A V^2 / d^2
    # With Johnson noise fluctuations delta_V ~ sqrt(4 k_B T R B)
    # For R ~ 1 MOhm, B ~ 1 Hz:
    delta_V = np.sqrt(4 * k_B * T * 1e6 * 1.0)
    F_patch = epsilon_0 * A_patch * 2 * V_patch * delta_V / d_gap**2
    tau_patch = F_patch * L_arm
    return tau_patch

def detector_readout_noise():
    """
    Optical lever / capacitive readout noise.
    Modern systems: ~1 nrad/sqrt(Hz) angular resolution.
    """
    theta_readout = 1e-9  # rad/sqrt(Hz)
    tau_readout = kappa * theta_readout
    return tau_readout

def total_noise_budget():
    """
    Combine all noise sources in quadrature (assuming uncorrelated).
    """
    tau_seismic = seismic_noise()
    tau_thermal = thermal_noise()
    tau_tilt = tilt_coupling_noise()
    tau_patch = casimir_patch_potential_noise()
    tau_readout = detector_readout_noise()
    
    tau_total = np.sqrt(
        tau_seismic**2 + 
        tau_thermal**2 + 
        tau_tilt**2 + 
        tau_patch**2 + 
        tau_readout**2
    )
    
    noise_breakdown = {
        'seismic': tau_seismic,
        'thermal': tau_thermal,
        'tilt': tau_tilt,
        'patch': tau_patch,
        'readout': tau_readout,
        'total': tau_total
    }
    
    return noise_breakdown

def load_geometric_cavendish_results():
    """
    Load torque predictions from geometric simulation sweep.
    """
    result_file = Path(__file__).parent.parent / "results" / "geometric_cavendish_sweep.json"
    if not result_file.exists():
        print(f"Warning: {result_file} not found. Run examples/geometric_cavendish.py first.")
        return None
    
    with open(result_file, 'r') as f:
        results = json.load(f)
    
    return results

def compute_snr(signal_torque, noise_torque, integration_time_s):
    """
    SNR = (signal / noise) * sqrt(T_int).
    """
    snr = (signal_torque / noise_torque) * np.sqrt(integration_time_s)
    return snr

def integration_time_for_snr(signal_torque, noise_torque, target_snr=5.0):
    """
    Required integration time to achieve target SNR.
    T_int = (target_snr * noise / signal)^2
    """
    t_int = (target_snr * noise_torque / signal_torque)**2
    return t_int

def main():
    print("=== Refined Cavendish Feasibility with Geometric Effects ===\n")
    
    # Get noise budget
    noise = total_noise_budget()
    print("Noise budget (per sqrt(Hz)):")
    for source, tau in noise.items():
        print(f"  {source:12s}: {tau:.3e} N·m/√Hz")
    print()
    
    # Load geometric results
    geometric_results = load_geometric_cavendish_results()
    if geometric_results is None:
        print("ERROR: Cannot proceed without geometric simulation data.")
        return
    
    print(f"Loaded {len(geometric_results)} geometric configurations.\n")
    
    # Compute SNR and integration times for each config
    target_snr = 5.0
    feasible_configs = []
    
    for result in geometric_results:
        tau_coherent = abs(result['tau_coherent'])
        tau_newtonian = abs(result['tau_newtonian'])
        
        # Signal is the *difference* between coherent and Newtonian
        signal = abs(tau_coherent - tau_newtonian)
        
        # Integration time for SNR=5
        t_int = integration_time_for_snr(signal, noise['total'], target_snr)
        
        result['signal_torque'] = signal
        result['snr_per_sqrt_sec'] = signal / noise['total']
        result['integration_time_hr'] = t_int / 3600
        
        # Extract coherence system name from coherent_position
        # (The JSON has Phi0 instead of explicit system name)
        phi0 = result['Phi0']
        if abs(phi0 - 3.65e6) < 1e5:
            coherence_name = 'rb87'
        elif abs(phi0 - 2.63e7) < 1e6:
            coherence_name = 'nb'
        elif abs(phi0 - 6.67e8) < 1e7:
            coherence_name = 'ybco'
        else:
            coherence_name = 'unknown'
        
        result['coherence'] = coherence_name
        
        # Position: check if offset or centered
        pos_z = result['coherent_position'][2]
        result['position'] = 'offset' if abs(pos_z + 0.08) < 0.01 else 'centered'
        
        # Flag as feasible if <24 hour integration
        if t_int < 24 * 3600:
            feasible_configs.append(result)
    
    # Sort by integration time
    feasible_configs.sort(key=lambda x: x['integration_time_hr'])
    
    print(f"Feasible configurations (T_int < 24 hr, SNR=5):")
    print(f"{'xi':>6s} {'Coherence':>10s} {'Position':>10s} {'ΔG/G':>10s} {'Signal (N·m)':>15s} {'T_int (hr)':>12s}")
    print("-" * 80)
    
    for config in feasible_configs[:10]:  # Show top 10
        print(f"{config['xi']:6.1f} {config['coherence']:>10s} {config['position']:>10s} "
              f"{config['delta_G_over_G']:10.2f} {config['signal_torque']:15.3e} "
              f"{config['integration_time_hr']:12.1f}")
    
    if len(feasible_configs) == 0:
        print("  (None - geometric signal below noise even at 24 hr)")
    
    print(f"\nTotal feasible: {len(feasible_configs)} / {len(geometric_results)}")
    
    # Generate contour plot: T_int vs (xi, system)
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    systems = ['rb87', 'nb', 'ybco']
    positions = ['offset', 'centered']
    
    for idx, position in enumerate(positions):
        ax = axes[idx]
        
        for system in systems:
            subset = [r for r in geometric_results 
                     if r['coherence'] == system and r['position'] == position]
            
            if len(subset) == 0:
                continue
            
            xi_vals = [r['xi'] for r in subset]
            t_int_vals = [r['integration_time_hr'] for r in subset]
            
            ax.plot(xi_vals, t_int_vals, 'o-', label=system, markersize=8)
        
        ax.axhline(24, color='red', linestyle='--', linewidth=2, label='24 hr limit')
        ax.axhline(1, color='green', linestyle=':', linewidth=1.5, label='1 hr')
        
        ax.set_xlabel('Coherence length ξ', fontsize=12)
        ax.set_ylabel('Integration time (hr) for SNR=5', fontsize=12)
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_title(f'Position: {position}', fontsize=13)
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_path = Path(__file__).parent / "figures" / "feasibility_integration_times.png"
    output_path.parent.mkdir(exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nFigure saved to: {output_path}")
    
    # Summary statistics
    print(f"\n=== Key Findings ===")
    if len(feasible_configs) > 0:
        best = feasible_configs[0]
        print(f"Best configuration:")
        print(f"  System: {best['coherence']}, ξ={best['xi']}, position={best['position']}")
        print(f"  ΔG/G: {best['delta_G_over_G']:.2f}")
        print(f"  Signal: {best['signal_torque']:.3e} N·m")
        print(f"  Integration time: {best['integration_time_hr']:.1f} hr")
        print(f"  SNR per second: {best['snr_per_sqrt_sec']:.2e}")
    else:
        print("No configurations feasible within 24 hr.")
        print("Geometric effects are strong but current noise floor is too high.")
        print("Recommendations:")
        print("  1. Lower temperature to 4K (reduces thermal noise by ~9×)")
        print("  2. Improve seismic isolation (active platform)")
        print("  3. Increase test mass (signal scales as m²)")

if __name__ == '__main__':
    main()
