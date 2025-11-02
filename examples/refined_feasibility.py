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
import sys
import matplotlib.pyplot as plt
from scipy.constants import k as k_B, hbar, c, G as G_newton
# Ensure repository root is on sys.path so `import examples...` works when running this file directly
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Physical constants (baseline geometry)
m_test_baseline = 10e-3  # Test mass 10 mg (typical for micro-torsion balance)
L_arm = 0.10    # Torsion bar arm length 10 cm
kappa = 1e-8    # Torsion constant ~10 nN·m/rad (ultra-soft fiber)
omega_0_baseline = np.sqrt(kappa / (m_test_baseline * L_arm**2))  # Natural frequency

# Noise budget presets
class NoiseProfile:
    """Experimental noise scenario with isolation and cryogenics."""
    def __init__(self, name: str, T: float, seismic_suppression: float, 
                 tilt_suppression: float, readout_improvement: float,
                 m_test_factor: float = 1.0):
        self.name = name
        self.T = T  # Temperature [K]
        self.seismic_suppression = seismic_suppression  # Factor reduction in seismic noise
        self.tilt_suppression = tilt_suppression  # Factor reduction in tilt noise
        self.readout_improvement = readout_improvement  # Factor improvement in readout
        self.m_test_factor = m_test_factor  # Test mass scaling factor
        
        # Derived params
        self.m_test = m_test_baseline * m_test_factor
        self.omega_0 = np.sqrt(kappa / (self.m_test * L_arm**2))

NOISE_PRESETS = {
    'room_temp_baseline': NoiseProfile(
        'Room temp (300K), basic isolation',
        T=300, seismic_suppression=1.0, tilt_suppression=1.0, 
        readout_improvement=1.0, m_test_factor=1.0
    ),
    'cryo_moderate': NoiseProfile(
        '4K cryo, 10× seismic, 10× tilt, 10× readout',
        T=4, seismic_suppression=10.0, tilt_suppression=10.0,
        readout_improvement=10.0, m_test_factor=1.0
    ),
    'cryo_advanced': NoiseProfile(
        '4K cryo, 30× seismic, 10× tilt, 10× readout',
        T=4, seismic_suppression=30.0, tilt_suppression=10.0,
        readout_improvement=10.0, m_test_factor=1.0
    ),
    'optimized': NoiseProfile(
        '4K cryo, 100× seismic, 10× tilt, 10× readout, 10× mass',
        T=4, seismic_suppression=100.0, tilt_suppression=10.0,
        readout_improvement=10.0, m_test_factor=10.0
    ),
}

# Noise sources (all in units of torque N·m for 1 Hz bandwidth)

def seismic_noise(profile: NoiseProfile, freq_Hz=1e-3):
    """
    Seismic acceleration noise: ~10^-9 g/sqrt(Hz) at 1 mHz.
    Couples to test mass as F_seismic = m * a_seismic.
    """
    a_seismic = 1e-9 * 9.81 / profile.seismic_suppression  # m/s^2/sqrt(Hz)
    F_seismic = profile.m_test * a_seismic
    tau_seismic = F_seismic * L_arm
    return tau_seismic

def thermal_noise(profile: NoiseProfile):
    """
    Brownian motion of torsion fiber: tau_th = sqrt(4 k_B T kappa / omega_0).
    This is the fundamental thermal limit.
    """
    # Fluctuation-dissipation theorem
    # For Q ~ 10^4 typical of good torsion fibers (improves at low T)
    Q_base = 1e4
    Q = Q_base * (300.0 / profile.T)  # Q improves ~linearly with T reduction
    tau_th = np.sqrt(4 * k_B * profile.T * kappa / (profile.omega_0 * Q))
    return tau_th

def tilt_coupling_noise(profile: NoiseProfile):
    """
    Earth tides + building tilt: ~1 μrad/sqrt(hr) at mHz.
    Couples to gravitational gradient.
    """
    tilt_rad = 1e-6 / (np.sqrt(3600) * profile.tilt_suppression)  # rad/sqrt(Hz)
    # Tilt creates spurious torque via gradient of local g
    # tau_tilt ~ m * g * L * theta
    g_local = 9.81
    tau_tilt = profile.m_test * g_local * L_arm * tilt_rad
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
    T_room = 300  # Patch noise not strongly T-dependent (uses room-temp electronics)
    delta_V = np.sqrt(4 * k_B * T_room * 1e6 * 1.0)
    F_patch = epsilon_0 * A_patch * 2 * V_patch * delta_V / d_gap**2
    tau_patch = F_patch * L_arm
    return tau_patch

def detector_readout_noise(profile: NoiseProfile):
    """
    Optical lever / capacitive readout noise.
    Modern systems: ~1 nrad/sqrt(Hz) angular resolution.
    Advanced interferometry: ~0.1 nrad/sqrt(Hz)
    """
    theta_readout = 1e-9 / profile.readout_improvement  # rad/sqrt(Hz)
    tau_readout = kappa * theta_readout
    return tau_readout

def total_noise_budget(profile: NoiseProfile):
    """
    Combine all noise sources in quadrature (assuming uncorrelated).
    """
    tau_seismic = seismic_noise(profile)
    tau_thermal = thermal_noise(profile)
    tau_tilt = tilt_coupling_noise(profile)
    tau_patch = casimir_patch_potential_noise()  # Not strongly T-dependent
    tau_readout = detector_readout_noise(profile)
    
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

def main(profile_name: str = 'room_temp_baseline'):
    profile = NOISE_PRESETS[profile_name]
    
    print("=== Refined Cavendish Feasibility with Geometric Effects ===")
    print(f"\nNoise profile: {profile.name}")
    print(f"  T = {profile.T} K")
    print(f"  Seismic suppression: {profile.seismic_suppression}×")
    print(f"  Tilt suppression: {profile.tilt_suppression}×")
    print(f"  Readout improvement: {profile.readout_improvement}×")
    print(f"  Test mass: {profile.m_test*1e3:.1f} mg\n")
    
    # Get noise budget
    noise = total_noise_budget(profile)
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
        signal_base = abs(tau_coherent - tau_newtonian)
        
        # Signal scales with test mass (torque ∝ m_test)
        signal = signal_base * profile.m_test_factor
        
        # Integration time for SNR=5
        t_int = integration_time_for_snr(signal, noise['total'], target_snr)
        
        result['signal_torque'] = signal
        result['signal_torque_base'] = signal_base
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
        print("  3. Increase test mass (signal scales as m_test)")


def sweep_noise_profiles(geometric_results, output_dir="figures"):
    """
    Sweep across all noise profiles and generate comparison plots.
    
    Parameters
    ----------
    geometric_results : list
        Geometric sweep results from geometric_cavendish_sweep.json
    output_dir : str
        Output directory for figures
    """
    print("\n" + "="*80)
    print("SWEEPING NOISE PROFILES")
    print("="*80)
    
    profile_names = list(NOISE_PRESETS.keys())
    
    # Store results for each profile
    all_results = {}
    
    for profile_name in profile_names:
        print(f"\n--- Profile: {profile_name} ---")
        profile = NOISE_PRESETS[profile_name]
        
        # Compute feasibility for this profile
        feasible_count = 0
        for result in geometric_results:
            tau_coherent = result['tau_coherent']
            tau_newtonian = result['tau_newtonian']
            signal_base = abs(tau_coherent - tau_newtonian)
            signal = signal_base * profile.m_test_factor
            
            noise = total_noise_budget(profile)
            target_snr = 5.0
            t_int = integration_time_for_snr(signal, noise['total'], target_snr)
            
            result_copy = result.copy()
            result_copy['signal_torque'] = signal
            result_copy['integration_time_hr'] = t_int / 3600
            
            if t_int < 24 * 3600:
                feasible_count += 1
        
        all_results[profile_name] = {
            'feasible_count': feasible_count,
            'total': len(geometric_results)
        }
        
        print(f"  Temperature: {profile.T} K")
        print(f"  Seismic suppression: {profile.seismic_suppression}×")
        print(f"  Tilt suppression: {profile.tilt_suppression}×")
        print(f"  Readout improvement: {profile.readout_improvement}×")
        print(f"  Test mass factor: {profile.m_test_factor}×")
        print(f"  Feasible configs: {feasible_count}/{len(geometric_results)}")
    
    # Bar plot comparison
    fig, ax = plt.subplots(figsize=(10, 6))
    
    x_pos = range(len(profile_names))
    feasible_counts = [all_results[p]['feasible_count'] for p in profile_names]
    
    bars = ax.bar(x_pos, feasible_counts, color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
    ax.axhline(18, color='gray', linestyle='--', linewidth=1, label='Total configs')
    
    ax.set_xlabel('Noise Profile', fontsize=12)
    ax.set_ylabel('Feasible Configurations (<24 hr)', fontsize=12)
    ax.set_title('Experimental Feasibility vs Noise Environment', fontsize=13)
    ax.set_xticks(x_pos)
    ax.set_xticklabels([p.replace('_', '\n') for p in profile_names], fontsize=10)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    
    # Add value labels on bars
    for bar, count in zip(bars, feasible_counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{count}', ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    output_path = Path(__file__).parent / output_dir / "noise_profile_sweep.png"
    output_path.parent.mkdir(exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\nComparison figure saved to: {output_path}")
    
    return all_results


def compare_optimized_vs_baseline(profile_name: str = 'cryo_moderate'):
    """
    Run geometry optimization and compare feasibility vs baseline.
    
    Generates a figure showing integration time improvements.
    """
    from examples.geometric_cavendish import optimize_geometry, run_geometric_cavendish
    
    profile = NOISE_PRESETS[profile_name]
    
    print("\n" + "="*80)
    print("OPTIMIZED vs BASELINE GEOMETRY COMPARISON")
    print("="*80)
    print(f"\nNoise profile: {profile.name}")
    
    # Load baseline results
    baseline_path = Path(__file__).parent.parent / "results" / "geometric_cavendish_sweep.json"
    if not baseline_path.exists():
        print(f"Error: {baseline_path} not found. Run geometric_cavendish.py first.")
        return
    
    with open(baseline_path) as f:
        baseline_results = json.load(f)
    
    # Pick a few representative configs to optimize
    test_configs = [
        {'xi': 100.0, 'Phi0': 6.67e8, 'name': 'YBCO ξ=100'},
        {'xi': 10.0, 'Phi0': 2.63e7, 'name': 'Nb ξ=10'},
        {'xi': 100.0, 'Phi0': 3.65e6, 'name': 'Rb87 ξ=100'},
    ]
    
    results_table = []
    noise = total_noise_budget(profile)
    
    for config in test_configs:
        print(f"\n--- Optimizing {config['name']} ---")
        
        # Find baseline
        baseline = None
        for r in baseline_results:
            if abs(r['xi'] - config['xi']) < 0.1 and abs(r['Phi0'] - config['Phi0']) < config['Phi0']*0.1:
                # Take offset position as baseline
                if abs(r['coherent_position'][2] + 0.08) < 0.01:
                    baseline = r
                    break
        
        if baseline is None:
            print(f"  Warning: No baseline found for {config['name']}")
            continue
        
        # Run optimization
        opt_result = optimize_geometry(xi=config['xi'], Phi0=config['Phi0'], verbose=False)
        
        # Extract optimized delta_tau
        delta_tau_baseline = abs(baseline['tau_coherent'] - baseline['tau_newtonian'])
        delta_tau_optimized = abs(opt_result['final_config']['delta_tau'])
        
        # Scale by mass factor
        signal_baseline = delta_tau_baseline * profile.m_test_factor
        signal_optimized = delta_tau_optimized * profile.m_test_factor
        
        # Integration times
        t_baseline_hr = integration_time_for_snr(signal_baseline, noise['total'], 5.0) / 3600
        t_optimized_hr = integration_time_for_snr(signal_optimized, noise['total'], 5.0) / 3600
        
        improvement = t_baseline_hr / t_optimized_hr if t_optimized_hr > 0 else 1.0
        
        results_table.append({
            'name': config['name'],
            'baseline_hr': t_baseline_hr,
            'optimized_hr': t_optimized_hr,
            'improvement': improvement
        })
        
        print(f"  Baseline T_int: {t_baseline_hr:.2f} hr")
        print(f"  Optimized T_int: {t_optimized_hr:.2f} hr")
        print(f"  Improvement: {improvement:.2f}×")
    
    # Plot comparison
    if len(results_table) > 0:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        names = [r['name'] for r in results_table]
        baseline_times = [r['baseline_hr'] for r in results_table]
        optimized_times = [r['optimized_hr'] for r in results_table]
        
        x = np.arange(len(names))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, baseline_times, width, label='Baseline', color='#1f77b4')
        bars2 = ax.bar(x + width/2, optimized_times, width, label='Optimized', color='#2ca02c')
        
        ax.axhline(24, color='red', linestyle='--', linewidth=2, label='24 hr limit')
        ax.axhline(1, color='orange', linestyle=':', linewidth=1.5, label='1 hr')
        
        ax.set_ylabel('Integration Time (hr) for SNR=5', fontsize=12)
        ax.set_title(f'Optimized vs Baseline Geometry\n({profile.name})', fontsize=13)
        ax.set_xticks(x)
        ax.set_xticklabels(names, fontsize=11)
        ax.legend(fontsize=11)
        ax.set_yscale('log')
        ax.grid(True, axis='y', alpha=0.3)
        
        # Add improvement labels
        for i, r in enumerate(results_table):
            ax.text(i, max(r['baseline_hr'], r['optimized_hr']) * 1.5,
                   f"{r['improvement']:.1f}×",
                   ha='center', fontsize=10, fontweight='bold', color='darkgreen')
        
        output_path = Path(__file__).parent / "figures" / "optimized_vs_baseline.png"
        output_path.parent.mkdir(exist_ok=True)
        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"\nComparison figure saved to: {output_path}")
    
    print("\n" + "="*80)


def run_convergence_study(
    grid_resolutions: list = [41, 61, 81, 101],
    test_configs: list = None
):
    """
    Run convergence study across multiple grid resolutions.
    
    Tests multiple representative configurations and generates:
    1. JSON file with convergence data
    2. CSV file with tabular results
    3. Matplotlib plot showing convergence trends
    
    Args:
        grid_resolutions: List of grid sizes to test
        test_configs: List of (xi, Phi0, name) tuples to test
    """
    import sys
    sys.path.insert(0, str(Path(__file__).parent))
    from geometric_cavendish import run_geometric_cavendish
    import csv
    from datetime import datetime
    
    print("\n" + "="*80)
    print("CONVERGENCE STUDY")
    print("="*80)
    print(f"Grid resolutions: {grid_resolutions}")
    
    if test_configs is None:
        # Default test configurations
        test_configs = [
            (100.0, 6.67e8, 'YBCO_xi100'),
            (10.0, 2.63e7, 'Nb_xi10'),
            (100.0, 3.65e6, 'Rb87_xi100'),
        ]
    
    all_results = []
    
    for xi, Phi0, name in test_configs:
        print(f"\n--- Testing {name} (ξ={xi}, Φ₀={Phi0:.2e}) ---")
        
        config_results = []
        
        for resolution in grid_resolutions:
            print(f"  Running {resolution}³ grid...", end=' ', flush=True)
            
            result = run_geometric_cavendish(
                xi=xi,
                Phi0=Phi0,
                geom_params={'coherent_position': (0.0, 0.0, -0.08)},
                grid_resolution=resolution,
                use_volume_average=True,
                verbose=False
            )
            
            config_results.append({
                'config_name': name,
                'xi': xi,
                'Phi0': Phi0,
                'resolution': resolution,
                'grid_points': resolution**3,
                'tau_newtonian': result['tau_newtonian'],
                'tau_coherent': result['tau_coherent'],
                'delta_tau': result['delta_tau'],
                'delta_tau_frac': result['delta_tau_frac'],
                'solve_time': result['solve_time_coherent'] + result['solve_time_newtonian']
            })
            
            print(f"Δτ = {result['delta_tau']:.4e} N·m, time = {config_results[-1]['solve_time']:.2f}s")
        
        # Compute convergence metrics
        for i in range(1, len(config_results)):
            coarse = config_results[i-1]
            fine = config_results[i]
            
            delta_tau_diff = abs(fine['delta_tau'] - coarse['delta_tau'])
            rel_error = delta_tau_diff / abs(fine['delta_tau']) if fine['delta_tau'] != 0 else 0
            
            config_results[i]['convergence_rel_error'] = rel_error
            config_results[i]['convergence_delta_tau_diff'] = delta_tau_diff
            
            # Estimate convergence order
            if i > 1:
                prev_error = config_results[i-1].get('convergence_rel_error', 0)
                if prev_error > 0 and rel_error > 0:
                    h_ratio = coarse['resolution'] / fine['resolution']
                    order = np.log(prev_error / rel_error) / np.log(h_ratio)
                    config_results[i]['convergence_order'] = order
        
        all_results.extend(config_results)
    
    # Save JSON results
    output_dir = Path(__file__).parent.parent / "results"
    output_dir.mkdir(exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    json_file = output_dir / f"convergence_{timestamp}.json"
    
    with open(json_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print(f"\n✅ JSON results saved to: {json_file}")
    
    # Save CSV results
    csv_file = output_dir / f"convergence_{timestamp}.csv"
    
    # Collect all possible field names from all results
    all_fieldnames = set()
    for result in all_results:
        all_fieldnames.update(result.keys())
    all_fieldnames = sorted(all_fieldnames)
    
    with open(csv_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=all_fieldnames)
        writer.writeheader()
        writer.writerows(all_results)
    
    print(f"✅ CSV results saved to: {csv_file}")
    
    # Generate convergence plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Plot 1: Δτ vs resolution
    ax1 = axes[0]
    for config_name in set(r['config_name'] for r in all_results):
        config_data = [r for r in all_results if r['config_name'] == config_name]
        resolutions = [r['resolution'] for r in config_data]
        delta_taus = [abs(r['delta_tau']) for r in config_data]
        ax1.plot(resolutions, delta_taus, 'o-', label=config_name, markersize=8)
    
    ax1.set_xlabel('Grid Resolution (N)', fontsize=12)
    ax1.set_ylabel('|Δτ| (N·m)', fontsize=12)
    ax1.set_title('Torque Signal vs Grid Resolution', fontsize=13, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # Plot 2: Relative convergence error
    ax2 = axes[1]
    for config_name in set(r['config_name'] for r in all_results):
        config_data = [r for r in all_results if r['config_name'] == config_name and 'convergence_rel_error' in r]
        if config_data:
            resolutions = [r['resolution'] for r in config_data]
            rel_errors = [r['convergence_rel_error'] for r in config_data]
            ax2.plot(resolutions, rel_errors, 's-', label=config_name, markersize=8)
    
    ax2.axhline(0.01, color='red', linestyle='--', linewidth=2, label='1% threshold')
    ax2.axhline(0.001, color='green', linestyle=':', linewidth=1.5, label='0.1% threshold')
    ax2.set_xlabel('Grid Resolution (N)', fontsize=12)
    ax2.set_ylabel('Relative Error: |Δτ_N - Δτ_{N-1}| / |Δτ_N|', fontsize=12)
    ax2.set_title('Grid Convergence', fontsize=13, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    
    # Plot 3: Computational time vs resolution
    ax3 = axes[2]
    for config_name in set(r['config_name'] for r in all_results):
        config_data = [r for r in all_results if r['config_name'] == config_name]
        grid_points = [r['grid_points'] for r in config_data]
        times = [r['solve_time'] for r in config_data]
        ax3.plot(grid_points, times, '^-', label=config_name, markersize=8)
    
    ax3.set_xlabel('Grid Points (N³)', fontsize=12)
    ax3.set_ylabel('Total Solve Time (s)', fontsize=12)
    ax3.set_title('Computational Cost', fontsize=13, fontweight='bold')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    
    plt.tight_layout()
    
    fig_dir = Path(__file__).parent / "figures"
    fig_dir.mkdir(exist_ok=True)
    fig_file = fig_dir / "convergence.png"
    
    plt.savefig(fig_file, dpi=300, bbox_inches='tight')
    print(f"✅ Convergence plot saved to: {fig_file}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("CONVERGENCE SUMMARY")
    print("="*80)
    
    for config_name in set(r['config_name'] for r in all_results):
        config_data = [r for r in all_results if r['config_name'] == config_name]
        finest = config_data[-1]
        
        print(f"\n{config_name}:")
        print(f"  Finest grid ({finest['resolution']}³): Δτ = {finest['delta_tau']:.4e} N·m")
        
        if 'convergence_rel_error' in finest:
            print(f"  Relative convergence: {finest['convergence_rel_error']:.4e}")
            if finest['convergence_rel_error'] < 0.01:
                print(f"  ✅ Well converged (<1%)")
            elif finest['convergence_rel_error'] < 0.05:
                print(f"  ⚠️  Moderately converged (<5%)")
            else:
                print(f"  ❌ Poor convergence (>{finest['convergence_rel_error']*100:.1f}%)")
        
        if 'convergence_order' in finest:
            print(f"  Estimated convergence order: {finest['convergence_order']:.2f}")
    
    print("\n" + "="*80)


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Compute experimental feasibility with parameterized noise models'
    )
    parser.add_argument(
        '--profile',
        type=str,
        default='room_temp_baseline',
        choices=list(NOISE_PRESETS.keys()),
        help='Noise profile to use (default: room_temp_baseline)'
    )
    parser.add_argument(
        '--sweep',
        action='store_true',
        help='Sweep across all noise profiles and compare feasibility'
    )
    parser.add_argument(
        '--optimize',
        action='store_true',
        help='Run geometry optimization and compare with baseline'
    )
    parser.add_argument(
        '--convergence',
        action='store_true',
        help='Run convergence study across multiple grid resolutions'
    )
    
    args = parser.parse_args()
    
    if args.convergence:
        run_convergence_study()
    elif args.optimize:
        compare_optimized_vs_baseline(profile_name=args.profile)
    elif args.sweep:
        # Load geometric results
        sweep_path = Path(__file__).parent.parent / "results" / "geometric_cavendish_sweep.json"
        if not sweep_path.exists():
            print(f"Error: {sweep_path} not found. Run geometric_cavendish.py first.")
            exit(1)
        
        with open(sweep_path) as f:
            geometric_results = json.load(f)
        
        sweep_noise_profiles(geometric_results)
    else:
        main(profile_name=args.profile)

