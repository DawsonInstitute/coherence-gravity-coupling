#!/usr/bin/env python3
"""
Overlay empirical lab constraints on coherence-gravity parameter space.

Reads data/empirical_constraints.json and plots exclusion regions on
(xi, Phi_0) space alongside theoretical predictions.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from pathlib import Path

# Load constraint data
constraint_file = Path(__file__).parent / "data" / "empirical_constraints.json"
with open(constraint_file, 'r') as f:
    constraints = json.load(f)

# Phi_0 values for known coherent systems (from phi_calibration.py)
PHI_CALIBRATIONS = {
    'rb87_bec': 3.65e6,      # m^-1
    'nb_superconductor': 2.63e7,
    'ybco_superconductor': 6.67e8
}

# Create figure
fig, ax = plt.subplots(figsize=(10, 8))

# Plot theoretical prediction contours (ΔG/G levels)
xi_range = np.logspace(-1, 3, 100)  # 0.1 to 1000
phi_range = np.logspace(5, 10, 100)  # 10^5 to 10^10 m^-1

Xi, Phi = np.meshgrid(xi_range, phi_range)

# Approximate ΔG/G ~ -Φ₀² ξ² / (m_Pl²) for negative coupling
# (This is path-average; geometric effects give -350% to +640%)
m_Pl_inv = 6.7e-11  # SI units: m^3 kg^-1 s^-2 → m^-1 via c=1 units
Delta_G_over_G = -(Phi**2) * (Xi**2) * (m_Pl_inv**2)

# Plot contours at experimental sensitivity levels
levels = [-1e-1, -1e-2, -1e-3, -1e-4, -1e-5, -1e-6, -1e-7]
contours = ax.contour(Xi, Phi, Delta_G_over_G, levels=levels, 
                      colors='gray', alpha=0.4, linewidths=0.5)
ax.clabel(contours, inline=True, fontsize=8, fmt='ΔG/G=%.0e')

# Plot calibrated systems as scatter points
for name, phi0 in PHI_CALIBRATIONS.items():
    # Typical coherence lengths for each system
    if 'bec' in name:
        xi_typical = [0.5, 1.0, 5.0]  # microns
    else:  # superconductors
        xi_typical = [10.0, 100.0, 300.0]  # nanometers
    
    for xi in xi_typical:
        ax.scatter(xi, phi0, s=100, marker='*', 
                  label=f'{name} (ξ={xi})' if xi == xi_typical[0] else "",
                  alpha=0.8, edgecolors='black', linewidths=0.5)

# Overlay experimental constraints as exclusion regions
for constraint in constraints['constraints']:
    if constraint['delta_G_over_G_limit'] is None:
        continue
    
    limit = abs(constraint['delta_G_over_G_limit'])
    length_scale = constraint['length_scale_m'][1]  # Upper end of range
    
    # Estimate corresponding (xi, Phi0) exclusion
    # If |ΔG/G| < limit, then Φ₀² ξ² < limit / m_Pl_inv²
    # For given xi, Phi0_max = sqrt(limit) / (xi * m_Pl_inv)
    
    xi_constraint = np.logspace(-1, 3, 50)
    phi_max = np.sqrt(limit) / (xi_constraint * m_Pl_inv)
    
    # Only plot if within our axis range
    mask = (phi_max > 1e5) & (phi_max < 1e10)
    if np.any(mask):
        label = constraint['experiment'][:20]  # Truncate for legend
        ax.fill_between(xi_constraint[mask], phi_max[mask], 1e10, 
                        alpha=0.15, label=label)

# Highlight "most promising" region (where geometric effects dominate)
# This is xi > 10 nm and Phi0 > 10^7 m^-1
promising_region = Rectangle((10, 1e7), 1000-10, 1e10-1e7, 
                             linewidth=2, edgecolor='green', 
                             facecolor='none', linestyle='--',
                             label='Geometric effects (|ΔG/G| > 1)')

ax.add_patch(promising_region)

# Formatting
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Coherence length ξ [dimensionless]', fontsize=12)
ax.set_ylabel('Order parameter Φ₀ [m⁻¹]', fontsize=12)
ax.set_title('Coherence-Gravity Parameter Space with Experimental Constraints', fontsize=14)
ax.legend(loc='upper right', fontsize=8, ncol=2)
ax.grid(True, alpha=0.3)

# Save figure
output_dir = Path(__file__).parent / "figures"
output_dir.mkdir(exist_ok=True)
output_path = output_dir / "parameter_space_with_constraints.png"
plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Figure saved to: {output_path}")

# Print summary statistics
print("\n=== Empirical Constraint Summary ===")
print(f"Total experiments surveyed: {len(constraints['constraints'])}")
print(f"\nStrongest limits by category:")
for constraint in sorted(constraints['constraints'], 
                         key=lambda x: x['delta_G_over_G_limit'] if x['delta_G_over_G_limit'] else 1.0):
    if constraint['delta_G_over_G_limit'] is not None:
        print(f"  {constraint['experiment']:40s}: |ΔG/G| < {abs(constraint['delta_G_over_G_limit']):.2e}")

print(f"\nRecommended next experiments:")
for rec in constraints['recommended_tests']:
    print(f"  - {rec['test_name']} (feasibility: {rec['feasibility']})")
    print(f"    Predicted signal: {rec['predicted_signal']}")
    print(f"    Required precision: {rec['required_precision']:.2e}")
    print()

plt.show()
