#!/usr/bin/env python3
"""
Generate plots for BSM parameter space paper (curvature_em_to_bsm.tex).

Creates:
1. ε_eff vs R (curvature) for fixed κ_R
2. g_aγγ vs R for fixed κ_R
3. Curvature amplification visualization
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add src to path
ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(ROOT / 'src'))

try:
    from analysis.bsm_bounds_from_kappa import (
        epsilon_equiv,
        axion_equiv_parametric,
        DEFAULT_ENVIRONMENTS
    )
except ImportError:
    print("Warning: Could not import BSM module, using fallback")
    epsilon_equiv = lambda k, R, C: C * k * R
    axion_equiv_parametric = lambda k, R, C, L: C * k * R / L
    DEFAULT_ENVIRONMENTS = {
        'lab_flat': type('obj', (), {'name': 'Lab (flat)', 'R_m2': 1e-30}),
        'earth_surface': type('obj', (), {'name': 'Earth surface', 'R_m2': 1e-26}),
        'leo': type('obj', (), {'name': 'Low Earth orbit', 'R_m2': 5e-27}),
        'magnetar_surface': type('obj', (), {'name': 'Magnetar surface', 'R_m2': 1e-6})
    }


def generate_epsilon_vs_R_plot(output_dir: Path):
    """Generate ε_eff vs curvature R plot."""
    print("Generating ε_eff vs R plot...")
    
    # Curvature range: 10^-30 to 10^-5 m^-2
    R_values = np.logspace(-30, -5, 100)
    
    # κ_R values to plot
    kappa_values = [1e-11, 5e17, 1e18]  # m^2
    kappa_labels = [r'$\kappa_R = 10^{-11}\,\mathrm{m}^2$',
                    r'$\kappa_R = 5 \times 10^{17}\,\mathrm{m}^2$ (lab limit)',
                    r'$\kappa_R = 10^{18}\,\mathrm{m}^2$']
    
    C_eps = 1.0  # Matching coefficient
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for kappa, label in zip(kappa_values, kappa_labels):
        eps_values = [epsilon_equiv(kappa, R, C_eps) for R in R_values]
        ax.loglog(R_values, eps_values, label=label, linewidth=2)
    
    # Mark environment points
    for env_key, env in DEFAULT_ENVIRONMENTS.items():
        for kappa in kappa_values:
            eps = epsilon_equiv(kappa, env.R_m2, C_eps)
            if eps > 0:
                ax.scatter(env.R_m2, eps, s=50, zorder=5)
    
    # Experimental limits
    ax.axhline(1e-3, color='red', linestyle='--', alpha=0.5, label='APEX limit')
    ax.axhline(1e-4, color='orange', linestyle='--', alpha=0.5, label='BaBar limit')
    
    ax.set_xlabel(r'Curvature $\mathcal{R}$ [m$^{-2}$]', fontsize=12)
    ax.set_ylabel(r'$\varepsilon_{\rm eff}$ (dark photon mixing)', fontsize=12)
    ax.set_title(r'Dark Photon Mixing vs Curvature ($C_\varepsilon=1$)', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    
    # Add annotations for environments
    env_labels = {
        'lab_flat': (1e-30, 1.5e-30, 'Lab'),
        'earth_surface': (1e-26, 1.5e-26, 'Earth'),
        'magnetar_surface': (1e-6, 1.5e-6, 'Magnetar')
    }
    for env_key, (x, y, text) in env_labels.items():
        ax.text(x, y, text, fontsize=9, alpha=0.7, rotation=45)
    
    plt.tight_layout()
    
    output_file = output_dir / 'epsilon_vs_curvature.pdf'
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    
    # Also save PNG
    fig.savefig(output_file.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def generate_axion_vs_R_plot(output_dir: Path):
    """Generate g_aγγ vs curvature R plot."""
    print("Generating g_aγγ vs R plot...")
    
    # Curvature range
    R_values = np.logspace(-30, -5, 100)
    
    # κ_R values
    kappa_values = [1e-11, 5e17]
    kappa_labels = [r'$\kappa_R = 10^{-11}\,\mathrm{m}^2$',
                    r'$\kappa_R = 5 \times 10^{17}\,\mathrm{m}^2$ (lab limit)']
    
    C_a = 1.0
    Lambda_GeV = 1e4  # 10 TeV
    
    fig, ax = plt.subplots(figsize=(8, 6))
    
    for kappa, label in zip(kappa_values, kappa_labels):
        g_values = [axion_equiv_parametric(kappa, R, C_a, Lambda_GeV) for R in R_values]
        ax.loglog(R_values, g_values, label=label, linewidth=2)
    
    # Experimental limits
    ax.axhline(1e-10, color='red', linestyle='--', alpha=0.5, label='CAST limit')
    ax.axhline(1e-14, color='orange', linestyle='--', alpha=0.5, label='ADMX sensitivity')
    
    ax.set_xlabel(r'Curvature $\mathcal{R}$ [m$^{-2}$]', fontsize=12)
    ax.set_ylabel(r'$g_{a\gamma\gamma}^{\rm equiv}$ [GeV$^{-1}$]', fontsize=12)
    ax.set_title(r'Axion Coupling Benchmark vs Curvature ($\Lambda=10$ TeV, $C_a=1$)', fontsize=14)
    ax.legend(fontsize=10, loc='best')
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    
    output_file = output_dir / 'axion_vs_curvature.pdf'
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    
    fig.savefig(output_file.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def generate_amplification_plot(output_dir: Path):
    """Generate curvature amplification visualization."""
    print("Generating curvature amplification plot...")
    
    # Environment curvatures
    environments = [
        ('Lab (flat)', 1e-30),
        ('Earth surface', 1e-26),
        ('Low Earth orbit', 5e-27),
        ('Magnetar surface', 1e-6)
    ]
    
    kappa_R = 1e-11  # m^2
    C_eps = 1.0
    
    env_names = [e[0] for e in environments]
    R_values = [e[1] for e in environments]
    eps_values = [epsilon_equiv(kappa_R, R, C_eps) for R in R_values]
    
    # Calculate amplification factors relative to lab
    amplification = [eps / eps_values[0] for eps in eps_values]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10))  # Stack vertically instead of horizontally
    
    # Plot 1: ε_eff by environment
    bars1 = ax1.bar(range(len(env_names)), eps_values, color=['blue', 'green', 'orange', 'red'], alpha=0.7)
    ax1.set_yscale('log')
    ax1.set_xticks(range(len(env_names)))
    ax1.set_xticklabels(env_names, rotation=45, ha='right')
    ax1.set_ylabel(r'$\varepsilon_{\rm eff}$ (dark photon mixing)', fontsize=11)
    ax1.set_title(r'Dark Photon Mixing Across Environments', fontsize=12)
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add value labels
    for i, (bar, val) in enumerate(zip(bars1, eps_values)):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height * 1.5,
                f'{val:.1e}', ha='center', va='bottom', fontsize=9)
    
    # Plot 2: Amplification factors
    bars2 = ax2.bar(range(len(env_names)), amplification, color=['blue', 'green', 'orange', 'red'], alpha=0.7)
    ax2.set_yscale('log')
    ax2.set_xticks(range(len(env_names)))
    ax2.set_xticklabels(env_names, rotation=45, ha='right')
    ax2.set_ylabel('Amplification Factor (relative to lab)', fontsize=11)
    ax2.set_title('Curvature Amplification Effect', fontsize=12)
    ax2.grid(True, alpha=0.3, axis='y')
    ax2.axhline(1, color='black', linestyle='--', alpha=0.5, linewidth=1)
    
    # Add amplification labels
    for i, (bar, amp) in enumerate(zip(bars2, amplification)):
        if amp > 1:
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height * 1.5,
                    f'{amp:.1e}×', ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    
    output_file = output_dir / 'curvature_amplification.pdf'
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Saved: {output_file}")
    
    fig.savefig(output_file.with_suffix('.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def main():
    """Generate all BSM plots."""
    print("=" * 70)
    print("BSM Parameter Space Plot Generation")
    print("=" * 70)
    
    # Create output directory
    output_dir = ROOT / 'papers' / 'kappaR_to_BSM' / 'figures'
    output_dir.mkdir(exist_ok=True)
    
    # Generate plots
    generate_epsilon_vs_R_plot(output_dir)
    generate_axion_vs_R_plot(output_dir)
    generate_amplification_plot(output_dir)
    
    print("\n" + "=" * 70)
    print("Plot generation complete!")
    print(f"Outputs in: {output_dir}")
    print("=" * 70)


if __name__ == '__main__':
    main()
