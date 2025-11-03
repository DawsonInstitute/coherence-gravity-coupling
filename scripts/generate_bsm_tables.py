#!/usr/bin/env python3
"""
Generate LaTeX tables for the BSM bounds paper from computed CSV data.
"""
from __future__ import annotations
import os
import csv

RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results', 'bsm_bounds')
PAPER_DIR = os.path.join(os.path.dirname(__file__), '..', 'papers', 'kappaR_to_BSM')


def read_epsilon_csv():
    """Read epsilon_equiv.csv and organize by environment and C_eps."""
    path = os.path.join(RESULTS_DIR, 'epsilon_equiv.csv')
    data = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            env = row['env']
            C_key = row['C_eps']
            kappa = row['kappa_m2']
            eps = float(row['epsilon'])
            if env not in data:
                data[env] = {}
            if C_key not in data[env]:
                data[env][C_key] = {}
            data[env][C_key][kappa] = eps
    return data


def generate_epsilon_table():
    """Generate LaTeX table for epsilon_eff values."""
    data = read_epsilon_csv()
    
    # Use a kappa value actually in the CSV
    kappa_val = "1.000e-11 m^2"  # Representative laboratory-scale bound
    
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Effective dark photon mixing $\varepsilon_{\rm eff}$ for $\kappa_R \sim 10^{-11}\,\mathrm{m}^2$ across curvature environments and matching coefficients.}")
    lines.append(r"\label{tab:epsilon}")
    lines.append(r"\begin{tabular}{lccc}")
    lines.append(r"\hline")
    lines.append(r"Environment & $\mathcal{R}\,[\mathrm{m}^{-2}]$ & $\varepsilon_{\rm eff}$ ($C_\varepsilon=1$) & $\varepsilon_{\rm eff}$ ($C_\varepsilon=10^{-2}$) \\")
    lines.append(r"\hline")
    
    env_labels = {
        'lab_flat': ('Lab (flat)', r'$10^{-30}$'),
        'earth_surface': ('Earth surface', r'$10^{-26}$'),
        'magnetar_surface': ('Magnetar', r'$10^{-6}$'),
    }
    
    for env_key, (label, R_str) in env_labels.items():
        if env_key in data:
            eps_1 = data[env_key].get('C_eps=1e+00', {}).get(kappa_val, 0)
            eps_2 = data[env_key].get('C_eps=1e-02', {}).get(kappa_val, 0)
            lines.append(f"{label} & {R_str} & {eps_1:.2e} & {eps_2:.2e} \\\\")
    
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    
    return "\n".join(lines)


def generate_axion_table():
    """Generate LaTeX table for axion g_equiv values."""
    path = os.path.join(RESULTS_DIR, 'g_axion_equiv.csv')
    
    data = {}
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            env = row['env']
            C_a = row['C_a']
            kappa = row['kappa_m2']
            g = float(row['g_equiv_GeVinv'])
            if env not in data:
                data[env] = {}
            if C_a not in data[env]:
                data[env][C_a] = {}
            data[env][C_a][kappa] = g
    
    kappa_val = "1.000e-11"  # Match the epsilon table
    
    lines = []
    lines.append(r"\begin{table}[htbp]")
    lines.append(r"\centering")
    lines.append(r"\caption{Parametric axion-equivalent coupling $g_{a\gamma\gamma}^{\rm equiv}$ for $\kappa_R \sim 10^{-11}\,\mathrm{m}^2$, $\Lambda=10\,\mathrm{TeV}$. \emph{Model-dependent; requires CP-odd portal.}}")
    lines.append(r"\label{tab:axion}")
    lines.append(r"\begin{tabular}{lccc}")
    lines.append(r"\hline")
    lines.append(r"Environment & $\mathcal{R}\,[\mathrm{m}^{-2}]$ & $g^{\rm equiv}_{a\gamma\gamma}\,[\mathrm{GeV}^{-1}]$ ($C_a=1$) & $g^{\rm equiv}_{a\gamma\gamma}\,[\mathrm{GeV}^{-1}]$ ($C_a=10^{-2}$) \\")
    lines.append(r"\hline")
    
    env_labels = {
        'lab_flat': ('Lab (flat)', r'$10^{-30}$'),
        'earth_surface': ('Earth surface', r'$10^{-26}$'),
        'magnetar_surface': ('Magnetar', r'$10^{-6}$'),
    }
    
    for env_key, (label, R_str) in env_labels.items():
        if env_key in data:
            g_1 = data[env_key].get('1e+00', {}).get(kappa_val, 0)
            g_2 = data[env_key].get('1e-02', {}).get(kappa_val, 0)
            lines.append(f"{label} & {R_str} & {g_1:.2e} & {g_2:.2e} \\\\")
    
    lines.append(r"\hline")
    lines.append(r"\end{tabular}")
    lines.append(r"\end{table}")
    
    return "\n".join(lines)


def main():
    os.makedirs(PAPER_DIR, exist_ok=True)
    
    eps_table = generate_epsilon_table()
    axion_table = generate_axion_table()
    
    # Write to separate .tex snippets
    with open(os.path.join(PAPER_DIR, 'table_epsilon.tex'), 'w') as f:
        f.write(eps_table)
    
    with open(os.path.join(PAPER_DIR, 'table_axion.tex'), 'w') as f:
        f.write(axion_table)
    
    print(f"Generated LaTeX tables in {PAPER_DIR}")
    print("  - table_epsilon.tex")
    print("  - table_axion.tex")


if __name__ == "__main__":
    main()
