#!/usr/bin/env python3
"""
Unified Report Generator for Coherence-Gravity Coupling Analysis

Consolidates all sweep results into publication-ready tables:
- CSV format for data analysis
- Markdown format for documentation
- LaTeX format for manuscripts

Usage:
    python generate_report.py --format csv
    python generate_report.py --format markdown --output report.md
    python generate_report.py --format latex --output tables.tex
    python generate_report.py --all  # Generate all formats

Author: GitHub Copilot (Claude Sonnet 4.5)
Date: October 2025
License: MIT
"""

import argparse
import json
import csv
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime
import glob


def load_latest_results(results_dir: Path, pattern: str) -> Dict[str, Any]:
    """Load the most recent result file matching pattern.
    
    Args:
        results_dir: Directory containing result files
        pattern: Glob pattern to match files (e.g., 'xi_sweep_*.json')
    
    Returns:
        Dict containing the result data, or None if not found
    """
    files = sorted(glob.glob(str(results_dir / pattern)), reverse=True)
    if not files:
        return None
    
    with open(files[0]) as f:
        return json.load(f)


def generate_csv_xi_sweep(results: Dict, output: Path):
    """Generate CSV table for ξ parameter sweep."""
    if not results or 'data' not in results:
        print(f"⚠️  No ξ sweep data found")
        return
    
    with open(output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['xi', 'delta_tau_Nm', 'delta_G_over_G', 'tau_coherent_Nm', 
                         'tau_newtonian_Nm', 'solve_time_s', 'elapsed_time_s'])
        
        for key, val in results['data'].items():
            if key.startswith('xi_'):
                writer.writerow([
                    val['xi'],
                    val['delta_tau'],
                    val['delta_G_over_G'],
                    val['tau_coherent'],
                    val['tau_newtonian'],
                    val['solve_time'],
                    val['elapsed_time']
                ])
    
    print(f"✅ CSV: {output}")


def generate_csv_materials(results: Dict, output: Path):
    """Generate CSV table for material comparison."""
    if not results or 'data' not in results:
        print(f"⚠️  No material comparison data found")
        return
    
    with open(output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['material', 'Phi0_m_inv', 'delta_tau_Nm', 'delta_G_over_G', 
                         'elapsed_time_s'])
        
        for key, val in results['data'].items():
            writer.writerow([
                val['name'],
                val['Phi0'],
                val['delta_tau'],
                val['delta_G_over_G'],
                val['elapsed_time']
            ])
    
    print(f"✅ CSV: {output}")


def generate_csv_curvature(results: Dict, output: Path, sweep_type: str):
    """Generate CSV table for curvature coupling sweeps."""
    if not results or 'data' not in results:
        print(f"⚠️  No {sweep_type} data found")
        return
    
    with open(output, 'w', newline='') as f:
        writer = csv.writer(f)
        
        # Determine header based on sweep type
        if 'B_' in list(results['data'].keys())[0]:
            writer.writerow(['B_T', 'R_m_minus_2', 'E_V_per_m', 'precision', 
                             'kappa_limit_m_sq', 'F_squared_T_sq'])
            param_key = 'B'
        elif 'R_' in list(results['data'].keys())[0]:
            writer.writerow(['R_m_minus_2', 'B_T', 'E_V_per_m', 'precision', 
                             'kappa_limit_m_sq', 'F_squared_T_sq'])
            param_key = 'R'
        else:  # precision sweep
            writer.writerow(['precision', 'B_T', 'R_m_minus_2', 'E_V_per_m', 
                             'kappa_limit_m_sq', 'F_squared_T_sq'])
            param_key = 'precision'
        
        for key, val in results['data'].items():
            writer.writerow([
                val[param_key],
                val['R'] if param_key != 'R' else val['B'],
                val['E'],
                val['precision'] if param_key != 'precision' else val['B'],
                val['kappa_limit'],
                val['F_squared']
            ])
    
    print(f"✅ CSV: {output}")


def generate_markdown_table(results: Dict, sweep_type: str) -> str:
    """Generate Markdown table from results."""
    if not results or 'data' not in results:
        return f"*No {sweep_type} data available*\n"
    
    if sweep_type == 'xi_sweep':
        lines = [
            "| ξ | Δτ [N·m] | ΔG/G | τ_coh [N·m] | τ_newt [N·m] | Solve Time [s] |",
            "|---|----------|------|-------------|--------------|----------------|"
        ]
        for key, val in results['data'].items():
            if key.startswith('xi_'):
                lines.append(
                    f"| {val['xi']:.1f} | {val['delta_tau']:.3e} | {val['delta_G_over_G']:.3e} | "
                    f"{val['tau_coherent']:.3e} | {val['tau_newtonian']:.3e} | {val['solve_time']:.2f} |"
                )
    
    elif sweep_type == 'materials':
        lines = [
            "| Material | Φ₀ [m⁻¹] | Δτ [N·m] | ΔG/G | Time [s] |",
            "|----------|----------|----------|------|----------|"
        ]
        for key, val in results['data'].items():
            lines.append(
                f"| {val['name']} | {val['Phi0']:.2e} | {val['delta_tau']:.3e} | "
                f"{val['delta_G_over_G']:.3e} | {val['elapsed_time']:.3f} |"
            )
    
    elif 'curvature' in sweep_type:
        # Determine parameter being swept
        first_key = list(results['data'].keys())[0]
        if 'B_' in first_key:
            lines = [
                "| B [T] | R [m⁻²] | Precision | κ_R limit [m²] | F² [T²] |",
                "|-------|---------|-----------|----------------|---------|"
            ]
            for key, val in results['data'].items():
                lines.append(
                    f"| {val['B']:.1f} | {val['R']:.2e} | {val['precision']:.1e} | "
                    f"{val['kappa_limit']:.2e} | {val['F_squared']:.1f} |"
                )
        elif 'R_' in first_key:
            lines = [
                "| R [m⁻²] | B [T] | Precision | κ_R limit [m²] | F² [T²] |",
                "|---------|-------|-----------|----------------|---------|"
            ]
            for key, val in results['data'].items():
                lines.append(
                    f"| {val['R']:.2e} | {val['B']:.1f} | {val['precision']:.1e} | "
                    f"{val['kappa_limit']:.2e} | {val['F_squared']:.1f} |"
                )
        else:  # precision sweep
            lines = [
                "| Precision | B [T] | R [m⁻²] | κ_R limit [m²] | F² [T²] |",
                "|-----------|-------|---------|----------------|---------|"
            ]
            for key, val in results['data'].items():
                lines.append(
                    f"| {val['precision']:.1e} | {val['B']:.1f} | {val['R']:.2e} | "
                    f"{val['kappa_limit']:.2e} | {val['F_squared']:.1f} |"
                )
    
    return '\n'.join(lines) + '\n'


def generate_latex_table(results: Dict, sweep_type: str) -> str:
    """Generate LaTeX table from results."""
    if not results or 'data' not in results:
        return f"% No {sweep_type} data available\n"
    
    if sweep_type == 'xi_sweep':
        lines = [
            "\\begin{table}[ht]",
            "\\centering",
            "\\caption{ξ Parameter Sweep Results}",
            "\\begin{tabular}{cccccc}",
            "\\hline",
            "$\\xi$ & $\\Delta\\tau$ [N·m] & $\\Delta G/G$ & $\\tau_{\\text{coh}}$ [N·m] & $\\tau_{\\text{newt}}$ [N·m] & Solve Time [s] \\\\",
            "\\hline"
        ]
        for key, val in results['data'].items():
            if key.startswith('xi_'):
                lines.append(
                    f"{val['xi']:.1f} & {val['delta_tau']:.2e} & {val['delta_G_over_G']:.2e} & "
                    f"{val['tau_coherent']:.2e} & {val['tau_newtonian']:.2e} & {val['solve_time']:.1f} \\\\"
                )
        lines.extend(["\\hline", "\\end{tabular}", "\\label{tab:xi_sweep}", "\\end{table}"])
    
    elif sweep_type == 'materials':
        lines = [
            "\\begin{table}[ht]",
            "\\centering",
            "\\caption{Material Comparison (ξ=100)}",
            "\\begin{tabular}{lcccc}",
            "\\hline",
            "Material & $\\Phi_0$ [m$^{-1}$] & $\\Delta\\tau$ [N·m] & $\\Delta G/G$ & Time [s] \\\\",
            "\\hline"
        ]
        for key, val in results['data'].items():
            mat_name = val['name'].replace('_', '\\_')
            lines.append(
                f"{mat_name} & {val['Phi0']:.2e} & {val['delta_tau']:.2e} & "
                f"{val['delta_G_over_G']:.2e} & {val['elapsed_time']:.2f} \\\\"
            )
        lines.extend(["\\hline", "\\end{tabular}", "\\label{tab:materials}", "\\end{table}"])
    
    elif 'curvature' in sweep_type:
        first_key = list(results['data'].keys())[0]
        if 'B_' in first_key:
            lines = [
                "\\begin{table}[ht]",
                "\\centering",
                "\\caption{Curvature-EM Coupling Exclusion Limits vs. Magnetic Field}",
                "\\begin{tabular}{ccccc}",
                "\\hline",
                "$B$ [T] & $R$ [m$^{-2}$] & Precision & $\\kappa_R$ limit [m$^2$] & $F^2$ [T$^2$] \\\\",
                "\\hline"
            ]
            for key, val in results['data'].items():
                lines.append(
                    f"{val['B']:.1f} & {val['R']:.1e} & {val['precision']:.1e} & "
                    f"{val['kappa_limit']:.2e} & {val['F_squared']:.1f} \\\\"
                )
        elif 'R_' in first_key:
            lines = [
                "\\begin{table}[ht]",
                "\\centering",
                "\\caption{Curvature-EM Coupling Exclusion Limits vs. Ricci Scalar}",
                "\\begin{tabular}{ccccc}",
                "\\hline",
                "$R$ [m$^{-2}$] & $B$ [T] & Precision & $\\kappa_R$ limit [m$^2$] & $F^2$ [T$^2$] \\\\",
                "\\hline"
            ]
            for key, val in results['data'].items():
                lines.append(
                    f"{val['R']:.1e} & {val['B']:.1f} & {val['precision']:.1e} & "
                    f"{val['kappa_limit']:.2e} & {val['F_squared']:.1f} \\\\"
                )
        else:  # precision sweep
            lines = [
                "\\begin{table}[ht]",
                "\\centering",
                "\\caption{Curvature-EM Coupling Exclusion Limits vs. Experimental Precision}",
                "\\begin{tabular}{ccccc}",
                "\\hline",
                "Precision & $B$ [T] & $R$ [m$^{-2}$] & $\\kappa_R$ limit [m$^2$] & $F^2$ [T$^2$] \\\\",
                "\\hline"
            ]
            for key, val in results['data'].items():
                lines.append(
                    f"{val['precision']:.1e} & {val['B']:.1f} & {val['R']:.1e} & "
                    f"{val['kappa_limit']:.2e} & {val['F_squared']:.1f} \\\\"
                )
        
        lines.extend(["\\hline", "\\end{tabular}", f"\\label{{tab:{sweep_type}}}", "\\end{table}"])
    
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(
        description="Generate unified reports from analysis results"
    )
    parser.add_argument('--format', choices=['csv', 'markdown', 'latex'], 
                        help='Output format')
    parser.add_argument('--output', type=Path, 
                        help='Output file path (default: auto-generated)')
    parser.add_argument('--all', action='store_true',
                        help='Generate all formats')
    parser.add_argument('--results-dir', type=Path, 
                        default=Path('results/analysis'),
                        help='Directory containing result files')
    
    args = parser.parse_args()
    
    if not args.format and not args.all:
        parser.print_help()
        return
    
    results_dir = args.results_dir
    if not results_dir.exists():
        print(f"❌ Results directory not found: {results_dir}")
        return
    
    # Load latest results
    xi_results = load_latest_results(results_dir, 'xi_sweep_*.json')
    mat_results = load_latest_results(results_dir, 'material_comparison_*.json')
    curv_B_results = load_latest_results(results_dir, 'curvature_limits_2*.json')
    curv_R_results = load_latest_results(results_dir, 'curvature_limits_vs_R_*.json')
    curv_P_results = load_latest_results(results_dir, 'curvature_limits_vs_precision_*.json')
    
    # Generate requested format(s)
    formats = ['csv', 'markdown', 'latex'] if args.all else [args.format]
    
    for fmt in formats:
        if fmt == 'csv':
            # Generate CSV files
            csv_dir = Path('results/reports/csv')
            csv_dir.mkdir(parents=True, exist_ok=True)
            
            generate_csv_xi_sweep(xi_results, csv_dir / 'xi_sweep.csv')
            generate_csv_materials(mat_results, csv_dir / 'materials.csv')
            generate_csv_curvature(curv_B_results, csv_dir / 'curvature_vs_B.csv', 'vs_B')
            generate_csv_curvature(curv_R_results, csv_dir / 'curvature_vs_R.csv', 'vs_R')
            generate_csv_curvature(curv_P_results, csv_dir / 'curvature_vs_precision.csv', 'vs_precision')
            
            print(f"\n📊 CSV files saved to: {csv_dir}")
        
        elif fmt == 'markdown':
            # Generate Markdown report
            md_path = args.output or Path('results/reports/analysis_report.md')
            md_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(md_path, 'w') as f:
                f.write("# Coherence-Gravity Coupling Analysis Report\n\n")
                f.write(f"**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                f.write("---\n\n")
                
                f.write("## ξ Parameter Sweep\n\n")
                f.write(generate_markdown_table(xi_results, 'xi_sweep'))
                f.write("\n")
                
                f.write("## Material Comparison\n\n")
                f.write(generate_markdown_table(mat_results, 'materials'))
                f.write("\n")
                
                f.write("## Curvature-EM Coupling: Exclusion Limits vs. B-field\n\n")
                f.write(generate_markdown_table(curv_B_results, 'curvature_B'))
                f.write("\n")
                
                if curv_R_results:
                    f.write("## Curvature-EM Coupling: Exclusion Limits vs. Ricci Scalar\n\n")
                    f.write(generate_markdown_table(curv_R_results, 'curvature_R'))
                    f.write("\n")
                
                if curv_P_results:
                    f.write("## Curvature-EM Coupling: Exclusion Limits vs. Precision\n\n")
                    f.write(generate_markdown_table(curv_P_results, 'curvature_precision'))
                    f.write("\n")
            
            print(f"\n📝 Markdown report saved to: {md_path}")
        
        elif fmt == 'latex':
            # Generate LaTeX tables
            tex_path = args.output or Path('results/reports/analysis_tables.tex')
            tex_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(tex_path, 'w') as f:
                f.write("% Coherence-Gravity Coupling Analysis Tables\n")
                f.write(f"% Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                f.write("% ξ Parameter Sweep\n")
                f.write(generate_latex_table(xi_results, 'xi_sweep'))
                f.write("\n")
                
                f.write("% Material Comparison\n")
                f.write(generate_latex_table(mat_results, 'materials'))
                f.write("\n")
                
                f.write("% Curvature-EM Coupling vs. B-field\n")
                f.write(generate_latex_table(curv_B_results, 'curvature_B'))
                f.write("\n")
                
                if curv_R_results:
                    f.write("% Curvature-EM Coupling vs. Ricci Scalar\n")
                    f.write(generate_latex_table(curv_R_results, 'curvature_R'))
                    f.write("\n")
                
                if curv_P_results:
                    f.write("% Curvature-EM Coupling vs. Precision\n")
                    f.write(generate_latex_table(curv_P_results, 'curvature_precision'))
                    f.write("\n")
            
            print(f"\n📄 LaTeX tables saved to: {tex_path}")
    
    print("\n✅ Report generation complete!")


if __name__ == '__main__':
    main()