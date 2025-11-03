"""
Robin Boundary Condition Parametric Study

Physical Motivation (Gorkavenko et al. arxiv:2409.04647):
    - Casimir energy E_ξ depends on vacuum coupling ξ
    - Robin BC parameter θ ∈ [-π/2, π/2] interpolates Dirichlet (θ=0) ↔ Neumann (θ=-π/2)
    - Boundary amplification: ξ sensitivity enhanced by BC modifications

Study Design:
    Sweep (mr_0, θ, F, ξ) parameter space:
        - mr_0: Mass × cavity size (dimensionless)
        - θ: Robin BC angle
        - F: External field strength
        - ξ: Vacuum coupling parameter
    
    Output: Contour plots of E_ξ(θ, ξ) showing boundary-amplified sensitivity

Implementation:
    1. Import robin_bc_poisson.py solver
    2. Define parameter grid
    3. Compute E_ξ for each (mr_0, θ, F, ξ)
    4. Generate contour plots with matplotlib

Usage:
    python scripts/robin_bc_parametric_sweep.py --output results/robin_bc_sweep.pdf

References:
    - Gorkavenko et al. (2024) arXiv:2409.04647 (Table 1)
    - Task: Run parametric study (task 12)
"""
from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import sys
import os

# Import robin_bc_poisson solver directly without going through __init__
coherence_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
robin_bc_path = os.path.join(coherence_path, 'src', 'solvers', 'robin_bc_poisson.py')

# Use importlib to load module without triggering __init__.py
import importlib.util
spec = importlib.util.spec_from_file_location("robin_bc_poisson", robin_bc_path)
robin_bc_module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(robin_bc_module)

solve_poisson_robin = robin_bc_module.solve_poisson_robin
compute_xi_dependent_energy = robin_bc_module.compute_xi_dependent_energy
validate_robin_bc_limits = robin_bc_module.validate_robin_bc_limits


def parametric_sweep_robin_bc(
    mr_0_values: np.ndarray,
    theta_values: np.ndarray,
    F_values: np.ndarray,
    xi_values: np.ndarray,
    grid_size: int = 50
) -> dict:
    """Perform 4D parameter sweep over (mr_0, θ, F, ξ).
    
    Args:
        mr_0_values: Dimensionless mass × radius
        theta_values: Robin BC parameter [radians]
        F_values: External field strength [arb units]
        xi_values: Vacuum coupling parameter
        grid_size: Spatial grid resolution
    
    Returns:
        dict with keys:
            'E_xi': 4D array of energies [mr_0, theta, F, xi]
            'parameters': Dict of parameter arrays
    """
    n_mr = len(mr_0_values)
    n_theta = len(theta_values)
    n_F = len(F_values)
    n_xi = len(xi_values)
    
    E_xi_grid = np.zeros((n_mr, n_theta, n_F, n_xi))
    
    print(f"Running sweep: {n_mr}×{n_theta}×{n_F}×{n_xi} = {n_mr*n_theta*n_F*n_xi} points")
    
    total = n_mr * n_theta * n_F * n_xi
    count = 0
    
    for i, mr_0 in enumerate(mr_0_values):
        for j, theta in enumerate(theta_values):
            for k, F in enumerate(F_values):
                for l, xi in enumerate(xi_values):
                    count += 1
                    if count % max(1, total // 20) == 0:
                        print(f"  Progress: {count}/{total} ({100*count/total:.1f}%)")
                    
                    # Compute E_ξ for this parameter point
                    r_0 = 1.0  # Fixed cavity size (units: m)
                    E = compute_xi_dependent_energy(
                        mr_0, r_0, theta, xi, F,
                        grid_size=grid_size
                    )
                    
                    E_xi_grid[i, j, k, l] = E
    
    return {
        'E_xi': E_xi_grid,
        'parameters': {
            'mr_0': mr_0_values,
            'theta': theta_values,
            'F': F_values,
            'xi': xi_values
        }
    }


def plot_E_xi_contours(
    sweep_result: dict,
    fixed_params: dict,
    output_path: Path
):
    """Generate contour plots of E_ξ(θ, ξ) for fixed mr_0, F.
    
    Args:
        sweep_result: Output from parametric_sweep_robin_bc
        fixed_params: {'mr_0_idx': 0, 'F_idx': 0} to select slices
        output_path: Save location for figure
    """
    E_xi = sweep_result['E_xi']
    params = sweep_result['parameters']
    
    i_mr = fixed_params.get('mr_0_idx', 0)
    k_F = fixed_params.get('F_idx', 0)
    
    # Extract 2D slice: E_ξ(θ, ξ)
    E_slice = E_xi[i_mr, :, k_F, :]  # Shape: (n_theta, n_xi)
    
    theta = params['theta']
    xi = params['xi']
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Contour plot
    ax = axes[0]
    THETA, XI = np.meshgrid(theta, xi, indexing='ij')
    contour = ax.contourf(THETA * 180/np.pi, XI, E_slice, levels=20, cmap='viridis')
    ax.set_xlabel('Robin BC angle θ [degrees]')
    ax.set_ylabel('Vacuum coupling ξ')
    ax.set_title(f'E_ξ(θ, ξ) | mr_0={params["mr_0"][i_mr]:.2f}, F={params["F"][k_F]:.1e}')
    plt.colorbar(contour, ax=ax, label='Energy [arb units]')
    
    # ξ-dependence curves for different θ
    ax = axes[1]
    theta_samples = [0, -np.pi/6, -np.pi/4, -np.pi/2]  # Dirichlet → Neumann
    colors = ['blue', 'green', 'orange', 'red']
    
    for theta_val, color in zip(theta_samples, colors):
        idx = np.argmin(np.abs(theta - theta_val))
        ax.plot(xi, E_slice[idx, :], color=color, marker='o', 
                label=f'θ = {theta_val*180/np.pi:.0f}°')
    
    ax.set_xlabel('Vacuum coupling ξ')
    ax.set_ylabel('Energy E_ξ')
    ax.set_title('ξ-dependence for different BC')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n✅ Saved figure: {output_path}")
    plt.close()


def validate_gorkavenko_scaling(sweep_result: dict) -> dict:
    """Check if E_ξ ∝ (1/4 - ξ) scaling holds.
    
    From Gorkavenko Table 1: E_ξ should show linear dependence on (1/4 - ξ).
    
    Returns:
        dict with validation metrics
    """
    E_xi = sweep_result['E_xi']
    xi = sweep_result['parameters']['xi']
    
    # Average over (mr_0, θ, F) to get mean E_ξ(ξ)
    E_mean = np.mean(E_xi, axis=(0, 1, 2))
    
    # Fit to E = a * (1/4 - ξ) + b
    quarter_minus_xi = 0.25 - xi
    
    # Linear regression
    A = np.vstack([quarter_minus_xi, np.ones(len(xi))]).T
    coeffs = np.linalg.lstsq(A, E_mean, rcond=None)[0]
    a, b = coeffs
    
    # Predicted vs actual
    E_predicted = a * quarter_minus_xi + b
    residuals = E_mean - E_predicted
    R_squared = 1 - np.sum(residuals**2) / np.sum((E_mean - np.mean(E_mean))**2)
    
    return {
        'scaling_coefficient': a,
        'offset': b,
        'R_squared': R_squared,
        'formula': f'E_ξ ≈ {a:.3e} × (1/4 - ξ) + {b:.3e}',
        'validation': '✅ Linear scaling confirmed' if R_squared > 0.95 else '⚠️ Deviations from linear'
    }


# CLI interface
def main():
    parser = argparse.ArgumentParser(description='Robin BC parametric sweep')
    parser.add_argument('--output', type=str, default='results/robin_bc_sweep.pdf',
                        help='Output figure path')
    parser.add_argument('--grid-size', type=int, default=30,
                        help='Spatial grid resolution')
    args = parser.parse_args()
    
    print("Robin BC Parametric Sweep")
    print("=" * 60)
    
    # Define parameter ranges
    mr_0_values = np.array([0.1, 0.5, 1.0])  # Light → moderate → heavy field
    theta_values = np.linspace(0, -np.pi/2, 8)  # Dirichlet → Neumann
    F_values = np.array([0.0, 1.0])  # No field / with field
    xi_values = np.array([0.0, 0.1, 0.25, 0.5, 1.0, 10.0])  # Range from Gorkavenko
    
    # Run sweep
    print("\nRunning parametric sweep...")
    result = parametric_sweep_robin_bc(
        mr_0_values, theta_values, F_values, xi_values,
        grid_size=args.grid_size
    )
    
    # Validate scaling
    print("\nValidating Gorkavenko (1/4 - ξ) scaling...")
    validation = validate_gorkavenko_scaling(result)
    print(f"  {validation['formula']}")
    print(f"  R² = {validation['R_squared']:.4f}")
    print(f"  {validation['validation']}")
    
    # Generate plots
    print("\nGenerating contour plots...")
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    plot_E_xi_contours(
        result,
        fixed_params={'mr_0_idx': 1, 'F_idx': 0},  # mr_0=0.5, F=0
        output_path=output_path
    )
    
    print("\n✅ Robin BC parametric study complete")
    print(f"   Results: {output_path}")


if __name__ == "__main__":
    main()
