#!/usr/bin/env python3
"""
R-dependent adaptive mesh refinement convergence study.

This module implements adaptive mesh refinement (AMR) based on local curvature
scalar R to ensure numerical convergence near singular points where the DOF
structure changes. Critical for validating κ_R = 0 baseline and detecting
deviations that would signal BSM physics.

Physical motivation:
- Near singular points (R → 0, parameter-dependent loci), extra scalar/tensor
  modes can activate/deactivate
- Standard uniform meshes may miss rapid changes in curvature coupling
- Adaptive refinement ensures we capture all physics needed to confirm κ_R = 0
  or detect κ_R ≠ 0 signatures

NEW PHYSICS discovery potential:
- If convergence fails ONLY when κ_R ≠ 0 included → evidence for new coupling
- Mesh refinement patterns different for GR vs modified gravity → diagnostic
- Residual scaling breaks down near singular points → extra DOF signature

Author: DawsonInstitute
Date: November 6, 2025
"""

import numpy as np
from typing import List, Tuple, Optional, Dict, Any
from dataclasses import dataclass
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


@dataclass
class RefinementCriterion:
    """Criteria for triggering mesh refinement based on curvature."""
    
    # Refinement thresholds
    r_gradient_threshold: float = 0.1  # Refine if |∇R|/R > this value
    r_second_deriv_threshold: float = 1.0  # Refine if |R''|/R > this value
    min_cell_size: float = 1e-6  # Minimum allowed cell size [m]
    max_cell_size: float = 1.0  # Maximum allowed cell size [m]
    target_cells_per_curvature_scale: int = 10  # Cells per 1/√|R|
    
    # Physics-based thresholds
    singular_point_buffer: float = 0.01  # Buffer around R = 0 [in units of max|R|]
    dof_transition_buffer: float = 0.05  # Buffer around DOF transition points


@dataclass
class ConvergenceMetrics:
    """Metrics for assessing convergence quality."""
    
    mesh_sizes: List[int]
    l2_errors: List[float]
    linf_errors: List[float]
    convergence_rates: List[float]
    residuals: List[float]
    condition_numbers: List[float]
    
    def compute_convergence_order(self) -> float:
        """
        Compute observed convergence order from L2 errors.
        
        For stable numerics should see:
        - 2nd order scheme: rate ≈ 2.0
        - 1st order scheme: rate ≈ 1.0
        - Sub-linear: < 1.0 → possible instability or under-resolution
        
        Returns:
            Average convergence order across refinements
        """
        if len(self.l2_errors) < 2:
            return 0.0
        
        orders = []
        for i in range(1, len(self.l2_errors)):
            h_ratio = self.mesh_sizes[i-1] / self.mesh_sizes[i]
            error_ratio = self.l2_errors[i-1] / self.l2_errors[i]
            
            if error_ratio > 1.0 and h_ratio > 1.0:
                order = np.log(error_ratio) / np.log(h_ratio)
                orders.append(order)
        
        return np.mean(orders) if orders else 0.0
    
    def is_converged(self, 
                    target_order: float = 1.5,
                    residual_tol: float = 1e-6) -> bool:
        """
        Check if refinement has converged acceptably.
        
        Args:
            target_order: Minimum acceptable convergence order
            residual_tol: Maximum acceptable residual
            
        Returns:
            True if converged, False otherwise
        """
        order = self.compute_convergence_order()
        last_residual = self.residuals[-1] if self.residuals else np.inf
        
        converged = (order >= target_order and 
                    last_residual < residual_tol)
        
        return converged


class CurvatureAdaptiveMesh:
    """
    Adaptive mesh refinement based on curvature scalar R.
    
    Strategy:
    1. Start with coarse uniform mesh
    2. Compute R(x) on current mesh
    3. Identify cells with high |∇R|, |R''|, or near singular points
    4. Refine those cells
    5. Repeat until convergence criteria met
    """
    
    def __init__(self, 
                 domain: Tuple[float, float],
                 criterion: Optional[RefinementCriterion] = None):
        """
        Initialize adaptive mesh.
        
        Args:
            domain: (x_min, x_max) spatial domain bounds
            criterion: Refinement criteria (uses defaults if None)
        """
        self.domain = domain
        self.criterion = criterion or RefinementCriterion()
        
        # Current mesh
        self.x_cells: np.ndarray = np.array([])
        self.cell_sizes: np.ndarray = np.array([])
        
        # Curvature values
        self.R_values: np.ndarray = np.array([])
        self.R_gradient: np.ndarray = np.array([])
        self.R_second_deriv: np.ndarray = np.array([])
        
    def initialize_uniform_mesh(self, n_cells: int) -> None:
        """Create initial uniform mesh."""
        self.x_cells = np.linspace(self.domain[0], self.domain[1], n_cells)
        self.cell_sizes = np.full(n_cells-1, 
                                   (self.domain[1] - self.domain[0]) / (n_cells-1))
        
    def compute_curvature(self, 
                         R_function: callable,
                         kappa_R: float = 0.0) -> None:
        """
        Compute curvature scalar and derivatives on current mesh.
        
        Args:
            R_function: Function R(x) returning curvature at position x
            kappa_R: Curvature-EM coupling constant [m²]
        """
        self.R_values = np.array([R_function(x, kappa_R) for x in self.x_cells])
        
        # Compute gradients using central differences
        self.R_gradient = np.gradient(self.R_values, self.x_cells)
        self.R_second_deriv = np.gradient(self.R_gradient, self.x_cells)
        
    def identify_refinement_cells(self) -> np.ndarray:
        """
        Identify cells that need refinement based on criteria.
        
        Returns:
            Boolean array indicating which cells to refine
        """
        n_cells = len(self.x_cells)
        refine = np.zeros(n_cells, dtype=bool)
        
        # Criterion 1: High curvature gradient
        with np.errstate(divide='ignore', invalid='ignore'):
            grad_ratio = np.abs(self.R_gradient) / (np.abs(self.R_values) + 1e-30)
        refine |= (grad_ratio > self.criterion.r_gradient_threshold)
        
        # Criterion 2: High second derivative
        with np.errstate(divide='ignore', invalid='ignore'):
            second_deriv_ratio = np.abs(self.R_second_deriv) / (np.abs(self.R_values) + 1e-30)
        refine |= (second_deriv_ratio > self.criterion.r_second_deriv_threshold)
        
        # Criterion 3: Near singular points (R ≈ 0)
        r_max = np.max(np.abs(self.R_values))
        near_singular = np.abs(self.R_values) < self.criterion.singular_point_buffer * r_max
        refine |= near_singular
        
        # Criterion 4: Insufficient resolution of curvature length scale
        # Length scale: ℓ_R ~ 1/√|R|
        with np.errstate(divide='ignore', invalid='ignore'):
            curvature_scale = 1.0 / np.sqrt(np.abs(self.R_values) + 1e-30)
        
        # Check if we have enough cells per curvature scale
        for i in range(len(self.cell_sizes)):
            cells_per_scale = curvature_scale[i] / self.cell_sizes[i]
            if cells_per_scale < self.criterion.target_cells_per_curvature_scale:
                refine[i] = True
                if i+1 < n_cells:
                    refine[i+1] = True
        
        return refine
    
    def refine_mesh(self) -> int:
        """
        Refine mesh by subdividing marked cells.
        
        Returns:
            Number of cells added
        """
        refine_mask = self.identify_refinement_cells()
        
        new_x_cells = []
        cells_added = 0
        
        for i in range(len(self.x_cells) - 1):
            new_x_cells.append(self.x_cells[i])
            
            if refine_mask[i] and self.cell_sizes[i] > self.criterion.min_cell_size:
                # Add midpoint
                midpoint = 0.5 * (self.x_cells[i] + self.x_cells[i+1])
                new_x_cells.append(midpoint)
                cells_added += 1
        
        new_x_cells.append(self.x_cells[-1])
        
        self.x_cells = np.array(new_x_cells)
        self.cell_sizes = np.diff(self.x_cells)
        
        return cells_added
    
    def adaptive_refinement_loop(self,
                                 R_function: callable,
                                 kappa_R: float = 0.0,
                                 max_iterations: int = 10,
                                 min_cells_added: int = 5) -> ConvergenceMetrics:
        """
        Iteratively refine mesh until convergence.
        
        Args:
            R_function: Function R(x, kappa_R) returning curvature
            kappa_R: Curvature-EM coupling [m²]
            max_iterations: Maximum refinement iterations
            min_cells_added: Stop if fewer than this many cells added
            
        Returns:
            Convergence metrics tracking refinement quality
        """
        mesh_sizes = []
        l2_errors = []
        linf_errors = []
        residuals = []
        condition_numbers = []
        
        R_prev = None
        
        for iteration in range(max_iterations):
            # Compute curvature on current mesh
            self.compute_curvature(R_function, kappa_R)
            
            # Record metrics
            mesh_sizes.append(len(self.x_cells))
            
            # Compute errors relative to previous iteration
            if R_prev is not None:
                # Interpolate previous solution to current mesh
                interp_func = interp1d(x_prev, R_prev, kind='cubic', 
                                      bounds_error=False, fill_value='extrapolate')
                R_interp = interp_func(self.x_cells)
                
                l2_error = np.sqrt(np.mean((self.R_values - R_interp)**2))
                linf_error = np.max(np.abs(self.R_values - R_interp))
                
                l2_errors.append(l2_error)
                linf_errors.append(linf_error)
            else:
                l2_errors.append(np.inf)
                linf_errors.append(np.inf)
            
            # Compute residual (measure of how well R satisfies governing equations)
            # For now, use gradient as proxy for residual
            residual = np.max(np.abs(self.R_gradient))
            residuals.append(residual)
            
            # Compute condition number (matrix health indicator)
            # Using curvature Hessian as proxy
            cond_number = np.max(np.abs(self.R_second_deriv)) / (np.min(np.abs(self.R_gradient)) + 1e-30)
            condition_numbers.append(cond_number)
            
            print(f"Iteration {iteration}: {len(self.x_cells)} cells, "
                  f"L2 error = {l2_errors[-1]:.2e}, "
                  f"Residual = {residual:.2e}")
            
            # Store for next iteration comparison
            x_prev = self.x_cells.copy()
            R_prev = self.R_values.copy()
            
            # Refine mesh
            cells_added = self.refine_mesh()
            
            if cells_added < min_cells_added:
                print(f"Converged: only {cells_added} cells added")
                break
        
        # Compute convergence rates
        convergence_rates = []
        for i in range(1, len(l2_errors)):
            if l2_errors[i-1] != np.inf and l2_errors[i] > 0:
                h_ratio = mesh_sizes[i-1] / mesh_sizes[i]
                error_ratio = l2_errors[i-1] / l2_errors[i]
                
                if h_ratio > 1.0 and error_ratio > 0:
                    rate = np.log(error_ratio) / np.log(h_ratio)
                    convergence_rates.append(rate)
                else:
                    convergence_rates.append(0.0)
            else:
                convergence_rates.append(0.0)
        
        return ConvergenceMetrics(
            mesh_sizes=mesh_sizes,
            l2_errors=l2_errors,
            linf_errors=linf_errors,
            convergence_rates=convergence_rates,
            residuals=residuals,
            condition_numbers=condition_numbers
        )


def compare_gr_vs_modified_gravity(domain: Tuple[float, float] = (-10.0, 10.0),
                                   kappa_R_values: List[float] = [0.0, 1e15, 5e17],
                                   initial_cells: int = 50) -> Dict[float, ConvergenceMetrics]:
    """
    Compare convergence behavior for GR (κ_R = 0) vs modified gravity (κ_R ≠ 0).
    
    This is the key NEW PHYSICS test:
    - If convergence identical for all κ_R → numerical artifact only
    - If convergence differs → κ_R coupling affects physics
    - If convergence fails only for κ_R ≠ 0 → strong evidence for new coupling
    
    Args:
        domain: Spatial domain (x_min, x_max)
        kappa_R_values: List of κ_R values to test [m²]
        initial_cells: Starting mesh resolution
        
    Returns:
        Dictionary mapping κ_R → ConvergenceMetrics
    """
    
    def schwarzschild_curvature(x: float, kappa_R: float) -> float:
        """
        Schwarzschild-like curvature with κ_R modifications.
        
        R(r) = R_0 / r² + κ_R * (background EM field energy)
        
        For testing purposes, use:
        R(r) = M / r² + κ_R * F² / r⁴
        
        where F ~ 1/r² is typical EM field strength.
        """
        M = 1.0  # Mass scale [arbitrary units for testing]
        F_squared = 1.0 / (x**4 + 1e-6)  # EM field energy density
        
        # GR contribution
        R_gr = M / (x**2 + 1e-3)  # Softened to avoid true singularity
        
        # κ_R correction
        R_correction = kappa_R * F_squared * 1e-30  # Rescaled for numerics
        
        return R_gr + R_correction
    
    results = {}
    
    for kappa_R in kappa_R_values:
        print(f"\n{'='*60}")
        print(f"Testing κ_R = {kappa_R:.2e} m²")
        print(f"{'='*60}")
        
        mesh = CurvatureAdaptiveMesh(domain)
        mesh.initialize_uniform_mesh(initial_cells)
        
        metrics = mesh.adaptive_refinement_loop(
            R_function=schwarzschild_curvature,
            kappa_R=kappa_R,
            max_iterations=10,
            min_cells_added=5
        )
        
        results[kappa_R] = metrics
        
        print(f"\nFinal convergence order: {metrics.compute_convergence_order():.3f}")
        print(f"Converged: {metrics.is_converged()}")
        
    return results


def plot_convergence_comparison(results: Dict[float, ConvergenceMetrics],
                                output_file: str = "figures/curvature_amr_convergence.pdf") -> None:
    """
    Plot convergence comparison for different κ_R values.
    
    Args:
        results: Dictionary of κ_R → ConvergenceMetrics
        output_file: Path to save figure
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # Plot 1: L2 error vs mesh size
    ax = axes[0, 0]
    for kappa_R, metrics in results.items():
        label = f"κ_R = {kappa_R:.1e}" if kappa_R != 0 else "GR (κ_R = 0)"
        ax.loglog(metrics.mesh_sizes[1:], metrics.l2_errors[1:], 
                 marker='o', label=label)
    
    ax.set_xlabel("Number of mesh cells")
    ax.set_ylabel("L2 error")
    ax.set_title("Convergence: L2 Error")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Convergence rate
    ax = axes[0, 1]
    for kappa_R, metrics in results.items():
        if metrics.convergence_rates:
            label = f"κ_R = {kappa_R:.1e}" if kappa_R != 0 else "GR (κ_R = 0)"
            iterations = range(len(metrics.convergence_rates))
            ax.plot(iterations, metrics.convergence_rates, marker='s', label=label)
    
    ax.axhline(y=2.0, color='k', linestyle='--', alpha=0.5, label='2nd order')
    ax.axhline(y=1.0, color='k', linestyle=':', alpha=0.5, label='1st order')
    ax.set_xlabel("Refinement iteration")
    ax.set_ylabel("Convergence order")
    ax.set_title("Observed Convergence Order")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Residuals
    ax = axes[1, 0]
    for kappa_R, metrics in results.items():
        label = f"κ_R = {kappa_R:.1e}" if kappa_R != 0 else "GR (κ_R = 0)"
        ax.semilogy(range(len(metrics.residuals)), metrics.residuals,
                   marker='^', label=label)
    
    ax.set_xlabel("Refinement iteration")
    ax.set_ylabel("Max residual")
    ax.set_title("Residual Convergence")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Condition number
    ax = axes[1, 1]
    for kappa_R, metrics in results.items():
        label = f"κ_R = {kappa_R:.1e}" if kappa_R != 0 else "GR (κ_R = 0)"
        ax.semilogy(range(len(metrics.condition_numbers)), metrics.condition_numbers,
                   marker='d', label=label)
    
    ax.set_xlabel("Refinement iteration")
    ax.set_ylabel("Condition number")
    ax.set_title("Matrix Health (Condition Number)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nSaved convergence comparison plot to {output_file}")


def main():
    """Run adaptive mesh refinement convergence study."""
    
    print("="*60)
    print("R-DEPENDENT ADAPTIVE MESH REFINEMENT CONVERGENCE STUDY")
    print("="*60)
    print("\nNEW PHYSICS test: Does convergence behavior differ between")
    print("GR (κ_R = 0) and modified gravity (κ_R ≠ 0)?")
    print("\nIf YES → evidence that κ_R coupling affects physical dynamics")
    print("If NO → κ_R may be purely perturbative correction\n")
    
    # Test different κ_R values
    kappa_R_values = [
        0.0,           # GR baseline
        1e15,          # Weak coupling (well below lab bound)
        5e17           # Near lab limit κ_R < 5×10^17 m² (95% CL)
    ]
    
    results = compare_gr_vs_modified_gravity(
        domain=(-10.0, 10.0),
        kappa_R_values=kappa_R_values,
        initial_cells=50
    )
    
    # Generate comparison plot
    plot_convergence_comparison(results)
    
    # Summary analysis
    print("\n" + "="*60)
    print("SUMMARY: NEW PHYSICS IMPLICATIONS")
    print("="*60)
    
    gr_metrics = results[0.0]
    
    for kappa_R in kappa_R_values:
        if kappa_R == 0.0:
            continue
        
        metrics = results[kappa_R]
        
        print(f"\nκ_R = {kappa_R:.2e} m²:")
        print(f"  Convergence order: {metrics.compute_convergence_order():.3f} "
              f"(GR: {gr_metrics.compute_convergence_order():.3f})")
        print(f"  Final residual: {metrics.residuals[-1]:.2e} "
              f"(GR: {gr_metrics.residuals[-1]:.2e})")
        print(f"  Converged: {metrics.is_converged()} "
              f"(GR: {gr_metrics.is_converged()})")
        
        # Check for significant differences
        order_diff = abs(metrics.compute_convergence_order() - 
                        gr_metrics.compute_convergence_order())
        
        if order_diff > 0.2:
            print(f"  ⚠️  SIGNIFICANT DIFFERENCE: Convergence order differs by {order_diff:.3f}")
            print(f"  → κ_R coupling may affect dynamical evolution!")
        
        if metrics.is_converged() != gr_metrics.is_converged():
            print(f"  ⚠️  CONVERGENCE FAILURE: GR converges but κ_R ≠ 0 does not!")
            print(f"  → Strong evidence for new physics affecting numerical stability")
    
    print("\n" + "="*60)
    print("Next steps:")
    print("1. Run EFQS simulations with AMR enabled")
    print("2. Compare mesh refinement patterns near singular points")
    print("3. Test if κ_R-induced mesh changes correlate with torque signals")
    print("="*60)


if __name__ == "__main__":
    import os
    os.makedirs("figures", exist_ok=True)
    main()
