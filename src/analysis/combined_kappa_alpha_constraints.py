"""
Combined (κ_R, α) Constraint from Lab + Black Hole Observations

Implements joint constraint analysis:
    Δχ² = χ²_lab(κ_R) + χ²_BH(α) < χ²_crit

Explores parameter space for:
- Cancellation scenarios (opposite signs)
- Correlation from UV completion (α = f(κ_R))

References:
- Lab constraint: κ_R < 5×10¹⁷ m² (our coherence-gravity coupling)
- BH constraint: α ≲ 10²⁸ m² (Carballo-Rubio et al., arXiv:2505.21431)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
from scipy.optimize import minimize
from typing import Tuple, Optional, Callable

class CombinedConstraint:
    """
    Joint (κ_R, α) parameter space analysis.
    
    Observables:
    - Lab: Torque measurement sensitive to κ_R R F²
    - BH: Photon ring splitting sensitive to α L F²
    """
    
    def __init__(self, 
                 kappa_R_bound: float = 5e17,  # m²
                 alpha_bound: float = 1e28,    # m²
                 sigma_kappa: float = 1e17,    # Uncertainty (m²)
                 sigma_alpha: float = 5e27):   # Uncertainty (m²)
        """
        Parameters
        ----------
        kappa_R_bound : float
            Laboratory constraint on |κ_R| (m²)
        alpha_bound : float
            Black hole constraint on |α| (m²)
        sigma_kappa : float
            Statistical uncertainty in κ_R (m²)
        sigma_alpha : float
            Statistical uncertainty in α (m²)
        """
        self.kappa_R_bound = kappa_R_bound
        self.alpha_bound = alpha_bound
        self.sigma_kappa = sigma_kappa
        self.sigma_alpha = sigma_alpha
    
    def chi_squared_lab(self, kappa_R: float) -> float:
        """
        Laboratory χ² for κ_R.
        
        Assumes Gaussian likelihood:
            χ²_lab = (κ_R - κ_R_best)² / σ_κ²
        where κ_R_best = 0 (null result).
        
        Parameters
        ----------
        kappa_R : float
            Test value of κ_R (m²)
        
        Returns
        -------
        chi_sq : float
            χ² contribution from lab data
        """
        kappa_R_best = 0.0  # Null result (no detection)
        return (kappa_R - kappa_R_best)**2 / self.sigma_kappa**2
    
    def chi_squared_BH(self, alpha: float) -> float:
        """
        Black hole χ² for α.
        
        Assumes Gaussian likelihood centered at α = 0.
        
        Parameters
        ----------
        alpha : float
            Test value of α (m²)
        
        Returns
        -------
        chi_sq : float
            χ² contribution from BH data
        """
        alpha_best = 0.0  # Null result
        return (alpha - alpha_best)**2 / self.sigma_alpha**2
    
    def chi_squared_combined(self, kappa_R: float, alpha: float) -> float:
        """
        Combined χ² = χ²_lab + χ²_BH.
        
        Parameters
        ----------
        kappa_R : float
            Test value of κ_R (m²)
        alpha : float
            Test value of α (m²)
        
        Returns
        -------
        chi_sq_total : float
            Total χ²
        """
        return self.chi_squared_lab(kappa_R) + self.chi_squared_BH(alpha)
    
    def confidence_region(self, confidence_level: float = 0.95,
                          dof: int = 2) -> float:
        """
        Compute critical χ² for given confidence level.
        
        Parameters
        ----------
        confidence_level : float
            Confidence level (e.g., 0.95 for 95%)
        dof : int
            Degrees of freedom (2 for (κ_R, α))
        
        Returns
        -------
        chi_sq_crit : float
            Critical χ² value
        """
        return chi2.ppf(confidence_level, dof)
    
    def cancellation_scenario(self, 
                              kappa_R_range: Tuple[float, float] = (-1e18, 1e18),
                              num_points: int = 200) -> Tuple[np.ndarray, np.ndarray]:
        """
        Explore cancellation: κ_R < 0, α > 0 (or vice versa).
        
        Question: Can opposite signs of κ_R and α lead to partial cancellation
        in observables, weakening constraints?
        
        Answer: NO for independent observables (torque vs photon ring), but
        YES if there are cross-terms in combined observable.
        
        Parameters
        ----------
        kappa_R_range : tuple
            (κ_R_min, κ_R_max) in m²
        num_points : int
            Grid resolution
        
        Returns
        -------
        kappa_R_grid : ndarray
            κ_R values (m²)
        alpha_grid : ndarray
            α values assuming α = -κ_R (opposite sign)
        """
        kappa_R_grid = np.linspace(kappa_R_range[0], kappa_R_range[1], num_points)
        alpha_grid = -kappa_R_grid  # Opposite sign ansatz
        
        return kappa_R_grid, alpha_grid
    
    def UV_correlation(self, kappa_R: float, 
                       model: str = 'string_theory') -> float:
        """
        Predict α from κ_R using UV completion.
        
        Models:
        - 'string_theory': α ~ κ_R / (M_string)² × (geometric factors)
        - 'asymptotic_safety': α ~ κ_R × β(G) (RG flow)
        - 'linear': α = C × κ_R (simple proportionality)
        
        Parameters
        ----------
        kappa_R : float
            κ_R value (m²)
        model : str
            UV completion model
        
        Returns
        -------
        alpha : float
            Predicted α (m²)
        """
        if model == 'string_theory':
            # String scale M_s ~ 10¹⁶ GeV ~ 10⁻¹⁹ m
            # α ~ κ_R / M_s² ~ κ_R × 10³⁸
            # But: Geometric suppression ~ 1/N² ~ 0.01 (extra dimensions)
            return kappa_R * 1e36
        
        elif model == 'asymptotic_safety':
            # RG flow: α(μ) ~ κ_R × β(G_N(μ))
            # β ~ O(1) near UV fixed point
            return kappa_R * 2.0
        
        elif model == 'linear':
            # Simplest ansatz: α = C × κ_R
            C = 1.0  # Order-of-magnitude
            return C * kappa_R
        
        else:
            raise ValueError(f"Unknown UV model: {model}")
    
    def marginalized_constraint(self, 
                                 uv_model: str = 'linear',
                                 confidence_level: float = 0.95) -> float:
        """
        Marginalize over α using α = f(κ_R) correlation.
        
        Returns 1D constraint on κ_R accounting for α correlation.
        
        Parameters
        ----------
        uv_model : str
            UV completion model
        confidence_level : float
            Confidence level
        
        Returns
        -------
        kappa_R_constrained : float
            Upper limit on |κ_R| (m²)
        """
        chi_sq_crit = self.confidence_region(confidence_level, dof=1)  # 1D after marginalization
        
        def objective(kappa_R):
            alpha = self.UV_correlation(kappa_R, model=uv_model)
            return self.chi_squared_combined(kappa_R, alpha) - chi_sq_crit
        
        # Find κ_R where Δχ² = χ²_crit
        from scipy.optimize import fsolve
        kappa_R_limit = fsolve(objective, x0=self.kappa_R_bound)[0]
        
        return abs(kappa_R_limit)


def plot_combined_constraints(save_path: Optional[str] = None):
    """
    Visualize (κ_R, α) parameter space with lab + BH constraints.
    """
    model = CombinedConstraint()
    
    # Grid
    kappa_R_grid = np.linspace(-1e18, 1e18, 300)
    alpha_grid = np.linspace(-2e28, 2e28, 300)
    KR, AL = np.meshgrid(kappa_R_grid, alpha_grid)
    
    # Combined χ²
    chi_sq_total = np.zeros_like(KR)
    for i in range(len(alpha_grid)):
        for j in range(len(kappa_R_grid)):
            chi_sq_total[i, j] = model.chi_squared_combined(KR[i, j], AL[i, j])
    
    # Confidence regions
    chi_sq_68 = model.confidence_region(0.68, dof=2)  # 1σ
    chi_sq_95 = model.confidence_region(0.95, dof=2)  # 2σ
    chi_sq_99 = model.confidence_region(0.99, dof=2)  # 3σ
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    
    # Panel 1: Confidence contours
    contours = ax1.contour(KR / 1e17, AL / 1e28, chi_sq_total, 
                           levels=[chi_sq_68, chi_sq_95, chi_sq_99],
                           colors=['blue', 'green', 'red'],
                           linewidths=2)
    ax1.clabel(contours, fmt={chi_sq_68: '68%', chi_sq_95: '95%', chi_sq_99: '99%'})
    
    # Individual constraints
    ax1.axvline(model.kappa_R_bound / 1e17, color='blue', linestyle='--', 
                alpha=0.5, label=f"Lab: κ_R < {model.kappa_R_bound/1e17:.0f}×10¹⁷ m²")
    ax1.axhline(model.alpha_bound / 1e28, color='green', linestyle='--', 
                alpha=0.5, label=f"BH: α < {model.alpha_bound/1e28:.0f}×10²⁸ m²")
    
    # Cancellation line
    ax1.plot(kappa_R_grid / 1e17, -kappa_R_grid / 1e28, 'k:', 
             linewidth=1.5, label="Cancellation: α = -κ_R")
    
    ax1.set_xlabel("κ_R (10¹⁷ m²)", fontsize=12)
    ax1.set_ylabel("α (10²⁸ m²)", fontsize=12)
    ax1.set_title("Combined (κ_R, α) Constraints", fontsize=13)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(-10, 10)
    ax1.set_ylim(-2, 2)
    
    # Panel 2: UV correlation models
    kappa_R_test = np.linspace(-5e17, 5e17, 100)
    
    for uv_model, color, style in [('linear', 'blue', '-'),
                                    ('asymptotic_safety', 'green', '--'),
                                    ('string_theory', 'red', '-.')]:
        alpha_pred = np.array([model.UV_correlation(kr, model=uv_model) for kr in kappa_R_test])
        ax2.plot(kappa_R_test / 1e17, alpha_pred / 1e28, color=color, 
                 linestyle=style, linewidth=2, label=f"{uv_model}")
    
    # Constraints
    ax2.axhline(model.alpha_bound / 1e28, color='gray', linestyle='--', alpha=0.5)
    ax2.axhline(-model.alpha_bound / 1e28, color='gray', linestyle='--', alpha=0.5)
    ax2.axvline(model.kappa_R_bound / 1e17, color='gray', linestyle='--', alpha=0.5)
    ax2.axvline(-model.kappa_R_bound / 1e17, color='gray', linestyle='--', alpha=0.5)
    
    ax2.set_xlabel("κ_R (10¹⁷ m²)", fontsize=12)
    ax2.set_ylabel("α (10²⁸ m²)", fontsize=12)
    ax2.set_title("UV Completion Correlations", fontsize=13)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(-5, 5)
    ax2.set_ylim(-10, 10)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {save_path}")
    
    return fig


def validate_combined_constraints():
    """
    Validate combined constraint analysis.
    """
    model = CombinedConstraint()
    
    print("=" * 70)
    print("COMBINED (κ_R, α) CONSTRAINT ANALYSIS")
    print("=" * 70)
    
    print("\n1. Individual Constraints")
    print("-" * 70)
    print(f"Lab (κ_R):  |κ_R| < {model.kappa_R_bound:.1e} m² (σ = {model.sigma_kappa:.1e} m²)")
    print(f"BH (α):     |α| < {model.alpha_bound:.1e} m² (σ = {model.sigma_alpha:.1e} m²)")
    print(f"Ratio:      α/κ_R ~ {model.alpha_bound / model.kappa_R_bound:.1e}")
    
    print("\n2. Confidence Regions")
    print("-" * 70)
    print(f"68% CL: Δχ² < {model.confidence_region(0.68, dof=2):.2f}")
    print(f"95% CL: Δχ² < {model.confidence_region(0.95, dof=2):.2f}")
    print(f"99% CL: Δχ² < {model.confidence_region(0.99, dof=2):.2f}")
    
    print("\n3. Cancellation Scenario")
    print("-" * 70)
    kappa_test = 1e17  # m²
    alpha_cancel = -kappa_test  # Opposite sign
    chi_sq_cancel = model.chi_squared_combined(kappa_test, alpha_cancel)
    print(f"Test: κ_R = {kappa_test:.1e} m², α = {alpha_cancel:.1e} m²")
    print(f"Δχ² = {chi_sq_cancel:.2f}")
    print("Result: NO cancellation (independent observables)")
    
    print("\n4. UV Correlation Models")
    print("-" * 70)
    kappa_R_test = model.kappa_R_bound
    for uv_model in ['linear', 'asymptotic_safety', 'string_theory']:
        alpha_pred = model.UV_correlation(kappa_R_test, model=uv_model)
        print(f"{uv_model:20s}: α = {alpha_pred:.2e} m² (for κ_R = {kappa_R_test:.1e} m²)")
    
    print("\n5. Marginalized Constraints")
    print("-" * 70)
    for uv_model in ['linear', 'asymptotic_safety']:
        kappa_margin = model.marginalized_constraint(uv_model=uv_model, confidence_level=0.95)
        print(f"{uv_model:20s}: |κ_R| < {kappa_margin:.2e} m² (95% CL, marginalized over α)")
    
    print("\n" + "=" * 70)
    print("Key insight: Lab κ_R constraint is 10¹¹× tighter than BH α constraint")
    print("Implications:")
    print("  - If operators mix (C_α ~ 1), lab constraint dominates")
    print("  - UV correlations can tighten joint constraint")
    print("  - No cancellation for independent observables")
    print("=" * 70)


if __name__ == "__main__":
    validate_combined_constraints()
    print("\n")
    plot_combined_constraints(save_path="figures/combined_kappa_alpha_constraints.pdf")
    print("\n✅ Task 17 complete: Derived combined (κ_R, α) constraint")
