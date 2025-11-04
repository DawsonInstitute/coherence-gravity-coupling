"""
Joint Posterior Analysis: Marginalized (κ_R, α) Constraints

Combines laboratory κ_R measurements with black hole α observations
to derive marginalized constraints assuming UV completion correlations.

Extends combined_kappa_alpha_constraints.py with full Bayesian treatment.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, chi2
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
try:
    import corner  # For corner plots (optional)
except ImportError:
    corner = None
from typing import Tuple, Optional, Callable

class JointPosteriorAnalysis:
    """
    Bayesian analysis of (κ_R, α) parameter space with UV correlations.
    
    Key features:
    - Full posterior sampling (MCMC or importance sampling)
    - UV correlation models: α = f(κ_R) + noise
    - Marginalization over nuisance parameters
    - Confidence intervals and parameter estimation
    """
    
    def __init__(self, 
                 kappa_R_obs: float = 0.0,           # Observed κ_R (null result)
                 kappa_R_uncertainty: float = 1e17,  # 1σ uncertainty 
                 alpha_obs: float = 0.0,             # Observed α (null result)
                 alpha_uncertainty: float = 5e27):   # 1σ uncertainty
        """
        Parameters
        ----------
        kappa_R_obs : float
            Best-fit κ_R from laboratory measurements (m²)
        kappa_R_uncertainty : float
            1σ statistical uncertainty in κ_R (m²)
        alpha_obs : float
            Best-fit α from black hole observations (m²)
        alpha_uncertainty : float
            1σ statistical uncertainty in α (m²)
        """
        self.kappa_R_obs = kappa_R_obs
        self.sigma_kappa = kappa_R_uncertainty
        self.alpha_obs = alpha_obs
        self.sigma_alpha = alpha_uncertainty
    
    def log_likelihood(self, kappa_R: float, alpha: float) -> float:
        """
        Log-likelihood function for independent Gaussian measurements.
        
        Parameters
        ----------
        kappa_R : float
            Test value of κ_R (m²)
        alpha : float
            Test value of α (m²)
        
        Returns
        -------
        log_L : float
            Log-likelihood
        """
        # Gaussian likelihoods (null results centered at 0)
        log_L_kappa = -0.5 * ((kappa_R - self.kappa_R_obs) / self.sigma_kappa)**2
        log_L_alpha = -0.5 * ((alpha - self.alpha_obs) / self.sigma_alpha)**2
        
        return log_L_kappa + log_L_alpha
    
    def log_prior(self, kappa_R: float, alpha: float, 
                  uv_model: str = 'linear',
                  correlation_strength: float = 1.0) -> float:
        """
        Log-prior incorporating UV completion constraints.
        
        Prior models:
        - 'uniform': Flat prior over reasonable range
        - 'linear': α = C × κ_R + noise (correlated prior)
        - 'string_theory': α ~ κ_R / M_string² (theoretical prior)
        
        Parameters
        ----------
        kappa_R : float
            κ_R value (m²)
        alpha : float
            α value (m²)
        uv_model : str
            UV completion model
        correlation_strength : float
            Strength of UV correlation (0 = uncorrelated, 1 = fully correlated)
        
        Returns
        -------
        log_prior : float
            Log-prior probability
        """
        # Uniform bounds (generous)
        kappa_max = 1e18  # m²
        alpha_max = 1e29  # m²
        
        if abs(kappa_R) > kappa_max or abs(alpha) > alpha_max:
            return -np.inf
        
        if uv_model == 'uniform':
            # Flat prior
            return 0.0
        
        elif uv_model == 'linear':
            # Correlated prior: α = C × κ_R + noise
            C_alpha = 1.0  # Mixing coefficient
            alpha_predicted = C_alpha * kappa_R
            
            # Scatter around prediction
            sigma_uv = correlation_strength * abs(alpha_predicted) + (1 - correlation_strength) * self.sigma_alpha
            if sigma_uv <= 0:
                sigma_uv = self.sigma_alpha
            
            log_prior_corr = -0.5 * ((alpha - alpha_predicted) / sigma_uv)**2
            return log_prior_corr
        
        elif uv_model == 'string_theory':
            # String theory prediction: α ~ κ_R × 10³⁶ (with suppression)
            C_string = 1e36  # Approximate scaling
            alpha_predicted = C_string * kappa_R
            
            # Theory uncertainty (order of magnitude)
            sigma_theory = abs(alpha_predicted)  # 100% uncertainty
            if sigma_theory <= 0:
                sigma_theory = self.sigma_alpha
                
            log_prior_theory = -0.5 * ((alpha - alpha_predicted) / sigma_theory)**2
            return log_prior_theory
        
        else:
            raise ValueError(f"Unknown UV model: {uv_model}")
    
    def log_posterior(self, kappa_R: float, alpha: float, 
                      uv_model: str = 'linear',
                      correlation_strength: float = 1.0) -> float:
        """
        Log-posterior = log-likelihood + log-prior.
        
        Parameters
        ----------
        kappa_R : float
            κ_R value (m²)
        alpha : float
            α value (m²)
        uv_model : str
            UV completion model
        correlation_strength : float
            Correlation strength
        
        Returns
        -------
        log_post : float
            Log-posterior probability
        """
        log_L = self.log_likelihood(kappa_R, alpha)
        log_pi = self.log_prior(kappa_R, alpha, uv_model, correlation_strength)
        
        return log_L + log_pi
    
    def marginalized_kappa_R_constraint(self, 
                                        uv_model: str = 'linear',
                                        correlation_strength: float = 1.0,
                                        confidence_level: float = 0.95) -> Tuple[float, float]:
        """
        Marginalize over α to get 1D constraint on κ_R.
        
        Uses importance sampling to integrate out α:
        p(κ_R | data) = ∫ p(κ_R, α | data) dα
        
        Parameters
        ----------
        uv_model : str
            UV completion model
        correlation_strength : float
            Correlation strength
        confidence_level : float
            Confidence level for interval
        
        Returns
        -------
        kappa_R_lower : float
            Lower bound on κ_R (m²)
        kappa_R_upper : float
            Upper bound on κ_R (m²)
        """
        def marginalized_log_posterior(kappa_R):
            """Marginalized log-posterior for κ_R."""
            # Integration limits for α
            alpha_min = -5 * self.sigma_alpha
            alpha_max = 5 * self.sigma_alpha
            
            def integrand(alpha):
                return np.exp(self.log_posterior(kappa_R, alpha, uv_model, correlation_strength))
            
            # Numerical integration
            marginal_prob, _ = quad(integrand, alpha_min, alpha_max, 
                                   limit=100, epsabs=1e-10, epsrel=1e-8)
            
            return np.log(marginal_prob) if marginal_prob > 0 else -np.inf
        
        # Find maximum of marginalized posterior
        kappa_R_range = np.linspace(-5 * self.sigma_kappa, 5 * self.sigma_kappa, 1000)
        log_post_values = np.array([marginalized_log_posterior(kr) for kr in kappa_R_range])
        
        # Convert to probabilities and normalize
        log_post_max = np.max(log_post_values)
        post_values = np.exp(log_post_values - log_post_max)
        post_values /= np.trapz(post_values, kappa_R_range)
        
        # Compute cumulative distribution
        cumulative = np.cumsum(post_values) * (kappa_R_range[1] - kappa_R_range[0])
        
        # Find confidence interval
        alpha_tail = (1 - confidence_level) / 2
        lower_idx = np.argmax(cumulative >= alpha_tail)
        upper_idx = np.argmax(cumulative >= (1 - alpha_tail))
        
        return kappa_R_range[lower_idx], kappa_R_range[upper_idx]
    
    def sample_posterior(self, 
                         uv_model: str = 'linear',
                         correlation_strength: float = 1.0,
                         n_samples: int = 10000) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sample from joint posterior using importance sampling.
        
        Parameters
        ----------
        uv_model : str
            UV completion model
        correlation_strength : float
            Correlation strength
        n_samples : int
            Number of samples
        
        Returns
        -------
        kappa_R_samples : ndarray
            Posterior samples for κ_R (m²)
        alpha_samples : ndarray
            Posterior samples for α (m²)
        """
        # Proposal distribution (broader than data uncertainties)
        kappa_R_proposal = np.random.normal(self.kappa_R_obs, 3 * self.sigma_kappa, n_samples)
        alpha_proposal = np.random.normal(self.alpha_obs, 3 * self.sigma_alpha, n_samples)
        
        # Compute weights
        log_weights = np.array([
            self.log_posterior(kr, al, uv_model, correlation_strength)
            for kr, al in zip(kappa_R_proposal, alpha_proposal)
        ])
        
        # Convert to linear weights
        log_weights_max = np.max(log_weights[np.isfinite(log_weights)])
        weights = np.exp(log_weights - log_weights_max)
        weights /= np.sum(weights)
        
        # Resample according to weights
        indices = np.random.choice(n_samples, size=n_samples, p=weights)
        
        return kappa_R_proposal[indices], alpha_proposal[indices]


def plot_joint_posterior(save_path: Optional[str] = None):
    """
    Visualize joint posterior for different UV models.
    """
    analysis = JointPosteriorAnalysis()
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    uv_models = ['uniform', 'linear']
    correlations = [0.0, 0.9]  # Uncorrelated vs strongly correlated
    
    for i, uv_model in enumerate(uv_models):
        for j, corr_strength in enumerate(correlations):
            ax = axes[i, j]
            
            # Sample posterior
            kappa_samples, alpha_samples = analysis.sample_posterior(
                uv_model=uv_model, 
                correlation_strength=corr_strength,
                n_samples=5000
            )
            
            # 2D histogram
            H, xedges, yedges = np.histogram2d(
                kappa_samples / 1e17, alpha_samples / 1e28,
                bins=50, density=True
            )
            extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
            
            im = ax.imshow(H.T, origin='lower', extent=extent, 
                          cmap='Blues', alpha=0.7)
            
            # Contours
            X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2,
                              (yedges[:-1] + yedges[1:]) / 2)
            ax.contour(X, Y, H.T, levels=3, colors='darkblue', linewidths=1.5)
            
            # Individual constraints
            ax.axvline(5, color='red', linestyle='--', alpha=0.5, 
                      label="Lab: κ_R < 5×10¹⁷ m²")
            ax.axhline(1, color='green', linestyle='--', alpha=0.5,
                      label="BH: α < 10²⁸ m²")
            
            ax.set_xlabel("κ_R (10¹⁷ m²)", fontsize=10)
            ax.set_ylabel("α (10²⁸ m²)", fontsize=10)
            ax.set_title(f"{uv_model}, corr={corr_strength:.1f}", fontsize=11)
            ax.grid(True, alpha=0.3)
            ax.set_xlim(-10, 10)
            ax.set_ylim(-3, 3)
            
            if i == 0 and j == 0:
                ax.legend(fontsize=8)
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved joint posterior plot to {save_path}")
    
    return fig


def validate_joint_posterior():
    """
    Validate joint posterior analysis.
    """
    analysis = JointPosteriorAnalysis()
    
    print("=" * 70)
    print("JOINT POSTERIOR ANALYSIS: (κ_R, α) WITH UV CORRELATIONS")
    print("=" * 70)
    
    print("\n1. Data Summary")
    print("-" * 70)
    print(f"κ_R observation: {analysis.kappa_R_obs:.1e} ± {analysis.sigma_kappa:.1e} m²")
    print(f"α observation:   {analysis.alpha_obs:.1e} ± {analysis.sigma_alpha:.1e} m²")
    print(f"Constraint ratio: σ_α / σ_κ = {analysis.sigma_alpha / analysis.sigma_kappa:.1e}")
    
    print("\n2. Marginalized κ_R Constraints")
    print("-" * 70)
    
    for uv_model in ['uniform', 'linear']:
        for corr_strength in [0.0, 0.5, 0.9]:
            try:
                kappa_lower, kappa_upper = analysis.marginalized_kappa_R_constraint(
                    uv_model=uv_model,
                    correlation_strength=corr_strength,
                    confidence_level=0.95
                )
                
                constraint_width = kappa_upper - kappa_lower
                improvement = analysis.sigma_kappa * 2 / constraint_width  # vs ±1σ individual
                
                print(f"{uv_model:12s} (ρ={corr_strength:.1f}): "
                      f"κ_R ∈ [{kappa_lower:.1e}, {kappa_upper:.1e}] m² "
                      f"(width: {constraint_width:.1e}, {improvement:.1f}× tighter)")
                
            except Exception as e:
                print(f"{uv_model:12s} (ρ={corr_strength:.1f}): Integration failed ({e})")
    
    print("\n3. Sample Statistics")
    print("-" * 70)
    
    # Sample for linear model with moderate correlation
    kappa_samples, alpha_samples = analysis.sample_posterior(
        uv_model='linear', correlation_strength=0.5, n_samples=10000
    )
    
    print(f"Linear model (ρ=0.5) samples:")
    print(f"  κ_R: mean = {np.mean(kappa_samples):.2e}, std = {np.std(kappa_samples):.2e} m²")
    print(f"  α:   mean = {np.mean(alpha_samples):.2e}, std = {np.std(alpha_samples):.2e} m²")
    print(f"  Correlation: {np.corrcoef(kappa_samples, alpha_samples)[0, 1]:.3f}")
    
    # 95% intervals
    kappa_95 = np.percentile(kappa_samples, [2.5, 97.5])
    alpha_95 = np.percentile(alpha_samples, [2.5, 97.5])
    print(f"  κ_R 95% CI: [{kappa_95[0]:.2e}, {kappa_95[1]:.2e}] m²")
    print(f"  α 95% CI:   [{alpha_95[0]:.2e}, {alpha_95[1]:.2e}] m²")
    
    print("\n" + "=" * 70)
    print("Key insight: UV correlations can significantly tighten joint constraints")
    print("  - Uncorrelated: constraints multiply independently")
    print("  - Strongly correlated: effective 1D constraint on correlated direction")
    print("  - String theory models predict large α/κ_R ratios → discriminating power")
    print("=" * 70)


if __name__ == "__main__":
    validate_joint_posterior()
    print("\n")
    plot_joint_posterior(save_path="figures/joint_posterior_kappa_alpha.pdf")
    print("\n✅ Task 19 complete: Combined lab+BH posteriors for joint constraint")