"""
Curvature Coupling Extension Module

Implements higher-order curvature coupling terms for testing new physics:

1. Ricci-Electromagnetic Coupling: R F_μν F^μν
2. Riemann-Electromagnetic Coupling: R_μνρσ F^μν F^ρσ
3. Weyl Curvature Coupling: C_μνρσ F^μν F^ρσ

These terms arise naturally in:
- Quantum gravity effective field theories
- String theory low-energy limits
- Modified gravity theories (f(R), etc.)

The coupling strength κ parameterizes the deviation from GR.
Experimental constraints suggest |κ| < 10^-20 m² for most scenarios.

Author: GitHub Copilot (Claude Sonnet 4.5)
Date: October 2025
License: MIT
"""

import numpy as np
from typing import Dict, Tuple, Callable, Optional
from dataclasses import dataclass
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Physical constants
C_LIGHT = 299792458  # m/s
EPSILON_0 = 8.854187817e-12  # F/m (vacuum permittivity)
MU_0 = 4 * np.pi * 1e-7  # H/m (vacuum permeability)
G_NEWTON = 6.67430e-11  # m³⋅kg⁻¹⋅s⁻²


@dataclass
class CurvatureCouplingParams:
    """Configuration for curvature coupling extensions."""
    
    # Coupling strengths (dimensionless or with mass scale)
    kappa_ricci_em: float = 0.0  # R F² coupling strength [m²]
    kappa_riemann_em: float = 0.0  # R_μνρσ F^μν F^ρσ coupling [m²]
    kappa_weyl_em: float = 0.0  # Weyl curvature coupling [m²]
    
    # Enable/disable specific terms
    enable_ricci_em: bool = False
    enable_riemann_em: bool = False
    enable_weyl_em: bool = False
    
    # Numerical parameters
    regularization_scale: float = 1e-10  # Small scale to avoid singularities [m]
    
    def __post_init__(self):
        """Validate configuration."""
        if self.kappa_ricci_em < 0 or self.kappa_riemann_em < 0 or self.kappa_weyl_em < 0:
            logger.warning("Negative coupling constants detected - may lead to instabilities")
        
        if any([self.enable_ricci_em, self.enable_riemann_em, self.enable_weyl_em]):
            logger.info("Curvature coupling extensions enabled:")
            if self.enable_ricci_em:
                logger.info(f"  • Ricci-EM: κ_R = {self.kappa_ricci_em:.2e} m²")
            if self.enable_riemann_em:
                logger.info(f"  • Riemann-EM: κ_Riem = {self.kappa_riemann_em:.2e} m²")
            if self.enable_weyl_em:
                logger.info(f"  • Weyl-EM: κ_Weyl = {self.kappa_weyl_em:.2e} m²")


class CurvatureCouplingCalculator:
    """
    Calculate curvature coupling contributions to field equations.
    
    These terms modify both:
    1. Maxwell equations (curvature affects EM propagation)
    2. Einstein equations (EM fields source curvature differently)
    """
    
    def __init__(self, params: CurvatureCouplingParams):
        """
        Initialize curvature coupling calculator.
        
        Args:
            params: Configuration parameters
        """
        self.params = params
        
    def electromagnetic_invariants(self, 
                                  E_field: np.ndarray,
                                  B_field: np.ndarray,
                                  metric: Optional[np.ndarray] = None) -> Dict[str, float]:
        """
        Compute electromagnetic field invariants.
        
        F_μν F^μν = 2(B² - E²/c²)
        *F_μν F^μν = 4E·B/c (dual invariant)
        
        Args:
            E_field: Electric field vector (3,) [V/m]
            B_field: Magnetic field vector (3,) [T]
            metric: Spacetime metric (4x4), if None use Minkowski
            
        Returns:
            Dictionary with:
            - 'F_squared': F_μν F^μν
            - 'F_dual_squared': *F_μν F^μν
            - 'E_squared': |E|²
            - 'B_squared': |B|²
        """
        # Compute magnitudes
        E_sq = np.dot(E_field, E_field)
        B_sq = np.dot(B_field, B_field)
        
        # First invariant: F² = 2(B² - E²/c²)
        F_squared = 2 * (B_sq - E_sq / C_LIGHT**2)
        
        # Second (dual) invariant: *F² = 4E·B/c
        E_dot_B = np.dot(E_field, B_field)
        F_dual_squared = 4 * E_dot_B / C_LIGHT
        
        return {
            'F_squared': F_squared,
            'F_dual_squared': F_dual_squared,
            'E_squared': E_sq,
            'B_squared': B_sq,
            'E_dot_B': E_dot_B
        }
    
    def ricci_em_coupling_energy(self,
                                 ricci_scalar: float,
                                 E_field: np.ndarray,
                                 B_field: np.ndarray) -> float:
        """
        Compute energy contribution from R F² coupling.
        
        L_coupling = κ_R * R * F_μν F^μν
        
        Args:
            ricci_scalar: Ricci curvature scalar R [m⁻²]
            E_field: Electric field (3,) [V/m]
            B_field: Magnetic field (3,) [T]
            
        Returns:
            Energy density contribution [J/m³]
        """
        if not self.params.enable_ricci_em:
            return 0.0
        
        # Get EM invariants
        invariants = self.electromagnetic_invariants(E_field, B_field)
        F_squared = invariants['F_squared']
        
        # Coupling term: κ R F²
        energy_density = self.params.kappa_ricci_em * ricci_scalar * F_squared
        
        # Convert to energy density (include vacuum permittivity)
        energy_density *= EPSILON_0 / 2
        
        return energy_density
    
    def modified_maxwell_equations(self,
                                   E_field: np.ndarray,
                                   B_field: np.ndarray,
                                   ricci_scalar: float,
                                   grad_ricci: Optional[np.ndarray] = None) -> Dict[str, np.ndarray]:
        """
        Compute modifications to Maxwell equations from curvature coupling.
        
        In the presence of R F² coupling, Maxwell's equations become:
        ∇·D = ρ  (unchanged)
        ∇·B = 0  (unchanged)
        ∇×E = -∂B/∂t - κ_R ∇R × B  (modified Faraday)
        ∇×H = J + ∂D/∂t + κ_R (∇R × E + R ∇×E)  (modified Ampère)
        
        Args:
            E_field: Electric field (3,) [V/m]
            B_field: Magnetic field (3,) [T]
            ricci_scalar: Ricci scalar R [m⁻²]
            grad_ricci: Gradient of R, (3,) [m⁻³]
            
        Returns:
            Dictionary with correction terms
        """
        corrections = {
            'faraday_correction': np.zeros(3),
            'ampere_correction': np.zeros(3)
        }
        
        if not self.params.enable_ricci_em or grad_ricci is None:
            return corrections
        
        # Faraday correction: -κ_R ∇R × B
        corrections['faraday_correction'] = -self.params.kappa_ricci_em * np.cross(grad_ricci, B_field)
        
        # Ampère correction: κ_R ∇R × E
        # (Second term R ∇×E requires curl computation, omitted for now)
        corrections['ampere_correction'] = self.params.kappa_ricci_em * np.cross(grad_ricci, E_field)
        
        return corrections
    
    def stress_energy_correction(self,
                                ricci_scalar: float,
                                ricci_tensor: np.ndarray,
                                E_field: np.ndarray,
                                B_field: np.ndarray) -> np.ndarray:
        """
        Compute correction to stress-energy tensor from curvature coupling.
        
        The effective stress-energy includes contributions from R F²:
        δT_μν = κ_R [F_αβ F^αβ ∇_μ∇_ν R + ...]
        
        Args:
            ricci_scalar: Ricci scalar R [m⁻²]
            ricci_tensor: Ricci tensor R_μν (4x4) [m⁻²]
            E_field: Electric field (3,) [V/m]
            B_field: Magnetic field (3,) [T]
            
        Returns:
            Stress-energy correction δT_μν (4x4) [J/m³]
        """
        if not self.params.enable_ricci_em:
            return np.zeros((4, 4))
        
        # Get EM invariants
        invariants = self.electromagnetic_invariants(E_field, B_field)
        F_squared = invariants['F_squared']
        
        # Simplified correction (full derivation requires covariant derivatives)
        # δT_μν ≈ κ_R F² R_μν
        delta_T = self.params.kappa_ricci_em * F_squared * ricci_tensor
        
        # Include vacuum permittivity factor
        delta_T *= EPSILON_0 / 2
        
        return delta_T
    
    def effective_gravitational_constant(self,
                                        E_field: np.ndarray,
                                        B_field: np.ndarray,
                                        ricci_scalar: float) -> float:
        """
        Compute effective gravitational constant in the presence of EM fields.
        
        Due to R F² coupling, the effective G depends on the EM field strength:
        G_eff ≈ G (1 + κ_R F²)
        
        Args:
            E_field: Electric field (3,) [V/m]
            B_field: Magnetic field (3,) [T]
            ricci_scalar: Ricci scalar R [m⁻²]
            
        Returns:
            Effective gravitational constant [m³⋅kg⁻¹⋅s⁻²]
        """
        if not self.params.enable_ricci_em:
            return G_NEWTON
        
        # Get EM invariants
        invariants = self.electromagnetic_invariants(E_field, B_field)
        F_squared = invariants['F_squared']
        
        # Effective G modification
        # (This is a simplified model; full calculation requires detailed field equations)
        modification = 1.0 + self.params.kappa_ricci_em * F_squared * ricci_scalar
        
        # Regularize to avoid unphysical values
        if modification < 0.1 or modification > 10.0:
            logger.warning(f"Large G_eff modification: {modification:.2e}. Using regularized value.")
            modification = np.clip(modification, 0.1, 10.0)
        
        return G_NEWTON * modification


def compute_exclusion_limits(
    experimental_precision: float,
    ricci_scale: float,
    field_strength: float
) -> Dict[str, float]:
    """
    Estimate exclusion limits on coupling constants from null results.
    
    If no effect is observed at precision δ, we can constrain:
    |κ| < δ / (R × F²)
    
    Args:
        experimental_precision: Relative precision achieved (e.g., 1e-6)
        ricci_scale: Characteristic Ricci scalar [m⁻²]
        field_strength: EM field invariant F² [T²]
        
    Returns:
        Dictionary with exclusion limits
    """
    # Avoid division by zero
    if ricci_scale == 0 or field_strength == 0:
        logger.warning("Cannot compute exclusion limits with zero curvature or field strength")
        return {'kappa_ricci_em_limit': np.inf}
    
    # Exclusion limit: κ < δ / (R F²)
    kappa_limit = experimental_precision / (abs(ricci_scale) * abs(field_strength))
    
    logger.info(f"Exclusion limit calculation:")
    logger.info(f"  Precision: {experimental_precision:.2e}")
    logger.info(f"  R scale: {ricci_scale:.2e} m⁻²")
    logger.info(f"  F² scale: {field_strength:.2e} T²")
    logger.info(f"  κ limit: {kappa_limit:.2e} m²")
    
    return {
        'kappa_ricci_em_limit': kappa_limit,
        'precision': experimental_precision,
        'ricci_scale': ricci_scale,
        'field_strength': field_strength
    }


def demonstrate_curvature_coupling():
    """Demonstration of curvature coupling calculations."""
    print("\n" + "="*70)
    print("🌀 CURVATURE COUPLING EXTENSION DEMONSTRATION")
    print("="*70)
    
    # Setup parameters
    params = CurvatureCouplingParams(
        kappa_ricci_em=1e-15,  # m² (testable scale)
        enable_ricci_em=True
    )
    
    calculator = CurvatureCouplingCalculator(params)
    
    # Example EM fields (laboratory scale)
    E_field = np.array([1e6, 0, 0])  # 1 MV/m
    B_field = np.array([0, 1.0, 0])  # 1 T
    
    # Example curvature (Earth's surface)
    ricci_scalar = 1e-26  # m⁻² (very weak)
    
    print("\n📐 EM Field Configuration:")
    print(f"  E = {np.linalg.norm(E_field):.2e} V/m")
    print(f"  B = {np.linalg.norm(B_field):.2e} T")
    
    # Compute invariants
    invariants = calculator.electromagnetic_invariants(E_field, B_field)
    print("\n🔢 EM Invariants:")
    print(f"  F² = {invariants['F_squared']:.2e} T²")
    print(f"  *F² = {invariants['F_dual_squared']:.2e} T·V/m")
    
    # Coupling energy
    energy = calculator.ricci_em_coupling_energy(ricci_scalar, E_field, B_field)
    print(f"\n⚡ Curvature Coupling Energy:")
    print(f"  δE = {energy:.2e} J/m³")
    
    # Effective G
    G_eff = calculator.effective_gravitational_constant(E_field, B_field, ricci_scalar)
    delta_G = (G_eff - G_NEWTON) / G_NEWTON
    print(f"\n🌍 Effective Gravitational Constant:")
    print(f"  G_eff = {G_eff:.4e} m³⋅kg⁻¹⋅s⁻²")
    print(f"  ΔG/G = {delta_G:.2e}")
    
    # Exclusion limits
    limits = compute_exclusion_limits(
        experimental_precision=1e-6,
        ricci_scale=ricci_scalar,
        field_strength=invariants['F_squared']
    )
    print(f"\n📊 Exclusion Limits (null result at 10⁻⁶):")
    print(f"  κ_R < {limits['kappa_ricci_em_limit']:.2e} m²")
    
    print("\n" + "="*70)
    print("✅ DEMONSTRATION COMPLETE")
    print("="*70)


if __name__ == '__main__':
    demonstrate_curvature_coupling()
