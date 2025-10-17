"""
Action principle for coherence-modulated gravity.

The modified Einstein-Hilbert action with non-minimally coupled coherence field:

S = ∫ d⁴x √(-g) [ R/(16πG) - (1/2)(∇Φ)² - V(Φ) - ξRΦ² + L_m ]

Key terms:
- R/(16πG): Einstein-Hilbert gravitational action
- -(1/2)(∇Φ)²: Coherence field kinetic energy
- -V(Φ): Self-interaction potential
- -ξRΦ²: NON-MINIMAL COUPLING (this is the key innovation!)
- L_m: Matter Lagrangian

The ξRΦ² term creates an effective field-dependent gravitational coupling:
G_eff(Φ) ≈ G(1 - 2ξΦ²) for weak coupling

This allows macroscopic quantum coherence to reduce the energy cost of curvature.
"""

import numpy as np
from dataclasses import dataclass
from typing import Callable, Dict, Tuple
import sympy as sp


@dataclass
class CoherenceGravityParams:
    """Parameters for coherence-gravity coupling theory."""
    
    G: float = 6.674e-11        # Gravitational constant (m³/kg/s²)
    c: float = 2.998e8          # Speed of light (m/s)
    hbar: float = 1.055e-34     # Reduced Planck constant (J·s)
    
    xi: float = 1.0             # Non-minimal coupling strength (dimensionless)
    alpha: float = 1.0          # Exponential suppression factor (for G_eff ansatz)
    
    # Potential parameters (V(Φ) = m²Φ²/2 + λΦ⁴/4)
    m_phi: float = 1.0          # Coherence field mass (kg)
    lambda_phi: float = 0.1     # Self-coupling strength
    
    def __post_init__(self):
        """Validate parameters."""
        if self.G <= 0:
            raise ValueError(f"G must be positive, got {self.G}")
        if self.c <= 0:
            raise ValueError(f"c must be positive, got {self.c}")


class CoherenceFieldPotential:
    """
    Catalog of coherence field potentials V(Φ).
    
    The potential determines the vacuum structure and self-interactions
    of the coherence field.
    """
    
    @staticmethod
    def quadratic(Phi: np.ndarray, m: float) -> np.ndarray:
        """
        Quadratic potential: V(Φ) = (1/2)m²Φ²
        
        This is the simplest case, giving a massive scalar field.
        
        Args:
            Phi: Coherence field value
            m: Mass parameter
            
        Returns:
            Potential energy density
        """
        return 0.5 * m**2 * Phi**2
    
    @staticmethod
    def quadratic_derivative(Phi: np.ndarray, m: float) -> np.ndarray:
        """Derivative dV/dΦ for quadratic potential."""
        return m**2 * Phi
    
    @staticmethod
    def quartic(Phi: np.ndarray, m: float, lambda_: float) -> np.ndarray:
        """
        Quartic potential: V(Φ) = (1/2)m²Φ² + (1/4)λΦ⁴
        
        The λ term provides self-interaction and can lead to spontaneous
        symmetry breaking if m² < 0.
        
        Args:
            Phi: Coherence field value
            m: Mass parameter
            lambda_: Self-coupling strength
            
        Returns:
            Potential energy density
        """
        return 0.5 * m**2 * Phi**2 + 0.25 * lambda_ * Phi**4
    
    @staticmethod
    def quartic_derivative(Phi: np.ndarray, m: float, lambda_: float) -> np.ndarray:
        """Derivative dV/dΦ for quartic potential."""
        return m**2 * Phi + lambda_ * Phi**3
    
    @staticmethod
    def double_well(Phi: np.ndarray, v: float, lambda_: float) -> np.ndarray:
        """
        Double-well potential: V(Φ) = (λ/4)(Φ² - v²)²
        
        This has degenerate vacua at Φ = ±v, enabling topological defects
        and domain walls.
        
        Args:
            Phi: Coherence field value
            v: Vacuum expectation value
            lambda_: Coupling strength
            
        Returns:
            Potential energy density
        """
        return 0.25 * lambda_ * (Phi**2 - v**2)**2
    
    @staticmethod
    def double_well_derivative(Phi: np.ndarray, v: float, lambda_: float) -> np.ndarray:
        """Derivative dV/dΦ for double-well potential."""
        return lambda_ * Phi * (Phi**2 - v**2)


class ActionFunctional:
    """
    Compute action functional and its variations for coherence-gravity system.
    
    The total action is:
    S = S_gravity + S_coherence + S_coupling + S_matter
    
    where:
    S_gravity = ∫ d⁴x √(-g) R/(16πG)
    S_coherence = ∫ d⁴x √(-g) [-(1/2)(∇Φ)² - V(Φ)]
    S_coupling = ∫ d⁴x √(-g) (-ξRΦ²)
    S_matter = ∫ d⁴x √(-g) L_m
    """
    
    def __init__(self, params: CoherenceGravityParams,
                 potential: Callable[[np.ndarray], np.ndarray],
                 potential_deriv: Callable[[np.ndarray], np.ndarray]):
        """
        Initialize action functional.
        
        Args:
            params: Physical parameters
            potential: V(Φ) function
            potential_deriv: dV/dΦ function
        """
        self.params = params
        self.V = potential
        self.dV_dPhi = potential_deriv
        
    def effective_gravitational_constant(self, Phi: np.ndarray) -> np.ndarray:
        """
        Compute effective gravitational coupling G_eff(Φ).
        
        In the non-minimal coupling formalism, the effective Newton constant is:
        G_eff(Φ) = G / (1 + 8πGξΦ²)  (exact)
        
        For weak coupling (8πGξΦ² << 1):
        G_eff(Φ) ≈ G(1 - 8πGξΦ²)
        
        Alternative ansatz (exponential suppression):
        G_eff(Φ) = G exp(-αΦ²)
        
        Args:
            Phi: Coherence field value
            
        Returns:
            Effective gravitational constant
        """
        # Weak-coupling approximation
        xi = self.params.xi
        G = self.params.G
        
        # Exact formula from conformal transformation
        denominator = 1.0 + 8 * np.pi * G * xi * Phi**2
        G_eff = G / denominator
        
        return G_eff
    
    def effective_coupling_ratio(self, Phi: np.ndarray) -> np.ndarray:
        """
        Compute G_eff / G ratio.
        
        This is the key quantity: how much is curvature "cheaper"?
        
        Args:
            Phi: Coherence field value
            
        Returns:
            G_eff(Φ) / G
        """
        G_eff = self.effective_gravitational_constant(Phi)
        return G_eff / self.params.G
    
    def einstein_hilbert_action(self, R: np.ndarray, sqrt_g: np.ndarray, 
                                 d4x: float) -> float:
        """
        Compute Einstein-Hilbert action: ∫ d⁴x √(-g) R/(16πG).
        
        Args:
            R: Ricci scalar at each point
            sqrt_g: √(-g) at each point
            d4x: Volume element
            
        Returns:
            Action value
        """
        integrand = sqrt_g * R / (16 * np.pi * self.params.G)
        return np.sum(integrand) * d4x
    
    def coherence_kinetic_action(self, grad_Phi: np.ndarray, sqrt_g: np.ndarray,
                                  d4x: float) -> float:
        """
        Compute coherence kinetic action: -∫ d⁴x √(-g) (1/2)(∇Φ)².
        
        Args:
            grad_Phi: Gradient of coherence field (magnitude squared)
            sqrt_g: √(-g) at each point
            d4x: Volume element
            
        Returns:
            Action value
        """
        integrand = -0.5 * sqrt_g * grad_Phi**2
        return np.sum(integrand) * d4x
    
    def coherence_potential_action(self, Phi: np.ndarray, sqrt_g: np.ndarray,
                                    d4x: float) -> float:
        """
        Compute potential action: -∫ d⁴x √(-g) V(Φ).
        
        Args:
            Phi: Coherence field
            sqrt_g: √(-g) at each point
            d4x: Volume element
            
        Returns:
            Action value
        """
        V_vals = self.V(Phi)
        integrand = -sqrt_g * V_vals
        return np.sum(integrand) * d4x
    
    def non_minimal_coupling_action(self, R: np.ndarray, Phi: np.ndarray,
                                     sqrt_g: np.ndarray, d4x: float) -> float:
        """
        Compute non-minimal coupling action: -∫ d⁴x √(-g) ξRΦ².
        
        THIS IS THE KEY TERM that makes coherence modify gravitational coupling!
        
        Args:
            R: Ricci scalar
            Phi: Coherence field
            sqrt_g: √(-g) at each point
            d4x: Volume element
            
        Returns:
            Action value
        """
        integrand = -self.params.xi * sqrt_g * R * Phi**2
        return np.sum(integrand) * d4x
    
    def total_action(self, R: np.ndarray, Phi: np.ndarray, grad_Phi: np.ndarray,
                     sqrt_g: np.ndarray, d4x: float) -> Dict[str, float]:
        """
        Compute total action and breakdown by component.
        
        Args:
            R: Ricci scalar
            Phi: Coherence field
            grad_Phi: ∇Φ magnitude squared
            sqrt_g: √(-g)
            d4x: Volume element
            
        Returns:
            Dictionary with total action and components
        """
        S_EH = self.einstein_hilbert_action(R, sqrt_g, d4x)
        S_kin = self.coherence_kinetic_action(grad_Phi, sqrt_g, d4x)
        S_pot = self.coherence_potential_action(Phi, sqrt_g, d4x)
        S_coupling = self.non_minimal_coupling_action(R, Phi, sqrt_g, d4x)
        
        return {
            'total': S_EH + S_kin + S_pot + S_coupling,
            'einstein_hilbert': S_EH,
            'coherence_kinetic': S_kin,
            'coherence_potential': S_pot,
            'non_minimal_coupling': S_coupling,
            'gravity_effective': S_EH + S_coupling  # Combined gravitational sector
        }


def symbolic_field_equations():
    """
    Derive field equations symbolically using SymPy.
    
    This provides exact symbolic expressions for the modified Einstein
    and coherence field equations.
    
    Returns:
        Dictionary with symbolic equations
    """
    # Define symbolic variables
    t, x, y, z = sp.symbols('t x y z', real=True)
    Phi = sp.Function('Phi')(t, x, y, z)
    R = sp.Symbol('R')  # Ricci scalar (symbolic)
    
    G_const, c, xi, m = sp.symbols('G c xi m', positive=True, real=True)
    lambda_ = sp.symbols('lambda', real=True)
    
    # Define Lagrangian density
    # Use symbolic placeholder for sqrt(-g) to avoid matrix operations
    sqrt_g = sp.Symbol('sqrt_g', positive=True)
    
    # Components (symbolic forms)
    L_gravity = R / (16 * sp.pi * G_const)
    L_coherence_kin = -sp.Rational(1, 2) * sp.Symbol('grad_Phi_sq')  # (∇Φ)²
    L_coherence_pot = -sp.Rational(1, 2) * m**2 * Phi**2 - lambda_ * Phi**4 / 4
    L_coupling = -xi * R * Phi**2
    
    # Total Lagrangian
    L_total = L_gravity + L_coherence_kin + L_coherence_pot + L_coupling
    
    # Field equation for Phi (from δS/δΦ = 0)
    # □Φ - dV/dΦ - 2ξRΦ = 0
    dV_dPhi = m**2 * Phi + lambda_ * Phi**3
    Phi_equation = sp.Symbol('square_Phi') - dV_dPhi - 2*xi*R*Phi
    
    return {
        'lagrangian': L_total,
        'coherence_equation': Phi_equation,
        'effective_G': G_const / (1 + 8*sp.pi*G_const*xi*Phi**2),
        'coupling_ratio': 1 / (1 + 8*sp.pi*G_const*xi*Phi**2)
    }
