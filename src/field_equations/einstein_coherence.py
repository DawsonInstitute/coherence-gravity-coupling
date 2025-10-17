"""
Modified Einstein equations with coherence field coupling.

From variation of the action:
S = ∫ d⁴x √(-g) [ R/(16πG) - (1/2)(∇Φ)² - V(Φ) - ξRΦ² + L_m ]

We get coupled field equations:

MODIFIED EINSTEIN EQUATION:
G_μν + ξ[2(∇_μ∇_ν - g_μν□)Φ² + 2Φ²G_μν - 4∇_μΦ∇_νΦ + 2g_μν(∇Φ)²] = 8πG T_μν

COHERENCE FIELD EQUATION:
□Φ - dV/dΦ - 2ξRΦ = 0

Key features:
1. Curvature sources coherence (2ξRΦ term)
2. Coherence back-reacts on curvature (ξ... terms in Einstein eq)
3. Creates effective field-dependent G
"""

import numpy as np
from typing import Dict, Tuple, Callable
from dataclasses import dataclass

from .action import CoherenceGravityParams, ActionFunctional
# Import weak-field solver and energy reduction from weak_field module
from .weak_field import WeakFieldSolver, compute_energy_cost_reduction


class ModifiedEinsteinEquations:
    """
    Solve modified Einstein equations with coherence coupling.
    
    The full tensor equation is complex, but we'll focus on special cases:
    1. Weak field limit (linearized)
    2. Spherically symmetric static solutions
    3. Cosmological (FLRW) backgrounds
    """
    
    def __init__(self, params: CoherenceGravityParams,
                 action: ActionFunctional):
        """
        Initialize modified Einstein equation solver.
        
        Args:
            params: Physical parameters
            action: Action functional for the theory
        """
        self.params = params
        self.action = action
        
    def modified_einstein_tensor(self, g_mu_nu: np.ndarray,
                                   Christoffel: np.ndarray,
                                   R_mu_nu: np.ndarray,
                                   R: float,
                                   Phi: float,
                                   grad_Phi: np.ndarray,
                                   Hessian_Phi: np.ndarray) -> np.ndarray:
        """
        Compute modified Einstein tensor including coherence corrections.
        
        G̃_μν = G_μν + ξ[2(∇_μ∇_ν - g_μν□)Φ² + 2Φ²G_μν 
                       - 4∇_μΦ∇_νΦ + 2g_μν(∇Φ)²]
        
        Args:
            g_mu_nu: Metric tensor (4x4)
            Christoffel: Christoffel symbols Γ^λ_μν
            R_mu_nu: Ricci tensor (4x4)
            R: Ricci scalar
            Phi: Coherence field value
            grad_Phi: Gradient ∇_μΦ (4-vector)
            Hessian_Phi: ∇_μ∇_νΦ² (4x4 matrix)
            
        Returns:
            Modified Einstein tensor G̃_μν (4x4)
        """
        # Standard Einstein tensor G_μν = R_μν - (1/2)g_μν R
        G_mu_nu = R_mu_nu - 0.5 * g_mu_nu * R
        
        # Coherence corrections (ξ terms)
        xi = self.params.xi
        
        # d'Alembertian operator on Φ²: □Φ² = g^μν ∇_μ∇_νΦ²
        g_inv = np.linalg.inv(g_mu_nu)
        dalembertian_Phi_sq = np.sum(g_inv * Hessian_Phi)
        
        # Term 1: 2(∇_μ∇_νΦ² - g_μν□Φ²)
        term1 = 2 * (Hessian_Phi - g_mu_nu * dalembertian_Phi_sq)
        
        # Term 2: 2Φ²G_μν
        term2 = 2 * Phi**2 * G_mu_nu
        
        # Term 3: -4∇_μΦ∇_νΦ
        grad_outer = np.outer(grad_Phi, grad_Phi)
        term3 = -4 * grad_outer
        
        # Term 4: 2g_μν(∇Φ)²
        grad_Phi_sq = np.dot(grad_Phi, np.dot(g_inv, grad_Phi))
        term4 = 2 * g_mu_nu * grad_Phi_sq
        
        # Modified Einstein tensor
        G_tilde = G_mu_nu + xi * (term1 + term2 + term3 + term4)
        
        return G_tilde
    
    def coherence_field_equation_rhs(self, Phi: float, R: float,
                                      dV_dPhi: Callable[[float], float]) -> float:
        """
        Right-hand side of coherence field equation.
        
        □Φ = dV/dΦ + 2ξRΦ
        
        This gives the source term for the d'Alembertian operator.
        
        Args:
            Phi: Coherence field value
            R: Ricci scalar
            dV_dPhi: Potential derivative function
            
        Returns:
            Source term for □Φ
        """
        xi = self.params.xi
        return dV_dPhi(Phi) + 2 * xi * R * Phi
    
    def effective_stress_energy(self, T_mu_nu_matter: np.ndarray,
                                 Phi: float,
                                 grad_Phi: np.ndarray,
                                 g_mu_nu: np.ndarray,
                                 R: float) -> np.ndarray:
        """
        Compute effective stress-energy tensor including coherence field.
        
        The coherence field contributes to the total stress-energy:
        T^(Phi)_μν = ∇_μΦ∇_νΦ - g_μν[(1/2)(∇Φ)² + V(Φ)]
        
        Args:
            T_mu_nu_matter: Matter stress-energy tensor
            Phi: Coherence field value
            grad_Phi: Gradient of Phi
            g_mu_nu: Metric tensor
            R: Ricci scalar (for potential energy)
            
        Returns:
            Total effective stress-energy tensor
        """
        # Coherence field stress-energy (scalar field form)
        grad_outer = np.outer(grad_Phi, grad_Phi)
        g_inv = np.linalg.inv(g_mu_nu)
        grad_sq = np.dot(grad_Phi, np.dot(g_inv, grad_Phi))
        
        V = self.action.V(np.array([Phi]))[0]
        
        T_Phi = grad_outer - g_mu_nu * (0.5 * grad_sq + V)
        
        # Total stress-energy
        T_total = T_mu_nu_matter + T_Phi
        
        return T_total
