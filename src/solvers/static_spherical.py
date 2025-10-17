"""
Spherically symmetric static solutions for coherence-modulated gravity.

Solves:
∇²Φ_grav = 4πG_eff(Φ₀)ρ(r)

in spherical coordinates for point mass and extended sources.
"""

import numpy as np
from typing import Tuple, Dict, Callable
import sys
from pathlib import Path

# Add parent to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from field_equations.action import CoherenceGravityParams
from .finite_difference import laplacian_spherical_1d
from .iterative import solve_poisson_1d


class StaticSphericalSolver:
    """
    Solve for static spherically symmetric gravitational potential
    with coherence-modulated coupling.
    """
    
    def __init__(self, params: CoherenceGravityParams):
        """
        Initialize solver.
        
        Args:
            params: Physical parameters including G, xi
        """
        self.params = params
    
    def effective_G(self, Phi0: float) -> float:
        """
        Compute effective gravitational constant.
        
        Args:
            Phi0: Background coherence field [m⁻¹]
            
        Returns:
            G_eff [m³/(kg·s²)]
        """
        G = self.params.G
        xi = self.params.xi
        return G / (1.0 + 8*np.pi*G*xi*Phi0**2)
    
    def point_mass_analytic(self, r: np.ndarray, M: float, Phi0: float) -> Dict[str, np.ndarray]:
        """
        Analytic solution for point mass.
        
        Φ_grav(r) = -G_eff M / r
        
        Args:
            r: Radial coordinates [m]
            M: Point mass [kg]
            Phi0: Background coherence [m⁻¹]
            
        Returns:
            Dictionary with potential, force, and diagnostics
        """
        G_eff = self.effective_G(Phi0)
        
        Phi_grav = -G_eff * M / r
        F = G_eff * M / r**2  # Radial force magnitude
        
        return {
            'r': r,
            'potential': Phi_grav,
            'force': F,
            'G_eff': G_eff,
            'M': M,
            'Phi0': Phi0,
            'suppression': G_eff / self.params.G
        }
    
    def point_mass_comparison(self, r: np.ndarray, M: float, Phi0: float) -> Dict[str, np.ndarray]:
        """
        Compare coherence-modulated vs standard GR for point mass.
        
        Args:
            r: Radial coordinates
            M: Point mass
            Phi0: Coherence background
            
        Returns:
            Comparison data
        """
        result = self.point_mass_analytic(r, M, Phi0)
        
        # Standard GR
        Phi_standard = -self.params.G * M / r
        F_standard = self.params.G * M / r**2
        
        return {
            **result,
            'potential_standard': Phi_standard,
            'force_standard': F_standard,
            'potential_ratio': result['potential'] / Phi_standard,
            'force_ratio': result['force'] / F_standard
        }
    
    def uniform_sphere_analytic(self, r: np.ndarray, M: float, R: float, Phi0: float) -> Dict:
        """
        Analytic solution for uniform density sphere.
        
        ρ = 3M/(4πR³) for r < R, 0 otherwise
        
        Interior (r < R): Φ = -G_eff M (3R² - r²) / (2R³)
        Exterior (r > R): Φ = -G_eff M / r
        
        Args:
            r: Radial coordinates
            M: Total mass
            R: Sphere radius
            Phi0: Coherence background
            
        Returns:
            Potential and force
        """
        G_eff = self.effective_G(Phi0)
        
        Phi = np.zeros_like(r)
        F = np.zeros_like(r)
        
        # Interior
        interior = r < R
        Phi[interior] = -G_eff * M * (3*R**2 - r[interior]**2) / (2*R**3)
        F[interior] = G_eff * M * r[interior] / R**3
        
        # Exterior
        exterior = r >= R
        Phi[exterior] = -G_eff * M / r[exterior]
        F[exterior] = G_eff * M / r[exterior]**2
        
        return {
            'r': r,
            'potential': Phi,
            'force': F,
            'G_eff': G_eff,
            'M': M,
            'R': R,
            'Phi0': Phi0
        }
    
    def solve_with_density(self, r: np.ndarray, rho: Callable[[np.ndarray], np.ndarray],
                           Phi0: float, method: str = 'gauss_seidel') -> Dict:
        """
        Numerically solve Poisson equation for arbitrary density profile.
        
        (1/r²) d/dr(r² dΦ/dr) = 4πG_eff ρ(r)
        
        Args:
            r: Radial grid
            rho: Density function ρ(r) [kg/m³]
            Phi0: Coherence background
            method: Iterative solver method
            
        Returns:
            Potential and diagnostics
        """
        G_eff = self.effective_G(Phi0)
        n = len(r)
        dr = r[1] - r[0]  # Assume uniform grid
        
        # Source term: 4πG_eff ρ(r)
        source = 4 * np.pi * G_eff * rho(r)
        
        # Initial guess (point mass approximation)
        M_total = np.trapz(4*np.pi*r**2 * rho(r), r)
        Phi = -G_eff * M_total / r
        
        # Solve using finite differences
        # Convert to Cartesian-like form for 1D solver
        # d²Φ/dr² + (2/r) dΦ/dr = 4πG_eff ρ
        
        # This requires custom treatment of the (2/r) term
        # For simplicity, use direct iteration on the spherical Laplacian
        
        max_iter = 1000
        tol = 1e-6
        
        for iteration in range(max_iter):
            Phi_old = Phi.copy()
            
            # Compute Laplacian
            lap = laplacian_spherical_1d(Phi, r)
            
            # Residual
            residual = lap - source
            
            # Update (damped for stability)
            omega = 0.5
            Phi = Phi - omega * dr**2 * residual
            
            # Apply BCs
            Phi[0] = Phi[1]  # Neumann at origin (regularity)
            Phi[-1] = -G_eff * M_total / r[-1]  # Asymptotic
            
            # Check convergence
            if np.max(np.abs(Phi - Phi_old)) < tol:
                break
        
        # Compute force: F = -dΦ/dr
        F = np.zeros_like(r)
        F[1:-1] = -(Phi[2:] - Phi[:-2]) / (2*dr)
        F[0] = 0  # Force vanishes at origin by symmetry
        F[-1] = G_eff * M_total / r[-1]**2
        
        return {
            'r': r,
            'potential': Phi,
            'force': F,
            'G_eff': G_eff,
            'density': rho(r),
            'iterations': iteration + 1,
            'Phi0': Phi0
        }
    
    def coherent_shell_example(self, r: np.ndarray, M_core: float, R_core: float,
                                Phi0_interior: float, Phi0_exterior: float) -> Dict:
        """
        Example: mass with coherent shell creating different G_eff in different regions.
        
        Args:
            r: Radial grid
            M_core: Core mass
            R_core: Core radius
            Phi0_interior: Coherence inside shell
            Phi0_exterior: Coherence outside shell
            
        Returns:
            Potential showing effect of coherence transition
        """
        G_int = self.effective_G(Phi0_interior)
        G_ext = self.effective_G(Phi0_exterior)
        
        Phi = np.zeros_like(r)
        
        # Match interior and exterior at R_core
        # Interior: Φ_int = -G_int M_core / R_core + C
        # Exterior: Φ_ext = -G_ext M_core / r
        
        # Continuity: Φ_int(R_core) = Φ_ext(R_core)
        # -G_int M_core / R_core + C = -G_ext M_core / R_core
        # C = (G_int - G_ext) M_core / R_core
        
        C = (G_int - G_ext) * M_core / R_core
        
        interior = r < R_core
        Phi[interior] = -G_int * M_core / R_core + C
        
        exterior = r >= R_core
        Phi[exterior] = -G_ext * M_core / r[exterior]
        
        return {
            'r': r,
            'potential': Phi,
            'G_int': G_int,
            'G_ext': G_ext,
            'Phi0_interior': Phi0_interior,
            'Phi0_exterior': Phi0_exterior,
            'R_core': R_core,
            'suppression_ratio': G_int / G_ext
        }


if __name__ == "__main__":
    # Test solver
    params = CoherenceGravityParams(xi=1.0)
    solver = StaticSphericalSolver(params)
    
    # Point mass
    r = np.linspace(1e-2, 10.0, 100)  # m
    M = 1e3  # kg
    Phi0 = 1e15  # m⁻¹
    
    result = solver.point_mass_comparison(r, M, Phi0)
    
    print("Point Mass Test:")
    print(f"G_eff/G = {result['suppression']:.6e}")
    print(f"Potential suppression: {result['potential_ratio'][50]:.6e}")
    print()
    
    # Coherent shell
    result_shell = solver.coherent_shell_example(r, M, 1.0, Phi0, 0.0)
    print("Coherent Shell Test:")
    print(f"G_interior/G_exterior = {result_shell['suppression_ratio']:.6e}")
