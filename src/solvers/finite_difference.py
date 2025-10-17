"""Finite difference schemes for PDE solvers."""

import numpy as np
from typing import Tuple


def laplacian_1d(u: np.ndarray, dx: float) -> np.ndarray:
    """
    Compute 1D Laplacian using 2nd-order central differences.
    
    ∇²u ≈ (u[i+1] - 2u[i] + u[i-1]) / dx²
    
    Args:
        u: 1D array of values
        dx: Grid spacing
        
    Returns:
        Laplacian at each interior point
    """
    n = len(u)
    laplacian = np.zeros(n)
    
    # Interior points
    laplacian[1:-1] = (u[2:] - 2*u[1:-1] + u[:-2]) / dx**2
    
    # Boundary points (use one-sided differences or set to zero)
    # For Dirichlet BCs, boundaries are fixed and Laplacian isn't needed
    
    return laplacian


def laplacian_3d(u: np.ndarray, dx: float, dy: float, dz: float) -> np.ndarray:
    """
    Compute 3D Laplacian using 2nd-order central differences.
    
    ∇²u = ∂²u/∂x² + ∂²u/∂y² + ∂²u/∂z²
    
    Args:
        u: 3D array of values (nx, ny, nz)
        dx, dy, dz: Grid spacings
        
    Returns:
        Laplacian at each interior point
    """
    nx, ny, nz = u.shape
    laplacian = np.zeros_like(u)
    
    # Interior points only
    laplacian[1:-1, 1:-1, 1:-1] = (
        (u[2:, 1:-1, 1:-1] - 2*u[1:-1, 1:-1, 1:-1] + u[:-2, 1:-1, 1:-1]) / dx**2 +
        (u[1:-1, 2:, 1:-1] - 2*u[1:-1, 1:-1, 1:-1] + u[1:-1, :-2, 1:-1]) / dy**2 +
        (u[1:-1, 1:-1, 2:] - 2*u[1:-1, 1:-1, 1:-1] + u[1:-1, 1:-1, :-2]) / dz**2
    )
    
    return laplacian


def gradient_1d(u: np.ndarray, dx: float) -> np.ndarray:
    """
    Compute 1D gradient using central differences.
    
    Args:
        u: 1D array
        dx: Grid spacing
        
    Returns:
        Gradient at each interior point
    """
    n = len(u)
    grad = np.zeros(n)
    
    # Central differences at interior points
    grad[1:-1] = (u[2:] - u[:-2]) / (2*dx)
    
    # Forward/backward at boundaries
    grad[0] = (u[1] - u[0]) / dx
    grad[-1] = (u[-1] - u[-2]) / dx
    
    return grad


def laplacian_spherical_1d(u: np.ndarray, r: np.ndarray) -> np.ndarray:
    """
    Compute Laplacian in spherical coordinates (radial part only).
    
    ∇²u = (1/r²) d/dr(r² du/dr)
    
    For spherically symmetric problems.
    
    Args:
        u: Radial function values
        r: Radial grid points
        
    Returns:
        Laplacian in spherical coordinates
    """
    n = len(r)
    laplacian = np.zeros(n)
    
    for i in range(1, n-1):
        dr_plus = r[i+1] - r[i]
        dr_minus = r[i] - r[i-1]
        dr_avg = 0.5 * (dr_plus + dr_minus)
        
        # First derivative: du/dr
        du_dr = (u[i+1] - u[i-1]) / (dr_plus + dr_minus)
        
        # Second derivative with r² factor
        # d/dr(r² du/dr) ≈ [r²_{i+1/2}(du/dr)_{i+1/2} - r²_{i-1/2}(du/dr)_{i-1/2}] / dr
        
        r_plus = 0.5 * (r[i+1] + r[i])
        r_minus = 0.5 * (r[i] + r[i-1])
        
        du_dr_plus = (u[i+1] - u[i]) / dr_plus
        du_dr_minus = (u[i] - u[i-1]) / dr_minus
        
        d_r2_du_dr = (r_plus**2 * du_dr_plus - r_minus**2 * du_dr_minus) / dr_avg
        
        laplacian[i] = d_r2_du_dr / r[i]**2
    
    return laplacian
