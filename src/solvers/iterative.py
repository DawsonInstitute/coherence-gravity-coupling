"""Iterative solvers for Poisson and elliptic PDEs."""

import numpy as np
from typing import Tuple, Optional, Callable


def gauss_seidel_1d(u: np.ndarray, source: np.ndarray, dx: float,
                    max_iter: int = 1000, tol: float = 1e-6) -> Tuple[np.ndarray, int]:
    """
    Solve 1D Poisson equation using Gauss-Seidel iteration.
    
    ∇²u = source
    
    With Dirichlet boundary conditions u[0] and u[-1] fixed.
    
    Args:
        u: Initial guess (boundaries should contain BCs)
        source: Right-hand side
        dx: Grid spacing
        max_iter: Maximum iterations
        tol: Convergence tolerance
        
    Returns:
        (solution, num_iterations)
    """
    u = u.copy()
    n = len(u)
    
    for iteration in range(max_iter):
        u_old = u.copy()
        
        # Update interior points
        for i in range(1, n-1):
            u[i] = 0.5 * (u[i+1] + u[i-1] - dx**2 * source[i])
        
        # Check convergence
        residual = np.max(np.abs(u - u_old))
        if residual < tol:
            return u, iteration + 1
    
    print(f"Warning: Gauss-Seidel did not converge after {max_iter} iterations")
    return u, max_iter


def jacobi_1d(u: np.ndarray, source: np.ndarray, dx: float,
              max_iter: int = 1000, tol: float = 1e-6) -> Tuple[np.ndarray, int]:
    """
    Solve 1D Poisson equation using Jacobi iteration.
    
    ∇²u = source
    
    Args:
        u: Initial guess
        source: Right-hand side
        dx: Grid spacing
        max_iter: Maximum iterations
        tol: Convergence tolerance
        
    Returns:
        (solution, num_iterations)
    """
    u = u.copy()
    n = len(u)
    u_new = np.zeros(n)
    
    for iteration in range(max_iter):
        u_new[0] = u[0]  # BC
        u_new[-1] = u[-1]  # BC
        
        # Update interior points (using OLD values for all neighbors)
        for i in range(1, n-1):
            u_new[i] = 0.5 * (u[i+1] + u[i-1] - dx**2 * source[i])
        
        # Check convergence
        residual = np.max(np.abs(u_new - u))
        if residual < tol:
            return u_new, iteration + 1
        
        u = u_new.copy()
    
    print(f"Warning: Jacobi did not converge after {max_iter} iterations")
    return u, max_iter


def sor_1d(u: np.ndarray, source: np.ndarray, dx: float,
           omega: float = 1.5, max_iter: int = 1000, tol: float = 1e-6) -> Tuple[np.ndarray, int]:
    """
    Solve 1D Poisson equation using Successive Over-Relaxation (SOR).
    
    SOR is Gauss-Seidel with relaxation parameter ω for faster convergence.
    
    Args:
        u: Initial guess
        source: Right-hand side
        dx: Grid spacing
        omega: Relaxation parameter (1.0 = Gauss-Seidel, optimal often 1.5-1.8)
        max_iter: Maximum iterations
        tol: Convergence tolerance
        
    Returns:
        (solution, num_iterations)
    """
    u = u.copy()
    n = len(u)
    
    for iteration in range(max_iter):
        u_old = u.copy()
        
        # Update interior points with over-relaxation
        for i in range(1, n-1):
            u_gs = 0.5 * (u[i+1] + u[i-1] - dx**2 * source[i])
            u[i] = omega * u_gs + (1 - omega) * u[i]
        
        # Check convergence
        residual = np.max(np.abs(u - u_old))
        if residual < tol:
            return u, iteration + 1
    
    print(f"Warning: SOR did not converge after {max_iter} iterations")
    return u, max_iter


def multigrid_v_cycle_1d(u: np.ndarray, source: np.ndarray, dx: float,
                         v1: int = 2, v2: int = 2) -> np.ndarray:
    """
    Single V-cycle of multigrid for 1D Poisson equation.
    
    Multigrid achieves optimal O(N) convergence by solving coarse-grid corrections.
    
    Args:
        u: Initial guess
        source: Right-hand side
        dx: Grid spacing
        v1: Pre-smoothing iterations
        v2: Post-smoothing iterations
        
    Returns:
        Improved solution
    """
    n = len(u)
    
    if n <= 3:
        # Coarsest level: direct solve
        u, _ = gauss_seidel_1d(u, source, dx, max_iter=100, tol=1e-10)
        return u
    
    # Pre-smoothing
    u, _ = gauss_seidel_1d(u, source, dx, max_iter=v1, tol=0)
    
    # Compute residual
    residual = np.zeros(n)
    for i in range(1, n-1):
        residual[i] = source[i] - (u[i+1] - 2*u[i] + u[i-1]) / dx**2
    
    # Restrict to coarse grid
    n_coarse = (n + 1) // 2
    residual_coarse = np.zeros(n_coarse)
    for i in range(1, n_coarse-1):
        residual_coarse[i] = 0.25 * (residual[2*i-1] + 2*residual[2*i] + residual[2*i+1])
    
    # Solve on coarse grid
    dx_coarse = 2 * dx
    error_coarse = np.zeros(n_coarse)
    error_coarse = multigrid_v_cycle_1d(error_coarse, residual_coarse, dx_coarse, v1, v2)
    
    # Interpolate correction back to fine grid
    error_fine = np.zeros(n)
    for i in range(1, n_coarse-1):
        error_fine[2*i] = error_coarse[i]
        error_fine[2*i-1] = 0.5 * (error_coarse[i-1] + error_coarse[i])
    if n % 2 == 0:
        error_fine[-2] = 0.5 * (error_coarse[-2] + error_coarse[-1])
    
    # Correct
    u = u + error_fine
    
    # Post-smoothing
    u, _ = gauss_seidel_1d(u, source, dx, max_iter=v2, tol=0)
    
    return u


def solve_poisson_1d(source: np.ndarray, dx: float,
                     u_left: float = 0.0, u_right: float = 0.0,
                     method: str = 'gauss_seidel',
                     max_iter: int = 1000, tol: float = 1e-6) -> Tuple[np.ndarray, int]:
    """
    Solve 1D Poisson equation with Dirichlet boundary conditions.
    
    ∇²u = source
    u[0] = u_left, u[-1] = u_right
    
    Args:
        source: Right-hand side
        dx: Grid spacing
        u_left, u_right: Boundary values
        method: 'jacobi', 'gauss_seidel', 'sor', or 'multigrid'
        max_iter: Maximum iterations
        tol: Convergence tolerance
        
    Returns:
        (solution, num_iterations)
    """
    n = len(source)
    u = np.zeros(n)
    u[0] = u_left
    u[-1] = u_right
    
    if method == 'jacobi':
        return jacobi_1d(u, source, dx, max_iter, tol)
    elif method == 'gauss_seidel':
        return gauss_seidel_1d(u, source, dx, max_iter, tol)
    elif method == 'sor':
        return sor_1d(u, source, dx, omega=1.7, max_iter=max_iter, tol=tol)
    elif method == 'multigrid':
        for _ in range(max_iter):
            u_old = u.copy()
            u = multigrid_v_cycle_1d(u, source, dx)
            if np.max(np.abs(u - u_old)) < tol:
                break
        return u, _
    else:
        raise ValueError(f"Unknown method: {method}")
