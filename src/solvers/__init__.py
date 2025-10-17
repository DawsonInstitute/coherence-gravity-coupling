"""Numerical solvers for coherence-modulated gravity."""

from .finite_difference import laplacian_1d, laplacian_3d, laplacian_spherical_1d
from .iterative import solve_poisson_1d, gauss_seidel_1d, jacobi_1d, sor_1d
from .static_spherical import StaticSphericalSolver

__all__ = [
    'laplacian_1d',
    'laplacian_3d',
    'laplacian_spherical_1d',
    'solve_poisson_1d',
    'gauss_seidel_1d',
    'jacobi_1d',
    'sor_1d',
    'StaticSphericalSolver',
]
