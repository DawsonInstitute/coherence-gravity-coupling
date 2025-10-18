"""
3D Poisson Solver for Spatially-Varying G_eff

Solves the modified gravitational Poisson equation:
    ∇·(G_eff(x)∇φ(x)) = 4πGρ(x)

where G_eff(x) = G / (1 + 8πGξΦ²(x))

Uses:
- 7-point finite difference Laplacian on cubic grid
- SciPy sparse iterative solver (CG or BiCGSTAB)
- Dirichlet boundary conditions φ(∞) = 0

Author: GitHub Copilot (Claude Sonnet 4.5)
License: MIT
"""

import numpy as np
from typing import Callable, Dict, Tuple, Optional
from dataclasses import dataclass
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import cg, bicgstab, LinearOperator
from pathlib import Path
import json

# Optional: PyAMG for algebraic multigrid preconditioning
try:
    import pyamg
    HAS_PYAMG = True
except ImportError:
    HAS_PYAMG = False

# scipy incomplete factorizations
try:
    from scipy.sparse.linalg import spilu
    HAS_SPILU = True
except ImportError:
    HAS_SPILU = False

# Physical constants
G_SI = 6.674e-11  # m³/(kg·s²)


@dataclass
class Grid3D:
    """3D cubic grid specification"""
    nx: int  # Grid points in x
    ny: int  # Grid points in y
    nz: int  # Grid points in z
    Lx: float  # Domain size in x [m]
    Ly: float  # Domain size in y [m]
    Lz: float  # Domain size in z [m]
    
    @property
    def dx(self) -> float:
        return self.Lx / (self.nx - 1)
    
    @property
    def dy(self) -> float:
        return self.Ly / (self.ny - 1)
    
    @property
    def dz(self) -> float:
        return self.Lz / (self.nz - 1)
    
    @property
    def total_points(self) -> int:
        return self.nx * self.ny * self.nz
    
    def index(self, i: int, j: int, k: int) -> int:
        """Convert (i,j,k) to flat index"""
        return i + self.nx * (j + self.ny * k)
    
    def coord(self, i: int, j: int, k: int) -> Tuple[float, float, float]:
        """Get physical coordinates of grid point"""
        x = -self.Lx/2 + i * self.dx
        y = -self.Ly/2 + j * self.dy
        z = -self.Lz/2 + k * self.dz
        return x, y, z


@dataclass
class PoissonSolution:
    """Solution to 3D Poisson equation"""
    grid: Grid3D
    phi: np.ndarray  # Potential field [m²/s²]
    G_eff: np.ndarray  # Effective G field [m³/(kg·s²)]
    rho: np.ndarray  # Mass density [kg/m³]
    Phi_coherence: np.ndarray  # Coherence field [m⁻¹]
    solver_info: Dict
    
    def save(self, filepath: Path):
        """Save solution to disk"""
        np.savez_compressed(
            filepath,
            phi=self.phi,
            G_eff=self.G_eff,
            rho=self.rho,
            Phi_coherence=self.Phi_coherence,
            grid_nx=self.grid.nx,
            grid_ny=self.grid.ny,
            grid_nz=self.grid.nz,
            grid_Lx=self.grid.Lx,
            grid_Ly=self.grid.Ly,
            grid_Lz=self.grid.Lz
        )
        
        # Save solver metadata
        json_path = filepath.with_suffix('.json')
        with open(json_path, 'w') as f:
            json.dump(self.solver_info, f, indent=2)
        
        print(f"   Saved solution: {filepath.name}")


class Poisson3DSolver:
    """
    3D Poisson solver for modified gravitational equation.
    
    Discretizes ∇·(G_eff(x)∇φ) = 4πGρ using 7-point stencil:
    
        (G_e·∂φ/∂x)|_{i+½} - (G_e·∂φ/∂x)|_{i-½}     (G_e·∂φ/∂y)|_{j+½} - (G_e·∂φ/∂y)|_{j-½}
        ───────────────────────────────────────  +  ───────────────────────────────────────
                        Δx²                                         Δy²
        
        (G_e·∂φ/∂z)|_{k+½} - (G_e·∂φ/∂z)|_{k-½}
      + ───────────────────────────────────────  =  4πGρ_{i,j,k}
                        Δz²
    
    where G_e = G_eff is evaluated at cell faces.
    """
    
    def __init__(self, grid: Grid3D, xi: float = 1.0):
        """
        Initialize solver.
        
        Args:
            grid: 3D grid specification
            xi: Non-minimal coupling strength (dimensionless)
        """
        self.grid = grid
        self.xi = xi
        
    def compute_G_eff(self, Phi_field: np.ndarray) -> np.ndarray:
        """
        Compute G_eff(x) from coherence field Φ(x).
        
        G_eff = G / (1 + 8πGξΦ²)
        
        Args:
            Phi_field: Coherence field Φ(x) [m⁻¹], shape (nx, ny, nz)
        
        Returns:
            G_eff field [m³/(kg·s²)], same shape
        """
        denominator = 1.0 + 8.0 * np.pi * G_SI * self.xi * Phi_field**2
        return G_SI / denominator
    
    def build_linear_system(
        self,
        rho: np.ndarray,
        G_eff: np.ndarray
    ) -> Tuple[csr_matrix, np.ndarray]:
        """
        Build sparse matrix A and RHS b for Ax = b.
        
        Uses 7-point stencil with G_eff evaluated at cell faces:
        - Face i+½: G_eff_{i+½,j,k} = 0.5*(G_eff_{i,j,k} + G_eff_{i+1,j,k})
        - Similar for j±½, k±½ faces
        
        Args:
            rho: Mass density [kg/m³], shape (nx, ny, nz)
            G_eff: Effective coupling [m³/(kg·s²)], shape (nx, ny, nz)
        
        Returns:
            A: Sparse coefficient matrix (N×N where N = nx*ny*nz)
            b: RHS vector (N,)
        """
        nx, ny, nz = self.grid.nx, self.grid.ny, self.grid.nz
        dx, dy, dz = self.grid.dx, self.grid.dy, self.grid.dz
        N = self.grid.total_points
        
        A = lil_matrix((N, N))
        b = np.zeros(N)
        
        print(f"   Building {N}×{N} sparse system...")
        
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    idx = self.grid.index(i, j, k)
                    
                    # Boundary conditions: φ = 0 at domain edges
                    if i == 0 or i == nx-1 or j == 0 or j == ny-1 or k == 0 or k == nz-1:
                        A[idx, idx] = 1.0
                        b[idx] = 0.0
                        continue
                    
                    # Interior points: 7-point stencil
                    # G_eff at cell faces (harmonic mean for better numerics)
                    G_xp = 2.0 / (1.0/G_eff[i,j,k] + 1.0/G_eff[i+1,j,k])  # i+½
                    G_xm = 2.0 / (1.0/G_eff[i,j,k] + 1.0/G_eff[i-1,j,k])  # i-½
                    G_yp = 2.0 / (1.0/G_eff[i,j,k] + 1.0/G_eff[i,j+1,k])  # j+½
                    G_ym = 2.0 / (1.0/G_eff[i,j,k] + 1.0/G_eff[i,j-1,k])  # j-½
                    G_zp = 2.0 / (1.0/G_eff[i,j,k] + 1.0/G_eff[i,j,k+1])  # k+½
                    G_zm = 2.0 / (1.0/G_eff[i,j,k] + 1.0/G_eff[i,j,k-1])  # k-½
                    
                    # Stencil coefficients
                    c_center = -(G_xp + G_xm)/dx**2 - (G_yp + G_ym)/dy**2 - (G_zp + G_zm)/dz**2
                    c_xp = G_xp / dx**2
                    c_xm = G_xm / dx**2
                    c_yp = G_yp / dy**2
                    c_ym = G_ym / dy**2
                    c_zp = G_zp / dz**2
                    c_zm = G_zm / dz**2
                    
                    # Fill matrix row
                    A[idx, idx] = c_center
                    A[idx, self.grid.index(i+1, j, k)] = c_xp
                    A[idx, self.grid.index(i-1, j, k)] = c_xm
                    A[idx, self.grid.index(i, j+1, k)] = c_yp
                    A[idx, self.grid.index(i, j-1, k)] = c_ym
                    A[idx, self.grid.index(i, j, k+1)] = c_zp
                    A[idx, self.grid.index(i, j, k-1)] = c_zm
                    
                    # RHS
                    b[idx] = 4.0 * np.pi * G_SI * rho[i,j,k]
        
        print(f"   Converting to CSR format...")
        A_csr = A.tocsr()
        nnz = A_csr.nnz
        sparsity = 100.0 * (1.0 - nnz / N**2)
        print(f"   Matrix sparsity: {sparsity:.2f}% (nnz={nnz:,})")
        
        return A_csr, b
    
    def solve(
        self,
        rho_func: Callable[[float, float, float], float],
        Phi_func: Callable[[float, float, float], float],
        method: str = 'cg',
        tol: float = 1e-8,
        maxiter: Optional[int] = None,
        preconditioner: str = 'none'
    ) -> PoissonSolution:
        """
        Solve modified Poisson equation.
        
        Args:
            rho_func: Mass density function ρ(x,y,z) [kg/m³]
            Phi_func: Coherence field function Φ(x,y,z) [m⁻¹]
            method: Solver method ('cg' or 'bicgstab')
            tol: Solver tolerance
            maxiter: Maximum iterations (None = auto)
            preconditioner: Preconditioner ('none', 'amg', 'ilu')
        
        Returns:
            PoissonSolution with φ, G_eff, ρ, Φ fields
        """
        print(f"\n3D Poisson Solver")
        print(f"   Grid: {self.grid.nx}×{self.grid.ny}×{self.grid.nz} = {self.grid.total_points:,} points")
        print(f"   Domain: [{-self.grid.Lx/2:.1f}, {self.grid.Lx/2:.1f}]³ m")
        print(f"   Resolution: Δx={self.grid.dx:.3f} m")
        print(f"   Coupling: ξ={self.xi}")
        
        # Evaluate fields on grid
        print(f"   Evaluating ρ(x) and Φ(x)...")
        rho = np.zeros((self.grid.nx, self.grid.ny, self.grid.nz))
        Phi = np.zeros_like(rho)
        
        for i in range(self.grid.nx):
            for j in range(self.grid.ny):
                for k in range(self.grid.nz):
                    x, y, z = self.grid.coord(i, j, k)
                    rho[i,j,k] = rho_func(x, y, z)
                    Phi[i,j,k] = Phi_func(x, y, z)
        
        # Compute G_eff field
        print(f"   Computing G_eff(x)...")
        G_eff = self.compute_G_eff(Phi)
        
        G_eff_min = G_eff.min()
        G_eff_max = G_eff.max()
        print(f"   G_eff range: [{G_eff_min/G_SI:.3e}, {G_eff_max/G_SI:.3e}] × G")
        
        # Build linear system
        A, b = self.build_linear_system(rho, G_eff)
        
        # Build preconditioner
        M = None
        precond_time = 0.0
        
        if preconditioner == 'amg' and HAS_PYAMG:
            print(f"   Building AMG preconditioner...")
            import time
            t0 = time.time()
            ml = pyamg.smoothed_aggregation_solver(A, max_levels=10, max_coarse=100)
            M = ml.aspreconditioner()
            precond_time = time.time() - t0
            print(f"   AMG setup time: {precond_time:.2f} s")
        elif preconditioner == 'ilu' and HAS_SPILU:
            print(f"   Building ILU preconditioner...")
            import time
            t0 = time.time()
            ilu = spilu(A.tocsc(), drop_tol=1e-4, fill_factor=10)
            M_x = lambda x: ilu.solve(x)
            M = LinearOperator(A.shape, M_x)
            precond_time = time.time() - t0
            print(f"   ILU setup time: {precond_time:.2f} s")
        elif preconditioner != 'none':
            print(f"   Warning: Preconditioner '{preconditioner}' not available, using none")
        
        # Solve
        precond_str = f" with {preconditioner.upper()} preconditioner" if preconditioner != 'none' else ""
        print(f"   Solving with {method.upper()}{precond_str}...")
        if maxiter is None:
            maxiter = 10 * self.grid.total_points
        
        import time
        t0 = time.time()
        if method == 'cg':
            phi_flat, info = cg(A, b, M=M, rtol=tol, maxiter=maxiter)
        elif method == 'bicgstab':
            phi_flat, info = bicgstab(A, b, M=M, rtol=tol, maxiter=maxiter)
        else:
            raise ValueError(f"Unknown method: {method}")
        solve_time = time.time() - t0
        
        if info == 0:
            print(f"   ✅ Converged in {solve_time:.2f} s!")
        elif info > 0:
            print(f"   ⚠️ Did not converge in {info} iterations ({solve_time:.2f} s)")
        else:
            print(f"   ❌ Solver error (code {info})")
        
        # Reshape solution
        phi = phi_flat.reshape((self.grid.nx, self.grid.ny, self.grid.nz))
        
        # Compute residual
        residual = np.linalg.norm(A @ phi_flat - b) / np.linalg.norm(b)
        print(f"   Relative residual: {residual:.3e}")
        
        # Solution statistics
        phi_min, phi_max = phi.min(), phi.max()
        print(f"   φ range: [{phi_min:.3e}, {phi_max:.3e}] m²/s²")
        
        solver_info = {
            'method': method,
            'preconditioner': preconditioner,
            'tolerance': tol,
            'maxiter': maxiter,
            'converged': (info == 0),
            'iterations': info if info > 0 else None,
            'solve_time': solve_time,
            'precond_time': precond_time,
            'residual': float(residual),
            'phi_min': float(phi_min),
            'phi_max': float(phi_max),
            'G_eff_min': float(G_eff_min),
            'G_eff_max': float(G_eff_max),
            'xi': self.xi
        }
        
        return PoissonSolution(
            grid=self.grid,
            phi=phi,
            G_eff=G_eff,
            rho=rho,
            Phi_coherence=Phi,
            solver_info=solver_info
        )


# ============================================================================
# Test Cases
# ============================================================================

def spherical_mass_coherent_shell(
    M: float = 1.0,  # Mass [kg]
    R_mass: float = 0.1,  # Mass radius [m]
    R_shell_inner: float = 0.2,  # Shell inner radius [m]
    R_shell_outer: float = 0.5,  # Shell outer radius [m]
    Phi0: float = 1e7  # Shell coherence [m⁻¹]
) -> Tuple[Callable, Callable, Dict]:
    """
    Create test case: point mass with spherical coherent shell.
    
    ρ(r) = (3M)/(4πR³) if r < R_mass, else 0
    Φ(r) = Φ₀ if R_in < r < R_out, else 0
    
    Returns:
        rho_func, Phi_func, metadata
    """
    rho0 = 3.0 * M / (4.0 * np.pi * R_mass**3)
    
    def rho_func(x: float, y: float, z: float) -> float:
        r = np.sqrt(x**2 + y**2 + z**2)
        return rho0 if r < R_mass else 0.0
    
    def Phi_func(x: float, y: float, z: float) -> float:
        r = np.sqrt(x**2 + y**2 + z**2)
        if R_shell_inner < r < R_shell_outer:
            return Phi0
        return 0.0
    
    metadata = {
        'M': M,
        'R_mass': R_mass,
        'R_shell_inner': R_shell_inner,
        'R_shell_outer': R_shell_outer,
        'Phi0': Phi0,
        'rho0': rho0
    }
    
    return rho_func, Phi_func, metadata


def newtonian_potential_sphere(r: float, M: float, R: float) -> float:
    """
    Analytic Newtonian potential for uniform sphere.
    
    φ(r) = -GM/R × (3/2 - r²/2R²)  for r < R
    φ(r) = -GM/r                    for r ≥ R
    """
    if r < R:
        return -G_SI * M / R * (1.5 - 0.5 * (r/R)**2)
    else:
        return -G_SI * M / r


# ============================================================================
# Main Test Script
# ============================================================================

if __name__ == '__main__':
    print("="*70)
    print("3D POISSON SOLVER - TEST CASE")
    print("Spherical mass with coherent shell")
    print("="*70)
    
    # Grid setup
    grid = Grid3D(
        nx=51, ny=51, nz=51,
        Lx=2.0, Ly=2.0, Lz=2.0  # -1 to +1 m cube
    )
    
    # Test case
    rho_func, Phi_func, meta = spherical_mass_coherent_shell(
        M=1.0,             # 1 kg mass
        R_mass=0.05,       # 5 cm radius
        R_shell_inner=0.15,  # Shell 15-40 cm
        R_shell_outer=0.40,
        Phi0=1e7           # Rb BEC scale
    )
    
    print(f"\nTest Configuration:")
    print(f"   Mass: {meta['M']} kg at r < {meta['R_mass']} m")
    print(f"   Coherent shell: Φ₀ = {meta['Phi0']:.2e} m⁻¹")
    print(f"   Shell region: {meta['R_shell_inner']} < r < {meta['R_shell_outer']} m")
    
    # Solve with ξ=100 (realistic coupling)
    solver = Poisson3DSolver(grid, xi=100.0)
    solution = solver.solve(rho_func, Phi_func, method='cg', tol=1e-8)
    
    # Extract radial profile along x-axis
    print(f"\n   Extracting radial profile...")
    j_center = grid.ny // 2
    k_center = grid.nz // 2
    
    r_vals = []
    phi_vals = []
    phi_newtonian = []
    G_eff_vals = []
    
    for i in range(grid.nx):
        x, y, z = grid.coord(i, j_center, k_center)
        r = abs(x)
        r_vals.append(r)
        phi_vals.append(solution.phi[i, j_center, k_center])
        phi_newtonian.append(newtonian_potential_sphere(r, meta['M'], meta['R_mass']))
        G_eff_vals.append(solution.G_eff[i, j_center, k_center])
    
    # Statistics
    print(f"\n   Radial Profile Analysis:")
    phi_vals_arr = np.array(phi_vals)
    phi_newton_arr = np.array(phi_newtonian)
    
    # Find deepest point
    idx_min = np.argmin(phi_vals_arr)
    print(f"   φ_min (modified): {phi_vals_arr[idx_min]:.6e} m²/s² at r={r_vals[idx_min]:.3f} m")
    print(f"   φ_min (Newtonian): {phi_newton_arr[idx_min]:.6e} m²/s²")
    
    ratio = phi_vals_arr[idx_min] / phi_newton_arr[idx_min]
    print(f"   Depth ratio: {ratio:.6f}")
    
    # Save solution
    output_dir = Path('results')
    output_dir.mkdir(exist_ok=True)
    solution.save(output_dir / 'poisson_3d_test.npz')
    
    # Save profile
    profile_data = {
        'r': r_vals,
        'phi_modified': phi_vals,
        'phi_newtonian': phi_newtonian,
        'G_eff': G_eff_vals,
        'metadata': meta,
        'solver_info': solution.solver_info
    }
    
    with open(output_dir / 'radial_profile.json', 'w') as f:
        # Convert numpy arrays to lists for JSON
        profile_json = {
            'r': [float(x) for x in r_vals],
            'phi_modified': [float(x) for x in phi_vals],
            'phi_newtonian': [float(x) for x in phi_newtonian],
            'G_eff_over_G': [float(x/G_SI) for x in G_eff_vals],
            'metadata': meta,
            'solver_info': solution.solver_info
        }
        json.dump(profile_json, f, indent=2)
    
    print(f"\n✅ Test complete!")
    print("="*70)
