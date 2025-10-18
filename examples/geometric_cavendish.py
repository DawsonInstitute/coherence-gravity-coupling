"""
Geometric Cavendish Simulation with 3D Poisson Solver

Models a realistic Cavendish torsion balance experiment with:
- Two source masses (spheres)
- Torsion bar with two test masses
- Coherent body (BEC or superconductor) positioned near sources
- Computes torque from spatial gradient of gravitational potential

This replaces the path-average approximation with full 3D field solution.

Author: GitHub Copilot (Claude Sonnet 4.5)
License: MIT
"""

import numpy as np
import json
from pathlib import Path
import sys
from typing import Dict, Tuple, Optional
import time

# Add parent for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.solvers.poisson_3d import Poisson3DSolver, Grid3D, PoissonSolution
from src.analysis.phi_calibration import get_all_calibrations

# Physical constants
G_SI = 6.674e-11  # m³/(kg·s²)

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches
    from mpl_toolkits.mplot3d import Axes3D
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available, plotting disabled")


class CavendishGeometry:
    """Defines the geometry of a Cavendish torsion balance experiment."""
    
    def __init__(
        self,
        # Source masses (large spheres)
        M_source: float = 1.0,  # kg
        R_source: float = 0.03,  # m (3 cm radius)
        source_separation: float = 0.20,  # m (20 cm between centers)
        
        # Test masses (small spheres on torsion bar)
        m_test: float = 0.010,  # kg (10 g)
        R_test: float = 0.01,  # m (1 cm radius)
        L_torsion: float = 0.10,  # m (10 cm arm length)
        
        # Coherent body (BEC or superconductor block)
        coherent_volume: Tuple[float, float, float] = (0.05, 0.05, 0.05),  # m (5 cm cube)
        coherent_position: Tuple[float, float, float] = (0.0, 0.0, -0.08),  # m
        Phi_coherent: float = 1e7,  # m⁻¹
        
        # Overall positioning
        source_offset_y: float = 0.10  # m (sources offset from torsion bar)
    ):
        """
        Initialize Cavendish geometry.
        
        Coordinate system:
        - Origin at center of torsion bar
        - x-axis: along torsion bar (test masses at ±L/2)
        - y-axis: toward source masses
        - z-axis: vertical (perpendicular to torsion plane)
        
        Standard configuration:
        - Test mass 1: (+L/2, 0, 0)
        - Test mass 2: (-L/2, 0, 0)
        - Source mass 1: (+sep/2, offset_y, 0)
        - Source mass 2: (-sep/2, offset_y, 0)
        - Coherent body: centered at coherent_position
        """
        self.M_source = M_source
        self.R_source = R_source
        self.source_separation = source_separation
        
        self.m_test = m_test
        self.R_test = R_test
        self.L_torsion = L_torsion
        
        self.coherent_volume = coherent_volume
        self.coherent_position = coherent_position
        self.Phi_coherent = Phi_coherent
        
        self.source_offset_y = source_offset_y
        
        # Derived positions
        self.test1_pos = np.array([L_torsion/2, 0.0, 0.0])
        self.test2_pos = np.array([-L_torsion/2, 0.0, 0.0])
        
        self.source1_pos = np.array([source_separation/2, source_offset_y, 0.0])
        self.source2_pos = np.array([-source_separation/2, source_offset_y, 0.0])
    
    def density_function(self, x: float, y: float, z: float) -> float:
        """
        Mass density field ρ(x,y,z) for all masses.
        
        Returns density in kg/m³.
        """
        r = np.array([x, y, z])
        rho = 0.0
        
        # Source mass 1 (uniform sphere)
        r1 = np.linalg.norm(r - self.source1_pos)
        if r1 < self.R_source:
            V1 = (4.0/3.0) * np.pi * self.R_source**3
            rho += self.M_source / V1
        
        # Source mass 2
        r2 = np.linalg.norm(r - self.source2_pos)
        if r2 < self.R_source:
            V2 = (4.0/3.0) * np.pi * self.R_source**3
            rho += self.M_source / V2
        
        # Test mass 1 (small)
        rt1 = np.linalg.norm(r - self.test1_pos)
        if rt1 < self.R_test:
            Vt = (4.0/3.0) * np.pi * self.R_test**3
            rho += self.m_test / Vt
        
        # Test mass 2
        rt2 = np.linalg.norm(r - self.test2_pos)
        if rt2 < self.R_test:
            Vt = (4.0/3.0) * np.pi * self.R_test**3
            rho += self.m_test / Vt
        
        return rho
    
    def coherence_function(self, x: float, y: float, z: float) -> float:
        """
        Coherence field Φ(x,y,z).
        
        Returns Φ in m⁻¹. Non-zero inside coherent body, zero elsewhere.
        """
        r = np.array([x, y, z])
        rc = np.array(self.coherent_position)
        
        # Check if inside coherent volume (rectangular box)
        dx, dy, dz = self.coherent_volume
        
        if (abs(r[0] - rc[0]) < dx/2 and 
            abs(r[1] - rc[1]) < dy/2 and 
            abs(r[2] - rc[2]) < dz/2):
            return self.Phi_coherent
        
        return 0.0
    
    def trilinear_interpolate(self, phi: np.ndarray, grid: Grid3D, pos: Tuple[float, float, float]) -> float:
        """
        Trilinear interpolation of potential φ at arbitrary position.
        
        Args:
            phi: 3D potential array
            grid: Grid3D object
            pos: (x, y, z) position in physical coordinates
            
        Returns:
            Interpolated φ value
        """
        x, y, z = pos
        
        # Convert to grid indices (fractional)
        i_frac = (x + grid.Lx/2) / grid.dx
        j_frac = (y + grid.Ly/2) / grid.dy
        k_frac = (z + grid.Lz/2) / grid.dz
        
        # Lower-left corner indices
        i0 = int(np.floor(i_frac))
        j0 = int(np.floor(j_frac))
        k0 = int(np.floor(k_frac))
        
        # Clamp to valid range
        i0 = max(0, min(grid.nx-2, i0))
        j0 = max(0, min(grid.ny-2, j0))
        k0 = max(0, min(grid.nz-2, k0))
        
        i1 = i0 + 1
        j1 = j0 + 1
        k1 = k0 + 1
        
        # Fractional distances within cell
        dx_frac = i_frac - i0
        dy_frac = j_frac - j0
        dz_frac = k_frac - k0
        
        # Trilinear interpolation
        # Along x at each (j, k) corner
        c00 = phi[i0,j0,k0]*(1-dx_frac) + phi[i1,j0,k0]*dx_frac
        c01 = phi[i0,j0,k1]*(1-dx_frac) + phi[i1,j0,k1]*dx_frac
        c10 = phi[i0,j1,k0]*(1-dx_frac) + phi[i1,j1,k0]*dx_frac
        c11 = phi[i0,j1,k1]*(1-dx_frac) + phi[i1,j1,k1]*dx_frac
        
        # Along y
        c0 = c00*(1-dy_frac) + c10*dy_frac
        c1 = c01*(1-dy_frac) + c11*dy_frac
        
        # Along z
        phi_interp = c0*(1-dz_frac) + c1*dz_frac
        
        return phi_interp
    
    def gradient_trilinear(self, phi: np.ndarray, grid: Grid3D, pos: Tuple[float, float, float]) -> np.ndarray:
        """
        Compute gradient ∇φ using trilinear interpolation with finite differences.
        
        Args:
            phi: 3D potential array
            grid: Grid3D object
            pos: (x, y, z) position in physical coordinates
            
        Returns:
            3D gradient vector [∂φ/∂x, ∂φ/∂y, ∂φ/∂z]
        """
        # Use small offset for numerical derivative
        h = min(grid.dx, grid.dy, grid.dz) * 0.5
        
        x, y, z = pos
        
        # Central differences using interpolation
        grad_x = (self.trilinear_interpolate(phi, grid, (x+h, y, z)) - 
                  self.trilinear_interpolate(phi, grid, (x-h, y, z))) / (2*h)
        
        grad_y = (self.trilinear_interpolate(phi, grid, (x, y+h, z)) - 
                  self.trilinear_interpolate(phi, grid, (x, y-h, z))) / (2*h)
        
        grad_z = (self.trilinear_interpolate(phi, grid, (x, y, z+h)) - 
                  self.trilinear_interpolate(phi, grid, (x, y, z-h))) / (2*h)
        
        return np.array([grad_x, grad_y, grad_z])
    
    def volume_average_force(self, phi: np.ndarray, grid: Grid3D, center: Tuple[float, float, float], radius: float, n_samples: int = 5) -> np.ndarray:
        """
        Compute volume-averaged force over a spherical test mass.
        
        Uses Simpson-like quadrature with trilinear-interpolated gradients.
        Reduces grid aliasing compared to point-sample.
        
        Args:
            phi: Potential field
            grid: Grid3D object
            center: Center of test mass
            radius: Radius of test mass sphere
            n_samples: Number of radial samples (odd for Simpson's rule)
            
        Returns:
            Average force vector F_avg = -m_test * <∇φ>_volume
        """
        if n_samples % 2 == 0:
            n_samples += 1  # Ensure odd for Simpson's rule
        
        # Radial samples from 0 to radius
        r_samples = np.linspace(0, radius, n_samples)
        
        # Angular samples (simplified: use theta-phi grid)
        n_theta = max(3, n_samples)
        n_phi = max(5, 2 * n_samples)
        
        theta = np.linspace(0, np.pi, n_theta)
        phi_angle = np.linspace(0, 2*np.pi, n_phi, endpoint=False)
        
        force_sum = np.zeros(3)
        weight_sum = 0.0
        
        for ir, r in enumerate(r_samples):
            # Simpson weights for radial integration
            if ir == 0 or ir == n_samples - 1:
                w_r = 1.0
            elif ir % 2 == 1:
                w_r = 4.0
            else:
                w_r = 2.0
            
            # Surface integral at this radius
            for th in theta:
                for ph in phi_angle:
                    # Spherical to Cartesian
                    x = center[0] + r * np.sin(th) * np.cos(ph)
                    y = center[1] + r * np.sin(th) * np.sin(ph)
                    z = center[2] + r * np.cos(th)
                    
                    pos = (x, y, z)
                    
                    # Get gradient at this point
                    try:
                        grad_phi_local = self.gradient_trilinear(phi, grid, pos)
                    except (IndexError, ValueError):
                        # Out of bounds, skip
                        continue
                    
                    # Volume element: r² sin(θ) dr dθ dφ
                    dtheta = np.pi / (n_theta - 1) if n_theta > 1 else np.pi
                    dphi = 2*np.pi / n_phi
                    dV_weight = r**2 * np.sin(th) * w_r * dtheta * dphi
                    
                    force_sum += -grad_phi_local * dV_weight
                    weight_sum += dV_weight
        
        # Normalize by total weight (approximate volume integral)
        if weight_sum > 0:
            force_avg = force_sum / weight_sum
        else:
            # Fallback to point sample
            grad_center = self.gradient_trilinear(phi, grid, center)
            force_avg = -grad_center
        
        # Scale by mass
        return self.m_test * force_avg
    
    def compute_torque(self, solution: PoissonSolution, use_interpolation: bool = True, use_volume_average: bool = False) -> Dict:
        """
        Compute gravitational torque on torsion bar from 3D potential field.
        
        Torque τ = r × F where F = -m ∇φ
        
        For torsion bar rotating about z-axis:
        τ_z = x * F_y - y * F_x
        
        where F = -m ∇φ evaluated at each test mass position.
        
        Args:
            solution: PoissonSolution with phi field
            use_interpolation: If True, use trilinear interpolation; otherwise nearest-grid
            use_volume_average: If True, volume-average force over test mass (reduces aliasing)
        """
        grid = solution.grid
        phi = solution.phi
        
        if use_volume_average:
            # Volume-averaged force (preferred for convergence)
            F1 = self.volume_average_force(phi, grid, self.test1_pos, self.R_test)
            F2 = self.volume_average_force(phi, grid, self.test2_pos, self.R_test)
            
            phi1 = self.trilinear_interpolate(phi, grid, self.test1_pos)
            phi2 = self.trilinear_interpolate(phi, grid, self.test2_pos)
        elif use_interpolation:
            # Use trilinear interpolation for gradient evaluation
            grad_phi1 = self.gradient_trilinear(phi, grid, self.test1_pos)
            grad_phi2 = self.gradient_trilinear(phi, grid, self.test2_pos)
            
            F1 = -self.m_test * grad_phi1
            F2 = -self.m_test * grad_phi2
            
            # Get interpolated phi values
            phi1 = self.trilinear_interpolate(phi, grid, self.test1_pos)
            phi2 = self.trilinear_interpolate(phi, grid, self.test2_pos)
        else:
            # Original nearest-grid method
            def nearest_index(pos):
                i = int((pos[0] + grid.Lx/2) / grid.dx)
                j = int((pos[1] + grid.Ly/2) / grid.dy)
                k = int((pos[2] + grid.Lz/2) / grid.dz)
                
                # Clamp to valid range
                i = max(1, min(grid.nx-2, i))
                j = max(1, min(grid.ny-2, j))
                k = max(1, min(grid.nz-2, k))
                
                return i, j, k
            
            i1, j1, k1 = nearest_index(self.test1_pos)
            i2, j2, k2 = nearest_index(self.test2_pos)
            
            # Compute force at test mass 1: F1 = -m1 * ∇φ
            grad_phi_x1 = (phi[i1+1,j1,k1] - phi[i1-1,j1,k1]) / (2*grid.dx)
            grad_phi_y1 = (phi[i1,j1+1,k1] - phi[i1,j1-1,k1]) / (2*grid.dy)
            grad_phi_z1 = (phi[i1,j1,k1+1] - phi[i1,j1,k1-1]) / (2*grid.dz)
            
            F1 = -self.m_test * np.array([grad_phi_x1, grad_phi_y1, grad_phi_z1])
            
            # Compute force at test mass 2
            grad_phi_x2 = (phi[i2+1,j2,k2] - phi[i2-1,j2,k2]) / (2*grid.dx)
            grad_phi_y2 = (phi[i2,j2+1,k2] - phi[i2,j2-1,k2]) / (2*grid.dy)
            grad_phi_z2 = (phi[i2,j2,k2+1] - phi[i2,j2,k2-1]) / (2*grid.dz)
            
            F2 = -self.m_test * np.array([grad_phi_x2, grad_phi_y2, grad_phi_z2])
            
            phi1 = phi[i1,j1,k1]
            phi2 = phi[i2,j2,k2]
        
        # Torque about z-axis (rotation axis)
        # τ = r × F, z-component only
        r1 = self.test1_pos
        r2 = self.test2_pos
        
        tau1_z = r1[0] * F1[1] - r1[1] * F1[0]
        tau2_z = r2[0] * F2[1] - r2[1] * F2[0]
        
        tau_total = tau1_z + tau2_z
        
        return {
            'torque_total': tau_total,
            'torque_mass1': tau1_z,
            'torque_mass2': tau2_z,
            'force_mass1': F1,
            'force_mass2': F2,
            'phi_mass1': phi1,
            'phi_mass2': phi2,
            'interpolation_used': use_interpolation
        }
    
    def to_dict(self):
        """Return dictionary representation of geometry parameters."""
        return {
            'M_source': self.M_source,
            'R_source': self.R_source,
            'source_separation': self.source_separation,
            'm_test': self.m_test,
            'R_test': self.R_test,
            'L_torsion': self.L_torsion,
            'coherent_volume': self.coherent_volume,
            'coherent_position': self.coherent_position,
            'Phi_coherent': self.Phi_coherent,
            'source_offset_y': self.source_offset_y,
        }


def run_geometric_cavendish(
    xi: float,
    Phi0: float,
    geom_params: Optional[Dict] = None,
    grid_resolution: int = 41,
    domain_size: float = 0.6,
    verbose: bool = True,
    use_interpolation: bool = True,
    use_volume_average: bool = False,
    grid_nx: Optional[int] = None,
    solver_method: str = 'cg',
    preconditioner: str = 'diagonal',
    cache: bool = False,
) -> Dict:
    """
    Run full geometric Cavendish simulation.
    
    Args:
        xi: Non-minimal coupling strength
        Phi0: Coherence field amplitude [m⁻¹]
        geom_params: Dictionary of geometry parameters to override defaults
        grid_resolution: Number of grid points per dimension
        domain_size: Total domain size [m]
        verbose: Print progress
        use_interpolation: Use trilinear interpolation for force
        use_volume_average: Use volume-averaged force (reduces aliasing)
        grid_nx: Alias for grid_resolution (for test compatibility)
        solver_method: Iterative solver ('cg' or 'bicgstab')
        preconditioner: Preconditioner type ('none', 'diagonal', 'amg', 'ilu')
        cache: Enable result caching (default: False)
    
    Returns:
        Dictionary with torque, ΔG/G, timing, etc.
    """
    # Handle optional parameter aliases
    if grid_nx is not None:
        grid_resolution = grid_nx
    
    # Initialize geometry parameters, allowing overrides
    base_geom_params = {
        'Phi_coherent': Phi0
    }
    if geom_params:
        base_geom_params.update(geom_params)
    
    # Check cache if enabled
    if cache:
        from src.utils.result_cache import get_cache
        cache_inst = get_cache()
        cache_key = cache_inst.compute_key(
            xi=xi,
            Phi0=Phi0,
            geom_params=base_geom_params,
            grid_resolution=grid_resolution,
            domain_size=domain_size,
            solver_method=solver_method,
            preconditioner=preconditioner
        )
        cached = cache_inst.load(cache_key)
        if cached is not None:
            if verbose:
                print(f"✅ Cache HIT: {cache_key}")
            return cached['result']
        elif verbose:
            print(f"⚠️  Cache MISS: {cache_key}")

    if verbose:
        print(f"\n{'='*70}")
        print(f"GEOMETRIC CAVENDISH SIMULATION")
        print(f"{'='*70}")
        print(f"   ξ = {xi}")
        print(f"   Φ₀ = {Phi0:.2e} m⁻¹")
        for key, value in base_geom_params.items():
            if key != 'Phi_coherent':
                print(f"   {key}: {value}")

    geom = CavendishGeometry(**base_geom_params)
    
    # Grid (needs to encompass all masses)
    grid = Grid3D(
        nx=grid_resolution,
        ny=grid_resolution,
        nz=grid_resolution,
        Lx=domain_size,
        Ly=domain_size,
        Lz=domain_size
    )
    
    # Solve with coherence
    if verbose:
        print(f"\n--- Solving with coherence (Φ₀ = {Phi0:.2e}) ---")
    
    t0 = time.time()
    solver = Poisson3DSolver(grid, xi=xi)
    solution_coherent = solver.solve(
        geom.density_function,
        geom.coherence_function,
        method=solver_method,
        preconditioner=preconditioner,
        tol=1e-8
    )
    t_coherent = time.time() - t0
    
    if verbose:
        print(f"   Solve time: {t_coherent:.2f} s")
    
    # Compute torque with specified method
    torque_result = geom.compute_torque(solution_coherent, use_interpolation=use_interpolation, use_volume_average=use_volume_average)
    tau_coherent = torque_result['torque_total']
    
    if verbose:
        print(f"   Torque: {tau_coherent:.6e} N·m")
    
    # Solve without coherence (Newtonian baseline)
    if verbose:
        print(f"\n--- Solving without coherence (Φ = 0) ---")
    
    def zero_coherence(x, y, z):
        return 0.0
    
    t0 = time.time()
    solution_newtonian = solver.solve(
        geom.density_function,
        zero_coherence,
        method=solver_method,
        preconditioner=preconditioner,
        tol=1e-8
    )
    t_newtonian = time.time() - t0
    
    if verbose:
        print(f"   Solve time: {t_newtonian:.2f} s")
    
    torque_newton = geom.compute_torque(solution_newtonian, use_interpolation=use_interpolation, use_volume_average=use_volume_average)
    tau_newtonian = torque_newton['torque_total']
    
    if verbose:
        print(f"   Torque: {tau_newtonian:.6e} N·m")
    
    # Compute fractional change
    delta_tau = tau_coherent - tau_newtonian
    delta_tau_frac = delta_tau / tau_newtonian if tau_newtonian != 0 else 0.0
    
    # Estimate effective ΔG/G from torque ratio
    # For small perturbations: Δτ/τ ≈ ΔG/G
    delta_G_over_G = delta_tau_frac
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"RESULTS")
        print(f"{'='*70}")
        print(f"   Torque (Newtonian): {tau_newtonian:.6e} N·m")
        print(f"   Torque (with coherence): {tau_coherent:.6e} N·m")
        print(f"   Δτ: {delta_tau:.6e} N·m")
        print(f"   Δτ/τ: {delta_tau_frac:.6e}")
        print(f"   ΔG_eff/G (estimated): {delta_G_over_G:.6e}")
        print(f"{'='*70}\n")
    
    result = {
        'xi': xi,
        'Phi0': Phi0,
        'geom_params': geom.to_dict(),
        'grid_resolution': grid_resolution,
        'tau_newtonian': float(tau_newtonian),
        'tau_coherent': float(tau_coherent),
        'delta_tau': float(delta_tau),
        'delta_tau_frac': float(delta_tau_frac),
        'delta_G_over_G': float(delta_G_over_G),
        'solve_time_coherent': t_coherent,
        'solve_time_newtonian': t_newtonian,
        'torque_details_coherent': {k: float(v) if np.isscalar(v) else v.tolist() 
                                     for k, v in torque_result.items()},
        'torque_details_newtonian': {k: float(v) if np.isscalar(v) else v.tolist() 
                                      for k, v in torque_newton.items()}
    }
    
    # Save to cache if enabled
    if cache:
        cache_inst.save(
            cache_key,
            result,
            phi_coherent=solution_coherent,
            phi_newtonian=solution_newtonian
        )
        if verbose:
            print(f"💾 Saved to cache: {cache_key}")
    
    return result


# ============================================================================
# Geometry Optimization Functions
# ============================================================================

def sweep_test_mass(
    m_test_range: np.ndarray,
    base_geom_params: Dict,
    xi: float = 100.0,
    Phi0: float = 6.67e8,
    fiber_stress_limit: float = 1e9,  # Pa (typical tungsten wire, ~1 GPa)
    verbose: bool = True
) -> Dict:
    """
    Sweep test mass to find optimal value maximizing |Δτ| while respecting fiber limits.
    
    Args:
        m_test_range: Array of test masses to sweep [kg]
        base_geom_params: Dictionary of baseline geometry parameters
        xi: Coherence coupling strength
        Phi0: Coherence field amplitude [m⁻¹]
        fiber_stress_limit: Maximum fiber stress [Pa]
        verbose: Print progress
        
    Returns:
        Dict with sweep results and optimal configuration
    """
    if verbose:
        print("\n" + "="*70)
        print("SWEEPING TEST MASS")
        print("="*70)
    
    results = []
    r_fiber = 25e-6  # m (typical torsion fiber radius)
    A_fiber = np.pi * r_fiber**2
    g = 9.81  # m/s²
    
    for m_test in m_test_range:
        # Check fiber stress constraint
        weight_force = m_test * g
        stress = weight_force / A_fiber
        
        if stress > fiber_stress_limit:
            if verbose:
                print(f"  m_test = {m_test*1e3:.2f} g: SKIPPED (stress {stress/1e6:.1f} MPa > limit)")
            continue
        
        # Run simulation with this test mass
        current_geom = base_geom_params.copy()
        current_geom['m_test'] = m_test
        
        result_sim = run_geometric_cavendish(
            xi=xi,
            Phi0=Phi0,
            geom_params=current_geom,
            grid_resolution=41,
            verbose=False
        )
        
        delta_tau = abs(result_sim['delta_tau'])
        
        results.append({
            'm_test': m_test,
            'm_test_mg': m_test * 1e6,  # Convert to mg
            'fiber_stress_MPa': stress / 1e6,
            'tau_newtonian': result_sim['tau_newtonian'],
            'tau_coherent': result_sim['tau_coherent'],
            'delta_tau': delta_tau,
            'delta_tau_frac': result_sim['delta_tau_frac']
        })
        
        if verbose:
            print(f"  m_test = {m_test*1e3:.2f} g: Δτ = {delta_tau:.3e} N·m "
                  f"(stress {stress/1e6:.1f} MPa)")
    
    # Find optimal (maximum |Δτ|)
    if len(results) == 0:
        raise ValueError("No valid test masses found within fiber stress limits")
    
    optimal_idx = np.argmax([r['delta_tau'] for r in results])
    optimal = results[optimal_idx]
    
    if verbose:
        print(f"\n✅ Optimal test mass: {optimal['m_test']*1e3:.2f} g")
        print(f"   Δτ = {optimal['delta_tau']:.3e} N·m")
        print(f"   Fiber stress: {optimal['fiber_stress_MPa']:.1f} MPa\n")
    
    return {
        'sweep_results': results,
        'optimal': optimal,
        'parameters': {
            'xi': xi,
            'Phi0': Phi0,
            'base_geom_params': base_geom_params,
            'fiber_stress_limit_MPa': fiber_stress_limit / 1e6
        }
    }


def sweep_source_mass(
    M_source_range: np.ndarray,
    base_geom_params: Dict,
    xi: float = 100.0,
    Phi0: float = 6.67e8,
    domain_size: float = 0.6,
    verbose: bool = True
) -> Dict:
    """
    Sweep source mass to maximize signal while keeping sources within domain.
    
    Args:
        M_source_range: Array of source masses to sweep [kg]
        base_geom_params: Dictionary of baseline geometry parameters
        xi: Coherence coupling strength
        Phi0: Coherence field amplitude [m⁻¹]
        domain_size: Computational domain size [m]
        verbose: Print progress
        
    Returns:
        Dict with sweep results and optimal configuration
    """
    if verbose:
        print("\n" + "="*70)
        print("SWEEPING SOURCE MASS")
        print("="*70)
    
    results = []
    
    # Get default values from a dummy instance
    default_geom = CavendishGeometry()
    M_default = default_geom.M_source
    R_source_default = default_geom.R_source
    
    for M_source in M_source_range:
        # Scale radius assuming constant density: R ∝ M^(1/3)
        R_source = R_source_default * (M_source / M_default)**(1/3)
        
        # Check if sources fit in domain with margin
        sep = base_geom_params.get('source_separation', default_geom.source_separation)
        max_extent = sep/2 + R_source
        
        if max_extent > domain_size/2 * 0.9:  # 90% of domain for safety
            if verbose:
                print(f"  M_source = {M_source:.2f} kg: SKIPPED (R={R_source*100:.1f}cm, too large for domain)")
            continue
        
        # Run simulation with this source mass and scaled radius
        current_geom = base_geom_params.copy()
        current_geom['M_source'] = M_source
        current_geom['R_source'] = R_source
        
        result_sim = run_geometric_cavendish(
            xi=xi,
            Phi0=Phi0,
            geom_params=current_geom,
            grid_resolution=41,
            domain_size=domain_size,
            verbose=False
        )
        
        delta_tau = abs(result_sim['delta_tau'])
        
        results.append({
            'M_source': M_source,
            'R_source': R_source,
            'tau_newtonian': result_sim['tau_newtonian'],
            'tau_coherent': result_sim['tau_coherent'],
            'delta_tau': delta_tau,
            'delta_tau_frac': result_sim['delta_tau_frac']
        })
        
        if verbose:
            print(f"  M_source = {M_source:.2f} kg (R={R_source*100:.1f}cm): Δτ = {delta_tau:.3e} N·m")
    
    if len(results) == 0:
        raise ValueError("No valid source masses found within domain constraints")
    
    optimal_idx = np.argmax([r['delta_tau'] for r in results])
    optimal = results[optimal_idx]
    
    if verbose:
        print(f"\n✅ Optimal source mass: {optimal['M_source']:.2f} kg")
        print(f"   Δτ = {optimal['delta_tau']:.3e} N·m\n")
    
    return {
        'sweep_results': results,
        'optimal': optimal,
        'parameters': {
            'xi': xi,
            'Phi0': Phi0,
            'base_geom_params': base_geom_params,
            'domain_size': domain_size
        }
    }


def sweep_coherent_position(
    y_range: np.ndarray,
    z_range: np.ndarray,
    base_geom_params: Dict,
    xi: float = 100.0,
    Phi0: float = 6.67e8,
    verbose: bool = True
) -> Dict:
    """
    Sweep coherent body position to find maximum |Δτ|.
    
    Args:
        y_range: Array of y positions to sweep [m]
        z_range: Array of z positions to sweep [m]
        base_geom_params: Dictionary of baseline geometry parameters
        xi: Coherence coupling strength
        Phi0: Coherence field amplitude [m⁻¹]
        verbose: Print progress
        
    Returns:
        Dict with sweep results and optimal position
    """
    if verbose:
        print("\n" + "="*70)
        print("SWEEPING COHERENT BODY POSITION")
        print("="*70)
    
    results = []
    
    # Fixed x=0 (centered on torsion bar)
    for y_pos in y_range:
        for z_pos in z_range:
            current_geom = base_geom_params.copy()
            current_geom['coherent_position'] = (0.0, y_pos, z_pos)
            
            # Run simulation
            result = run_geometric_cavendish(
                xi=xi,
                Phi0=Phi0,
                geom_params=current_geom,
                grid_resolution=41,
                verbose=False
            )
            
            results.append({
                'position': current_geom['coherent_position'],
                'y': y_pos,
                'z': z_pos,
                'tau_newtonian': result['tau_newtonian'],
                'tau_coherent': result['tau_coherent'],
                'delta_tau': result['delta_tau'],
                'delta_tau_frac': result['delta_tau_frac'],
                'delta_G_over_G': result['delta_G_over_G']
            })
            
            if verbose:
                print(f"  Position (0, {y_pos:.3f}, {z_pos:.3f}): "
                      f"Δτ = {result['delta_tau']:.3e} N·m, "
                      f"ΔG/G = {result['delta_G_over_G']:.3f}")
    
    optimal_idx = np.argmax([abs(r['delta_tau']) for r in results])
    optimal = results[optimal_idx]
    
    if verbose:
        print(f"\n✅ Optimal position: (0, {optimal['y']:.3f}, {optimal['z']:.3f}) m")
        print(f"   Δτ = {optimal['delta_tau']:.3e} N·m")
        print(f"   ΔG/G = {optimal['delta_G_over_G']:.3f}\n")
    
    return {
        'sweep_results': results,
        'optimal': optimal,
        'parameters': {
            'xi': xi,
            'Phi0': Phi0,
            'base_geom_params': base_geom_params,
            'y_range': y_range.tolist(),
            'z_range': z_range.tolist()
        }
    }


def convergence_test(
    grid_resolutions: list = [41, 61, 81],
    geom_params: Optional[Dict] = None,
    xi: float = 100.0,
    Phi0: float = 6.67e8,
    use_interpolation: bool = True,
    verbose: bool = True
) -> Dict:
    """
    Test grid convergence: compare torque results across different resolutions.
    
    Args:
        grid_resolutions: List of grid sizes to test (e.g., [41, 61, 81])
        geom_params: Dictionary of geometry parameters
        xi: Coherence coupling strength
        Phi0: Coherence field amplitude [m⁻¹]
        use_interpolation: Use trilinear interpolation (recommended)
        verbose: Print progress
        
    Returns:
        Dict with convergence results
    """
    if verbose:
        print("\n" + "="*70)
        print("GRID CONVERGENCE TEST")
        print(f"Resolutions: {grid_resolutions}")
        print(f"Interpolation: {use_interpolation}")
        print("="*70)
    
    results = []
    
    for resolution in grid_resolutions:
        if verbose:
            print(f"\n--- Testing {resolution}³ grid ---")
        
        result_full = run_geometric_cavendish(
            xi=xi,
            Phi0=Phi0,
            geom_params=geom_params,
            grid_resolution=resolution,
            verbose=False
        )
        
        results.append({
            'resolution': resolution,
            'grid_points': resolution**3,
            'tau_newtonian': result_full['tau_newtonian'],
            'tau_coherent': result_full['tau_coherent'],
            'delta_tau': result_full['delta_tau'],
            'delta_tau_abs': abs(result_full['delta_tau'])
        })
        
        if verbose:
            print(f"  τ_N = {result_full['tau_newtonian']:.6e} N·m")
            print(f"  τ_coh = {result_full['tau_coherent']:.6e} N·m")
            print(f"  Δτ = {result_full['delta_tau']:.6e} N·m")
    
    # Compute convergence rates
    convergence_data = []
    for i in range(1, len(results)):
        coarse = results[i-1]
        fine = results[i]
        
        # Relative difference in Δτ
        rel_diff = abs(fine['delta_tau'] - coarse['delta_tau']) / abs(fine['delta_tau']) if fine['delta_tau'] != 0 else 0
        
        # Grid refinement ratio
        h_ratio = coarse['resolution'] / fine['resolution']
        
        # Estimate convergence order: log(error_ratio) / log(h_ratio)
        if i > 1:
            prev_fine = results[i-1]
            prev_coarse = results[i-2]
            prev_diff = abs(prev_fine['delta_tau'] - prev_coarse['delta_tau']) / abs(prev_fine['delta_tau']) if prev_fine['delta_tau'] != 0 else 0
            if prev_diff > 0 and rel_diff > 0:
                convergence_order = np.log(prev_diff / rel_diff) / np.log(h_ratio)
            else:
                convergence_order = None
        else:
            convergence_order = None
        
        convergence_data.append({
            'from_resolution': coarse['resolution'],
            'to_resolution': fine['resolution'],
            'relative_difference': rel_diff,
            'convergence_order': convergence_order
        })
        
        if verbose:
            print(f"\n{coarse['resolution']}³ → {fine['resolution']}³:")
            print(f"  Relative Δτ difference: {rel_diff:.4e}")
            if convergence_order is not None:
                print(f"  Estimated convergence order: {convergence_order:.2f}")
    
    # Summary
    finest = results[-1]
    if verbose:
        print("\n" + "="*70)
        print("CONVERGENCE SUMMARY")
        print("="*70)
        print(f"Finest grid ({finest['resolution']}³):")
        print(f"  Δτ = {finest['delta_tau']:.6e} N·m")
        if len(convergence_data) > 0:
            final_conv = convergence_data[-1]
            print(f"  Relative convergence: {final_conv['relative_difference']:.4e}")
            if final_conv['relative_difference'] < 0.01:
                print("  ✅ Converged (<1% change)")
            else:
                print(f"  ⚠️  Not fully converged (>{final_conv['relative_difference']*100:.2f}% change)")
    
    return {
        'resolution_sweep': results,
        'convergence_data': convergence_data,
        'parameters': {
            'xi': xi,
            'Phi0': Phi0,
            'geom_params': geom_params,
            'use_interpolation': use_interpolation
        }
    }


def optimize_geometry(
    xi: float = 100.0,
    Phi0: float = 6.67e8,
    verbose: bool = True
) -> Dict:
    """
    Combined geometry optimization: sweep coherent position, test mass, and source mass.
    
    Performs sequential optimization:
    1. Sweep coherent position (fastest, biggest impact)
    2. Sweep test mass (constrained by fiber stress)
    3. Sweep source mass (constrained by domain size)
    
    Args:
        xi: Coherence coupling strength
        Phi0: Coherence field amplitude [m⁻¹]
        verbose: Print progress
        
    Returns:
        Dict with all sweep results and final optimized configuration
    """
    if verbose:
        print("\n" + "="*70)
        print("COMBINED GEOMETRY OPTIMIZATION")
        print(f"ξ = {xi}, Φ₀ = {Phi0:.2e} m⁻¹")
        print("="*70)
    
    base_geom = {}

    # Step 1: Optimize coherent position
    position_sweep = sweep_coherent_position(
        y_range=np.linspace(-0.05, 0.10, 4),
        z_range=np.linspace(-0.12, 0.0, 4),
        base_geom_params=base_geom,
        xi=xi,
        Phi0=Phi0,
        verbose=verbose
    )
    base_geom['coherent_position'] = position_sweep['optimal']['position']
    
    # Step 2: Optimize test mass
    mass_sweep = sweep_test_mass(
        m_test_range=np.linspace(0.005, 0.020, 5),  # 5-20 g
        base_geom_params=base_geom,
        xi=xi,
        Phi0=Phi0,
        verbose=verbose
    )
    base_geom['m_test'] = mass_sweep['optimal']['m_test']
    
    # Step 3: Optimize source mass
    source_sweep = sweep_source_mass(
        M_source_range=np.linspace(0.5, 2.0, 4),  # 0.5-2 kg
        base_geom_params=base_geom,
        xi=xi,
        Phi0=Phi0,
        verbose=verbose
    )
    base_geom['M_source'] = source_sweep['optimal']['M_source']
    base_geom['R_source'] = source_sweep['optimal']['R_source']

    # Final configuration
    final_config = {
        'xi': xi,
        'Phi0': Phi0,
        'geom_params': base_geom,
        'delta_tau': source_sweep['optimal']['delta_tau']
    }
    
    if verbose:
        print("\n" + "="*70)
        print("FINAL OPTIMIZED CONFIGURATION")
        print("="*70)
        for key, value in base_geom.items():
            print(f"  {key}: {value}")
        print(f"  Δτ (optimized): {final_config['delta_tau']:.3e} N·m")
        print("="*70)
    
    return {
        'position_sweep': position_sweep,
        'mass_sweep': mass_sweep,
        'source_sweep': source_sweep,
        'final_config': final_config
    }


def parameter_sweep(
    xi_values: list = [1.0, 10.0, 100.0],
    calibration_names: list = ['rb87_bec', 'nb_cavity', 'ybco_cuprate'],
    coherent_positions: list = [(0.0, 0.0, -0.08), (0.0, 0.05, 0.0)],
    output_dir: Path = Path('results')
) -> Dict:
    """
    Sweep over (ξ, Φ₀, geometry) parameter space.
    
    Args:
        xi_values: List of coupling values
        calibration_names: List of calibration preset names
        coherent_positions: List of coherent body positions
        output_dir: Output directory
    
    Returns:
        Dictionary with all results
    """
    print(f"\n{'='*70}")
    print(f"GEOMETRIC CAVENDISH PARAMETER SWEEP")
    print(f"{'='*70}")
    
    calibrations = get_all_calibrations()
    results = []
    
    total_runs = len(xi_values) * len(calibration_names) * len(coherent_positions)
    run_count = 0
    
    for xi in xi_values:
        for cal_name in calibration_names:
            if cal_name not in calibrations:
                print(f"Warning: {cal_name} not in calibrations, skipping")
                continue
            
            Phi0 = calibrations[cal_name].Phi0
            
            for coh_pos in coherent_positions:
                run_count += 1
                print(f"\n[{run_count}/{total_runs}] ξ={xi}, {cal_name}, pos={coh_pos}")
                
                geom_p = {'coherent_position': coh_pos}
                
                result = run_geometric_cavendish(
                    xi=xi,
                    Phi0=Phi0,
                    geom_params=geom_p,
                    grid_resolution=41,
                    verbose=False
                )
                
                result['calibration_name'] = cal_name
                result['run_id'] = run_count
                
                print(f"   → ΔG/G = {result['delta_G_over_G']:.4e}, "
                      f"Δτ/τ = {result['delta_tau_frac']:.4e}")
                
                results.append(result)
    
    # Save results
    output_dir.mkdir(exist_ok=True)
    output_file = output_dir / 'geometric_cavendish_sweep.json'
    
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n✅ Sweep complete! Results saved to {output_file}")
    
    return {
        'results': results,
        'summary': {
            'total_runs': total_runs,
            'xi_values': xi_values,
            'calibrations': calibration_names,
            'positions': coherent_positions
        }
    }


def plot_results(results: Dict, output_dir: Path = Path('results')):
    """Generate plots from parameter sweep."""
    if not HAS_MATPLOTLIB:
        print("Matplotlib not available, skipping plots")
        return
    
    sweep_results = results['results']
    
    # Extract data
    xi_vals = [r['xi'] for r in sweep_results]
    Phi_vals = [r['Phi0'] for r in sweep_results]
    delta_G = [r['delta_G_over_G'] for r in sweep_results]
    cal_names = [r['calibration_name'] for r in sweep_results]
    
    # Plot 1: ΔG/G vs ξ for each calibration
    fig, ax = plt.subplots(figsize=(10, 6))
    
    unique_cals = list(set(cal_names))
    colors = plt.cm.viridis(np.linspace(0, 1, len(unique_cals)))
    
    for cal, color in zip(unique_cals, colors):
        mask = [c == cal for c in cal_names]
        xi_cal = [xi_vals[i] for i, m in enumerate(mask) if m]
        dG_cal = [delta_G[i] for i, m in enumerate(mask) if m]
        
        ax.scatter(xi_cal, dG_cal, label=cal, s=100, alpha=0.7, color=color)
    
    ax.set_xlabel('Coupling ξ', fontsize=12)
    ax.set_ylabel('ΔG_eff/G', fontsize=12)
    ax.set_title('Geometric Cavendish: ΔG/G vs Coupling', fontsize=14, fontweight='bold')
    ax.set_xscale('log')
    ax.axhline(0, color='k', linestyle='--', alpha=0.3)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'geometric_cavendish_delta_G.png', dpi=200)
    print(f"   Saved: geometric_cavendish_delta_G.png")
    plt.close()


# ============================================================================
# Main Script
# ============================================================================

if __name__ == '__main__':
    # Single test run
    print("Running single geometric Cavendish test...")
    
    result = run_geometric_cavendish(
        xi=100.0,
        Phi0=3.65e6,  # Rb BEC
        geom_params={
            'coherent_volume': (0.05, 0.05, 0.05),
            'coherent_position': (0.0, 0.0, -0.08),
        },
        grid_resolution=41,
        verbose=True
    )
    
    # Parameter sweep
    print("\n" + "="*70)
    print("Running parameter sweep...")
    print("="*70)
    
    sweep_results = parameter_sweep(
        xi_values=[1.0, 10.0, 100.0],
        calibration_names=['rb87_bec', 'nb_cavity', 'ybco_cuprate'],
        coherent_positions=[
            (0.0, 0.0, -0.08),  # Below torsion bar
            (0.0, 0.05, 0.0)    # Between source and test masses
        ],
        output_dir=Path('results')
    )
    
    # Plot
    plot_results(sweep_results, output_dir=Path('results'))
    
    print("\n✅ Geometric Cavendish simulation complete!")
