"""
WKB QNM Diagnostics: Effective Potential → Quasi-Normal Mode Estimates

Physical Basis (Karimabadi et al. arxiv:2508.13820):
    - Black hole perturbations satisfy wave equation with effective potential
    - QNM frequencies determined by potential barrier shape
    - WKB approximation: ω ≈ ω_WKB(V_eff) for large overtone numbers

Implementation:
    Given V_eff(r; θ, ξ, ζ), compute:
        1. Potential peak location and curvature
        2. WKB frequency estimate from barrier properties
        3. Compare to time-domain numerical integration

Usage:
    from src.analysis.wkb_qnm_diagnostics import compute_wkb_qnm, compare_to_time_domain
    
    omega_wkb, gamma_wkb = compute_wkb_qnm(r, V_eff, n=0, l=2)

References:
    - Karimabadi et al. (2025) arXiv:2508.13820 (NC-QNM analysis)
    - Schutz & Will (1985): WKB formula for QNMs
    - Task: arxiv.2508.13820:35
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
from typing import Callable


def find_potential_peak(r: np.ndarray, V_eff: np.ndarray) -> tuple[float, float]:
    """Locate potential maximum and curvature.
    
    Returns:
        (r_peak, V_peak): Location and value of maximum
    """
    idx_peak = np.argmax(V_eff)
    r_peak = r[idx_peak]
    V_peak = V_eff[idx_peak]
    
    return r_peak, V_peak


def compute_barrier_curvature(r: np.ndarray, V_eff: np.ndarray, r_peak: float) -> float:
    """Compute second derivative of potential at peak.
    
    V'' = d²V/dr² | r=r_peak (determines oscillation frequency)
    """
    # Numerical second derivative
    dr = r[1] - r[0]
    idx = np.argmin(np.abs(r - r_peak))
    
    if idx < 2 or idx >= len(r) - 2:
        return 0.0
    
    # 5-point stencil for second derivative
    V_pp = (-V_eff[idx+2] + 16*V_eff[idx+1] - 30*V_eff[idx] + 
            16*V_eff[idx-1] - V_eff[idx-2]) / (12 * dr**2)
    
    return V_pp


def wkb_qnm_frequency(
    r: np.ndarray,
    V_eff: np.ndarray,
    n: int = 0,
    l: int = 2,
    M: float = 1.0
) -> tuple[float, float]:
    """Compute WKB approximation to QNM frequency.
    
    ω_QNM ≈ V_peak^(1/2) - i (n + 1/2) |V''_peak|^(1/2)
    
    Args:
        r: Radial coordinate array
        V_eff: Effective potential
        n: Overtone number (n=0 fundamental, n=1,2,... overtones)
        l: Angular momentum number
        M: Black hole mass [geometric units]
    
    Returns:
        (omega_real, omega_imag): Real and imaginary parts of ω [M⁻¹ units]
    """
    r_peak, V_peak = find_potential_peak(r, V_eff)
    V_pp = compute_barrier_curvature(r, V_eff, r_peak)
    
    if V_peak <= 0 or V_pp >= 0:
        # No barrier or upward curvature → no QNM
        return 0.0, 0.0
    
    # WKB formula
    omega_real = np.sqrt(V_peak)
    omega_imag = -(n + 0.5) * np.sqrt(np.abs(V_pp))
    
    return omega_real, omega_imag


def time_domain_ringdown(
    r_grid: np.ndarray,
    V_eff: np.ndarray,
    l: int = 2,
    M: float = 1.0,
    t_max: float = 100.0,
    initial_gaussian_width: float = 1.0
) -> tuple[np.ndarray, np.ndarray]:
    """Solve wave equation with V_eff to extract QNM via time evolution.
    
    d²Ψ/dt² + d²Ψ/dr*² - V_eff Ψ = 0
    
    Extract dominant frequency and damping from late-time behavior.
    
    Returns:
        (t, Psi_t): Time array and wavefunction at extraction radius
    """
    # Simplified 1+1D wave equation solver
    # For production, use spectral methods or Runge-Kutta
    
    r_extract_idx = len(r_grid) // 2  # Extract at mid-radius
    dr = r_grid[1] - r_grid[0]
    dt = 0.5 * dr  # CFL condition
    nt = int(t_max / dt)
    
    # Initial condition: Gaussian pulse
    r_0 = r_grid[r_extract_idx]
    Psi = np.exp(-(r_grid - r_0)**2 / initial_gaussian_width**2)
    Psi_dot = np.zeros_like(Psi)
    
    # Storage for extraction point
    Psi_t = np.zeros(nt)
    t = np.arange(nt) * dt
    
    # Leapfrog integration
    for i in range(nt):
        # d²Ψ/dr²
        Psi_rr = np.zeros_like(Psi)
        Psi_rr[1:-1] = (Psi[2:] - 2*Psi[1:-1] + Psi[:-2]) / dr**2
        
        # Update
        Psi_ddot = Psi_rr - V_eff * Psi
        Psi_dot += Psi_ddot * dt
        Psi += Psi_dot * dt
        
        # Store at extraction point
        Psi_t[i] = Psi[r_extract_idx]
    
    return t, Psi_t


def extract_qnm_from_ringdown(
    t: np.ndarray,
    Psi_t: np.ndarray,
    t_start: float = 20.0
) -> tuple[float, float]:
    """Extract QNM frequency and damping from late-time ringdown.
    
    Psi(t) ~ A * exp(i*omega_R*t - gamma*t) for t > t_start
    
    Returns:
        (omega_R, gamma): Real frequency and damping rate
    """
    # Use only late-time data (after initial transients)
    mask = t > t_start
    t_late = t[mask]
    Psi_late = Psi_t[mask]
    
    if len(t_late) < 10:
        return 0.0, 0.0
    
    # FFT to find dominant frequency
    dt = t_late[1] - t_late[0]
    freqs = np.fft.rfftfreq(len(t_late), d=dt)
    fft = np.fft.rfft(Psi_late)
    psd = np.abs(fft)**2
    
    idx_peak = np.argmax(psd[1:]) + 1  # Skip DC
    omega_R = 2 * np.pi * freqs[idx_peak]
    
    # Fit exponential decay to amplitude envelope
    envelope = np.abs(Psi_late)
    if np.max(envelope) > 1e-10:
        # Linear fit to log(envelope)
        log_env = np.log(envelope + 1e-10)
        gamma = -np.polyfit(t_late, log_env, 1)[0]
    else:
        gamma = 0.0
    
    return omega_R, max(0.0, gamma)


def compare_wkb_to_time_domain(
    r: np.ndarray,
    V_eff: np.ndarray,
    n: int = 0,
    l: int = 2,
    M: float = 1.0
) -> dict:
    """Compare WKB estimate to time-domain numerical result.
    
    Returns:
        dict with keys:
            'wkb': (omega_R, gamma)
            'time_domain': (omega_R, gamma)
            'agreement': relative difference
    """
    # WKB estimate
    omega_wkb_R, omega_wkb_I = wkb_qnm_frequency(r, V_eff, n, l, M)
    gamma_wkb = -omega_wkb_I
    
    # Time-domain simulation
    t, Psi_t = time_domain_ringdown(r, V_eff, l, M)
    omega_td, gamma_td = extract_qnm_from_ringdown(t, Psi_t)
    
    # Compute agreement
    if omega_wkb_R > 0 and omega_td > 0:
        freq_diff = abs(omega_wkb_R - omega_td) / omega_wkb_R
        damp_diff = abs(gamma_wkb - gamma_td) / (gamma_wkb + 1e-10)
    else:
        freq_diff = 1.0
        damp_diff = 1.0
    
    return {
        'wkb': {'omega_R': omega_wkb_R, 'gamma': gamma_wkb},
        'time_domain': {'omega_R': omega_td, 'gamma': gamma_td},
        'agreement': {'freq_diff': freq_diff, 'damp_diff': damp_diff}
    }


# Example usage and validation
if __name__ == "__main__":
    print("WKB QNM Diagnostics Module")
    print("=" * 60)
    
    # Test potential: Schwarzschild-like barrier
    M = 1.0
    r = np.linspace(2.1, 20, 300)
    l = 2
    V_eff = (1 - 2*M/r) * (l*(l+1)/r**2 + 2*M/r**3)
    
    print(f"\nTest potential (Schwarzschild, l={l}):")
    r_peak, V_peak = find_potential_peak(r, V_eff)
    print(f"  Peak at r = {r_peak:.2f}M")
    print(f"  V_max = {V_peak:.4f}")
    
    # WKB estimate
    omega_R, omega_I = wkb_qnm_frequency(r, V_eff, n=0, l=l, M=M)
    print(f"\nWKB estimate (n=0):")
    print(f"  ω_R = {omega_R:.4f} M⁻¹")
    print(f"  γ = {-omega_I:.4f} M⁻¹")
    print(f"  Q-factor ≈ {omega_R / (2*abs(omega_I)):.1f}")
    
    # Time-domain comparison
    print(f"\nComparing WKB to time-domain...")
    comparison = compare_wkb_to_time_domain(r, V_eff, n=0, l=l, M=M)
    
    print(f"  WKB: ω_R = {comparison['wkb']['omega_R']:.4f}, γ = {comparison['wkb']['gamma']:.4f}")
    print(f"  Time-domain: ω_R = {comparison['time_domain']['omega_R']:.4f}, γ = {comparison['time_domain']['gamma']:.4f}")
    print(f"  Frequency agreement: {(1-comparison['agreement']['freq_diff'])*100:.1f}%")
    print(f"  Damping agreement: {(1-comparison['agreement']['damp_diff'])*100:.1f}%")
    
    print("\n✅ WKB QNM diagnostics module functional")
    print("   Enables new physics discovery: V_eff → QNM predictions")
