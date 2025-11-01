"""
Tests for Curvature Coupling Extension Module

Validates:
1. EM invariant calculations
2. Ricci-EM coupling energy
3. Modified Maxwell equations
4. Stress-energy corrections
5. Effective G calculations
6. Exclusion limit estimates

Author: GitHub Copilot (Claude Sonnet 4.5)
Date: October 2025
License: MIT
"""

import pytest
import numpy as np
from pathlib import Path
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.field_equations.curvature_coupling import (
    CurvatureCouplingParams,
    CurvatureCouplingCalculator,
    compute_exclusion_limits,
    C_LIGHT,
    EPSILON_0,
    G_NEWTON
)


class TestEMInvariants:
    """Test electromagnetic field invariant calculations."""
    
    def test_zero_fields(self):
        """Test that zero fields give zero invariants."""
        params = CurvatureCouplingParams()
        calc = CurvatureCouplingCalculator(params)
        
        E = np.zeros(3)
        B = np.zeros(3)
        
        inv = calc.electromagnetic_invariants(E, B)
        
        assert inv['F_squared'] == 0.0
        assert inv['F_dual_squared'] == 0.0
        assert inv['E_squared'] == 0.0
        assert inv['B_squared'] == 0.0
    
    def test_pure_electric_field(self):
        """Test invariants for pure electric field."""
        params = CurvatureCouplingParams()
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1.0, 0.0, 0.0])  # 1 V/m
        B = np.zeros(3)
        
        inv = calc.electromagnetic_invariants(E, B)
        
        # F² = 2(B² - E²/c²) = -2E²/c²
        expected_F_sq = -2 * 1.0 / C_LIGHT**2
        
        assert np.isclose(inv['F_squared'], expected_F_sq)
        assert inv['F_dual_squared'] == 0.0  # No E·B
        assert inv['E_squared'] == 1.0
        assert inv['B_squared'] == 0.0
    
    def test_pure_magnetic_field(self):
        """Test invariants for pure magnetic field."""
        params = CurvatureCouplingParams()
        calc = CurvatureCouplingCalculator(params)
        
        E = np.zeros(3)
        B = np.array([0.0, 1.0, 0.0])  # 1 T
        
        inv = calc.electromagnetic_invariants(E, B)
        
        # F² = 2B²
        expected_F_sq = 2 * 1.0
        
        assert np.isclose(inv['F_squared'], expected_F_sq)
        assert inv['F_dual_squared'] == 0.0
        assert inv['E_squared'] == 0.0
        assert inv['B_squared'] == 1.0
    
    def test_parallel_fields(self):
        """Test invariants for parallel E and B."""
        params = CurvatureCouplingParams()
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1.0, 0.0, 0.0])
        B = np.array([1.0, 0.0, 0.0])
        
        inv = calc.electromagnetic_invariants(E, B)
        
        # Dual invariant: *F² = 4E·B/c
        expected_dual = 4 * 1.0 / C_LIGHT
        
        assert np.isclose(inv['F_dual_squared'], expected_dual)
        assert inv['E_dot_B'] == 1.0


class TestRicciEMCoupling:
    """Test Ricci-EM coupling calculations."""
    
    def test_disabled_coupling_gives_zero(self):
        """Test that disabled coupling returns zero energy."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e-10,
            enable_ricci_em=False
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1e6, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        
        energy = calc.ricci_em_coupling_energy(R, E, B)
        assert energy == 0.0
    
    def test_coupling_energy_scaling(self):
        """Test that coupling energy scales correctly with parameters."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e-15,
            enable_ricci_em=True
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1e6, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        
        energy1 = calc.ricci_em_coupling_energy(R, E, B)
        
        # Double the Ricci scalar
        energy2 = calc.ricci_em_coupling_energy(2*R, E, B)
        
        # Energy should double
        assert np.isclose(energy2, 2*energy1)
    
    def test_coupling_energy_sign(self):
        """Test energy sign consistency."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e-15,
            enable_ricci_em=True
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.zeros(3)
        B = np.array([0, 1.0, 0])  # Pure magnetic -> F² > 0
        R = 1e-20  # Positive curvature
        
        energy = calc.ricci_em_coupling_energy(R, E, B)
        
        # With positive κ, R, and F², energy should be positive
        assert energy > 0


class TestMaxwellCorrections:
    """Test modified Maxwell equations."""
    
    def test_no_correction_when_disabled(self):
        """Test that corrections are zero when coupling is disabled."""
        params = CurvatureCouplingParams(enable_ricci_em=False)
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1.0, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        grad_R = np.array([1e-22, 0, 0])
        
        corr = calc.modified_maxwell_equations(E, B, R, grad_R)
        
        assert np.allclose(corr['faraday_correction'], 0)
        assert np.allclose(corr['ampere_correction'], 0)
    
    def test_no_correction_without_gradient(self):
        """Test that corrections are zero without gradient."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e-15,
            enable_ricci_em=True
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1.0, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        
        corr = calc.modified_maxwell_equations(E, B, R, grad_ricci=None)
        
        assert np.allclose(corr['faraday_correction'], 0)
        assert np.allclose(corr['ampere_correction'], 0)
    
    def test_faraday_correction_perpendicular(self):
        """Test Faraday correction is perpendicular to both ∇R and B."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e-15,
            enable_ricci_em=True
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1.0, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        grad_R = np.array([1e-22, 0, 0])
        
        corr = calc.modified_maxwell_equations(E, B, R, grad_R)
        faraday_corr = corr['faraday_correction']
        
        # Should be perpendicular to both ∇R and B
        assert np.isclose(np.dot(faraday_corr, grad_R), 0, atol=1e-30)
        assert np.isclose(np.dot(faraday_corr, B), 0, atol=1e-30)


class TestStressEnergyCorrection:
    """Test stress-energy tensor corrections."""
    
    def test_correction_shape(self):
        """Test that correction has correct shape."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e-15,
            enable_ricci_em=True
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1e6, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        R_mu_nu = np.diag([R/4, R/4, R/4, R/4])
        
        delta_T = calc.stress_energy_correction(R, R_mu_nu, E, B)
        
        assert delta_T.shape == (4, 4)
    
    def test_correction_scales_with_kappa(self):
        """Test that correction scales with coupling constant."""
        params1 = CurvatureCouplingParams(
            kappa_ricci_em=1e-15,
            enable_ricci_em=True
        )
        calc1 = CurvatureCouplingCalculator(params1)
        
        params2 = CurvatureCouplingParams(
            kappa_ricci_em=2e-15,
            enable_ricci_em=True
        )
        calc2 = CurvatureCouplingCalculator(params2)
        
        E = np.array([1e6, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        R_mu_nu = np.diag([R/4, R/4, R/4, R/4])
        
        delta_T1 = calc1.stress_energy_correction(R, R_mu_nu, E, B)
        delta_T2 = calc2.stress_energy_correction(R, R_mu_nu, E, B)
        
        # Should scale linearly with κ
        assert np.allclose(delta_T2, 2*delta_T1)


class TestEffectiveG:
    """Test effective gravitational constant."""
    
    def test_no_modification_when_disabled(self):
        """Test G_eff = G when coupling is disabled."""
        params = CurvatureCouplingParams(enable_ricci_em=False)
        calc = CurvatureCouplingCalculator(params)
        
        E = np.array([1e6, 0, 0])
        B = np.array([0, 1.0, 0])
        R = 1e-20
        
        G_eff = calc.effective_gravitational_constant(E, B, R)
        
        assert G_eff == G_NEWTON
    
    def test_modification_increases_with_fields(self):
        """Test that G_eff increases with stronger fields."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e-5,  # Very large κ to see effect
            enable_ricci_em=True
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.zeros(3)
        B1 = np.array([0, 10.0, 0])  # 10 T
        B2 = np.array([0, 20.0, 0])  # 20 T
        R = 1e-10  # Strong curvature
        
        G_eff1 = calc.effective_gravitational_constant(E, B1, R)
        G_eff2 = calc.effective_gravitational_constant(E, B2, R)
        
        # Stronger field should give larger modification
        # Since modification = 1 + κRF², and F² ∝ B², we expect:
        # (G_eff2 - G) / (G_eff1 - G) ≈ (B2/B1)² = 4
        delta1 = G_eff1 - G_NEWTON
        delta2 = G_eff2 - G_NEWTON
        
        assert delta2 > delta1  # Basic monotonicity
        assert np.isclose(delta2 / delta1, 4.0, rtol=0.1)  # Quadratic scaling
    
    def test_regularization_prevents_negative_g(self):
        """Test that regularization prevents unphysical G_eff."""
        params = CurvatureCouplingParams(
            kappa_ricci_em=1e10,  # Absurdly large κ
            enable_ricci_em=True
        )
        calc = CurvatureCouplingCalculator(params)
        
        E = np.zeros(3)
        B = np.array([0, 100.0, 0])  # Very strong field
        R = 1e-10  # Strong curvature
        
        G_eff = calc.effective_gravitational_constant(E, B, R)
        
        # Should be regularized to reasonable range
        assert G_eff > 0
        assert G_eff < 100 * G_NEWTON


class TestExclusionLimits:
    """Test exclusion limit calculations."""
    
    def test_limit_scales_inversely(self):
        """Test that limit scales inversely with precision."""
        R = 1e-20
        F_sq = 1.0
        
        limit1 = compute_exclusion_limits(1e-6, R, F_sq)
        limit2 = compute_exclusion_limits(1e-7, R, F_sq)
        
        # Tighter precision gives tighter (smaller) limit
        assert limit2['kappa_ricci_em_limit'] < limit1['kappa_ricci_em_limit']
    
    def test_limit_with_zero_inputs(self):
        """Test that zero inputs return infinite limit."""
        limit = compute_exclusion_limits(1e-6, 0.0, 1.0)
        assert limit['kappa_ricci_em_limit'] == np.inf
        
        limit = compute_exclusion_limits(1e-6, 1e-20, 0.0)
        assert limit['kappa_ricci_em_limit'] == np.inf
    
    def test_limit_order_of_magnitude(self):
        """Test that exclusion limits have reasonable magnitude."""
        # Typical experimental setup
        precision = 1e-6
        R = 1e-26  # Earth curvature scale
        F_sq = 1.0  # 1 T field
        
        limit = compute_exclusion_limits(precision, R, F_sq)
        
        # Limit = precision / (R × F²) = 1e-6 / (1e-26 × 1) = 1e20 m²
        # This is actually a very weak constraint (very large allowed κ)
        assert limit['kappa_ricci_em_limit'] > 1e15


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
