"""
Tests for action principle and field equations.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from field_equations.action import (
    CoherenceGravityParams,
    CoherenceFieldPotential,
    ActionFunctional,
    symbolic_field_equations
)
from field_equations.einstein_coherence import (
    compute_energy_cost_reduction,
    WeakFieldSolver
)


def test_params_creation():
    """Test parameter dataclass creation and validation."""
    params = CoherenceGravityParams()
    
    assert params.G == 6.674e-11
    assert params.c == 2.998e8
    assert params.xi == 1.0
    
    # Test invalid G
    with pytest.raises(ValueError):
        CoherenceGravityParams(G=-1.0)


def test_potentials():
    """Test coherence field potentials."""
    Phi = np.array([0.0, 1.0, 2.0, -1.0])
    m = 1.0
    lambda_ = 0.1
    
    # Quadratic
    V_quad = CoherenceFieldPotential.quadratic(Phi, m)
    assert V_quad[0] == 0.0
    assert V_quad[1] == 0.5
    
    # Quartic
    V_quart = CoherenceFieldPotential.quartic(Phi, m, lambda_)
    assert V_quart[0] == 0.0
    assert V_quart[1] > V_quad[1]  # Quartic adds positive term


def test_effective_g():
    """Test effective gravitational constant calculation."""
    params = CoherenceGravityParams(xi=1.0)
    pot = lambda x: CoherenceFieldPotential.quadratic(x, 1.0)
    pot_deriv = lambda x: CoherenceFieldPotential.quadratic_derivative(x, 1.0)
    
    action = ActionFunctional(params, pot, pot_deriv)
    
    # Test G_eff(Φ=0) = G
    Phi0 = np.array([0.0])
    G_eff = action.effective_gravitational_constant(Phi0)
    assert np.isclose(G_eff[0], params.G)
    
    # Test G_eff(Φ>0) < G (suppression)
    Phi_large = np.array([1e15])
    G_eff_large = action.effective_gravitational_constant(Phi_large)
    assert G_eff_large[0] < params.G
    
    # Test ratio
    ratio = action.effective_coupling_ratio(Phi_large)
    assert ratio[0] < 1.0


def test_energy_cost_reduction():
    """Test energy cost reduction calculation."""
    G = 6.674e-11
    xi = 1.0
    
    # Small coherence: minimal reduction
    # For G_eff/G > 0.9, need 8πGξΦ₀² < 1/9
    # With G ≈ 6.674e-11 and ξ=1: Φ₀ < ~8e3
    Phi0_small = 1e3
    result_small = compute_energy_cost_reduction(Phi0_small, xi, G)
    assert result_small['G_eff_ratio'] < 1.0
    assert result_small['G_eff_ratio'] > 0.9  # Small reduction
    
    # Large coherence: significant reduction
    Phi0_large = 1e20
    result_large = compute_energy_cost_reduction(Phi0_large, xi, G)
    assert result_large['G_eff_ratio'] < result_small['G_eff_ratio']
    
    # Check that required Phi0 for target reduction makes sense
    assert result_large['Phi0_for_1e6_reduction'] > 0


def test_weak_field_solver():
    """Test weak-field linearized solver."""
    params = CoherenceGravityParams(xi=1.0)
    pot = lambda x: CoherenceFieldPotential.quadratic(x, 1.0)
    pot_deriv = lambda x: CoherenceFieldPotential.quadratic_derivative(x, 1.0)
    
    action = ActionFunctional(params, pot, pot_deriv)
    solver = WeakFieldSolver(params, action)
    
    # Test linearized equation with coherence background
    rho = np.array([1e3])  # kg/m³
    Phi0 = 1e15
    phi = np.zeros_like(rho)
    
    result = solver.linearized_einstein_equation(rho, Phi0, phi)
    
    assert 'G_effective' in result
    assert 'coupling_ratio' in result
    assert result['coupling_ratio'] < 1.0  # Suppression


def test_symbolic_equations():
    """Test symbolic field equation derivation."""
    eqs = symbolic_field_equations()
    
    assert 'lagrangian' in eqs
    assert 'coherence_equation' in eqs
    assert 'effective_G' in eqs
    assert 'coupling_ratio' in eqs


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
