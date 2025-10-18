"""
Tests for the refactored geometry parameterization.
"""

import pytest
import numpy as np
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from examples.geometric_cavendish import CavendishGeometry, run_geometric_cavendish

def test_parameter_override():
    """
    Test that passing geom_params to run_geometric_cavendish correctly
    overrides the default geometry.
    """
    # Non-default parameters
    geom_params = {
        'm_test': 0.025,
        'M_source': 5.0,
        'L_torsion': 0.15
    }
    
    result = run_geometric_cavendish(
        xi=10.0,
        Phi0=1e7,
        geom_params=geom_params,
        verbose=False
    )
    
    # Check that the returned geometry in the results matches the override
    result_geom = result['geom_params']
    assert abs(result_geom['m_test'] - 0.025) < 1e-9
    assert abs(result_geom['M_source'] - 5.0) < 1e-9
    assert abs(result_geom['L_torsion'] - 0.15) < 1e-9

def test_near_linear_scaling_small_perturbation():
    """
    Assert that for a small change in mass, the torque scales almost linearly.
    This validates the old assumption for a limited regime and serves as a
    good sanity check on the full simulation.
    """
    
    # Baseline
    base_params = {'m_test': 0.010}
    res_base = run_geometric_cavendish(
        xi=10.0, Phi0=1e7, geom_params=base_params, verbose=False
    )
    tau_base = res_base['tau_coherent']

    # Perturbed (1% mass increase)
    perturbed_params = {'m_test': 0.010 * 1.01}
    res_perturbed = run_geometric_cavendish(
        xi=10.0, Phi0=1e7, geom_params=perturbed_params, verbose=False
    )
    tau_perturbed = res_perturbed['tau_coherent']
    
    # Linear scaling prediction
    tau_predicted = tau_base * 1.01

    # Check that the actual result is very close to the linear prediction
    relative_diff = abs(tau_perturbed - tau_predicted) / abs(tau_predicted)
    
    print(f"\n--- Small Perturbation Test ---")
    print(f"Base torque: {tau_base:.6e}")
    print(f"Perturbed torque (actual): {tau_perturbed:.6e}")
    print(f"Perturbed torque (predicted by linear scaling): {tau_predicted:.6e}")
    print(f"Relative difference: {relative_diff:.3%}")

    assert relative_diff < 0.005, \
        f"Torque should scale near-linearly for small mass changes, but difference was {relative_diff:.3%}"

def test_sweep_functions_run():
    """
    Test that the refactored sweep functions execute without error.
    This is a smoke test.
    """
    from examples.geometric_cavendish import sweep_test_mass, sweep_source_mass

    base_geom = {'coherent_position': (0.0, 0.0, -0.08)}
    
    # Test mass sweep
    m_sweep_result = sweep_test_mass(
        m_test_range=np.array([0.01, 0.015]),
        base_geom_params=base_geom,
        verbose=False
    )
    assert len(m_sweep_result['sweep_results']) == 2
    assert 'optimal' in m_sweep_result

    # Source mass sweep
    s_sweep_result = sweep_source_mass(
        M_source_range=np.array([1.0, 1.5]),
        base_geom_params=base_geom,
        verbose=False
    )
    assert len(s_sweep_result['sweep_results']) == 2
    assert 'optimal' in s_sweep_result

if __name__ == '__main__':
    pytest.main([__file__, '-v'])
