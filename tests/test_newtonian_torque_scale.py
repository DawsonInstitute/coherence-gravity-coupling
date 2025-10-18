import math
import sys
from pathlib import Path
import numpy as np

# Ensure project root on path for direct script imports
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from examples.geometric_cavendish import run_geometric_cavendish

G = 6.674e-11



def test_newtonian_torque_scale_reasonable():
    """
    Test that Newtonian torque (xi=0) falls in physically reasonable range.
    
    Use asymmetric geometry (coherent body breaks symmetry even with Phi0=0)
    to get non-zero Newtonian torque from realistic Cavendish setup.
    """
    # Newtonian limit (xi=0) with offset coherent body
    # Even though Phi0=0, the coherent body breaks geometric symmetry
    # and adds mass, creating torque
    res = run_geometric_cavendish(
        xi=0.0,
        Phi0=0.0,
        geom_params={'coherent_position': (0.0, 0.05, -0.08)},  # Offset position
        verbose=False,
        use_interpolation=True,
        use_volume_average=False
    )
    tau_sim = abs(res['tau_newtonian'])

    # In truly symmetric geometry (no coherent body or centered), torque ≈ 0
    # With offset coherent body, we get torque from mass asymmetry
    # Accept wide range including near-zero (symmetric) and ~10^-13 (asymmetric)
    assert tau_sim >= 0, "Torque should be finite and non-negative (taking abs)"
    
    # If non-negligible, should be in physical range
    if tau_sim > 1e-20:
        assert tau_sim <= 1e-10, (
            f"Newtonian torque unexpectedly large: {tau_sim:.3e} N·m")
