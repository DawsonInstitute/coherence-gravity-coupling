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
    # Newtonian limit (xi=0, Phi0=0)
    res = run_geometric_cavendish(xi=0.0, Phi0=0.0, verbose=False)
    tau_sim = abs(res['tau_newtonian'])

    # Assert the torque magnitude sits within broad physically reasonable bounds
    # for kilogram–gram at 0.1–0.2 m separations with this geometry.
    assert 1e-14 <= tau_sim <= 1e-11, (
        f"Newtonian torque magnitude out of expected range: {tau_sim:.3e} N·m")
