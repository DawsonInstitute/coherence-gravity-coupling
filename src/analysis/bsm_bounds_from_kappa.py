"""
BSM coupling bounds from curvature–EM coupling κ_R.

Scope:
- Provide a conservative, unit-consistent mapping from κ_R (R·F^2 operator) to
  an "equivalent" dark photon kinetic mixing ε in a background with curvature scale R.
- For axions, provide only a model-dependent parametric mapping with explicit portal
  assumptions; do not overinterpret CP-even ↔ CP-odd operator equivalences.

Conventions:
- Natural units: ℏ = c = 1 for intermediate computations.
- κ_R has units of length^2 in existing code (m^2). In natural units, [κ_R] = GeV^-2.
- Ricci scalar or characteristic curvature scale R has units of length^-2; in natural units, [R] = GeV^2.

Key mapping (robust):
    ε_eff ≈ C_ε · (κ_R · R)
This is dimensionless and captures the size of Maxwell equation modifications in curved space
relative to standard kinetic terms, up to an O(1) matching coefficient C_ε.

Axion note:
    The axion operator is CP-odd: (g_{aγγ}/4) a F·F~ (dimension-5), while κ_R couples to
    CP-even F^2. Any mapping requires additional CP-odd portal assumptions (e.g., a·R or
    curvature-induced axion background). We expose a parametric helper that includes explicit
    model coefficients and a reference mass scale Λ to keep dimensions consistent.

"""
from __future__ import annotations

from typing import Dict, Iterable, Tuple
import math


# Unit conversions
M_TO_GEV_INV = 5.067730716156395e15       # 1 m ≈ 5.07e15 GeV^-1
M2_TO_GEV2_INV = M_TO_GEV_INV ** 2        # 1 m^-2 ≈ (5.07e15)^2 GeV^2
M2_TO_GEV2 = M_TO_GEV_INV ** 2            # same factor for converting m^-2 -> GeV^2
M2_TO_GEVM2 = (1.0 / M_TO_GEV_INV) ** 2   # 1 m^2 ≈ (1/5.07e15)^2 GeV^-2


class CurvatureEnvironment:
    def __init__(self, name: str, R_m2: float, note: str = ""):
        self.name = name
        self.R_m2 = R_m2
        self.note = note
    @property
    def R_GeV2(self) -> float:
        return self.R_m2 * M2_TO_GEV2


DEFAULT_ENVIRONMENTS: Dict[str, CurvatureEnvironment] = {
    # These are order-of-magnitude benchmarks; refine as needed with specific geometries.
    # Lab table-top (near flat): effectively zero; use a tiny regulator to show scaling
    "lab_flat": CurvatureEnvironment("lab_flat", R_m2=1e-30, note="Near-flat lab background"),
    # Earth surface tidal scale (very rough): ~ 1e-27 to 1e-26 m^-2 (order-of-magnitude)
    "earth_surface": CurvatureEnvironment("earth_surface", R_m2=1e-26, note="Order-of-magnitude"),
    # Low Earth Orbit: similar order
    "leo": CurvatureEnvironment("leo", R_m2=5e-27, note="Order-of-magnitude"),
    # Magnetar surface (as used in prior code): ~1e-6 m^-2
    "magnetar_surface": CurvatureEnvironment("magnetar_surface", R_m2=1e-6, note="Strong curvature benchmark"),
}


def kappa_m2_to_GeV2_inv(kappa_m2: float) -> float:
    """Convert κ_R in m^2 to GeV^-2."""
    return kappa_m2 * M2_TO_GEVM2


def curvature_m2_to_GeV2(R_m2: float) -> float:
    """Convert curvature scale R in m^-2 to GeV^2."""
    return R_m2 * M2_TO_GEV2


def epsilon_equiv(kappa_m2: float, R_m2: float, C_eps: float = 1.0) -> float:
    """Equivalent dark photon kinetic mixing from κ_R at curvature scale R.

    ε_eff ≈ C_ε · (κ_R · R)

    Args:
        kappa_m2: κ_R in m^2 (as used elsewhere in repository)
        R_m2: curvature scale in m^-2
        C_eps: O(1) matching coefficient capturing model dependence

    Returns:
        epsilon_eff (dimensionless)
    """
    return C_eps * (kappa_m2 * R_m2)


def epsilon_sweep(
    kappa_values_m2: Iterable[float],
    envs: Dict[str, CurvatureEnvironment] = DEFAULT_ENVIRONMENTS,
    C_eps_values: Iterable[float] = (1.0, 1e-2, 1e-4)
) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Compute ε_eff across κ_R values, environments, and C_ε choices.

    Returns nested dict: env -> Cε -> {kappa_str: epsilon}
    """
    out: Dict[str, Dict[str, Dict[str, float]]] = {}
    for env_name, env in envs.items():
        out_env: Dict[str, Dict[str, float]] = {}
        for C in C_eps_values:
            table: Dict[str, float] = {}
            for kappa in kappa_values_m2:
                e = epsilon_equiv(kappa, env.R_m2, C)
                table[f"{kappa:.3e} m^2"] = e
            out_env[f"C_eps={C:.0e}"] = table
        out[env_name] = out_env
    return out


def axion_equiv_parametric(
    kappa_m2: float,
    R_m2: float,
    C_a: float = 1.0,
    Lambda_GeV: float = 1.0e4,
) -> float:
    """Parametric axion-equivalent coupling g_{aγγ}^{equiv} [GeV^-1].

    Warning: Model-dependent! We present a portal-inspired scaling:
        g_equiv ≈ C_a · (κ_R · R) / Λ
    where Λ is a UV suppression scale (GeV). This keeps dimensions consistent
    and vanishes in flat space (R→0). The CP structure mismatch remains; this is
    only a benchmarking device to compare orders of magnitude with axion limits.

    Args:
        kappa_m2: κ_R in m^2
        R_m2: curvature scale in m^-2
        C_a: dimensionless O(1) matching coefficient
        Lambda_GeV: UV scale suppressing the CP-odd portal (GeV)

    Returns:
        g_equiv [GeV^-1]
    """
    # Dimensionless core factor
    dimless = kappa_m2 * R_m2
    return (C_a * dimless) / max(Lambda_GeV, 1e-30)


def pretty_si(value: float) -> str:
    """Pretty SI string for dimensionless small numbers."""
    if value == 0:
        return "0"
    exp = int(math.floor(math.log10(abs(value))))
    mant = value / (10 ** exp)
    return f"{mant:.2f}e{exp:+d}"


def demo_printout() -> None:
    """Small demo to print ε_equiv and g_equiv tables for common κ_R bounds."""
    kappa_grid = [1e-19, 1e-11, 1e+0, 1e+3]  # m^2 (illustrative spread)
    print("[Dark photon: ε_equiv ≈ Cε (κ_R·R)]")
    eps_tables = epsilon_sweep(kappa_grid)
    for env_name, env_tab in eps_tables.items():
        print(f"\nEnvironment: {env_name} (R ≈ {envs_to_si(env_name)})")
        for C_key, row in env_tab.items():
            print(f"  {C_key}:")
            for kappa_str, eps in row.items():
                print(f"    κ_R={kappa_str}: ε ≈ {pretty_si(eps)}")

    print("\n[Axion (illustrative, model-dependent): g_equiv ≈ Ca (κ_R·R)/Λ]")
    for env_name, env in DEFAULT_ENVIRONMENTS.items():
        g = axion_equiv_parametric(1e-11, env.R_m2, C_a=1.0, Lambda_GeV=1e4)
        print(f"  Env {env_name}: κ_R=1e-11 m^2 → g_equiv ≈ {g:.2e} GeV^-1 (Λ=10 TeV)")


def envs_to_si(env_name: str) -> str:
    env = DEFAULT_ENVIRONMENTS[env_name]
    return f"{env.R_m2:.1e} m^-2"


if __name__ == "__main__":
    demo_printout()
