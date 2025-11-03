# EFQS next steps (in light of null_results.tex)

This note translates the laboratory nulls into a concrete plan for the Extreme-Field QED Simulator (EFQS) gravitational-coupling plugin and scripts.

## Summary of null-result implications
- Laboratory $R\sim10^{-26}\,\mathrm{m^{-2}}$ implies weak sensitivity to curvature–EM operators; treat our bound $\kappa_R < 5\times10^{17}\,\mathrm{m^2}}$ (at $B=10$\,T) as a baseline prior.
- No observed $\xi$-dependence across materials in present geometries indicates we are at numerical floors for $41^3$; production runs should use $\ge 61^3$ and exploit geometry/boundary amplifiers.

## Immediate engineering tasks
1) API parity for scripts (simulate, sweep, experiments)
- Add in `efqs.gravitational_coupling` the expected functions:
  - `quadrupole_moment(positions, energy_elements_J)` (wrapper over compute_quadrupole)
  - `strain_far_field(Q_t, dt, R, use_tt=True, use_spectral=True)`
  - `radiated_power_from_quadrupole(Q_t, dt, use_spectral=True)`
  - `dominant_frequency(series, dt, component=None)`
  - `stress_energy_from_fields(E, B, include_qed=False)` (return `T00 = (ε0 E^2 + B^2/μ0)/2` with optional QED corrections)

2) ✅ Frequency-domain diagnostics
- Report peak frequency, amplitude at peak, and -3 dB bandwidth for selected strain components.

3) ✅ Boundary-condition amplifier
- Optional TT-projection direction (via `line_of_sight` parameter in `strain_far_field`) implemented
- Robin/impedance boundary-condition hooks: To be implemented as separate module (see tasks 10-15 for Gorkavenko et al. related implementation)

## Validation
- Re-enable `scripts/simulate_gravity_coupling.py` quickstart; verify nonzero but tiny $h_\mathrm{rms}$ and sensible dominant frequency.
- Run `scripts/run_sweep.py` on a short grid of parameters; check monotonic trends and spectra.

## Implementation status
EFQS repository confirmed present in workspace at `../extreme-field-qed-simulator`. API wrappers and validation tasks ready to execute.

## Stretch goals
- Add a simple TT-projection selector by line-of-sight vector and unit tests for projection identities.
- Provide an astrophysical recast notebook that propagates our lab $\kappa_R$ prior to compact-object QNM parameter bounds (see Karimabadi et al.).
