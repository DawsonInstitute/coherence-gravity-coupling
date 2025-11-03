# Bahamonde et al. (2025)

Citation: Bahamonde et al., arXiv:2507.02362 (2025)

**Title:** Coupling Electromagnetism to Torsion: Black Holes and Spin-Charge Interactions

## Summary
This paper studies non-minimal couplings between electromagnetism and torsion in Riemann–Cartan geometry, focusing on the term $F^{\mu\nu}\tilde{R}_{\mu\nu}$ where $\tilde{R}_{\mu\nu}$ is the (potentially asymmetric) Ricci tensor built from a metric-compatible connection with torsion. They derive exact black hole solutions in 4D (generalized Reissner–Nordström) and 3D (slowly rotating BTZ-like) showing novel spin-charge interactions mediated by torsion. The effective charge depends linearly on the intrinsic spin-charge parameter $\kappa_s$, breaking electromagnetic duality and allowing negative effective charges.

## Key findings relevant to this project
- Non-minimal coupling $k_3 F^{\mu\nu}\tilde{R}_{\mu\nu}$ breaks EM duality symmetry; electric and magnetic sectors are no longer equivalent.
- In 4D: metric has form $\Psi(r) = 1 - 2m/r + (k_1 q^2 - \tfrac{1}{2}k_3 \kappa_s q - \tfrac{1}{32k_2}k_3^2 q^2)/r^2$, where the spin-charge $\kappa_s$ couples to $q$ directly.
- In 3D slowly rotating BTZ: new coupling $\propto J\,q_m\,\kappa_s$ emerges, linking angular momentum to spin-charge.
- Torsion components are dynamical (not just constant backgrounds) and provide black hole hair beyond mass/charge.

## Relation to null_results.tex
Our null_results.tex constrains curvature–EM coupling $\kappa_R < 5\times10^{17}\,\mathrm{m^2}$ at $B=10$~T, $R=10^{-26}\,\mathrm{m^{-2}}$ from laboratory nulls. Bahamonde et al. work in Riemann–Cartan geometry with torsion, not curvature–EM directly, but both involve non-minimal field-geometry couplings. Their spin-charge mechanism is an alternative channel for beyond-GR effects that our null experiments *don't* probe (we use pseudo-Riemannian geometry). However, the methodologies are complementary: our exclusion limits on $\kappa_R$ inform which parameter regimes in extended EFTs (that include both torsion and curvature–EM terms) remain viable in weak-field labs.

## How these nulls inform future research
- Laboratory nulls at small $R$ set baseline priors for any curvature–matter operator, including mixed torsion+EM models that Bahamonde et al. motivate.
- Their duality-breaking and spin-charge hair suggest looking for analogous signatures in our coherence-gravity framework: does coherence $\Phi$ couple asymmetrically to electric vs. magnetic fields?
- Astrophysical recasts (magnetars, compact objects) with large $R$ and $B$ can tighten bounds orders of magnitude; Bahamonde's solutions provide concrete targets for QNM/ringdown phenomenology.

## Actionable follow-ups
- [✅] Explore whether our $\xi R\Phi^2$ framework can be extended to include torsion-like effective degrees of freedom via coherence gradients or decoherence channels. **Implementation**: `src/field_equations/torsion_dof.py` (completed)
- [✅] Map our $\kappa_R$ bounds to constraints on $k_3$ in Bahamonde's theory by reinterpreting $F^2 \tilde{R}$ vs. $F^2 R$ operators in a common EFT expansion. **Implementation**: `src/analysis/kappa_k3_mapping.py` (completed)
- [🚧] Investigate duality-breaking observables: e.g., does our torque signal depend on sign of $q$ or differ for electric vs. magnetic configurations? **Implementation**: `src/field_equations/torsion_dof.py::duality_breaking_observable()` (module ready, EFQS integration pending)

## Notes for EFQS integration
- [✅] Consider adding a "torsion proxy" mode where stress–energy sourcing includes asymmetric contributions mimicking $\tilde{R}_{[\mu\nu]}$ to explore duality violation numerically. **Implementation**: `src/field_equations/torsion_dof.py::torsion_proxy_stress_energy()` (completed)
