# Gorkavenko et al. (2025)

Citation: Gorkavenko et al., arXiv:2509.06815 (2025)

**Title:** Impact of space-time curvature coupling on the vacuum energy induced by a magnetic topological defect in flat space-time of arbitrary dimension

## Summary
This paper investigates vacuum polarization of a quantized charged massive scalar field with non-minimal curvature coupling $\xi$ near a magnetic topological defect (modeled as an impenetrable tube with magnetic flux) in flat spacetime of arbitrary dimension. The key finding: under generalized Robin boundary conditions at the tube surface, the total induced vacuum energy depends on $\xi$ even in flat spacetime, unlike the special Dirichlet/Neumann cases where it is $\xi$-independent. The energy splits as $E = E_{\mathrm{can}} + E_\xi$ with $E_\xi \propto (1/4 - \xi)$. For small $mr_0$ (thin tubes) and higher dimensions, $E_\xi$ can dominate $E_{\mathrm{can}}$.

## Key findings relevant to this project
- **Flat-space $\xi$ sensitivity via boundaries:** Robin BCs reintroduce $\xi$-dependence in the total energy, providing an independent probe of the curvature coupling constant without requiring macroscopic curvature.
- **Scaling:** $E_\xi$ grows with dimension and for thin tubes; explicit formulas given for $(d+1)$-dimensional spacetime with fractional flux parameter $F = \Phi/\Phi_0$.
- **Physical mechanism:** Vacuum energy shifts arise from Aharonov–Bohm phase + Robin impedance matching; the canonical and $\xi$-dependent contributions are renormalized and separately calculable.

## Relation to null_results.tex
Our null_results.tex finds **no $\xi$ dependence** in torque measurements across materials and $\xi \in \{50,100\}$ at terrestrial $R \sim 10^{-26}\,\mathrm{m^{-2}}$, consistent with weak curvature suppressing direct $\xi R\Phi^2$ effects. Gorkavenko et al. show that **boundary engineering** can bypass this suppression: even at $R\approx0$, Robin conditions on cavity surfaces couple to $\xi$ via vacuum polarization. This suggests:

- Our nulls reflect geometry choice (free-space torque balance with Dirichlet-like or simple boundaries), not a fundamental absence of $\xi$ sensitivity.
- **Actionable insight:** Future experiments should incorporate tunable impedance boundaries (Robin parameter $\theta$) or thin-structure resonators ($mr_0 \ll 1$) to amplify $E_\xi$ and access the $(1/4 - \xi)$ term experimentally.

## How these nulls inform future research
- **Boundary-amplified $\xi$ detection:** Design tabletop experiments with microwave/optical cavities, tunable coatings (variable $\theta$), and magnetic flux inserts to measure $E_\xi$; bypass need for large $R$.
- **Numerical validation:** Extend our PDE solver to support Robin BCs; compute $E_\xi(\xi, mr_0, \theta, F)$ and compare to Gorkavenko's analytic/numerical results in $(3+1)$D.
- **Dimensional scaling:** Our $(3+1)$D nulls are consistent; explore whether lower-dimensional analogs (graphene-like $(2+1)$D systems with disclinations) enhance sensitivity as Gorkavenko predicts.

## Actionable follow-ups
- [🚧] Implement Robin BC solver module with parameter $\theta \in [-\pi/2, \pi/2]$; validate against known Dirichlet ($\theta=0$) and Neumann ($\theta=-\pi/2$) limits. **Implementation**: `src/solvers/robin_bc_poisson.py`
- [⏳] Parametric sweep: vary $mr_0 \in [10^{-3}, 1]$, $\theta$, $F$, $\xi \in [0, 1/2]$; produce contour plots of $|E_\xi|/|E_{\mathrm{can}}|$ to identify high-sensitivity regimes.
- [⏳] Materials feasibility: identify coatings/surfaces realizing target $\theta$ values in cryogenic torsion balance or superconducting cavity geometries.
- [⏳] Cross-check: reproduce Gorkavenko's Table 1 ($(3+1)$D, $mr_0=0.01$, $F=0.25$) numerically and extend to parameter ranges relevant to our experimental proposals.

## Notes for EFQS integration
- [🚧] Add boundary-condition hooks: EFQS stress–energy sourcing should support Robin-impedance surfaces when computing quadrupole moments. **Integration point**: `src/solvers/robin_bc_poisson.py` ready for EFQS coupling
- [⏳] Test: verify $\xi$-independence recovers for $\theta = 0, -\pi/2$ (Dirichlet/Neumann).
