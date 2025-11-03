# Hell & Lüst (2025)

Citation: Hell & Lüst, arXiv:2509.20217 (2025)

**Title:** Aspects of non-minimally coupled curvature with power laws

## Summary
This work systematically classifies theories with power-law curvature terms ($R^\ell$, $\sigma^n R^m$ where $\sigma$ is a scalar) including non-minimal couplings. The authors analyze degrees of freedom (DOF) structure on cosmological (FLRW) backgrounds by mapping to Einstein frame and studying singular points in the Jordan frame. Key result: the number and nature of propagating scalar modes depend sensitively on model parameters $(\ell, m, n)$ and background values of $R$ and $\sigma$; some cases exhibit vanishing or strongly coupled modes at leading order near special points (e.g., $R=0$, $m/\ell = 0,1$).

## Key findings relevant to this project
- **Parameter-dependent DOF:** For power-law models, the scalar mode count and their kinetic structure change across parameter space; some models have no propagating scalar at $R=0$ (flat background), others activate modes only away from singular loci.
- **Frame dependence:** Singular points in Jordan frame (e.g., $R=0$ for certain $\ell,m$) can be regular in Einstein frame or vice versa; constraints derived in one frame don't always translate naively.
- **Small-$R$ regimes:** Laboratory backgrounds ($R \sim 10^{-26}\,\mathrm{m^{-2}}$) sit near the $R \to 0$ limit for many power-law models, where extra DOF may decouple or become non-dynamical.

## Relation to null_results.tex
Our null_results.tex observes **no detectable signal** from $\xi$ variation ($\xi \in \{50,100\}$) or materials at $R \sim 10^{-26}\,\mathrm{m^{-2}}$. Hell & Lüst's DOF analysis provides a theoretical context:

- **Consistency with decoupling:** If the effective scalar mode in a non-minimally coupled theory (like our $\xi R\Phi^2$ or extended power-law models) has vanishing kinetic term or becomes non-propagating at small $R$, laboratory nulls are expected.
- **Parameter priors:** Our constraints help identify which $(\ell, m, n, \xi)$ combinations remain viable: models predicting strong scalar-mode activity at $R \approx 0$ are disfavored unless they evade our torque/curvature observable.
- **Singular-point proximity:** Terrestrial $R$ is extremely close to the $R=0$ singular locus; Hell & Lüst show this can suppress or eliminate leading-order scalar dynamics in certain power-law classes.

## How these nulls inform future research
- **DOF-guided parameter scans:** Prioritize $(\ell, m, \xi)$ regions where DOF analysis predicts non-trivial but weakly coupled scalars at small $R$; use our nulls to calibrate expected signal floors.
- **Frame-invariant observables:** Translate our $\kappa_R$ bounds and $\xi$ insensitivity into Einstein-frame statements; check whether frame choice affects interpretation of "null" vs. "excluded."
- **Dynamical transitions:** When extending to time-dependent $\Phi(t)$ or varying $R(t)$ (e.g., cosmological simulations), monitor for crossing singular points where DOF structure changes; implement diagnostics (condition number, residual stability) to detect mode activation/deactivation.

## Actionable follow-ups
- [⏳] Tabulate DOF predictions from Hell & Lüst for standard power-law extensions of our $\xi R\Phi^2$ model; map to expected torque signal scalings.
- [⏳] Implement numerical health checks: flag runs approaching singular points (e.g., $R \to 0$, parameter-dependent loci); output warnings if DOF structure is ambiguous.
- [⏳] Document our Jordan vs. Einstein frame conventions; cross-reference Hell & Lüst's Table 1–3 classifications to ensure consistency in reporting constraints.
- [⏳] Extend convergence study to include $R$-dependent mesh refinement near singular points; verify stability of nulls under resolution changes.

## Notes for EFQS integration
- [⏳] Add a "DOF mode" selector: users specify $(\ell, m, n)$ for power-law curvature terms; EFQS warns if chosen parameters place simulation near singular point or imply decoupled scalars.
- [⏳] Cross-validate: reproduce Hell & Lüst's flat-FLRW mode counts numerically for toy models; ensure our stress–energy → quadrupole pipeline respects frame transformations.
