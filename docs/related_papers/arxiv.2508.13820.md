# Karimabadi, Mahdavian Yekta & Alavi (2025)

Citation: Karimabadi, Mahdavian Yekta & Alavi, arXiv:2508.13820 (2025)

**Title:** Effects of non-minimal scalar field couplings with curvature tensors on perturbations in non-commutative Schwarzschild spacetimes

## Summary
This paper compares quasi-normal modes (QNMs) of scalar field perturbations in non-commutative Schwarzschild (NC-Sch) black holes under two distinct non-minimal curvature couplings: (i) **Scalar model** — field couples to Ricci scalar $R$, yielding modified Klein–Gordon with curvature-dependent mass; (ii) **Tensor model** — field derivatives couple to Einstein tensor $G^{\mu\nu}$. Using WKB-6 approximation and time-domain integration, they find QNM spectra agree at low overtones but diverge at high $n$; instabilities appear for large NC parameter $\theta$ or coupling $\zeta$. Near-horizon analytic QNMs confirm consistency and support area quantization.

## Key findings relevant to this project
- **Coupling sensitivity in strong fields:** QNM frequencies (real and imaginary parts) depend measurably on coupling type and NC parameter; effects are pronounced at high overtones and large $\theta$.
- **Time-domain ringdown:** Decay rates controlled by near-horizon physics; oscillation frequencies reflect coupling nature. Both models show damped oscillations confirming stability for small $\theta$, instabilities for large $\theta/\zeta$.
- **Near-horizon universality:** Analytic QNMs in Schwarzschild limit ($\theta \to 0$) recover standard results; deviations scale with $\theta$ and coupling strength.
- **Area quantization:** Consistent with Bekenstein–Hod framework under both coupling scenarios.

## Relation to null_results.tex
Our null_results.tex derives $\kappa_R < 5\times10^{17}\,\mathrm{m^2}$ at $B=10$~T, $R=10^{-26}\,\mathrm{m^{-2}}$ from laboratory nulls with no $\xi$ dependence detected. Karimabadi et al. work in **strong-field regime** (BH horizons, $R \sim M^{-2}$) orders of magnitude larger. The connection:

- **Weak vs. strong curvature:** Laboratory $R \sim 10^{-26}$ suppresses non-minimal coupling effects (our nulls); BH $R \sim 10^{10}\,\mathrm{m^{-2}}$ (solar mass) amplifies them (Karimabadi's QNM shifts).
- **Astrophysical reinterpretation:** Our $\kappa_R$ bound, when propagated to BH environments with $R \sim 10^{10}\,\mathrm{m^{-2}}$ and $B \sim 10^8$~T (magnetars), tightens by $\sim 10^{22}$ orders (as null_results.tex notes). Karimabadi's QNM framework provides observational targets: gravitational ringdown from BH mergers could probe $\kappa_R$ via frequency/damping anomalies.
- **Frequency diagnostics:** They implement dominant-frequency extraction and spectral analysis — precisely the tools our EFQS needs to bridge lab nulls and astrophysical signals.

## How these nulls inform future research
- **QNM-constrained EFT:** Use Karimabadi's QNM($\theta$, $\zeta$) formulas + GW observations (LIGO/Virgo ringdown) to bound non-minimal couplings in strong fields; combine with our lab $\kappa_R$ prior for multi-scale EFT validation.
- **Frequency-domain pipeline:** Implement their WKB-6 + time-domain methods in EFQS; validate against lab-scale analogs (cavity resonances, LC circuits with effective curvature proxies).
- **Instability windows:** Their instability onset at large $\theta/\zeta$ guides parameter-space exclusions; cross-check whether our coherence $\Phi$ or EM configurations approach analogous thresholds.

## Actionable follow-ups
- [ ] Add `dominant_frequency(h_t, dt)` utility to EFQS gravitational_coupling.py; output peak $f$, amplitude, −3 dB bandwidth per Karimabadi's Sec. 4.
- [ ] Notebook: map our $\kappa_R < 5\times10^{17}$ to BH QNM constraints assuming $R \sim (2M)^{-2}$, $B \sim 10^8$~T; estimate detectable frequency shift $\Delta \omega/\omega$ for LIGO-band sources.
- [ ] Reproduce Karimabadi's Fig. 3 (time-domain ringdown) for toy NC-Sch potential in flat-space limit; validate our spectral_derivative and damping fits.
- [ ] Explore "laboratory QNMs": can engineered potentials (superconducting cavities, BEC traps) mimic NC-Sch effective potentials to test coupling scalings in tabletop?

## Notes for EFQS integration
- Implement WKB effective-potential diagnostics: given $V_{\mathrm{eff}}(r; \theta, \xi, \zeta)$, compute QNM estimates and compare to time-domain runs.
- Add "astrophysical recast mode": user inputs $M_{\mathrm{BH}}$, $B_{\mathrm{BH}}$; code scales lab $\kappa_R$ bounds and predicts $\Delta f$ in GW ringdown.
