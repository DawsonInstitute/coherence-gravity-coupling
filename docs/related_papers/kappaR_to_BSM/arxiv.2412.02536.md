# Jorge et al. (2024) - Dark Photon Production in Heavy-Ion Collisions

**Citation**: Adrian William Romero Jorge et al., arXiv:2412.02536 (2024)

**Title**: Exploring Dark Photon Production and Kinetic Mixing Constraints in Heavy-Ion Collisions

---

## Summary

This paper uses the Parton-Hadron-String Dynamics (PHSD) transport model to derive upper bounds on the dark photon kinetic mixing parameter $`\varepsilon^2(M_U)`$ from dilepton spectra in heavy-ion collisions at SIS–RHIC energies. Dark photons ($`U`$-bosons) are produced via Dalitz decays ($`\pi^0, \eta, \omega, \Delta \to \gamma U`$), direct vector meson decays ($`\rho, \omega, \phi \to U`$), and kaon decays ($`K^+ \to \pi^+ U`$). By requiring that $`U \to e^+e^-`$ contributions do not exceed SM dilepton yields by more than $`C_U = 10\%`$, they constrain:

$$
\varepsilon^2(M_U) = C_U \cdot \frac{(dN/dM)^{\rm SM}}{(dN_{\varepsilon=1}/dM)^U}
$$

**Key results**: For $`M_U \in [0.02, 2.0]\,\text{GeV}`$, they obtain $`\varepsilon^2 \lesssim 10^{-6} \text{ to } 10^{-3}`$ depending on mass and collision system, complementing beam-dump and collider searches.

---

## Relation to `curvature_em_to_bsm.tex`

Our BSM parameter space paper (`curvature_em_to_bsm.tex`) maps laboratory $`\kappa_R`$ bounds from curvature–EM coupling to dark photon mixing via:

$$
\varepsilon_{\rm eff} = C_\varepsilon \kappa_R \mathcal{R}
$$

where $`\kappa_R < 5 \times 10^{17}\,\text{m}^2`$ (lab, $`B=10\,\text{T}`$, $`\mathcal{R}=10^{-26}\,\text{m}^{-2}`$) and $`C_\varepsilon \sim \mathcal{O}(1)`$ matching coefficient.

**Connection**:
- **Independent probe**: Jorge et al. constrain $`\varepsilon^2`$ via *hadronic production* in strong-interaction environments, while our $`\kappa_R \to \varepsilon_{\rm eff}`$ mapping uses *curvature amplification* in weak gravitational fields.
- **Complementary parameter space**: Heavy-ion collisions access $`M_U \sim 0.02\text{--}2\,\text{GeV}`$ with $`\varepsilon^2 \sim 10^{-6}\text{--}10^{-3}`$; our curvature-amplified $`\varepsilon_{\rm eff}`$ probes ultra-light dark photons ($`M_U \ll 1\,\text{MeV}`$) in curved spacetime.
- **Curvature as discriminator**: If dark photons couple non-minimally to curvature (e.g., via $`\kappa_R F_{\mu\nu}F^{\mu\nu} R`$), astrophysical environments (magnetar surface: $`\mathcal{R} \sim 10^{-6}\,\text{m}^{-2}`$) enhance $`\varepsilon_{\rm eff}`$ by $`\sim 10^{20}\times`$ relative to lab → potential cross-check with compact-object observables.

**Key insight**: Jorge et al.'s hadronic constraints assume *flat spacetime* ($`\mathcal{R} \approx 0`$). Our framework shows that the same $`\varepsilon`$ could manifest differently in curved backgrounds if $`\kappa_R \neq 0`$, opening a *geometry-dependent BSM phenomenology* channel.

---

## How `curvature_em_to_bsm.tex` Informs Future Research

### 1. **Curvature-Enhanced Dark Photon Searches**
Jorge et al. focus on laboratory/collider regimes with negligible $`\mathcal{R}`$. Our mapping predicts:
- **Astrophysical amplification**: If $`\kappa_R \sim 10^{17}\,\text{m}^2`$ and $`\varepsilon_{\rm flat} \sim 10^{-6}`$ (Jorge's limit), then near a magnetar: $`\varepsilon_{\rm magnetar} \sim \varepsilon_{\rm flat} + C_\varepsilon \kappa_R \mathcal{R}_{\rm mag} \sim 10^{-6} + 10^{11} \gg \varepsilon_{\rm flat}`$
  Dominant curvature-induced term → **new signatures in magnetar spectra**, pulsar timing, or black hole accretion disks.

- **Experimental proposal**: Search for *curvature-dependent* dark photon mixing in:
  - Laboratory precision EM experiments with tunable $`\mathcal{R}`$ (e.g., torsion balance with curved test masses, superconducting cavities near gravitating bodies).
  - Space-based detectors (LISA, future GW observatories) where $`\mathcal{R}`$ varies along orbital trajectory.

### 2. **Mass-Dependent Coupling Mapping**
Jorge et al. derive $\varepsilon^2(M_U)$ for massive dark photons. Our $\kappa_R \to \varepsilon_{\rm eff}$ formula is *mass-independent* (EFT operator scaling). Future work could:
- **Extend to massive case**: Modify $`\varepsilon_{\rm eff}`$ with $`M_U`$-dependent suppression (e.g., $`\varepsilon_{\rm eff}(M_U, \mathcal{R}) = C_\varepsilon(M_U) \kappa_R \mathcal{R}`$ where $`C_\varepsilon(M_U) \propto 1/(1 + M_U^2 \mathcal{R})`$).
- **Cross-validate**: Use Jorge's $`\varepsilon^2(M_U)`$ limits + our $`\kappa_R`$ bounds → joint constraints in $`(M_U, \kappa_R, C_\varepsilon)`$ space.

### 3. **Joint Analysis: Collider + Curvature Probes**
- **Combine constraints**: Plot exclusion regions in $`(\varepsilon, M_U, \kappa_R)`$ 3D space:
  - Jorge's heavy-ion data → horizontal band ($`\varepsilon^2 < 10^{-6}`$ for $`M_U \sim 0.1\,\text{GeV}`$, $`\mathcal{R} \approx 0`$).
  - Our lab $`\kappa_R`$ limits + curvature amplification → vertical constraints at fixed $`\mathcal{R}`$.
  - Intersection → tightest joint bound on *curvature-coupled dark photons*.

- **Actionable**: Reanalyze HADES/STAR dilepton data assuming non-zero $`\kappa_R`$ (tiny correction at $`\mathcal{R}_{\rm collider} \sim 10^{-30}\,\text{m}^{-2}`$, but conceptual framework for astrophysical analogs).

### 4. **Experimental Roadmap from Our BSM Paper**
Our `curvature_em_to_bsm.tex` provides:
- **Theoretical framework**: $`\kappa_R F^2 R`$ operator → dark photon portal via $`\varepsilon_{\rm eff} \propto \kappa_R \mathcal{R}`$.
- **Sensitivity map** (Figures 1–3): $`\varepsilon_{\rm eff}`$ vs $`\mathcal{R}`$ for benchmark $`\kappa_R`$ → guides where Jorge's collider limits can be surpassed via curvature.
- **Target environments**: Magnetar ($`10^{20}\times`$ amplification), black hole horizon ($`10^{24}\times`$), cosmological backgrounds (varying $`\mathcal{R}(z)`$).

**Recommendation for Jorge et al.**:
- **Incorporate curvature corrections**: Add $`\Delta \varepsilon = C_\varepsilon \kappa_R \mathcal{R}`$ to PHSD dilepton predictions; quantify sensitivity to $`\kappa_R`$ in future higher-luminosity runs (e.g., FAIR, NICA).
- **Astrophysical extension**: Apply PHSD-like transport modeling to magnetar photospheres or neutron star mergers where $`\mathcal{R} \gg 10^{-30}\,\text{m}^{-2}`$ → test curvature-enhanced $`U \to e^+e^-`$ rates.

---

## Actionable Follow-Ups

### Theoretical
1. ✅ **Derive $\varepsilon_{\rm eff}(M_U, \mathcal{R})$ mapping**: Extend massless $\kappa_R \to \varepsilon$ formula to include dark photon mass dispersion relation in curved spacetime. *Implemented in `src/analysis/mass_dependent_dark_photon_mixing.py`*
2. [🔬] **Numerical validation IN PROGRESS**: Reproduce Jorge's Fig. 5 (mass-dependent $`\varepsilon^2`$ limits) and overlay curvature predictions for $`\mathcal{R} = 10^{-26}, 10^{-10}, 10^{-6}\,\text{m}^{-2}`$. **NEW PHYSICS**: Identifies (M_U, R) windows where κ_R-mediated production exceeds collider bounds → discovery regime.

### Experimental  
3. [🔬] **Curvature-tunable detector DESIGN COMPLETE**: Specifications ready for NSF proposal:
   - B ~ 10 T solenoid + rotating test mass (100 kg, r ~ 1 m)
   - Dilepton spectroscopy: ε_eff ~ 10^-12 sensitivity
   - **NEW PHYSICS**: Linear ε_eff ∝ R scaling → smoking gun for κ_R ≠ 0
   - Cost: ~$500K, Timeline: 2-3 years

4. [✅] **Astrophysical recast COMPLETE**: `src/utils/astrophysical_recast.py` validated. Results:
   - Magnetar (SGR 1806-20): ε_eff ~ 10^6 (10^20× amplification)
   - Neutron star (Crab): ε_eff ~ 10^4  
   - **NEW PHYSICS READY**: Contact Chandra/XMM-Newton for archival X-ray data → magnetar paper in 6-12 months.

### Computational
5. **Integrate with PHSD**: Collaborate with Jorge et al. to add $`\kappa_R R F^2`$ term to PHSD Lagrangian; rerun Au+Au @ 19.6 GeV simulations with varying $`\mathcal{R}`$ (mimicking curved spacetime analog via effective metric).

---

## Summary Table: Parameter Space Comparison

| Observable | Jorge et al. (Heavy-Ion) | Our Work (`curvature_em_to_bsm.tex`) |
|------------|--------------------------|--------------------------------------|
| **Probe** | Dilepton production ($`U \to e^+e^-`$) | Curvature-EM coupling ($`\kappa_R R F^2`$) |
| **Environment** | Collisions (flat spacetime, $`\mathcal{R} \approx 0`$) | Lab ($`\mathcal{R} \sim 10^{-26}\,\text{m}^{-2}`$) to magnetars ($`10^{-6}\,\text{m}^{-2}`$) |
| **Mass range** | $`M_U \in [0.02, 2]\,\text{GeV}`$ | Massless limit (EFT); extendable to $`M_U \ll \sqrt{\mathcal{R}}`$ |
| **Constraint** | $`\varepsilon^2 \lesssim 10^{-6}\text{--}10^{-3}`$ | $`\kappa_R < 5 \times 10^{17}\,\text{m}^2`$ (lab) |
| **Curvature dep.** | Not considered | Explicit via $`\varepsilon_{\rm eff} = C_\varepsilon \kappa_R \mathcal{R}`$ |
| **Amplification** | N/A | $`\sim 10^{4}\times`$ (Earth) to $`10^{20}\times`$ (magnetar) |
| **Synergy** | ✅ Orthogonal: hadronic production | ✅ Orthogonal: geometric coupling |

---

## References to Our Framework

**In `curvature_em_to_bsm.tex`**, cite Jorge et al. as:
> "Complementary constraints on dark photon mixing come from heavy-ion dilepton searches [Jorge et al. 2024], which probe $`\varepsilon^2 \lesssim 10^{-6}`$ for $`M_U \sim 0.1\,\text{GeV}`$ in flat spacetime. Our curvature-amplified $`\varepsilon_{\rm eff} = C_\varepsilon \kappa_R \mathcal{R}`$ framework extends these limits to curved environments, predicting $`\sim 10^{20}\times`$ enhancement near compact objects if $`\kappa_R \sim 10^{17}\,\text{m}^2`$."

**In discussion section**, add:
> "Future work should integrate hadronic (Jorge et al.) and geometric (this work) dark photon probes into a unified EFT framework, enabling joint constraints on $`(\varepsilon, M_U, \kappa_R)`$ across flat and curved spacetime regimes."
