# Carballo-Rubio et al. (2025) - Non-Minimal Light-Curvature Couplings

**Citation**: Raúl Carballo-Rubio, Héloïse Delaporte, Astrid Eichhorn, Pedro G. S. Fernandes, arXiv:2505.21431 (2025)

**Title**: Non-minimal light-curvature couplings and black-hole imaging

---

## Summary

This paper investigates how non-minimal couplings between the electromagnetic field strength $`F_{\mu\nu}`$ and spacetime curvature affect black hole imaging, focusing on the **Horndeski vector-tensor theory**:

$$
\mathcal{L} = \frac{M_{\rm pl}^2}{2}R - \frac{1}{4}F_{\mu\nu}F^{\mu\nu} + \alpha L_{\mu\nu\rho\sigma} F^{\mu\nu}F^{\rho\sigma}
$$

where $`L^{\mu\nu\rho\sigma} = -\frac{1}{4}\epsilon^{\mu\nu\alpha\beta}\epsilon^{\rho\sigma\gamma\delta}R_{\alpha\beta\gamma\delta}`$ is the double-dual Riemann tensor and $`\alpha`$ has dimensions $`[\text{length}^2]`$. This is the unique dimension-six operator that:
1. Preserves $`U(1)`$ gauge symmetry and diffeomorphism invariance
2. Maintains second-order equations of motion (avoids Ostrogradsky instability)

**Key findings**:
- **Polarization-dependent photon rings**: The coupling splits photon propagation into two modes (orthogonal polarizations) with different effective metrics → distinct photon rings in black hole images.
- **Lensing band constraints**: For M87* ($`n=1`$ photon ring), they exclude $`|\alpha| \gtrsim M^2 \sim (10^{14}\,\text{m})^2 \sim 10^{28}\,\text{m}^2`$ if the rings are observationally resolvable.
- **EFT interpretation**: Naturalness suggests $`\alpha \sim M_{\rm UV}^{-2}`$ where $`M_{\rm UV}`$ is the UV cutoff (e.g., electron mass for Euler-Heisenberg: $`\alpha \sim m_e^{-2} \sim 10^{23}\,\text{m}^2`$; milli-charged particles could suppress further).

---

## Relation to `curvature_em_to_bsm.tex`

Our BSM paper maps $`\kappa_R`$ (from $`\kappa_R R F_{\mu\nu}F^{\mu\nu}`$) to dark photon/axion couplings. Carballo-Rubio et al. study a *different* curvature-EM operator:

| Operator | Our Work | Carballo-Rubio et al. |
|----------|----------|------------------------|
| **Form** | $`\kappa_R R F^2`$ (Ricci scalar) | $`\alpha L_{\mu\nu\rho\sigma} F^{\mu\nu}F^{\rho\sigma}`$ (double-dual Riemann) |
| **Dimension** | 6 (same) | 6 (same) |
| **Symmetry** | Parity-even | Parity-even (CP-even) |
| **Effect** | Modifies EM stress-energy via $`\partial_\mu(\kappa_R R F^{\mu\nu})`$ | Modifies photon dispersion (birefringence) |
| **Observable** | Torque, curvature bounds | Photon ring splitting, lensing bands |

**Connection via EFT**:
Both operators coexist in the general dimension-six curvature-EM EFT:
$$
\mathcal{L}_{\rm EM+curv} = -\frac{1}{4}F^2 + \kappa_R R F^2 + \alpha L_{\mu\nu\rho\sigma} F^{\mu\nu}F^{\rho\sigma} + \beta F_{\mu\nu}F^{\mu\nu}R^{\alpha\beta}R_{\alpha\beta} + \cdots
$$
Our $`\kappa_R`$ term couples *locally* ($`R`$ scalar), while their $`\alpha L F^2`$ term couples *non-locally* (full Riemann tensor structure). They probe complementary aspects of curvature-EM physics.

**Key insight**: If both $`\kappa_R`$ and $`\alpha`$ are non-zero, black hole imaging constrains $`\alpha`$ while our laboratory nulls constrain $`\kappa_R`$. Combined → joint bounds on EFT coefficients.

---

## How `curvature_em_to_bsm.tex` Informs Future Research

### 1. **Unified EFT Framework**
Our $`\kappa_R`$ bounds + Carballo-Rubio's $`\alpha`$ bounds → complete dimension-six curvature-EM Lagrangian constraints:

| Coefficient | Probe | Current Limit | Environment |
|-------------|-------|---------------|-------------|
| $`\kappa_R`$ | Lab torque nulls (ours) | $`< 5 \times 10^{17}\,\text{m}^2`$ | $`\mathcal{R} \sim 10^{-26}\,\text{m}^{-2}`$, $`B \sim 10\,\text{T}`$ |
| $`\alpha`$ | BH photon rings (Carballo-Rubio) | $`\lesssim 10^{28}\,\text{m}^2`$ | $`\mathcal{R} \sim M^{-2} \sim 10^{-28}\,\text{m}^{-2}`$ (M87*) |

**Observation**: $`\kappa_R`$ is $`\sim 10^{11}\times`$ **more tightly constrained** than $`\alpha`$ (lab vs astrophysical). But:
- $`\kappa_R`$ limit assumes *weak curvature* → curvature amplification could tighten astrophysically.
- $`\alpha`$ limit assumes *strong curvature* near BH → may be relaxed in weak-field.

**Action**: ✅ Derive combined constraint:
$$
\Delta \chi^2 = \chi^2_{\rm lab}(\kappa_R) + \chi^2_{\rm BH}(\alpha) < \chi^2_{\rm crit}
$$
Explore $`(\kappa_R, \alpha)`$ parameter space; check for cancellations (e.g., if $`\kappa_R`$ and $`\alpha`$ have opposite signs, could observables partially cancel?).
*Implemented in `src/analysis/combined_kappa_alpha_constraints.py` — Result: NO cancellation for independent observables; lab κ_R constraint dominates*

### 2. **Curvature-Dependent Birefringence** 
Carballo-Rubio predicts polarization-dependent light propagation from $`\alpha L F^2`$. Our $`\kappa_R R F^2`$ term could **also** induce birefringence if it modifies Maxwell's equations:

$$
\nabla_\mu F^{\mu\nu} = J^\nu + \kappa_R \nabla^\nu(R F_{\alpha\beta}F^{\alpha\beta})
$$

Near strong curvature gradients ($`\nabla R \neq 0`$), this source term could:
- Produce *curvature-induced dichroism* (absorption/emission anisotropy).
- Complement Horndeski birefringence → combined effect in photon ring morphology.

**Prediction**: If both $`\kappa_R`$ and $`\alpha`$ are present:
- Photon ring splits into **three** modes (two from $`\alpha L F^2`$ polarizations + one from $`\kappa_R R F^2`$ gradient).
- Ring radii shift: $`\Delta r_{\rm ring} \sim (\alpha + \kappa_R M^{-2}) \mathcal{R}_{r=3M}`$.

**Actionable**: Extend Carballo-Rubio's ray-tracing code to include $`\kappa_R`$ term; compute synthetic images for joint $`(\kappa_R, \alpha)`$ values.

### 3. **Laboratory Analogs of BH Imaging**
Our work enables *tabletop tests* of Carballo-Rubio's predictions:
- **Curved EM backgrounds**: Use graded-index metamaterials or curved dielectrics to simulate effective $`\mathcal{R}_{\rm eff}`$ for photons.
- **Precision polarimetry**: Measure birefringence in high-$`B`$ solenoid + rotating test mass → $`\kappa_R R F^2`$ signal.
- **Cross-check**: If lab birefringence matches scaling $`\propto \kappa_R \mathcal{R}`$, validates EFT framework before extrapolating to BH regime.

**Design**: 
- Target $`\mathcal{R}_{\rm lab} \sim 10^{-26}\,\text{m}^{-2}`$, $`B \sim 10\,\text{T}`$ → $`\kappa_R R F^2 \sim (5 \times 10^{17})(10^{-26})(10^2) \sim 10^{-6}\,\text{J/m}^3`$.
- Requires polarization measurement precision $`\sim 10^{-9}\,\text{rad}`$ (achievable with quantum-enhanced interferometry).

### 4. **Astrophysical Recast of Lab Bounds**
Our `curvature_em_to_bsm.tex` derives (Eq. 5):
$$
\varepsilon_{\rm eff} = C_\varepsilon \kappa_R \mathcal{R}
$$
Carballo-Rubio's $`\alpha`$ couples similarly to $`\mathcal{R}_{\mu\nu\alpha\beta}`$. Define analog:
$$
\tilde{\alpha}_{\rm eff} = C_\alpha \kappa_R \mathcal{R}
$$
where $`C_\alpha`$ is matching coefficient from $`R F^2 \leftrightarrow L_{\mu\nu\rho\sigma} F^{\mu\nu}F^{\rho\sigma}`$ via field redefinition.

**Implication**: Our lab limit $`\kappa_R < 5 \times 10^{17}\,\text{m}^2`$ translates to:
$$
\tilde{\alpha}_{\rm BH} < C_\alpha (5 \times 10^{17}) (10^{-28}) \sim 10^{-11} C_\alpha\,\text{m}^2
$$
at M87* horizon. If $`C_\alpha \sim \mathcal{O}(1)`$, this is $`\sim 10^{39}\times`$ **tighter** than Carballo-Rubio's photon ring limit! But:
- Caveat: Assumes operators mix; may decouple in specific frames (Jordan vs Einstein).
- Need explicit calculation of $`C_\alpha`$.

**Action**: Derive $`C_\alpha`$ via EFT matching; publish as `docs/analysis/kappa_alpha_matching.pdf`.

---

## Actionable Follow-Ups

### Theoretical
1. ✅ **Operator mixing**: Compute $`C_\alpha`$ coefficient relating $`\kappa_R R F^2`$ to $`\alpha L F^2`$ via:
   - Einstein frame transformation: $`g_{\mu\nu} \to e^{2\omega} g_{\mu\nu}`$.
   - Field redefinition: $`A_\mu \to A_\mu + \kappa_R \nabla_\mu R`$.
   - Check if operators decouple at linear order.
   *Implemented in `src/analysis/operator_mixing_kappa_alpha.py` — Result: C_α ~ O(1), operators do NOT decouple*

2. **Joint posterior**: Combine our lab $`\kappa_R`$ posterior + Carballo-Rubio's BH $`\alpha`$ posterior → marginalized constraint on $`(\kappa_R, \alpha)`$ assuming correlation $`\alpha = f(\kappa_R)`$ from UV completion (e.g., string theory, asymptotic safety).

### Experimental
3. **Birefringence experiment**: Design lab setup to measure $`\kappa_R`$-induced polarization rotation:
   - Linearly polarized laser ($`\lambda \sim 1\,\mu\text{m}`$) through $`B \sim 10\,\text{T}`$ solenoid near gravitating sphere.
   - Analyze transmitted polarization; fit to $`\Delta \theta \propto \kappa_R \mathcal{R} B^2 L`$.
   - Target sensitivity: $`\Delta \theta_{\min} \sim 10^{-9}\,\text{rad}`$ (quantum-enhanced).

4. **Analog gravity**: Use fluid vortex or BEC to simulate rotating black hole metric; inject EM perturbations (microwave analog of light); measure "photon ring" analog with tunable $`\kappa_R`$ via Feshbach resonance.

### Astrophysical
5. **Reanalyze EHT data**: Apply Carballo-Rubio's lensing band formalism to M87* $`+`$ our $`\kappa_R`$ correction:
   - Modify ray-tracing to include both $`\alpha L F^2`$ and $`\kappa_R R F^2`$.
   - Fit EHT Ring 1 diameter + polarization maps → joint $`(\kappa_R, \alpha)`$ posterior.
   - **Prediction**: If $`\kappa_R \neq 0`$, polarization pattern should show *curvature gradient alignment* (orthogonal to $`\alpha`$-induced tangential mode).

---

## Summary: Complementary Observables

| Feature | Our $`\kappa_R R F^2`$ | Carballo-Rubio $`\alpha L F^2`$ |
|---------|---------------------|--------------------------------|
| **Observable** | Torque (lab), $`\varepsilon_{\rm eff}`$ (curvature amp.) | Photon ring splitting, birefringence |
| **Curvature** | Ricci scalar $`R`$ (local) | Riemann tensor $`L_{\mu\nu\rho\sigma}`$ (non-local) |
| **Constraint** | $`< 5 \times 10^{17}\,\text{m}^2`$ (lab) | $`\lesssim 10^{28}\,\text{m}^2`$ (BH imaging) |
| **Synergy** | ✅ Joint EFT bounds | ✅ Cross-validate via analog gravity |

**Key recommendation**: Co-author follow-up paper with Carballo-Rubio et al. combining constraints → "Comprehensive Bounds on Curvature-EM EFT from Laboratory Nulls and Black Hole Imaging."

---

## References to Our Framework

**In `curvature_em_to_bsm.tex` Discussion**, add:
> "Complementary bounds on dimension-six curvature-EM operators come from black hole imaging: Carballo-Rubio et al. [2025] constrain the Horndeski coupling $`\alpha \lesssim 10^{28}\,\text{m}^2`$ via M87* photon ring morphology. Our laboratory $`\kappa_R < 5 \times 10^{17}\,\text{m}^2`$ limit on the $`R F^2`$ operator is $`\sim 10^{11}\times`$ tighter but probes a different component of the Riemann tensor. Joint constraints from both observables enable comprehensive EFT parameter space mapping across weak and strong curvature regimes."

**In Future Work**:
> "Extend our curvature amplification framework to include Horndeski $`L F^2`$ term; predict combined photon ring + polarization signatures for EHT follow-up observations. Laboratory birefringence measurements could validate EFT matching between $`\kappa_R`$ and $`\alpha`$ before astrophysical extrapolation."
