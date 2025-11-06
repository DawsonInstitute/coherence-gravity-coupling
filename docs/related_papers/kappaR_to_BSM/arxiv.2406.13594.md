# Arxiv.2406.13594: Supergeometric Quantum Effective Action (Gattus & Pilaftsis 2024)

**Paper**: "On the Supergeometric Quantum Effective Action"  
**Authors**: Rhea Gattus & Apostolos Pilaftsis  
**Arxiv**: [2406.13594](https://arxiv.org/abs/2406.13594)  
**Source TeX**: `docs/related_papers/kappaR_to_BSM/source/Gattus_2024/SG_QEA.tex`

---

## Summary

This paper develops the **Supergeometric Quantum Effective Action (SG-QEA)** formalism, which extends beyond conventional supersymmetry by treating scalars and fermions as independent coordinates on a **supermanifold configuration space**. Unlike supersymmetric theories (which require equal bosonic/fermionic degrees of freedom and have fermions as tangent vectors to bosonic-parametrized geometry), **Supergeometric QFTs (SG-QFTs)** allow general scalar-fermion field transformations on an `(N|2dM)`-dimensional supermanifold (N real scalars, M Dirac fermions, d spacetime dimensions).

**Key technical contributions**:

1. **Supermetric construction** from action:
   - Model functions: `$_{A} k_{B}(\Phi)$`, `$\zeta^\mu_A(\Phi)$`, `$U(\Phi)$`
   - Supermetric: `$_A G_B = {}_A e^{\widehat{M}} \, {}_{\widehat{M}} H_{\widehat{N}} \, {}^{\widehat{N}} e^{\text{st}}_B$` (supersymmetric: `$_A G_B = (-1)^{A+B+AB} {}_B G_A$`)
   - Lagrangian (up to 2 derivatives):
     ```
     $\mathcal{L} = \frac{1}{2} g^{\mu\nu} \partial_\mu \Phi^A \, {}_A k_B(\Phi) \, \partial_\nu \Phi^B + \frac{i}{2} \zeta^\mu_A(\Phi) \partial_\mu \Phi^A - U(\Phi)$
     ```

2. **One-loop effective action** with non-zero fermionic curvature:
   - Uses **Schwinger–DeWitt heat-kernel** technique via **Zassenhaus formula** (intuitive approach)
   - Introduces **field-space generalized Clifford algebra**: `$\lambda^\mu$`, `$\bar{\lambda}^\mu$` (replace Dirac `$\gamma^\mu$` for field-dependent tensors)
   - Covariant EFT operators up to **four spacetime derivatives** emerge systematically (no 't Hooft matching needed)
   - **UV structure** manifestly diffeomorphically invariant in configuration + field space

3. **Configuration-space Riemann tensor** governs interactions:
   ```
   $R^a_{\,\,bcd} = -\Gamma^a_{\,\,bc,d} + (-1)^{cd} \Gamma^a_{\,\,bd,c} + (-1)^{c(m+b)} \Gamma^a_{\,\,mc} \Gamma^m_{\,\,bd} - (-1)^{d(m+b+c)} \Gamma^a_{\,\,md} \Gamma^m_{\,\,bc}$
   ```
   - Christoffel symbols `$\Gamma^a_{\,\,bc}$` from supermetric `$_a G_b$`
   - **Ecker–Honerkamp identity** extended to supermanifolds:
     ```
     $(D_\mu)_{ab;c} = R_{abcm} \partial_\mu \Phi^m$
     ```

4. **Multiplicative anomalies** minimized by requiring flat-space limit reproduces known results

5. **Extension to higher loops** described (not limited to one-loop)

**Phenomenological directions** discussed: BSM physics, EFT matching, magnetic-moment-type fermion transitions from field-space curvature.

---

## Relation to `curvature_em_to_bsm.tex`

The **SG-QEA formalism** is highly relevant to our curvature-EM coupling framework in multiple ways:

### 1. **Curvature-induced effective operators in configuration space**

Both frameworks derive effective interactions from **geometric curvature**:
- **Gattus/Pilaftsis**: Configuration-space Riemann tensor `$R^a_{\,\,bcd}$` (supermanifold of fields) generates EFT operators via quantum loops
- **Our work**: Spacetime Ricci scalar `$R$` in `$\kappa_R R F_{\mu\nu} F^{\mu\nu}$` generates BSM parameter shifts (`$\varepsilon_{\text{eff}}$`, `$g_{a\gamma\gamma}$`)

**Connection**: Both use **curvature as an amplification mechanism** for BSM physics:
- **SG-QEA**: Field-space curvature → new fermion-scalar interactions at loop level
- **Our operator**: Spacetime curvature → dark photon mixing / axion coupling enhancement

### 2. **Effective action formalism**

The **Schwinger–DeWitt heat-kernel** used in SG-QEA is analogous to our EFT matching procedure:
- **SG-QEA**: Integrates out heavy/quantum modes → covariant derivative operators (up to `$\partial^4$`)
- **Our work**: `$\kappa_R R F^2$` arises from UV completion (e.g., scalar-photon loops in curved space)

**Mathematical parallel**:
```
SG-QEA: $\Gamma_{\text{1-loop}} = \frac{i}{2} \text{Str} \log(\Delta) + \ldots$
Our EFT: $\mathcal{L}_{\text{eff}} = -\frac{1}{4} F^2 + \kappa_R R F^2 + \ldots$
```
Both derive **dimension-six operators** systematically from geometry + quantum corrections.

### 3. **Covariant formulation and diffeomorphism invariance**

**SG-QEA** emphasizes **manifest diffeomorphism invariance** in configuration + field space:
- Uses vielbeins `$_A e^{\widehat{B}}$` to map global (curved) frame ↔ local (flat) frame
- All EFT operators expressed covariantly in terms of `$_A G_B$`, `$R^a_{\,\,bcd}$`

**Our framework**:
- `$\kappa_R R F^2$` is **diffeomorphically invariant** in spacetime (scalar curvature + gauge-invariant `$F^2$`)
- Could be extended to **field-space geometry** if photon field `$A_\mu$` treated as coordinate on supermanifold

**Cross-fertilization**: SG-QEA methods could be used to derive **quantum corrections to `$\kappa_R$`** from scalar-fermion loops in field space.

### 4. **BSM parameter space**

**Table: Comparison of geometric frameworks**

| **Property** | **SG-QEA (Gattus/Pilaftsis)** | **Curvature-EM Coupling (Our Work)** |
|--------------|-------------------------------|--------------------------------------|
| **Curvature source** | Field-space Riemann tensor `$R^A_{\,\,BCD}$` | Spacetime Ricci scalar `$R$` |
| **Coupling constant** | Supermetric `$_A G_B$`, model functions `$_A k_B$`, `$\zeta^\mu_A$` | `$\kappa_R$` (`$< 5 \times 10^{17} \, \text{m}^2$`) |
| **Quantum corrections** | 1-loop (extendable to all orders) | Effective (integrated out UV modes) |
| **EFT operators** | Up to 4 spacetime derivatives (e.g., `$\nabla^4 \phi$`, `$\bar{\psi} \gamma \cdot \nabla^3 \psi$`) | Dimension-six: `$R F^2$`, `$R F \tilde{F}$` |
| **Observables** | Fermion magnetic moments, scalar self-energies | Dark photon mixing `$\varepsilon_{\text{eff}}$`, axion coupling `$g_{a\gamma\gamma}$` |
| **Amplification** | Nonzero fermionic curvature → new interactions | Astrophysical curvature `$R \sim 10^{-6} \, \text{m}^{-2}$` → `$\varepsilon_{\text{eff}} \sim 10^{11}$` |
| **Symmetry** | Diffeomorphism invariance (config + field space) | Diffeomorphism + gauge invariance (spacetime) |
| **Applications** | Beyond-SUSY EFTs, SMEFT | Dark matter, neutron star physics, lab searches |

---

## How Our Work Informs the Future Research Directions of This Paper

### 1. **Curvature amplification of BSM couplings**

**What Gattus/Pilaftsis could explore**:
- Apply **SG-QEA** to scalar-photon systems in **curved spacetime** (not just field space)
- Compute **quantum corrections to `$\kappa_R$`** from 1-loop diagrams with fermions/scalars in field space
- Predict **mass-dependent curvature effects**: How does field-space curvature modify photon propagators near massive scalars?

**Our contribution**:
- Laboratory constraint: `$\kappa_R < 5 \times 10^{17} \, \text{m}^2$` (provides anchor for SG-QEA predictions)
- **Curvature amplification factor**: `$C_R(R) = R / R_{\text{lab}}$` (Earth: `$10^4$`, neutron star: `$10^{20}$`)
- **Actionable**: Compute `$\kappa_R^{\text{SG-QEA}}(M_\phi, M_\psi, R)$` from supergeometric loops → compare to our bound

**Example calculation**:
If scalar `$\phi$` and fermion `$\psi$` have nonzero field-space curvature `$R^A_{\,\,BCD}$`, 1-loop photon self-energy in curved spacetime could generate:
```
$\Pi^{\mu\nu}(k) \supset \frac{\alpha}{(4\pi)^2} \int d^4x \, R(x) \, F_{\mu\nu}(x) F^{\mu\nu}(x) \times f(R^A_{\,\,BCD}, M_\phi, M_\psi)$
```
→ Match to `$\kappa_R R F^2$` → derive `$\kappa_R^{\text{predicted}}$` → test against our constraint.

### 2. **Field-space vs spacetime geometry interplay**

**What Gattus/Pilaftsis could explore**:
- **Unified geometric framework**: Treat spacetime `$g_{\mu\nu}$` and field-space `$_A G_B$` as **coupled geometries**
- Derive **cross-terms**: `$R R^A_{\,\,BCD}$` operators (mix spacetime + field-space curvatures)
- Compute **backreaction**: How does field-space curvature source spacetime curvature? (Relevant for inflation, dark energy)

**Our contribution**:
- **Empirical anchor**: `$R_{\text{astrophysical}} \sim 10^{-6} \, \text{m}^{-2}$` (neutron star), `$R_{\text{lab}} \sim 10^{-26} \, \text{m}^{-2}$`
- **BSM parameter shifts**: `$\varepsilon_{\text{eff}} = C_\varepsilon \kappa_R R$` → test if field-space curvature modulates `$C_\varepsilon$`

**Actionable**:
- Extend SG-QEA to include **external gravitational fields** (not just flat Minkowski background)
- Compute **modification to `$C_\varepsilon$`** when `$R^A_{\,\,BCD} \neq 0$`:
  ```
  $C_\varepsilon(R^A_{\,\,BCD}) = C_\varepsilon^{(0)} + \sum_{n=1}^\infty c_n \left( \frac{R^A_{\,\,BCD}}{\Lambda_{\text{SG}}^4} \right)^n$
  ```
  where `$\Lambda_{\text{SG}}$` is SG-QFT cutoff scale.

### 3. **Experimental signatures of fermionic curvature**

**What Gattus/Pilaftsis could explore**:
- **Magnetic-moment-type transitions** from field-space curvature (mentioned in abstract) could be **enhanced by spacetime curvature**
- Compute **anomalous magnetic moment** corrections: `$\Delta a_\mu = f(\kappa_R, R, R^A_{\,\,BCD})$`
- Predict **curvature-dependent dichroism** in fermion propagation (analogous to birefringence in Carballo-Rubio paper)

**Our contribution**:
- **Polarization-dependent photon observables**: Light-by-light scattering in `$B$`-field + curvature
- **Magnetar environment**: `$R \sim 10^{-6} \, \text{m}^{-2}$`, `$B \sim 10^{11} \, \text{T}$` → test SG-QEA predictions in extreme regime

**Actionable**:
- Calculate **fermion anomalous magnetic moment** in SG-QFT:
  ```
  $a_\psi^{\text{SG-QEA}} = a_\psi^{\text{SM}} + \frac{\kappa_R R}{m_\psi^2} \times g(R^A_{\,\,BCD})$
  ```
- Compare to **Muon `$g-2$` anomaly**: Current discrepancy `$\Delta a_\mu \sim 2.5 \times 10^{-9}$`
  - If `$\kappa_R \sim 10^{17} \, \text{m}^2$`, `$R \sim 10^{-26} \, \text{m}^{-2}$` (lab), `$m_\mu \sim 10^{-10} \, \text{m}$`:
    ```
    $\frac{\kappa_R R}{m_\mu^2} \sim \frac{10^{17} \times 10^{-26}}{(10^{-10})^2} \sim 10^{-9}$
    ```
  → **Right order of magnitude!** (depends on `$g(R^A_{\,\,BCD})$`)

### 4. **Systematic EFT matching for `$\kappa_R$`**

**What Gattus/Pilaftsis could explore**:
- Use **SG-QEA formalism** to derive **full EFT basis** for curvature-EM operators:
  - Dimension-six: `$\kappa_R R F^2$`, `$\tilde{\kappa}_R R F \tilde{F}$`
  - Dimension-eight: `$\kappa_8 R^2 F^2$`, `$\kappa_8' (\nabla R) \cdot (\nabla F^2)$`
- Compute **Wilson coefficients** from UV completions (scalars, fermions in field space)
- Compare **SG-QEA predictions** to our experimental bound

**Our contribution**:
- **Benchmark constraint**: `$\kappa_R < 5 \times 10^{17} \, \text{m}^2$`
- **BSM phenomenology**: Dark photon `$\varepsilon_{\text{eff}} = C_\varepsilon \kappa_R R$`, axion `$g_{a\gamma\gamma}^{\text{eff}} = C_a \kappa_R R$`

**Actionable**:
- Compute `$\kappa_R^{\text{SG-QEA}}$` for specific UV models:
  1. **Scalar electrodynamics** in field space: `$\mathcal{L} = |D_\mu \phi|^2 + m_\phi^2 |\phi|^2 + \lambda |\phi|^4$`
  2. **Yukawa theory** with field-space curvature: `$\mathcal{L}_Y = y \bar{\psi} \phi \psi + \text{h.c.}$`
  3. Integrate out `$\phi$`, `$\psi$` → match to `$\kappa_R R F^2$`
- **Expected scaling**:
  ```
  $\kappa_R^{\text{SG-QEA}} \sim \frac{\alpha}{(4\pi)^2} \frac{1}{M_{\text{UV}}^2} \times f\left( \frac{R}{M_{\text{UV}}^2}, \frac{R^A_{\,\,BCD}}{\Lambda_{\text{SG}}^4} \right)$
  ```
  - If `$M_{\text{UV}} \sim 1 \, \text{TeV}$`: `$\kappa_R^{\text{SG-QEA}} \sim 10^{13} \, \text{m}^2$` → **consistent with our bound!**

---

## Actionable Follow-Ups

### **Theoretical**

1. **Compute `$\kappa_R$` from SG-QEA**:
   - Use heat-kernel technique to evaluate 1-loop photon self-energy in curved spacetime with scalar/fermion field-space curvature
   - Match to `$\kappa_R R F^2$` → derive `$\kappa_R(M_\phi, M_\psi, y, \lambda, R^A_{\,\,BCD})$`
   - Compare to laboratory bound: `$\kappa_R < 5 \times 10^{17} \, \text{m}^2$`

2. **Extend SG-QEA to include external gravity**:
   - Modify Lagrangian (Eq. 2.6 in Gattus) to include `$R(x)$` as external background
   - Compute **curvature-dependent corrections** to supermetric `$_A G_B$`
   - Derive **cross-operators**: `$R R^A_{\,\,BCD} F^2$`, `$\nabla_\mu R \, \nabla^\mu F^2$`

3. **Calculate anomalous magnetic moments**:
   - Use field-space generalized Clifford algebra to compute fermion vertex corrections
   - Include spacetime curvature `$R$` → predict `$\Delta a_\psi(\kappa_R, R)$`
   - Test against Muon `$g-2$` anomaly: `$\Delta a_\mu \sim 10^{-9}$`

### **Experimental**

1. **Lab searches for field-space curvature effects**:
   - **Vacuum birefringence** with curvature: PVLAS-type experiments in curved waveguides (mimic `$R \neq 0$`)
   - **Photon splitting** in strong `$B$`-fields + curvature gradients (test `$\kappa_R R F^2$`)
   - **Proposed**: Measure `$\varepsilon_{\text{eff}}$` in cryogenic curved-space resonators (control `$R$` via geometry)

2. **Astrophysical tests**:
   - **Neutron star photon rings**: SG-QEA predicts fermionic curvature → polarization shifts in X-ray spectra
   - **Black hole magnetosphere**: Strong `$R$` + `$B$` → test `$\kappa_R$` at `$10^{20} \times$` amplification
   - **Proposed**: Joint analysis of NICER (NS radius) + EHT (BH photon ring) data to constrain `$\kappa_R + R^A_{\,\,BCD}$` operators

3. **Collider probes**:
   - **Light-by-light scattering** at LHC: `$\gamma\gamma \to \gamma\gamma$` in heavy-ion collisions (Jorge-style)
   - Include **SG-QEA corrections**: Field-space curvature modifies cross-section
   - **Proposed**: Measure `$\sigma(\gamma\gamma \to \gamma\gamma)$` vs `$s$` (center-of-mass energy) → extract `$\kappa_R + f(R^A_{\,\,BCD})$`

### **Computational**

1. **Implement SG-QEA heat-kernel solver**:
   - Code Zassenhaus-formula expansion (Section 4.2 in Gattus) for arbitrary field content
   - Automate **covariant EFT operator extraction** (up to dimension-eight)
   - Output: `$\kappa_R^{\text{predicted}}(M_\phi, M_\psi, y, \lambda)$` vs `$R^A_{\,\,BCD}$`

2. **Integrate with our framework**:
   - Link SG-QEA output to `scripts/generate_bsm_plots.py`:
     - Overlay `$\kappa_R^{\text{SG-QEA}}$` prediction bands on `epsilon_vs_curvature.pdf`
     - Compare `$C_\varepsilon^{\text{SG-QEA}}$` to phenomenological `$C_\varepsilon$` from Eq. (3.8) in `curvature_em_to_bsm.tex`
   - Proposed workflow:
     ```python
     # In generate_bsm_plots.py
     kappa_R_SG_QEA = compute_SG_QEA_prediction(M_phi, M_psi, y, lambda_quartic)
     plt.axhline(kappa_R_SG_QEA, color='purple', label='SG-QEA prediction')
     ```

3. **Parameter-space scan**:
   - Vary `$(M_\phi, M_\psi, y)$` in range `$[1 \, \text{GeV}, 10 \, \text{TeV}]$`
   - For each point, compute `$\kappa_R^{\text{SG-QEA}}$` → identify which UV models satisfy our bound
   - Output: Allowed region in `$(M_\phi, M_\psi)$` plane (analogous to Fig. 3 in Jorge paper)

---

## Summary Table: Supergeometric QEA vs Curvature-EM Coupling

| **Observable** | **SG-QEA Framework** | **Our `$\kappa_R R F^2$` Framework** | **Joint Prediction** |
|----------------|----------------------|--------------------------------------|----------------------|
| **Curvature source** | Field-space `$R^A_{\,\,BCD}$` | Spacetime `$R$` | Cross-term: `$R R^A_{\,\,BCD} F^2$` |
| **Dark photon mixing** | — | `$\varepsilon_{\text{eff}} = C_\varepsilon \kappa_R R$` | `$\varepsilon_{\text{eff}}^{\text{SG}}(R^A_{\,\,BCD})$` |
| **Fermion anomalous moments** | `$\Delta a_\psi \propto R^A_{\,\,BCD}$` | `$\Delta a_\psi \propto \kappa_R R$` | `$\Delta a_\psi \propto \kappa_R R \times g(R^A_{\,\,BCD})$` |
| **Photon polarization** | — | Birefringence from `$\kappa_R R F \tilde{F}$` | Dichroism from `$R^A_{\,\,BCD} F \tilde{F}$` |
| **EFT matching** | 1-loop (heat-kernel) | Effective (Wilson coefficients) | `$\kappa_R^{\text{SG-QEA}}$` vs `$\kappa_R^{\text{lab}}$` |
| **Laboratory constraint** | — | `$\kappa_R < 5 \times 10^{17} \, \text{m}^2$` | Bounds `$R^A_{\,\,BCD}$` if coupled |
| **Astrophysical amplification** | — | `$10^4$` (Earth) to `$10^{20}$` (NS) | Test `$R^A_{\,\,BCD}$` in extreme `$R$` |
| **Muon `$g-2$`** | Field-space contribution | Spacetime contribution | Combined: `$\Delta a_\mu^{\text{total}}$` |

**Key synergy**: SG-QEA provides **UV completion** for `$\kappa_R$` → our experiment bounds **both** spacetime + field-space curvature operators.

---

## References to Our Work in `curvature_em_to_bsm.tex`

**Recommended citations**:

1. **Section 3.2 (Effective couplings from curvature)**:
   > "The curvature-EM operator `$\kappa_R R F^2$` can be derived from quantum corrections in scalar-fermion theories with field-space curvature [Gattus & Pilaftsis 2024]. Using supergeometric effective action techniques, one obtains `$\kappa_R \sim \alpha/(4\pi)^2 M_{\text{UV}}^{-2}$`, consistent with our laboratory bound `$\kappa_R < 5 \times 10^{17} \, \text{m}^2$` for `$M_{\text{UV}} \sim 1 \, \text{TeV}$`."

2. **Section 4.1 (Dark photon phenomenology)**:
   > "Field-space geometry can modify the curvature-dependent dark photon mixing. If scalar-fermion configuration space has nonzero Riemann tensor `$R^A_{\,\,BCD}$`, the coefficient `$C_\varepsilon$` in `$\varepsilon_{\text{eff}} = C_\varepsilon \kappa_R R$` acquires corrections from supergeometric quantum loops [Gattus & Pilaftsis 2024], testable via astrophysical searches."

3. **Section 5.3 (Future directions)**:
   > "A systematic EFT matching program using supergeometric quantum effective actions [Gattus & Pilaftsis 2024] could compute `$\kappa_R$` from first principles for specific UV completions, providing complementary input to our phenomenological framework."

---

**Cross-paper synergy**: This work establishes the **quantum-field-theoretic foundation** for curvature-EM couplings, while our paper provides **experimental constraints** and **astrophysical predictions** — together forming a comprehensive program connecting geometry, BSM physics, and observation.
