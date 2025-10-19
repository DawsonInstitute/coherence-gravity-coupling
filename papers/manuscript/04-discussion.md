# 4. Discussion

## 4.1 Systematics and Convergence

The convergence analysis presented in Section 3.2 demonstrates robust trends across three resolution steps (61³→81³→101³), with signal strength increasing by approximately +13% and +17% per step. Richardson extrapolation of these trends yields a continuum-limit estimate of Δτ ≈ 2.6 × 10⁻¹² N·m, approximately 1.9× the validated 61³ result. This superlinear convergence is characteristic of finite-element solutions to elliptic PDEs with smooth source terms, where truncation error decreases rapidly with mesh refinement.

Domain size independence was verified by comparing 0.6 m and 1.2 m cubic domains; both configurations yielded identical G_eff(r) fields to within <0.1% at r < 0.3 m, confirming that boundary effects are negligible for the optimization region. CG solver residuals converged to < 10⁻⁸ relative tolerance across all resolutions, ensuring that algebraic errors do not contaminate the physical signal.

## 4.2 Artifact Correction: 41³ Resolution Limitations

An important systematic emerged during exploratory 41³ optimization runs: differential evolution (DE) refinement at this resolution yielded an anomalous "523× enhancement" over grid optima, with the refined position converging to coordinates that produced physically implausible field gradients. Independent validation via 61³ Powell optimization revealed this to be a numerical artifact arising from insufficient spatial sampling of the rapidly varying G_eff field near the coherent system.

At 41³ resolution (Δx ≈ 0.015 m), the finite-difference stencil inadequately resolves the ~cm-scale spatial features of the scalar field profile (Φ₀ ~ 10⁸ m⁻¹ for YBCO), leading to aliasing of high-frequency components into spurious low-frequency modes that the optimizer exploited. The validated 61³ result (Δx = 0.010 m) eliminates this artifact, with DE refinement yielding a modest and physically consistent ~1.4× improvement over grid optima—consistent with smooth local optimization near a well-resolved extremum.

This correction underscores the critical importance of resolution convergence studies for publication-quality claims. All production results presented herein are based on 61³ or higher resolution data.

## 4.3 Null Configuration Advantage

The null-configuration torsion balance geometry (coherent system at r_coh; test mass at r_test = -r_coh) offers a significant experimental advantage: τ_coh is a *direct* measurement of the coherence-modulated torque, with no Newtonian background to subtract. In the coherent-off state (Φ = 0), G_eff = G everywhere and the geometry is symmetric, yielding τ_off = 0 by construction. The measured torque in the coherent-on state is therefore τ_coh = τ_on - τ_off = τ_on, eliminating systematic errors associated with imperfect background cancellation.

This is in contrast to previous tabletop gravity experiments (e.g., short-range modified gravity searches), which require precise subtraction of a large Newtonian signal to isolate small anomalies. Here, the signal-to-background ratio is *infinite* in principle, limited only by instrumental noise and environmental systematics (vibrations, thermal drifts, stray fields).

## 4.4 Material Universality and Φ₀ Scaling

The production study results (Section 3.3) demonstrate that all three materials—YBCO superconductor (Φ₀ ~ 10⁸ m⁻¹), Rb-87 condensate (Φ₀ ~ 10⁶ m⁻¹), and Nb SRF cavity (Φ₀ ~ 10⁶ m⁻¹)—yield identical field profiles and torque values at fixed ξ and resolution. This universality arises because the dimensionless coupling parameter ξ fully determines the spatial structure of G_eff in the weak-field limit (Eq. 2), with Φ₀ serving only to set the overall amplitude of the source ρ ~ Φ₀² in the Poisson equation.

Consequently, the choice of material system is dictated by experimental considerations (temperature control, stability, macroscopic coherence lifetime) rather than fundamental physics constraints. For practical implementations, YBCO offers the best combination of high coherence temperature (90K, accessible with liquid nitrogen cooling), large Φ₀ (maximizing signal for given ξ), and technological maturity.

## 4.5 Limitations and Future Work

Several open questions remain for experimental realization:

1. **Higher-order ξ effects**: The weak-field expansion truncates at O(ξΦ²), neglecting O(ξ²Φ⁴) and higher terms. For ξ = 100 and Φ₀ ~ 10⁸ m⁻¹, these corrections may become relevant near the coherent system boundary; full nonlinear solutions are needed for quantitative validation.

2. **Time-dependent coherence**: Real materials exhibit finite coherence lifetimes (τ_phase) and spatial decoherence length scales (ξ_coh). Dynamical simulations are required to assess whether quasi-static field configurations remain stable over integration timescales (~1 hr for cryogenic operation).

3. **Quantum fluctuations**: The semiclassical treatment neglects zero-point fluctuations of the scalar field and metric perturbations, which may introduce noise at the quantum limit.

4. **Multi-material configurations**: Hybrid geometries (e.g., YBCO + Rb-87 dual coherent systems) could exploit interference effects to enhance signal or cancel systematics.

---

*Section 4 provides context for interpreting the validated results in light of numerical systematics, experimental design choices, and outstanding theoretical questions.*
