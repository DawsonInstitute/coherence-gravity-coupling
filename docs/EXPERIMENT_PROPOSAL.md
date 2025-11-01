# Experimental Proposal: Tabletop Constraints on Curvature–EM Coupling (R·F²)

Date: 2025-10-31

## 1. Objective
Tighten laboratory bounds on curvature–EM anomalous coupling $\kappa_R$ using a cryogenic torsion balance near a high-field magnetic region and controlled curvature proxy.

## 2. Concept of Operation
- Measure torque on a torsion balance while modulating an EM invariant region ($F^2 = 2(B^2 - E^2/c^2)$) in the presence of a known curvature proxy (effective Ricci scalar $R$ from calibrated mass distribution and geometry).
- Null experiment: detect no torque change at the modulation frequency → set upper bound on $\kappa_R$ via $\kappa_R < \delta / (R |F^2|)$.

## 3. Apparatus
- Torsion balance with fiber torsional constant $\kappa_{\text{fiber}} \sim 10^{-8}$ N·m/rad; Q > 10⁴ at 4 K.
- Test mass pair (e.g., 10–50 g tungsten) at radius 5–10 cm.
- Superconducting solenoid or permanent magnet pair achieving B = 1–10 T in the interaction region.
- Curvature proxy mass near field region to induce controlled R (modeled with finite elements; target $R \sim 10^{-26}$–$10^{-22}$ m⁻²).
- Cryostat enabling 4 K operation; low-vibration feedthroughs; μ-metal shielding.
- Readout: interferometric or capacitive angle sensor (≤ 1 nrad/√Hz).

## 4. Modulation & Readout
- On/Off or AC modulation of B-field at torsional resonance (0.1–5 Hz), phase-locked to lock-in detection.
- Synchronous demodulation at 1f; optional 2f to probe quadratic systematics.
- Integration: 0.5–24 h per configuration.

## 5. Noise & Systematics
- Thermal torque noise: $S_\tau^{1/2} \approx \sqrt{4 k_B T \kappa_{\text{fiber}} / (\omega Q)}$.
- Seismic/tilt coupling: active isolation (10–100×). Monitor tilt with dual-axis sensors.
- Magnetic cross-talk: Faraday/μ-metal shielding; null runs with field reversed.
- Eddy-current damping: choose low-loss materials; verify phase at modulation.
- Gravitational gradients: perform Newtonian nulling runs; swap masses to average.

## 6. Operating Points
- Baseline: B = 1.0 T, R = 1e-26 m⁻², precision δ = 1e-6 (fractional torque sensitivity).
- Aggressive: B = 10 T, R = 1e-22 m⁻², δ = 1e-8.
- Expected bounds (from analysis pipeline):
  - vs B (R=1e-26, δ=1e-6): $\kappa_R < 5\times10^{19}$ m² at 1 T; $5\times10^{17}$ m² at 10 T.
  - vs R (B=1 T, δ=1e-6): $\kappa_R < 5\times10^{23}\,\to\,5\times10^{15}$ m² for R from $10^{-30}$ to $10^{-22}$ m⁻².
  - vs δ (B=1 T, R=1e-26): $\kappa_R < 5\times10^{21}\,\to\,5\times10^{15}$ m² for δ from $10^{-4}$ to $10^{-10}$.

## 7. Measurement Plan
1. Commission torsion balance at room temp; measure transfer function.
2. Cool to 77 K (LN₂), then 4 K (LHe); verify Q, resonance, noise floor.
3. Run baseline modulation at B = 1 T, R = 1e-26 m⁻² for 2 h; compute bound.
4. Increase B and/or R; repeat; map $(B, R)$ plane for δ achieved.
5. Perform reversed-field and mass-swap nulls to isolate magnetic/gravitational systematics.

## 8. Deliverables
- κ-bounds table across (B, R, δ) with uncertainties (from lock-in variance and systematics).
- Open dataset: CSVs + plots + lab log → DOI (Zenodo). Link in preprint.
- Methods write-up enabling reproduction by independent labs.

## 9. Timeline & Budget (ROM)
- Design & procurement: 6–8 weeks
- Integration & commissioning: 4–6 weeks
- Data runs: 4 weeks
- Analysis & publication: 2 weeks
- Budget: $150–$300k (torsion balance, cryostat, isolation, magnet)

## 10. References
- See `papers/null_results_preprint.md` and `results/reports/` for current pipeline outputs and tables.
