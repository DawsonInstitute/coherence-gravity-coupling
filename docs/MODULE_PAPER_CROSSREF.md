# Module-Paper Cross-Reference Map

**Purpose**: Track which code modules address which research paper tasks

Last Updated: 2025-01-XX

---

## arXiv:2507.02362 (Bahamonde et al. - Torsion-EM Coupling)

| Task | Module | Status | Notes |
|------|--------|--------|-------|
| 25 | `src/field_equations/torsion_dof.py` | ✅ | Coherence gradients as torsion proxy |
| 26 | `src/analysis/kappa_k3_mapping.py` | ✅ | κ_R → k_3 EFT mapping (3 scenarios) |
| 27 | `src/field_equations/torsion_dof.py::duality_breaking_observable()` | 🚧 | Module ready, EFQS integration pending |
| 30 | `extreme-field-qed-simulator/src/efqs/coherence_gravity_efqs_integration.py` | ✅ | Torsion-proxy stress in EFQS pipeline |

**Key Results**:
- Conservative k_3 bound: < 10^19 m² (from lab κ_R < 5×10^17 m²)
- Magnetar amplification: → 10^-19 m² (38 orders of magnitude)
- Duality-breaking: ∫ E·(∇×A_torsion) ≠ 0 for asymmetric Φ

---

## arXiv:2508.13820 (Karimabadi et al. - NC-Schwarzschild QNMs)

| Task | Module | Status | Notes |
|------|--------|--------|-------|
| 32 | `src/analysis/qnm_laboratory.py` | ✅ | BEC trap potential matching |
| 33 | `src/analysis/karimabadi_qnm.py` | ✅ | Fig. 3 reproduction: δω/ω ∝ θ |
| 35 | `src/analysis/wkb_qnm_diagnostics.py` | ✅ | WKB extraction from V_eff |

**Key Results**:
- QNM shift scaling: δω/ω ≈ 5.15e-2 × θ (linear for small θ)
- WKB vs time-domain: ~10% frequency agreement (damping needs work)
- Laboratory QNM matching: BEC trap fidelity ~0.7 for r ∈ [2.5M, 10M]

---

## arXiv:2509.06815 (Hell & Lüst - Bumblebee FLRW)

| Task | Module | Status | Notes |
|------|--------|--------|-------|
| Cross-validate | `src/analysis/hell_lust_flrw.py` | ✅ | Fig. 2 mode growth rates |

**Key Results**:
- Standard FLRW: λ ≈ 0.000-0.005 (unexpected - needs review)
- Lorentz violation (β=0.1): Similar range (perturbative effect)
- Mode detail (k=0.1 Mpc⁻¹): Final δ(a=1) = 1.0e-5, λ = 0.000

---

## arXiv:2409.04647 (Gorkavenko et al. - Robin BC Vacuum Energy)

| Task | Module | Status | Notes |
|------|--------|--------|-------|
| 27 | `src/solvers/robin_bc_poisson.py` | ✅ | Robin BC solver with θ parameter |
| 28 | `src/solvers/robin_bc_poisson.py::compute_xi_dependent_energy()` | ✅ | E_ξ ∝ (1/4 - ξ) validated |
| 29-30 | `scripts/robin_bc_parametric_sweep.py` | 🚧 | Script created, import issues |
| 33 | `src/solvers/robin_bc_poisson.py::validate_robin_bc_limits()` | ✅ | Dirichlet/Neumann limits checked |
| 34 | Pending | ⏳ | Astrophysical boundary condition mapping |

**Key Results**:
- θ=0 (Dirichlet): E_ξ(0) = 2.97e-7 J, E_ξ(0.25) = 0, E_ξ(0.5) = -2.97e-7 J
- Scaling: E_ξ ∝ (1/4 - ξ) with slope 2.97e-7 J
- Robin interpolation: θ ∈ [0, -π/2] → Dirichlet ↔ Neumann

---

## Null Results Paper (Zenodo 17504852)

| Experimental Constraint | Module | Value | Notes |
|-------------------------|--------|-------|-------|
| κ_R upper bound | `src/analysis/kappa_k3_mapping.py` | < 5×10^17 m² | From B=10 T, R=10^-26 m^-2 nulls |
| Mapping to k_3 | `src/analysis/kappa_k3_mapping.py` | < 10^19 m² | Conservative EFT scenario |
| Astrophysical recast | `coherence_gravity_coupling/notebooks/astrophysical_recast_qnm.ipynb` | 10^38× improvement | Magnetar B~10^15 G |

---

## EFQS Integration Summary

### Core EFQS Modules
- `extreme-field-qed-simulator/src/efqs/gravitational_coupling.py`: Quadrupole pipeline, TT projection, frequency diagnostics
- `extreme-field-qed-simulator/src/efqs/source_geometries.py`: EM field sources (fixed divide-by-zero bugs)
- `extreme-field-qed-simulator/scripts/run_experiments.py`: Experimental pipeline (added frequency output)

### Integration Layer
- `extreme-field-qed-simulator/src/efqs/coherence_gravity_efqs_integration.py`: Hooks for:
  - Torsion-proxy stress addition
  - κ_R → k_3 constraint mapping  
  - Duality-breaking observable evaluation

### Status
- ✅ Module interfaces defined
- ✅ Standalone testing complete
- 🚧 Full EFQS workflow integration pending (requires YAML config examples)

---

## Module Dependency Graph

```
EFQS Pipeline
    ├── source_geometries.py (EM fields)
    │   └── E(x,t), B(x,t)
    ├── gravitational_coupling.py
    │   ├── stress_energy_from_fields(E,B) → T00
    │   ├── quadrupole_moment(T00) → Q_ij(t)
    │   ├── strain_far_field(Q_ij) → h_TT
    │   └── dominant_frequency(h_TT) → spectrum
    └── coherence_gravity_efqs_integration.py
        ├── add_torsion_proxy_stress(T00, Φ)
        │   └── → torsion_dof.py
        ├── compute_k3_constraints(κ_R)
        │   └── → kappa_k3_mapping.py
        └── evaluate_duality_breaking(E, B, Φ)
            └── → torsion_dof.py

Analysis Modules (Standalone)
    ├── hell_lust_flrw.py (FLRW mode validation)
    ├── karimabadi_qnm.py (NC-QNM cross-check)
    ├── qnm_laboratory.py (tabletop BEC matching)
    ├── wkb_qnm_diagnostics.py (V_eff → ω extraction)
    └── robin_bc_poisson.py (boundary amplification)
        └── → robin_bc_parametric_sweep.py
```

---

## Next Steps

### High Priority
1. ⏳ Create YAML config examples using integration hooks
2. ⏳ Run full EFQS pipeline with torsion/duality/k_3 outputs
3. ⏳ Fix `robin_bc_parametric_sweep.py` import issues (sympy dependency)
4. ⏳ Add extended DOF predictions (axion-photon, dark photon) to kappa_k3_mapping

### Medium Priority
5. ⏳ Implement EFQS health checks (energy conservation, mode orthogonality)
6. ⏳ Materials feasibility study for Robin BC (θ → dielectric properties)
7. ⏳ Gorkavenko Table 1 numerical reproduction

### Documentation
8. ✅ Frame conventions guide (`docs/conventions/frames.md`)
9. ⏳ README roadmap section
10. ⏳ API documentation for integration module

---

## Validation Checklist

- [x] Torsion-DOF module: Antisymmetric norm = 0 for symmetric Φ ✓
- [x] κ_R→k_3 mapping: Three EFT scenarios output reasonable bounds ✓
- [x] Robin BC solver: Validates Dirichlet/Neumann limits ✓
- [x] Robin BC solver: E_ξ ∝ (1/4 - ξ) scaling confirmed ✓
- [x] Karimabadi QNM: Linear θ-dependence δω/ω = 5.24e-2 × θ ✓
- [x] Hell & Lüst FLRW: Module functional (growth rates need review) ⚠️
- [x] WKB QNM: Functional (time-domain accuracy ~10%, needs refinement) ⚠️
- [x] EFQS integration: All three hooks tested standalone ✓
- [ ] EFQS integration: Full pipeline run with real experiment ⏳
- [ ] Robin BC sweep: Contour plots generated ⏳

---

## Contact & Contribution

**Repository**: `coherence-gravity-coupling/` (main), `extreme-field-qed-simulator/` (EFQS)  
**Documentation**: `docs/EFQS_NEXT_STEPS.md`, `docs/conventions/frames.md`  
**Issues**: Open GitHub issue or see `README.md`

---

*Last updated: See git log for this file*
