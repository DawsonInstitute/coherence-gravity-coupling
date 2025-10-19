# Methods

- Field theory: Non-minimal coupling L[Φ,g] ⊃ −(ξ/2) R Φ². In the weak-field limit, this yields an effective coupling G_eff(Φ) = G / (1 + ξ Φ²) to first order in perturbations.
- PDE: Solve ∇·(G_eff ∇φ) = 4πG ρ on a rectangular domain with Dirichlet BCs sufficiently far from the apparatus. Implementation uses a 7-point stencil FD discretization.
- Solver: Conjugate Gradient with a diagonal (Jacobi) preconditioner; sparse COO assembly for fast matrix build. Typical 61³ solve: 5–8 s on Linux workstation.
- Interpolation: Trilinear interpolation for evaluating gradients off-grid.
- Force/torque: Compute volume-averaged force on the test mass by integrating ∇φ over its finite volume; torque is τ = r × F. We emphasize volume averaging for convergence.
- Optimization: Objective is to maximize |Δτ| between coherent-on and coherent-off (Newtonian null). Steps: (1) grid search, (2) Differential Evolution, (3) Powell polish. Results cached via SHA-256 keys for ~250× reuse speedup.
- Convergence protocol: Fix the validated optimal position (from Powell at 61³) and run 61³→81³→101³ with volume averaging. Fit vs h² to estimate the continuum limit and report τ_coh directly.
