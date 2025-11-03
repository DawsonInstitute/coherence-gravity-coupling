#!/usr/bin/env python3
"""
DOF Mode Selector for Power-Law Curvature Models

Implements DOF (degrees of freedom) mode selection for theories with power-law
curvature terms R^ℓ σ^n R^m, following Hell & Lüst (2025) classification.

Based on: Hell & Lüst, arXiv:2509.20217 (2025)
Context: arxiv.2509.20217.md:34 actionable follow-up

Key Features:
- User-specified (ℓ, m, n) parameters for power-law curvature models
- DOF structure classification (number of propagating scalars)
- Singular point detection and warnings
- Frame-dependence analysis (Jordan vs Einstein)
- Decoupling diagnostics for small-R regimes

Mathematical Framework:
General action: S = ∫ d⁴x √(-g) [f(R, σ, ∇σ)]
Power-law case: f ∝ R^ℓ + σ^n R^m
DOF count depends on (ℓ, m, n) and background R, σ values
"""

from __future__ import annotations
import numpy as np
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass
from enum import Enum
import logging

logger = logging.getLogger(__name__)


class Frame(Enum):
    """Reference frame for DOF analysis."""
    JORDAN = "Jordan"
    EINSTEIN = "Einstein"


class DOFStructure(Enum):
    """Degrees of freedom structure classification."""
    ZERO_SCALARS = "0 propagating scalars"
    ONE_SCALAR = "1 propagating scalar (standard)"
    TWO_SCALARS = "2 propagating scalars (extended)"
    STRONGLY_COUPLED = "Strongly coupled (IR breakdown)"
    SINGULAR = "Singular point (DOF undefined)"


@dataclass
class PowerLawParameters:
    """Parameters for power-law curvature model R^ℓ σ^n R^m."""
    ell: int  # Pure R power (R^ℓ term)
    m: int    # Mixed R power (σ^n R^m term)
    n: int    # Scalar power (σ^n R^m term)
    
    def __post_init__(self):
        """Validate parameters."""
        if self.ell < 0 or self.m < 0 or self.n < 0:
            raise ValueError("Power-law exponents must be non-negative")
    
    def __repr__(self) -> str:
        return f"(ℓ={self.ell}, m={self.m}, n={self.n})"


@dataclass
class BackgroundState:
    """Background curvature and field values."""
    R: float  # Ricci scalar [m^-2] or dimensionless
    sigma: float  # Coherence field value (dimensionless)
    
    def is_near_singular(self, tolerance: float = 1e-10) -> bool:
        """Check if background is near singular point (R→0 or σ→0)."""
        return abs(self.R) < tolerance or abs(self.sigma) < tolerance


@dataclass
class DOFDiagnostics:
    """DOF analysis diagnostics."""
    dof_structure: DOFStructure
    frame: Frame
    num_scalars: int
    is_singular: bool
    is_decoupled: bool
    warnings: List[str]
    recommendations: List[str]
    
    # Technical details
    kinetic_matrix_det: Optional[float] = None
    effective_mass_squared: Optional[float] = None
    coupling_strength: Optional[float] = None


class DOFModeSelector:
    """
    DOF mode selector for power-law curvature theories.
    
    Implements classification scheme from Hell & Lüst (2025) to identify
    number and nature of propagating scalar modes based on model parameters
    and background state.
    """
    
    def __init__(
        self,
        params: PowerLawParameters,
        frame: Frame = Frame.JORDAN
    ):
        """
        Initialize DOF mode selector.
        
        Args:
            params: Power-law model parameters (ℓ, m, n)
            frame: Reference frame for analysis (Jordan or Einstein)
        """
        self.params = params
        self.frame = frame
        self.diagnostic_cache: Dict[Tuple, DOFDiagnostics] = {}
        
    def analyze_dof_structure(
        self,
        background: BackgroundState,
        verbose: bool = True
    ) -> DOFDiagnostics:
        """
        Analyze DOF structure for given background state.
        
        Args:
            background: Background curvature and field values
            verbose: Print warnings and recommendations
            
        Returns:
            DOF diagnostics with classification and warnings
        """
        # Check cache
        cache_key = (self.params.ell, self.params.m, self.params.n,
                     background.R, background.sigma, self.frame)
        if cache_key in self.diagnostic_cache:
            return self.diagnostic_cache[cache_key]
        
        warnings_list = []
        recommendations = []
        is_singular = False
        is_decoupled = False
        
        # Check for singular point
        if background.is_near_singular():
            warnings_list.append(
                "Background near singular point (R→0 or σ→0). "
                "DOF structure may change discontinuously."
            )
            is_singular = True
        
        # Classify DOF structure based on (ℓ, m, n)
        dof_structure, num_scalars = self._classify_dof(background)
        
        # Check for decoupling in small-R regime
        if abs(background.R) < 1e-26:  # Laboratory R ~ 10^-26 m^-2
            is_decoupled = self._check_small_r_decoupling(background)
            if is_decoupled:
                warnings_list.append(
                    f"Small-R regime (R={background.R:.2e} m⁻²). "
                    "Scalar modes may decouple or become non-dynamical."
                )
                recommendations.append(
                    "Consider increasing curvature scale or using "
                    "astrophysical background for observable effects."
                )
        
        # Frame-dependent checks
        if self.frame == Frame.JORDAN:
            self._check_jordan_frame_issues(background, warnings_list, recommendations)
        else:
            self._check_einstein_frame_issues(background, warnings_list, recommendations)
        
        # Compute kinetic matrix determinant (simplified)
        kinetic_det = self._compute_kinetic_determinant(background)
        
        # Estimate effective mass²
        m_eff_sq = self._estimate_effective_mass_squared(background)
        
        # Coupling strength
        coupling = self._estimate_coupling_strength(background)
        
        diagnostics = DOFDiagnostics(
            dof_structure=dof_structure,
            frame=self.frame,
            num_scalars=num_scalars,
            is_singular=is_singular,
            is_decoupled=is_decoupled,
            warnings=warnings_list,
            recommendations=recommendations,
            kinetic_matrix_det=kinetic_det,
            effective_mass_squared=m_eff_sq,
            coupling_strength=coupling
        )
        
        # Cache result
        self.diagnostic_cache[cache_key] = diagnostics
        
        if verbose:
            self._print_diagnostics(diagnostics)
        
        return diagnostics
    
    def _classify_dof(
        self,
        background: BackgroundState
    ) -> Tuple[DOFStructure, int]:
        """
        Classify DOF structure based on (ℓ, m, n) and background.
        
        Implements simplified version of Hell & Lüst classification:
        - (1, 0, 0): f(R) theory → 0 or 1 scalar depending on R
        - (2, 1, 1): f(R) + σR → 1 scalar (standard)
        - (1, 2, 2): Extended models → 2 scalars possible
        """
        ell, m, n = self.params.ell, self.params.m, self.params.n
        R, sigma = background.R, background.sigma
        
        # Pure f(R) theories (n = 0)
        if n == 0 and m == 0:
            if ell == 1:
                # GR: no extra DOF
                return DOFStructure.ZERO_SCALARS, 0
            elif ell == 2 and abs(R) < 1e-30:
                # f(R) = R² near R=0: mode decouples
                return DOFStructure.ZERO_SCALARS, 0
            else:
                # Generic f(R): 1 scalar (scalaron)
                return DOFStructure.ONE_SCALAR, 1
        
        # Non-minimal coupling σ^n R^m
        if ell == 1 and m == 1 and n == 2:
            # Standard ξRΦ² (our model)
            if abs(R) > 1e-30 and abs(sigma) > 1e-10:
                return DOFStructure.ONE_SCALAR, 1
            else:
                # Near singular point: decoupling
                return DOFStructure.ZERO_SCALARS, 0
        
        # Extended models with multiple scalars
        if n >= 2 and m >= 2:
            if abs(R * sigma**n) > 1e-20:
                return DOFStructure.TWO_SCALARS, 2
            else:
                return DOFStructure.ZERO_SCALARS, 0
        
        # Strong coupling regime (IR breakdown)
        if m >= 3 or n >= 3:
            return DOFStructure.STRONGLY_COUPLED, 1
        
        # Default: 1 scalar
        return DOFStructure.ONE_SCALAR, 1
    
    def _check_small_r_decoupling(
        self,
        background: BackgroundState
    ) -> bool:
        """Check if scalar modes decouple in small-R limit."""
        ell, m = self.params.ell, self.params.m
        R = background.R
        
        # For R^ℓ theories: mode decouples if R → 0 and ℓ ≥ 2
        if ell >= 2 and abs(R) < 1e-26:
            return True
        
        # For σ^n R^m: mode decouples if R^m → 0
        if m >= 1 and abs(R**m) < 1e-30:
            return True
        
        return False
    
    def _check_jordan_frame_issues(
        self,
        background: BackgroundState,
        warnings: List[str],
        recommendations: List[str]
    ) -> None:
        """Check for Jordan frame specific issues."""
        # Singular points in Jordan frame
        if self.params.ell >= 2 and abs(background.R) < 1e-28:
            warnings.append(
                "Jordan frame singular at R→0 for ℓ≥2. "
                "Consider Einstein frame analysis."
            )
    
    def _check_einstein_frame_issues(
        self,
        background: BackgroundState,
        warnings: List[str],
        recommendations: List[str]
    ) -> None:
        """Check for Einstein frame specific issues."""
        # Conformal transformation singularities
        if self.params.n >= 2 and abs(background.sigma) < 1e-10:
            warnings.append(
                "Einstein frame may be singular at σ→0 for n≥2. "
                "Verify conformal factor regularity."
            )
    
    def _compute_kinetic_determinant(
        self,
        background: BackgroundState
    ) -> Optional[float]:
        """
        Compute determinant of kinetic matrix (simplified).
        
        For 2-scalar system: det(K) determines propagation.
        det(K) < 0 → ghost, det(K) = 0 → strong coupling
        """
        # Simplified: use R and σ as proxies
        R, sigma = background.R, background.sigma
        ell, m, n = self.params.ell, self.params.m, self.params.n
        
        # Kinetic matrix approximately:
        # K ~ [[∂²f/∂R², ∂²f/∂R∂σ], [∂²f/∂σ∂R, ∂²f/∂σ²]]
        # For f ~ R^ℓ + σ^n R^m:
        K11 = ell * (ell - 1) * R**(ell - 2) if ell >= 2 else 0.0
        K12 = n * m * sigma**(n-1) * R**(m-1) if n >= 1 and m >= 1 else 0.0
        K22 = n * (n - 1) * sigma**(n-2) * R**m if n >= 2 else 0.0
        
        det_K = K11 * K22 - K12**2
        
        return float(det_K)
    
    def _estimate_effective_mass_squared(
        self,
        background: BackgroundState
    ) -> Optional[float]:
        """
        Estimate effective scalar mass² (simplified).
        
        For non-minimal coupling: m²_eff ~ 2ξ|R|
        """
        R = background.R
        m, n = self.params.m, self.params.n
        
        if m == 1 and n == 2:
            # Standard ξRΦ² case
            xi_effective = 50.0  # Example coupling
            m_eff_sq = 2 * xi_effective * abs(R)
        else:
            # Generic estimate
            m_eff_sq = abs(R)
        
        return float(m_eff_sq)
    
    def _estimate_coupling_strength(
        self,
        background: BackgroundState
    ) -> Optional[float]:
        """Estimate dimensionless coupling strength."""
        R, sigma = background.R, background.sigma
        n, m = self.params.n, self.params.m
        
        # Coupling ~ σ^n R^m (dimensionless)
        coupling = abs(sigma**n * R**m) if n > 0 and m > 0 else 0.0
        
        return float(coupling)
    
    def _print_diagnostics(self, diagnostics: DOFDiagnostics) -> None:
        """Print diagnostic information."""
        print("\n" + "=" * 70)
        print(f"DOF Mode Analysis: {self.params}")
        print(f"Frame: {self.frame.value}")
        print("=" * 70)
        
        print(f"\nDOF Structure: {diagnostics.dof_structure.value}")
        print(f"Number of propagating scalars: {diagnostics.num_scalars}")
        
        if diagnostics.is_singular:
            print("\n⚠️  SINGULAR POINT detected")
        if diagnostics.is_decoupled:
            print("\n⚠️  DECOUPLING regime detected")
        
        if diagnostics.kinetic_matrix_det is not None:
            print(f"\nKinetic matrix determinant: {diagnostics.kinetic_matrix_det:.2e}")
            if diagnostics.kinetic_matrix_det < 0:
                print("  ⚠️  Negative determinant → Ghost instability")
            elif abs(diagnostics.kinetic_matrix_det) < 1e-20:
                print("  ⚠️  Near-zero determinant → Strong coupling")
        
        if diagnostics.effective_mass_squared is not None:
            print(f"Effective mass²: {diagnostics.effective_mass_squared:.2e}")
        
        if diagnostics.coupling_strength is not None:
            print(f"Coupling strength: {diagnostics.coupling_strength:.2e}")
        
        if diagnostics.warnings:
            print("\n⚠️  WARNINGS:")
            for i, warning in enumerate(diagnostics.warnings, 1):
                print(f"  {i}. {warning}")
        
        if diagnostics.recommendations:
            print("\n💡 RECOMMENDATIONS:")
            for i, rec in enumerate(diagnostics.recommendations, 1):
                print(f"  {i}. {rec}")
        
        print("=" * 70 + "\n")


def demo_dof_mode_selector():
    """Demonstrate DOF mode selector for various model parameters."""
    print("=" * 70)
    print("DOF Mode Selector Demonstration")
    print("=" * 70)
    
    # Test cases
    test_cases = [
        # (ℓ, m, n), R, σ, description
        (PowerLawParameters(1, 0, 0), 1e-26, 1.0, "GR (baseline)"),
        (PowerLawParameters(2, 0, 0), 1e-26, 1.0, "f(R)=R² small-R"),
        (PowerLawParameters(1, 1, 2), 1e-26, 1.0, "ξRΦ² (our model) lab"),
        (PowerLawParameters(1, 1, 2), 1e-6, 1.0, "ξRΦ² magnetar"),
        (PowerLawParameters(2, 2, 2), 1e-26, 0.5, "Extended 2-scalar"),
        (PowerLawParameters(1, 1, 2), 1e-35, 1.0, "Near R→0 singular"),
    ]
    
    for params, R, sigma, description in test_cases:
        print(f"\nTest: {description}")
        print("-" * 70)
        
        background = BackgroundState(R=R, sigma=sigma)
        selector = DOFModeSelector(params, frame=Frame.JORDAN)
        diagnostics = selector.analyze_dof_structure(background, verbose=True)
    
    print("\n" + "=" * 70)
    print("Demo complete.")
    print("=" * 70)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    demo_dof_mode_selector()
