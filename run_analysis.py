"""
Main runner for coherence-gravity coupling analysis.

Executes the complete parameter space scan and produces summary report.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from analysis.g_eff_scanner import run_parameter_scan

if __name__ == "__main__":
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  COHERENCE-MODULATED GRAVITY COUPLING - Phase D Analysis        ║")
    print("║                                                                  ║")
    print("║  Exploring field-dependent G_eff via ξRΦ² non-minimal coupling  ║")
    print("║  Can macroscopic coherence reduce curvature energy cost?        ║")
    print("╚══════════════════════════════════════════════════════════════════╝")
    print()
    
    grid, curves, benchmarks = run_parameter_scan()
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print("\nKey Question: Can we reach G_eff << G with realizable coherence?")
    print()
    print("Next Steps:")
    print("1. Review results/parameter_scan.json")
    print("2. Identify most promising (ξ, Φ₀) configurations")
    print("3. Assess experimental feasibility")
    print("4. If realizable: Design proof-of-concept measurement")
    print("5. If gap remains: Explore coherence amplification mechanisms")
    print()
