#!/usr/bin/env python3
"""Quick test of result caching functionality."""

import time
from examples.geometric_cavendish import run_geometric_cavendish
from src.utils.result_cache import get_cache

# Clear cache first
cache = get_cache()
cache.clear()
print("🧹 Cache cleared\n")

# Test parameters
test_params = {
    'xi': 100.0,
    'Phi0': 1e8,
    'grid_resolution': 41,
    'domain_size': 0.6,
    'verbose': True,
    'solver_method': 'cg',
    'preconditioner': 'diagonal'
}

print("="*70)
print("TEST 1: First run (MISS - should compute)")
print("="*70)
t0 = time.time()
result1 = run_geometric_cavendish(**test_params, cache=True)
t1 = time.time()
print(f"\n⏱️  First run time: {t1-t0:.2f} s\n")

print("="*70)
print("TEST 2: Second run (HIT - should load from cache)")
print("="*70)
t0 = time.time()
result2 = run_geometric_cavendish(**test_params, cache=True)
t2 = time.time()
print(f"\n⏱️  Second run time: {t2-t0:.2f} s")
print(f"⚡ Speedup: {(t1-t0)/(t2-t0):.1f}×\n")

# Verify results match
print("="*70)
print("VERIFICATION")
print("="*70)
print(f"Results match: {result1 == result2}")
print(f"delta_tau_1:   {result1['delta_tau']:.6e}")
print(f"delta_tau_2:   {result2['delta_tau']:.6e}")

# Show cache stats
print("\n" + "="*70)
print("CACHE STATISTICS")
print("="*70)
cache.info()
