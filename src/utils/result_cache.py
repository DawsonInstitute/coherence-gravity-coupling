#!/usr/bin/env python3
"""
Result caching for coherence-gravity simulations.

Provides fast result retrieval for repeated configurations by hashing
input parameters and storing/loading φ fields.

Author: GitHub Copilot (Claude Sonnet 4.5)
"""

import hashlib
import json
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Any
import time


class ResultCache:
    """Cache for simulation results."""
    
    def __init__(self, cache_dir: str = "results/cache"):
        """
        Initialize result cache.
        
        Args:
            cache_dir: Directory for cache storage (default: results/cache)
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._hits = 0
        self._misses = 0
        self.stats = {'hits': 0, 'misses': 0}
    
    def compute_key(
        self,
        xi: float,
        Phi0: float,
        geom_params: Dict,
        grid_resolution: int,
        domain_size: float,
        solver_method: str = 'cg',
        preconditioner: str = 'diagonal'
    ) -> str:
        """
        Compute cache key from simulation parameters.
        
        Args:
            xi: Non-minimal coupling
            Phi0: Coherence field amplitude
            geom_params: Geometry parameters dict
            grid_resolution: Grid points per dimension
            domain_size: Domain size [m]
            solver_method: Solver type
            preconditioner: Preconditioner type
        
        Returns:
            16-character hex hash
        """
        # Sort keys for consistent hashing
        cache_data = {
            'xi': float(xi),
            'Phi0': float(Phi0),
            'geom_params': {k: float(v) if isinstance(v, (int, float)) else v 
                           for k, v in sorted(geom_params.items())},
            'grid_resolution': int(grid_resolution),
            'domain_size': float(domain_size),
            'solver_method': str(solver_method),
            'preconditioner': str(preconditioner)
        }
        
        # JSON with sorted keys for reproducibility
        json_str = json.dumps(cache_data, sort_keys=True)
        hash_obj = hashlib.sha256(json_str.encode())
        return hash_obj.hexdigest()[:16]
    
    def save(
        self,
        key: str,
        result: Dict,
        phi_coherent: np.ndarray,
        phi_newtonian: np.ndarray,
        metadata: Optional[Dict] = None
    ) -> None:
        """
        Save result to cache.
        
        Args:
            key: Cache key from compute_key()
            result: Result dictionary to cache
            phi_coherent: Coherence field solution
            phi_newtonian: Newtonian field solution
            metadata: Optional additional metadata
        """
        cache_file = self.cache_dir / f"{key}.npz"
        meta_file = self.cache_dir / f"{key}.json"
        
        # Save compressed fields
        np.savez_compressed(
            cache_file,
            phi_coherent=phi_coherent,
            phi_newtonian=phi_newtonian
        )
        
        # Save metadata with result
        meta_data = {
            'timestamp': time.time(),
            'result': result,
            'custom': metadata or {}
        }
        
        with open(meta_file, 'w') as f:
            json.dump(meta_data, f, indent=2)
    
    def load(self, key: str) -> Optional[Dict]:
        """
        Load cached result.
        
        Args:
            key: Cache key from compute_key()
        
        Returns:
            Dict with 'result', 'phi_coherent', 'phi_newtonian', 'metadata' or None if not found
        """
        cache_file = self.cache_dir / f"{key}.npz"
        meta_file = self.cache_dir / f"{key}.json"
        
        if not cache_file.exists() or not meta_file.exists():
            self._misses += 1
            return None
        
        try:
            # Load compressed fields
            data = np.load(cache_file, allow_pickle=True)
            
            # Load metadata
            with open(meta_file, 'r') as f:
                metadata = json.load(f)
            
            self._hits += 1
            return {
                'result': metadata['result'],
                'phi_coherent': data['phi_coherent'],
                'phi_newtonian': data['phi_newtonian'],
                'metadata': metadata
            }
        except Exception as e:
            print(f"Warning: Failed to load cache {key}: {e}")
            self._misses += 1
            return None
    
    def clear(self):
        """Clear all cached results."""
        import shutil
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.stats = {'hits': 0, 'misses': 0}
    
    def info(self) -> Dict[str, Any]:
        """Get and print cache statistics."""
        cache_files = list(self.cache_dir.glob("*.npz"))
        total_size = sum(f.stat().st_size for f in cache_files)
        
        stats = {
            'cache_dir': str(self.cache_dir),
            'num_entries': len(cache_files),
            'total_size_mb': total_size / 1024 / 1024,
            'hits': self._hits,
            'misses': self._misses,
            'hit_rate': self._hits / (self._hits + self._misses) 
                       if (self._hits + self._misses) > 0 else 0.0
        }
        
        print(f"\n{'='*60}")
        print("CACHE STATISTICS")
        print(f"{'='*60}")
        print(f"Cache directory: {stats['cache_dir']}")
        print(f"Entries:         {stats['num_entries']}")
        print(f"Total size:      {stats['total_size_mb']:.2f} MB")
        print(f"Hits:            {stats['hits']}")
        print(f"Misses:          {stats['misses']}")
        print(f"Hit rate:        {stats['hit_rate']:.1%}")
        print(f"{'='*60}\n")
        
        return stats


# Global cache instance
_cache = ResultCache()


def get_cache() -> ResultCache:
    """Get global cache instance."""
    return _cache
