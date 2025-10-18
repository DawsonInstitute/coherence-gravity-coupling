"""
Publication-Quality Plotting Utilities

Provides consistent styling, figure helpers, and specialized plots for:
- Parameter sweeps (xi, materials)
- Optimization convergence traces
- Grid search landscape visualization

Author: GitHub Copilot (Claude Sonnet 4.5)
License: MIT
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import json


# Publication-quality style defaults
def set_publication_style():
    """
    Apply publication-quality matplotlib style settings.
    
    Features:
    - LaTeX-style fonts
    - Larger font sizes for readability
    - High-DPI output
    - Clean grid and spines
    """
    rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 14,
        'figure.dpi': 100,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.linewidth': 1.0,
        'grid.linewidth': 0.5,
        'lines.linewidth': 1.5,
        'lines.markersize': 6,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'axes.axisbelow': True,
    })


def save_figure(fig: plt.Figure, output_path: Path, formats: List[str] = ['png', 'pdf']):
    """
    Save figure in multiple formats for publication.
    
    Args:
        fig: Matplotlib figure object
        output_path: Base path (without extension)
        formats: List of file extensions to save
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    for fmt in formats:
        filepath = output_path.with_suffix(f'.{fmt}')
        fig.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"   Saved: {filepath}")


def plot_parameter_sweep(
    results: List[Dict[str, Any]],
    param_name: str,
    output_path: Optional[Path] = None,
    title: Optional[str] = None
) -> plt.Figure:
    """
    Plot results of a parameter sweep.
    
    Args:
        results: List of result dictionaries with keys: param_value, delta_tau, compute_time
        param_name: Name of swept parameter (e.g., 'xi', 'Phi0')
        output_path: Optional path to save figure
        title: Optional plot title
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    # Extract data
    param_values = [r.get('param_value', r.get(param_name, 0)) for r in results]
    delta_tau = [r['delta_tau'] for r in results]
    compute_times = [r.get('compute_time', 0) for r in results]
    
    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    
    # Plot 1: Signal vs parameter
    ax1 = axes[0]
    ax1.plot(param_values, np.abs(delta_tau), 'o-', color='#2E86AB', label='|Δτ|')
    ax1.set_xlabel(f'{param_name}')
    ax1.set_ylabel('|Δτ| [N·m]')
    ax1.set_yscale('log')
    if title:
        ax1.set_title(title)
    else:
        ax1.set_title(f'Signal Strength vs {param_name}')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Plot 2: Compute time vs parameter
    ax2 = axes[1]
    cache_hits = [r.get('cache_hit', False) for r in results]
    colors = ['#A23B72' if hit else '#F18F01' for hit in cache_hits]
    ax2.scatter(param_values, compute_times, c=colors, s=80, alpha=0.7)
    ax2.set_xlabel(f'{param_name}')
    ax2.set_ylabel('Compute Time [s]')
    ax2.set_yscale('log')
    ax2.set_title('Computation Time')
    ax2.grid(True, alpha=0.3)
    
    # Add legend for cache hits
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#F18F01', alpha=0.7, label='Computed'),
        Patch(facecolor='#A23B72', alpha=0.7, label='Cache hit')
    ]
    ax2.legend(handles=legend_elements, loc='best')
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig


def plot_material_comparison(
    results: List[Dict[str, Any]],
    output_path: Optional[Path] = None
) -> plt.Figure:
    """
    Plot comparison of different materials/coherence profiles.
    
    Args:
        results: List of result dictionaries with keys: material, delta_tau, config
        output_path: Optional path to save figure
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    # Extract data
    materials = [r['material'] for r in results]
    delta_tau = [r['delta_tau'] for r in results]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Bar plot
    x = np.arange(len(materials))
    bars = ax.bar(x, np.abs(delta_tau), color='#2E86AB', alpha=0.7, edgecolor='black')
    
    ax.set_xticks(x)
    ax.set_xticklabels(materials, rotation=45, ha='right')
    ax.set_ylabel('|Δτ| [N·m]')
    ax.set_title('Material Comparison: Signal Strength')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, delta_tau)):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height*1.1,
                f'{np.abs(val):.2e}',
                ha='center', va='bottom', fontsize=8, rotation=0)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig


def plot_optimization_trace(
    trace: List[Dict[str, Any]],
    output_path: Optional[Path] = None,
    title: Optional[str] = None
) -> plt.Figure:
    """
    Plot optimization convergence trace.
    
    Args:
        trace: List of iteration records with keys: iteration, objective, position, best_so_far
        output_path: Optional path to save figure
        title: Optional plot title
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    # Extract data
    iterations = [t['iteration'] for t in trace]
    objectives = [np.abs(t['objective']) for t in trace]
    best_so_far = [np.abs(t.get('best_so_far', t['objective'])) for t in trace]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot both traces
    ax.plot(iterations, objectives, 'o-', color='#F18F01', alpha=0.5, 
            label='Current iteration', markersize=4)
    ax.plot(iterations, best_so_far, 's-', color='#2E86AB', linewidth=2,
            label='Best so far', markersize=5)
    
    ax.set_xlabel('Iteration')
    ax.set_ylabel('|Objective| [N·m]')
    ax.set_yscale('log')
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Optimization Convergence')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Add final improvement annotation
    if len(best_so_far) > 1:
        improvement = best_so_far[-1] / best_so_far[0]
        ax.text(0.98, 0.02, f'Final improvement: {improvement:.2f}×',
                transform=ax.transAxes, ha='right', va='bottom',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig


def plot_landscape_slice(
    grid_results: Dict[str, Any],
    slice_axis: str = 'z',
    slice_value: float = 0.0,
    output_path: Optional[Path] = None
) -> plt.Figure:
    """
    Plot 2D slice of 3D optimization landscape.
    
    Args:
        grid_results: Dictionary with grid_points and objectives arrays
        slice_axis: Axis to slice ('x', 'y', or 'z')
        slice_value: Value along slice axis
        output_path: Optional path to save figure
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    points = np.array(grid_results['grid_points'])
    objectives = np.array(grid_results['objectives'])
    
    # Determine slice indices
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    slice_idx = axis_map[slice_axis]
    other_axes = [i for i in range(3) if i != slice_idx]
    
    # Filter points near slice
    mask = np.abs(points[:, slice_idx] - slice_value) < 1e-6
    sliced_points = points[mask]
    sliced_objectives = objectives[mask]
    
    if len(sliced_points) == 0:
        # Find nearest slice
        unique_vals = np.unique(points[:, slice_idx])
        slice_value = unique_vals[np.argmin(np.abs(unique_vals - slice_value))]
        mask = np.abs(points[:, slice_idx] - slice_value) < 1e-6
        sliced_points = points[mask]
        sliced_objectives = objectives[mask]
    
    # Create 2D grid
    x_vals = sliced_points[:, other_axes[0]]
    y_vals = sliced_points[:, other_axes[1]]
    z_vals = np.abs(sliced_objectives)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plot with color map
    scatter = ax.scatter(x_vals, y_vals, c=z_vals, cmap='viridis',
                         s=100, edgecolors='black', linewidths=0.5)
    
    # Find and mark optimum
    opt_idx = np.argmax(z_vals)
    ax.scatter(x_vals[opt_idx], y_vals[opt_idx], 
              marker='*', s=500, color='red', edgecolors='white', linewidths=2,
              label=f'Optimum: {z_vals[opt_idx]:.2e}')
    
    axis_names = ['x', 'y', 'z']
    ax.set_xlabel(f'{axis_names[other_axes[0]]} [m]')
    ax.set_ylabel(f'{axis_names[other_axes[1]]} [m]')
    ax.set_title(f'Signal Landscape Slice ({slice_axis}={slice_value:.3f} m)')
    
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label('|Δτ| [N·m]')
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig


def plot_landscape_3d(
    grid_results: Dict[str, Any],
    output_path: Optional[Path] = None,
    view_angles: Tuple[int, int] = (30, 45)
) -> plt.Figure:
    """
    Plot 3D scatter of optimization landscape.
    
    Args:
        grid_results: Dictionary with grid_points and objectives arrays
        output_path: Optional path to save figure
        view_angles: (elevation, azimuth) viewing angles in degrees
    
    Returns:
        Matplotlib figure object
    """
    set_publication_style()
    
    from mpl_toolkits.mplot3d import Axes3D
    
    points = np.array(grid_results['grid_points'])
    objectives = np.abs(np.array(grid_results['objectives']))
    
    # Create figure
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Color by objective value
    scatter = ax.scatter(points[:, 0], points[:, 1], points[:, 2],
                        c=objectives, cmap='viridis', s=80,
                        edgecolors='black', linewidths=0.5, alpha=0.8)
    
    # Mark optimum
    opt_idx = np.argmax(objectives)
    ax.scatter(*points[opt_idx], marker='*', s=500, color='red',
              edgecolors='white', linewidths=2,
              label=f'Optimum: {objectives[opt_idx]:.2e}')
    
    ax.set_xlabel('x [m]')
    ax.set_ylabel('y [m]')
    ax.set_zlabel('z [m]')
    ax.set_title('3D Signal Landscape')
    
    ax.view_init(elev=view_angles[0], azim=view_angles[1])
    
    cbar = plt.colorbar(scatter, ax=ax, pad=0.1, shrink=0.8)
    cbar.set_label('|Δτ| [N·m]')
    
    ax.legend()
    
    plt.tight_layout()
    
    if output_path:
        save_figure(fig, output_path)
    
    return fig


def load_and_plot_sweep(json_path: Path, output_dir: Optional[Path] = None):
    """
    Load sweep results from JSON and generate plots.
    
    Args:
        json_path: Path to sweep results JSON file
        output_dir: Directory to save plots (defaults to same as JSON)
    """
    with open(json_path) as f:
        data = json.load(f)
    
    if output_dir is None:
        output_dir = json_path.parent
    
    sweep_type = data.get('sweep_type', 'unknown')
    
    if sweep_type == 'xi':
        param_name = 'xi'
        results = data['results']
        output_path = output_dir / f"{json_path.stem}_plot"
        plot_parameter_sweep(results, param_name, output_path, 
                           title=f"Coherence Coupling Sweep (ξ)")
    elif sweep_type == 'materials':
        results = data['results']
        output_path = output_dir / f"{json_path.stem}_plot"
        plot_material_comparison(results, output_path)
    
    print(f"✅ Generated plots from {json_path}")


if __name__ == '__main__':
    # Example usage and tests
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate plots from saved results')
    parser.add_argument('json_file', type=Path, help='Path to results JSON file')
    parser.add_argument('--output-dir', type=Path, help='Output directory for plots')
    
    args = parser.parse_args()
    
    load_and_plot_sweep(args.json_file, args.output_dir)
