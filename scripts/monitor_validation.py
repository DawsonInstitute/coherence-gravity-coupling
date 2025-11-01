#!/usr/bin/env python3
"""
Monitor validation runs and parse results when complete.

Checks validation_61_rb87.log and validation_61_nb.log for completion,
then extracts optimized position and delta_tau values.
"""
import time
import re
import subprocess
import sys
from pathlib import Path


def check_process_running(log_file):
    """Check if the validation process for this log is still running."""
    result = subprocess.run(
        ["pgrep", "-f", f"optimize_geometry.*{log_file[:-4]}"],
        capture_output=True,
        text=True
    )
    return bool(result.stdout.strip())


def parse_final_result(log_file):
    """Parse the final optimized result from a completed log."""
    if not Path(log_file).exists():
        return None
    
    with open(log_file) as f:
        content = f.read()
    
    # Look for the final "Optimization result" section
    result_match = re.search(
        r"Optimization result:.*?"
        r"Position: \[([-\d.e+, ]+)\].*?"
        r"Δτ = ([-\d.e+]+) N·m",
        content,
        re.DOTALL
    )
    
    if result_match:
        position_str = result_match.group(1)
        delta_tau = float(result_match.group(2))
        position = [float(x.strip()) for x in position_str.split(',')]
        return {
            'position': position,
            'delta_tau': delta_tau
        }
    
    return None


def get_log_stats(log_file):
    """Get basic statistics from log file."""
    if not Path(log_file).exists():
        return None
    
    lines = Path(log_file).read_text().splitlines()
    
    # Count solver iterations (look for "Solve time:" lines)
    solver_runs = len([l for l in lines if 'Solve time:' in l])
    
    # Get last torque value if available
    last_torque = None
    for line in reversed(lines):
        if 'Torque:' in line:
            match = re.search(r'Torque: ([-\d.e+]+) N·m', line)
            if match:
                last_torque = float(match.group(1))
                break
    
    return {
        'total_lines': len(lines),
        'solver_runs': solver_runs,
        'last_torque': last_torque
    }


def monitor_validation(interval=30, max_wait=3600):
    """Monitor validation runs until completion."""
    
    logs = {
        'Rb-87': 'validation_61_rb87.log',
        'Nb': 'validation_61_nb.log'
    }
    
    print("=" * 70)
    print("VALIDATION MONITORING")
    print("=" * 70)
    print(f"Checking every {interval} seconds (max {max_wait/60:.0f} min)\n")
    
    start_time = time.time()
    completed = {name: False for name in logs}
    
    while time.time() - start_time < max_wait:
        all_complete = True
        
        print(f"\n[{time.strftime('%H:%M:%S')}] Status check:")
        print("-" * 70)
        
        for name, log_file in logs.items():
            running = check_process_running(log_file)
            stats = get_log_stats(log_file)
            
            if running:
                all_complete = False
                status = f"🔄 RUNNING"
                if stats:
                    status += f" | {stats['solver_runs']} solves"
                    if stats['last_torque']:
                        status += f" | τ={stats['last_torque']:.3e} N·m"
            else:
                if not completed[name]:
                    completed[name] = True
                    print(f"\n✅ {name} validation COMPLETED!")
                status = "✅ COMPLETE"
            
            print(f"  {name:8s}: {status}")
        
        if all_complete:
            print("\n" + "=" * 70)
            print("ALL VALIDATIONS COMPLETE!")
            print("=" * 70 + "\n")
            break
        
        time.sleep(interval)
    
    # Parse final results
    print("\nFINAL RESULTS:")
    print("=" * 70)
    
    for name, log_file in logs.items():
        result = parse_final_result(log_file)
        if result:
            print(f"\n{name} Validation:")
            print(f"  Position: [{result['position'][0]:.6f}, "
                  f"{result['position'][1]:.6f}, {result['position'][2]:.6f}] m")
            print(f"  Δτ = {result['delta_tau']:.6e} N·m")
        else:
            print(f"\n{name}: Could not parse result (check {log_file})")
    
    # Compare with 41³ DE result
    print("\n" + "-" * 70)
    print("COMPARISON WITH 41³ DE REFINEMENT:")
    print("-" * 70)
    print("41³ result: Position ~[0.0234, 0.0216, -0.0224] m")
    print("41³ result: Δτ ~ -6.174e-10 N·m (for Rb-87 & Nb)")
    print("\nEnhancement factor: 523× over 41³ grid search optimum")
    print("\nValidation objective: Confirm if 61³ reproduces similar Δτ")
    print("  - If yes: Geometric resonance effect is REAL")
    print("  - If no:  Numerical artifact in 41³ DE optimization")
    print("=" * 70)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Monitor validation runs')
    parser.add_argument('--interval', type=int, default=30,
                       help='Check interval in seconds (default: 30)')
    parser.add_argument('--max-wait', type=int, default=3600,
                       help='Maximum wait time in seconds (default: 3600)')
    
    args = parser.parse_args()
    
    try:
        monitor_validation(interval=args.interval, max_wait=args.max_wait)
    except KeyboardInterrupt:
        print("\n\nMonitoring interrupted by user.")
        sys.exit(0)
