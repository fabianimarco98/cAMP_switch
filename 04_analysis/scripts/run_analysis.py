#!/usr/bin/env python3
"""
Master analysis script for MD trajectories.

Runs all standard analyses and generates summary report.
"""

import os
import sys
import subprocess
from pathlib import Path


def run_command(cmd, description):
    """Run a shell command and handle errors."""
    print(f"\n=== {description} ===")
    print(f"Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=True, check=True, 
                              capture_output=True, text=True)
        print(result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error: {e}")
        print(e.stderr)
        return False


def check_files():
    """Check if required input files exist."""
    required = [
        '../03_md_simulations/setup/md.tpr',
        '../03_md_simulations/setup/md.xtc'
    ]
    
    missing = [f for f in required if not os.path.exists(f)]
    
    if missing:
        print("Error: Missing required files:")
        for f in missing:
            print(f"  - {f}")
        return False
    
    return True


def main():
    print("=== cAMP Binder MD Analysis Pipeline ===")
    
    # Check input files
    if not check_files():
        print("\nPlease ensure MD simulation has completed.")
        sys.exit(1)
    
    # Create output directories
    os.makedirs('results', exist_ok=True)
    os.makedirs('figures', exist_ok=True)
    
    # List of analysis scripts to run
    analyses = [
        ('python scripts/calculate_rmsd.py', 'RMSD Analysis'),
        ('python scripts/calculate_rmsf.py', 'RMSF Analysis'),
        ('python scripts/analyze_contacts.py', 'Contact Analysis'),
        ('python scripts/analyze_hbonds.py', 'Hydrogen Bond Analysis'),
    ]
    
    # Run each analysis
    success_count = 0
    for cmd, desc in analyses:
        if run_command(cmd, desc):
            success_count += 1
    
    print(f"\n=== Analysis Complete ===")
    print(f"Successfully completed {success_count}/{len(analyses)} analyses")
    
    # Generate report if all analyses succeeded
    if success_count == len(analyses):
        print("\nGenerating summary report...")
        if run_command('python scripts/generate_report.py', 'Report Generation'):
            print("\nFull analysis report available in results/REPORT.md")
    
    print("\nOutput files:")
    print("  - Results: 04_analysis/results/")
    print("  - Figures: 04_analysis/figures/")


if __name__ == '__main__':
    main()
