#!/usr/bin/env python3
"""
Rank RFDiffusion designs by confidence scores and structural quality metrics.

Usage:
    python rank_designs.py --input_dir ../outputs/ --output_csv ../outputs/rankings.csv
"""

import argparse
import os
import glob
from pathlib import Path
import pandas as pd
import numpy as np


def extract_plddt_from_pdb(pdb_file):
    """Extract mean pLDDT score from PDB file (if available in B-factor column)."""
    plddts = []
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    bfactor = float(line[60:66].strip())
                    plddts.append(bfactor)
                except (ValueError, IndexError):
                    pass
    
    if plddts:
        return np.mean(plddts)
    return None


def count_residues(pdb_file):
    """Count number of residues in PDB file."""
    residues = set()
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM'):
                try:
                    res_num = int(line[22:26].strip())
                    chain = line[21]
                    residues.add((chain, res_num))
                except (ValueError, IndexError):
                    pass
    return len(residues)


def rank_designs(input_dir, output_csv, top_n=None):
    """
    Rank RFDiffusion designs by quality metrics.
    
    Args:
        input_dir: Directory containing PDB files
        output_csv: Output CSV file path
        top_n: Number of top designs to save separately (optional)
    """
    # Find all PDB files
    pdb_files = glob.glob(os.path.join(input_dir, '*.pdb'))
    
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return
    
    print(f"Found {len(pdb_files)} PDB files")
    
    # Collect metrics for each design
    results = []
    for pdb_file in pdb_files:
        design_name = Path(pdb_file).stem
        
        # Extract metrics
        mean_plddt = extract_plddt_from_pdb(pdb_file)
        num_residues = count_residues(pdb_file)
        
        results.append({
            'design_name': design_name,
            'pdb_file': pdb_file,
            'mean_plddt': mean_plddt if mean_plddt else 0,
            'num_residues': num_residues
        })
    
    # Create DataFrame and sort by pLDDT
    df = pd.DataFrame(results)
    df = df.sort_values('mean_plddt', ascending=False)
    df['rank'] = range(1, len(df) + 1)
    
    # Reorder columns
    df = df[['rank', 'design_name', 'mean_plddt', 'num_residues', 'pdb_file']]
    
    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Rankings saved to {output_csv}")
    
    # Print top designs
    print("\nTop 10 designs:")
    print(df.head(10).to_string(index=False))
    
    # Optionally save top N designs to separate directory
    if top_n:
        top_dir = os.path.join(input_dir, 'top_designs')
        os.makedirs(top_dir, exist_ok=True)
        
        for idx, row in df.head(top_n).iterrows():
            src = row['pdb_file']
            dst = os.path.join(top_dir, f"top_{row['rank']:02d}_{row['design_name']}.pdb")
            import shutil
            shutil.copy(src, dst)
        
        print(f"\nTop {top_n} designs copied to {top_dir}")


def main():
    parser = argparse.ArgumentParser(
        description='Rank RFDiffusion designs by confidence scores'
    )
    parser.add_argument(
        '--input_dir',
        type=str,
        required=True,
        help='Directory containing PDB files from RFDiffusion'
    )
    parser.add_argument(
        '--output_csv',
        type=str,
        required=True,
        help='Output CSV file for rankings'
    )
    parser.add_argument(
        '--top_n',
        type=int,
        default=None,
        help='Number of top designs to save separately'
    )
    
    args = parser.parse_args()
    rank_designs(args.input_dir, args.output_csv, args.top_n)


if __name__ == '__main__':
    main()
