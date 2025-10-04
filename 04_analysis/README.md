# Analysis and Visualization

This directory contains scripts and results for analyzing MD simulations and evaluating designs.

## Directory Structure

- `scripts/` - Analysis scripts (Python, R, or other)
- `results/` - Numerical results and data tables
- `figures/` - Generated plots and visualizations

## Analysis Components

1. **Structural Stability Analysis**
   - RMSD (Root Mean Square Deviation) over time
   - RMSF (Root Mean Square Fluctuation) per residue
   - Radius of gyration
   - Secondary structure evolution

2. **Binding Analysis**
   - Protein-cAMP interaction energies
   - Hydrogen bond analysis
   - Contact maps and interaction frequencies
   - Binding pose clustering

3. **Free Energy Calculations**
   - MM-PBSA/MM-GBSA calculations
   - Binding free energy decomposition
   - Key residue contributions

4. **Comparative Analysis**
   - Compare multiple designs
   - Rank candidates by binding affinity and stability
   - Identify design-performance relationships

## Usage

See `docs/04_analysis_guide.md` for detailed instructions.

## Expected Outputs

- Summary tables of key metrics
- Time series plots (RMSD, energy, etc.)
- Structural visualizations
- Final design rankings and recommendations
