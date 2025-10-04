# Molecular Dynamics Simulations - GROMACS

This directory contains files and scripts for running MD simulations with GROMACS.

## Directory Structure

- `setup/` - System preparation (topology, solvation, etc.)
- `simulations/` - MD parameter files (MDP files) and run scripts
- `trajectories/` - Output trajectory files and checkpoints
- `scripts/` - Workflow automation scripts

## Workflow

1. System Preparation:
   - Generate topology files for protein and cAMP complex
   - Solvate the system in a water box
   - Add ions for neutralization
   - Energy minimization

2. Equilibration:
   - NVT equilibration (constant volume)
   - NPT equilibration (constant pressure)
   - Gradual release of restraints

3. Production MD:
   - Run production simulations (typically 100-500 ns)
   - Monitor stability and convergence
   - Save trajectories for analysis

4. Binding Analysis:
   - Calculate binding free energies (MM-PBSA or similar)
   - Analyze protein-ligand interactions
   - Identify stable binding modes

## Usage

See `docs/03_md_protocol.md` for detailed instructions.

## Expected Outputs

- Trajectory files (XTC format)
- Energy files (EDR format)
- Structure files at key timepoints
- Binding metrics and interaction maps
