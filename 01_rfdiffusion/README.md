# RFDiffusion - Structure Generation

This directory contains files and scripts for generating protein backbone structures using RFDiffusion.

## Directory Structure

- `inputs/` - Configuration files and input parameters for RFDiffusion
- `outputs/` - Generated protein backbone structures (PDB files)
- `scripts/` - Scripts to run RFDiffusion and process outputs

## Workflow

1. Prepare input configuration specifying:
   - Target ligand (cAMP) structure
   - Binding site constraints
   - Scaffold length and architecture preferences

2. Run RFDiffusion to generate candidate binders:
   - Multiple designs with varied architectures
   - Filter by geometric complementarity to cAMP

3. Post-process outputs:
   - Rank designs by confidence scores
   - Filter by structural quality metrics
   - Select top candidates for sequence design

## Usage

See `docs/01_rfdiffusion_setup.md` for detailed instructions.

## Expected Outputs

- Backbone PDB files for top designs
- Confidence scores and metrics
- Visualization of binding poses
