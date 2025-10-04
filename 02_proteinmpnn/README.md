# ProteinMPNN - Sequence Design

This directory contains files and scripts for designing amino acid sequences using ProteinMPNN.

## Directory Structure

- `inputs/` - Backbone structures from RFDiffusion
- `outputs/` - Designed sequences (PDB and FASTA files)
- `scripts/` - Scripts to run ProteinMPNN and process sequences

## Workflow

1. Prepare backbone structures:
   - Import top designs from RFDiffusion
   - Define fixed positions (if any)
   - Specify design regions

2. Run ProteinMPNN:
   - Generate multiple sequence variants per backbone
   - Optimize for binding interface residues
   - Consider structural stability

3. Post-process sequences:
   - Rank by predicted confidence
   - Filter by sequence diversity
   - Check for common motifs and stability

## Usage

See `docs/02_proteinmpnn_setup.md` for detailed instructions.

## Expected Outputs

- FASTA files with designed sequences
- PDB files with full atom models
- Sequence diversity metrics
- Confidence scores per design
