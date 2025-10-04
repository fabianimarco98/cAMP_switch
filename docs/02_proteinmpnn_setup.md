# ProteinMPNN Setup Guide

## Overview
ProteinMPNN is a neural network for protein sequence design. Given a backbone structure, it predicts amino acid sequences that are likely to fold into that structure.

## Installation

### Option 1: Using Conda (Recommended)
```bash
# Clone ProteinMPNN repository
git clone https://github.com/dauparas/ProteinMPNN.git
cd ProteinMPNN

# Create conda environment
conda create -n proteinmpnn python=3.9
conda activate proteinmpnn

# Install dependencies
pip install torch torchvision torchaudio
pip install biopython numpy
```

### Option 2: Using Existing Environment
```bash
conda activate your_env
pip install git+https://github.com/dauparas/ProteinMPNN.git
```

## Preparing Input Files

### 1. Backbone Structures
Copy selected designs from RFDiffusion:
```bash
cp 01_rfdiffusion/outputs/cAMP_binder_top*.pdb 02_proteinmpnn/inputs/
```

### 2. Chain Specification (Optional)
Create a JSON file to specify which chains to design:
```json
{
  "cAMP_binder_top1.pdb": {
    "designed_chains": ["A"],
    "fixed_chains": []
  }
}
```

## Running ProteinMPNN

### Basic Command
```bash
conda activate proteinmpnn

# Run ProteinMPNN for single structure
python /path/to/ProteinMPNN/protein_mpnn_run.py \
  --pdb_path 02_proteinmpnn/inputs/cAMP_binder_top1.pdb \
  --out_folder 02_proteinmpnn/outputs/ \
  --num_seq_per_target 10 \
  --sampling_temp 0.1
```

### Batch Processing
```bash
# Design sequences for all input structures
python /path/to/ProteinMPNN/protein_mpnn_run.py \
  --pdb_path_chains 02_proteinmpnn/inputs/chain_definitions.jsonl \
  --out_folder 02_proteinmpnn/outputs/ \
  --num_seq_per_target 10 \
  --sampling_temp 0.1 \
  --batch_size 1
```

### Advanced Options
```bash
# With specific residue positions fixed
python /path/to/ProteinMPNN/protein_mpnn_run.py \
  --pdb_path 02_proteinmpnn/inputs/cAMP_binder_top1.pdb \
  --out_folder 02_proteinmpnn/outputs/ \
  --num_seq_per_target 20 \
  --sampling_temp 0.1 \
  --fixed_positions "A10 A11 A12" \
  --omit_AAs "C" \
  --bias_AA_jsonl 02_proteinmpnn/inputs/aa_bias.jsonl
```

## Output Files

ProteinMPNN generates:
1. **FASTA files**: Designed sequences
2. **Score files**: Sequence recovery and confidence metrics
3. **PDB files**: Full-atom models with designed sequences

Example output structure:
```
02_proteinmpnn/outputs/
├── seqs/
│   ├── cAMP_binder_top1.fa
│   └── cAMP_binder_top2.fa
├── scores/
│   ├── cAMP_binder_top1_scores.csv
│   └── cAMP_binder_top2_scores.csv
└── backbones/
    └── cAMP_binder_top1.pdb
```

## Post-Processing

### 1. Sequence Analysis
```python
# Example script: 02_proteinmpnn/scripts/analyze_sequences.py
from Bio import SeqIO
import pandas as pd

# Load sequences
sequences = list(SeqIO.parse("02_proteinmpnn/outputs/seqs/cAMP_binder_top1.fa", "fasta"))

# Calculate diversity metrics
# Check for common motifs
# Filter by confidence scores
```

### 2. Building Full-Atom Models
Use a structure prediction tool to validate designs:
```bash
# Option 1: AlphaFold2
# Option 2: ESMFold
# Option 3: Simple threading onto backbone
```

### 3. Ranking Sequences
Consider:
- Sequence recovery score
- Predicted confidence (pLDDT-like metric)
- Diversity from training set
- Presence of stabilizing motifs

## Quality Control

### Metrics to Check
1. **Sequence Recovery**: How well sequences match natural distributions
2. **Hydrophobic Core**: Proper burial of hydrophobic residues
3. **Surface Properties**: Appropriate charge distribution
4. **Binding Interface**: Residues at cAMP binding site

### Filtering Criteria
```python
# Example filtering
def filter_sequences(scores_df):
    filtered = scores_df[
        (scores_df['score'] < -2.0) &  # Good confidence
        (scores_df['seq_recovery'] > 0.3) &  # Reasonable recovery
        (scores_df['global_score'] < -1.0)  # Overall quality
    ]
    return filtered
```

## Troubleshooting

### Common Issues
1. **Low sequence recovery**: Adjust sampling temperature (try 0.05-0.3)
2. **Limited diversity**: Increase number of sequences or temperature
3. **Unrealistic sequences**: Check backbone quality, consider constraints

## Next Steps

After sequence design:
1. Select top 5-10 sequence variants per backbone
2. Build full-atom models (if not already done)
3. Prepare systems for MD simulations
4. Copy to `03_md_simulations/setup/`

## References
- Dauparas et al., "Robust deep learning-based protein sequence design using ProteinMPNN", Science (2022)
- ProteinMPNN GitHub: https://github.com/dauparas/ProteinMPNN
