# RFDiffusion Setup Guide

## Overview
RFDiffusion is a method for structure generation of proteins with defined shapes and functions. For this project, we use it to generate protein scaffolds that can bind cAMP.

## Installation

### Option 1: Using Conda (Recommended)
```bash
# Clone RFDiffusion repository
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion

# Create conda environment
conda env create -f env/SE3nv.yml
conda activate SE3nv

# Install additional dependencies
pip install -e .

# Download model weights
mkdir models
cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
```

### Option 2: Using Docker
```bash
docker pull ghcr.io/rosettacommons/rfdiffusion:latest
```

## Preparing Input Files

### 1. cAMP Structure
- Obtain cAMP structure from PDB or generate using molecular modeling tools
- Save as PDB file in `data/cAMP.pdb`
- Ensure proper atom naming and formatting

### 2. Configuration File
Create a configuration file specifying:
- Ligand coordinates
- Binding site definition
- Scaffold constraints (length, secondary structure preferences)
- Number of designs to generate

Example configuration (`01_rfdiffusion/inputs/config.yaml`):
```yaml
# RFDiffusion configuration for cAMP binder design
inference:
  output_prefix: "01_rfdiffusion/outputs/cAMP_binder"
  num_designs: 100
  
contigmap:
  contigs: ["A1-100"]  # Design 100-residue protein
  
ligand:
  ligand_path: "data/cAMP.pdb"
  ligand_chain: "X"
  
ppi:
  hotspot_res: []  # Define if specific residues should be near ligand
```

## Running RFDiffusion

### Basic Command
```bash
conda activate SE3nv

# Run RFDiffusion for cAMP binder design
python /path/to/RFdiffusion/scripts/run_inference.py \
  --config_path 01_rfdiffusion/inputs/config.yaml \
  --output_dir 01_rfdiffusion/outputs/
```

### Advanced Options
```bash
# With specific architectural constraints
python /path/to/RFdiffusion/scripts/run_inference.py \
  inference.output_prefix=01_rfdiffusion/outputs/cAMP_binder \
  inference.num_designs=100 \
  'contigmap.contigs=[A10-40/0 A10-40]' \
  'ppi.hotspot_res=[A30,A31,A32]' \
  inference.input_pdb=data/cAMP.pdb
```

## Output Processing

### 1. Extracting Designs
Generated designs will be saved as:
- `cAMP_binder_0.pdb`, `cAMP_binder_1.pdb`, etc.

### 2. Ranking and Filtering
```bash
# Sort by confidence score (usually in PDB header)
python 01_rfdiffusion/scripts/rank_designs.py \
  --input_dir 01_rfdiffusion/outputs/ \
  --output_csv 01_rfdiffusion/outputs/rankings.csv
```

### 3. Visual Inspection
Use PyMOL or ChimeraX to visually inspect top designs:
```bash
pymol 01_rfdiffusion/outputs/cAMP_binder_top10.pdb data/cAMP.pdb
```

## Quality Metrics

Evaluate designs based on:
1. **Confidence Score**: Model's confidence in the design
2. **pLDDT**: Per-residue confidence metric
3. **Binding Pocket Geometry**: Shape complementarity to cAMP
4. **Designability**: Predicted stability of the scaffold

## Troubleshooting

### Common Issues
1. **CUDA out of memory**: Reduce batch size or scaffold length
2. **Poor binding geometry**: Adjust hotspot residues or contig constraints
3. **Low diversity**: Increase number of designs or adjust temperature

## Next Steps

After generating and filtering designs:
1. Select top 10-20 designs for sequence optimization
2. Copy selected PDB files to `02_proteinmpnn/inputs/`
3. Proceed to ProteinMPNN stage

## References
- Watson et al., "De novo design of protein structure and function with RFdiffusion", Nature (2023)
- RFDiffusion GitHub: https://github.com/RosettaCommons/RFdiffusion
