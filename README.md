# cAMP_switch

## Project Overview
This project aims to de novo design a protein binder for cyclic adenosine monophosphate (cAMP) using computational protein design methods. The workflow integrates several state-of-the-art tools for protein structure generation, sequence optimization, molecular dynamics simulations, and analysis.

## Workflow

The design pipeline consists of four main stages:

### 1. Structure Generation with RFDiffusion
- Generate initial protein backbone structures capable of binding cAMP
- RFDiffusion uses a diffusion-based generative model to create novel protein structures
- Output: Backbone coordinates (PDB files)

### 2. Sequence Design with ProteinMPNN
- Fill the generated backbones with optimal amino acid sequences
- ProteinMPNN predicts sequences that are compatible with the designed structures
- Output: Fully sequenced protein structures (PDB/FASTA files)

### 3. Molecular Dynamics Simulations with GROMACS
- Validate structural stability and binding interactions
- Run MD simulations to assess dynamic behavior
- Evaluate binding affinity and conformational changes
- Output: Trajectory files, energy profiles, binding metrics

### 4. Final Analysis
- Analyze MD trajectories
- Calculate binding free energies
- Evaluate structural stability (RMSD, RMSF)
- Identify key residues for binding
- Generate visualizations and reports

## Project Structure

```
cAMP_switch/
├── 01_rfdiffusion/          # RFDiffusion structure generation
│   ├── inputs/              # Input files and configurations
│   ├── outputs/             # Generated structures
│   └── scripts/             # RFDiffusion scripts
├── 02_proteinmpnn/          # ProteinMPNN sequence design
│   ├── inputs/              # Backbone structures from RFDiffusion
│   ├── outputs/             # Designed sequences
│   └── scripts/             # ProteinMPNN scripts
├── 03_md_simulations/       # GROMACS MD simulations
│   ├── setup/               # System preparation scripts
│   ├── simulations/         # MD run files and parameters
│   ├── trajectories/        # Output trajectories
│   └── scripts/             # GROMACS workflow scripts
├── 04_analysis/             # Analysis and visualization
│   ├── scripts/             # Analysis scripts
│   ├── results/             # Analysis outputs
│   └── figures/             # Generated plots and visualizations
├── data/                    # Shared data (cAMP structure, etc.)
└── docs/                    # Documentation
```

## Getting Started

### Prerequisites
- RFDiffusion (https://github.com/RosettaCommons/RFdiffusion)
- ProteinMPNN (https://github.com/dauparas/ProteinMPNN)
- GROMACS (version 2020 or later)
- Python 3.8+
- Required Python packages (see requirements.txt)

### Installation
Detailed installation instructions for each component can be found in the `docs/` directory.

## Usage

Each stage of the pipeline is designed to be run sequentially:

1. **RFDiffusion**: Generate candidate binder structures
   ```bash
   cd 01_rfdiffusion/scripts
   # Run RFDiffusion scripts (see docs for details)
   ```

2. **ProteinMPNN**: Design sequences for the generated structures
   ```bash
   cd 02_proteinmpnn/scripts
   # Run ProteinMPNN scripts (see docs for details)
   ```

3. **MD Simulations**: Validate designs through molecular dynamics
   ```bash
   cd 03_md_simulations/scripts
   # Run GROMACS simulation workflow (see docs for details)
   ```

4. **Analysis**: Analyze results and identify best candidates
   ```bash
   cd 04_analysis/scripts
   # Run analysis scripts (see docs for details)
   ```

## Documentation

Detailed documentation for each stage can be found in:
- [RFDiffusion Setup](docs/01_rfdiffusion_setup.md)
- [ProteinMPNN Setup](docs/02_proteinmpnn_setup.md)
- [MD Simulation Protocol](docs/03_md_protocol.md)
- [Analysis Guide](docs/04_analysis_guide.md)

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this workflow, please cite the original tools:
- RFDiffusion: Watson et al., Nature, 2023
- ProteinMPNN: Dauparas et al., Science, 2022
- GROMACS: Abraham et al., SoftwareX, 2015

## Contact

For questions or issues, please open an issue on GitHub.