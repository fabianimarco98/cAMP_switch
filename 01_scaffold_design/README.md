# Phase 1: De Novo Scaffold Design

This directory contains all the necessary files and scripts for the initial design of protein scaffolds using RFDiffusion. The goal of this phase is to generate a diverse set of protein backbones that possess a pre-defined binding pocket suitable for cAMP.

## Objective

The primary objective of this phase is to generate protein scaffolds featuring two distinct conformational states:
1.  A **"closed" state**, designed by conditioning the diffusion model on the geometry of a cAMP molecule within a target binding motif.
2.  An **"open" state**, generated using inpainting to create a plausible ligand-free conformation.

## Methodology

We employ **RFDiffusion**, specifically using the "inpainting" and "motif scaffolding" protocols. The design process is conditioned on a set of spatial constraints derived from an idealized cAMP binding pocket. Key interacting residues and their geometries are defined as a motif, which the model builds a stable scaffold around.

## Directory Structure
/01_scaffold_design
│
├── README.md             # This file
│
├── inputs/
│   ├── constraints/      # Files defining geometric constraints for the design
│   └── motifs/           # PDB files of the target binding motif (e.g., key residues)
│
├── scripts/
│   └── run_rfdiffusion.sh # Main script to execute the RFDiffusion jobs
│
└── outputs/
└── pdb/              # Directory where generated PDB scaffolds are saved

## Usage

To run this design phase, execute the main script from within the `/scripts` directory.

1.  **Navigate to the scripts directory:**
    ```bash
    cd 01_scaffold_design/scripts
    ```

2.  **Inspect and modify the script:**
    Open `run_rfdiffusion.sh` and ensure all paths to the RFDiffusion installation, input files, and output directories are correctly set. Adjust model parameters such as the number of designs (`n_designs`) and diffusion timesteps as needed.

3.  **Execute the script:**
    ```bash
    bash run_rfdiffusion.sh
    ```
    This process can be computationally intensive and should ideally be run on a high-performance computing (HPC) cluster with GPU resources.

## Outputs

The script will generate a specified number of protein scaffold designs in PDB format, which will be saved in the `outputs/pdb/` directory.

Each successful design will have a corresponding PDB file (e.g., `design_001.pdb`). These structures are the direct input for the next phase of the project, **Phase 02: Sequence Design**. The top candidates should be selected based on model confidence scores and visual inspection of their structural integrity and pocket geometry.
