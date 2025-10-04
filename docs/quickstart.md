# Quick Start Guide

This guide will help you get started with the cAMP binder design project.

## Prerequisites

Before starting, ensure you have:
- Linux/Unix system (or WSL on Windows)
- Python 3.8+
- Conda package manager
- NVIDIA GPU with CUDA (recommended)

## Step 1: Installation

1. Clone the repository:
```bash
git clone https://github.com/fabianimarco98/cAMP_switch.git
cd cAMP_switch
```

2. Install dependencies:
```bash
# Create conda environment
conda create -n camp_binder python=3.9
conda activate camp_binder

# Install Python packages
pip install -r requirements.txt
```

3. Install external tools (see `docs/installation.md` for details):
   - RFDiffusion
   - ProteinMPNN
   - GROMACS

## Step 2: Prepare Input Data

1. Obtain cAMP structure:
```bash
# Download from PDB or generate using molecular modeling
# Save as data/cAMP.pdb
```

2. Verify structure:
```bash
# Use PyMOL or similar to visualize
pymol data/cAMP.pdb
```

## Step 3: Run RFDiffusion

1. Configure RFDiffusion:
```bash
cd 01_rfdiffusion/inputs
cp config_example.yaml config.yaml
# Edit config.yaml with your parameters
```

2. Run structure generation:
```bash
conda activate SE3nv  # RFDiffusion environment

# Run RFDiffusion (adjust path to your installation)
python $RFDIFFUSION_HOME/scripts/run_inference.py \
  inference.output_prefix=01_rfdiffusion/outputs/cAMP_binder \
  inference.num_designs=100 \
  'contigmap.contigs=[A1-100]' \
  inference.input_pdb=data/cAMP.pdb
```

3. Rank and select designs:
```bash
cd 01_rfdiffusion
python scripts/rank_designs.py \
  --input_dir outputs/ \
  --output_csv outputs/rankings.csv \
  --top_n 10
```

## Step 4: Run ProteinMPNN

1. Copy selected backbones:
```bash
cp 01_rfdiffusion/outputs/top_designs/*.pdb 02_proteinmpnn/inputs/
```

2. Run sequence design:
```bash
conda activate camp_binder

cd 02_proteinmpnn
for pdb in inputs/*.pdb; do
  python $PROTEINMPNN_HOME/protein_mpnn_run.py \
    --pdb_path "$pdb" \
    --out_folder outputs/ \
    --num_seq_per_target 10 \
    --sampling_temp 0.1
done
```

3. Analyze sequences:
```bash
python scripts/analyze_sequences.py \
  --fasta outputs/seqs/top_01*.fa \
  --output_csv outputs/sequence_analysis.csv
```

## Step 5: Setup MD Simulations

1. Prepare system:
```bash
cd 03_md_simulations/setup

# Copy best design from ProteinMPNN
cp ../../02_proteinmpnn/outputs/backbones/best_design.pdb protein.pdb

# Run setup script
bash ../scripts/setup_system.sh protein.pdb ../../../data/cAMP.pdb
```

2. Run MD workflow:
```bash
# This will run energy minimization and equilibration
bash ../scripts/run_md_workflow.sh 4  # Use 4 CPU cores

# Start production MD (long running)
gmx mdrun -v -deffnm md -nt 4 -nb gpu
```

## Step 6: Analyze Results

1. Wait for MD to complete (or stop at desired time)

2. Run analysis pipeline:
```bash
cd 04_analysis
python scripts/run_analysis.py
```

3. Review results:
```bash
# View plots
ls -lh figures/
ls -lh results/

# Read summary report
cat results/REPORT.md
```

## Expected Timeline

- **RFDiffusion**: 1-4 hours (depending on number of designs)
- **ProteinMPNN**: 10-30 minutes per design
- **MD Setup**: 30 minutes
- **MD Simulation**: 1-24 hours per 100 ns (depending on system size and hardware)
- **Analysis**: 1-2 hours

## Troubleshooting

### Common Issues

1. **CUDA out of memory (RFDiffusion)**
   - Reduce number of designs
   - Use smaller scaffold size

2. **GROMACS errors during setup**
   - Check protein structure quality
   - Verify ligand parameters
   - See `docs/03_md_protocol.md`

3. **Analysis scripts fail**
   - Ensure MD completed successfully
   - Check file paths in scripts
   - Verify Python packages installed

### Getting Help

- Check documentation in `docs/` directory
- Review error messages carefully
- Open an issue on GitHub with:
  - Error message
  - Command that failed
  - System information

## Next Steps

After completing initial workflow:
1. Compare multiple designs
2. Refine top candidates
3. Perform longer MD simulations
4. Calculate binding free energies
5. Plan experimental validation

## Best Practices

1. **Keep organized**: Document all parameters and decisions
2. **Version control**: Use git to track changes
3. **Backup data**: Save important results regularly
4. **Monitor resources**: Check CPU/GPU usage and disk space
5. **Iterate**: Use insights from analysis to improve designs

## Additional Resources

- [Full Installation Guide](docs/installation.md)
- [RFDiffusion Details](docs/01_rfdiffusion_setup.md)
- [ProteinMPNN Details](docs/02_proteinmpnn_setup.md)
- [MD Protocol](docs/03_md_protocol.md)
- [Analysis Guide](docs/04_analysis_guide.md)

Happy designing! 🧬
