# Installation Guide

This guide covers the installation of all software required for the cAMP binder design project.

## Table of Contents
1. [System Requirements](#system-requirements)
2. [Python Environment Setup](#python-environment-setup)
3. [RFDiffusion Installation](#rfdiffusion-installation)
4. [ProteinMPNN Installation](#proteinmpnn-installation)
5. [GROMACS Installation](#gromacs-installation)
6. [Verification](#verification)

## System Requirements

### Hardware
- **CPU**: Modern multi-core processor (8+ cores recommended)
- **RAM**: Minimum 16 GB, 32+ GB recommended
- **GPU**: NVIDIA GPU with CUDA support (for RFDiffusion and GROMACS)
  - Minimum 8 GB VRAM for RFDiffusion
  - CUDA 11.0 or later
- **Storage**: 50+ GB free space for software and data

### Operating System
- Linux (Ubuntu 20.04+ recommended)
- macOS (with limitations on GPU acceleration)
- Windows with WSL2 (for Linux compatibility)

## Python Environment Setup

We recommend using conda for environment management:

```bash
# Install Miniconda (if not already installed)
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Create base environment for the project
conda create -n camp_binder python=3.9
conda activate camp_binder

# Install Python packages
cd /path/to/cAMP_switch
pip install -r requirements.txt
```

## RFDiffusion Installation

### Prerequisites
```bash
conda activate camp_binder
conda install -c conda-forge cudatoolkit=11.3 cudnn=8.2
```

### Installation Steps
```bash
# Clone repository
cd ~/software  # or your preferred location
git clone https://github.com/RosettaCommons/RFdiffusion.git
cd RFdiffusion

# Create dedicated environment
conda env create -f env/SE3nv.yml
conda activate SE3nv

# Install RFDiffusion
pip install -e .

# Download model weights
mkdir models && cd models
wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
cd ../..
```

### Verify Installation
```bash
conda activate SE3nv
python -c "import rfdiffusion; print('RFDiffusion installed successfully')"
```

## ProteinMPNN Installation

```bash
# Clone repository
cd ~/software
git clone https://github.com/dauparas/ProteinMPNN.git
cd ProteinMPNN

# Install in base environment
conda activate camp_binder
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Test installation
python protein_mpnn_run.py --help
```

## GROMACS Installation

### Option 1: Conda Installation (Easiest)
```bash
conda activate camp_binder
conda install -c conda-forge gromacs
```

### Option 2: From Source (For GPU Support)
```bash
# Install dependencies
sudo apt-get update
sudo apt-get install -y cmake build-essential libfftw3-dev

# Download GROMACS
cd ~/software
wget http://ftp.gromacs.org/gromacs/gromacs-2023.tar.gz
tar xfz gromacs-2023.tar.gz
cd gromacs-2023

# Configure with GPU support
mkdir build && cd build
cmake .. \
    -DGMX_BUILD_OWN_FFTW=ON \
    -DGMX_GPU=CUDA \
    -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
    -DGMX_SIMD=AVX2_256

# Build (this will take 15-30 minutes)
make -j8

# Install
sudo make install

# Source GROMACS environment
source /usr/local/gromacs/bin/GMXRC
```

### Verify GROMACS Installation
```bash
gmx --version
```

Expected output should show GROMACS version and GPU support (if compiled with CUDA).

## Additional Tools

### PyMOL (for visualization)
```bash
conda activate camp_binder
conda install -c conda-forge pymol-open-source
```

### Jupyter Lab (for interactive analysis)
```bash
conda activate camp_binder
conda install -c conda-forge jupyterlab
```

## Environment Configuration

### Setting Environment Variables
Add to your `~/.bashrc` or `~/.zshrc`:

```bash
# GROMACS
source /usr/local/gromacs/bin/GMXRC  # if installed from source

# RFDiffusion
export RFDIFFUSION_HOME=~/software/RFdiffusion

# ProteinMPNN
export PROTEINMPNN_HOME=~/software/ProteinMPNN

# CUDA (if applicable)
export CUDA_HOME=/usr/local/cuda
export PATH=$CUDA_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
```

Then reload:
```bash
source ~/.bashrc
```

## Verification

Run the verification script:

```bash
cd /path/to/cAMP_switch
python << EOF
import sys
print("Python version:", sys.version)

try:
    import numpy
    print("✓ NumPy:", numpy.__version__)
except ImportError:
    print("✗ NumPy not found")

try:
    import pandas
    print("✓ Pandas:", pandas.__version__)
except ImportError:
    print("✗ Pandas not found")

try:
    import MDAnalysis
    print("✓ MDAnalysis:", MDAnalysis.__version__)
except ImportError:
    print("✗ MDAnalysis not found")

try:
    import torch
    print("✓ PyTorch:", torch.__version__)
    if torch.cuda.is_available():
        print("  CUDA available:", torch.cuda.get_device_name(0))
    else:
        print("  CUDA not available")
except ImportError:
    print("✗ PyTorch not found")
EOF
```

Check GROMACS:
```bash
gmx --version
```

## Troubleshooting

### CUDA Issues
If you encounter CUDA-related errors:
```bash
# Check NVIDIA driver
nvidia-smi

# Verify CUDA installation
nvcc --version

# Update environment variables
export CUDA_HOME=/usr/local/cuda
export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH
```

### Import Errors
If Python packages fail to import:
```bash
# Ensure correct environment is activated
conda activate camp_binder

# Reinstall problematic package
pip install --upgrade --force-reinstall <package_name>
```

### GROMACS GPU Issues
If GROMACS doesn't detect GPU:
```bash
# Recompile with correct CUDA path
cmake .. -DGMX_GPU=CUDA -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda
```

## Next Steps

After successful installation:
1. Prepare input data (cAMP structure)
2. Start with RFDiffusion to generate binder scaffolds
3. Follow the workflow in main README.md

For detailed usage instructions, see the documentation in `docs/` directory.

## Getting Help

- **RFDiffusion**: https://github.com/RosettaCommons/RFdiffusion/issues
- **ProteinMPNN**: https://github.com/dauparas/ProteinMPNN/issues
- **GROMACS**: https://gromacs.org/documentation.html
- **Project issues**: Open an issue on the project GitHub repository
