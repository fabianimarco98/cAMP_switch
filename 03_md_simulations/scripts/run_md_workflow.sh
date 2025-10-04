#!/bin/bash
#
# Run complete MD simulation workflow
# This script runs energy minimization, equilibration, and production MD
#

set -e

echo "=== GROMACS MD Simulation Workflow ==="
echo ""

# Check if system is set up
if [ ! -f "ionized.gro" ] || [ ! -f "topol.top" ]; then
    echo "Error: System not set up. Run setup_system.sh first."
    exit 1
fi

# Configuration
NPROC="${1:-4}"  # Number of processors
USE_GPU="${2:-auto}"  # GPU usage: auto, yes, no

echo "Configuration:"
echo "  Processors: $NPROC"
echo "  GPU: $USE_GPU"
echo ""

# Energy Minimization
echo "=== Step 1: Energy Minimization ==="
if [ ! -f "em.gro" ]; then
    echo "Running energy minimization..."
    gmx grompp -f ../simulations/em.mdp -c ionized.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em -nt "$NPROC"
    
    echo "Checking minimization..."
    echo "Potential" | gmx energy -f em.edr -o potential.xvg
    echo "Energy minimization complete!"
else
    echo "Energy minimization already done (em.gro exists)"
fi

# NVT Equilibration
echo ""
echo "=== Step 2: NVT Equilibration ==="
if [ ! -f "nvt.gro" ]; then
    echo "Running NVT equilibration..."
    gmx grompp -f ../simulations/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -v -deffnm nvt -nt "$NPROC" -nb "$USE_GPU"
    
    echo "Checking temperature..."
    echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg
    echo "NVT equilibration complete!"
else
    echo "NVT equilibration already done (nvt.gro exists)"
fi

# NPT Equilibration
echo ""
echo "=== Step 3: NPT Equilibration ==="
if [ ! -f "npt.gro" ]; then
    echo "Running NPT equilibration..."
    gmx grompp -f ../simulations/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    gmx mdrun -v -deffnm npt -nt "$NPROC" -nb "$USE_GPU"
    
    echo "Checking pressure and density..."
    echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
    echo "Density" | gmx energy -f npt.edr -o density.xvg
    echo "NPT equilibration complete!"
else
    echo "NPT equilibration already done (npt.gro exists)"
fi

# Production MD
echo ""
echo "=== Step 4: Production MD ==="
echo "Preparing production MD..."
gmx grompp -f ../simulations/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

echo ""
echo "Ready to run production MD!"
echo "To start production run:"
echo "  gmx mdrun -v -deffnm md -nt $NPROC -nb $USE_GPU"
echo ""
echo "For long runs, consider using:"
echo "  gmx mdrun -v -deffnm md -nt $NPROC -nb $USE_GPU -maxh 24"
echo "  (stops after 24 hours, can be continued with -cpi md.cpt -append)"
echo ""
