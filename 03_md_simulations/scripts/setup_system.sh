#!/bin/bash
#
# Setup script for MD simulation system
# This script prepares a protein-ligand complex for GROMACS MD simulation
#

set -e  # Exit on error

echo "=== GROMACS System Setup Script ==="
echo ""

# Configuration
PROTEIN_PDB="${1:-protein.pdb}"
LIGAND_PDB="${2:-../../../data/cAMP.pdb}"
FORCE_FIELD="${3:-amber99sb-ildn}"
WATER_MODEL="${4:-tip3p}"
BOX_DISTANCE="1.0"  # nm from protein to box edge

echo "Input files:"
echo "  Protein: $PROTEIN_PDB"
echo "  Ligand: $LIGAND_PDB"
echo "  Force field: $FORCE_FIELD"
echo "  Water model: $WATER_MODEL"
echo ""

# Check if input files exist
if [ ! -f "$PROTEIN_PDB" ]; then
    echo "Error: Protein PDB file not found: $PROTEIN_PDB"
    exit 1
fi

if [ ! -f "$LIGAND_PDB" ]; then
    echo "Warning: Ligand PDB file not found: $LIGAND_PDB"
    echo "Proceeding with protein-only setup"
    LIGAND_PDB=""
fi

# Step 1: Generate protein topology
echo "Step 1: Generating protein topology..."
gmx pdb2gmx -f "$PROTEIN_PDB" -o protein_processed.gro -water "$WATER_MODEL" -ff "$FORCE_FIELD" -ignh

# Step 2: Combine protein and ligand if ligand is provided
if [ -n "$LIGAND_PDB" ]; then
    echo "Step 2: Combining protein and ligand..."
    # This is a placeholder - actual ligand parameterization is complex
    # and typically requires tools like ACPYPE, CGenFF, or ATB
    echo "  Note: Ligand topology must be prepared separately using ACPYPE or similar"
    echo "  For now, using protein-only system"
    cp protein_processed.gro complex.gro
else
    cp protein_processed.gro complex.gro
fi

# Step 3: Define simulation box
echo "Step 3: Defining simulation box..."
gmx editconf -f complex.gro -o box.gro -c -d "$BOX_DISTANCE" -bt cubic

# Step 4: Solvate the system
echo "Step 4: Solvating system..."
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top

# Step 5: Add ions
echo "Step 5: Adding ions for neutralization..."

# Create minimal ions.mdp if it doesn't exist
if [ ! -f "../simulations/ions.mdp" ]; then
    echo "Creating minimal ions.mdp..."
    cat > ions.mdp << EOF
; Ions MDP
integrator  = steep
emtol       = 1000.0
emstep      = 0.01
nsteps      = 5000
nstlist     = 1
cutoff-scheme   = Verlet
ns_type         = grid
coulombtype     = PME
rcoulomb        = 1.0
rvdw            = 1.0
pbc             = xyz
EOF
    MDP_FILE="ions.mdp"
else
    MDP_FILE="../simulations/ions.mdp"
fi

gmx grompp -f "$MDP_FILE" -c solvated.gro -p topol.top -o ions.tpr -maxwarn 1

echo "SOL" | gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral

echo ""
echo "=== Setup Complete ==="
echo ""
echo "Output files:"
echo "  - ionized.gro: Final system structure"
echo "  - topol.top: System topology"
echo ""
echo "Next steps:"
echo "  1. Run energy minimization: gmx grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr"
echo "  2. Continue with equilibration (see docs/03_md_protocol.md)"
echo ""
