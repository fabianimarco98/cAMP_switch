# GROMACS MD Simulation Protocol

## Overview
This guide describes the protocol for running molecular dynamics simulations of cAMP-protein complexes using GROMACS.

## Installation

### GROMACS Installation
```bash
# Option 1: Using conda
conda install -c conda-forge gromacs

# Option 2: From source (for GPU acceleration)
wget http://ftp.gromacs.org/gromacs/gromacs-2023.tar.gz
tar xfz gromacs-2023.tar.gz
cd gromacs-2023
mkdir build
cd build
cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA
make
sudo make install
source /usr/local/gromacs/bin/GMXRC
```

### Additional Tools
```bash
# For analysis and visualization
pip install MDAnalysis matplotlib numpy scipy
conda install -c conda-forge pymol-open-source
```

## Workflow Overview

The MD simulation workflow consists of:
1. System preparation
2. Energy minimization
3. NVT equilibration
4. NPT equilibration
5. Production MD
6. Analysis

## 1. System Preparation

### Step 1.1: Prepare Protein-Ligand Complex
```bash
cd 03_md_simulations/setup

# Copy designed structure
cp ../../02_proteinmpnn/outputs/best_design.pdb protein.pdb

# Prepare cAMP ligand parameters
# Option A: Use existing force field parameters (e.g., AMBER GAFF)
# Option B: Generate parameters using ACPYPE or ATB server
```

### Step 1.2: Generate Topology
```bash
# Generate protein topology
gmx pdb2gmx -f protein.pdb -o protein_processed.gro -water tip3p -ff amber99sb-ildn

# If cAMP parameters are in separate file, combine topologies
# Edit topol.top to include cAMP parameters
```

### Step 1.3: Define Simulation Box
```bash
# Create cubic box with 1.0 nm minimum distance to edges
gmx editconf -f protein_processed.gro -o box.gro -c -d 1.0 -bt cubic
```

### Step 1.4: Solvate the System
```bash
# Add water molecules
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
```

### Step 1.5: Add Ions
```bash
# Prepare for ion addition
gmx grompp -f ../simulations/ions.mdp -c solvated.gro -p topol.top -o ions.tpr

# Add ions to neutralize system
gmx genion -s ions.tpr -o ionized.gro -p topol.top -pname NA -nname CL -neutral
```

## 2. Energy Minimization

### MDP Parameters (ions.mdp and em.mdp)
Create `03_md_simulations/simulations/em.mdp`:
```
; Energy Minimization Parameters
integrator  = steep         ; Steepest descent
emtol       = 1000.0        ; Stop when max force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Initial step size
nsteps      = 50000         ; Maximum number of steps

; Output control
nstxout     = 1000          ; Write coordinates every 1000 steps
nstvout     = 1000          ; Write velocities every 1000 steps
nstlog      = 100           ; Update log file every 100 steps
nstenergy   = 100           ; Write energies every 100 steps

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; van der Waals
vdwtype         = Cut-off
rvdw            = 1.0

; PBC
pbc             = xyz
```

### Run Energy Minimization
```bash
cd 03_md_simulations/setup

# Prepare energy minimization
gmx grompp -f ../simulations/em.mdp -c ionized.gro -p topol.top -o em.tpr

# Run minimization
gmx mdrun -v -deffnm em

# Check results
gmx energy -f em.edr -o potential.xvg
# Select "Potential" when prompted
```

## 3. NVT Equilibration

### NVT MDP Parameters
Create `03_md_simulations/simulations/nvt.mdp`:
```
; NVT Equilibration
integrator  = md            ; Leap-frog integrator
nsteps      = 50000         ; 100 ps
dt          = 0.002         ; 2 fs

; Output control
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

; Bond parameters
continuation    = no
constraint_algorithm = lincs
constraints     = h-bonds
lincs_iter      = 1
lincs_order     = 4

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; van der Waals
vdwtype         = Cut-off
rvdw            = 1.0

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = 300     300

; PBC
pbc             = xyz

; Velocity generation
gen_vel         = yes
gen_temp        = 300
gen_seed        = -1
```

### Run NVT Equilibration
```bash
# Prepare NVT
gmx grompp -f ../simulations/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr

# Run NVT equilibration
gmx mdrun -v -deffnm nvt

# Check temperature convergence
gmx energy -f nvt.edr -o temperature.xvg
# Select "Temperature" when prompted
```

## 4. NPT Equilibration

### NPT MDP Parameters
Create `03_md_simulations/simulations/npt.mdp`:
```
; NPT Equilibration
integrator  = md
nsteps      = 50000         ; 100 ps
dt          = 0.002

; Output control
nstxout     = 500
nstvout     = 500
nstenergy   = 500
nstlog      = 500

; Bond parameters
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics and VdW (same as NVT)
coulombtype     = PME
rcoulomb        = 1.0
vdwtype         = Cut-off
rvdw            = 1.0

; Temperature coupling (same as NVT)
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = 300     300

; Pressure coupling
pcoupl          = Berendsen
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

; PBC
pbc             = xyz

; Velocity generation
gen_vel         = no
```

### Run NPT Equilibration
```bash
# Prepare NPT
gmx grompp -f ../simulations/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr

# Run NPT equilibration
gmx mdrun -v -deffnm npt

# Check pressure and density convergence
gmx energy -f npt.edr -o pressure.xvg
gmx energy -f npt.edr -o density.xvg
```

## 5. Production MD

### Production MDP Parameters
Create `03_md_simulations/simulations/md.mdp`:
```
; Production MD
integrator  = md
nsteps      = 50000000      ; 100 ns (adjust as needed)
dt          = 0.002

; Output control
nstxout-compressed  = 5000  ; Save compressed coordinates every 10 ps
compressed-x-grps   = System
nstvout     = 0             ; Don't save velocities (save space)
nstfout     = 0             ; Don't save forces
nstenergy   = 5000
nstlog      = 5000

; Bond parameters
continuation    = yes
constraint_algorithm = lincs
constraints     = h-bonds

; Neighbor searching
cutoff-scheme   = Verlet
ns_type         = grid
nstlist         = 10
rlist           = 1.0

; Electrostatics
coulombtype     = PME
rcoulomb        = 1.0

; van der Waals
vdwtype         = Cut-off
rvdw            = 1.0

; Temperature coupling
tcoupl          = V-rescale
tc-grps         = Protein Non-Protein
tau_t           = 0.1     0.1
ref_t           = 300     300

; Pressure coupling
pcoupl          = Parrinello-Rahman
pcoupltype      = isotropic
tau_p           = 2.0
ref_p           = 1.0
compressibility = 4.5e-5

; PBC
pbc             = xyz

; Velocity generation
gen_vel         = no
```

### Run Production MD
```bash
# Prepare production run
gmx grompp -f ../simulations/md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

# Run production MD (GPU acceleration recommended)
gmx mdrun -v -deffnm md -nb gpu

# For long simulations, consider using:
gmx mdrun -v -deffnm md -nb gpu -maxh 24  # Stop after 24 hours
# Then continue with:
gmx mdrun -v -deffnm md -nb gpu -cpi md.cpt -append
```

## 6. Basic Analysis

### RMSD Analysis
```bash
cd 03_md_simulations/trajectories

# Calculate RMSD
gmx rms -s ../setup/md.tpr -f ../setup/md.xtc -o rmsd.xvg -tu ns
# Select backbone for both groups
```

### RMSF Analysis
```bash
# Calculate RMSF
gmx rmsf -s ../setup/md.tpr -f ../setup/md.xtc -o rmsf.xvg -res
```

### Hydrogen Bonds
```bash
# Analyze protein-ligand hydrogen bonds
gmx hbond -s ../setup/md.tpr -f ../setup/md.xtc -num hbnum.xvg
```

### Radius of Gyration
```bash
gmx gyrate -s ../setup/md.tpr -f ../setup/md.xtc -o gyrate.xvg
```

## Troubleshooting

### Common Issues
1. **System explodes during minimization**: Check topology, especially ligand parameters
2. **Temperature/Pressure fluctuations**: Normal in small systems; check coupling parameters
3. **Slow performance**: Use GPU acceleration, optimize nstlist
4. **LINCS warnings**: Reduce time step or check constraints

## Best Practices

1. **Always visualize the system** before starting long simulations
2. **Monitor initial equilibration** closely
3. **Save checkpoints frequently** for long runs
4. **Keep raw data** until analysis is complete
5. **Document all parameters** used

## Next Steps

After MD simulations:
1. Copy trajectory files to `04_analysis/` for detailed analysis
2. Run comprehensive binding analysis
3. Calculate free energies
4. Compare multiple designs

## References
- GROMACS Manual: https://manual.gromacs.org/
- Best practices: https://www.mdtutorials.com/gmx/
- Force field selection: Robustelli et al., PNAS (2018)
