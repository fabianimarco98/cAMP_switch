# Analysis Guide

## Overview
This guide covers the analysis of MD simulation trajectories and evaluation of designed cAMP binders.

## Required Software

```bash
# Python packages
pip install MDAnalysis matplotlib seaborn numpy scipy pandas
pip install biopython prody

# For visualization
conda install -c conda-forge pymol-open-source nglview

# For binding free energy calculations
# gmx_mmpbsa or g_mmpbsa (see installation below)
```

## Analysis Workflow

### 1. Trajectory Preprocessing

#### Remove Periodic Boundary Artifacts
```bash
cd 04_analysis

# Center protein and handle PBC
gmx trjconv -s ../03_md_simulations/setup/md.tpr \
  -f ../03_md_simulations/setup/md.xtc \
  -o md_centered.xtc \
  -center -pbc mol -ur compact
# Select "Protein" for centering and "System" for output
```

#### Extract Specific Frames
```bash
# Extract every 10th frame
gmx trjconv -s ../03_md_simulations/setup/md.tpr \
  -f md_centered.xtc \
  -o md_reduced.xtc \
  -skip 10

# Extract last 50 ns for analysis
gmx trjconv -s ../03_md_simulations/setup/md.tpr \
  -f md_centered.xtc \
  -o md_last50ns.xtc \
  -b 50000
```

### 2. Structural Stability Analysis

#### RMSD Analysis Script
Create `04_analysis/scripts/calculate_rmsd.py`:
```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt
import numpy as np

# Load trajectory
u = mda.Universe('../03_md_simulations/setup/md.tpr', 
                 'md_centered.xtc')

# Select protein backbone
protein = u.select_atoms('protein and name CA')

# Calculate RMSD
R = rms.RMSD(protein, protein, select='backbone')
R.run()

# Plot RMSD
plt.figure(figsize=(10, 6))
plt.plot(R.rmsd[:, 1], R.rmsd[:, 2])
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (Å)')
plt.title('Backbone RMSD over time')
plt.savefig('results/rmsd.png', dpi=300)
plt.close()

# Save data
np.savetxt('results/rmsd.txt', R.rmsd[:, [1, 2]], 
           header='Time(ps) RMSD(A)')
```

#### RMSF Analysis Script
Create `04_analysis/scripts/calculate_rmsf.py`:
```python
import MDAnalysis as mda
from MDAnalysis.analysis import rms
import matplotlib.pyplot as plt

u = mda.Universe('../03_md_simulations/setup/md.tpr', 
                 'md_centered.xtc')

# Select C-alpha atoms
ca = u.select_atoms('protein and name CA')

# Calculate RMSF
R = rms.RMSF(ca).run()

# Plot RMSF
plt.figure(figsize=(12, 6))
plt.plot(ca.resids, R.rmsf)
plt.xlabel('Residue Number')
plt.ylabel('RMSF (Å)')
plt.title('Per-residue flexibility')
plt.savefig('results/rmsf.png', dpi=300)
plt.close()

# Save data
np.savetxt('results/rmsf.txt', 
           np.column_stack([ca.resids, R.rmsf]),
           header='Residue RMSF(A)')
```

### 3. Binding Analysis

#### Contact Analysis
Create `04_analysis/scripts/analyze_contacts.py`:
```python
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import matplotlib.pyplot as plt
import numpy as np

u = mda.Universe('../03_md_simulations/setup/md.tpr', 
                 'md_centered.xtc')

# Define selections
protein = u.select_atoms('protein')
ligand = u.select_atoms('resname cAMP')  # Adjust resname as needed

# Calculate contacts
ca = contacts.Contacts(u, 
                       select=('protein', 'resname cAMP'),
                       refgroup=(protein, ligand),
                       radius=4.5)
ca.run()

# Plot contact timeline
plt.figure(figsize=(10, 6))
plt.plot(ca.timeseries[:, 0], ca.timeseries[:, 1])
plt.xlabel('Time (ps)')
plt.ylabel('Number of Contacts')
plt.title('Protein-cAMP Contacts')
plt.savefig('results/contacts.png', dpi=300)
plt.close()
```

#### Hydrogen Bond Analysis
Create `04_analysis/scripts/analyze_hbonds.py`:
```python
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
import matplotlib.pyplot as plt

u = mda.Universe('../03_md_simulations/setup/md.tpr', 
                 'md_centered.xtc')

# Analyze hydrogen bonds between protein and ligand
hbonds = HydrogenBondAnalysis(
    universe=u,
    donors_sel='protein',
    hydrogens_sel='protein',
    acceptors_sel='resname cAMP'
)
hbonds.run()

# Get results
print(f"Total hydrogen bonds found: {len(hbonds.results.hbonds)}")

# Analyze persistence
# Plot histogram of H-bond lifetimes
plt.figure(figsize=(10, 6))
plt.hist(hbonds.count_by_time(), bins=50)
plt.xlabel('Number of H-bonds')
plt.ylabel('Frequency')
plt.title('Hydrogen Bond Distribution')
plt.savefig('results/hbond_distribution.png', dpi=300)
plt.close()
```

### 4. Binding Free Energy Calculations

#### Using gmx_MMPBSA
```bash
# Install gmx_MMPBSA
pip install gmx_MMPBSA

# Create input file
cat > mmpbsa.in << EOF
&general
startframe=500,
endframe=1000,
interval=5,
/

&gb
igb=5,
saltcon=0.150,
/
EOF

# Run MM-PBSA calculation
gmx_MMPBSA -O -i mmpbsa.in \
  -cs ../03_md_simulations/setup/md.tpr \
  -ct md_centered.xtc \
  -ci ../03_md_simulations/setup/md.tpr \
  -cg 1 13 \
  -o results/FINAL_RESULTS_MMPBSA.dat
```

#### Energy Decomposition
Create `04_analysis/scripts/decompose_energy.py`:
```python
import pandas as pd
import matplotlib.pyplot as plt

# Parse MM-PBSA results (adjust based on output format)
# This is a template - actual parsing depends on gmx_MMPBSA output

results = {
    'VDWAALS': -45.2,
    'EEL': -78.3,
    'EGB': 65.4,
    'ESURF': -3.2,
    'DELTA_G': -61.3
}

# Plot energy components
plt.figure(figsize=(10, 6))
components = list(results.keys())[:-1]
energies = [results[k] for k in components]

plt.bar(components, energies)
plt.xlabel('Energy Component')
plt.ylabel('Energy (kcal/mol)')
plt.title('Binding Free Energy Decomposition')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('results/energy_decomposition.png', dpi=300)
plt.close()
```

### 5. Comparative Analysis

#### Compare Multiple Designs
Create `04_analysis/scripts/compare_designs.py`:
```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data for multiple designs
designs = {
    'Design_1': {
        'rmsd_mean': 2.3,
        'rmsd_std': 0.4,
        'binding_energy': -61.3,
        'contacts_mean': 12.5,
        'hbonds_mean': 3.2
    },
    'Design_2': {
        'rmsd_mean': 3.1,
        'rmsd_std': 0.6,
        'binding_energy': -45.7,
        'contacts_mean': 10.2,
        'hbonds_mean': 2.5
    },
    # Add more designs...
}

df = pd.DataFrame(designs).T

# Create comparison plots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# RMSD comparison
axes[0, 0].bar(df.index, df['rmsd_mean'], yerr=df['rmsd_std'])
axes[0, 0].set_ylabel('RMSD (Å)')
axes[0, 0].set_title('Structural Stability')
axes[0, 0].tick_params(axis='x', rotation=45)

# Binding energy
axes[0, 1].bar(df.index, df['binding_energy'])
axes[0, 1].set_ylabel('Binding Energy (kcal/mol)')
axes[0, 1].set_title('Binding Affinity')
axes[0, 1].tick_params(axis='x', rotation=45)

# Contacts
axes[1, 0].bar(df.index, df['contacts_mean'])
axes[1, 0].set_ylabel('Number of Contacts')
axes[1, 0].set_title('Protein-Ligand Contacts')
axes[1, 0].tick_params(axis='x', rotation=45)

# Hydrogen bonds
axes[1, 1].bar(df.index, df['hbonds_mean'])
axes[1, 1].set_ylabel('Number of H-bonds')
axes[1, 1].set_title('Hydrogen Bonds')
axes[1, 1].tick_params(axis='x', rotation=45)

plt.tight_layout()
plt.savefig('results/design_comparison.png', dpi=300)
plt.close()

# Save summary table
df.to_csv('results/design_summary.csv')
```

### 6. Visualization

#### Create PyMOL Visualization Script
Create `04_analysis/scripts/visualize_binding.py`:
```python
from pymol import cmd

# Load structure
cmd.load('../03_md_simulations/setup/md.gro', 'protein')
cmd.load('../data/cAMP.pdb', 'ligand')

# Styling
cmd.hide('everything')
cmd.show('cartoon', 'protein')
cmd.show('sticks', 'ligand')
cmd.color('cyan', 'protein')
cmd.color('yellow', 'ligand')

# Show binding site residues
cmd.select('binding_site', 'protein within 5 of ligand')
cmd.show('sticks', 'binding_site')
cmd.color('green', 'binding_site')

# Add hydrogen bonds
cmd.distance('hbonds', 'protein', 'ligand', 3.5, mode=2)

# Set view
cmd.zoom('ligand', 5)

# Save image
cmd.png('figures/binding_site.png', width=1200, height=1200, dpi=300, ray=1)

# Save session
cmd.save('figures/binding_visualization.pse')
```

### 7. Generate Final Report

#### Summary Report Script
Create `04_analysis/scripts/generate_report.py`:
```python
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime

# Compile all results
report_data = {
    'Date': datetime.now().strftime('%Y-%m-%d'),
    'Design': 'cAMP_binder_v1',
    'RMSD_mean': 2.3,
    'RMSD_std': 0.4,
    'RMSF_max': 4.5,
    'Binding_Energy': -61.3,
    'Contacts_mean': 12.5,
    'HBonds_mean': 3.2,
    'Simulation_time': '100 ns'
}

# Create markdown report
report = f"""
# cAMP Binder Analysis Report

**Generated:** {report_data['Date']}
**Design:** {report_data['Design']}

## Structural Stability
- Mean RMSD: {report_data['RMSD_mean']:.2f} ± {report_data['RMSD_std']:.2f} Å
- Maximum RMSF: {report_data['RMSF_max']:.2f} Å

## Binding Metrics
- Binding Free Energy: {report_data['Binding_Energy']:.2f} kcal/mol
- Average Contacts: {report_data['Contacts_mean']:.1f}
- Average H-bonds: {report_data['HBonds_mean']:.1f}

## Conclusion
[Add interpretation and recommendations]

## Figures
- RMSD plot: results/rmsd.png
- RMSF plot: results/rmsf.png
- Contact analysis: results/contacts.png
- Binding site visualization: figures/binding_site.png
"""

with open('results/REPORT.md', 'w') as f:
    f.write(report)

print("Report generated: results/REPORT.md")
```

## Running the Complete Analysis

Create master script `04_analysis/scripts/run_full_analysis.sh`:
```bash
#!/bin/bash

echo "Starting analysis pipeline..."

# Preprocess trajectory
bash preprocess_trajectory.sh

# Run Python analyses
python calculate_rmsd.py
python calculate_rmsf.py
python analyze_contacts.py
python analyze_hbonds.py

# Run MM-PBSA (if applicable)
# bash run_mmpbsa.sh

# Generate visualizations
python visualize_binding.py

# Compare designs (if multiple)
python compare_designs.py

# Generate final report
python generate_report.py

echo "Analysis complete! Check results/ and figures/ directories."
```

## Quality Metrics Summary

Track these key metrics:
1. **Stability**: RMSD < 3 Å (good), < 2 Å (excellent)
2. **Flexibility**: RMSF < 2 Å for binding site residues
3. **Binding**: ΔG < -5 kcal/mol (weak), < -10 kcal/mol (strong)
4. **Interactions**: > 3 persistent H-bonds, > 10 contacts

## References
- MDAnalysis: Michaud-Agrawal et al., J. Comput. Chem. (2011)
- gmx_MMPBSA: Valdés-Tresanco et al., J. Chem. Theory Comput. (2021)
- Best practices: Gowers et al., SciPy (2016)
