# Data Directory

This directory contains shared data files used across the pipeline.

## Expected files:

### cAMP.pdb
Structure of cyclic AMP (cAMP) ligand.
- Can be obtained from PDB (e.g., from protein-cAMP complexes)
- Or generated using molecular modeling tools
- Should be properly formatted with correct atom naming

### Other potential data:
- Reference structures
- Force field parameters for cAMP
- Template scaffolds
- Experimental data for validation

## Preparing cAMP Structure

If you need to generate cAMP structure:

1. **From PDB database:**
   - Search for structures containing cAMP
   - Extract cAMP coordinates
   - Clean up and save as separate PDB

2. **Using molecular modeling:**
   ```bash
   # Using Open Babel
   obabel -ismi -osdf -O cAMP.sdf
   obabel cAMP.sdf -O cAMP.pdb
   ```

3. **Parameterization for MD:**
   - Use ACPYPE for AMBER force field
   - Or ATB (Automated Topology Builder) for GROMOS
   - Or CGenFF for CHARMM

See `docs/03_md_protocol.md` for details on ligand parameterization.
