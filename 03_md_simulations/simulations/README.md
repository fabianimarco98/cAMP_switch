# This directory contains GROMACS MDP parameter files

MDP (Molecular Dynamics Parameters) files control simulation settings.

## Standard workflow files:

### ions.mdp
Minimal parameters for adding ions to the system.

### em.mdp
Energy minimization parameters.

### nvt.mdp
NVT (constant Number, Volume, Temperature) equilibration.

### npt.mdp
NPT (constant Number, Pressure, Temperature) equilibration.

### md.mdp
Production MD simulation parameters.

## Usage

These files are referenced by `grompp` when preparing simulations:
```bash
gmx grompp -f em.mdp -c input.gro -p topol.top -o em.tpr
```

## Customization

Modify these files to adjust:
- Simulation length (nsteps)
- Temperature and pressure (ref_t, ref_p)
- Output frequency (nstxout, nstenergy)
- Integrator and timestep (integrator, dt)

See `docs/03_md_protocol.md` for detailed parameter descriptions.
