# This directory will contain outputs from RFDiffusion

Generated PDB files will be saved here.

## Expected files after running RFDiffusion:
- `cAMP_binder_0.pdb`, `cAMP_binder_1.pdb`, etc. - Generated designs
- `rankings.csv` - Ranked designs by confidence (after running rank_designs.py)
- `top_designs/` - Subdirectory with top-ranked designs (optional)

Files are excluded from git tracking (see .gitignore).
Keep only the best designs for the next pipeline stage.
