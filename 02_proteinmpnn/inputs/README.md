# This directory contains input files for ProteinMPNN

Copy your top-ranked backbone structures from RFDiffusion here.

## Files to include:
- Backbone PDB files from RFDiffusion (top designs)
- `chain_definitions.jsonl` - Optional file specifying which chains to design

Example chain_definitions.jsonl:
```json
{"cAMP_binder_top1.pdb": {"designed_chains": ["A"], "fixed_chains": []}}
{"cAMP_binder_top2.pdb": {"designed_chains": ["A"], "fixed_chains": []}}
```

See `docs/02_proteinmpnn_setup.md` for detailed instructions.
