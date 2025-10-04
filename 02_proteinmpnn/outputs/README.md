# This directory will contain outputs from ProteinMPNN

Generated sequences and models will be saved here.

## Expected files after running ProteinMPNN:
- `seqs/` - FASTA files with designed sequences
- `scores/` - CSV files with confidence scores
- `backbones/` - PDB files with sequences threaded onto backbones

Directory structure:
```
outputs/
├── seqs/
│   ├── design1.fa
│   └── design2.fa
├── scores/
│   ├── design1_scores.csv
│   └── design2_scores.csv
└── backbones/
    ├── design1.pdb
    └── design2.pdb
```

Files are excluded from git tracking (see .gitignore).
Select best sequences for MD simulations.
