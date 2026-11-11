Output of the DMRG runs (`scripts/submit.sh` writes here), with the structure written by `observables_calc.jl`:

```
data/
└── results_nu1/
    └── Lx44/
        ├── fidelity/          fidelity_L44_J1.00_Jp0.10_U1.200_maxdim2000_cutoff1e-10.txt
        ├── Lz/                Lz_L44_...
        ├── Lz_corr/           Lz_corrL44_...
        ├── compressibility/   comp_L44_...
```

The analysis notebook in `../analysis` reads the files from here.
