# DMRG for the π-flux Bose-Hubbard ladder

DMRG simulations of interacting bosons on a two-leg ladder with a π flux per plaquette and alternating hopping. The code computes the ground state and the observables used to map the phase diagram: loop currents and their correlations, fidelity susceptibility, charge gap and entanglement entropy.

The DMRG is done with [ITensor](https://itensor.org/) (`ITensors.jl` / `ITensorMPS.jl`). This repository contains the model, the observables, the simulation workflow and the analysis.

<!-- TODO: add the reference of the paper with these results -->

## Model

```
 2 — 4 — 6 — 8 — ...      upper leg, Peierls phase e^{iφ} on the horizontal bonds
 |   |   |   |
 1 — 3 — 5 — 7 — ...      lower leg
   J   J'  J   J'
```

$$\hat H = -\sum_{\langle ij\rangle} J_{ij}\, e^{i\theta_{ij}}\,\hat b_i^\dagger \hat b_j + \frac{U}{2}\sum_i \hat n_i^2 ,$$

where:

- the hopping alternates between $J$ (strong plaquettes) and $J'$ (weak plaquettes) along the legs;
- the rungs have hopping $J$;
- there is a flux $\phi = \pi$ per plaquette.

At fixed particle number, writing the interaction as $\hat n^2$ instead of $\hat n(\hat n - 1)$ only shifts the energy by a constant, which cancels in the charge gap.

## What the code computes

For each value of $U$, `observables_calc.jl` does the following:

1. **Selects a chirality.** It runs DMRG at flux $\pi + \delta\phi$, which breaks the degeneracy between the two chiral states and selects one of them.
2. **Finds the ground state.** It runs DMRG at flux $\pi$, starting from that state, increasing the bond dimension step by step.
3. **Measures the observables:**
   - the loop current on each strong plaquette, $\langle\hat L_z\rangle$, together with the current on each bond;
   - the correlations $\langle\hat L_{z,0}\hat L_{z,i}\rangle$;
   - the fidelities $|\langle\psi(U)|\psi(U\pm\delta)\rangle|$, which give the fidelity susceptibility;
   - optionally, the ground state energies with $N\pm1$ particles, for the charge gap $E(N+1)+E(N-1)-2E(N)$;
   - optionally, the entanglement entropy profile.

## Structure

```
├── observables_calc.jl      # main script: ground state + observables
├── input.jl                 # parameters and their default values
├── src/
│   ├── hamiltonian_piflux_ladder.jl   # Hamiltonian and initial product states
│   ├── operators_piflux_ladder.jl     # loop current operators
│   ├── loop_current_calc.jl
│   ├── loop_current_correlator.jl
│   ├── fidelity_calc.jl
│   ├── compressibility_calc.jl        # energies with N +- 1 particles
│   └── entanglement_entropy.jl
├── scripts/
│   ├── submit.sh            # scan in U
│   └── excited_states.jl    # ground state and first excited states
├── test/
│   └── check_ed.jl          # DMRG vs exact diagonalization on a small ladder
├── analysis/
│   └── dmrg_analysis.ipynb  # fidelity susceptibility, loop current, charge gap, correlations, entropy
└── data/                    # results read by the analysis notebook
```

## Usage

The code needs the Julia packages `ITensors`, `ITensorMPS` and `ArgParse`. From the repository folder, check that DMRG reproduces exact diagonalization on a small ladder:

```bash
julia test/check_ed.jl
```

Run one point of the phase diagram:

```bash
julia observables_calc.jl --L 20 --U 1.0 --nu 0.25 --outdir data/results_nu0.25/Lx20
```

Or run a scan in $U$. The parameters are at the top of the script:

```bash
bash scripts/submit.sh $(which julia)
```

All the options are listed in `input.jl`. For example, `--compressibility true` also computes the charge gap, and `--entropy true` computes the entanglement entropy.

The results are written in `data/results_nu{nu}/Lx{L}/`, where the analysis notebook `analysis/dmrg_analysis.ipynb` reads them. The analysis is done in Python (NumPy, SciPy, Matplotlib).

## Contact

Maria Lanaro, University of Padova
