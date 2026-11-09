# excited_states.jl
#
# Ground state and first two excited states of a short ladder with DMRG:
# each new state is found with a penalty (weight) on the overlap with the
# previous ones.
# Usage:  julia scripts/excited_states.jl

using ITensorMPS, ITensors, Printf, LinearAlgebra

include("../src/hamiltonian_piflux_ladder.jl")

Nthreads = 1
BLAS.set_num_threads(Nthreads)

# Parameters
Lx = 10;
Ns = 2*Lx;
J1 = 1.;
J2 = 0.1;
U = 0.8;
phi = pi;
nu = 0.5;
ncut = 3;
nsweeps = 30;
maxdim = 600;
error = 1E-10;

sites = siteinds("Boson",Ns,dim=ncut+1,conserve_qns=true);
ampo = Hamiltonian(J1,J2,U,Lx,phi)
H = MPO(ampo, sites);

state = initial_state(Ns, nu)
psi0 = productMPS(sites,state);

print("Finding ground state Lx = $(Lx), J1 = $(J1), J2 = $(J2), U = $(U)\n\n")

sweeps = Sweeps(nsweeps)
setmaxdim!(sweeps,10,100,200,maxdim)
setcutoff!(sweeps,error)
energy1, psi1 = dmrg(H,psi0,sweeps)

print("\nFirst excited state\n")
energy2, psi2 = dmrg(H,[psi1],psi0,sweeps;weight = J1)

print("\nSecond excited state\n")
energy3, psi3 = dmrg(H,[psi1, psi2],psi0,sweeps;weight = J1)

print("\nCheck orthogonality\n")
@show inner(psi1, psi2)
@show inner(psi1, psi3)
@show inner(psi2, psi3)

@printf("\nE0 = %.10f\nE1 = %.10f\nE2 = %.10f\n", energy1, energy2, energy3)
print("\n Simulation done. \n")
