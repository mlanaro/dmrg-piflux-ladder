# check_ed.jl
#
# Check of DMRG against exact diagonalization on a small ladder (L = 4 rungs,
# 8 sites), at the fillings nu = 1/4 and nu = 1.
# The ED is written independently of the ITensor OpSum: the Hamiltonian is built
# directly in the basis of occupation numbers with fixed total number of bosons,
# following the same definitions of hamiltonian_piflux_ladder.jl.
# Usage (from the repository folder):  julia test/check_ed.jl

using ITensorMPS, ITensors, LinearAlgebra, Printf

include("../src/hamiltonian_piflux_ladder.jl")


function ed_energies(J, Jp, U, L, phi, ncut, Npart)
    Ns = 2*L
    d = ncut + 1

    # basis: occupations (n_1, ..., n_Ns) with n_i <= ncut and sum = Npart
    basis = [collect(Tuple(c)) .- 1 for c in CartesianIndices(ntuple(_ -> d, Ns)) if sum(Tuple(c)) - Ns == Npart]
    index = Dict(b => i for (i, b) in enumerate(basis))

    # hopping terms c b_i^dag b_j (the Hermitian conjugate is added below)
    hops = Tuple{ComplexF64, Int, Int}[]
    for j in 1:(L-1)
        Jj = iseven(j) ? Jp : J
        push!(hops, (-Jj, 2*j-1, 2*j+1))                  # lower leg
        push!(hops, (-Jj * exp(im*phi), 2*j, 2*j+2))      # upper leg, with the flux
    end
    for j in 1:L
        push!(hops, (-J, 2*j-1, 2*j))                     # rungs
    end

    H = zeros(ComplexF64, length(basis), length(basis))
    for (a, s) in enumerate(basis)
        H[a, a] += 0.5 * U * sum(s .^ 2)
        for (c, i, j) in hops
            for (cc, ii, jj) in ((c, i, j), (conj(c), j, i))
                (s[jj] == 0 || s[ii] == ncut) && continue
                t = copy(s)
                amp = sqrt(t[jj]); t[jj] -= 1
                amp *= sqrt(t[ii] + 1); t[ii] += 1
                H[index[t], a] += cc * amp
            end
        end
    end
    return length(basis), eigvals(Hermitian(H))
end


L = 4
Ns = 2*L
ncut = 2          # at most 2 bosons per site, the same for ED and DMRG
J = 1.0
Jp = 0.1
phi = pi

for nu in [0.25, 1.0]
    state = initial_state(Ns, nu)
    Npart = sum(parse.(Int, state))
    println("\nnu = $nu: $Npart bosons on $Ns sites")

    for U in [0.5, 1.0, 2.0, 5.0]
        dim, E_ed = ed_energies(J, Jp, U, L, phi, ncut, Npart)

        # The two chiral states can be split by a very small gap, so DMRG needs
        # noise at the beginning, many sweeps and an accurate local eigensolver.
        sites = siteinds("Boson", Ns; dim=ncut+1, conserve_qns=true)
        H = MPO(Hamiltonian(J, Jp, U, L, phi), sites)
        psi0 = productMPS(sites, state)
        E_dmrg, psi = dmrg(H, psi0; nsweeps=100, maxdim=[10, 20, 50, 100, 200], cutoff=1e-12,
                           noise=[1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 0.0],
                           eigsolve_krylovdim=10, eigsolve_maxiter=4, outputlevel=0)

        @printf("  U = %.1f   E0_ED = %.12f   E_DMRG = %.12f   difference = %.1e   (E1-E0 = %.4f, ED dim = %d)\n",
                U, E_ed[1], E_dmrg, abs(E_ed[1] - E_dmrg), E_ed[2] - E_ed[1], dim)
    end
end