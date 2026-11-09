function compute_compressibility(
    psi, J, Jp, U, L, phi, deltaPhi, state, energy, Ns, ncut, qns;
    nsweeps0=10, nsweeps=10, maxdim=1000, delta_maxdim=200, max_maxdim=2000,
    cutoff=1e-10, noise=0.0, observer=nothing,
    output_dir="", ID="")

    # Ground state energies with one particle more (p) and one less (h).
    # The charge gap is E(N+1) + E(N-1) - 2E(N).
    # As for N particles, we first converge with flux phi + deltaPhi, to select
    # one chirality, and then at flux phi.

    # add a particle on the first site that can take one more, remove one from the first occupied site
    state_p, state_h = copy(state), copy(state)
    i_p = findfirst(s -> parse(Int, s) < ncut, state)
    i_h = findfirst(s -> parse(Int, s) > 0, state)
    state_p[i_p] = string(parse(Int, state[i_p]) + 1)
    state_h[i_h] = string(parse(Int, state[i_h]) - 1)
    println("state p", state_p, '\n')
    println("state h", state_h, '\n')

    sites_p = siteinds("Boson", Ns; dim=ncut+1, conserve_qns=qns)
    sites_h = siteinds("Boson", Ns; dim=ncut+1, conserve_qns=qns)
    psi0_p = productMPS(sites_p, state_p)
    psi0_h = productMPS(sites_h, state_h)

    ampo = Hamiltonian(J, Jp, U, L, phi+deltaPhi)
    H_p = MPO(ampo, sites_p);
    H_h = MPO(ampo, sites_h);

    # a fresh copy of the observer for each DMRG run
    new_obs() = observer === nothing ? nothing : deepcopy(observer)

    energy_p, psi00_p = dmrg(H_p, psi0_p; nsweeps=nsweeps0, maxdim=maxdim, cutoff=cutoff, noise=noise, observer=new_obs())
    energy_h, psi00_h = dmrg(H_h, psi0_h; nsweeps=nsweeps0, maxdim=maxdim, cutoff=cutoff, noise=noise, observer=new_obs())

    # Build dmrg sweep schedule
    function build_maxdim(psi00)
        mdim = []
        current = maximum(linkdims(psi00)) + delta_maxdim
        while current <= max_maxdim
            push!(mdim, current, current)
            current += delta_maxdim
        end;
        push!(mdim, max_maxdim)
        mdim
    end;
    maxdim_dmrg_p = build_maxdim(psi00_p)
    maxdim_dmrg_h = build_maxdim(psi00_h)

    ampo2 = Hamiltonian(J, Jp, U, L, phi)
    H2_p = MPO(ampo2, sites_p)
    H2_h = MPO(ampo2, sites_h)

    println("DMRG for n+1 particles")
    energy_p, psi_p = dmrg(H2_p, psi00_p; nsweeps=nsweeps, maxdim=maxdim_dmrg_p, cutoff=cutoff, noise=noise, observer=new_obs())

    println("DMRG for n-1 particles")
    energy_h, psi_h = dmrg(H2_h, psi00_h; nsweeps=nsweeps, maxdim=maxdim_dmrg_h, cutoff=cutoff, noise=noise, observer=new_obs())

    # check that the three states really have N, N+1 and N-1 bosons
    N0 = sum(expect(psi, "N"))
    N_p = sum(expect(psi_p, "N"))
    N_h = sum(expect(psi_h, "N"))
    @printf("Number of bosons: N = %.6f   N+1 state: %.6f   N-1 state: %.6f\n", N0, N_p, N_h)

    compressibility_dir = joinpath(output_dir, "compressibility")
    mkpath(compressibility_dir)
    filename = joinpath(compressibility_dir, "comp_" * ID * ".txt")

    open(filename, "w") do file
        write(file, "# E Ep Eh\n")
        writedlm(file, [energy energy_p energy_h])
    end;

    println("Energies for the charge gap saved to $filename")
end;
