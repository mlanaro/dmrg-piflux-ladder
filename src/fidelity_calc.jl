function fidelity(sites, psi, psi_DMRG, J1, J2, U, L, delta, phi; nsweeps=2000, max_maxdim=1600, delta_maxdim=[100], cutoff=[1e-10], observer=nothing)

    # Ground state at U + delta, starting from psi, and its overlap with psi_DMRG

    ampo = Hamiltonian(J1,J2,U+delta,L,phi)
    H = MPO(ampo, sites)

    maxdim_psi00 = maximum(linkdims(psi))
    maxdim_dmrg = []
    current = maxdim_psi00
    for Δ in delta_maxdim
        next_val = current + Δ
        if next_val > max_maxdim
            break
        end
        push!(maxdim_dmrg, next_val)
        push!(maxdim_dmrg, next_val)
        current = next_val
    end
    push!(maxdim_dmrg, max_maxdim)

    @show maxdim_psi00
    @show maxdim_dmrg

    # a fresh copy of the observer for each DMRG run, so that the energy
    # convergence check does not use the energies of a previous run
    obs = observer === nothing ? nothing : deepcopy(observer)
    energy, psi1 = dmrg(H, psi; nsweeps=nsweeps, maxdim=maxdim_dmrg, cutoff=cutoff, observer=obs)
    overlap = abs(inner(psi1,psi_DMRG))
    return overlap, energy, psi1

end;


function compute_fidelity(sites, psi, J, Jp, U, L, phi;
                          nsweeps=2000, max_maxdim=500, delta_maxdim=[20], cutoff=[1e-11], delta0=0.005,
                          observer=nothing,
                          psi_DMRG=psi, output_dir="", ID="")

    # Fidelities |<psi(U)|psi(U +- delta)>| and |<psi(U+delta)|psi(U-delta)>|.
    # The fidelity susceptibility is then
    #   chi_F = -(F_+ + F_- - 2) / (L delta^2)

    delta_list = [delta0]
    overlap_list_plus = zeros(length(delta_list))
    overlap_list_minus = zeros(length(delta_list))
    overlap_list_sim = zeros(length(delta_list))
    energy_list_plus = zeros(length(delta_list))
    energy_list_minus = zeros(length(delta_list))

    psi0_plus = psi
    psi0_minus = psi

    for i in eachindex(delta_list)
        local_delta_plus = delta_list[i]
        local_delta_minus = -delta_list[i]

        println("Computing fidelity for +δ, -δ = $local_delta_plus, $local_delta_minus")

        overlap_plus, energy_plus, psi_plus = fidelity(
            sites, psi0_plus, psi_DMRG, J, Jp, U, L, local_delta_plus, phi;
            nsweeps=nsweeps, max_maxdim=max_maxdim, delta_maxdim=delta_maxdim,
            cutoff=cutoff, observer=observer)

        overlap_minus, energy_minus, psi_minus = fidelity(
            sites, psi0_minus, psi_DMRG, J, Jp, U, L, local_delta_minus, phi;
            nsweeps=nsweeps, max_maxdim=max_maxdim, delta_maxdim=delta_maxdim,
            cutoff=cutoff, observer=observer)

        overlap_list_plus[i] = overlap_plus
        overlap_list_minus[i] = overlap_minus
        overlap_list_sim[i] = abs(inner(psi_plus, psi_minus))
        energy_list_plus[i] = energy_plus
        energy_list_minus[i] = energy_minus

        psi0_plus = psi_plus
        psi0_minus = psi_minus
    end;

    fid_dir = joinpath(output_dir, "fidelity")
    mkpath(fid_dir)
    filename = joinpath(fid_dir, "fidelity_" * ID * ".txt")

    open(filename, "w") do file
        write(file, "# delta_list overlap_plus overlap_minus overlap_sim GS_energy_plus GS_energy_minus\n")
        writedlm(file, [delta_list overlap_list_plus overlap_list_minus overlap_list_sim energy_list_plus energy_list_minus])
    end;

    println("Fidelity results saved to $filename")
end;
