# observables_calc.jl
#
# Ground state of the Bose-Hubbard pi-flux ladder with DMRG and its observables.
# Usage:
#   julia observables_calc.jl --U 1.0
# All the parameters and their default values are in input.jl: any of them can
# be changed from the command line, e.g. --L 40 --nu 1.
# The results are saved in outdir/results_nu{nu}/Lx{L}/ (outdir = data by default).

using ITensorMPS, ITensors, DelimitedFiles, Printf, LinearAlgebra

include("src/hamiltonian_piflux_ladder.jl")
include("src/operators_piflux_ladder.jl")
include("src/loop_current_calc.jl")
include("src/loop_current_correlator.jl")
include("src/fidelity_calc.jl")
include("src/compressibility_calc.jl")
include("src/entanglement_entropy.jl")

let

    include("input.jl")
    args = parse_args(settings) # parsing args in input.jl and ARGS by default

    L = args["L"]
    qns = args["qns"]
    ncut = args["ncut"]
    nu = args["nu"]

    J = args["J"]
    Jp = args["Jp"]
    U = args["U"]
    phi = 2*pi*args["phi"]
    deltaPhi = args["deltaPhi"]

    if isodd(L)
        error("L must be even: the ladder has to end with a full J plaquette")
    end

    println("DMRG on a ladder of 2x$L sites with state per sites: $ncut+1")
    Ns = 2*L

    sites = siteinds("Boson", Ns; dim=ncut+1, conserve_qns=qns)

    # parameters for DMRG sweeps
    nsweeps0 = args["nsweeps0"]
    nsweeps = args["nsweep"]
    print("Number of sweeps first DMRG: $(nsweeps0)\n")
    print("Number of sweeps second DMRG: $(nsweeps)\n")

    maxdim = args["maxdim"]
    delta_maxdim = args["delta_maxdim"]
    max_maxdim = args["max_maxdim"]
    cutoff = args["cutoff"]
    noise = args["noise"]
    etol_cutoff = 1e-11   # cutoff of the first DMRG, also used in the file names
    ID = @sprintf("L%d_J%.2f_Jp%.2f_U%.3f_maxdim%d_cutoff%.0e", L, J, Jp, U, max_maxdim, etol_cutoff)

    enConv = args["enConv"]
    energy_tol = args["etol"]

    # results in outdir/results_nu{nu}/Lx{L}/
    nu_str = isinteger(nu) ? string(Int(nu)) : string(nu)
    output_dir = joinpath(args["outdir"], "results_nu$(nu_str)", "Lx$(L)")
    mkpath(output_dir)
    println("Results will be saved in $output_dir")
    println("Parameters: L = $L, nu = $nu, J = $J, Jp = $Jp, U = $U, max_maxdim = $max_maxdim\n")

    # ---------------------------------------------------------------
    # 1) DMRG with flux phi + deltaPhi: the small perturbation breaks the
    #    degeneracy between the two chiralities and selects one of them
    # ---------------------------------------------------------------
    ampo = Hamiltonian(J, Jp, U, L, phi+deltaPhi)
    H = MPO(ampo, sites);

    state = initial_state(Ns, nu)
    println("state ", state, '\n')
    psi0 = productMPS(sites, state)

    println("Finding ground state L = $(L), J = $(J), Jp = $(Jp), U = $(U)")

    maxdim_sweeps0 = [10, 50, 50, 50, 100, 100, 150, 150, 200, 200, 300, 300, 400, 400, 500]
    cutoff_0 = [etol_cutoff]
    println("maxdim0: $maxdim_sweeps0")
    println("cutoff0: $cutoff_0")
    obs0 = enConv ? DMRGObserver(; energy_tol=1e-5, minsweeps=length(cutoff_0), energy_type=Float64) : nothing

    energy, psi00 = dmrg(H, psi0; nsweeps=nsweeps0, maxdim=maxdim_sweeps0, cutoff=cutoff_0, noise=noise, observer=obs0)
    println("Ground-state energy symmetry broken: $energy")

    println("Compute loop current \n")
    compute_loop_current(sites, psi00, L, J, Jp, U; output_dir=joinpath(output_dir, "Lz0"), ID=ID)
    if args["entropy"]
        compute_entropy(psi00; output_dir=joinpath(output_dir, "entropy0"), ID=ID)
    end

    # ---------------------------------------------------------------
    # 2) DMRG at flux phi, starting from the state with broken symmetry
    # ---------------------------------------------------------------
    maxdim_psi00 = maximum(linkdims(psi00))
    maxdim_dmrg = []
    current = maxdim_psi00
    for Δ in delta_maxdim
        next_val = current + Δ
        if next_val > max_maxdim
            break
        end
        push!(maxdim_dmrg, next_val)
        current = next_val
    end
    push!(maxdim_dmrg, max_maxdim)

    @show maxdim_psi00
    @show maxdim_dmrg
    println("cutoff DMRG: $cutoff")
    obs = enConv ? DMRGObserver(; energy_tol=energy_tol, minsweeps=length(noise), energy_type=Float64) : nothing
    ampo2 = Hamiltonian(J, Jp, U, L, phi)
    H2 = MPO(ampo2, sites);
    energy, psi = dmrg(H2, psi00; nsweeps=nsweeps, maxdim=maxdim_dmrg, cutoff=cutoff, noise=noise, observer=obs)
    println("Ground-state energy: $energy")

    # ---------------------------------------------------------------
    # 3) Observables
    # ---------------------------------------------------------------
    println("Compute loop current \n")
    compute_loop_current(sites, psi, L, J, Jp, U; output_dir=output_dir, ID=ID)

    if args["correlator"]
        println("Compute loop current correlator \n")
        compute_4body_correlator(sites, psi, L, J, Jp, U; output_dir=output_dir, ID=ID)
    end

    if args["entropy"]
        println("Compute entanglement entropy \n")
        compute_entropy(psi; output_dir=output_dir, ID=ID)
    end

    if args["fidelity"]
        println("Compute fidelity \n")
        obs_fidelity = enConv ? DMRGObserver(; energy_tol=energy_tol, minsweeps=length(noise), energy_type=Float64) : nothing
        compute_fidelity(sites, psi, J, Jp, U, L, phi;
                         nsweeps=2000, max_maxdim=max_maxdim,
                         delta_maxdim=delta_maxdim, cutoff=[1e-11],
                         observer=obs_fidelity, psi_DMRG=psi, output_dir=output_dir, ID=ID)
    end

    if args["compressibility"]
        println("Compute compressibility \n")
        obs_compressibility = enConv ? DMRGObserver(; energy_tol=1E-5, minsweeps=2, energy_type=Float64) : nothing
        compute_compressibility(
            psi, J, Jp, U, L, phi, deltaPhi, state, energy, Ns, ncut, qns;
            nsweeps=nsweeps, maxdim=maxdim, cutoff=cutoff, noise=noise, observer=obs_compressibility,
            output_dir=output_dir, ID=ID)
    end

end;

println("\nSimulation done. \n")