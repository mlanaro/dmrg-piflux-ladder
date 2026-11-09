function compute_4body_correlator(sites, psi, L, J, Jp, U; output_dir="", ID="")

    # Correlations <Lz_j Lz_k> between the loop current of the central plaquette j
    # and the one of the plaquettes k on its right.
    # The loop current is Hermitian, so <psi|Lz_j Lz_k|psi> = <Lz_j psi|Lz_k psi>:
    # inner(Lzj, psi, Lzk, psi) computes it exactly, without truncating
    # the MPS after applying the operators.

    Ns = 2*L
    if isodd(L ÷ 2)
        index_central_plaq = (Ns ÷ 2 - 1)
    else
        index_central_plaq = (Ns ÷ 2 - 3)
    end

    println("index_central_plaq = $index_central_plaq")
    Lzvals_total = ComplexF64[]
    k_vector = collect(index_central_plaq : 4 : Ns - 3)
    for k in k_vector
        ampo_j, ampo_k = loop_current_correlator(index_central_plaq, k)
        Lzj = MPO(ampo_j, sites)
        Lzk = MPO(ampo_k, sites)

        val = inner(Lzj, psi, Lzk, psi)
        push!(Lzvals_total, val)
    end
    lz_dir = joinpath(output_dir, "Lz_corr")
    mkpath(lz_dir)
    filename = joinpath(lz_dir, "Lz_corr" * ID * ".txt")

    open(filename, "w") do file
        write(file, "# Loop current correlations with central plaquette $index_central_plaq\n")
        write(file, "# Format: k_plaquette_index\tRe(Lz)\tIm(Lz)\n")
        for (k, val) in zip(k_vector, Lzvals_total)
            write(file, "$(k)\t$(real(val))\t$(imag(val))\n")
        end;
    end;

    println("Loop current correlators saved to $filename")

end;
