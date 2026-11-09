function compute_loop_current(sites, psi, L, J, Jp, U; output_dir="", ID="")
    Ns = 2*L
    Lz_list_total = [MPO(local_loop_current(j)[1], sites) for j in 1:4:Ns]
    Lz_list_1 = [MPO(local_loop_current(j)[2], sites) for j in 1:4:Ns]
    Lz_list_2 = [MPO(local_loop_current(j)[3], sites) for j in 1:4:Ns]
    Lz_list_3 = [MPO(local_loop_current(j)[4], sites) for j in 1:4:Ns]
    Lz_list_4 = [MPO(local_loop_current(j)[5], sites) for j in 1:4:Ns]

    Lzvals_total = [inner(psi', Lzj, psi) for Lzj in Lz_list_total]
    Lzvals_1 = [inner(psi', Lzj, psi) for Lzj in Lz_list_1]
    Lzvals_2 = [inner(psi', Lzj, psi) for Lzj in Lz_list_2]
    Lzvals_3 = [inner(psi', Lzj, psi) for Lzj in Lz_list_3]
    Lzvals_4 = [inner(psi', Lzj, psi) for Lzj in Lz_list_4]

    lz_dir = joinpath(output_dir, "Lz")
    mkpath(lz_dir)
    filename = joinpath(lz_dir, "Lz_" * ID * ".txt")

    open(filename, "w") do file
        write(file, "# Loop current expectation values per plaquette\n")
        write(file, "# Format: plaquette_index\tRe(Lz)\tIm(Lz)\n")
        for (i, val) in enumerate(Lzvals_total)
            write(file, "$(i)\t$(real(val))\t$(imag(val))\n")
        end;
        write(file, "# Loop current expectation values Horizontal hopping 1 -> 3\n")
        for (i, val) in enumerate(Lzvals_1)
            write(file, "$(i)\t$(real(val))\t$(imag(val))\n")
        end;
        write(file, "# Loop current expectation values Vertical hopping 3 -> 4\n")
        for (i, val) in enumerate(Lzvals_2)
            write(file, "$(i)\t$(real(val))\t$(imag(val))\n")
        end;
        write(file, "# Loop current expectation values Horizontal hopping 4 -> 2\n")
        for (i, val) in enumerate(Lzvals_3)
            write(file, "$(i)\t$(real(val))\t$(imag(val))\n")
        end;
        write(file, "# Loop current expectation values Vertical hopping 2 -> 1\n")
        for (i, val) in enumerate(Lzvals_4)
            write(file, "$(i)\t$(real(val))\t$(imag(val))\n")
        end;
    end;

    println("Loop currents saved to $filename")
end;
