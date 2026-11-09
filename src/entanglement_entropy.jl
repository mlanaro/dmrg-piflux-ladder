function entanglement_entropy(psi::MPS)

    # Von Neumann entropy S(b) for a cut between site b and b+1

    N = length(psi)
    S = zeros(Float64, N - 1)
    for b in 2:(N - 1)
        psi = orthogonalize(psi, b)
        lidx = linkind(psi, b - 1)
        sidx = siteind(psi, b)
        _, Svals = svd(psi[b], (lidx, sidx))

        s = Float64[]
        for (i, val) in enumerate(values(Svals))
            push!(s, abs(val))
        end

        p = s .^ 2
        p ./= sum(p)
        S[b] = -sum(p .* log.(p .+ 1e-12))
    end
    return S
end


function compute_entropy(psi; output_dir="", ID="")
    S = entanglement_entropy(psi)
    entropy_dir = joinpath(output_dir, "entropy")
    mkpath(entropy_dir)
    filename = joinpath(entropy_dir, "entropy_" * ID * ".txt")
    open(filename, "w") do file
        write(file, "# site entropy\n")
        for i in 1:length(S)
            @printf(file, "%d %.12f\n", i, S[i])
        end
    end
    println("Entanglement entropy saved to $filename")
end
