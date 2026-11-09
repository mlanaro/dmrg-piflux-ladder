function Hamiltonian(J1, J2, U, L, phi)

    # Bose-Hubbard model on a two-leg ladder with alternating hopping J1 (strong
    # plaquettes) and J2 (weak plaquettes), and a flux phi per plaquette.
    # Site labelling:
    # 2 4 6 ..
    # 1 3 5 ..
    #
    # The interaction is U/2 n^2 instead of U/2 n(n-1): at fixed particle number
    # the difference is the constant -U N/2, which cancels in the charge gap.

    ampo = OpSum()

    # Horizontal hopping lower
    for j = 1:(L-1)
        if j%2 == 0
            J = J2
        else
            J = J1
        end
        ampo += -J, "Adag", 2*j-1, "A", 2*j+1
        ampo += -J, "A", 2*j-1, "Adag", 2*j+1
    end

    # Horizontal hopping upper, with the Peierls phase
    for j = 1:(L-1)
        if j % 2 == 0
            J = J2
        else
            J = J1
        end
        ampo += -J * exp(im * phi), "Adag", 2*j, "A", 2*j+2
        ampo += -J * exp(-im * phi), "A", 2*j, "Adag", 2*j+2
    end

    # Vertical hopping
    for j = 1:L
        ampo += -J1, "Adag", 2*j-1, "A", 2*j
        ampo += -J1, "A", 2*j-1, "Adag", 2*j
    end

    # Interaction
    for j = 1:2*L
        ampo += 0.5*U, "N * N", j
    end

    return ampo

end;


function initial_state(Ns, nu)

    # Product state with filling nu (bosons per site), repeating the same
    # pattern on every plaquette of 4 sites (1 = occupied, 0 = empty).
    patterns = Dict(
        0.25 => ["0", "0", "0", "1"],   # 0001-0001
        0.5  => ["1", "0", "0", "1"],   # 1001-1001
        0.75 => ["1", "1", "0", "1"],   # 1101-1101
        1.0  => ["1", "1", "1", "1"],
    )
    if !haskey(patterns, nu)
        error("Initial state not defined for nu = $nu")
    end
    return [patterns[nu][mod1(i, 4)] for i in 1:Ns]

end;
