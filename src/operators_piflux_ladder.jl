function local_loop_current(j)
    
    # Calculates the loop current in the j-th strong-bond plaquettes (J1) 
    # Current is defined counter-clockwise: 1 -> 3 -> 4 -> 2 -> 1

    # 2 4 6 ..
    # 1 3 5 ..

    # j must be 1, 5, 9, ...

    if mod(j-1, 4) != 0
        print("Plaquette index ill defined.")
    end

    ampo = OpSum()
    ampo1 = OpSum()
    ampo2 = OpSum()
    ampo3 = OpSum()
    ampo4 = OpSum()

    # Horizontal hopping 1 -> 3
    ampo1 += im, "Adag", j, "A", j+2
    ampo1 += -im, "A", j, "Adag", j+2

    # Vertical hopping 3 -> 4
    ampo2 += im, "Adag", j+2, "A", j+3
    ampo2 += -im, "A", j+2, "Adag", j+3    

    # Horizontal hopping 4 -> 2
    ampo3 += -im, "Adag", j+3, "A", j+1
    ampo3 += im, "A", j+3, "Adag", j+1

    # Vertical hopping 2 -> 1
    ampo4 += im, "Adag", j+1, "A", j
    ampo4 += -im, "A", j+1, "Adag", j
    
    ampo = ampo1 + ampo2 + ampo3 + ampo4
    return ampo, ampo1, ampo2, ampo3, ampo4

end;


function loop_current_correlator(j,k)
    
    # Calculates the correlators of two loop currents on different plaquettes
    # Current is defined counter-clockwise: 1 -> 3 -> 4 -> 2 -> 1
    # j is the lower left site of the central plaquette
    # k is the lower left site of the displaced plaquette

    # 2 4 6 ..
    # 1 3 5 ..

    # j,k must be 1, 5, 9, ...
    

    if mod(j-1, 4) != 0
        print("Plaquette index ill defined.")
    end

    if mod(k-1, 4) != 0
        print("Plaquette index ill defined.")
    end

    ampo_j = OpSum()
    ampo_k = OpSum()
    ampo = OpSum()
    
    #Central plaquette
    # Horizontal hopping 1 -> 3
    ampo_j += im, "Adag", j, "A", j+2
    ampo_j += -im, "A", j, "Adag", j+2

    # Vertical hopping 3 -> 4
    ampo_j += im, "Adag", j+2, "A", j+3
    ampo_j += -im, "A", j+2, "Adag", j+3    

    # Horizontal hopping 4 -> 2
    ampo_j += -im, "Adag", j+3, "A", j+1
    ampo_j += im, "A", j+3, "Adag", j+1

    # Vertical hopping 2 -> 1
    ampo_j += im, "Adag", j+1, "A", j
    ampo_j += -im, "A", j+1, "Adag", j


    # Displaced plaquette

    # Horizontal hopping 1 -> 3
    ampo_k += im, "Adag", k, "A", k+2
    ampo_k += -im, "A", k, "Adag", k+2

    # Vertical hopping 3 -> 4
    ampo_k += im, "Adag", k+2, "A", k+3
    ampo_k += -im, "A", k+2, "Adag", k+3    

    # Horizontal hopping 4 -> 2
    ampo_k += -im, "Adag", k+3, "A", k+1
    ampo_k += im, "A", k+3, "Adag", k+1

    # Vertical hopping 2 -> 1
    ampo_k += im, "Adag", k+1, "A", k
    ampo_k += -im, "A", k+1, "Adag", k
    
    return ampo_j, ampo_k 

end;
