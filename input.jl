using ArgParse

settings = ArgParseSettings()
@add_arg_table! settings begin

# ---- Model ----
  "--L", "-L"
  help = "Ladder length (number of rungs, must be even)"
  arg_type = Int
  default = 106

  "--J"
  help = "Hopping amplitude J (strong plaquettes)"
  arg_type = Float64
  default = 1.0

  "--Jp"
  help = "Hopping amplitude J' (weak plaquettes)"
  arg_type = Float64
  default = 0.1

  "--U"
  help = "Interaction strength"
  arg_type = Float64
  default = 0.5

  "--nu"
  help = "Filling (bosons per site): 0.25, 0.5, 0.75 or 1"
  arg_type = Float64
  default = 0.25

  "--phi"
  help = "Flux per plaquette in units of the flux quantum (0.5 = pi flux)"
  arg_type = Float64
  default = 0.5

  "--deltaPhi"
  help = "Flux perturbation used in the first DMRG to select one chirality"
  arg_type = Float64
  default = 0.001

  "--ncut"
  help = "Maximum number of bosons per site"
  arg_type = Int
  default = 5

  "--qns"
  help = "Quantum number conservation"
  arg_type = Bool
  default = true

  "--outdir"
  help = "Output directory: results are saved in outdir/results_nu{nu}/Lx{L}"
  arg_type = String
  default = "data"

# ---- Observables to compute ----
  "--correlator"
  help = "Compute the loop current correlations"
  arg_type = Bool
  default = true

  "--fidelity"
  help = "Compute the fidelity"
  arg_type = Bool
  default = true

  "--compressibility"
  help = "Compute the energies with N+1 and N-1 particles (charge gap)"
  arg_type = Bool
  default = true

  "--entropy"
  help = "Compute the entanglement entropy"
  arg_type = Bool
  default = false

# ---- DMRG ----
  "--nsweeps0"
  help = "Number of sweeps for the first DMRG with explicitly broken symmetry"
  arg_type = Int64
  default = 100

  "--nsweep", "-n"
  help = "Number of sweeps of DMRG"
  arg_type = Int64
  default = 500

  "--max_maxdim"
  help = "Max bond dimension last sweep in DMRG"
  arg_type = Int64
  default = 500

  "--maxdim"
  help = "Maximum bond dimension in DMRG (used for the charge gap)"
  nargs = '+'
  arg_type = Int64
  default = [50, 100, 150, 150, 150, 150, 200, 200, 200, 200, 200, 300, 300, 300, 300, 300, 300, 400, 400, 500, 500]

  "--delta_maxdim"
  help = "Max bond dimension steps in DMRG sweeps"
  nargs = '+'
  arg_type = Int64
  default = [20]

  "--cutoff", "-c"
  help = "Truncation cutoff in DMRG"
  nargs = '+'
  arg_type = Float64
  default = [1.e-8, 1e-8, 1e-9, 1e-9, 1e-10]

  "--noise"
  help = "Magnitude of noise to add in DMRG"
  nargs = '+'
  arg_type = Float64
  default = [1e-2, 1e-2, 1e-4, 1e-4, 1e-8, 0]

  "--enConv"
  help = "Stop DMRG when the energy is converged"
  arg_type = Bool
  default = true

  "--etol"
  help = "Tolerance on energy"
  arg_type = Float64
  default = 1E-8

end