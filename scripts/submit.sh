#!/bin/bash
# Scan in U: one run of observables_calc.jl for each value of U.
# All the other parameters (L, nu, J, Jp, bond dimension, observables, ...)
# are taken from input.jl.
#
# Usage:  bash scripts/submit.sh /path/to/julia

export JULIA_NUM_THREADS=12 # only when one process is running
julia_exe=$1   # first argument is Julia executable path

repo_dir="$(cd "$(dirname "$0")/.." && pwd)"

# Values of U
Ustart=0.4
Ustop=0.9
Usteps=11

Ulist=$(python3 - <<END
import numpy as np
Ulist = np.linspace($Ustart, $Ustop, $Usteps)
print(' '.join(map(str, Ulist)))
END
)

Uarray=($Ulist)
echo "Ulist: $Ulist"

# log files of the runs
logdir="$repo_dir/data/logs"
mkdir -p "$logdir"
echo "Log files in: $logdir"

max_procs=1  # number of parallel Julia processes

for U in "${Uarray[@]}"
do
    logfile="$logdir/output_U${U}_$(date '+%Y%m%d_%H%M%S').txt"

    {
        echo "Started at: $(date '+%Y-%m-%d %H:%M:%S')";

        "$julia_exe" "$repo_dir/observables_calc.jl" --U $U --outdir "$repo_dir/data";

        echo "Finished at: $(date '+%Y-%m-%d %H:%M:%S')";
    } > "$logfile" 2>&1 &

    # Wait if we reached max parallel jobs
    while [ $(jobs -r | wc -l) -ge $max_procs ]; do
        sleep 1
    done
done

wait
echo "All jobs finished."