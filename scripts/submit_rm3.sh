#!/bin/bash
#PBS -N MPI_job
#PBS -l nodes=1:ppn=69
#PBS -j oe
#PBS -q al9_mpi
#PBS -l walltime=240:00:00

# Match the existing cluster environment before enabling strict shell checks.
if [[ -r /storage/local/home/theorm3/silvestrini/.bashrc_ui02 ]]; then
    source /storage/local/home/theorm3/silvestrini/.bashrc_ui02
fi
set -euo pipefail
: "${PBS_NODEFILE:?Submit through PBS so CPU slots are known}"

NCPU=$(awk 'NF { n++ } END { print n+0 }' "$PBS_NODEFILE")
NNODES=$(awk 'NF { hosts[$1]=1 } END { for (h in hosts) n++; print n+0 }' "$PBS_NODEFILE")
if [[ "$NCPU" -lt 1 ]]; then
    echo "PBS_NODEFILE contains no CPU slots" >&2
    exit 1
fi

WORKDIR=${GOLDENMODES_WORKDIR:-${PBS_O_WORKDIR:?Set PBS_O_WORKDIR or GOLDENMODES_WORKDIR}}
MPIEXEC=${GOLDENMODES_MPIEXEC:-/storage/theor/Users/silvestrini/Al9/bin/mpiexec}
BAT_LIBDIR=${GOLDENMODES_BAT_LIBDIR:-/storage/theor/Users/silvestrini/Al9/BAT_parallel/lib}
cd "$WORKDIR"
export LD_LIBRARY_PATH="$BAT_LIBDIR:${LD_LIBRARY_PATH:-}"
export GOLDENMODES_MPI_EXPECTED_RANKS="$NCPU"
# Each PBS CPU slot runs one MPI process. Avoid nested BLAS/OpenMP threads.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

# PBS slots may be SMT/hyperthread CPUs. Bind to individual hardware threads,
# without permitting two ranks to share the same processing unit by force.
MPI_PLACEMENT=(--use-hwthread-cpus --map-by hwthread --bind-to hwthread --report-bindings)

if [[ "$#" -gt 0 ]]; then
    arguments=("$@")
else
    # ARGS is a whitespace-separated argument list, as in the existing qsub jobs.
    # No eval or wildcard expansion is performed. Use positional arguments for
    # paths with spaces when invoking this script directly for a local test.
    read -r -a arguments <<< "${ARGS:-}"
fi
echo "Launching $NCPU MPI processes on $NNODES node(s)"
echo "MPI placement: ${MPI_PLACEMENT[*]}"
if command -v lscpu >/dev/null 2>&1; then
    LC_ALL=C lscpu | awk -F: '/^CPU\(s\):|^On-line CPU|^Thread\(s\) per core:|^Core\(s\) per socket:|^Socket\(s\):/ {print}'
fi
if [[ -r /proc/self/status ]]; then
    awk '/^Cpus_allowed_list:/ {print "Launcher " $0}' /proc/self/status
fi
exec "$MPIEXEC" -machinefile "$PBS_NODEFILE" -n "$NCPU" "${MPI_PLACEMENT[@]}" "$WORKDIR/mainGM" "${arguments[@]}"
