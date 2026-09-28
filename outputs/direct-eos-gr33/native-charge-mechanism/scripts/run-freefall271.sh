# Phase271 runner: bash run-freefall271.sh <runtime> <mode> <readout folder> <label>
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd "$1" || exit 1
tr -d '\r' < $S/phase271-freefall.py > .phase271-freefall.py
taskset -c 9 python3 .phase251-reader-launch.py .phase271-freefall.py "$2" "$3" "$4" $5 2>&1 | grep -v 'screen size'
