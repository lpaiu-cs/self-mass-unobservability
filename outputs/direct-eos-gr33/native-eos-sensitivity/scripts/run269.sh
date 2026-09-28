# Phase269 stage 1 runner: bash run269.sh <args for phase269-eos.py> (runs in the original runtime, one build per process).
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
tr -d '\r' < $S/phase269-eos.py > .phase269-eos.py
taskset -c 7 python3 .phase269-eos.py "$@"
