# Phase273: restoring timescales and long-wavelength force ratio of the declared background (4x runtime modules).
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-refined268-runtime || exit 1
tr -d '\r' < $S/phase273-modes.py > .phase273-modes.py
taskset -c 3,5 python3 .phase273-modes.py .phase273-modes.json > .phase273-modes.log 2>&1; echo "exit=$?"; grep -v 'screen size' .phase273-modes.log | head -80
