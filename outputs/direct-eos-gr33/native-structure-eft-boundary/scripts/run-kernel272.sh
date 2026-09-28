# Phase272: face-density kernel on the 4x grid (three face ranges and the full interior charge in parallel), then analysis.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-refined268-runtime || exit 1
tr -d '\r' < $S/phase272-kernel.py > .phase272-kernel.py
O=readout268-quad64-work
run() { taskset -c $1 python3 .phase251-reader-launch.py .phase272-kernel.py "${@:2}" > .phase272-$2-$4.log 2>&1; }
run 7 faces $O quad64 0 14 & run 9 faces $O quad64 15 29 & run 11 faces $O quad64 30 43 & run 13 full $O quad64 &
wait
python3 .phase251-reader-launch.py .phase272-kernel.py analyze $O quad64 > .phase272-analyze.log 2>&1
echo "exit=$?" > .phase272.done
