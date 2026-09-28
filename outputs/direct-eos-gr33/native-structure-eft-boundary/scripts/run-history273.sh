# Phase273: charge history q(t) of the 4x solution for the masks all/state/geometry (128 times each, in parallel).
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-refined268-runtime || exit 1
tr -d '\r' < $S/phase273-history.py > .phase273-history.py
i=0; for mask in all state geometry; do
  taskset -c $((13 + i)) python3 .phase251-reader-launch.py .phase273-history.py readout268-quad64-work $mask 128 > .phase273-history-$mask.log 2>&1 &
  i=$((i + 1))
done
wait; echo done > .phase273-history.done
