# Phase268 readouts: 4x compact charge at t61, t62, t63 and T as each accepted segment appears (phase-267 readout path),
# plus depth bands at T (every refined sub-cell of original cells 8-15, the thin cells, atmosphere quarters, boundary).
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
R=/home/lpaiu/work/native-refined268-runtime; R2=/home/lpaiu/work/native-refined267-runtime
cd $R || exit 1
for f in phase267-audited phase267-readout phase267-depth; do tr -d '\r' < $S/$f.py > .$f.py; done
W=primary268-quad64-work; LOGS=.phase268-quad64-logs
readout() {  # $1 = macro (61..64), $2 = core
  N=$1; C=$2; O=readout268-quad$N-work; L=.phase268-readout-quad$N
  python3 -c "import numpy as np; a=float(np.load('$W/sweep-1/photons/seg-$N.npz')['actual_step_edges'][-1]); b=float(np.load('$R2/primary267-refined64-work/sweep-1/photons/seg-$N.npz')['actual_step_edges'][-1]); print(repr(a), repr(b), a == b)" > $L-time.txt
  r() { taskset -c $C python3 .phase267-audited.py $L-$1-record.json .phase251-reader-launch.py .phase267-readout.py "$@" > $L-$1.stdout.log 2> $L-$1.stderr.log || { echo "failed $1" > $L.done; return 1; }; }
  r endpoints $O $W/sweep-1/photons/seg-$N.npz captures:$W/captures && r source $O && r field $O || return 1
  if [ $N = 64 ]; then
    bands="cells:0-7"; for c in $(seq 8 42); do bands="$bands cells:$c-$c"; done; bands="$bands cells:43-170 cells:171-298 cells:299-426 cells:427-554 boundary all"
    taskset -c $C python3 .phase251-reader-launch.py .phase267-depth.py $O bands $bands > $L-depth.stdout.log 2> $L-depth.stderr.log || { echo "failed depth" > $L.done; return 1; }
  fi
  echo ok > $L.done
}
core=7
for N in 61 62 63 64; do
  until [ -f $W/seg-$N-driver.json ] && grep -q '"passed": true' $W/seg-$N-driver.json; do
    [ -f $LOGS/stopped.txt ] && { echo "run stopped before seg-$N" > .phase268-readouts.done; wait; exit 1; }
    sleep 60
  done
  echo "$(date '+%F %T') readout $N start (core $core)" >> .phase268-readouts.log
  ( readout $N $core; echo "$(date '+%F %T') readout $N end $(cat .phase268-readout-quad$N.done)" >> .phase268-readouts.log ) &
  core=$((core + 2))
done
wait
echo done > .phase268-readouts.done
