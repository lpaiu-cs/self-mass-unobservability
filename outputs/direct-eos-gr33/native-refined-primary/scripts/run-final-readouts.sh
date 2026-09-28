# Phase267 final readouts: refined charge at t62, t63, T (+ depth bands at T); original charge at t62, t63.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
R=/home/lpaiu/work/native-refined267-runtime; G=/home/lpaiu/work/native-retained-tail-runtime
for d in $R $G; do for f in phase267-audited phase267-readout phase267-depth; do tr -d '\r' < $S/$f.py > $d/.$f.py; done; done
W=$R/primary267-refined64-work
refined() {  # $1 = segment label end (62|63|64), $2 = core
  N=$1; C=$2; cd $R; O=readout267-refined$N-work; L=.phase267-readout-refined$N
  r() { taskset -c $C python3 .phase267-audited.py $L-$1-record.json .phase251-reader-launch.py .phase267-readout.py "$@" > $L-$1.stdout.log 2> $L-$1.stderr.log || { echo "failed $1" > $L.done; return 1; }; }
  r endpoints $O primary267-refined64-work/sweep-1/photons/seg-$N.npz captures:primary267-refined64-work/captures && r source $O && r field $O || return 1
  if [ $N = 64 ]; then
    bands="cells:0-7"; for c in $(seq 8 26); do bands="$bands cells:$c-$c"; done; bands="$bands cells:27-154 cells:155-282 cells:283-410 cells:411-538 boundary all"
    taskset -c $C python3 .phase251-reader-launch.py .phase267-depth.py $O bands $bands > $L-depth.stdout.log 2> $L-depth.stderr.log
  fi
  echo ok > $L.done
}
refined 63 9 & refined 64 11 &
( cd $G
  for n in 63; do
    t=$(cd $R && python3 -c "import numpy as np; print(repr(float(np.load('primary267-refined64-work/sweep-1/photons/seg-$n.npz')['actual_step_edges'][-1])))")
    taskset -c 13 python3 .phase251-reader-launch.py .phase267-readout.py field readout267-identity2-work $t > .phase267-t$n-field.stdout.log 2> .phase267-t$n-field.stderr.log
  done; echo ok > .phase267-t63.done ) &
wait
echo done > $R/.phase267-final-readouts.done
