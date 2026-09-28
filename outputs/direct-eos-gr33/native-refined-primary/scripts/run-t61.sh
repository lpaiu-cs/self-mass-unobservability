# Phase267: compact charge at t61 (= end of the accepted refined macro 61) on both grids, with depth bands.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
R=/home/lpaiu/work/native-refined267-runtime; G=/home/lpaiu/work/native-retained-tail-runtime
for d in $R $G; do for f in phase267-audited phase267-readout phase267-depth; do tr -d '\r' < $S/$f.py > $d/.$f.py; done; done
[ -f $R/.phase251-reader-launch.py ] || cp $G/.phase251-reader-launch.py $R/
T61=$(cd $R && python3 -c "import numpy as np; print(repr(float(np.load('primary267-refined64-work/sweep-1/photons/seg-61.npz')['actual_step_edges'][-1])))")
echo "t61=$T61" > $R/.phase267-t61.txt
# original grid, same time, from the full accepted source (retarded field depends only on the past)
( cd $G
  taskset -c 9 python3 .phase251-reader-launch.py .phase267-readout.py field readout267-identity2-work $T61 > .phase267-t61-field.stdout.log 2> .phase267-t61-field.stderr.log
  bands=""; for c in $(seq 0 18); do bands="$bands cells:$c-$c"; done; bands="$bands cells:19-146 cells:147-274 cells:275-402 cells:403-530 boundary all"
  PHASE267_TIME=$T61 taskset -c 9 python3 .phase251-reader-launch.py .phase267-depth.py readout267-identity2-work t61 $bands > .phase267-t61-depth.stdout.log 2> .phase267-t61-depth.stderr.log
  echo "exit=$?" > .phase267-t61.done ) &
# refined grid: accepted history to macro 61 with its own captures
( cd $R
  O=readout267-refined61-work; L=.phase267-readout-refined61; W=primary267-refined64-work
  run() { taskset -c 7 python3 .phase267-audited.py $L-$1-record.json .phase251-reader-launch.py .phase267-readout.py "$@" > $L-$1.stdout.log 2> $L-$1.stderr.log || { echo "failed $1" > $L.done; exit 1; }; }
  run endpoints $O $W/sweep-1/photons/seg-61.npz captures:$W/captures && run source $O && run field $O || exit 1
  bands="cells:0-7"; for c in $(seq 8 26); do bands="$bands cells:$c-$c"; done; bands="$bands cells:27-154 cells:155-282 cells:283-410 cells:411-538 boundary all"
  taskset -c 7 python3 .phase251-reader-launch.py .phase267-depth.py $O bands $bands > $L-depth.stdout.log 2> $L-depth.stderr.log
  echo "exit=$?" > $L.done ) &
wait
echo done > $R/.phase267-t61-all.done
