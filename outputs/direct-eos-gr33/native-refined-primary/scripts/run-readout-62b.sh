S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-refined267-runtime
tr -d '\r' < $S/phase267-readout.py > .phase267-readout.py
mv readout267-refined62-work readout267-refined62-failed-exact-inverse; for f in .phase267-readout-refined62*; do mv $f ${f/refined62/refined62-failed-exact-inverse}; done
O=readout267-refined62-work; L=.phase267-readout-refined62
r() { taskset -c 7 python3 .phase267-audited.py $L-$1-record.json .phase251-reader-launch.py .phase267-readout.py "$@" > $L-$1.stdout.log 2> $L-$1.stderr.log || { echo "failed $1" > $L.done; exit 1; }; }
r endpoints $O primary267-refined64-work/sweep-1/photons/seg-62.npz captures:primary267-refined64-work/captures
r source $O
r field $O
echo ok > $L.done
