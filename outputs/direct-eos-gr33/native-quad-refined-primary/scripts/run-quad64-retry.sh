# Phase268: preserve attempt 1 of the 4x seg-64 and rerun it with driver v6 (long-double corrections 12 -> 24 under the
# 2026-09-24 resource policy; every acceptance gate and the registered last-segment rule unchanged), then read T with depth bands.
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-refined268-runtime || exit 1
W=primary268-quad64-work; L=.phase268-quad64-logs
test -d $W/attempt1-seg-64 && { echo "attempt1 already preserved"; exit 1; }
mkdir -p $W/attempt1-seg-64/captures
for f in seg-64-fallbacks.json seg-64-driver.json seg-64-nonlinear-exceptions.json integer-64.json polish.json rejected-joint-stage.npz rejected-joint-stage.json; do [ -f $W/$f ] && mv $W/$f $W/attempt1-seg-64/; done
for f in $W/sweep-1/photons/seg-64*; do [ -f "$f" ] && mv "$f" $W/attempt1-seg-64/; done
for i in $(seq 252 299); do f=$(printf '%s/captures/captured-64-%03d.npz' $W $i); [ -f $f ] && mv $f $W/attempt1-seg-64/captures/; done
for f in seg-64.stdout.log seg-64.stderr.log exists-seg-64.json stopped.txt; do [ -f $L/$f ] && mv $L/$f $L/attempt1-$f; done
echo "$(date '+%F %T') attempt1 of seg-64 stopped at stage 252: joint solve diverged (relative 1.6e3); 12 long-double corrections reached physical moments 7.5e-15 but material components 1.1e-12 > 1e-13 (falling about tenfold per correction); production polish 3.6e-6; preserved in $W/attempt1-seg-64" >> $L/status.txt
tr -d '\r' < $S/phase268-driver.py > .phase268-driver.py
echo "$(date '+%F %T') driver v6 installed (resource policy 2026-09-24: long-double corrections 12 -> 24 for seg-64; every gate unchanged)" >> $L/status.txt
echo "$(date '+%F %T') start seg-64 restart=seg-63" >> $L/status.txt
PHASE268_LD_CORRECTIONS=24 PHASE267_VECTOR_EXCEPTION=inf PHASE267_NONLINEAR_EXCEPTION=inf PHASE267_POLISH_CAP=60 PHASE267_EXISTS=$L/exists-seg-64.json \
  /usr/bin/time -v taskset -c 5 python3 .phase268-driver.py $W seg-64 64 seg-63 - 10800 > $L/seg-64.stdout.log 2> $L/seg-64.stderr.log
code=$?
echo "$(date '+%F %T') end seg-64 exit=$code" >> $L/status.txt
if [ $code -ne 0 ]; then echo "stopped at seg-64 exit=$code" > $L/stopped.txt; exit $code; fi
echo "$(date '+%F %T') complete" > $L/complete.txt
# T readout with depth bands (same procedure as run-readouts268.sh; that watcher exited when attempt 1 stopped)
R2=/home/lpaiu/work/native-refined267-runtime; N=64; C=13; O=readout268-quad$N-work; Q=.phase268-readout-quad$N
echo "$(date '+%F %T') readout $N start (core $C)" >> .phase268-readouts.log
python3 -c "import numpy as np; a=float(np.load('$W/sweep-1/photons/seg-$N.npz')['actual_step_edges'][-1]); b=float(np.load('$R2/primary267-refined64-work/sweep-1/photons/seg-$N.npz')['actual_step_edges'][-1]); print(repr(a), repr(b), a == b)" > $Q-time.txt
r() { taskset -c $C python3 .phase267-audited.py $Q-$1-record.json .phase251-reader-launch.py .phase267-readout.py "$@" > $Q-$1.stdout.log 2> $Q-$1.stderr.log || { echo "failed $1" > $Q.done; echo "$(date '+%F %T') readout $N failed $1" >> .phase268-readouts.log; exit 1; }; }
r endpoints $O $W/sweep-1/photons/seg-$N.npz captures:$W/captures
r source $O
r field $O
bands="cells:0-7"; for c in $(seq 8 42); do bands="$bands cells:$c-$c"; done; bands="$bands cells:43-170 cells:171-298 cells:299-426 cells:427-554 boundary all"
taskset -c $C python3 .phase251-reader-launch.py .phase267-depth.py $O bands $bands > $Q-depth.stdout.log 2> $Q-depth.stderr.log || { echo "failed depth" > $Q.done; echo "$(date '+%F %T') readout $N failed depth" >> .phase268-readouts.log; exit 1; }
echo ok > $Q.done; echo "$(date '+%F %T') readout $N end ok" >> .phase268-readouts.log
