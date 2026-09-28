# Phase267 refined 64-clock primary: restartable segments of 4 macro steps (skip the begin=60 boundary).
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-refined267-runtime
[ -f .phase267-driver-run.py ] || tr -d '\r' < $S/phase267-driver.py > .phase267-driver-run.py
W=primary267-refined64-work; mkdir -p .phase267-refined64-logs
prev=-
for lim in 4 8 12 16 20 24 28 32 36 40 44 48 52 56 61 62 63 64; do
  label=$(printf 'seg-%02d' $lim)
  if [ -f $W/$label-driver.json ] && grep -q '"passed": true' $W/$label-driver.json; then prev=$label; continue; fi
  echo "$(date '+%F %T') start $label restart=$prev" >> .phase267-refined64-logs/status.txt
  EXC=; NLEXC=; CAP=; case $label in seg-62|seg-63|seg-64) EXC=inf; NLEXC=inf; CAP=60;; esac  # user approval 3 (2026-09-27): last segment, physical/material gates, vector recorded only
  PHASE267_VECTOR_EXCEPTION=$EXC PHASE267_NONLINEAR_EXCEPTION=$NLEXC PHASE267_POLISH_CAP=$CAP PHASE267_EXISTS=.phase267-refined64-logs/exists-$label.json /usr/bin/time -v taskset -c 5 python3 .phase267-driver-run.py $W $label $lim $prev - 10800 \
    > .phase267-refined64-logs/$label.stdout.log 2> .phase267-refined64-logs/$label.stderr.log
  code=$?
  echo "$(date '+%F %T') end $label exit=$code" >> .phase267-refined64-logs/status.txt
  if [ $code -ne 0 ]; then echo "stopped at $label exit=$code" > .phase267-refined64-logs/stopped.txt; exit $code; fi
  prev=$label
done
echo "$(date '+%F %T') complete" > .phase267-refined64-logs/complete.txt
