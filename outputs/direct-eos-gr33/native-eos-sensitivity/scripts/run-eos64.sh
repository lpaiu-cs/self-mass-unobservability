# Phase269 stage 2: 2x primary with the PL-off native EOS and level libraries, restartable segments (phase-268 plan).
# Macro 0-61 under every original gate; segments 62-64 with the registered last-segment rule (physical moments and
# material components 1e-13, vectors recorded, polish cap 60 s). Long-double corrections 24 throughout (resource policy).
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
cd /home/lpaiu/work/native-eos269-runtime || exit 1
W=primary269-eos64-work; L=.phase269-eos64-logs; mkdir -p $L
prev=-
for lim in 4 8 12 16 20 24 28 32 36 40 44 48 52 56 61 62 63 64; do
  label=$(printf 'seg-%02d' $lim)
  if [ -f $W/$label-driver.json ] && grep -q '"passed": true' $W/$label-driver.json; then prev=$label; continue; fi
  echo "$(date '+%F %T') start $label restart=$prev" >> $L/status.txt
  EXC=; NLEXC=; CAP=; case $label in seg-62|seg-63|seg-64) EXC=inf; NLEXC=inf; CAP=60;; esac
  PHASE268_LD_CORRECTIONS=24 PHASE267_VECTOR_EXCEPTION=$EXC PHASE267_NONLINEAR_EXCEPTION=$NLEXC PHASE267_POLISH_CAP=$CAP PHASE267_EXISTS=$L/exists-$label.json \
    /usr/bin/time -v taskset -c 5 python3 .phase268-driver.py $W $label $lim $prev - 10800 > $L/$label.stdout.log 2> $L/$label.stderr.log
  code=$?
  echo "$(date '+%F %T') end $label exit=$code" >> $L/status.txt
  if [ $code -ne 0 ]; then echo "stopped at $label exit=$code" > $L/stopped.txt; exit $code; fi
  prev=$label
done
echo "$(date '+%F %T') complete" > $L/complete.txt
