# Phase268 preparation, resumed at the Born step (Born requires a fresh phase-155 folder; the input copy had already
# placed the grid-independent plan JSON there, which phase 267 placed by hand after Born). The JSON is moved aside
# during Born and restored bitwise afterwards.
set -u
NEW=/home/lpaiu/work/native-refined268-runtime
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
LOG=$NEW/.phase268-prepare.log
step() { echo "$(date '+%F %T') $*" >> $LOG; }
fail() { step "FAILED: $*"; echo failed > $NEW/.phase268-prepare.done; exit 1; }
cd $NEW; rm -f .phase268-prepare.done; mv .phase268-born.log .phase268-born-attempt1.log
P=native-incident-drive155-work/material-precision-plan.json; H=$(sha256sum $P | cut -d' ' -f1)
step "born and placeholders start (plan JSON set aside, sha256 $H)"
mv $P .phase268-plan-aside.json && rmdir native-incident-drive155-work || fail set-aside
python3 .phase267-born.py native-incident-drive155-work > .phase268-born.log 2>&1 || fail "born (see .phase268-born.log)"
mv .phase268-plan-aside.json $P || fail restore
test "$(sha256sum $P | cut -d' ' -f1)" = "$H" || fail restore-hash
step "born done; plan JSON restored bitwise"
python3 .phase267-placeholders.py >> $LOG 2>&1 || fail placeholders
step "smoke start"
PHASE267_EXISTS=.phase268-exists-smoke.json /usr/bin/time -v taskset -c 3 python3 .phase267-driver.py primary268-smoke-work smoke-01 1 - - 1800 > .phase268-smoke.stdout.log 2> .phase268-smoke.stderr.log || fail "smoke (see .phase268-smoke.stderr.log)"
python3 .phase268-checks.py exists >> $LOG 2>&1 || fail exists-check
tail -1 .phase268-smoke.stdout.log | cut -c1-400 >> $LOG
step "prepared"
echo ok > $NEW/.phase268-prepare.done
