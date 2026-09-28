# Phase268: 4x interior grid runtime and every input the final driven model uses (phase-267 procedures).
set -u
S=/mnt/c/Users/lpaiu/AppData/Local/Temp/claude/E--lab-self-mass-unobservability--claude-worktrees-eft-massive-objects-gravity-f20672/c8fbf92c-e1e4-431b-9665-aff6bfed2aa3/scratchpad
OLD=/home/lpaiu/work/native-retained-tail-runtime; NEW=/home/lpaiu/work/native-refined268-runtime
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
test ! -e $NEW || { echo "exists"; exit 1; }
mkdir -p $NEW/verification $NEW/outputs/direct-eos-gr33; LOG=$NEW/.phase268-prepare.log
step() { echo "$(date '+%F %T') $*" >> $LOG; }
fail() { step "FAILED: $*"; echo failed > $NEW/.phase268-prepare.done; exit 1; }
(cd $OLD/verification && find . -maxdepth 1 -name '*.py' -print0 | xargs -0 -I{} cp -p {} $NEW/verification/)
cp -p $OLD/.phase251-reader-launch.py $NEW/
(cd $OLD/verification && sha256sum *.py) > /tmp/o268.sha; (cd $NEW/verification && sha256sum *.py) > /tmp/n268.sha; cmp -s /tmp/o268.sha /tmp/n268.sha || fail code
step "code copied identical: $(ls $NEW/verification/*.py | wc -l) modules"
cd $NEW
for f in phase268-bank phase268-thermal phase268-initial phase268-initial-audit phase267-template phase267-undriven phase267-p150-setup phase267-p150 phase267-born phase267-placeholders phase267-driver phase268-checks; do tr -d '\r' < $S/$f.py > .$f.py; done
python3 .phase268-checks.py inputs >> $LOG 2>&1 || fail inputs
step "bank start"
/usr/bin/time -f "%e s %M KB" taskset -c 3 python3 .phase268-bank.py quad outputs/direct-eos-gr33/def-native-boundary-layer > .phase268-bank.log 2>&1 || fail "bank (see .phase268-bank.log)"
python3 .phase268-checks.py bank >> $LOG 2>&1 || fail bank-check
step "thermal start"
/usr/bin/time -f "%e s %M KB" taskset -c 3 python3 .phase268-thermal.py outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined > .phase268-thermal.log 2>&1 || fail "thermal (see .phase268-thermal.log)"
python3 .phase268-checks.py thermal >> $LOG 2>&1 || fail thermal-check
F=outputs/direct-eos-gr33/def-native-initial-constraints/finite-volume
step "initial start"
/usr/bin/time -f "%e s %M KB" taskset -c 3 python3 .phase268-initial.py $F > .phase268-initial.log 2>&1 || fail "initial (see .phase268-initial.log)"
/usr/bin/time -f "%e s %M KB" taskset -c 3 python3 .phase268-initial-audit.py $F > .phase268-initial-audit.log 2>&1 || fail "initial-audit (see .phase268-initial-audit.log)"
python3 .phase268-checks.py initial >> $LOG 2>&1 || fail initial-check
step "template and undriven start"
python3 .phase267-template.py outputs/direct-eos-gr33/def-native-anisotropic-gr/source-128.npz >> $LOG 2>&1 || fail template
T=outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/gr/source-128.npz; mkdir -p $(dirname $T); cp outputs/direct-eos-gr33/def-native-anisotropic-gr/source-128.npz $T
/usr/bin/time -f "%e s %M KB" taskset -c 3 python3 .phase267-undriven.py outputs/direct-eos-gr33/native-retained-completion 128 - > .phase268-undriven.log 2>&1 || fail "undriven (see .phase268-undriven.log)"
grep -v 'screen size' .phase268-undriven.log | tail -2 | cut -c1-600 >> $LOG
step "phase-150 bank start"
EARLY=outputs/direct-eos-gr33/def-native-conservative-rates/thermal-refined/updated-gr-return/lapse/corrected/metric-128-g8.npz
R=$(python3 .phase268-checks.py ratio) || fail ratio
python3 .phase267-p150-setup.py retained-native-return150-work $EARLY $R >> $LOG 2>&1 || fail p150-setup
RETAINED_NATIVE_OUTPUT=retained-native-return150-work /usr/bin/time -f "%e s %M KB" taskset -c 3 python3 .phase267-p150.py .phase268-p150-audit.json bank > .phase268-p150.log 2>&1 || fail "p150 bank (see .phase268-p150.log)"
step "photon bank points: $(ls retained-native-return150-work/photons/bank-128/point-*.npz | wc -l)"
step "born and placeholders start"
python3 .phase267-born.py native-incident-drive155-work > .phase268-born.log 2>&1 || fail "born (see .phase268-born.log)"
python3 .phase267-placeholders.py >> $LOG 2>&1 || fail placeholders
step "smoke start"
PHASE267_EXISTS=.phase268-exists-smoke.json /usr/bin/time -v taskset -c 3 python3 .phase267-driver.py primary268-smoke-work smoke-01 1 - - 1800 > .phase268-smoke.stdout.log 2> .phase268-smoke.stderr.log || fail "smoke (see .phase268-smoke.stderr.log)"
python3 .phase268-checks.py exists >> $LOG 2>&1 || fail exists-check
tail -1 .phase268-smoke.stdout.log | cut -c1-400 >> $LOG
step "prepared"
echo ok > $NEW/.phase268-prepare.done
