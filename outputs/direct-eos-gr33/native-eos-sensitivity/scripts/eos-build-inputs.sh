# Inputs for the phase-269 library rebuilds (read-only).
ls -la /home/lpaiu/work/direct-eos-gr33/native-cold-population/ | cut -c30-
ls /home/lpaiu/work/direct-eos-gr33/native-cold-population/stable/
ls -la /mnt/e/lab/self-mass-unobservability/outputs/direct-eos-gr33/gr-radiation-eos-split/gas-bridge.f90 /mnt/e/lab/self-mass-unobservability/outputs/direct-eos-gr33/def-native-cold-population/diagnostic-bridge.f90 2>&1 | cut -c30-
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'EOF'
import json, sys
sys.path.insert(0, 'verification')
import def_native_cold_population as cp, def_photon_eos_populations as pp
r = json.load(open(cp.OUT/'stable-build.json'))
link = r['receipts'][1]['command']; print('LINK objects', len(link)); print([w for w in link if 'mod_free_eos' in w or 'mod_eos_calc' in w or 'excitation' in w])
print('levels SOURCE', pp.SOURCE)
b = json.load(open(pp.OUT/'levels-repaired/build.json')); print('LEVELS CMD', b['command'])
EOF
