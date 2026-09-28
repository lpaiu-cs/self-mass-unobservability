# Print the native EOS lineage paths and the exact stable-build commands (read-only).
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'EOF'
import json, sys
sys.path.insert(0, 'verification')
import def_native_ion_closure as ic, def_native_cold_population as cp, def_native_hydrogen_exchange as hx
print('ion closure SOURCE', ic.SOURCE); print('ion closure BUILD', ic.BUILD); print('ion closure CACHE', ic.CACHE)
print('cold CACHE', cp.CACHE, 'OUT', cp.OUT); print('LEVELS', hx.LEVELS, 'ATOMIC', hx.ATOMIC)
r = json.load(open(cp.OUT/'stable-build.json'))
for x in r['receipts']: print('CMD', ' '.join(x['command'])[:900])
EOF
grep -n "ifmodified.eq.11" -A4 $(python3 -c "import sys; sys.path.insert(0,'verification'); import def_native_ion_closure as ic; print(ic.SOURCE)")/mod_free_eos.f90 | head -8
grep -n "calculate planck-larkin occupation probabilities" $(python3 -c "import sys; sys.path.insert(0,'verification'); import def_native_ion_closure as ic; print(ic.SOURCE)")/*.f90 | head
