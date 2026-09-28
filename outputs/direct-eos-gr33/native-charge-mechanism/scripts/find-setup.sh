# Locate the retarded-field setup used by the readout and how it consumes the 244 source arrays (read-only).
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'EOF'
import inspect, sys
sys.path.insert(0, 'verification')
import extend_retarded_history as ext, read_full_captured_history as cap
rch = cap.prior
print('rch module', rch.__name__, 'KEYS', rch.KEYS)
setup = ext.previous.prior.Response.setup
print('setup defined in', inspect.getsourcefile(setup), 'line', inspect.getsourcelines(setup)[1])
src = inspect.getsource(setup)
print(src[:6000])
EOF
