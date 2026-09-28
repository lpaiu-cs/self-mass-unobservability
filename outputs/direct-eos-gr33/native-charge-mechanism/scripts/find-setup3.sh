# The original source-coefficient assembly (state readout + pulse x geometry) used by the readout (read-only).
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'EOF'
import inspect, sys
sys.path.insert(0, 'verification')
import read_complete_radau_history as rch
p = rch.prior
print('=== prior.coefficients in', inspect.getsourcefile(p.coefficients))
print(inspect.getsource(p.coefficients)[:5000])
EOF
