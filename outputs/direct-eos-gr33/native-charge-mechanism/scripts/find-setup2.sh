# How the 244 source components become the scalar source (read-only): bridge.coefficients and base gr Response.setup.
cd /home/lpaiu/work/native-retained-tail-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'EOF'
import inspect, sys
sys.path.insert(0, 'verification')
import apply_actual_stage_gr as a
bridge = a.bridge
print('=== bridge.coefficients in', inspect.getsourcefile(bridge.coefficients))
print(inspect.getsource(bridge.coefficients)[:3500])
setup = bridge.base.gr.Response.setup
print('=== base gr Response.setup in', inspect.getsourcefile(setup), 'line', inspect.getsourcelines(setup)[1])
print(inspect.getsource(setup)[:7000])
EOF
