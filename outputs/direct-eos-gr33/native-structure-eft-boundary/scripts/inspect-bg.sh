# What the declared background sampler provides (keys, radial coverage, core values) - read-only.
cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'EOF'
import inspect, sys
import numpy as np
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
bg = bl.chem.prior.Background()
print('Background class', type(bg).__module__, [k for k in dir(bg) if not k.startswith('_')][:40])
print('R', getattr(bg, 'R', None))
for attr in ['d', 'env']:
    v = getattr(bg, attr, None)
    if isinstance(v, dict): print(attr, {k: np.shape(x) for k, x in v.items()})
r = np.array([1e7, 1e8, 1e9, 3e9, 5e9, 6.5e9, 6.85e9, 6.9e9])
try:
    s = bg.sample(r); print('sample keys', list(s.keys()))
    for k, v in s.items(): print(k, np.array2string(np.asarray(v), precision=4))
except Exception as exc: print('sample failed', repr(exc)[:300])
print(inspect.getsource(type(bg).sample)[:2500])
EOF
