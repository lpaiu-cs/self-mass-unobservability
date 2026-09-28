cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
python3 - <<'PY'
import sys
import numpy as np
sys.path.insert(0, 'verification')
import def_native_boundary_layer as bl
bg = bl.chem.prior.Background(); d = bg.d
for k, v in d.items():
    a = np.asarray(v)
    if a.dtype.kind in 'fi' and a.size > 1: print(k, a.shape, '%.4e .. %.4e' % (a.min(), a.max()))
    else: print(k, a.shape, str(a)[:80])
PY
