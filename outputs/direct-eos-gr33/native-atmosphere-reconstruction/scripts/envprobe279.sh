cd /home/lpaiu/work/native-refined268-runtime || exit 1
export PYTHONPATH=/home/lpaiu/work/nutimo_pilot/request13_deps:verification
timeout 600 python3 - <<'PY'
import sys, time, numpy as np
sys.path.insert(0, 'verification')
t0 = time.monotonic()
import def_native_radiative_envelope as env
E = env.Envelope()
print('Envelope ready %.1fs' % (time.monotonic() - t0))
import gr_two_carrier_evolution as two
names = getattr(two.e.g.c, 'NAMES', None)
print('X nonzero:', [(i, float(v)) for i, v in enumerate(E.X) if v > 1e-12][:10], 'A', two.e.g.c.A[:6], 'Z', two.e.g.c.Z[:6])
print('cx', E.cx, 'arad', E.arad, 'R', E.R, 'Lref', E.Lref)
for lT, lP in [(np.log(1.4e4), np.log(1e3)), (np.log(1.7e4), np.log(2e4)), (np.log(2.5e4), np.log(1.1e5)), (np.log(4e4), np.log(6e5))]:
    a = E.eos(1, float(lP), float(lT), E.X); rho = a[0]
    kap = two.opacity_parts(E.opacity, (np.log(rho), lT, E.X))
    print('T %.0f P %.2e -> rho %.3e  a[1]/P %.9f  kappa parts %s' % (np.exp(lT), np.exp(lP), rho, a[1]/np.exp(lP), np.array2string(np.atleast_1d(np.asarray(kap, dtype=object))[:3], precision=4)))
PY
