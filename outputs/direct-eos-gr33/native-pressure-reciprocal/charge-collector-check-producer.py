import json,resource,signal,time
from pathlib import Path
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));signal.alarm(15)
started=time.monotonic()
import read_native_pressure_charge as r
source=r.active_geometry_source(r.compact_owner.run.base.old.aligned.base.source)
compile(source,'active-source','exec')
fake=r.SimpleNamespace(t=r.np.array([0.,1.,2.]),metric={'delta_u':r.np.full((3,2),99.)},
    fields=lambda t:(None,None,r.np.array([[t,t+1],[0.,0.],[2*t,3*t]])))
ns=dict(np=r.np,m=fake,AMP=.25);exec(r.GEOMETRY.replace('\n        ','\n'),ns)
assert r.np.array_equal(ns['phi'],.25*r.np.array([[0,1],[1,2],[2,3]]))
assert r.np.array_equal(ns['lam'],.25*r.np.array([[0,0],[2,3],[4,6]]))
assert "phi=m.metric" not in source and 'field=fields[k]' in source
value=dict(classification='Proven',passed=True,seconds=time.monotonic()-started,
    scope='Generated collector compiles and queries each actual clock value instead of a stale constructor cache. Synthetic nonzero time-dependent regression; no accepted175 GR execution or physical charge result.')
Path('native-pressure-reciprocal175-work/charge-collector-code-check.json').write_text(json.dumps(value,indent=2))
print(json.dumps(value))
