from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np,json
import gr_flux_coupled_evolution as f
import gr_conservative_evolution as p
import gr_coupled_evolution as e
d=np.load(p.OUT/'initial.npz');rows=[]
with ProcessPoolExecutor(max_workers=1,initializer=e.worker_init) as pool:
 for cls in [e.Star,f.FluxStar]:
  star=cls(5735,pool);star.base=d['base'].copy();star.qscale=d['qscale'].copy()
  star.material_cache={e.material_key(row):aux for row,aux in zip(zip(star.base[:,0],star.base[:,1],star.base[:,4:]),d['aux'])}
  rate,z=star.rhs(np.zeros_like(star.base))
  assert np.all(star.base[:,2]==0) and np.all(rate[:,4:]==0)
  er=z['rho']*(z['rest']+z['u']+z['aux'][:,9]);et=z['rho']*z['aux'][:,10]
  actual=er*rate[:,0]+et*rate[:,1]+2*z['Q']*rate[:,2]
  expected=-e.C*star.divergence(star.faces(z['N']*z['S']/z['a'],odd=True))
  error=(actual-expected)/et
  row=dict(implementation=cls.__name__,maximum_initial_energy_rate_defect_over_heat_capacity=float(np.max(abs(error))),relative_rate_score=float(np.max(abs(actual-expected)/(et+abs(expected)))),worst_original_cell=int(star.indices[np.argmax(abs(error))]))
  rows.append(row)
assert rows[0]['relative_rate_score']>1e-5
assert rows[1]['relative_rate_score']<1e-10
out=dict(classification='Counterexample candidate',passed=True,rows=rows,initial_v_zero_only=True,bindings={str(q.relative_to(e.ROOT)):e.digest(q) for q in [p.OUT/'initial-manifest.json',f.OUT/'initial-manifest.json',Path(f.__file__),Path(e.__file__)]})
path=f.OUT/'initial-energy-identity.json'
if path.exists(): assert json.loads(path.read_text())==out
else: e.write(path,out)
print(json.dumps(out))
