"""A stored failing evolved state must lose the baryon chain defect."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
import numpy as np
import gr_coupled_evolution as e
import gr_flux_coupled_evolution as heat
import gr_baryon_flux_evolution as fixed

initial = np.load(heat.OUT/'initial.npz')
end = np.load(heat.OUT/'path-8/step-0008.npz')
aux = np.load(heat.OUT/'budget-8.npz')['aux']
y = initial['base']+end['delta']
records = []
with ProcessPoolExecutor(max_workers=1, initializer=e.worker_init) as pool:
    for cls in [heat.FluxStar, fixed.BaryonStar]:
        star = cls(5735, pool)
        star.base, star.qscale = initial['base'].copy(), initial['qscale'].copy()
        star.material_cache = {e.material_key(row): value for row, value in
            zip(zip(y[:, 0], y[:, 1], y[:, 4:]), aux)}
        rate, z = star.rhs(end['delta'])
        v,W,rho,P,Q,w,a,N = [z[k] for k in ['v','W','rho','P','Q','w','a','N']]
        er = rho*(z['rest']+z['u']+aux[:, 9])
        et = rho*aux[:, 10]
        pr,pt = P*aux[:, 5],P*aux[:, 6]
        E_t = (W*W*(er+pr*v*v)*rate[:, 0]+W*W*(et+pt*v*v)*rate[:, 1]
            +2*W**4*(v*w+Q*(1+v*v))*rate[:, 2]+2*v*W*W*star.qscale*rate[:, 3])
        expected_E = -e.C*star.divergence(star.faces(N*z['S']/a, odd=True))
        B_t = a*(z['D']*rate[:, 0]+rho*W**3*v*rate[:, 2])+z['at']*z['D']
        expected_B = -e.C*star.divergence(star.faces(N*z['D']*v, odd=True))
        # No omitted chemical time derivative in this selected surface row.
        assert np.all(rate[-1, 4:] == 0)
        records.append(dict(implementation=cls.__name__,
            maximum_relative_baryon_rate_defect=float(np.max(abs((B_t-expected_B)/(a*rho)))),
            outer_energy_rate_defect_over_heat_capacity=float((E_t[-1]-expected_E[-1])/et[-1]),
            outer_composition_time_derivative_zero=True))
assert records[0]['maximum_relative_baryon_rate_defect'] > 1e-8
assert records[1]['maximum_relative_baryon_rate_defect'] < 1e-15
assert abs(records[0]['outer_energy_rate_defect_over_heat_capacity']) > 1
assert abs(records[1]['outer_energy_rate_defect_over_heat_capacity']) < 1e-4
inputs = [Path(__file__), Path(fixed.__file__), Path(heat.__file__),
          heat.OUT/'initial-manifest.json', heat.OUT/'path-8/manifest.json', heat.OUT/'budget-8-manifest.json']
result = dict(classification='Counterexample candidate', passed=True, rows=records,
    scope='One actual evolved state; baryon rates at all nodes and energy rate at its zero-composition-rate surface. No whole-path or continuous energy certificate.',
    bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in inputs})
path = fixed.OUT/'baryon-chain-check.json'
if path.exists(): assert json.loads(path.read_text()) == result
else: e.write(path, result)
print(json.dumps(result), flush=True)
