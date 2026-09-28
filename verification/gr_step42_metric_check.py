"""Counterexample candidate: isolate omitted live GR metric derivatives."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json

import numpy as np
import gr_step42_analysis as test

full, m, e, ld = test.full, test.m, test.e, test.ld
OUT = m.BASE/'step42-metric-check'
METRICS = ['m','mf','a','da','N','nur']


class LiveMetricTangent(m.original.CompatibleTangent):
    def evaluate(self, delta):
        # Reuse the established EOS/opacity tangent and projected composition.
        # A second cached evaluation retains its actual GR metric response.
        super().evaluate(delta)
        last,older,h,coefficients = self.composition_context
        projected = m.species(self,delta,last,older,h,coefficients,self.linearization)
        candidate = delta.copy()
        candidate[:,5:] = self.anchor[:,5:]+(projected-self.projected_anchor)
        return m.prior.ConservativeStar.evaluate(self,candidate)


class FrozenNativeMetric(m.original.CompatibleStar):
    def evaluate(self,delta):
        z = super().evaluate(delta)
        for key in METRICS:
            z[key] = self.fixed_metric[key]
        self.finish_moments(delta,z,z['dU'][:,1])
        return z


def run():
    assert not OUT.exists()
    plan = full.bindings()
    for rel,digest in json.loads((test.OUT/'manifest.json').read_text())['sha256'].items():
        assert e.digest(e.ROOT/rel) == digest,rel
    OUT.mkdir()
    with ProcessPoolExecutor(max_workers=15,initializer=e.worker_init) as pool:
        star = m.initialize(pool)
        times = full.time_nodes(plan,1)
        older,previous = [full.restored(star,full.OUT/'path-1',k,times[k]) for k in (40,41)]
        delta,z = test.failed_state(star)
        h = times[42]-times[41]
        coefficients = m.weights(h,times[41]-times[40])
        value,z = m.residual(star,delta,previous,older,h,coefficients)
        assert float(np.max(abs(value)/m.ATOL)) == 18.768935482897138
        cp = np.load(test.OUT/'refreshed_composition.npz')
        change = np.zeros_like(delta)
        change[:,:5] = cp['direction']
        star.composition_context = (previous,older,h,coefficients)
        star.composition_data = dict(coefficients=cp['composition_coefficients'],inverse=cp['composition_inverse'])
        tangent = m.tangent(star,delta,z)
        tangent.__class__ = LiveMetricTangent
        frozen = FrozenNativeMetric.__new__(FrozenNativeMetric)
        frozen.__dict__ = dict(star.__dict__,fixed_metric=z,material_cache=dict(star.material_cache))
        native,held,live = [],[],[]
        for sign in (1,-1):
            candidate = delta+sign*change
            candidate[:,5:] = m.species(star,candidate,previous,older,h,coefficients,z)
            native.append(m.residual(star,candidate,previous,older,h,coefficients)[0])
            held.append(m.residual(frozen,candidate,previous,older,h,coefficients)[0])
            live.append(m.residual(tangent,delta+sign*change,previous,older,h,coefficients)[0])
        actions={name:(rows[0]-rows[1])/2 for name,rows in [('native',native),('native_fixed_metric',held),('live_metric_tangent',live)]}
        old = cp['full_tangent_action']
        errors={name:(np.max(abs(v[:,:5]-actions['native'][:,:5]),axis=0)/m.ATOL[:5]).astype(float).tolist()
                for name,v in [('old_tangent',old),*actions.items()]}
        energies={name:float(v[3142,1]) for name,v in [('old_tangent',old),*actions.items()]}
        result=dict(classification='Counterexample candidate',errors_from_native=errors,limiting_energy_action=energies,
            scope='One unchanged failed native state and direction; not repaired stage or full time convergence.')
        e.write(OUT/'result.json',result)
        np.savez_compressed(OUT/'actions.npz',**actions,old_tangent=old)
        print('METRIC DERIVATIVE ISOLATION',json.dumps(result),flush=True)
    inputs=[Path(__file__),test.OUT/'manifest.json',*OUT.iterdir()]
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in inputs if p.is_file()}))


if __name__ == '__main__':
    run()
