"""Counterexample candidate: isolate a donor switch from native EOS response."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json

import numpy as np
import gr_step42_analysis as test

full,m,e,ld = test.full,test.m,test.e,test.ld
OUT = m.BASE/'step42-upwind-check'


class FixedDonorNative(m.original.CompatibleStar):
    def finish_moments(self,delta,z,dE):
        super().finish_moments(delta,z,dE)
        m.original.finish_fluxes(self,z,self.fixed_donor)


def run():
    assert not OUT.exists()
    plan=full.bindings()
    OUT.mkdir()
    with ProcessPoolExecutor(max_workers=15,initializer=e.worker_init) as pool:
        star=m.initialize(pool)
        times=full.time_nodes(plan,1)
        older,previous=[full.restored(star,full.OUT/'path-1',k,times[k]) for k in (40,41)]
        delta,z=test.failed_state(star)
        h=times[42]-times[41]
        coefficients=m.weights(h,times[41]-times[40])
        value,z=m.residual(star,delta,previous,older,h,coefficients)
        cp=np.load(test.OUT/'refreshed_composition.npz')
        change=np.zeros_like(delta);change[:,:5]=cp['direction']
        star.composition_context=(previous,older,h,coefficients)
        star.composition_data=dict(coefficients=cp['composition_coefficients'],inverse=cp['composition_inverse'])
        tangent=m.tangent(star,delta,z)
        anchor_donor=star.faces(z['N']*z['v']/z['a'],odd=True)
        frozen=FixedDonorNative.__new__(FixedDonorNative)
        frozen.__dict__=dict(star.__dict__,fixed_donor=anchor_donor,material_cache=dict(star.material_cache))
        rows=[]
        for fraction in (ld(1),ld(1)/128,ld(1)/16384):
            candidate=delta+fraction*change
            candidate[:,5:]=m.species(star,candidate,previous,older,h,coefficients,z)
            native,native_z=m.residual(star,candidate,previous,older,h,coefficients)
            held,held_z=m.residual(frozen,candidate,previous,older,h,coefficients)
            surrogate,surrogate_z=m.residual(tangent,delta+fraction*change,previous,older,h,coefficients)
            donor=star.faces(native_z['N']*native_z['v']/native_z['a'],odd=True)
            changed=np.flatnonzero((donor>=0)!=(anchor_donor>=0))
            item=dict(fraction=float(fraction),native_score=float(np.max(abs(native)/m.ATOL)),
                fixed_donor_native_score=float(np.max(abs(held)/m.ATOL)),
                tangent_score=float(np.max(abs(surrogate)/m.ATOL)),switched_faces=changed.tolist(),
                limiting_cell_energy=[float(a[3142,1]) for a in (value,native,held,surrogate)],
                native_vs_tangent_energy_parts=dict(
                    stored_energy=float(coefficients[0]*(native_z['dU'][3142,1]-surrogate_z['dU'][3142,1])/star.heat0[3142]),
                    transport=float(h*e.C*star.divergence(native_z['fluxes'][1]-surrogate_z['fluxes'][1])[3142]/star.heat0[3142])),
                native_aux_energy=float(native_z['aux'][3142,2]),tangent_aux_energy=float(surrogate_z['aux'][3142,2]),
                cell_direction=change[3142,:5].astype(float).tolist(),
                adjacent_anchor_donor=anchor_donor[3142:3144].astype(float).tolist(),
                adjacent_trial_donor=donor[3142:3144].astype(float).tolist())
            rows.append(item)
            print('DONOR ISOLATION',json.dumps(item),flush=True)
        result=dict(classification='Counterexample candidate',rows=rows,
            initial_score=float(np.max(abs(value)/m.ATOL)),scope='Same failed state, native EOS and BDF step. Only donor selection differs in the native counterfactual.')
        e.write(OUT/'result.json',result)
    files=[Path(__file__),test.OUT/'manifest.json',*OUT.iterdir()]
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files if p.is_file()}))


if __name__=='__main__':
    run()
