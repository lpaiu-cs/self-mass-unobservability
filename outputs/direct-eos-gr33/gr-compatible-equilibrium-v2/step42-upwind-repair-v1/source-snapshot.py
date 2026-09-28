"""Counterexample candidate: pivot the iteration donor after a failed line search.

The native upwind flux and all 31 acceptance tolerances remain unchanged.
Only the sparse iteration matrix can select the donor predicted by its trial.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json

import numpy as np
import gr_step42_analysis as test

full,m,e,ld = test.full,test.m,test.e,test.ld
OUT = m.BASE/'step42-upwind-repair'


class PredictedDonorTangent(m.original.CompatibleTangent):
    def finish_moments(self,delta,z,dE):
        super().finish_moments(delta,z,dE)
        m.original.finish_fluxes(self,z,self.selected_donor)


def correction(star,delta,z,previous,older,h,coefficients):
    value,z=m.residual(star,delta,previous,older,h,coefficients)
    norm=float(np.max(abs(value)/m.ATOL))
    donor=star.faces(z['N']*z['v']/z['a'],odd=True)
    trials=[]
    for pivot in range(8):
        model=m.tangent(star,delta,z)
        model.__class__=PredictedDonorTangent
        model.selected_donor=donor
        matrix=m.jacobian(model,delta,previous,older,h,coefficients)
        direction=-m.splu(matrix).solve(np.asarray(value[:,:5]/m.SCALE,float).ravel())
        change=direction.reshape(star.n,5).astype(ld)*m.SCALE
        fraction=min(1.,.1/max(float(abs(change[:,:2]).max()),1e-300),.01/max(float(abs(change[:,2]).max()),1e-300))
        candidate=delta.copy()
        candidate[:,:5]+=fraction*change
        candidate[:,5:]=m.species(star,candidate,previous,older,h,coefficients,z)
        native,state=m.residual(star,candidate,previous,older,h,coefficients)
        score=float(np.max(abs(native)/m.ATOL))
        predicted=star.faces(state['N']*state['v']/state['a'],odd=True)
        switched=int(np.count_nonzero((donor>=0)!=(predicted>=0)))
        trials.append(dict(pivot=pivot,native_score=score,changed_donors=switched,fraction=fraction))
        print('NATIVE DONOR PIVOT',json.dumps(trials[-1]),flush=True)
        if score<=1 or score<norm*(1-1e-4*fraction):
            return candidate,native,state,trials
        if not switched:
            break
        donor=predicted
    raise RuntimeError(('Predicted donor correction failed',norm,trials))


def finish(star,delta,z,previous,older,h,coefficients,log,next_iteration):
    for iteration in range(next_iteration,24):
        delta,value,z,trials=correction(star,delta,z,previous,older,h,coefficients)
        star.last_delta,star.last_state=delta.copy(),z
        norm=float(np.max(abs(value)/m.ATOL))
        log(dict(iteration=iteration,residual_norm=norm,
            maximum_absolute_residual=np.max(abs(value),axis=0).astype(float).tolist(),donor_pivots=trials))
        if norm<=1:
            return delta,z
    raise RuntimeError(('Native donor repair exhausted original 24-iteration limit',norm))


def stage(star,previous,older,h,coefficients,log):
    records=[]
    def record(row):
        records.append(row)
        log(row)
    try:
        return m.stage(star,previous,older,h,coefficients,record)
    except RuntimeError as error:
        if not error.args or not isinstance(error.args[0],tuple) or error.args[0][0]!='Conservative native line search failed':
            raise
        return finish(star,star.last_delta,star.last_state,previous,older,h,coefficients,record,len(records))


def check():
    assert not OUT.exists()
    plan=full.bindings()
    assert m.prior.symbolic()['passed']
    OUT.mkdir()
    with ProcessPoolExecutor(max_workers=15,initializer=e.worker_init) as pool:
        star=m.initialize(pool)
        times=full.time_nodes(plan,1)
        older,previous=[full.restored(star,full.OUT/'path-1',k,times[k]) for k in (40,41)]
        delta,z=test.failed_state(star)
        h=times[42]-times[41]
        coefficients=m.weights(h,times[41]-times[40])
        value,z=m.residual(star,delta,previous,older,h,coefficients)
        initial=float(np.max(abs(value)/m.ATOL))
        assert initial==18.768935482897138
        with np.load(test.OUT/'lagged_composition.npz') as cp:
            star.composition_data=dict(coefficients=cp['composition_coefficients'],inverse=cp['composition_inverse'])
        star.composition_context=(previous,older,h,coefficients)
        records=[]
        delta,z=finish(star,delta,z,previous,older,h,coefficients,records.append,14)
        value,z=m.residual(star,delta,previous,older,h,coefficients)
        assert float(np.max(abs(value)/m.ATOL))<=1
        with np.load(full.OUT/'path-1/step-0041.npz') as last, np.load(full.OUT/'path-1/step-0040.npz') as before:
            c0,c1,c2=coefficients
            area=4*np.pi*star.rf**2
            energy=(-c1*last['integrated_energy_flux']-c2*before['integrated_energy_flux']+h*e.C*area*z['fluxes'][1])/c0
            baryon=(-c1*last['integrated_baryon_flux']-c2*before['integrated_baryon_flux']+h*e.C*area*z['fluxes'][0])/c0
        budget,energy_defect,baryon_defect=m.budget(star,z,energy,baryon)
        cone=m.parent.cones(z)
        assert m.budget_passed(budget,plan) and cone['sampled_cone_inside_light_cone']
        m.save_state(OUT/'step-0042.npz',delta,z,times[42],energy,baryon,energy_defect,baryon_defect)
        result=dict(classification='Counterexample candidate',passed=True,initial_score=initial,
            final_score=float(np.max(abs(value)/m.ATOL)),iterations=records,budget=budget,cone=cone,
            symbolic_passed=True,same_native_equations=True,same_time_step=True,same_tolerances=True,
            full_duration_completed=False,scope='One repaired native step from the exact preserved failed iterate; whole duration and time convergence remain open.')
        e.write(OUT/'result.json',result)
        print('REPAIRED NATIVE STEP',json.dumps(result),flush=True)
    files=[Path(__file__),test.OUT/'manifest.json',*OUT.iterdir()]
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files if p.is_file()}))


if __name__=='__main__':
    check()
