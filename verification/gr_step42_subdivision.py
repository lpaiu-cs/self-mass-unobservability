"""Counterexample candidate: solve the failed interval as two native BDF steps.

Keep the last two accepted states and accumulated flux histories. This changes
only the time grid; no repaired-donor counterfactual becomes an accepted state.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
import time

import numpy as np
import gr_compatible_full_duration as full

m,e,ld=full.cached,full.e,full.ld
OUT=m.BASE/'step42-subdivision'


def run():
    assert not OUT.exists()
    plan=full.bindings()
    times=full.time_nodes(plan,1)
    nodes=np.linspace(times[41],times[42],3,dtype=ld)
    files=[Path(__file__),full.OUT/'plan.json',full.OUT/'path-1/failure.json',
        full.OUT/'path-1/step-0040.npz',full.OUT/'path-1/step-0041.npz']
    OUT.mkdir()
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in files},
        times_seconds=list(map(str,nodes)),rule='Two native variable-step BDF2 half-steps from the same accepted 40/41 history. Original native solver, equations, 24 iterations, eight backtracks, tolerances and gates. Preserve failures.',
        full_duration=False))
    began=time.perf_counter()
    rows=[]
    with ProcessPoolExecutor(max_workers=15,initializer=e.worker_init) as pool:
        star=m.initialize(pool)
        older,previous=[full.restored(star,full.OUT/'path-1',k,times[k]) for k in (40,41)]
        star.operator_step,star.operator_time=41,times[41]
        cp40,cp41=[np.load(full.OUT/'path-1'/f'step-{k:04d}.npz') for k in (40,41)]
        energy,old_energy=cp41['integrated_energy_flux'].copy(),cp40['integrated_energy_flux'].copy()
        baryon,old_baryon=cp41['integrated_baryon_flux'].copy(),cp40['integrated_baryon_flux'].copy()
        old_h=times[41]-times[40]
        for step,(a,b) in enumerate(zip(nodes[:-1],nodes[1:]),42):
            h=b-a;coefficients=m.weights(h,old_h)
            def record(row):
                with (OUT/'iterations.jsonl').open('a') as stream:
                    stream.write(json.dumps(dict(step=step,**row))+'\n')
                print('NATIVE HALF STEP',step,row['iteration'],row['residual_norm'],flush=True)
            try:
                delta,z=m.stage(star,previous,older,h,coefficients,record)
            except Exception as error:
                if hasattr(star,'last_delta'):
                    np.savez_compressed(OUT/'failed-iterate.npz',delta=star.last_delta,aux=star.last_state['aux'],time_seconds=b)
                e.write(OUT/'failure.json',dict(classification='Counterexample candidate',step=step,reason=repr(error)))
                raise
            c0,c1,c2=coefficients;area=4*np.pi*star.rf**2
            energy,old_energy=(-c1*energy-c2*old_energy+h*e.C*area*z['fluxes'][1])/c0,energy
            baryon,old_baryon=(-c1*baryon-c2*old_baryon+h*e.C*area*z['fluxes'][0])/c0,baryon
            budget,energy_defect,baryon_defect=m.budget(star,z,energy,baryon)
            cone=m.parent.cones(z)
            native,_=m.residual(star,delta,previous,older,h,coefficients)
            assert np.max(abs(native)/m.ATOL)<=1 and m.budget_passed(budget,plan) and cone['sampled_cone_inside_light_cone']
            m.save_state(OUT/f'step-{step:04d}.npz',delta,z,b,energy,baryon,energy_defect,baryon_defect)
            rows.append(dict(step=step,time_seconds=str(b),native_score=float(np.max(abs(native)/m.ATOL)),budget=budget,cone=cone,elapsed_seconds=time.perf_counter()-began))
            e.write(OUT/'progress.json',dict(classification='Counterexample candidate',rows=rows))
            print('ACCEPTED NATIVE HALF STEP',json.dumps(rows[-1]),flush=True)
            older,previous=previous,(delta,z)
            old_h=h
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,
        scope='Original failed interval completed with two half-steps; full duration and three-grid convergence remain open.')
    e.write(OUT/'result.json',result)
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in [Path(__file__),*OUT.iterdir()] if p.is_file()}))


if __name__=='__main__':
    run()
