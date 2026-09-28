"""Counterexample candidate: physical matter work, not scalar scheme defects.

The predecessor is frozen and failed. Reuse its undriven path; attempt only
the driven native step. A scalar energy gate failure remains a failure.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import time
import numpy as np
import def_spherical_coupled as prior

e,m,ld=prior.e,prior.m,prior.ld
OUT=e.g.OUT/'def-spherical-physical-exchange'


def vacuum_check():
    s=SimpleNamespace(n=4,r=np.array([1,3,5,7],dtype=ld)*ld('1e6'),
        rf=np.array([0,2,4,6,8],dtype=ld)*ld('1e6'),h=ld('1e-5'),amplitude=ld('1e-6'))
    s.volume=4*np.pi/3*np.diff(s.rf**3);s.area=4*np.pi*s.rf**2
    s.partial_volume=4*np.pi/3*(s.r**3-s.rf[:-1]**3)
    s.distance=np.r_[ld(1),np.diff(s.r),s.rf[-1]-s.r[-1]];s.face_weight=s.area*s.distance
    zero=np.zeros(s.n,dtype=ld);zf=np.zeros(s.n+1,dtype=ld)
    s.b0=np.ones(s.n,dtype=ld);s.bf0=np.ones(s.n+1,dtype=ld)
    s.reference=dict(E=zero,m=zero,mf=zf,a=s.b0)
    s.divergence=lambda f:np.diff(s.area*f)/s.volume
    z=prior.metric(s,zero,zero,zero,zf);s.scalar_reference=dict(z,trace=zero)
    for _ in range(12):
        psi,Pi,Phi=prior.wave(s,z,zero)
        z=prior.metric(s,zero,zero,Pi,Phi)
    work=prior.scalar_work(s,z,zero,psi,Pi,Phi)
    assert work['scalar_residual']<2e-17
    assert np.max(abs(work['geometric_work']))>0
    return dict(classification='Counterexample candidate',
        vacuum_physical_matter_work=0,predecessor_nonzero_vacuum_work=True,
        maximum_spurious_energy_density_increment=float(np.max(abs(s.h*work['geometric_work']))),
        integrated_spurious_energy_erg=float(s.h*np.sum(s.volume*work['geometric_work'])),
        scalar_energy_erg=float(np.sum(s.volume*z['Ephi'])),
        scalar_residual=work['scalar_residual'],old_scheme_relative_energy_residual=work['scalar_relative_balance'])


def physical_work(star,z,trace,psi,Pi,Phi):
    measured=prior.scalar_work(star,z,trace,psi,Pi,Phi)
    old=star.reference
    b=(z['b']+star.scalar_reference['b'])/2
    bf=(z['bf']+star.scalar_reference['bf'])/2
    logH=(z['logH']+star.scalar_reference['logH'])/2
    scalarE=(Pi*Pi*b+prior.allocation(star.face_weight*Phi*Phi*bf)/star.volume)/(32*np.pi*e.GRAV)
    scalarS=-(Pi/2)*(Phi[:-1]+Phi[1:])/4*b/(4*np.pi*e.GRAV)
    Sm=(z['S']+old['S'])/2
    EmRm=(z['E']+z['R']+old['E']+old['R'])/2
    geometric=4*np.pi*e.GRAV*e.C*star.r*np.exp(logH)*(2*scalarE*Sm-scalarS*EmRm)
    physical=measured['direct_work']+geometric
    defect=z['Ephi']+star.h*e.C*star.divergence(measured['scalar_flux'])-star.h*physical
    denom=np.maximum(z['Ephi']+abs(star.h*e.C*star.divergence(measured['scalar_flux']))+abs(star.h*physical),ld('1e-100'))
    return dict(measured,scalar_work=physical,geometric_work=geometric,
        rejected_discretization_work=measured['scalar_work']-physical,
        scalar_balance=defect,scalar_relative_balance=float(np.max(abs(defect)/denom)))


class PhysicalStar(prior.CoupledStar):
    evaluate=FunctionType(prior.CoupledStar.evaluate.__code__,dict(vars(prior),scalar_work=physical_work))


def prepare():
    assert not OUT.exists()
    failed=json.loads((prior.OUT/'failure.json').read_text());assert failed['path']==1
    control=json.loads((prior.OUT/'path-0/result.json').read_text());assert control['passed']
    paths=[Path(__file__),Path(prior.__file__),prior.OUT/'plan.json',prior.OUT/'failure.json',
        prior.OUT/'path-0/result.json',prior.OUT/'path-0/endpoint.npz',prior.OUT/'path-1/failed-iterate.npz']
    OUT.mkdir()
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in paths},
        predecessor_plan=json.loads((prior.OUT/'plan.json').read_text()),vacuum_negative_control=vacuum_check(),
        correction='Use the derived continuous scalar-matter geometric cross work with midpoint scalar and averaged matter. Never turn the scalar scheme energy defect into material heating. Preserve and apply the unchanged scalar and matter acceptance gates.',
        decision='Does removing the false matter source recover the native joint solve, and what scalar/metric conservation defect remains? A native solve alone cannot pass the milestone.',
        reuse='Do not recompute the identical undriven path. One corrected driven attempt only.',
        budget=dict(workers=8,hard_timeout_seconds=600,maximum_runs=1,expected_wall_seconds=[140,450],
            basis='Undriven native path measured 142.2s. Failed paired attempt 393.5s. This one driven correction remains within the original 1200s overall wall budget; stop at 24 iterations or 600s.')))
    print('FROZEN physical exchange correction',json.dumps(vacuum_check()),flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert e.digest(e.ROOT/rel)==sha,rel
    assert not (OUT/'result.json').exists() and not (OUT/'failure.json').exists()
    start=time.monotonic();p=plan['predecessor_plan'];records=[]
    with ProcessPoolExecutor(max_workers=8,initializer=prior.imported.initial.original.worker_init) as pool:
        star=prior.initialize(pool,ld(p['boundary_amplitudes'][1]),ld(p['h_seconds']))
        star.__class__=PhysicalStar
        zero=np.zeros_like(star.base);amp=star.amplitude
        star.amplitude=ld(0);z0=star.evaluate(zero);star.amplitude=amp
        def log(row):
            records.append(row)
            with (OUT/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
            print('PHYSICAL JOINT',row['iteration'],row['residual_norm'],flush=True)
        try:
            delta,z,value=prior.stage(star,(zero,z0),log)
            np.savez_compressed(OUT/'endpoint.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
            energy=star.h*e.C*star.area*z['fluxes'][1]
            baryon=star.h*e.C*star.area*z['fluxes'][0]
            corrected=dict(z,dU=z['dU'].copy());corrected['dU'][:,1]+=star.h*z['scalar_work']
            budget,_,_=m.budget(star,corrected,energy,baryon)
            import inspect
            source=inspect.getsource(prior.two.cones).replace('def cones(', 'def frame_cones(').replace('e.TAU',"z['tau_cond']")
            ns=dict(vars(prior.two));exec(compile(source,__file__,'exec'),ns);cone=ns['frame_cones'](z)
            gates=dict(native=np.max(abs(value)/prior.ATOL)<=1,matter_budget=m.budget_passed(budget,p),
                cone=cone['sampled_cone_inside_light_cone'],scalar_residual=z['scalar_residual']<p['scalar_residual_tolerance'],
                scalar_budget=z['scalar_relative_balance']<p['scalar_relative_budget_tolerance'])
            gates={k:bool(v) for k,v in gates.items()}
            control=np.load(prior.OUT/'path-0/endpoint.npz')['delta']
            result=dict(classification='Counterexample candidate',passed=all(gates.values()),gates=gates,
                iterations=len(records),native_residual_norm=float(np.max(abs(value)/prior.ATOL)),budget=budget,cone=cone,
                scalar_residual=z['scalar_residual'],scalar_relative_balance=z['scalar_relative_balance'],
                maximum_discarded_false_heat_over_initial_heat_capacity=float(np.max(abs(star.h*z['rejected_discretization_work']/star.heat0))),
                scalar_energy_erg=float(np.sum(z['Ephi']*star.volume)),
                direct_matter_exchange_erg=float(-star.h*np.sum(z['direct_work']*star.volume)),
                geometric_matter_exchange_erg=float(-star.h*np.sum(z['geometric_work']*star.volume)),
                total_scalar_energy_defect_erg=float(np.sum(z['scalar_balance']*star.volume)),
                maximum_driven_minus_undriven=np.max(abs(delta-control),axis=0).astype(float).tolist(),
                native_calls=star.pool.evaluations,seconds=time.monotonic()-start,
                simultaneous_native_equations_solved=bool(gates['native'] and gates['scalar_residual']),
                validated_coupled_evolution=False,observational_closure=False)
            e.write(OUT/'result.json',result)
            print('PHYSICAL JOINT VERDICT',json.dumps(result),flush=True)
        except Exception as error:
            if hasattr(star,'last_state'):
                np.savez_compressed(OUT/'failed-iterate.npz',delta=star.last_delta,**{k:v for k,v in star.last_state.items() if isinstance(v,np.ndarray)})
            e.write(OUT/'failure.json',dict(classification='Counterexample candidate',error=repr(error),seconds=time.monotonic()-start))
            raise
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['vacuum_check','prepare','run'])
    print(globals()[parser.parse_args().action]())
