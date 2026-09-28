"""Counterexample candidate: energy-projected midpoint scalar integration.

The projection belongs to the numerical scalar equation, never matter heat.
It is not an unmodified midpoint solve or a continuum convergence certificate.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType
import argparse
import inspect
import json
import time
import numpy as np
from scipy.linalg import solve_banded
import sympy as sp
import def_spherical_coupled as original
import def_spherical_physical_exchange as physical

e,m,ld=original.e,original.m,original.ld
OUT=e.g.OUT/'def-spherical-balanced-wave'


def correction(star,z,trace,psi,Pi,Phi):
    geometric=original.scalar_work(star,z,trace,psi,Pi,Phi)['geometric_work']
    target=physical.physical_work(star,z,trace,psi,Pi,Phi)['geometric_work']
    coefficients=original.scalar_coefficients(star,z,trace)
    k,b=coefficients[4],coefficients[6]
    pm=Pi/2
    result=np.zeros(star.n,dtype=ld)
    nonzero=abs(pm)>ld('1e-200')
    result[nonzero]=star.h**2*np.pi*e.GRAV*e.C*k[nonzero]*(target[nonzero]-geometric[nonzero])/(b[nonzero]*pm[nonzero])
    return result


def wave(star,z,trace):
    lower,diagonal,upper,rhs,k,*_=original.scalar_coefficients(star,z,trace)
    psi=z['psi']
    Pi=2*psi/(star.h*e.C*k)
    Phi=np.r_[ld(0),np.diff(psi)/np.diff(star.r),(star.amplitude-psi[-1])/star.distance[-1]]
    rhs=rhs+correction(star,z,trace,psi,Pi,Phi)
    band=np.zeros((3,star.n));band[1]=diagonal
    band[0,1:],band[2,:-1]=upper,lower
    mid=solve_banded((1,1),band,np.asarray(rhs,float)).astype(ld)
    for _ in range(3):
        defect=rhs-diagonal*mid
        defect[1:]-=lower*mid[:-1];defect[:-1]-=upper*mid[1:]
        mid+=solve_banded((1,1),band,np.asarray(defect,float)).astype(ld)
    psi=2*mid;Pi=4*mid/(star.h*e.C*k)
    Phi=np.r_[ld(0),np.diff(psi)/np.diff(star.r),(star.amplitude-psi[-1])/star.distance[-1]]
    return psi,Pi,Phi


def work(star,z,trace,psi,Pi,Phi):
    value=physical.physical_work(star,z,trace,psi,Pi,Phi)
    lower,diagonal,upper,rhs,k,*_=original.scalar_coefficients(star,z,trace)
    projected=correction(star,z,trace,psi,Pi,Phi)
    mid=psi/2;defect=diagonal*mid-rhs-projected
    defect[1:]+=lower*mid[:-1];defect[:-1]+=upper*mid[1:]
    scale=max(abs(star.amplitude),ld('1e-100'))
    value['unprojected_scalar_residual']=value['scalar_residual']
    value['scalar_residual']=float(max(np.max(abs(defect))/scale,np.max(abs(psi-star.h*e.C*k*Pi/2))/scale))
    value['scalar_projection']=projected
    value['maximum_projection_over_boundary_amplitude']=float(np.max(abs(projected))/scale)
    return value


class BalancedStar(original.CoupledStar):
    evaluate=FunctionType(original.CoupledStar.evaluate.__code__,dict(vars(original),wave=wave,scalar_work=work))


class BalancedTangent(original.Tangent):
    def evaluate(self,delta):
        z=super().evaluate(delta)
        # Differentiate the physical exchange with the native material tangent.
        # The field/metric block remains frozen only in this iteration matrix.
        z.update(physical.physical_work(self,z,z['trace'],z['psi'],z['Pi'],z['Phi']))
        self.current=z
        return z


def tangent(star,delta,z):
    model=original.tangent(star,delta,z);model.__class__=BalancedTangent
    return model


stage=FunctionType(original.stage.__code__,dict(vars(original),tangent=tangent))


def check():
    h,g,c,k,b,p,d=sp.symbols('h g c k b p d',nonzero=True)
    corrected_rhs=h*h*sp.pi*g*c*k*d/(b*p)
    extra_Pi_rate=4*corrected_rhs/(h*h*c*k)
    assert sp.cancel(b*p*extra_Pi_rate/(4*sp.pi*g)-d)==0
    # Reuse the exact predecessor vacuum mesh, changing only the wave update.
    source=inspect.getsource(physical.vacuum_check).replace('def vacuum_check(', 'def vacuum_balanced(')
    source=source.replace('for _ in range(12):','for _ in range(12):\n        z.update(psi=locals().get("psi",zero),E=zero,R=zero,S=zero)')
    source=source.replace('s.reference=dict(E=zero,m=zero,mf=zf,a=s.b0)','s.reference=dict(E=zero,R=zero,S=zero,m=zero,mf=zf,a=s.b0)')
    source=source.replace('prior.wave(s,z,zero)','wave(s,z,zero)').replace('work=prior.scalar_work(s,z,zero,psi,Pi,Phi)','z.update(psi=psi,E=zero,R=zero,S=zero)\n    measured=work(s,z,zero,psi,Pi,Phi)\n    assert measured["scalar_residual"]<2e-17\n    assert measured["scalar_relative_balance"]<2e-15,measured["scalar_relative_balance"]\n    return measured')
    source=source[:source.index("    assert work['scalar_residual']")]
    ns=dict(vars(physical),wave=wave,work=work);exec(compile(source,__file__,'exec'),ns)
    value=ns['vacuum_balanced']()
    assert np.all(value['scalar_work']==0)
    return dict(classification='Counterexample candidate',passed=True,vacuum_material_work_exactly_zero=True,
        scalar_residual=value['scalar_residual'],scalar_relative_balance=value['scalar_relative_balance'],
        projection_over_amplitude=value['maximum_projection_over_boundary_amplitude'],
        identity=dict(classification='Proven',passed=True,scope='The added scalar numerical increment supplies exactly the difference between physical cross work and the staggered geometric product remainder. This does not prove PDE consistency or continuum convergence.'))


def prepare():
    assert not OUT.exists()
    checked=check()
    paths=[Path(__file__),Path(original.__file__),Path(physical.__file__),original.OUT/'plan.json',
        original.OUT/'path-0/result.json',original.OUT/'path-0/endpoint.npz',physical.OUT/'plan.json']
    p=json.loads((original.OUT/'plan.json').read_text());OUT.mkdir()
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in paths},
        original_plan=p,check=checked,
        method='Energy-projected midpoint scalar equation: add h^2*pi*g*c*k*(physical_geometric_work-discrete_geometric_remainder)/(b*Pi_mid) to the scalar midpoint RHS. Zero at exactly zero fields; extended precision residual refinement. Physical matter work unchanged. The scalar numerical scheme changes explicitly; original failed equations remain frozen.',
        gates='Same matter, cone and finite energy gates. The scalar algebra residual now tests the explicitly modified update; also record the unprojected wave defect and projection size. Neither small defect nor a passed finite budget proves continuum convergence.',
        budget=dict(workers=8,hard_timeout_seconds=450,maximum_runs=1,expected_wall_seconds=[140,400],
            cumulative_cap_seconds=1200,basis='Original attempt 393.5s, physical-only correction capped at 600s; launch only if completed elapsed totals plus 450s remain below 1200s. Reuse accepted zero-drive state; do not refine or extend time.'),
        decision='Can a numerical scalar energy correction eliminate the false vacuum matter source and retain a native simultaneous driven solve under unchanged finite matter/energy thresholds?',
        continuum_consistency_proved=False,time_convergence_verified=False,observational_closure=False))
    print('FROZEN balanced scalar update',json.dumps(checked),flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text());p=plan['original_plan']
    for rel,sha in plan['bindings'].items():assert e.digest(e.ROOT/rel)==sha,rel
    assert not (OUT/'result.json').exists() and not (OUT/'failure.json').exists()
    elapsed=json.loads((original.OUT/'failure.json').read_text())['seconds']
    previous=physical.OUT/('result.json' if (physical.OUT/'result.json').exists() else 'failure.json')
    elapsed+=json.loads(previous.read_text())['seconds']
    assert elapsed+plan['budget']['hard_timeout_seconds']<plan['budget']['cumulative_cap_seconds'],elapsed
    start=time.monotonic();records=[]
    with ProcessPoolExecutor(max_workers=8,initializer=original.imported.initial.original.worker_init) as pool:
        star=original.initialize(pool,ld(p['boundary_amplitudes'][1]),ld(p['h_seconds']));star.__class__=BalancedStar
        zero=np.zeros_like(star.base);amp=star.amplitude;star.amplitude=ld(0)
        z0=star.evaluate(zero);star.amplitude=amp
        def log(row):
            records.append(row)
            with (OUT/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
            print('BALANCED JOINT',row['iteration'],row['residual_norm'],flush=True)
        try:
            delta,z,value=stage(star,(zero,z0),log)
            np.savez_compressed(OUT/'endpoint.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
            corrected=dict(z,dU=z['dU'].copy());corrected['dU'][:,1]+=star.h*z['scalar_work']
            budget,_,_=m.budget(star,corrected,star.h*e.C*star.area*z['fluxes'][1],star.h*e.C*star.area*z['fluxes'][0])
            source=inspect.getsource(original.two.cones).replace('def cones(', 'def frame_cones(').replace('e.TAU',"z['tau_cond']")
            ns=dict(vars(original.two));exec(compile(source,__file__,'exec'),ns);cone=ns['frame_cones'](z)
            gates=dict(native=np.max(abs(value)/original.ATOL)<=1,matter_budget=m.budget_passed(budget,p),
                cone=cone['sampled_cone_inside_light_cone'],scalar_residual=z['scalar_residual']<p['scalar_residual_tolerance'],
                scalar_budget=z['scalar_relative_balance']<p['scalar_relative_budget_tolerance'])
            gates={k:bool(v) for k,v in gates.items()}
            control=np.load(original.OUT/'path-0/endpoint.npz')['delta']
            total=z['dEm']+z['Ephi']+star.h*e.C*star.divergence(z['fluxes'][1]+z['scalar_flux'])
            allowance=star.heat0*ld('1e-9')+ld('2e-15')*(z['Ephi']+abs(star.h*e.C*star.divergence(z['scalar_flux'])))
            gates['total_energy']=bool(np.max(abs(total)/allowance)<=1)
            result=dict(classification='Counterexample candidate',passed=all(gates.values()),gates=gates,
                iterations=len(records),native_residual_norm=float(np.max(abs(value)/original.ATOL)),budget=budget,cone=cone,
                scalar_residual=z['scalar_residual'],scalar_relative_balance=z['scalar_relative_balance'],
                unprojected_scalar_residual=z['unprojected_scalar_residual'],
                scalar_projection_over_amplitude=z['maximum_projection_over_boundary_amplitude'],
                total_energy_allowance_ratio=float(np.max(abs(total)/allowance)),
                scalar_energy_erg=float(np.sum(z['Ephi']*star.volume)),
                direct_matter_exchange_erg=float(-star.h*np.sum(z['direct_work']*star.volume)),
                geometric_matter_exchange_erg=float(-star.h*np.sum(z['geometric_work']*star.volume)),
                maximum_driven_minus_undriven=np.max(abs(delta-control),axis=0).astype(float).tolist(),
                native_calls=star.pool.evaluations,seconds=time.monotonic()-start,
                simultaneous_discrete_step_solved=True,continuum_consistency_proved=False,
                time_convergence_verified=False,companion_matched=False,observational_closure=False)
            e.write(OUT/'result.json',result)
            print('BALANCED JOINT VERDICT',json.dumps(result),flush=True)
        except Exception as error:
            if hasattr(star,'last_state'):
                np.savez_compressed(OUT/'failed-iterate.npz',delta=star.last_delta,**{k:v for k,v in star.last_state.items() if isinstance(v,np.ndarray)})
            e.write(OUT/'failure.json',dict(classification='Counterexample candidate',error=repr(error),seconds=time.monotonic()-start))
            raise
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['check','prepare','run'])
    print(globals()[parser.parse_args().action]())
