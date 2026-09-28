"""Time-centered conservative transport with the exponential scalar carrier.

The stiff two-carrier heat closure remains backward Euler. No claim of a
globally second-order method is made. Frozen Phase44/45 code is not mutated.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import inspect
import json
import numpy as np
import sympy as sp
import def_exponential_coupled as x
import def_exponential_coupled_run as prior_run
import def_exponential_coupled_audit as prior_audit

s,e,ld=x.s,x.e,x.ld
OUT=e.g.OUT/'def-centered-matter'


def species(star,candidate,previous,older,h,coefficients,metric):
    p=previous[1];f=s.m.fluxes(star,p)[0]
    X=star.base[:,5:]+previous[0][:,5:]
    flux=f[:,None]*s.m.prior.upwind(X,f)
    B=p['B']-h*e.C/2*star.divergence(f)
    BX=p['B'][:,None]*X-h*e.C/2*np.diff(star.area[:,None]*flux,axis=0)/star.volume[:,None]
    assert B.min()>0 and BX.min()>=-1e-10
    effective=previous[0].copy();effective[:,5:]=BX/B[:,None]-star.base[:,5:]
    state=dict(p,B=B)
    return s.m.species(star,candidate,(effective,state),(effective,state),h/2,(ld(1),ld(-1),ld(0)),metric)


def residual(star,delta,previous,older,h,coefficients):
    value,z=s.residual(star,delta,previous,older,h,coefficients)
    p=previous[1];old=s.m.fluxes(star,p);new=z['fluxes']
    value[:,0]+=h*e.C/2*star.divergence(old[0]-new[0])/star.B0
    value[:,1]+=h*e.C/2*star.divergence(z['Hface']*(old[1]-new[1]))/(z['H']*star.heat0)
    def force(a):
        result=-z['at']*a['S']-e.C*a['N']*a['nur']*a['E']
        result+=e.C*a['N']*a['P']*star.area_difference_over_volume
        result+=e.C*a['N']*star.beta*a['psi']*a['trace']*(a['Phi'][:-1]+a['Phi'][1:])/2
        return result
    value[:,2]=(z['dU'][:,2]-p['dU'][:,2]+h*e.C*star.divergence((old[2]+new[2])/2)-h*(force(p)+force(z))/2)/star.momentum_scale
    XP=star.base[:,5:]+previous[0][:,5:];XN=star.base[:,5:]+delta[:,5:]
    fp=old[0][:,None]*s.m.prior.upwind(XP,old[0]);fn=new[0][:,None]*s.m.prior.upwind(XN,new[0])
    value[:,5:]+=h*e.C/2*np.diff(star.area[:,None]*(fp-fn),axis=0)/(star.volume*star.B0)[:,None]
    return value,z


_m=dict(vars(s.m));_m['species']=species;material=SimpleNamespace(**_m)


class Tangent(s.Tangent):
    tangent_aux=FunctionType(s.old.CoupledStar.tangent_aux.__code__,dict(s.old.CoupledStar.tangent_aux.__globals__,m=material))


def tangent(star,delta,z):
    model=Tangent.__new__(Tangent);model.__dict__=dict(star.__dict__)
    model.anchor,model.linearization=delta.copy(),z
    model.projected_anchor=species(star,delta,*star.composition_context,z)
    return model


jacobian=FunctionType(s.m.prior.jacobian.__code__,dict(s.m.prior.jacobian.__globals__,residual=residual))
stage=FunctionType(s.stage.__code__,dict(s.stage.__globals__,m=material,residual=residual,tangent=tangent,jacobian=jacobian))
_s=dict(vars(s));_s.update(residual=residual,stage=stage);solver=SimpleNamespace(**_s)


def symbolic():
    old,new,h,fp,fn=sp.symbols('old new h fp fn')
    assert sp.expand(new-old+h*(fp+fn)/2-(new-(old-h*fp/2)+h*fn/2))==0
    g0,g1,g2=sp.symbols('g0 g1 g2')
    assert sp.expand((g1-g0)+(g2-g1)-(g2-g0))==0
    return dict(classification='Proven',passed=True,
        scope='Explicit half-flux plus implicit half-flux is the same conservative trapezoidal update; shared species faces telescope. No second-order claim for the backward-Euler heat closure or paired lapse.')


def controls():
    return dict(classification='Counterexample candidate',passed=True,exponential=x.controls(),symbolic=symbolic())


# Reuse the tested pilot and result protocol, but bind the actual new residual
# and composition iteration. This does not mutate any imported module.
source=inspect.getsource(x.pilot).replace('files=[Path(__file__),','files=[Path(x.__file__),Path(__file__),')
source=source.replace("checkpoint='2ba5f1ea'","checkpoint='c840ee93'")
source=source.replace('First native reciprocal coupling of the exponential scalar carrier and boundary energy reservoir; no defect assigned to material heat.',
    'Time-center shared baryon, energy, momentum and isotope transport plus momentum force. Retain the exponential reciprocal scalar and both stiff backward-Euler heat laws.')
namespace=dict(vars(x),x=x,OUT=OUT,s=solver,__file__=__file__,controls=controls)
exec(compile(source,__file__,'exec'),namespace);pilot=namespace['pilot']

_x=dict(vars(x));_x['OUT']=OUT
runtime=dict(vars(prior_run),x=SimpleNamespace(**_x),s=solver,OUT=OUT/'production',__file__=__file__)
for name in ['bindings','run_path','compare','run','prepare']:
    runtime[name]=FunctionType(getattr(prior_run,name).__code__,runtime)


def prepare():
    runtime['prepare']()
    path=OUT/'production/plan.json';plan=json.loads(path.read_text())
    plan['bindings'][Path(prior_run.__file__).relative_to(s.ROOT).as_posix()]=e.digest(Path(prior_run.__file__))
    plan['method']='Time-centered matter transport/force and the same conservative species faces; exponential reciprocal scalar with boundary reservoir. Both stiff heat closures remain backward Euler; no global second-order claim.'
    plan['comparison']='Same failed Phase45 input, grids, cutoff readouts and original gates; only temporal material quadrature changes. Two new-method pilot steps reused, no old-method state reused.'
    e.write(path,plan)


def audit():
    target=OUT/'production/saved-audit.json';assert not target.exists()
    globals_=dict(prior_audit.audit.__globals__,run=SimpleNamespace(**runtime),x=x,s=solver,__file__=__file__)
    result=FunctionType(prior_audit.audit.__code__,globals_)();result['centered_species_symbolic']=symbolic()
    e.write(target,result);print(json.dumps(result));assert result['passed']


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['pilot','prepare','run','audit'])
    action=parser.parse_args().action
    runtime['run']() if action=='run' else globals()[action]()
