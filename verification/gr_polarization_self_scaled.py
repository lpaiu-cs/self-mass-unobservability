"""Resolve the occupied momentum region at large Q, without a response floor."""
from types import FunctionType
import json, shutil, sys
import numpy as np
import mpmath as mp
from scipy.integrate import quad_vec
from scipy.special import expit
import gr_polarization_self_refinement as previous

original=previous.original;g=original.g;OUT=g.OUT/'gr-polarization-self-scaled'
BASE=original.finite.module()


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=original.finite.previous.original.previous.scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();previous.bindings()
    error=g.ROOT/'outputs/gr-polarization-self-refinement33-run.log'
    assert 'assert np.all(r>0)' in error.read_text() and not (previous.OUT/'states.npz').exists()
    shutil.copy2(error,OUT/'original-error.log')
    paths=[g.ROOT/'verification/gr_polarization_self_scaled.py',previous.OUT/'plan.json',
        previous.OUT/'order-256.npz',previous.OUT/'symbolic.json',previous.OUT/'tail-bounds.json',OUT/'original-error.log']
    save('plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        original_failure='The 512-point run reached its highest-wavenumber nodes, where the old u=p/Q quadrature returned zero for a localized occupied momentum distribution. The strict positive-response assertion stopped it before publishing an order512 or final result.',
        change='When s=max_i(sqrt(beta_i*t_i*(2+beta_i*t_i))/Q_i)<1/16, t_i=max(eta_i,0)+1, map the low-u interval by u=s*w. Integrate w on [0,16], then u on [16*s,1] and [1,infinity], assigning one third of the original tolerance to each. No discarded interval or positive floor. Otherwise reuse the original response unchanged. Recompute 256/512 points with unchanged science gates.',
        response_grid_absolute_gate=1e-8,high_Q_control_relative_gate=1e-7,
        controls='Replay all41678 old response grid points. At the highest order512 wavenumber, preserve old zero responses and compare new positive responses to 60-digit cold-function Fermi convolution at the five registered control positions.',
        physical_EOS_certified=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return previous.bindings()


def response(eta,beta,Q,S0,tolerance):
    if np.all(Q==0):return BASE['response'](eta,beta,Q,S0,tolerance)
    t=np.maximum(eta,0)+1;s=float(np.max(np.sqrt(beta*t*(2+beta*t))/Q))
    if s>=1/16:return BASE['response'](eta,beta,Q,S0,tolerance)
    def f(u):
        p=Q*u;gamma=np.sqrt(1+p*p);t=p*p/(beta*(gamma+1))
        return BASE['positive_kernel'](u,Q)*expit(eta-t)/S0
    options=dict(epsabs=tolerance/3,epsrel=tolerance/3,norm='max',limit=1600)
    return (quad_vec(lambda w:s*f(s*w),0,16,**options)[0]
        +quad_vec(f,16*s,1,**options)[0]+quad_vec(f,1,np.inf,**options)[0])


def audit():
    plan=bindings();own=json.loads((OUT/'plan.json').read_text())
    a0=dict(np.load(original.finite.OUT/'states.npz'));eta,beta,S0=[a0[k] for k in ['eta','beta','S0']]
    errors=[]
    for j in range(len(a0['k_over_kappa'])):
        value=response(eta,beta,a0['Q'][:,j],S0,plan['inner_tolerance'])
        errors.append(float(np.max(abs(value-a0['response_over_k0'][:,j]))))
    end=np.arctan(float(plan['z_cut']));nodes,_=np.polynomial.legendre.leggauss(512);z=float(np.tan((nodes[-1]+1)*end/2))
    a=a0['Q'][:,int(np.flatnonzero(a0['k_over_kappa']==1)[0])];Q=a*z
    old=BASE['response'](eta,beta,Q,S0,plan['inner_tolerance']);new=response(eta,beta,Q,S0,plan['inner_tolerance'])
    assert np.any(old==0) and np.all(new>0)
    controls=[];mp.mp.dps=60
    for i in plan['control_positions']:
        ref=BASE['convolution'](eta[i],beta[i],Q[i],S0[i]);score=abs(new[i]/ref-1)
        controls.append(dict(cell=int(a0['cells'][i]),Q=float(Q[i]),old=float(old[i]),new=float(new[i]),reference=ref,relative_difference=float(score),passed=score<own['high_Q_control_relative_gate']))
    passed=max(errors)<own['response_grid_absolute_gate'] and all(v['passed'] for v in controls)
    save('response-audit.json',dict(classification='Counterexample candidate',passed=passed,
        previous_grid_maximum_absolute_difference=max(errors),grid_points=a0['response_over_k0'].size,
        largest_k_over_kappa=z,old_zero_states=int(np.sum(old==0)),new_zero_states=int(np.sum(new==0)),controls=controls,
        scope='Actual old quadrature failure reproduced at fixed high Q, resolved by exact change of variable and independently checked. No floor or epsilon shift.'))
    assert passed


def run():
    audit()
    # Reuse the original orchestration, replacing only its response function.
    model=dict(BASE,response=response)
    class BoundModel:
        OUT=original.finite.OUT
        @staticmethod
        def module():return model
    ns=dict(original.__dict__,OUT=OUT,save=save,bindings=bindings,verify=verify,finite=BoundModel)
    for name in ['symbolic','tail_bounds','run']:ns[name]=FunctionType(original.__dict__[name].__code__,ns)
    ns['run']()


def verify():
    ns=dict(original.__dict__,OUT=OUT,bindings=bindings)
    FunctionType(original.verify.__code__,ns)()
    assert json.loads((OUT/'response-audit.json').read_text())['passed']
    for name in ['symbolic.json','tail-bounds.json']:assert (original.OUT/name).read_bytes()==(OUT/name).read_bytes(),name
    print('PASS occupied-momentum scaling; original failures and all analytic bounds retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
