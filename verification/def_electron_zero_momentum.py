"""Exploit exact isotropy at p1=0 and sample the partner Fermi sea directly.

The near-Fermi logistic proposal used for the bulk bracket is inefficient
for a deep-hole loss rate. Use two-dimensional deterministic quadrature,
with the same energy domain and fewer total collision evaluations.
"""
from pathlib import Path
import json
import time
import numpy as np
import sympy as sp
from scipy.special import expit
from numpy.polynomial.legendre import leggauss
import def_electron_screened_fermi as resolved

exchange=resolved.exchange
OUT=exchange.OUT/'zero-momentum'


def integrate(state,order):
    gx,gw=leggauss(order);u=(gx+1)/2;uw=gw/2
    # Split the two branches of the screened angular proposal at 1/2.
    angle=np.r_[u/2,.5+u/2];angle_weight=np.tile(uw/2,2)
    a,b=np.meshgrid(u,angle,indexing='ij');weight=np.outer(uw,angle_weight).ravel()
    samples=np.c_[a.ravel(),a.ravel(),np.full(a.size,.5),b.ravel(),np.full(a.size,.37)]
    total=0.;maximum=dict(energy=0.,momentum=0.,detailed_balance_log=0.)
    for start in range(0,len(samples),4096):
        value,_,checks=exchange.events(samples[start:start+4096],state,fixed_p=0.)
        total+=weight[start:start+4096]@value
        for key in maximum:maximum[key]=max(maximum[key],checks[key])
    return float(total*exchange.RATE),maximum


def main():
    assert not OUT.exists();OUT.mkdir();h=exchange.h;source=resolved.OUT
    plan=json.loads((source/'plan.json').read_text())
    paths=[Path(__file__),Path(resolved.__file__),exchange.OUT/'collision-readout/plan.json',exchange.OUT/'collision-readout/result.json',source/'result.json']
    exchange.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        preserved_failure='The 5D near-Fermi importance readout failed the 3% gate at cell 5734 (3.138%). Its deep-hole loss rate under-samples partner energies above the Fermi surface. Preserve its raw rates as unaccepted, not an accurate p=0 result.',
        repair='At p1=0, rotational invariance removes the two orientation integrals. Directly integrate the partner kinetic energy x2 in [0,eta+40] and the screened scattering-angle proposal, including the exact replacement Jacobian. Same kernel and energy domain; no physical relaxation floor.',
        cells=plan['cells'],orders=[64,128],gates=dict(relative_quadrature=.03,isotropy_relative=2e-12),
        budget=dict(hard_seconds=60,event_evaluations=3*(2*64**2+2*128**2),CPU_workers=1,native_calls=0,stellar_steps=0,automatic_expansion=False)))
    # Generalize only the importance proposal for the partner. The two first
    # sample coordinates are equal; p1 itself is overwritten by fixed_p=0.
    def uniform_energy(u,state):
        top=state['eta']+40;x=top*u;f=expit(state['eta']-x)
        return x,f,top*f*(1-f)
    exchange.energy=uniform_energy
    f1,f2,f3,f4,L=sp.symbols('f1 f2 f3 f4 L',positive=True)
    assert sp.simplify(L*f2*(1-f2)*(1-f3)*(1-f4)/((1-f1)*(1-f2))-L*f2*(1-f3)*(1-f4)/(1-f1))==0
    rows=[];start=time.monotonic()
    for index in plan['cells']:
        state=exchange.equilibrium(index);d=np.load(source/f'polarization-{index}.npz')
        exchange.amplitude=resolved.leading.amplitude(state,d['phase'],d['polarization'])
        u=np.linspace(.03,.97,32);base=np.c_[u,u,np.full(32,.5),np.linspace(.01,.99,32),np.full(32,.37)]
        a,_,_=exchange.events(base,state,fixed_p=0.)
        alternate=base.copy();alternate[:,2]=np.linspace(.1,.9,32);alternate[:,4]=np.linspace(.05,.95,32)
        b,_,_=exchange.events(alternate,state,fixed_p=0.)
        isotropy=float(np.max(abs(a/b-1)));assert isotropy<2e-12
        coarse,_=integrate(state,64);fine,checks=integrate(state,128)
        rows.append(dict(cell=index,loss_rate_per_second=fine,relative_quadrature=abs(coarse/fine-1),isotropy_relative=isotropy,conservation=checks))
        print('CELL',index,'LOSS0',fine,'QUADRATURE',rows[-1]['relative_quadrature'],flush=True)
    result=dict(classification='Counterexample candidate',rows=rows,
        numerical_gates_passed=all(r['relative_quadrature']<.03 for r in rows),
        first_readout_loss_rates_superseded=True,continuum_inverse_moment_certified=False,
        elapsed_seconds=time.monotonic()-start,full_dynamic_charge_solved=False)
    exchange.write(OUT/'result.json',result)


if __name__=='__main__':main()
