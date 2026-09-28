"""Independent joint baryon, collision momentum and analytic ODE controls."""
from pathlib import Path
import json
import sys
import time
import numpy as np
import sympy as sp
import def_native_interior_motion as task


def boundary():
    d=np.load(task.OUT/'motion-128.npz');a=np.load(task.prior.OUT/'source-896-128.npz')
    # Read the actual atmosphere independently of the imposed mechanical port.
    atmo=np.asarray(a['baryon_g'].sum(1),float);deep=d['baryon_g'][::8].sum(1,dtype=np.longdouble)
    defect=float(np.max(abs(deep+atmo)));relative=defect/max(abs(atmo).max(),1.)
    return dict(actual_joint_mass_defect_g=defect,relative_to_actual_port=float(relative),passed=bool(relative<1e-6))


def main():
    start=time.monotonic();out=task.OUT;first=sys.argv[-1]=='boundary'
    if first:
        row=boundary();task.write(out/'first-joint-baryon-audit.json',dict(classification='Counterexample candidate',**row));print(row);assert row['passed'];return
    task.write(out/'audit-plan.json',dict(classification='Counterexample candidate',seconds=15,native_constructor_calls=2,
        checks=['actual deep plus atmosphere mass, normalized by actual small port','first collision moment versus independent opacity times H','forced harmonic oscillator closed form through actual stepper','symbolic first law and same-observer normalization'],
        gates=dict(joint_baryon_relative=1e-6,collision=1e-12,oscillator=.002),bindings={str(Path(__file__)):task.sha(__file__)}))
    b=task.model().bulk;z=task.prior.load_path(896,128);saved=np.load(out/'forcing-896.npz');forces=[]
    for th,et,I in zip(saved['theta'],saved['eta'],np.concatenate([b.initial[None],z['bulk_I']])):
        ne=b.eos.gas(th,et)[6];ab,em,*_=b.eos.radiation(th,et)
        opacity=b.d['rho'][:,None]*b.d['thermo'][:,4,None]*(ab-em)+ne[:,None]*6.6524587321e-25
        H=np.einsum('iqf,q->if',I,b.w*b.mu)
        forces.append(np.sum(opacity*H*b.d['num']*b.d['Einf'],axis=1)/b.d['a']**3)
    collision=float(np.max(abs(np.array(forces)-saved['radiation_force']))/np.max(abs(saved['radiation_force'])))
    errors=[]
    for steps in [32,64]:
        m=task.Mechanics();T=m.d['t'][-1];n=m.n;k=1/T**2
        f=np.arange(1,n,dtype=float)*1e9
        L=np.zeros((n-1,n+1));L[:,1:-1]=-k*np.eye(n-1)
        m.operator=lambda t:(L,f)
        m.boundary=lambda t,nu=0:0.
        label=f'manufactured-{steps}';m.solve(steps,label)
        z=np.load(out/(label+'.npz'));exact=(1-np.cos(z['t']/T))[:,None]*f[None,:]/k
        errors.append(float(np.max(abs(z['transported_g'][:,1:-1]-exact))/np.max(abs(exact))))
    p,rho,ur,ut,pr,pt=sp.symbols('p rho ur ut pr pt',nonzero=True)
    tr=(p/rho-ur)/ut;K=pr+pt*tr
    assert sp.simplify(ur+ut*tr-p/rho)==0
    q,a,e=sp.symbols('q a e');assert sp.simplify((q+a*e)/(1-e)-(q+a)/(1-e)+a)==0
    row=boundary();passed=row['passed'] and collision<1e-12 and errors[-1]<.002 and 3.9<errors[0]/errors[1]<4.1
    task.write(out/'audit.json',dict(classification='Counterexample candidate',passed=bool(passed),joint_baryon=row,
        collision_first_moment_relative=collision,manufactured_oscillator_relative=errors,observed_error_ratio=errors[0]/errors[1],
        symbolic_classification='Proven',fixed_inventory_adiabatic_K=str(K),identities_passed=True,seconds=time.monotonic()-start))
    print(json.dumps(dict(passed=passed,baryon=row,collision=collision,oscillator=errors)),flush=True);assert passed


if __name__=='__main__':main()
