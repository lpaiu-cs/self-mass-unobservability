"""Native-state and independent algebra controls for the repaired GR source."""
from pathlib import Path
import json
import signal
import time
import numpy as np
import sympy as s
import def_native_metric_release as task
import def_native_metric_charge as charge


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic();signal.alarm(30)
    task.write(out/'audit-plan.json',dict(classification='Counterexample candidate',seconds=30,native_calls=12,
        checks=['actual corrected boundary cells and dilute gas against native EOS','lab-frame E-Pr boost identity',
        'general GR source vacuum reduction','retarded mass-source Green integral for a manufactured linear history','initial face acceleration after the final local repair'],
        native_relative_gate=.002,algebra_scope='The scalar wave and GR constraints are linearized around the declared static isotropic background. No chemical kinetics or nonlinear feedback certificate.'))
    r,m,N,Phi,A,alpha,E,P=s.symbols('r m N Phi A alpha E P',nonzero=True)
    b=1-2*m/r;C2=N*N*b
    KJ=2*C2*Phi/(r*b*b)*(1+4*s.pi*r*r*A**4*(P-E))-8*s.pi*N*N*alpha*A**4*(E-3*P)/b
    assert s.simplify(KJ.subs({E:0,P:0})-2*N*N*Phi/(r*b))==0
    fcoef=KJ*r*r*b*Phi+16*s.pi*N*N*alpha*r*r*A**4*Phi*(P-E)-4*s.pi*N*N*r*A**4*(s.Symbol('beta')+4*alpha**2)*(E-3*P)
    V=C2*(2*m/(r*r*b)+4*s.pi*r*A**4*(P-E)/b)/r-fcoef/r
    assert s.simplify(V.subs({E:0,P:0})-2*N*N*(m/r**3-Phi**2))==0
    rho,u,p,v,cx=s.symbols('rho u p v cx')
    labE=(rho*(cx+u)+p)/(1-v*v)-p
    labPr=(rho*(cx+u)+p)*v*v/(1-v*v)+p
    assert s.simplify(labE-labPr-(rho*(cx+u)-p))==0
    fan=task.old.bank.task.prior.Fan(call_cap=12,reuse=True);rows=[];controls=[]
    for n in [896,1792]:
        flow=task.Flow(n);base=flow.base;d=np.load(out/f'cells-{n}.npz')
        active=np.flatnonzero(d['rho']>base.eos.rho0*base.eos.floor)
        ids=np.unique([0,1,int(active[np.argmin(abs(np.log(d['rho'][active]/base.eos.rho0)+5))]),int(active[-1])])
        samples=[]
        for i in ids:
            raw=fan.call(float(np.log(d['rho'][i])),float(np.log(d['T'][i])))
            rr=d['rho'][i]/flow.eos.rho0;vv=d['velocity_cm_s'][i]/task.C;root=np.sqrt(1-vv*vv);W=1/root;wm=vv*vv/(root*(1+root))
            pp=d['pressure'][i]/(flow.eos.rho0*task.C**2)
            tau=(d['U'][2,i]-(base.a[i]-base.a0)*flow.eos.cx*d['U'][0,i])/base.a[i]
            actual_u=(tau-flow.eos.cx*d['U'][0,i]*wm-pp*W*W*vv*vv)/(rr*W*W)*task.C**2
            errors=[float(abs(d['pressure'][i]/raw[1]-1)),float(abs(actual_u/raw[2]-1))]
            samples.append(dict(cell=int(i),relative=errors))
        rate,_,_=flow.rhs(flow.initial,0.)
        pp,uu,*_=flow.eos(flow.initial[0],flow.seed);hh=flow.eos.cx+uu+pp/np.maximum(flow.initial[0],flow.eos.floor)
        accel=rate[1,:8]/(flow.initial[0,:8]*hh[:8])*task.C
        rows.append(dict(cells=n,native_controls=samples,initial_acceleration_cm_s2=accel.tolist()))
        source=np.load(out/'charge'/f'mass-source-{n}.npz');model=charge.adapter(n);geom=charge.prior.Geometry(model)
        _,delay,a,B,_=geom(model.x);coefficient=source['scalar_source_coefficient_cm3']
        end=float(source['source_seconds'][-1]);query=end/2
        poly=charge.prior.polynomial(np.array([0.,end]),np.array([np.zeros(n),np.full(n,end)])).antiderivative()
        numerical=task.C/2*np.sum(model.dx*B/a*coefficient*charge.prior.paired(poly,query+delay))
        exact=task.C/4*np.sum(model.dx*B/a*coefficient*np.maximum(query+delay,0)**2)
        controls.append(float(abs(numerical/exact-1)))
    passed=max(max(x['relative']) for row in rows for x in row['native_controls'])<.002 and max(controls)<1e-12
    result=dict(classification='Counterexample candidate',passed=bool(passed),rows=rows,
        symbolic=dict(classification='Proven',vacuum_scalar_metric_reduction=True,boost_invariant_stress=True,scope='Identities in the declared linearized model.'),
        manufactured_Green_relative=controls,native_calls=fan.calls,seconds=time.monotonic()-start,
        whole_GR_solution_certified=False,full_goal_complete=False)
    task.write(out/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert passed


if __name__=='__main__':main()
