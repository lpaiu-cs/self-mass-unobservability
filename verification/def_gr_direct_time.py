"""Full coupled Cayley propagation with exact heat-pole forcing moments.

No modal reduction. The homogeneous midpoint map has unit-modulus multipliers
for positive stiffness; the original ramp and microscopic startup are integrated
in each forcing moment. Step-size comparisons remain mandatory.
"""
from pathlib import Path
import argparse
import json
import time
import signal
import resource
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import def_gr_hierarchical as space

task=space.task;OUT=task.OUT/'direct-time';FIELDS=task.FIELDS;write=task.write


class Stepper:
    def __init__(self,model,steps):
        self.model=model;self.steps=steps;self.dt=h=1/steps
        M=model.M;A=(M+h*h/4*model.K).tocsc()
        self.scale=np.sqrt(A.diagonal());D=diags(1/self.scale)
        self.lu=splu((D@A@D).tocsc());self.extended=A.astype(np.longdouble);self.error=0.
        self.mass_scale=np.sqrt(M.diagonal());D=diags(1/self.mass_scale)
        mass=splu((D@M@D).tocsc())
        heat=model.heat;tc=model.original.radiation.geometry.tc;lam=heat.rates*tc
        assert lam.min()*h>700
        moments=np.zeros((len(heat.edges),4))
        for k in range(4):moments[heat.face_ids,k]=np.sum(heat.amplitude*tc/lam**k,axis=1)
        self.force=model.load@moments
        self.acceleration=mass.solve(self.force/self.mass_scale[:,None])/self.mass_scale[:,None]
        self.q=np.zeros(model.size);self.v=self.q.copy();self.n=0
        self.omitted_exponential_bound=float(np.exp(-700))

    def step(self):
        h=self.dt;t=self.n*h;F=self.force;B=self.acceleration
        # I1 = integral_0^h (h-s) F(t+s) ds. Beyond the first step all
        # exponential tails are smaller than exp(-700) of their coefficients.
        integral=F[:,0]*(t*h*h/2+h**3/6)-F[:,1]*(h*h/2)
        correction=B[:,0]*(h*h/6)
        if self.n==0:
            integral+=h*F[:,2]-F[:,3]
            correction+=-B[:,2]+2/h*B[:,3]
        rhs=self.model.M@(2*self.q+h*self.v)+integral
        r=self.lu.solve(rhs/self.scale)/self.scale
        extended_rhs=rhs.astype(np.longdouble)
        for _ in range(2):
            defect=extended_rhs-self.extended@r.astype(np.longdouble)
            r+=self.lu.solve(np.asarray(defect/self.scale,float))/self.scale
        defect=extended_rhs-self.extended@r.astype(np.longdouble)
        error=float(np.max(abs(defect)/(abs(self.extended)@abs(r)+abs(rhs)+1e-100)))
        self.error=max(self.error,error);assert self.error<1e-9
        newq=r-self.q
        self.v=2/h*(newq-self.q)-self.v+correction
        self.q=newq;self.n+=1


def evolve(model,steps,label):
    start=time.monotonic();p=Stepper(model,steps);rows=[];velocities=[];scalars=[]
    r=model.original.native;N,a=model.original.radiation.geometry.metric(r)
    speed=a/N*r*task.fem.base.task.h.gr.C;w=model.original.weights
    for n in range(steps+1):
        if n%(steps//64)==0:
            t=n/steps;flux,energy=model.heat.faces(t)
            q=p.q-model.H@energy;qt=p.v-model.H@(flux*model.original.radiation.geometry.tc)
            cv=speed*(model.nativeV[0]@qt);cf=model.nativeV[1]@q
            rows.append(dict(tau=t,velocity_mass_RMS_m_s=float(np.sqrt(w@(cv*cv))),scalar_mass_RMS=float(np.sqrt(w@(cf*cf))),
                **{key:float(np.sqrt(w[mask]@cv[mask]**2/w[mask].sum())) for key,mask in zip(FIELDS[2:],model.original.masks)}))
            velocities.append(cv);scalars.append(cf)
        if n<steps:p.step()
    np.savez_compressed(OUT/(label+'.npz'),grid=model.grid,cells=model.cells,indices=model.indices,q=q,qt=qt,
        native_radius=r,native_velocity=velocities,native_scalar=scalars,weights=w,masks=model.original.masks,heat_energy=energy,heat_flux=flux)
    result=dict(classification='Counterexample candidate',degree=model.degree,steps=steps,history=rows,linear_residual=p.error,
        heat_telescoping=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100)),
        seconds=time.monotonic()-start,omitted_exponential_coefficient_factor=p.omitted_exponential_bound)
    write(OUT/(label+'.json'),result);print('DIRECT TIME',label,result['seconds'],rows[-1],flush=True);return result


def control():
    import sympy as s
    x,h,t,l=s.symbols('x h t l',positive=True)
    source=t+x-1/l+s.exp(-l*(t+x))/l
    I0=t*h+h*h/2-h/l+s.exp(-l*t)*(1-s.exp(-l*h))/l**2
    I1=t*h*h/2+h**3/6-h*h/(2*l)+s.exp(-l*t)*(h/l**2-(1-s.exp(-l*h))/l**3)
    assert s.simplify(s.integrate(source,(x,0,h))-I0)==0
    assert s.simplify(s.integrate((h-x)*source,(x,0,h))-I1)==0
    errors=[]
    for steps in [128,256,512]:
        dt=1/steps;q=0.;v=0.;omega=17.
        for n in range(steps):
            r=(2*q+dt*v+(n*dt)*dt*dt/2+dt**3/6)/(1+dt*dt*omega**2/4)
            newq=r-q;v=2/dt*(newq-q)-v+dt*dt/6;q=newq
        errors.append(abs(q-(1-np.sin(omega)/omega)/omega**2))
    orders=np.log2(np.array(errors[:-1])/errors[1:]);assert min(orders)>1.9
    return dict(classification='Proven',forcing_moments_checked=True,
        identity='The homogeneous update is the implicit midpoint map. Exact I0 and I1 retain the original heat-pole startup; velocity correction is M^-1*(I0-2*I1/h).',
        numerical_classification='Counterexample candidate',oscillator_errors=errors,orders=orders.tolist(),
        scope='Forcing identities and oscillator control only; actual four histories must pass unchanged relative2% and order1.5.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='95d99987',
        claim='Remove mixed-field reduced eigendecomposition as a possible error source by directly evolving the entire same-space coupled matrix.',
        method='Fixed512/1024/2048 Cayley steps, exact forcing moments, two extended residual corrections per step. Original degrees4,2,1, source, horizon and65 readouts. One sparse real factorization reused per path.',
        gates=dict(time_relative=.02,time_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Measure128 steps on degree4. Execute only if forecast fits300s. Degree4 first; stop on any time gate failure. Then degrees2,1. Run the three original contrasts only after spatial acceptance; no automatic extra steps, degree or cells.',
        budget=dict(hard_seconds=300,pilot_steps=128,maximum_steps=2048,CPU_threads=1,memory_GB=3,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(space.__file__),Path(task.__file__),task.OUT/'hierarchical/result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());start=time.monotonic();m=space.Model(4,task.BANK/'fine-bank.npz');p=Stepper(m,2048);setup=time.monotonic()-start
    start=time.monotonic()
    for _ in range(128):p.step()
    elapsed=time.monotonic()-start
    # Main degrees4,2,1 weighted conservatively as2 degree4 paths; contrasts3.
    forecast=1.4*(6*setup+(2*3584+3*2048)/128*elapsed)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,steps128_seconds=elapsed,
        forecast_seconds=forecast,linear_residual=p.error,assumptions='Main degrees cost2 times degree4 and3 contrasts cost degree4; six setups plus40 percent allowance. Other degrees and late steps unmeasured.'))
    print('TIME FORECAST',forecast,flush=True);signal.alarm(0)


def run():
    plan=json.loads((OUT/'plan.json').read_text());assert not (OUT/'result.json').exists()
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<300
    signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic();cases={};temporal={};comparisons={}
    time_pass=True;spatial_pass=False
    for degree in [4,2,1]:
        m=space.Model(degree,task.BANK/'fine-bank.npz');rows={n:evolve(m,n,f'p{degree}-{n}') for n in [512,1024,2048]};del m
        cases[str(degree)]=rows[2048];temporal[str(degree)]={}
        for field in FIELDS:
            a,b,c=[np.array([r[field] for r in rows[n]['history']]) for n in [512,1024,2048]]
            norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
            temporal[str(degree)][field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
        time_pass=all(v['last']<.02 and v['order']>1.5 for d in temporal.values() for v in d.values())
        if not time_pass:break
    if len(cases)==3:
        for field in FIELDS:
            a,b,c=[np.array([r[field] for r in cases[str(n)]['history']]) for n in [1,2,4]]
            norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
            comparisons[field]=dict(previous=float(d1),last=float(d2),decreased=bool(d2<d1))
        spatial_pass=all(v['last']<.02 and v['decreased'] for v in comparisons.values())
    if time_pass and spatial_pass:
        for label,bank,outer,quadrature in [('coefficient','coarse-bank.npz',2,6),('outer','fine-bank.npz',3,6),('quadrature','fine-bank.npz',2,8)]:
            m=space.Model(4,task.BANK/bank,outer,quadrature);cases[label]=evolve(m,2048,label);del m
            for field in FIELDS:
                c=np.array([r[field] for r in cases['4']['history']]);d=np.array([r[field] for r in cases[label]['history']])
                comparisons[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
    passed=time_pass and spatial_pass and all(v.get('coefficient',1)<.02 and v.get('outer',1)<.002 and v.get('quadrature',1)<.002 for v in comparisons.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,passed=passed,time_passed=time_pass,
        spatial_passed=spatial_pass,temporal=temporal,comparisons=comparisons,paths=list(cases),seconds=time.monotonic()-start,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('TIME RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
