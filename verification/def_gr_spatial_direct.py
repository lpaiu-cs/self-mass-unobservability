"""Unreduced canonical GR resolvents with exact inertial-source subtraction."""
from pathlib import Path
import argparse
import json
import signal
import time
import resource
import numpy as np
from scipy.sparse.linalg import splu
import def_gr_spatial_repair as task

OUT=task.OUT/'direct';write=task.write;digest=task.digest;TIMES=task.modes.TIMES;FIELDS=task.FIELDS


class Problem:
    def __init__(self,degree=4,bank=None,outer=2,nquad=6):
        self.model=task.Model(degree,bank or task.BANK/'fine-bank.npz',outer,nquad);m=self.model
        self.M=m.M;self.K=m.K;self.ms=np.sqrt(self.M.diagonal())
        from scipy.sparse import diags
        D=diags(1/self.ms);self.mass=splu((D@self.M@D).tocsc());self.error=0.
        r=m.original.native;N,a=m.original.radiation.geometry.metric(r)
        self.speed=a/N*r*task.fem.base.task.h.gr.C

    def mass_solve(self,x):return self.mass.solve(x/self.ms)/self.ms

    def transform(self,s):
        m=self.model;heat=m.heat;tc=m.original.radiation.geometry.tc;lam=heat.rates*tc
        E=np.zeros(len(heat.edges),complex);E[heat.face_ids]=np.sum(heat.amplitude*tc*lam/(s*s*(s+lam)),axis=1)
        rhs=m.load@E;A=(self.K+s*s*self.M).tocsc();scale=np.asarray(abs(A).sum(1)).ravel()
        lu=splu(A.multiply((1/scale)[:,None]).tocsc());u=lu.solve(rhs/scale);extended=A.astype(np.clongdouble)
        for _ in range(2):u+=lu.solve(np.asarray((rhs.astype(np.clongdouble)-extended@u.astype(np.clongdouble))/scale,complex))
        error=float(np.max(abs(rhs.astype(np.clongdouble)-extended@u.astype(np.clongdouble))/(abs(extended)@abs(u)+abs(rhs)+1e-100)))
        self.error=max(self.error,error);assert error<1e-9
        inertial=self.mass_solve(rhs.real).astype(complex)+1j*self.mass_solve(rhs.imag)
        inertial/=s*s
        return u,inertial

    def inertial(self,t):
        m=self.model;heat=m.heat;tc=m.original.radiation.geometry.tc;lam=heat.rates*tc
        # Only t=0 or t>=1/64 is evaluated; no cancellation-prone small lambda*t.
        E=np.zeros(len(heat.edges));V=E.copy()
        E[heat.face_ids]=np.sum(heat.amplitude*tc*(t**3/6-t*t/(2*lam)+t/lam**2+np.expm1(-lam*t)/lam**3),axis=1)
        V[heat.face_ids]=np.sum(heat.amplitude*tc*(t*t/2-t/lam-np.expm1(-lam*t)/lam**2),axis=1)
        return self.mass_solve(m.load@E),self.mass_solve(m.load@V)

    def run(self,label,sigma=6,pilot=None):
        began=time.monotonic();m=self.model;nt=len(TIMES);nr=len(m.original.native)
        sums=np.zeros((nt,2,nr),complex);end=np.zeros((2,m.size),complex)
        history={};period=4.;prefactor=2/period*np.exp(sigma*TIMES)
        for k in range(2049):
            s=sigma+2j*np.pi*k/period
            u,b=pilot[k] if pilot is not None and k<len(pilot) else self.transform(s)
            z=u-b;native=np.array([self.speed*(m.nativeV[0]@(s*z)),m.nativeV[1]@z])
            phase=np.exp(2j*np.pi*k*TIMES/period);weight=.5 if k==0 else 1.
            sums+=weight*phase[:,None,None]*native[None,:,:]
            end+=weight*phase[-1]*np.array([z,s*z])
            if k in [512,1024,2048]:
                values=prefactor[:,None,None]*sums.real;rows=[]
                for j,t in enumerate(TIMES):
                    bu,bv=self.inertial(t);flux,energy=m.heat.faces(t)
                    values[j,0]+=self.speed*(m.nativeV[0]@(bv-m.H@(flux*m.original.radiation.geometry.tc)))
                    values[j,1]+=m.nativeV[1]@bu
                    if j==0:values[j]=0.
                    cv,cf=values[j];w=m.original.weights
                    rows.append(dict(tau=float(t),velocity_mass_RMS_m_s=float(np.sqrt(w@(cv*cv))),scalar_mass_RMS=float(np.sqrt(w@(cf*cf))),
                        **{name:float(np.sqrt(w[mask]@cv[mask]**2/w[mask].sum())) for name,mask in zip(FIELDS[2:],m.original.masks)}))
                bu,bv=self.inertial(1.);flux,energy=m.heat.faces(1.)
                q=prefactor[-1]*end[0].real+bu-m.H@energy;qt=prefactor[-1]*end[1].real+bv-m.H@(flux*m.original.radiation.geometry.tc)
                balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
                name=f'{label}-{k}'
                np.savez_compressed(OUT/(name+'.npz'),grid=m.grid,cells=m.cells,indices=m.indices,q=q,qt=qt,
                    native_radius=m.original.native,native_velocity=values[:,0],native_scalar=values[:,1],
                    weights=m.original.weights,masks=m.original.masks,heat_energy=energy,heat_flux=flux)
                row=dict(classification='Counterexample candidate',degree=m.degree,cutoff=k,sigma=sigma,history=rows,
                    linear_residual=self.error,heat_telescoping=balance,seconds=time.monotonic()-began)
                write(OUT/(name+'.json'),row);history[k]=row;print('DIRECT',name,row['seconds'],rows[-1],flush=True)
        return history


def control():
    import sympy as s
    t,l,z=s.symbols('t l z',positive=True)
    b=t**3/6-t*t/(2*l)+t/l**2+(s.exp(-l*t)-1)/l**3
    assert s.simplify(s.diff(b,t,2)-t-(s.exp(-l*t)-1)/l)==0
    assert b.subs(t,0)==0 and s.diff(b,t).subs(t,0)==0
    errors={}
    for omega in [3,50,300,800,2000]:
        vals=[]
        for cutoff in [512,1024,2048]:
            zz=6+2j*np.pi*np.arange(cutoff+1)/4
            transform=1/(zz*zz*(zz*zz+omega**2))-1/zz**4
            values=task.fem.base.invert((zz*transform)[:,None],6,4)[:,0]+TIMES**2/2
            exact=(1-np.cos(omega*TIMES))/omega**2
            vals.append(float(np.max(abs(values[1:]-exact[1:]))/max(abs(exact).max(),1e-100)))
        errors[str(omega)]=vals
    return dict(classification='Proven',identity_checked=True,numerical_classification='Counterexample candidate',
        oscillator_relative_velocity_errors=errors,formula='Subtract M^-1*F_hat/s^2 from the resolvent and add its exact zero-initial-data time integral. This is an inverse-transform identity, not a physical filter.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='95d99987',
        claim='Evaluate the same polynomial GR model without a reduced eigenbasis, preserving the weak scalar component in direct coupled linear solves.',
        method='Bromwich Fourier line sigma6, period4, nested cutoffs512/1024/2048, with exact inertial-source subtraction. Reuse the first16 resolvents. Same degrees1/2/4, source, horizon and readouts.',
        gates=dict(time_relative=.02,time_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,alias=.0002,linear_residual=1e-9),
        decision='Measure each degree on16 solves. Run degree4 first; stop on any time gate failure. Only then complete degrees2,1. Only if spatial gates pass, run coefficient, outer3R, quadrature8 and sigma8 contrasts. No cutoff or degree expansion.',
        budget=dict(hard_seconds=300,maximum_models=7,CPU_threads=1,memory_GB=3,new_EOS_calls=0),
        bindings={str(p):digest(p) for p in [Path(__file__),Path(task.__file__),task.OUT/'balanced/result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());measurements={}
    for degree in [4,2,1]:
        start=time.monotonic();p=Problem(degree);setup=time.monotonic()-start;start=time.monotonic()
        samples=np.array([p.transform(6+2j*np.pi*k/4) for k in range(16)])
        elapsed=time.monotonic()-start;np.save(OUT/f'pilot-p{degree}.npy',samples)
        measurements[str(degree)]=dict(setup_seconds=setup,solve16_seconds=elapsed,linear_residual=p.error)
        del p
    # Five p4-sized models include the main and four conditional contrasts.
    forecast=1.35*sum((5 if degree=='4' else 1)*(v['setup_seconds']+2049/16*v['solve16_seconds']+5) for degree,v in measurements.items())
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',measurements=measurements,
        forecast_seconds=forecast,assumption='Five p4 models plus p2 and p1,16-solve measurements with35 percent allowance and5s per inverse transform. High-frequency, quadrature and exterior costs unmeasured.'))
    print('DIRECT FORECAST',forecast,flush=True);signal.alarm(0)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<300
    signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
    cases={};temporal={};time_pass=True
    for degree in [4,2,1]:
        problem=Problem(degree);rows=problem.run(f'p{degree}',pilot=np.load(OUT/f'pilot-p{degree}.npy'))
        cases[str(degree)]=rows[2048];temporal[str(degree)]={}
        for field in FIELDS:
            a,b,c=[np.array([h[field] for h in rows[n]['history']]) for n in [512,1024,2048]]
            norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
            temporal[str(degree)][field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)))
        del problem
        time_pass=all(v['last']<.02 and v['order']>1.5 for p in temporal.values() for v in p.values())
        if not time_pass:break
    spatial={};spatial_pass=False
    if len(cases)==3:
        for field in FIELDS:
            a,b,c=[np.array([h[field] for h in cases[str(n)]['history']]) for n in [1,2,4]]
            norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
            spatial[field]=dict(previous=float(d1),last=float(d2),decreased=bool(d2<d1))
        spatial_pass=all(v['last']<.02 and v['decreased'] for v in spatial.values())
    if time_pass and spatial_pass:
        for label,bank,outer,nquad,sigma in [('coefficient','coarse-bank.npz',2,6,6),('outer','fine-bank.npz',3,6,6),('quadrature','fine-bank.npz',2,8,6),('alias','fine-bank.npz',2,6,8)]:
            problem=Problem(4,task.BANK/bank,outer,nquad);cases[label]=problem.run(label,sigma)[2048];del problem
            for field in FIELDS:
                c=np.array([h[field] for h in cases['4']['history']]);d=np.array([h[field] for h in cases[label]['history']])
                spatial[field][label]=float(np.max(abs(c-d))/max(abs(c).max(),1e-100))
    passed=time_pass and spatial_pass and all(v.get('coefficient',1)<.02 and v.get('outer',1)<.002 and v.get('quadrature',1)<.002 and v.get('alias',1)<.0002 for v in spatial.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,passed=passed,time_passed=time_pass,
        spatial_passed=spatial_pass,temporal=temporal,spatial=spatial,paths=list(cases),seconds=time.monotonic()-start,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('DIRECT RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
