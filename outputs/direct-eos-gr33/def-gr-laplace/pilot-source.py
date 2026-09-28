"""Exact-source resolvent of the unchanged GR descriptor equation.

Counterexample candidate. Fourier inversion is a numerical approximation;
its cutoff/alias error must pass the original component thresholds.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
import sympy as sp
from scipy.sparse.linalg import splu
import def_radial_conduction_reconnect as prior

task=prior.old.task
OUT=prior.OUT.parent/'def-gr-laplace'
write=prior.write
TIMES=np.arange(65)/64
FIELDS=['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']


def invert(samples,sigma,period,times=TIMES):
    # Unwindowed Bromwich trapezoid. No modal filtering or time averaging.
    k=np.arange(len(samples));phase=np.exp(2j*np.pi*times[:,None]*k/period)
    phase[:,0]*=.5
    return 2/period*np.exp(sigma*times)[:,None]*(phase@samples).real


def control():
    t,s,rate=sp.symbols('t s rate',positive=True)
    energy=t+(sp.exp(-rate*t)-1)/rate
    assert sp.simplify(sp.laplace_transform(energy,t,s,noconds=True)-rate/(s*s*(s+rate)))==0
    errors=[]
    for cutoff in [128,256,512]:
        z=6+2j*np.pi*np.arange(cutoff+1)/4
        transform=1/(z*z*(z*z+9))
        y=invert(transform[:,None],6,4)[:,0]
        v=invert((z*transform)[:,None],6,4)[:,0]
        exact=(TIMES-np.sin(3*TIMES)/3)/9
        exact_v=(1-np.cos(3*TIMES))/9
        errors.append(float(max(np.max(abs(y[1:]-exact[1:])),np.max(abs(v[1:]-exact_v[1:])))))
    assert errors[-1]<1e-6 and errors[-1]<errors[0]/8
    return dict(classification='Proven',passed=True,energy_transform='tc*amplitude*lambda/(s^2*(s+lambda))',
        GR_transform='(K-s^2 D) u_hat = F_hat + K H_hat, u(0)=u_prime(0)=0; y=u-H',
        numerical_classification='Counterexample candidate',oscillator_errors=errors,
        scope='The transform identity is exact. Finite Fourier cutoff and periodic aliasing are not proved small by this identity.')


class Problem:
    def __init__(self,bank,outer=2):
        ray=dict(np.load(task.coupled.OUT/'fine-rays.npz'))
        self.radiation=task.Radiation(ray,bank,False);self.heat=self.radiation.heat
        self.bg=task.coupled.Background(self.radiation,outer)
        fn,_=task.coupled.reactive.symbolic()
        self.K,self.D,self.base=task.coupled.assemble(self.bg,fn)
        self.original_faces=self.heat.faces
        self.error=0.;self.outer=outer
        self.native=self.bg.native['radius_cm'][::-1]/(100*self.radiation.geometry.R)
        self.weights=self.bg.native['dm'][::-1];self.weights/=self.weights.sum()
        r=self.bg.grid;N,a=self.radiation.geometry.metric(r)
        self.speed=a/N*r*task.h.gr.C
        faces=self.bg.native['faces_cm']/(100*self.radiation.geometry.R)
        self.masks=[abs(self.native-faces[i])<=538000/self.radiation.geometry.R for i in [4123,1723]]

    def transform(self,s):
        heat=self.heat;lam=heat.rates*self.radiation.geometry.tc
        energy=np.zeros(len(heat.edges),complex)
        energy[heat.face_ids]=np.sum(heat.amplitude*self.radiation.geometry.tc*lam/(s*s*(s+lam)),axis=1)
        pieces=[]
        try:
            # Reuse the exact real-valued source assembly for both components.
            for component in [energy.real,energy.imag]:
                heat.faces=lambda t:(np.zeros(len(component)),component)
                pieces.append(self.base(0.)+self.K@heat.lift(0.,self.bg.nodes))
        finally:
            heat.faces=self.original_faces
        rhs=pieces[0]+1j*pieces[1]
        matrix=(self.K-s*s*self.D).tocsc();scale=np.asarray(abs(matrix).sum(1)).ravel()
        lu=splu(matrix.multiply((1/scale)[:,None]).tocsc())
        answer=lu.solve(rhs/scale)
        extended=matrix.astype(np.clongdouble);r=rhs.astype(np.clongdouble)
        for _ in range(2):answer+=lu.solve(np.asarray((r-extended@answer)/scale,complex))
        residual=float(np.max(abs(r-extended@answer)/(abs(extended)@abs(answer)+abs(r)+1e-100)))
        self.error=max(self.error,residual)
        assert residual<1e-9 and np.all(np.isfinite(answer))
        # Store Eulerian scalar before Fourier cancellation. This is an exact
        # output transformation, not a changed spatial operator or scalar ODE.
        answer=answer.reshape(-1,4)
        answer[:,2]-=self.bg.grid*self.bg.nodes['v']*answer[:,0]
        return answer.ravel()

    def series(self,samples,sigma,period,cutoff,label):
        samples=samples[:cutoff+1];s=sigma+2j*np.pi*np.arange(cutoff+1)/period
        u=invert(samples,sigma,period).reshape(65,-1,4)
        v=invert(s[:,None]*samples,sigma,period).reshape(65,-1,4)
        raw_initial=dict(displacement=float(abs(u[0]).max()),velocity=float(abs(v[0]).max()))
        # The initial state is prescribed, not inferred from a truncated Fourier
        # series at its endpoint. All positive-time samples remain untouched.
        u[0]=0;v[0]=0;history=[]
        r=self.bg.grid
        for j,t in enumerate(TIMES):
            lift=self.heat.lift(t,self.bg.nodes).reshape(-1,4)
            lift_v=self.heat.lift(t,self.bg.nodes,True).reshape(-1,4)
            z=u[j,:,0]-lift[:,0];vz=v[j,:,0]-lift_v[:,0]
            cv=np.interp(self.native,r,self.speed*vz);cf=np.interp(self.native,r,u[j,:,2])
            history.append(dict(tau=float(t),velocity_mass_RMS_m_s=float(np.sqrt(self.weights@(cv*cv))),
                scalar_mass_RMS=float(np.sqrt(self.weights@(cf*cf))),
                **{name:float(np.sqrt(self.weights[m]@cv[m]**2/self.weights[m].sum())) for name,m in zip(FIELDS[2:],self.masks)}))
        q=u[-1].copy();p=v[-1].copy()
        q[:,2]+=r*self.bg.nodes['v']*q[:,0];p[:,2]+=r*self.bg.nodes['v']*p[:,0]
        q-=self.heat.lift(1.,self.bg.nodes).reshape(-1,4)
        p-=self.heat.lift(1.,self.bg.nodes,True).reshape(-1,4)
        flux,energy=self.heat.faces(1.)
        balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
        row=dict(classification='Counterexample candidate',cutoff=cutoff,sigma=sigma,period=period,outer=self.outer,
            linear_residual=self.error,heat_telescoping=balance,raw_Fourier_initial=raw_initial,history=history)
        np.savez_compressed(OUT/(label+'.npz'),grid=r,response=q,velocity=p,Eulerian_scalar=u[-1,:,2],
            heat_luminosity_erg_s=flux,heat_cumulative_energy_erg=energy)
        write(OUT/(label+'.json'),row)
        print('PATH',label,history[-1],flush=True)
        return row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(task.__file__),Path(task.coupled.__file__),prior.OUT/'fine-bank.npz',
           prior.OUT/'coarse-bank.npz',prior.OUT/'solver-source.py',prior.OUT/'result.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='f5c11912',
        claim='Propagate the original failed GR response by the exact-source resolvent without changing K,D, the heat source, initial state, grid, duration or outputs.',
        method='Unwindowed Fourier inversion of the Bromwich line, sigma=6, period=4; cutoffs128/256/512. Reuse nested transform samples. Same65 physical sample times. These are frequency quadratures, not renamed16/32/64 time steps.',
        decision='First measure16 actual resolvents. Then evaluate one512-frequency response; stop if any original component fails cutoff relative2% or order1.5. Only if all pass, perform fixed coefficient, outer and sigma8 contrasts. No automatic frequency/cutoff expansion.',
        gates=dict(relative=.02,order=1.5,coefficient=.02,outer=.002,alias=.0002,linear_residual=1e-9,heat_balance=2e-13),
        budget=dict(pilot_resolvents=16,main_resolvents=513,conditional_contrast_resolvents=1539,
                    hard_seconds=240,CPU_threads=1,memory_GB=3,native_EOS_calls=0),
        forecast='To be determined from16 actual resolvents before production. A new numerical representation, not a full nonlinear GR or photon completion claim.',
        bindings={str(p):prior.digest(p) for p in paths}))
    write(OUT/'control.json',control())
    signal.alarm(60);start=time.monotonic();problem=Problem(prior.OUT/'fine-bank.npz')
    setup=time.monotonic()-start;start=time.monotonic()
    samples=np.array([problem.transform(6+2j*np.pi*k/4) for k in range(16)])
    seconds=time.monotonic()-start
    np.save(OUT/'pilot.npy',samples)
    forecast=setup*4+seconds/16*2052*1.6+20
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',setup_seconds=setup,
        resolvent_seconds=seconds,resolvents=16,forecast_total_seconds=forecast,
        linear_residual=problem.error,assumption='Four setups plus2052 solves at1.6 times measured mean, plus20s inversion; high-frequency and exterior costs unmeasured.'))
    print('PILOT',seconds,'FORECAST',forecast,flush=True);signal.alarm(0)


def run():
    assert not (OUT/'result.json').exists()
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert prior.digest(Path(p))==h,p
    pilot=json.loads((OUT/'pilot.json').read_text());assert pilot['forecast_total_seconds']<240
    signal.alarm(240);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)))
    start=time.monotonic();problem=Problem(prior.OUT/'fine-bank.npz');samples=list(np.load(OUT/'pilot.npy'))
    for k in range(len(samples),513):samples.append(problem.transform(6+2j*np.pi*k/4))
    samples=np.array(samples);np.save(OUT/'fine-transform.npy',samples)
    cases={str(n):problem.series(samples,6,4,n,'fine-'+str(n)) for n in [128,256,512]}
    comparisons={}
    for field in FIELDS:
        series=lambda key:np.array([x[field] for x in cases[key]['history']])
        a,b,c=[series(str(n)) for n in [128,256,512]];norm=max(abs(c).max(),1e-100)
        e1=float(max(abs(a-b))/norm);e2=float(max(abs(b-c))/norm)
        comparisons[field]=dict(cutoff_previous=e1,cutoff_last=e2,order=float(np.log2(e1/e2)))
    passed=all(x['cutoff_last']<.02 and x['order']>1.5 for x in comparisons.values())
    if passed:
        for key,bank,outer,sigma in [('coefficient',prior.OUT/'coarse-bank.npz',2,6),('outer',prior.OUT/'fine-bank.npz',3,6),('alias',prior.OUT/'fine-bank.npz',2,8)]:
            p=Problem(bank,outer);z=np.array([p.transform(sigma+2j*np.pi*k/4) for k in range(513)])
            np.save(OUT/(key+'-transform.npy'),z);row=p.series(z,sigma,4,512,key);cases[key]=row
            for field in FIELDS:
                c=np.array([v[field] for v in cases['512']['history']]);d=np.array([v[field] for v in row['history']])
                comparisons[field][key]=float(max(abs(c-d))/max(abs(c).max(),1e-100))
        passed=all(v['coefficient']<.02 and v['outer']<.002 and v['alias']<.0002 for v in comparisons.values())
    passed=passed and max(v['linear_residual'] for v in cases.values())<1e-9 and max(v['heat_telescoping'] for v in cases.values())<2e-13
    result=dict(classification='Counterexample candidate',passed=passed,comparisons=comparisons,
        actual_GR_descriptor_propagated=True,original_same_input_failure_resolved=passed,
        seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        full_temperature_feedback=False,photon_transport=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
