"""Refine original Gauss stages, including stage transform and state updates."""
from pathlib import Path
import argparse
import json
import signal
import time
import resource
import numpy as np
import def_gr_gauss_time as gauss

common=gauss.common;task=gauss.task;space=gauss.space;write=task.write
OUT=task.OUT.parent/'def-gr-gauss-refined';ld=np.longdouble
c=np.array([ld('.5')-np.sqrt(ld(3))/6,ld('.5')+np.sqrt(ld(3))/6])
ai=np.array([[3,-3+2*np.sqrt(ld(3))],[-3-2*np.sqrt(ld(3)),3]],dtype=ld);B=ai@ai


class Stepper(gauss.Stepper):
    def __init__(self,model,steps):
        super().__init__(model,steps)
        self.M=model.M.astype(ld);self.K=model.K.astype(ld)
        self.q=self.q.astype(ld);self.v=self.v.astype(ld);self.force=self.force.astype(ld)
        self.h=ld(1)/steps

    def solve(self,rhs):
        r=(self.inverse@rhs)[0];z=self.lu.solve(np.asarray(r/self.scale,complex))
        return (self.vectors@np.array([z,z.conjugate()])).real.astype(ld)

    def step(self):
        h=self.h;t=self.n*h;F=self.force
        forcing=(t+h*c[:,None])*F[:,0]-F[:,1]
        if self.n==0:forcing+=(4-6*c[:,None])*F[:,2]/h+(12*c[:,None]-6)*F[:,3]/h**2
        history=self.q[None,:]+h*c[:,None]*self.v
        rhs=forcing+B@np.array([self.M@q for q in history])/h**2
        def product(z):return np.array([self.K@q for q in z])+B@np.array([self.M@q for q in z])/h**2
        stages=self.solve(rhs)
        for _ in range(3):stages+=self.solve(rhs-product(stages))
        defect=rhs-product(stages)
        denominator=np.array([abs(self.K)@abs(q) for q in stages])+abs(B)@np.array([abs(self.M)@abs(q) for q in stages])/h**2+abs(rhs)+1e-100
        self.error=max(self.error,float(np.max(abs(defect)/denominator)));assert self.error<1e-9
        stagev=ai@(stages-self.q[None,:])/h
        self.q+=h*np.mean(stagev,axis=0)
        self.v+=np.mean(ai@(stagev-self.v[None,:]),axis=0);self.n+=1


def evolve(model,n,label):
    gauss.previous.Stepper=Stepper;gauss.previous.OUT=OUT
    return gauss.previous.evolve(model,n,label)


def control():
    from scipy.sparse import csc_matrix
    from types import SimpleNamespace
    from scipy.linalg import eig
    K=np.array([[9.,1e-11],[1e-11,25.]]);scale=np.array([1.,2e-12])
    heat=SimpleNamespace(rates=np.array([[1e12]]),amplitude=np.array([[1.]]),edges=np.array([0.]),face_ids=np.array([0]))
    model=SimpleNamespace(M=csc_matrix(np.eye(2)),K=csc_matrix(K),size=2,heat=heat,load=csc_matrix(scale[:,None]),original=SimpleNamespace(radiation=SimpleNamespace(geometry=SimpleNamespace(tc=1.))))
    k,V=eig(K*scale[None,:]/scale[:,None]);k=k.real;V=V.real;b=np.linalg.solve(V,np.ones(2));q=[];v=[]
    for ki in k:
        S,C,T,co=common.modes.functions(ki,1.);lam=1e12
        q.append(T-C/lam+S/(lam*lam+ki)-co/(lam*(lam*lam+ki)))
        v.append(C-S/lam+co/(lam*lam+ki)+ki*S/(lam*(lam*lam+ki)))
    exactq=scale*(V@(b*q));exactv=scale*(V@(b*v));errors=[]
    for n in [16,32,64]:
        step=Stepper(model,n)
        for _ in range(n):step.step()
        errors.append(np.maximum(abs(step.q-exactq),abs(step.v-exactv))/scale)
    errors=np.array(errors,float);orders=np.log2(errors[:-1]/errors[1:]);assert orders.min()>3.8
    return dict(**gauss.control(),actual_refined_weak_field_errors=errors.tolist(),actual_refined_orders=orders.tolist())


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='caf5c5b',
        claim='Isolate rounding outside the transformed complex solve before attributing Gauss failure to unresolved physical waves.',
        evidence='Saved Gauss weak scalar and new-boundary differences are about1e-9 and nonmonotone; SDIRK extended original stages recover2nd order for those fields but not global/old velocity. Refine original full Gauss stages and retain extended accepted states.',
        method='Identical Gauss tableau, same512/1024/2048 steps and original source moments. Complex LU only preconditions the extended real original-stage residual; three corrections, no physical changes, clipping or damping.',
        gates=dict(time_relative=.02,time_order=1.5,linear_residual=1e-9,heat_balance=2e-13,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002),
        decision='Pilot64 steps. Only if first triple forecast fits360s run degree4 triple; stop any failed component. Fixed degree2/1 and original contrasts remain conditional under900s total. No step,degree,period or source expansion.',
        budget=dict(pilot_seconds=90,first_case_seconds=360,total_seconds=900,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(gauss.__file__),Path(gauss.previous.__file__),Path(space.__file__),Path(task.__file__),task.OUT.parent/'def-gr-stiff-time/stage-result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');step=Stepper(model,2048);setup=time.monotonic()-start
    start=time.monotonic()
    for _ in range(64):step.step()
    elapsed=time.monotonic()-start;forecast=1.4*(setup+3584/64*elapsed+6)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,steps64_seconds=elapsed,first_case_forecast_seconds=forecast,linear_residual=step.error,
        assumption='Scale64 measured early steps to3584 with40 percent margin,one setup and6s readout allowance; other spaces and late steps unmeasured.'))
    signal.alarm(0);print('REFINED FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<360
    signal.alarm(360);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');rows=[evolve(model,n,f'p4-{n}') for n in [512,1024,2048]]
    cmp=common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    contrasts={}
    for n,row in zip([512,1024,2048],rows):
        saved=json.loads((common.OUT/'gauss'/f'p4-{n}.json').read_text());contrasts[str(n)]={}
        for key in task.FIELDS:
            a=np.array([r[key] for r in saved['history']]);b=np.array([r[key] for r in row['history']])
            contrasts[str(n)][key]=float(max(abs(a-b))/max(abs(b).max(),1e-100))
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,comparisons=cmp,original_arithmetic_differences=contrasts,
        seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('REFINED RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
