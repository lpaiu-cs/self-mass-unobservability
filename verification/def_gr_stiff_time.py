"""SDIRK2 of the full unchanged GR matrix, with original heat moments.

Numerical damping is tested by the original four-component refinement gates;
no physical viscosity, source smoothing, modal cutoff or changed readout.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu
import def_gr_direct_time as previous
import def_gr_inverse_projection as common

task=common.task;space=common.space;OUT=task.OUT.parent/'def-gr-stiff-time';write=task.write
ld=np.longdouble;gamma=1-1/np.sqrt(ld(2))
A=np.array([[gamma,0],[1-gamma,gamma]],dtype=ld);weights=A[-1];nodes=A.sum(1)
# Two original forcing integrals: integral f and integral (h-s)f.
moment_matrix=np.array([weights,weights@A]);det=np.linalg.det(np.asarray(moment_matrix,float))
moment_inverse=np.array([[moment_matrix[1,1],-moment_matrix[0,1]],[-moment_matrix[1,0],moment_matrix[0,0]]])
moment_inverse/=moment_matrix[0,0]*moment_matrix[1,1]-moment_matrix[0,1]*moment_matrix[1,0]


class Stepper(previous.Stepper):
    def __init__(self,model,steps):
        super().__init__(model,steps)
        self.h=ld(1)/steps;self.alpha=gamma*self.h
        self.mass=model.M.astype(ld)
        # Residuals use the original extended stage equation; the double LU is
        # a preconditioner. Store the accepted state in extended precision too.
        self.extended=self.mass+self.alpha**2*model.K.astype(ld)
        self.scale=np.sqrt(np.asarray(self.extended.diagonal(),float));D=diags(1/self.scale)
        self.lu=splu((D@self.extended.astype(float)@D).tocsc())
        self.q=self.q.astype(ld);self.v=self.v.astype(ld);self.force=self.force.astype(ld)

    def implicit(self,q,v,forcing):
        alpha=self.alpha;rhs=self.mass@(q+alpha*v)+alpha**2*forcing
        z=(self.lu.solve(np.asarray(rhs/self.scale,float))/self.scale).astype(ld)
        for _ in range(2):
            defect=rhs-self.extended@z
            z+=(self.lu.solve(np.asarray(defect/self.scale,float))/self.scale).astype(ld)
        defect=rhs-self.extended@z
        error=float(np.max(abs(defect)/(abs(self.extended)@abs(z)+abs(rhs)+1e-100)))
        self.error=max(self.error,error);assert self.error<1e-9
        return z,(z-q)/alpha

    def step(self):
        h=self.h;t=self.n*h;F=self.force
        forces=(t+h*nodes[:,None])*F[:,0]-F[:,1]
        if self.n==0:
            forces+=moment_inverse@np.array([F[:,2]/h,F[:,2]/h-F[:,3]/h**2])
        q1,v1=self.implicit(self.q,self.v,forces[0])
        qbase=self.q+(1-gamma)*h*v1
        vbase=self.v+(1-gamma)/gamma*(v1-self.v)
        self.q,self.v=self.implicit(qbase,vbase,forces[1]);self.n+=1


def evolve(model,steps,label):
    previous.Stepper=Stepper;previous.OUT=OUT
    return previous.evolve(model,steps,label)


def control():
    import sympy as s
    from scipy.sparse import csc_matrix
    from types import SimpleNamespace
    g=1-1/s.sqrt(2);aa=s.Matrix([[g,0],[1-g,g]]);bb=aa[-1,:]
    z,y=s.symbols('z y',real=True)
    R=s.simplify(1+z*(bb*(s.eye(2)-z*aa).inv()*s.ones(2,1))[0])
    assert s.simplify(R-(1+(1-2*g)*z)/(1-g*z)**2)==0
    assert s.simplify((bb*aa*s.ones(2,1))[0]-s.Rational(1,2))==0
    assert s.limit(R,z,s.oo)==0
    assert s.simplify(abs((1-g*s.I*y)**2)**2-abs(1+(1-2*g)*s.I*y)**2-g**4*y**4)==0
    h,I0,I1=s.symbols('h I0 I1');mom=bb.col_join(bb*aa)
    ff=mom.inv()*s.Matrix([I0/h,I1/h**2])
    assert s.simplify((h*bb*ff)[0]-I0)==0 and s.simplify((h*h*bb*aa*ff)[0]-I1)==0
    # Test the actual implementation with a weak coupled field, and a heat pole
    # too fast for the time mesh. All fields retain both coupling directions.
    M=csc_matrix(np.eye(2));K=csc_matrix([[9.,1e-11],[1e-11,25.]])
    heat=SimpleNamespace(rates=np.array([[1e12]]),amplitude=np.array([[1.]]),edges=np.array([0.]),face_ids=np.array([0]))
    model=SimpleNamespace(M=M,K=K,size=2,heat=heat,load=csc_matrix([[1.],[2e-12]]),original=SimpleNamespace(radiation=SimpleNamespace(geometry=SimpleNamespace(tc=1.))))
    from scipy.linalg import eig
    scale=np.array([1.,2e-12]);k,V=eig(K.toarray()*scale[None,:]/scale[:,None]);k=k.real;V=V.real
    b=np.linalg.solve(V,np.ones(2));q=[];v=[]
    for ki in k:
        S,C,T,co=common.modes.functions(ki,1.);lam=1e12
        q.append(T-C/lam+S/(lam*lam+ki)-co/(lam*(lam*lam+ki)))
        v.append(C-S/lam+co/(lam*lam+ki)+ki*S/(lam*(lam*lam+ki)))
    exactq=scale*(V@(b*q));exactv=scale*(V@(b*v));errors=[]
    for n in [64,128,256]:
        step=Stepper(model,n)
        for _ in range(n):step.step()
        errors.append(np.maximum(abs(step.q-exactq),abs(step.v-exactv))/scale)
    errors=np.array(errors,float);orders=np.log2(errors[:-1]/errors[1:]);assert orders.min()>1.8
    return dict(classification='Proven',tableau_order_and_stiff_limit_checked=True,source_moments_checked=True,
        numerical_classification='Counterexample candidate',weak_coupled_oscillator_errors=errors.tolist(),orders=orders.tolist(),
        scope='The same second-order L-stable tableau already used in the repository, with full coupled weak-field control; actual GR and spatial acceptance remain separate.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='caf5c5b',
        claim='Resolve unresolved wave remnants and state-update roundoff together using full-matrix L-stable second-order propagation, without changing the physical generator or input.',
        retained='Same original input, zero initial state,512/1024/2048 steps,65 output times,degree4/2/1 spaces,readout definitions and unchanged gates. Exact two forcing moments for microscopic startup. No new EOS calls.',
        method='SDIRK2 reused from repository tableau; one real factorization per path, original extended stage residuals and extended state storage. Numerical damping tends to zero with dt and must pass all four refinement tests and comparison with saved Gauss.',
        gates=dict(time_relative=.02,time_order=1.5,spatial_relative=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13,saved_Gauss_difference=.02),
        decision='Measure128 steps, forecast first degree4 triple; only run within300s. Stop if any propagation gate fails. If passed, fixed degrees2/1 triples and conditional coefficient/outer/quadrature single paths under900s total. No step,degree,basis or period expansion.',
        budget=dict(pilot_seconds=90,first_case_seconds=300,total_seconds=900,CPU_threads=1,memory_GB=4,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(previous.__file__),Path(space.__file__),Path(task.__file__),task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz',common.OUT/'gauss/stage-result.json']}))
    write(OUT/'control.json',control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');p=Stepper(model,2048);setup=time.monotonic()-start
    start=time.monotonic()
    for _ in range(128):p.step()
    elapsed=time.monotonic()-start;forecast=1.4*(setup+3584/128*elapsed+6)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,steps128_seconds=elapsed,
        first_case_forecast_seconds=forecast,linear_residual=p.error,dofs=model.size,
        assumptions='Measured128 early steps scaled to3584,one model setup plus6s readout allowance and40 percent margin; other degrees and late steps unmeasured. Remaining paths conditional.'))
    signal.alarm(0);print('STIFF FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<300
    signal.alarm(300);resource.setrlimit(resource.RLIMIT_AS,(int(4e9),int(4e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');rows=[evolve(model,n,f'p4-{n}') for n in [512,1024,2048]]
    cmp=common.comparison(rows);gauss=json.loads((common.OUT/'gauss/p4-2048.json').read_text());contrast={}
    for key in task.FIELDS:
        a=np.array([r[key] for r in gauss['history']]);b=np.array([r[key] for r in rows[-1]['history']])
        contrast[key]=float(max(abs(a-b))/max(abs(b).max(),1e-100))
    passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values()) and max(contrast.values())<.02
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
        comparisons=cmp,saved_Gauss_differences=contrast,seconds=time.monotonic()-start,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('STIFF RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
