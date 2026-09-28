"""Separate stage-solve roundoff from unresolved GR temporal dynamics.

Same saved transport input, K/D equations, initial state, grid and horizon.
Counterexample candidate; no opacity or whole-star completion claim.
"""
from pathlib import Path
from functools import lru_cache
import argparse
import json
import resource
import signal
import time
import numpy as np
import sympy as sp
from scipy.sparse.linalg import splu
import def_radial_conduction_reconnect as prior

OUT=prior.OUT.parent/'def-gr-transient-radau';write=prior.write;digest=prior.digest


@lru_cache(maxsize=2)
def tableau(stages):
    x=sp.symbols('x')
    nodes=([sp.Rational(1,3),sp.Integer(1)] if stages==2 else [(4-sp.sqrt(6))/10,(4+sp.sqrt(6))/10,sp.Integer(1)])
    columns=[sp.prod((x-nodes[k])/(nodes[j]-nodes[k]) for k in range(stages) if k!=j) for j in range(stages)]
    A=sp.Matrix([[sp.integrate(q,(x,0,c)) for q in columns] for c in nodes]);inverse=sp.simplify(A.inv())
    convert=lambda M:np.array([[str(sp.N(v,40)) for v in row] for row in M.tolist()],dtype=np.longdouble)
    return convert(inverse),np.array([str(sp.N(c,40)) for c in nodes],dtype=np.longdouble),A


class Step:
    stages=3
    def __init__(self,K,D,dt):
        self.K=K.astype(np.longdouble);self.D=D.astype(np.longdouble);self.dt=np.longdouble(dt)
        self.AI,self.c,_=tableau(self.stages);self.B=self.AI@self.AI
        values,V=np.linalg.eig(np.asarray(self.B,float));self.V=V;self.inverse=np.linalg.inv(V)
        self.solvers=[]
        for value in values:
            matrix=(K.astype(complex)-value*D/dt**2).tocsc();scale=np.asarray(abs(matrix).sum(1)).ravel()
            self.solvers.append((splu(matrix.multiply((1/scale)[:,None]).tocsc()),scale))
        self.error=0.;self.correction=0.

    def advance(self,y,v,t,forcing):
        dt=self.dt;history=np.asarray(y,np.longdouble)[None,:]+dt*self.c[:,None]*np.asarray(v,np.longdouble)
        rhs=np.array([forcing(float(t+dt*c)) for c in self.c],np.longdouble)-self.B@np.array([self.D@q for q in history])/dt**2
        def solve(r):
            transformed=self.inverse@r
            z=np.array([lu.solve(np.asarray(q/scale,complex)) for q,(lu,scale) in zip(transformed,self.solvers)])
            return np.asarray((self.V@z).real,np.longdouble)
        def product(z):return np.array([self.K@q for q in z])-self.B@np.array([self.D@q for q in z])/dt**2
        stages=solve(rhs)
        for _ in range(3):
            correction=solve(rhs-product(stages));stages+=correction
        residual=rhs-product(stages)
        scale=np.array([abs(self.K)@abs(q) for q in stages])+abs(self.B)@np.array([abs(self.D)@abs(q) for q in stages])/dt**2+abs(rhs)+1e-100
        self.error=max(self.error,float(np.max(abs(residual)/scale)))
        self.correction=max(self.correction,float(np.linalg.norm(correction)/max(np.linalg.norm(stages),1e-100)))
        assert self.error<1e-9
        return np.asarray(stages[-1],float),np.asarray(self.AI[-1]@(stages-y[None,:])/dt,float)


class PrecisionStep(Step):
    stages=2


def solve(steps,bank,outer,label,precision=False):
    ns=dict(vars(prior.old),OUT=OUT,Step=PrecisionStep if precision else Step)
    source=(prior.OUT/'solver-source.py').read_text();exec(compile(source,str(prior.OUT/'solver-source.py'),'exec'),ns)
    return ns['solve'](steps,bank,False,outer,label)


def control():
    from scipy.sparse import csc_matrix
    K=csc_matrix([[-9.,0.],[-1.,1.]]);D=csc_matrix([[1.,0.],[0.,0.]])
    errors=[];exact=(1-np.sin(3)/3)/9;exact_v=(1-np.cos(3))/9
    for n in [4,8,16]:
        step=Step(K,D,1/n);y=np.zeros(2);v=y.copy()
        for j in range(n):y,v=step.advance(y,v,j/n,lambda t:np.array([-t,0.]))
        assert abs(y[1]-y[0])<1e-13
        errors.append(float(np.hypot(3*(y[0]-exact),v[0]-exact_v)))
    orders=np.log2(np.array(errors[:-1])/errors[1:]);assert orders.min()>4.7
    _,_,A=tableau(3);z=sp.symbols('z')
    R=sp.simplify(1+z*(A[-1,:]*(sp.eye(3)-z*A).inv()*sp.ones(3,1))[0])
    expected=(1+2*z/5+z*z/20)/(1-3*z/5+3*z*z/20-z**3/60)
    assert sp.simplify(R-expected)==0 and sp.limit(R,z,sp.oo)==0
    return dict(classification='Proven',symbolic_passed=True,stability_function=str(expected),
        numerical_control_classification='Counterexample candidate',errors=errors,orders=orders.tolist(),
        scope='Fifth-order L-stable three-stage Radau on a known constrained oscillator; no automatic GR accuracy claim.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(prior.__file__),Path(prior.old.__file__),prior.OUT/'solver-source.py',prior.OUT/'fine-bank.npz',prior.OUT/'coarse-bank.npz',prior.OUT/'result.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='2d1f5299',
        claim='Separate floating-point stage-solve error from low-order temporal resolution, then test one higher-order same-grid method on the actual failed GR component.',
        decision='First run one same64-step two-stage path using extended-precision original-stage residuals and iterative refinement. If readout differences from saved64 exceed1e-4 of their own history norms, stop to repair arithmetic. Otherwise run only the three-stage fifth-order Radau candidate on the original16/32/64 and two input/boundary contrasts.',
        retained='Same frozen-source K,D,material/scalar/metric,source,initial state,5864 interior grid nodes,period and acceptance criteria. No smoothing, cutoff removal, extra time steps, changing readouts or relaxed gates.',
        budget=dict(precision_paths=1,candidate_paths=5,hard_seconds=180,cpu_threads=1,memory_GB=3,new_native_EOS_calls=0,new_collision_states=0),
        forecast='Phase71 five paths25.48s. Additional64 precision path and three-stage complex factorizations/stage residuals forecast50-130s, hard180s; changed method cost unmeasured. Stop on arithmetic branch or numerical failure; no next method/refinement automatically.',
        gates=dict(time_relative=.02,time_order=1.5,coefficient_relative=.02,outer_relative=.002,linear_residual=1e-9,heat_balance=2e-13),
        bindings={str(p):digest(p) for p in paths}))
    write(OUT/'control.json',control())


def run():
    assert not (OUT/'result.json').exists();signal.alarm(180);begin=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)))
    for p,h in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert digest(Path(p))==h,p
    row=solve(64,prior.OUT/'fine-bank.npz',2,'precision-64',True)
    original=json.loads((prior.OUT/'heat-64.json').read_text());arithmetic={}
    fields=list(json.loads((prior.OUT/'result.json').read_text())['comparisons'])
    for field in fields:
        a=np.array([r[field] for r in original['history']]);b=np.array([r[field] for r in row['history']])
        arithmetic[field]=float(max(abs(a-b))/max(abs(a).max(),1e-100))
    write(OUT/'precision.json',dict(classification='Counterexample candidate',differences=arithmetic,seconds=time.monotonic()-begin))
    assert max(arithmetic.values())<1e-4,'Arithmetic changes the failed response; stop higher-order branch'
    cases={}
    for n in [16,32,64]:cases[str(n)]=solve(n,prior.OUT/'fine-bank.npz',2,'heat-'+str(n))
    cases['coarse']=solve(64,prior.OUT/'coarse-bank.npz',2,'coefficient-64')
    cases['outer']=solve(64,prior.OUT/'fine-bank.npz',3,'outer-64')
    comparisons={}
    for field in fields:
        series=lambda name:np.array([x[field] for x in cases[name]['history']])
        a,b,c=series('16'),series('32'),series('64');norm=max(abs(c).max(),1e-100)
        e1=float(max(abs(a-b[::2]))/norm);e2=float(max(abs(b-c[::2]))/norm)
        comparisons[field]=dict(time_previous=e1,time_last=e2,order=float(np.log2(e1/e2)),
            coefficients=float(max(abs(c-series('coarse')))/norm),outer=float(max(abs(c-series('outer')))/norm))
    passed=all(x['time_last']<.02 and x['order']>1.5 and x['coefficients']<.02 and x['outer']<.002 for x in comparisons.values())
    balance=max(r['heat_telescoping'] for r in cases.values());residual=max(r['linear_residual'] for r in cases.values())
    passed=passed and balance<2e-13 and residual<1e-9
    result=dict(classification='Counterexample candidate',actual_GR_fluid_scalar_metric_evolved=True,passed=passed,
        comparisons=comparisons,arithmetic_differences=arithmetic,seconds=time.monotonic()-begin,
        max_heat_balance=balance,max_linear_residual=residual,
        memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_same_input_numerical_failure_resolved=passed,original_complete_input_failure_resolved=False,
        full_temperature_feedback=False,photon_transport=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('TRANSIENT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
