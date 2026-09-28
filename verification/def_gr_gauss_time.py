"""Fourth-order full-matrix Gauss propagation; original heat moments retained."""
from pathlib import Path
import argparse
import time
import json
import signal
import resource
import numpy as np
from scipy.sparse.linalg import splu
import def_gr_direct_time as previous
import def_gr_inverse_projection as common

task=common.task;space=common.space;write=task.write;OUT=common.OUT/'gauss'
c=np.array([.5-np.sqrt(3)/6,.5+np.sqrt(3)/6]);a=np.array([[.25,.25-np.sqrt(3)/6],[.25+np.sqrt(3)/6,.25]])
ai=np.linalg.inv(a);b=ai@ai


class Stepper(previous.Stepper):
    def __init__(self,model,steps):
        super().__init__(model,steps)
        values,self.vectors=np.linalg.eig(b);self.inverse=np.linalg.inv(self.vectors)
        assert np.allclose(values[1],values[0].conjugate())
        self.matrix=(model.K+values[0]*model.M/self.dt**2).tocsc()
        self.scale=np.asarray(abs(self.matrix).sum(1)).ravel()
        self.lu=splu(self.matrix.multiply((1/self.scale)[:,None]).tocsc());self.extended=self.matrix.astype(np.clongdouble)

    def step(self):
        h=self.dt;t=self.n*h;F=self.force
        forcing=(t+h*c[:,None])*F[:,0]-F[:,1]
        if self.n==0:
            # Exact I0 and I1 of the microscopic pole pulse fix its two stage
            # load coefficients. This is a quadrature rule for the original
            # source, not a new physical source or a clipped startup.
            forcing+=(4-6*c[:,None])*F[:,2]/h+(12*c[:,None]-6)*F[:,3]/h**2
        history=self.q[None,:]+h*c[:,None]*self.v
        rhs=forcing+b@np.array([self.model.M@q for q in history])/h**2
        r=(self.inverse@rhs)[0];z=self.lu.solve(r/self.scale)
        for _ in range(2):
            defect=r.astype(np.clongdouble)-self.extended@z.astype(np.clongdouble)
            z+=self.lu.solve(np.asarray(defect/self.scale,complex))
        defect=r.astype(np.clongdouble)-self.extended@z.astype(np.clongdouble)
        error=float(np.max(abs(defect)/(abs(self.extended)@abs(z)+abs(r)+1e-100)))
        self.error=max(self.error,error);assert self.error<1e-9
        stages=(self.vectors@np.array([z,z.conjugate()])).real
        stagev=ai@(stages-self.q[None,:])/h
        newq=self.q+h*np.mean(stagev,axis=0)
        self.v+=np.array([.5,.5])@ai@(stagev-self.v[None,:]);self.q=newq;self.n+=1


def evolve(model,steps,label):
    previous.Stepper=Stepper;previous.OUT=OUT
    return previous.evolve(model,steps,label)


def control():
    import sympy as s
    z=s.symbols('z');aa=s.Matrix([[s.Rational(1,4),s.Rational(1,4)-s.sqrt(3)/6],
        [s.Rational(1,4)+s.sqrt(3)/6,s.Rational(1,4)]])
    R=1+z*(s.Matrix([[s.Rational(1,2),s.Rational(1,2)]])*(s.eye(2)-z*aa).inv()*s.ones(2,1))[0]
    assert s.simplify(R-(1+z/2+z*z/12)/(1-z/2+z*z/12))==0
    w=s.symbols('w',real=True);assert s.simplify(R.subs(z,s.I*w)*R.subs(z,-s.I*w)-1)==0
    h,I0,I1=s.symbols('h I0 I1');cc=[s.Rational(1,2)-s.sqrt(3)/6,s.Rational(1,2)+s.sqrt(3)/6]
    ff=s.Matrix([I0/h+(6*I0/h-12*I1/h**2)*(v-s.Rational(1,2)) for v in cc]);bb=s.Matrix([[s.Rational(1,2),s.Rational(1,2)]])
    assert s.simplify((h*bb*ff)[0]-I0)==0 and s.simplify((h*h*bb*aa*ff)[0]-I1)==0
    errors=[];omega=17.;lam=1e12
    for steps in [64,128,256]:
        dt=1/steps;y=np.zeros(2);L=np.array([[0.,1.],[-omega**2,0.]])
        matrix=np.eye(4)-dt*np.kron(a,L)
        for n in range(steps):
            f=n*dt+c*dt-1/lam
            if n==0:f+=(4-6*c)/(dt*lam**2)+(12*c-6)/(dt*dt*lam**3)
            stages=np.linalg.solve(matrix,np.tile(y,2)+dt*np.kron(a,np.eye(2))@np.column_stack([np.zeros(2),f]).ravel()).reshape(2,2)
            y+=dt*np.mean(stages@L.T+np.column_stack([np.zeros(2),f]),axis=0)
        S,C,F,co=common.modes.functions(omega**2,1.)
        exact=np.array([F-C/lam+S/(lam*lam+omega**2)-co/(lam*(lam*lam+omega**2)),
            C-S/lam+co/(lam*lam+omega**2)+omega**2*S/(lam*(lam*lam+omega**2))])
        errors.append(float(np.linalg.norm((y-exact)*[omega,1.])))
    orders=np.log2(np.array(errors[:-1])/errors[1:]);assert min(orders)>3.9
    return dict(classification='Proven',stability_and_moment_identities_checked=True,stability_function='(1+z/2+z^2/12)/(1-z/2+z^2/12)',
        numerical_classification='Counterexample candidate',fast_pole_oscillator_errors=errors,orders=orders.tolist(),
        scope='Unit-modulus homogeneous stability and exact first two source moments; high-mode startup quadrature and full GR still require actual unchanged convergence gates.')


def prepare():
    assert not OUT.exists();OUT.mkdir();signal.alarm(90)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1260c415',
        claim='Reduce full-matrix wave phase error from midpoint to fourth-order Gauss without adding steps or mixing separately accepted components.',
        method='Same512/1024/2048 steps, degree4/2/1 spaces and original heat input. Exact first two microscopic forcing moments feed Gauss stage loads; all later linear ramp loads are exact. Full coupled sparse solve, no modal reduction or added damping.',
        gates=dict(time_relative=.02,time_order=1.5,spatial=.02,spatial_decrease=True,coefficient=.02,outer=.002,quadrature=.002,linear_residual=1e-9,heat_balance=2e-13),
        decision='Measure128 steps, then degree4 only if first-case forecast fits180s. Stop on any time failure. Only if all four pass perform fixed degree2/1 paths and conditional coefficient/outer/quadrature contrasts under600s total cap. No step or degree expansion.',
        budget=dict(pilot_cap_seconds=90,first_case_cap_seconds=180,total_cap_seconds=600,CPU_threads=1,memory_GB=3,new_EOS_calls=0),
        bindings={str(p):task.digest(p) for p in [Path(__file__),Path(previous.__file__),Path(space.__file__),Path(task.__file__),common.OUT/'fields/stage-result.json',task.BANK/'fine-bank.npz',task.BANK/'coarse-bank.npz']}))
    write(OUT/'control.json',control());start=time.monotonic();model=space.Model(4,task.BANK/'fine-bank.npz');p=Stepper(model,2048);setup=time.monotonic()-start
    start=time.monotonic()
    for _ in range(128):p.step()
    elapsed=time.monotonic()-start;forecast=1.4*(setup+3584/128*elapsed+4)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',setup_seconds=setup,steps128_seconds=elapsed,
        first_case_forecast_seconds=forecast,linear_residual=p.error,
        assumption='Same three paths,128 measured steps scaled to3584 with40 percent allowance. Remaining cases conditional and unmeasured. Complex direct factorization costs more than midpoint but changes no grid or step counts.'))
    signal.alarm(0);print('GAUSS FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'stage-result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert task.digest(Path(p))==h,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['first_case_forecast_seconds']<180
    signal.alarm(180);resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)));start=time.monotonic()
    model=space.Model(4,task.BANK/'fine-bank.npz');rows=[evolve(model,n,f'p4-{n}') for n in [512,1024,2048]]
    cmp=common.comparison(rows);passed=all(v['last']<.02 and v['order']>1.5 for v in cmp.values())
    result=dict(classification='Counterexample candidate',actual_same_input_evolved=True,propagation_passed=passed,
        comparisons=cmp,seconds=time.monotonic()-start,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_failure_resolved=False,full_dynamic_charge_solved=False)
    write(OUT/'stage-result.json',result);signal.alarm(0);print('GAUSS RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
