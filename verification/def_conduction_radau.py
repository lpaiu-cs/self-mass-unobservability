"""Same coupled equations with L-stable two-stage Radau IIA time stepping.

Preserves the failed Newmark heat-component test. No spatial/time refinement
or changed acceptance criterion; one conjugate sparse factorization per path.
"""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import splu
import def_conductive_cauchy as task

OUT=task.OUT/'radau'
A=np.array([[5.,-1.],[9.,3.]])/12
AI=np.linalg.inv(A);B=AI@AI


class Step:
    def __init__(self,K,D,dt):
        self.K=K;self.D=D;self.dt=dt;values,V=np.linalg.eig(B)
        assert np.allclose(values[1],values[0].conjugate()) and np.allclose(V[:,1],V[:,0].conjugate())
        self.V=V;self.inverse=np.linalg.inv(V);self.matrix=K-values[0]*D/dt**2
        self.scale=np.asarray(abs(self.matrix).sum(1)).ravel();self.lu=splu(self.matrix.multiply((1/self.scale)[:,None]).tocsc())
        self.extended=self.matrix.astype(np.clongdouble);self.error=0.

    def advance(self,y,v,t,forcing):
        dt=self.dt;history=y[None,:]+dt*np.array([1/3,1])[:,None]*v
        rhs=np.array([forcing(t+dt/3),forcing(t+dt)])-B@np.array([self.D@q for q in history])/dt**2
        r=(self.inverse@rhs)[0];z=self.lu.solve(r/self.scale)
        for _ in range(2):
            defect=r.astype(np.clongdouble)-self.extended@z.astype(np.clongdouble)
            z+=self.lu.solve(np.asarray(defect/self.scale,complex))
        error=np.max(abs(self.matrix@z-r)/(abs(self.matrix)@abs(z)+abs(r)+1e-100));self.error=max(self.error,float(error))
        stages=(self.V@np.array([z,z.conjugate()])).real
        assert self.error<1e-9 and np.all(np.isfinite(stages))
        return stages[1],AI[1]@(stages-y[None,:])/dt


def control():
    errors=[];position_errors=[];discrete_errors=[];omega=3.
    K=csc_matrix([[-omega**2,0.],[-1.,1.]]);D=csc_matrix([[1.,0.],[0.,0.]])
    for n in [16,32,64]:
        step=Step(K,D,1/n);y=np.zeros(2);v=y.copy()
        for j in range(n):y,v=step.advance(y,v,j/n,lambda t:np.array([-t,0.]))
        exact=(1-np.sin(omega)/omega)/omega**2;exact_v=(1-np.cos(omega))/omega**2
        position_errors.append(float(abs(y[0]/exact-1)))
        error=np.hypot(omega*(y[0]-exact),v[0]-exact_v)/np.hypot(omega*exact,exact_v)
        z=-1j*omega/n;R=(1+z/3)/(1-2*z/3+z*z/6);W=R**n*(-1j/omega**2)
        discrete_errors.append(float(abs(omega*(y[0]-1/omega**2)+1j*(v[0]-1/omega**2)-W)))
        assert abs(y[1]-y[0])<1e-14
        errors.append(float(error))
    orders=np.log2(np.array(errors[:-1])/errors[1:]);assert np.min(orders)>2.8 and max(discrete_errors)<1e-12
    # Stability function of two-stage Radau IIA, derived directly from A.
    z=task.sp.symbols('z');aa=task.sp.Matrix([[task.sp.Rational(5,12),-task.sp.Rational(1,12)],[task.sp.Rational(3,4),task.sp.Rational(1,4)]])
    R=1+z*(aa[-1,:]*(task.sp.eye(2)-z*aa).inv()*task.sp.ones(2,1))[0]
    assert task.sp.simplify(R-(1+z/3)/(1-2*z/3+z*z/6))==0 and task.sp.limit(R,z,task.sp.oo)==0
    return dict(classification='Counterexample candidate',passed=True,constrained_oscillator_errors=errors,orders=orders.tolist(),position_errors=position_errors,exact_RK_discrete_errors=discrete_errors,
        stability_function='(1+z/3)/(1-2*z/3+z^2/6); tends to zero for large negative z. Newmark average acceleration has unit modulus on the imaginary axis.',
        limitation='L-stability does not prove accuracy of unresolved physical modes; the original heat-component refinement criterion is still applied.')


def solve(steps,bank,neutrino=False,outer=2,label=None):
    start=time.monotonic();ray=dict(np.load(task.coupled.OUT/'fine-rays.npz'));radiation=task.Radiation(ray,bank,neutrino)
    bg=task.coupled.Background(radiation,outer);fn,_=task.coupled.reactive.symbolic();K,D,base=task.coupled.assemble(bg,fn)
    dt=1/steps;step=Step(K,D,dt);y=np.zeros(K.shape[0]);v=y.copy();history=[]
    native=bg.native['radius_cm'][::-1]/(100*radiation.geometry.R);weights=bg.native['dm'][::-1];weights/=weights.sum()
    r=bg.grid;N,a=radiation.geometry.metric(r)
    def output(t):
        q=(y-radiation.heat.lift(t,bg.nodes)).reshape(-1,4);p=(v-radiation.heat.lift(t,bg.nodes,True)).reshape(-1,4)
        speed=np.interp(native,r,a/N*r*p[:,0]*task.h.gr.C);scalar=np.interp(native,r,q[:,2]-r*q[:,0]*bg.nodes['v'])
        return dict(tau=t,velocity_mass_RMS_m_s=float(np.sqrt(weights@(speed*speed))),scalar_mass_RMS=float(np.sqrt(weights@(scalar*scalar)))),q,p
    forcing=lambda t:base(t)+K@radiation.heat.lift(t,bg.nodes)
    history.append(output(0.)[0]);began=time.monotonic()
    for j in range(steps):
        y,v=step.advance(y,v,j*dt,forcing);history.append(output((j+1)*dt)[0])
    step_seconds=time.monotonic()-began;_,q,p=output(1.);flux,energy=radiation.heat.faces(1.)
    balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    row=dict(classification='Counterexample candidate',steps=steps,neutrino=neutrino,outer=outer,linear_residual=step.error,
        seconds=time.monotonic()-start,step_seconds=step_seconds,heat_telescoping=balance,history=history)
    if label:
        np.savez_compressed(OUT/(label+'.npz'),grid=bg.grid,response=q,velocity=p,heat_luminosity_erg_s=flux,heat_cumulative_energy_erg=energy)
        task.core.ex.write(OUT/(label+'.json'),row)
    print('RADAU',label,steps,'SECONDS',row['seconds'],'END',history[-1],flush=True)
    return row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(task.__file__),task.OUT/'plan.json',task.OUT/'component/result.json',task.core.OUT/'bank.npz',task.core.OUT/'coarse-bank.npz']
    task.core.ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',bindings={p.relative_to(task.h.ROOT).as_posix():task.h.digest(p) for p in paths},
        preserved_failure='Newmark heat-only velocity time history: 32/64 relative difference 0.021965831943863916 and order 0.6167685881501919; scalar passed. Combined-source success does not repair this failure.',
        change='Two-stage Radau IIA on exactly the same K,D,source and analytic microscopic lift. Third order and L-stable; no changed spatial mesh, source, time counts, boundary or accuracy gate.',
        paths=['heat-16','heat-32','heat-64','coarse-bank-64','outer-64','combined-64'],
        gates=dict(time_relative=.02,time_order=1.5,coefficient_relative=.02,outer_relative=.002,linear_residual=1e-9,heat_balance=2e-13),
        budget=dict(pilot_steps=8,production_steps=304,paths=6,hard_seconds=180,CPU_workers=1,BLAS_threads=1,automatic_expansion=False)))
    task.core.ex.write(OUT/'control.json',control());row=solve(8,task.core.OUT/'bank.npz',label='pilot')
    setup=row['seconds']-row['step_seconds'];forecast=(6*setup+304/8*row['step_seconds'])*1.5
    task.core.ex.write(OUT/'pilot-budget.json',dict(setup_seconds=setup,step_seconds=row['step_seconds'],forecast_seconds=forecast,
        assumption='Six preparations plus 304 measured steps, with 50 percent overhead allowance; same 180 second cap.'))
    print('FORECAST',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text())
    for p,s in plan['bindings'].items():assert task.h.digest(task.h.ROOT/p)==s,p
    assert json.loads((OUT/'pilot-budget.json').read_text())['forecast_seconds']<170
    start=time.monotonic();cases={}
    for steps in [16,32,64]:cases[str(steps)]=solve(steps,task.core.OUT/'bank.npz',label=f'heat-{steps}')
    cases['coarse']=solve(64,task.core.OUT/'coarse-bank.npz',label='coarse-bank-64')
    cases['outer']=solve(64,task.core.OUT/'bank.npz',outer=3,label='outer-64')
    cases['combined']=solve(64,task.core.OUT/'bank.npz',True,label='combined-64')
    comparisons={}
    for field in ['velocity_mass_RMS_m_s','scalar_mass_RMS']:
        series=lambda name:np.array([x[field] for x in cases[name]['history']])
        a,b,c=series('16'),series('32'),series('64');norm=max(abs(c).max(),1e-100)
        d1=np.max(abs(a-b[::2]))/norm;d2=np.max(abs(b-c[::2]))/norm
        comparisons[field]=dict(time_previous=float(d1),time_last=float(d2),order=float(np.log2(d1/d2)),
            coefficients=float(np.max(abs(c-series('coarse')))/norm),outer=float(np.max(abs(c-series('outer')))/norm))
    passed=all(x['time_last']<.02 and x['order']>1.5 and x['coefficients']<.02 and x['outer']<.002 for x in comparisons.values())
    passed=passed and max(x['heat_telescoping'] for x in cases.values())<2e-13
    result=dict(classification='Counterexample candidate',passed=passed,comparisons=comparisons,seconds=time.monotonic()-start,
        heat_endpoint=cases['64']['history'][-1],combined_endpoint=cases['combined']['history'][-1],
        physical_core_interface_closed=False,full_temperature_feedback=False,spatial_error_certified=False,photon_transport=False,full_dynamic_charge_solved=False)
    task.core.ex.write(OUT/'result.json',result);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
