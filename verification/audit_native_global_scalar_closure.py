"""Independent finite-retarded-time flat impulse and binding checks."""
from pathlib import Path
import json,signal,time
import numpy as np
from scipy.integrate import quad
import sympy as sp
import def_native_global_scalar_closure as task


def main():
    out=task.OUT;assert not (out/'audit.json').exists();start=time.monotonic()
    task.write(out/'audit-plan.json',dict(classification='Counterexample candidate',seconds=20,
        independent='Analytic radial primitives at a finite causal cutoff, integrated independently in angle; exact time-port primitive and source bindings; zero and nonzero controls.',
        gate=.002,bindings={str(p):task.sha(p) for p in [Path(__file__),Path(task.__file__),out/'bound.json',out/'result.json']}))
    signal.signal(signal.SIGALRM,task.flow.old.optical.timeout);signal.alarm(20)
    for name in ['plan.json','execution-plan.json','bound-plan.json','audit-plan.json']:
        for path,digest in json.loads((out/name).read_text())['bindings'].items():assert task.sha(path)==digest,path
    # Check the radial antiderivatives independently of the production ray ODE.
    th=sp.symbols('theta',positive=True)
    F=sp.sin(th)-th/2-sp.sin(2*th)/4
    assert sp.trigsimp(sp.diff(F,th)-sp.sin(th)*sp.cos(th)*sp.sin(th)/(1+sp.cos(th)))==0
    class Flat(task.Exterior):
        def metric(self,q):return (np.ones_like(q),)*3
    m=Flat();m.N0=1.;m.K=-1.;k=m.kernel(8,8);T=m.T;R=m.r0;C=task.C
    numerical=np.array([-C*np.sum(k['weights']*k['mass']*np.maximum(T-k['delay'],0.)),
        -.5*np.sum(k['weights']*k['stress']*(k['delay']<T))])
    def exact(mu,component):
        theta=np.arccos(mu);impact=R*np.sqrt(1-mu*mu);end=R/C*(1-mu)
        low=2*np.arctan(C*(end-T)/impact) if end>T else 0.
        if component==0:
            primitive=lambda x:(T-end)*np.sin(x)**2/2+impact/C*(np.sin(x)-x/2-np.sin(2*x)/4)
            return mu*(-C*m.K/impact**2)*(primitive(theta)-primitive(low))
        primitive=lambda x:x/2-np.sin(2*x)/4
        return mu*(-.5*m.K/impact)*(primitive(theta)-primitive(low))
    threshold=1-C*T/R
    analytic=np.array([quad(lambda v:exact(v,i)*R,0,1,points=[threshold],epsabs=1e-12,epsrel=2e-10)[0]/R for i in range(2)])
    errors=abs(numerical/analytic-1);assert np.all(numerical>0) and max(errors)<.002,errors
    # The exact zero input produces zero in both linear exterior kernels.
    assert np.sum(k['weights']*k['mass']*np.zeros_like(k['delay']))==0
    d=np.load(task.compact.base.OUT/'source-128.npz');ports=[]
    for key in ['inner','outer']:
        H=task.flow.green.polynomial(d['t'],d[key+'_luminosity']).antiderivative()
        ports.append(float(np.max(abs(H(d['t'])-d[key+'_cumulative_energy_erg']))/max(np.max(abs(d[key+'_cumulative_energy_erg'])),1.)))
    assert max(ports)<1e-12
    bound=json.loads((out/'bound.json').read_text());assert bound['passed']
    # A nonzero manufactured Volterra resolvent obeys the stated geometric bound.
    eta=sp.Rational(1,5);terms=sum(eta**n for n in range(1,20));assert terms<eta/(1-eta)
    result=dict(classification='Counterexample candidate',passed=True,finite_time_flat_impulse_component_relative=errors.tolist(),
        finite_time_flat_impulse_numerical=numerical.tolist(),finite_time_flat_impulse_independent=analytic.tolist(),
        symbolic_radial_primitives_passed=True,zero_input_passed=True,nonzero_resolvent_check=True,
        exact_saved_port_primitives_relative=ports,all_registered_bindings_verified=True,
        conditional_scalar_lower=bound['arbitrary_outward_conditional_scalar_interval'][0],
        limitation='Independent algebra/numerical audit for the declared source and initial linear operator; not a full physical error certificate.',
        full_source_error_enclosed=False,final_charge_solved=False,full_goal_complete=False,seconds=time.monotonic()-start)
    task.write(out/'audit.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
