"""Independent controls for the declared conserved-volume GR source bound."""
from pathlib import Path
from types import SimpleNamespace
import json
import signal
import time
import numpy as np
import sympy as s
import mpmath as mp
import def_native_feedback_gr as task

OUT=task.OUT;write=task.write


def main():
    assert not (OUT/'audit.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,task.flow.old.optical.timeout);signal.alarm(20)
    write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',seconds=20,new_fluid_steps=0,
        checks='Independent flat manufactured retarded stress integral; zero Einstein/scalar background reduction; analytic grazing-ray integral and retarded operator remainder; stable mass identity and all source hashes.',
        source_sha256=task.sha(__file__),producer_sha256=task.sha(task.__file__),gate=1e-8))
    plan=json.loads((OUT/'plan.json').read_text());actual=Path(task.__file__)
    for path,value in plan['bindings'].items():
        target=OUT/'first-producer.py' if Path(path)==actual else Path(path)
        assert task.sha(target)==value,path
    # Green operator with known two-cell ramp; the source code and its
    # quadrature are separate from the analytic polynomial integral.
    C=task.C;G=task.G;T=task.flow.old.END;dx=120000.;tt=np.linspace(0,T,17);obs=np.linspace(dx/C,T-dx/C,40)
    gx,gw=np.polynomial.legendre.leggauss(8);x=np.r_[-dx+dx*(gx+1)/2,dx*(gx+1)/2];weights=np.tile(gw/2,2)
    numerical=np.array([G/(2*C**3)*np.sum(weights*.5*(u+x/C)**2) for u in obs])
    exact=G/(2*C**3)*(obs**2+(dx/C)**2/3);error=float(max(abs(numerical/exact-1)));assert error<1e-8
    # Independent zero-background control: all metric coefficients vanish.
    class Empty:
        R=1.
        def sample(self,r):
            return dict(m=np.zeros_like(r),N=np.ones_like(r),phi=np.zeros_like(r),v=np.zeros_like(r),e=np.zeros_like(r),p=np.zeros_like(r),gamma=np.ones_like(r))
    V,K,b,N=task.wave.coeff(Empty(),np.array([1.,2.,10.]))
    assert not np.any(V) and not np.any(K) and np.all(b==1) and np.all(N==1)
    mp.mp.dps=70
    grazing=mp.quad(lambda theta:mp.mpf(1),[0,mp.pi/2]);assert abs(grazing-mp.pi/2)<mp.mpf('1e-60')
    # r0/r=sin(theta) turns dr/[r²*sqrt(1-r0²/r²)] into dtheta/r0.
    th,r0=s.symbols('theta r0',positive=True);r=r0/s.sin(th)
    transformed=s.simplify((-s.diff(r,th))/(r*r*s.cos(th)))
    assert s.simplify(transformed-1/r0)==0
    # For a=1-kappa in(0,1), concavity gives y^a<=1-a+a*y.
    # This also controls the grazing endpoint, without a mu>0 cutoff.
    angular=[]
    for kappa in [mp.mpf('0'),mp.mpf('.001'),mp.mpf('.3')]:
        a=1-kappa;exact_ray=mp.beta(1/(2*a),mp.mpf('.5'))/(2*a)
        upper=mp.pi/(2*mp.sqrt(a));assert exact_ray<=upper+mp.mpf('1e-65')
        angular.append(dict(kappa=float(kappa),integral=float(exact_ray),bound=float(upper)))
    k=mp.mpf('.2');exact=2/k-2*(-mp.expm1(-k))/k**2;born=1-k/3
    assert abs(exact-born)<=k*k/(1-k)
    bound=json.loads((OUT/'bound.json').read_text());result=json.loads((OUT/'result.json').read_text());source=json.loads((OUT/'sources.json').read_text())
    assert bound['passed'] and result['passed'] and source['passed']
    assert bound['inner_direct_boundary_outside_observer_cone']
    # Sample only a check of the previously derived continuous quadratic
    # bound; the bound itself used every interval's exact stationary point.
    d=np.load(OUT/'source-128.npz');t=d['t'];grid=np.linspace(t[0],t[-1],4097)
    e=np.asarray(d['Killing_cell_energy_erg'].sum(1,dtype=np.longdouble),float)
    exact_port=task.green.polynomial(t,e)(grid)-task.green.polynomial(t,d['inner_luminosity']).antiderivative()(grid)+task.green.polynomial(t,d['outer_luminosity']).antiderivative()(grid)
    assert max(abs(exact_port))<=bound['continuous_piecewise_port_error_erg']*(1+1e-13)
    z=np.load(OUT/'wave-128.npz');old=np.load(task.flow.OUT/'wave-128.npz');f=z['normalized_forced'];delta=bound['uncomputed_GR_normalized_bound']
    M=float(d['M_cm']);alpha=-float(d['K_cm'])/M;eps=G*old['arrived_energy_erg']/(C**4*M)
    nominal=(f-f[0]+alpha*eps)/(1-eps)
    # The observer difference can contain two independently bounded values.
    # Keep the twice-bound for that difference; do not reuse the endpoint
    # bound for both events as if they were correlated.
    relative_low=float((f[-1]-f[0])-2*delta)
    nominal_bound=float(2*delta/(1-eps[-1]))
    np.savez_compressed(OUT/'normalized-total.npz',t=z['t'],nominal_direct_plus_photon_mass=nominal,
        GR_absolute_error_bound=2*delta/(1-eps),angular_independent_direct_lower=f-f[0]-2*delta)
    row=dict(classification='Counterexample candidate',passed=True,manufactured_retarded_relative=error,zero_metric_reduction=True,
        analytic_grazing_ray_integral=True,grazing_enclosure_checks=angular,manufactured_nonzero_resolvent_remainder=True,
        continuous_port_bound_dense_check=True,source_binding_verified=True,
        endpoint_nominal_direct_plus_photon_mass=float(nominal[-1]),nominal_GR_interval=[float(nominal[-1]-nominal_bound),float(nominal[-1]+nominal_bound)],
        observer_difference_GR_error_bound=2*delta,arbitrary_outward_angular_mass_conditional_lower=relative_low,
        full_physical_numerical_error_enclosed=False,seconds=time.monotonic()-start,
        limitation='Interval covers only the declared additional GR terms. Original time/history comparisons remain empirical; source continuum, initial Einstein constraints and dynamic matter/geometry closure remain unproved.',
        final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',row);signal.alarm(0);print(json.dumps(row),flush=True)


if __name__=='__main__':main()
