"""Independent shooting and core-data binding for the momentary comparison."""
from pathlib import Path
import json, resource, time
import numpy as np
from scipy.integrate import solve_ivp
import def_retained_static_response as p


def main():
    start=time.monotonic();cpu=time.process_time();p.native.deadline(45)
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    plan=p.OUT/'independent-plan.json';assert not plan.exists()
    p.write(plan,dict(classification='Counterexample candidate',seconds=45,
        claim='Check the whole-star homogeneous static susceptibility with an independent outward ODE shoot, and bind the actual restored core Gamma table.',
        gates=dict(shooting=1e-6,static_mass_work=1e-12),
        scope='Independent discretization of the same frozen first-variation coefficients, not physical EOS or hydrostatic/thermal equilibrium certification.',
        bindings={str(f):p.sha(f) for f in [Path(__file__),Path(p.__file__),p.OUT/'result.json',p.OUT/'audit.json']}))
    model=p.construct();s=p.Static(model,8)
    core_file=p.gr.flow.initial.previous.chem.prior.star.OUT/'background.npz'
    core=np.load(core_file)
    for key in ['radius_cm','gamma1']:assert np.array_equal(core[key],model.core[key]),key
    np.savez_compressed(p.OUT/'core-coefficients.npz',radius_cm=model.core['radius_cm'],gamma1=model.core['gamma1'])
    np.savez_compressed(p.OUT/'static-coefficients.npz',radius_cm=s.r,a=s.a,W=s.W,
        exact_exterior_derivatives=s.der,asurf=s.asurf,R=s.R,den=s.den,qs=s.qs,bsurf=s.bsurf,M=s.M,K=s.K)
    x=np.r_[0.,np.asarray(s.r.ravel()/s.R,float),1.]
    a=np.r_[float(s.a.flat[0]),np.asarray(s.a.ravel(),float),float(s.asurf)]
    W=np.r_[0.,np.asarray(s.W.ravel(),float),float(s.W.flat[-1])]
    def rhs(xi,y):
        aa=np.interp(xi,x,a);ww=np.interp(xi,x,W)
        return [y[1]/(aa*xi*xi),ww*y[0]]
    first=x[1];D0=W[1]*first/3
    sol=solve_ivp(rhs,[first,1],[1.,D0],method='DOP853',rtol=3e-11,atol=[2e-13,1e-15])
    assert sol.success
    f,D=sol.y[:,-1];incident=s.den*f+s.der[2,1]*D/s.asurf
    du=s.bsurf*s.qs*f/incident;dv=D/(s.asurf*incident)
    dM,dK,_=s.R*(s.der@np.array([du,dv],np.longdouble))
    charge=float(-dK/s.M+s.K*dM/s.M**2)
    record=next(r for r in p.read(p.OUT/'result.json')['rows'] if r['case']=='unit_static_incident' and r['order']==8)
    err=abs(charge/record['normalized_charge']-1)
    work=float(abs(dM/s.K-1));assert err<1e-6 and work<1e-12,(err,work)
    np.savez_compressed(p.OUT/'independent-shooting.npz',x=sol.t,states=sol.y)
    result=dict(classification='Counterexample candidate',passed=True,shooting_charge=charge,
        shooting_relative=err,static_mass_work_relative=work,ODE_evaluations=sol.nfev,
        core_source_path=str(core_file),core_source_sha256=p.sha(core_file),
        response_clock_zero_is_saved_precision_only=True,background_time_convergence_recomputed=False,
        contraction_scope='Floating arithmetic row-sum estimate for the saved finite operator, not an outward-rounded continuum certificate.',
        seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        full_goal_complete=False)
    p.write(p.OUT/'independent.json',result)
    print(json.dumps(result))


if __name__=='__main__':main()
