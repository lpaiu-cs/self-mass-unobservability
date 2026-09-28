"""Moving-emitter photon contribution on the fixed native vacuum metric.

Counterexample candidate: apply the measured material surface history to
retarded photon moments and their coupled GR response. This isolates the
kinematic contribution; full moving stress/lapse junction closure remains open.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline
from types import SimpleNamespace
import def_native_balanced_flux as task

OUT=task.OUT
write=task.write


class Rays:
    def __init__(self,metric,r,angles=48,order=32):
        self.metric=metric;self.r=np.asarray(r)
        x,w=np.polynomial.legendre.leggauss(angles);self.mu0=mu0=(x+1)/2
        self.mu,self.weights,self.delay=metric.rays(self.r,angles,order)
        m,n,b,c,v=metric.metric(self.r)
        ms,ns,bs,cs,vs=[float(a[0]) for a in metric.metric(np.array([1.]))]
        self.cs=cs;ks=1-ms/bs-vs**2/2;self.ks=ks
        k=1-m/(self.r*b)-self.r**2*v**2/2
        f=self.r/(c*k)
        span=((self.r-1)*(self.r+1))[:,None]/(np.sqrt(((self.r-1)*(self.r+1))[:,None]+mu0**2)+mu0)
        y,wq=np.polynomial.legendre.leggauss(order)
        z=mu0[None,:,None]+span[:,:,None]*(y+1)/2
        rad=np.sqrt(z*z+(1-mu0**2)[None,:,None])
        mm,nn,bb,cc,vv=[a.reshape(rad.shape) for a in metric.metric(rad.ravel())]
        mr=rad**2*bb*vv**2/2;cr=2*mm/(rad**2*bb)
        kk=1-mm/(rad*bb)-rad**2*vv**2/2
        kr=(mm-rad*mr)/(rad**2*bb**2)+rad*vv**2+rad**2*cr*vv**2
        fp=rad/(cc*kk)*(1/rad-cr-kr/kk)
        mulocal=np.sqrt(1-(1-mu0**2)[None,:,None]/ns**2*nn**2/rad**2)
        integral=span/2*np.sum(fp*z/(rad*mulocal)*wq,axis=2)
        # Integration by parts removes the grazing 1/mu0 cancellation.
        self.DR=ks*(-f[:,None]/self.mu+integral)
        self.DB=-(1/cs+mu0*self.DR)/ks
        self.A=mu0/cs-self.DR;self.B=-self.DB/cs
        self.Cz=-(1-self.mu**2)/self.mu*ks
        self.Cv=(1-self.mu**2)/self.mu*mu0/cs
        self.n,self.b,self.c,self.v=n,b,c,v

    def moments(self,t,history):
        u=t-self.delay;active=u>0;u=np.clip(u,0,history.x[-1])
        z=history(u);zd=history(u,1);zdd=history(u,2)
        cumulative=(self.A*z+self.B*zd)*active
        flux=(self.A*zd+self.B*zdd)*active
        dmu=(self.Cz*z+self.Cv*zd)*active
        energy_pressure=((1/self.mu-self.mu)*flux-(1+1/self.mu**2)*dmu)@self.weights
        pressure=(self.mu*flux+dmu)@self.weights
        return cumulative@self.weights,energy_pressure,pressure

    def force(self,t,history,f0,R,test):
        h,ep,_=self.moments(t,history);factor=task.old.G/(task.old.C**4*R)*f0
        J=-factor*h*np.sqrt(self.b)/self.n
        EP=factor*ep/(4*np.pi*self.r**2*self.n**2)
        return test@(2*self.c*self.v*J/self.b**2-4*np.pi*self.r**3*self.c*self.v*EP/self.b)


def controls():
    import sympy as s
    r,C,k,mu,I,z,zd=s.symbols('r C k mu I z zd',nonzero=True)
    DR=-1/(C*mu)+k*I;DB=-mu*I
    assert s.simplify(DB+(1/C+mu*DR)/k)==0
    assert s.integrate(2*mu**2,(mu,0,1))==s.Rational(2,3)
    # Exact flat finite-amplitude rays provide an independent first-variation
    # check of retardation, aberration, pressure and accumulated energy.
    flat=SimpleNamespace(Nb=1.)
    flat.metric=lambda r:(np.zeros_like(r),np.ones_like(r),np.ones_like(r),np.ones_like(r),np.zeros_like(r))
    flat.rays=lambda r,a,o:task.old.photons.Exterior.rays(flat,r,a,o)
    rays=Rays(flat,np.array([1.02,1.3,2.]),32,64);t=2.5
    def zeta(u):return .02*u**2+.003*u**3
    def speed(u):return .04*u+.009*u**2
    def acc(u):return .04+.018*u
    hist=CubicHermiteSpline([0.,3.],zeta(np.array([0.,3.])),speed(np.array([0.,3.])))
    analytic=np.array(rays.moments(t,hist))
    def exact(eps):
        u=t-rays.delay;mu0=rays.mu0
        for _ in range(10):
            R=1+eps*zeta(u);beta=eps*speed(u);bdot=eps*acc(u)
            lab=(mu0+beta)/(1+beta*mu0);labdot=(1-mu0**2)/(1+beta*mu0)**2*bdot
            impact=R**2*(1-lab**2);root=np.sqrt(rays.r[:,None]**2-impact)
            delay=root-R*lab
            idot=2*R*beta*(1-lab**2)-2*R**2*lab*labdot
            ddot=-idot/(2*root)-beta*lab-R*labdot
            u-=(u+delay-t)/(1+ddot)
        mu=root/rays.r[:,None];flux=(1+mu0*beta)/(1+ddot)
        return np.array([(u+eps*mu0*zeta(u))@rays.weights,
            ((1/mu-mu)*flux)@rays.weights,(mu*flux)@rays.weights])
    fd=(exact(1e-5)-exact(-1e-5))/(2e-5)
    relative=np.max(abs(fd-analytic),axis=1)/np.max(abs(analytic),axis=1)
    assert max(relative)<2e-6,relative
    # Independent moving lower-limit and impact-parameter quadrature in the
    # actual curved background. It does not use DR/DB or their derivation.
    bg=task.prior.Background();rr=np.array([1.001,1.01,1.1]);curve=Rays(bg,rr,16,64)
    def travel(R,beta):
        mu0=curve.mu0;lab=(mu0+beta)/(1+mu0*beta)
        ns=bg.metric(np.array([R]))[1][0];impact=R**2*(1-lab**2)/ns**2
        lo=R*lab;hi=np.sqrt(rr[:,None]**2-R**2*(1-lab**2));y,w=np.polynomial.legendre.leggauss(64)
        z=lo[None,:,None]+(hi-lo)[:,:,None]*(y+1)/2
        rad=np.sqrt(z*z+R**2*(1-lab**2)[None,:,None])
        _,n,_,c,_=[a.reshape(rad.shape) for a in bg.metric(rad.ravel())]
        mu=np.sqrt(1-impact[None,:,None]*n*n/(rad*rad))
        return (hi-lo)/2*np.sum(z/(rad*c*mu)*w,axis=2)
    eps=1e-7
    rd=(travel(1+eps,0)-travel(1-eps,0))/(2*eps)
    bd=(travel(1,eps)-travel(1,-eps))/(2*eps)
    errors=[float(np.max(abs(a-b))/np.max(abs(b))) for a,b in [(rd,curve.DR),(bd,curve.DB)]]
    assert max(errors)<2e-5,errors
    return dict(classification='Proven',symbolic_passed=True,
        identity='First-order impact change is delta ln B = (1-nu_s_prime)*zeta-mu0*zeta_dot/Cs. Retarded cumulative emission includes the angular pressure-work term f0*mu0*zeta/Cs.',
        numerical_classification='Counterexample candidate',flat_moment_relative=relative.tolist(),curved_delay_relative=errors,
        scope='Fixed vacuum metric, constant reference comoving emitted power. No evolving metric, area/lapse luminosity feedback, full stress junction or final charge certificate.')


def prepare():
    assert not (OUT/'motion-plan.json').exists();result=json.loads((OUT/'result.json').read_text());assert result['time_passed']
    start=time.monotonic();check=controls()
    write(OUT/'motion-plan.json',dict(classification='Counterexample candidate',
        claim='Apply the corrected actual surface motion to the moving-ray photon energy and pressure, and evolve its additive coupled matter/scalar/metric/thermal response.',
        scope='First kinematic correction on the fixed curved background, driven by the saved baseline surface. This is NOT the full moving-surface solution: radiation stress junction, area/lapse feedback and dynamic ray metric remain to be closed.',
        source='Phase96 p4-64 physical displacement and direct velocity; normalized reference hemispheric luminosity. No prescribed synthetic motion in production.',
        paths=['motion-p4-8 pilot','motion-p4-32','motion-p4-64'],
        gates=dict(temperature=.02,velocity=.03,scalar=.03,source_history_relative=.02),
        prior_seconds=result['seconds']+3,controls_seconds=time.monotonic()-start,total_phase_budget_seconds=180,
        decision='Keep the same degree, grid, horizon and gates. Stop on budget or failed controls. A small isolated contribution is not a bound on omitted physical mechanisms.',
        controls=check,bindings={str(p.relative_to(task.old.ROOT)):task.old.photons.digest(p) for p in [Path(__file__),Path(task.__file__),OUT/'p4-32.npz',OUT/'p4-64.npz',OUT/'result.json']}))


def run():
    assert not (OUT/'motion-result.json').exists();spec=json.loads((OUT/'motion-plan.json').read_text())
    for p,h in spec['bindings'].items():assert task.old.photons.digest(task.old.ROOT/p)==h,p
    previous=spec['prior_seconds']+spec['controls_seconds'];signal.alarm(int(180-previous));start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    model=task.Model();rays=Rays(model.bg,model.points['r'][model.outside]);a=np.load(OUT/'p4-64.npz');b=np.load(OUT/'p4-32.npz')
    hist=CubicHermiteSpline(a['emission_times'],a['surface'][:,0],a['surface'][:,2])
    coarse=CubicHermiteSpline(b['emission_times'],b['surface'][:,0],b['surface'][:,2])
    f0=float(model.f0[-1]);R=model.bg.R;end=a['emission_times'][-1]
    moment=np.array([rays.moments(t,hist) for t in a['emission_times']])
    moment32=np.array([rays.moments(t,coarse) for t in a['emission_times']])
    error=np.max(abs(moment-moment32),axis=(0,2))/np.maximum(np.max(abs(moment),axis=(0,2)),1e-100)
    write(OUT/'motion-source-control.json',dict(classification='Counterexample candidate',source_history_relative=error.tolist(),passed=bool(max(error)<.02)))
    assert max(error)<.02,error
    # Solve only this additive correction in the unchanged coupled operator.
    # Keeping f0 sets luminosity normalization, not a second baseline drive.
    model.F0*=0;model.TE0*=0;model.LE0*=0;model.J0*=0
    model.photon_force=lambda t:rays.force(t,hist,f0,R,model.photon_test)
    horizon=float(end*model.bg.tc);pilot=model.evolve(horizon,8,'motion-p4-8')
    forecast=previous+time.monotonic()-start+1.5*(pilot['seconds']*96/8+5)
    write(OUT/'motion-pilot-budget.json',dict(classification='Counterexample candidate',forecast_total_seconds=forecast,
        assumption='Measured same additive-forcing eight-step path scaled to96 more steps,50percent margin, previous baseline/control costs included.'))
    assert forecast<180,'Stop: total registered phase budget would be exceeded'
    for n in [32,64]:model.evolve(horizon,n,f'motion-p4-{n}')
    c=np.load(OUT/'motion-p4-32.npz');d=np.load(OUT/'motion-p4-64.npz');errors={};relative={}
    for f in ['temperature','velocity','scalar']:
        errors[f]=float(np.max(abs(c[f]-d[f][::2]))/max(np.max(abs(d[f])),1e-100))
        relative[f]=float(np.max(abs(d[f]))/max(np.max(abs(a[f])),1e-100))
    ps=model.bg.sample(np.array([1.]));bs=1-2*ps['m'][0];Cs=ps['N'][0]*np.sqrt(bs)
    work=2*f0/(3*Cs)*a['surface'][:,0];emitted=f0*a['emission_times']
    env=np.load(task.old.prior.OUT/'final-envelope.npz');gaswork=work*float(env['Pgas'][-1]/env['Prad'][-1])
    np.savez_compressed(OUT/'motion-moments.npz',r=rays.r,t=a['emission_times'],moments=moment,photon_work_erg=work,external_gas_work_erg=gaswork,
        reference_emitted_energy_erg=emitted,source_history_relative=error)
    write(OUT/'motion-result.json',dict(classification='Counterexample candidate',time_passed=all(errors[f]<spec['gates'][f] for f in errors),
        time_relative=errors,additive_over_baseline=relative,photon_work_endpoint_erg=float(work[-1]),external_gas_work_endpoint_erg=float(gaswork[-1]),
        photon_work_over_emitted=float(work[-1]/emitted[-1]),seconds=time.monotonic()-start,accounted_total_seconds=previous+time.monotonic()-start,
        actual_surface_history_applied=True,additive_coupled_GR_evolved=True,moving_surface_fully_solved=False,final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
