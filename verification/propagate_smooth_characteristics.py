"""Integrate identical exterior null variations without metric time jets.

Counterexample candidate. k=h+mu^2*lambda removes the known lambda_t from
the explicit ray ODE. The metric, source knots, horizon and gates stay fixed.
"""
from pathlib import Path
import inspect,json,resource,sys,textwrap,time
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import propagate_returned_moments as previous

b=previous.base; p=b.p; C=p.C; G=p.G; LD=p.LD
OUT=Path('native-smooth-photon263-work'); OLD=previous.OUT
b.OUT=p.OUT=OUT
read,write,sha=b.read,b.write,b.sha


class Metric(b.Metric):
    def at(self,t,r):
        z=super().at(t,r); t,r=np.broadcast_arrays(t,r); shape=t.shape; d=self.d
        t=t.ravel(); r=r.ravel(); N=z['N'].ravel(); a=z['a'].ravel()
        anisotropy=np.zeros(len(t),LD)
        for k,mu0 in enumerate(self.mu):
            v=np.clip((r-d.r0)/(self.reach[k]-d.r0),0,1)
            age=d.T*self.inverse[k].sol(v)[0]; age=np.where(r>self.reach[k],2*d.T,age)
            direction=np.sqrt(np.maximum(1-(1-mu0*mu0)*(N*d.r0/(d.z0['lapse'][0]*r))**2,0))
            L=self.emission.value(t-age,np.full(len(t),self.bins[k],int))
            anisotropy+=self.aw[k]*L*(1-direction*direction)/direction
        extra=LD(G)/LD(C)**4*anisotropy/(C*r*a*a)
        z['lamp']=z['nup']-2*z['lam']/(r.reshape(shape)*z['b'])+np.asarray(extra,float).reshape(shape)
        return z


class Photons(p.Photons):
    def integrate(self,age,te,owner,scale):
        n=len(te); r0,mu0=self.ray(np.zeros(n),owner)
        launch=self.metric.at(te,r0)
        initial=np.r_[np.zeros(n),mu0*mu0*launch['lam']/scale]
        def rhs(s,y):
            rr,mu=self.ray(age*s,owner); z=self.metric.at(te+age*s,rr)
            mdot=(1-mu*mu)*C*z['a']*(1/rr-z['nr']); speed=C*z['a']*mu
            h=y[n:]-mu*mu*z['lam']/scale
            delay=-z['zeta']/scale-(1-mu*mu)*h/(mu*mu)
            kdot=(-speed*z['nup']+2*mu*mdot*z['lam']+mu*mu*speed*z['lamp'])/scale
            return np.r_[age*delay,age*kdot]
        sol=solve_ivp(rhs,[0.,1.],initial,method='DOP853',rtol=2e-8,atol=2e-11)
        assert sol.success,sol.message
        r,mu=self.ray(age,owner); z=self.metric.at(te+age,r)
        h=sol.y[n:,-1]-mu*mu*z['lam']/scale
        # The work below follows the exact integrated kinematic identity.
        # It is not an independent numerical work-integral check. That check
        # remains the unchanged direct-jet pilot and the comparison below.
        work=(z['nu']-launch['nu'])/scale+h
        result=type('Solution',(),{})()
        result.success=True; result.message=sol.message; result.nfev=sol.nfev
        result.y=np.r_[sol.y[:n,-1],h,work][:,None]
        return result


def install():
    # Reuse the packet, stress and angular bookkeeping exactly; replace only
    # the triangular ODE solve and its now non-independent work-test label.
    source=textwrap.dedent(inspect.getsource(p.Photons.propagate))
    old="sol=solve_ivp(rhs,[0.,1.],np.zeros(3*n),method='DOP853',rtol=2e-8,atol=2e-11)"
    assert source.count(old)==1
    source=source.replace(old,'sol=self.integrate(age,te,owner,scale)')
    ns=dict(p.Photons.propagate.__globals__);exec(compile(source,__file__,'exec'),ns)
    Photons.propagate=ns['propagate'];p.Metric=Metric;b.Moments=previous.Moments


def symbolic():
    v,mu,md,lam,lr,lt,nr=sp.symbols('v mu md lam lr lt nr')
    hdot=-v*nr-mu*mu*lt
    kdot=hdot+2*mu*md*lam+mu*mu*(lt+v*lr)
    assert sp.expand(kdot-(-v*nr+2*mu*md*lam+mu*mu*v*lr))==0
    return dict(classification='Proven',passed=True,
        transform='k=h+mu^2*lambda; kdot=-v*nu_prime+2*mu*mudot*lambda+mu^2*v*lambda_prime; k_launch=mu_launch^2*lambda_launch.',
        scope='Exact change of variables for the declared linear nonstiff null characteristic ODE. Numerical equivalence needs a direct-jet control. Work reconstructed from k is not an independent energy-integral test.')


def prepare():
    assert not OUT.exists(); OUT.mkdir();old=read(OLD/'execution-plan.json')
    files=[Path(__file__),Path(previous.__file__),Path(b.__file__),Path(p.__file__),OLD/'check.json',OLD/'arithmetic-check.json']
    bindings=dict(old['bindings']);bindings.update({str(q):sha(q) for q in files})
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Finish the same returned-metric photon propagation by analytically removing time-jet jumps from the integrated state, without loosening the physical or ODE tolerances.',
        decision='Compare a completed actual direct-jet cohort, then use only a verified equivalent smoother solve for the full-period stress. If reference is absent or the gate fails, do not promote production.',
        reuse='Identical rays, metric,575source knots, accepted high/low histories and angular bins. k=h+mu^2*lambda with nonzero launch condition; original direct pilot remains unchanged.',
        gates=dict(actual_direct_comparison=.002,derivative=.002,rtol=2e-8,atol=2e-11,angular_invariant=1e-10,quadrature=.002),
        budget=dict(check_seconds=600,pilot_seconds=1800,each_path_seconds=14400,virtual_GiB=16,CPU_threads=1),
        work='Computed work is the integrated identity and must not be advertised as independent work integration. The old direct-jet pilot supplies that independent check.',
        limits='Same outgoing scalar continuation. Full scalar reciprocal/background operator, source-clock controls, coupled boundary application and final charge remain open.',
        bindings=bindings))
    write(OUT/'symbolic.json',symbolic())


def check():
    install(); p.initialize(); m=Metric(8); d=m.d
    t=np.array([.23,.48,.71,.91])*d.T;r=d.inverse(C*d.T*np.array([.01,.04,.06,.13]));dr=1e-7*C*d.T
    z=m.at(t,r);derivative=(m.at(t,r+dr)['lam']-m.at(t,r-dr)['lam'])/(2*dr)
    relative=float(np.max(abs(derivative-z['lamp']))/max(np.max(abs(z['lamp'])),1e-290));assert relative<.002,relative
    # An independent varying lambda, varying mu manufactured characteristic.
    def values(t):return .3+.1*t,.1,np.sin(3*t),3*np.cos(3*t),.02*t
    def direct(t,y):
        mu,md,lam,lt,nr=values(t);return [-nr-mu*mu*lt]
    def transformed(t,y):
        mu,md,lam,lt,nr=values(t);return [-nr+2*mu*md*lam]
    old=solve_ivp(direct,[0,1],[0],rtol=1e-11,atol=1e-13,method='DOP853')
    new=solve_ivp(transformed,[0,1],[0],rtol=1e-11,atol=1e-13,method='DOP853')
    mu,_,lam,_,_=values(1);error=float(abs(old.y[0,-1]-(new.y[0,-1]-mu*mu*lam)))
    assert error<1e-10,error
    write(OUT/'check.json',dict(classification='Counterexample candidate',passed=True,lambda_radial_relative=relative,
        manufactured_transformed_absolute=error,actual_direct_comparison_pending=True))


def pilot():
    assert read(OUT/'check.json')['passed']; install();p.initialize();m=Photons(8,8);rows=[]
    for cell in [0,7,15]:
        z,row=m.propagate(m.d.T,8,[cell]);row['work_identity_independent']=False
        np.savez_compressed(OUT/f'pilot-{cell}.npz',**z);write(OUT/f'pilot-{cell}.json',row);rows.append(row)
        write(OUT/'pilot-progress.json',dict(rows=rows))
    upper=2*max(r['seconds']/r['packets'] for r in rows)*8*32*sum(range(1,17))+120
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',rows=rows,forecast_upper_seconds_per_full_path=upper,
        cost_eligible=upper<14400,actual_direct_comparison_pending=True,eligible=False))


if __name__=='__main__':
    action=sys.argv[1];start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);p.incident.native.deadline(1800 if action=='pilot' else 600)
    try:
        if action!='prepare':
            for q,h in read(OUT/'plan.json')['bindings'].items():assert sha(q)==h,q
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
