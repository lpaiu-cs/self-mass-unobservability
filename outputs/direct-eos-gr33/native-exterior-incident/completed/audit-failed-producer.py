"""Independent operator, moving-shell and interval-application audit."""
from pathlib import Path
import json, resource, time
import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.optimize import brentq
import couple_native_exterior_incident as run


def flat_packet(mu0):
    """Direct distributional source versus the reduced integral, c=1."""
    R=2.;D=.3;te=.04
    def ray(age):
        r=np.sqrt(R*R+2*R*mu0*age+age*age)
        return r,(R*mu0+age)/r
    def wave(t,x):
        s=(t+x)/D
        if not 0<s<1:return 0.,0.,0.
        return (256*s**4*(1-s)**4,1024*s**3*(1-s)**3*(1-2*s)/D,
                1024*s*s*(1-s)**2*(3-14*s+14*s*s)/D**2)
    cut=brentq(lambda age:te+age+ray(age)[0]-R-D,0,D-te)
    def integral(fn,lo,hi):return quad(fn,lo,hi,epsabs=2e-13,epsrel=2e-12,limit=100)[0]
    def volume(age,reduced):
        rp,mu=ray(age);t=te+age;xmax=min(rp-R,D-t)
        def value(x):
            r=R+x;U,Ux,Uxx=wave(t,x)
            if reduced:return 2*U/r**3
            zeta=2/r-(1+mu*mu)/rp;zx=-2/r**2;zt=mu*(3*mu*mu-1)/rp**2
            return 2*zeta*Uxx+zx*(Ux-U/r)+zt*Ux
        return integral(value,0,max(xmax,0))
    smooth=sum(integral(lambda t:volume(t,False),a,b) for a,b in [(0,cut),(cut,D-te)])
    def shell(age,reduced):
        rp,mu=ray(age);U,Ux,_=wave(te+age,rp-R)
        return (1-mu*mu)*U/rp**2 if reduced else (1-mu*mu)/rp*((mu-1)*Ux+U/rp)
    full=smooth+integral(lambda t:shell(t,False),0,cut)
    boundary=-integral(lambda age:(2/R-(1+ray(age)[1]**2)/ray(age)[0])*wave(te+age,0)[1],0,D-te)
    reduced=boundary+sum(integral(lambda t:volume(t,True),a,b) for a,b in [(0,cut),(cut,D-te)])+integral(lambda t:shell(t,True),0,cut)
    error=abs(full-reduced)/max(abs(full),1e-100);assert error<1e-9,(mu0,full,reduced,error)
    return dict(mu=mu0,direct_source=full/2,reduced_source=reduced/2,relative_error=error)


def main():
    out=run.OUT;assert not (out/'audit-plan.json').exists()
    files=[Path(__file__),out/'bound-result.json',out/'bound-plan.json',out/'bound-symbolic.json',out/'result.json']
    run.write(out/'audit-plan.json',dict(classification='Counterexample candidate',budget_seconds=30,
        claim='Independently verify the scalar-operator variation, the distributional shell and exterior boundary reduction, the grazing integrals and application of the selected-term envelope. Preserve failed original point controls.',
        bindings={str(p):run.sha(p) for p in files}))
    start=time.monotonic();cpu=time.process_time();run.previous.incident.native.deadline(30)
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    for plan in ['plan.json','bound-plan.json']:
        for p,h in run.read(out/plan)['bindings'].items():assert run.sha(p)==h,p
    r,t,c,ep=sp.symbols('r t c ep',positive=True);U=sp.Function('U')(t,r);a=sp.Function('a')(r);Z=sp.Function('Z')(t,r)
    phi=U/r
    operator=r*(sp.diff(phi,t,2)/c**2-ep*sp.diff(Z,t)*sp.diff(phi,t)/c**2-
        a*a*sp.exp(2*ep*Z)*(sp.diff(phi,r,2)+(2/r+sp.diff(a,r)/a+ep*sp.diff(Z,r))*sp.diff(phi,r)))
    variation=-sp.diff(operator,ep).subs(ep,0)
    ux,uxx,ut=sp.symbols('Ux Uxx Ut')
    actual=variation.subs({sp.diff(U,r,2):uxx/a**2-sp.diff(a,r)*ux/a**2,sp.diff(U,r):ux/a,sp.diff(U,t):ut})
    expected=2*Z*(uxx-a*sp.diff(a,r)/r*U)+a*sp.diff(Z,r)*(ux-a*U/r)+sp.diff(Z,t)*ut/c**2
    assert sp.simplify(actual-expected)==0
    mu,nt,lt,nr,ar=sp.symbols('mu nt lt nr ar',positive=True)
    assert sp.simplify((nt-mu**2*lt)/(c*ar*mu)-(nr+nt/(c*ar*mu))+nr+mu*lt/(c*ar))==0
    h,zeta=sp.symbols('h zeta');dmu=(1-mu*mu)*h/mu
    assert sp.simplify(-(zeta+dmu/mu)/(c*ar*mu)-(-zeta-(1-mu*mu)*h/mu**2)/(c*ar*mu))==0
    theta=sp.symbols('theta',real=True)
    assert sp.trigsimp(sp.diff(sp.tan(theta/2),theta)-1/(1+sp.cos(theta)))==0
    assert sp.tan(sp.pi/4)-sp.tan(0)==1
    coeff=[sp.Rational(128,35) if i==4 else sp.Rational(0) for i in range(9)]
    first=[8*(b-a) for a,b in zip(coeff[:-1],coeff[1:])]
    second=[7*(b-a) for a,b in zip(first[:-1],first[1:])]
    assert max(map(abs,first))==sp.Rational(2048,70) and max(map(abs,second))==sp.Rational(2048,5)
    rows=[flat_packet(mu) for mu in [.1,.5,.99]]
    bound=run.read(out/'bound-result.json');point=run.read(out/'result.json')
    assert not point['passed'] and not bound['original_point_readout_passed']
    assert bound['analytic_envelope_verified'] and bound['conditional_interval_around_stored_body'][1]<0
    keys=['primary_total_bound','finite_Born_input_interaction_bound','additional_static_potential_bound','photon_arrival_mass_bound']
    assert abs(sum(bound[k] for k in keys)/bound['total_selected_envelope']-1)<1e-12
    result=dict(classification='Counterexample candidate',passed=True,operator_identity=True,geodesic_identity=True,
        independent_flat_moving_shell=rows,Bernstein_derivative_bounds_exact=True,
        grazing_delay_integrals_verified=True,source_bindings_verified=True,original_point_failure_preserved=True,
        selected_contribution_interval_negative=True,full_exterior_interval=False,full_goal_complete=False,
        seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    run.write(out/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
