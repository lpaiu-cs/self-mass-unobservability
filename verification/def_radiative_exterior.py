"""Outgoing dynamic scalar/metric vacuum boundary on a finite Just background.

No matter approximation is made in this vacuum module; the interior fluid and
physical atmosphere remain to be connected. exp(-i omega t), r_star(R)=0.
"""
from pathlib import Path
import argparse
import json
import time
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
import def_normalized_charge as q

e,s=q.e,q.s
OUT=e.g.OUT/'def-radiative-exterior'


def symbolic():
    r,m,p,z,zr=sp.symbols('r m p z zr',nonzero=True)
    b=1-2*m/r;mp_=r*r*b*p*p/2
    nup=m/(r*r*b)+r*p*p/2
    lamp=(mp_/r-m/r**2)/b
    pp=-(2/r+nup-lamp)*p
    assert sp.simplify(nup-lamp-2*m/(r*r*b))==0
    dlambda=r*p*z
    dnuprime=p*z/b+r*p*zr
    dlambdaprime=(p+r*pp)*z+r*p*zr
    assert sp.simplify((dnuprime-dlambdaprime)*p-2*p*p*z/b)==0
    # Differentiate the radial mass constraint using dm=r^2 b p z.
    dm=r*r*b*p*z
    dmr=sp.diff(dm,r)+sp.diff(dm,m)*mp_+sp.diff(dm,p)*pp+sp.diff(dm,z)*zr
    assert sp.simplify(dmr-(-r*p*p*dm+r*r*b*p*zr))==0
    return dict(classification='Proven',passed=True,
        background='b=1-2m/r; F=N sqrt(b); m_prime=r^2 b Phi^2/2; (ln F)_prime=2m/(r^2 b).',
        metric_constraint='For nonzero frequency dm=r^2 b Phi dphi and dlambda=r Phi dphi; the homogeneous time-independent mass perturbation is outside this radiative sector.',
        radial_equation='dphi_rr+(2/r+2m/(r^2 b))*dphi_r+(omega^2/F^2+2 Phi^2/b)*dphi=0.',
        wave_potential='u=r dphi; u_rstar_rstar+(omega^2-V)u=0; V=2 N^2(m/r^3-Phi^2).',
        scope='Linear vacuum scalar plus metric perturbation around a static massless-scalar background; no cold-fluid assumption is imported.')


def mul(a,b,n):return np.convolve(a,b)[:n]


def inverse(a,n):
    z=np.zeros(n,dtype=np.result_type(a,complex));z[0]=1/a[0]
    for j in range(1,n):z[j]=-sum(a[k]*z[j-k] for k in range(1,min(j+1,len(a))))/a[0]
    return z


def series(M,K,omega,order=8):
    """Outgoing h after u=exp(i omega r_star)h, expanded at x=R/r=0."""
    n=order+3;mass=np.zeros(n,complex);F=np.zeros(n,complex);mass[0]=M;F[0]=1
    for j in range(order+1):
        b=np.r_[1,-2*mass[:-1]];ib=inverse(b,n);iF=inverse(F,n)
        mass[j+1]=(-K*K/2*mul(b,mul(iF,iF,n),n))[j]/(j+1)
        F[j+1]=(-2*mul(mul(mass,F,n),ib,n))[j]/(j+1)
    b=np.r_[1,-2*mass[:-1]];ib=inverse(b,n);iF=inverse(F,n)
    B=-2j*omega*iF;B[1]+=2;B[2:]-=2*mul(mass,ib,n)[:-2]
    C=np.zeros(n,complex);C[1:]-=2*mul(mass,ib,n)[:-1]
    C[2:]+=2*K*K*mul(ib,mul(iF,iF,n),n)[:-2]
    h=np.zeros(n,complex);h[0]=1
    for j in range(order+(not bool(omega))):
        dh=np.arange(1,n)*h[1:];ddh=np.arange(1,n-1)*dh[1:]
        coeff=mul(B,dh,n)[j]+mul(C,h,n)[j]
        if j>=2:coeff+=ddh[j-2]
        if omega:h[j+1]=-coeff/(B[0]*(j+1))
        else:
            # At omega=0 the x^j equation fixes h_j, j>=1.
            if j:h[j]=-coeff/(j*(j+1))
    return mass,F,h[:order+1]


def background(mu,flux,tol):
    mp.mp.dps=60;ext=q.exterior.exact(mp.mpf(str(mu)),mp.mpf(str(flux)))
    M=float(mp.mpf(str(mu))+mp.mpf(str(flux))**2*ext[0])
    K=float(mp.mpf(str(flux))*ext[1])
    def rhs(x,y):
        m,F=y;b=1-2*m*x
        return [-.5*b*K*K/(F*F),-2*m*F/b]
    sol=solve_ivp(rhs,(0,1),[M,1],rtol=tol,atol=tol*1e-3,dense_output=True,method='DOP853')
    assert sol.success and abs(sol.y[0,-1]-mu)<max(1e-14,20*tol*mu)
    assert abs(sol.y[1,-1]-float(ext[1]))<20*tol
    return M,K,sol


def outgoing(mu,flux,omega,extent=40,tol=2e-12):
    M,K,bg=background(mu,flux,tol)
    _,_,h=series(M,K,omega)
    x0=omega/extent if omega else 1e-5
    h0=np.polynomial.polynomial.polyval(x0,h)
    dh0=np.polynomial.polynomial.polyval(x0,np.arange(1,len(h))*h[1:])
    def rhs(x,y):
        m,F=bg.sol(x);b=1-2*m*x
        B=2/x-2*m/b-2j*omega/(F*x*x)
        C=-2*m/(b*x)+2*K*K/(b*F*F)
        return [y[1],-B*y[1]-C*y[0]]
    sol=solve_ivp(rhs,(x0,1),[h0,dh0],rtol=tol,atol=tol*1e-4,method='DOP853')
    assert sol.success
    hr,hx=sol.y[:,-1];F=bg.y[1,-1]
    Z=1j*omega/F-hx/hr-1
    flux_current=float(np.imag(np.conj(hr)*(1j*omega*hr-F*hx)))
    return dict(impedance=Z,amplitude_transfer=1/hr,h=hr,
        current=flux_current,current_relative_error=abs(flux_current/omega-1) if omega else 0,
        rhs_calls=sol.nfev,background_rhs_calls=bg.nfev,cutoff_x=x0)


def static_exact(mu,flux):
    mp.mp.dps=60;m=mp.mpf(str(mu));p=mp.mpf(str(flux));b=1-2*m
    G=lambda a,c:q.exterior.exact(a,c)[2]
    charge=lambda a,c:-c*q.exterior.exact(a,c)[1]
    dphi=-(G(m,p)+p*mp.diff(lambda t:G(m,t),p))/(1+b*p*p*mp.diff(lambda t:G(t,p),m))
    dm=b*p*dphi
    dcharge=mp.diff(lambda t:charge(m,t),p)+mp.diff(lambda t:charge(t,p),m)*dm
    mass=lambda a,c:a+c*c*q.exterior.exact(a,c)[0]
    mass_change=mp.diff(lambda t:mass(t,p),m)*dm+mp.diff(lambda t:mass(m,t),p)
    assert abs(mass_change)<mp.mpf('1e-45')
    return dict(impedance=float(1/dphi),amplitude_transfer=float(dcharge/dphi),fixed_ADM_mass_verified=True)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(q.exterior.__file__),q.OUT/'result.json',q.OUT/'finite.npz',s.OUT/'companion-benchmark.json']
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',symbolic=symbolic(),
        bindings={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in files},
        source='https://arxiv.org/html/1802.07847v2 ; full metric/scalar radial perturbations, vacuum reduction derived independently here. Cold barotropic interior is not imported.',
        parameters='Phase41 finite-background surface mu and scalar flux, plus manufactured mu=0.1, q=0.03 and flat-space controls.',
        frequency_set='0 and first three declared orbital harmonics for the actual surface; 0,0.05,0.2 on manufactured strong vacuum; 0.001,0.2 flat.',
        boundary='Outgoing exp(-i omega t + i omega r_star); r_star(R)=0. Eighth-order asymptotic h with u=exp(i omega r_star)h. Return R*dphi_prime/dphi and asymptotic u/(R*dphi_surface).',
        comparisons=[dict(extent=40,tol=2e-12),dict(extent=80,tol=2e-13)],
        gates=dict(absolute_impedance_and_transfer_difference=2e-9,current_relative=2e-8,static_exact=2e-10,flat_exact=2e-10),
        budget=dict(hard_timeout_seconds=90,maximum_frequency_solves=18,native_calls=0,new_fluid_steps=0,automatic_expansion=False),
        decision='Install dynamic radiative boundary only if static Just, flat outgoing, conserved flux and outer/tolerance comparisons pass. No full dynamic-charge completion without coupled fluid, background and benchmark response.'))


def run():
    began=time.monotonic();plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert e.digest(s.ROOT/rel)==sha,rel
    assert not (OUT/'result.json').exists()
    saved=np.load(q.OUT/'finite.npz');R=float(np.load(s.OUT/'initial.npz')['rf'][-1]);mu=float(saved['mf'][-1]/R)
    prior=json.loads((q.OUT/'result.json').read_text());finite=next(x for x in prior['rows'] if x['name']=='finite_background')
    flux=finite['surface_flux_q'];omega=json.loads((s.OUT/'companion-benchmark.json').read_text())['leading_drive']['omega_R_over_c']
    cases=[('stellar',mu,flux,w) for w in [0,omega,2*omega,3*omega]]
    cases += [('strong',.1,.03,w) for w in [0,.05,.2]]+[('flat',0,0,w) for w in [.001,.2]]
    rows=[]
    for name,m,p,w in cases:
        values=[]
        for setting in plan['comparisons']:
            if name=='flat':
                # The same solver supports the flat limit with a minimal empty metric.
                original=globals()['background']
                def flat(*_):
                    class Metric:
                        y=np.array([[0.,0.],[1.,1.]])
                        nfev=0
                        def sol(self,x):return np.array([0.,1.])
                    return 0.,0.,Metric()
                globals()['background']=flat
                try:value=outgoing(m,p,w,**setting)
                finally:globals()['background']=original
            else:value=outgoing(m,p,w,**setting)
            values.append(value)
        difference=max(abs(values[0][k]-values[1][k]) for k in ['impedance','amplitude_transfer'])
        control=static_exact(m,p) if not w else None
        error=max(abs(values[-1][k]-control[k]) for k in ['impedance','amplitude_transfer']) if control else 0.
        flat_error=max(abs(values[-1]['impedance']-(1j*w-1)),abs(values[-1]['amplitude_transfer']-1)) if name=='flat' else 0.
        passed=difference<2e-9 and max(v['current_relative_error'] for v in values)<2e-8 and error<2e-10 and flat_error<2e-10
        row=dict(case=name,mu=m,scalar_flux=p,omega_R_over_c=w,refinement_difference=float(difference),static_exact_error=float(error),flat_error=float(flat_error),passed=bool(passed),static_control=control,
            solves=[{k:[float(v.real),float(v.imag)] if isinstance(v,complex) else v for k,v in z.items()} for z in values])
        rows.append(row);e.write(OUT/'progress.json',dict(completed=len(rows),total=len(cases),last=row));print('VACUUM',name,w,passed,difference,error,flush=True)
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows,seconds=time.monotonic()-began,
        coupled_interior_solved=False,physical_surface_solved=False,full_dynamic_charge_solved=False,symbolic=symbolic())
    e.write(OUT/'result.json',result)
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in OUT.iterdir() if p.is_file() and p.name not in ['manifest.json','run.log']}))
    print('FINAL',result['passed'],result['seconds'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
