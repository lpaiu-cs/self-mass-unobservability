"""Exact Just exterior and independently integrated scalar metric backreaction."""
import json, shutil, sys
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp
import sympy as sp
import gr_scalar_piecewise_control as previous

g=previous.g;OUT=g.OUT/'gr-scalar-nonlinear-exterior'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();previous.verify()
    shutil.copy2(g.ROOT/'outputs/def1996-scalar.pdf',OUT/'def1996.pdf')
    paths=[g.ROOT/'verification/gr_scalar_nonlinear_exterior.py',OUT/'def1996.pdf',
        previous.OUT/'manifest.json',g.OUT/'initial-state-17-4.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='575fae7',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        primary_source=dict(url='https://arxiv.org/pdf/gr-qc/9602056',DOI='10.1103/PhysRevD.54.1474',
            classification='Imported from prior work',equations=['3.6a-e','3.8-3.10','3.13a-h'],
            notation='The source nu is twice log lapse; retain this factor. Its mass is the Einstein ADM mass. Its effective charge ratio is not automatically a derivative along a fixed-baryon stellar sequence.'),
        surface_fluxes=[0.,-1e-8,1e-8,-1e-4,1e-4,-.1,.1],tolerances=[1e-11,2e-13],
        comparison_relative_tolerance=1e-10,points=129,
        input='Use the current star R and surface m/R only as scales for manufactured vacuum boundary controls. q=R*phi_prime(R) is prescribed and phi(R)=0. These are not finite-amplitude interior solutions, actual background estimates or observed source charges.',
        method='Compare the exact source Just matching map in 80-digit arithmetic with an independent compactified vacuum IVP x=R/r from x=1 to x=0, retaining scalar mass/lapse increments separately. No finite outer-radius truncation.',
        scope='Resolve the nonlinear exterior mapping prerequisite only. The finite-pressure photosphere is not matched to vacuum as a physical complete star. No atmosphere, same-baryon nonlinear interior, physical EOS, temporal response or observational inference is claimed.'))


def exact(mu,q):
    mu=mp.mpf(mu);q=mp.mpf(q);f=1-2*mu
    assert 0<mu<mp.mpf('.5')
    if not q:
        # Continuous coefficients of q^2 and q about the Schwarzschild solution.
        return [-(f*f)*mp.log(f)/(4*mu),f,-f*mp.log(f)/(2*mu),f/(2*mu)+f*f*mp.log(f)/(4*mu*mu)]
    d=q*q+2*mu/f;alpha=2*q/d;Q1=mp.sqrt(1+alpha*alpha)
    arg=Q1*d/(d+2);assert abs(arg)<1
    nu=-2*mp.atanh(arg)/Q1
    mass=d*mp.sqrt(f)*mp.exp(nu/2)/2
    return [(mass-mu)/(q*q),alpha*mass/q,-alpha*nu/(2*q),(mp.log(f)-nu)/(2*q*q)]


def symbolic():
    x,R,mu,q,v,z=sp.symbols('x R mu q v z',positive=True)
    m=mu+q*q*v;f=1-2*m*x;f0=1-2*mu*x
    vp=-f*z*z/2;zp=2*m*z/f;dp=-z;wp=-v/(f*f0)-z*z*x/2
    mass_r=-q*q*x*x*vp
    assert sp.simplify(mass_r-f*q*q*z*z*x*x/2)==0
    psi_r=-(q/R**2)*x*x*(zp*x*x+2*x*z)
    expected=-2*x/R*(1-m*x)/f*(q*z*x*x/R)
    assert sp.simplify(psi_r-expected)==0
    nu_r=-x*x/R*(-mu/f0+q*q*wp)
    assert sp.simplify(nu_r-(m*x*x/(R*f)+q*q*z*z*x**3/(2*R)))==0
    assert sp.simplify(-q*dp*x*x/R-q*z*x*x/R)==0
    # Independent first asymptotic order in the physical Jordan metric.
    h,M,charge,alpha,Ainf=sp.symbols('h M charge alpha Ainf',real=True)
    conformal=1+2*alpha*charge*h
    tt=sp.series(conformal*(1-2*M*h),h,0,2).removeO()
    rr=sp.series(conformal*(1+2*M*h),h,0,2).removeO()
    assert sp.expand(tt-(1-2*(M-alpha*charge)*h))==0
    assert sp.expand(rr-(1+2*(M+alpha*charge)*h))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        compactification='m/R=mu+q^2*v; r^2*phi_prime/R=q*z; phi-phi_s=q*d; logN-logN_s=0.5*log((1-2mu*x)/(1-2mu))+q^2*w. Then v_prime=-f*z^2/2, z_prime=2*(mu+q^2*v)*z/f, d_prime=-z, w_prime=-v/(f*f0)-x*z^2/2.',
        exterior='Scalar vacuum energy is nonnegative: dm/dr=0.5*r*(r-2m)*phi_prime^2. For nonzero flux the exterior mass is not constant.',
        physical_frame='For phi=phi_infinity+charge/r+O(r^-2), normalize T=A_infinity*t and Jordan areal radius R_J=A(phi)*r. The 1/R_J coefficient lengths are A_infinity*(M_E-alpha_infinity*charge) in -g_TT and A_infinity*(M_E+alpha_infinity*charge) in g_RJRJ. They must not be silently identified with each other or M_E.',
        limitation='These algebraic vacuum/asymptotic statements do not determine a self-consistent stellar charge or actual orbit.'))


def run():
    mp.mp.dps=80;plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    symbolic();state=dict(np.load(g.OUT/'initial-state-17-4.npz'))
    R=float(state['radius_faces_m'][0]);mu=float(state['mass_faces_geom'][0]/R);rows=[]
    grid=np.linspace(1,0,plan['points'])
    for index,q in enumerate(plan['surface_fluxes']):
        expected=exact(mu,q);cases=[]
        for tol in plan['tolerances']:
            def rhs(x,y):
                v,z,d,w=y;m=mu+q*q*v;f=1-2*m*x;f0=1-2*mu*x
                assert f>0 and f0>0
                return [-f*z*z/2,2*m*z/f,-z,-v/(f*f0)-x*z*z/2]
            sol=solve_ivp(rhs,(1,0),[0,1,0,0],method='DOP853',rtol=tol,
                atol=tol*1e-3,max_step=.025,t_eval=grid)
            assert sol.success
            errors=[float(abs(mp.mpf(float(a))-b)/max(abs(b),mp.mpf('1e-40'))) for a,b in zip(sol.y[:,-1],expected)]
            passed=max(errors)<plan['comparison_relative_tolerance']
            np.savez_compressed(OUT/f'case-{index}-{tol}.npz',x=sol.t,scaled_variables=sol.y)
            cases.append(dict(tolerance=tol,maximum_relative_difference=max(errors),components=errors,passed=passed))
        mass=mp.mpf(mu)+mp.mpf(q)**2*expected[0]
        phi=mp.mpf(q)*expected[2];charge=-mp.mpf(R)*mp.mpf(q)*expected[1]
        rows.append(dict(surface_q=q,exact_scaled_values=[mp.nstr(v,70) for v in expected],cases=cases,
            exact_mass_increment_geom_m=mp.nstr(mp.mpf(R)*mp.mpf(q)**2*expected[0],70),
            exact_Einstein_ADM_mass_geom_m=mp.nstr(mp.mpf(R)*mass,70),
            exact_phi_infinity=mp.nstr(phi,70),exact_scalar_tail_m=mp.nstr(charge,70),
            exact_exterior_charge_over_mass=mp.nstr(-charge/(mp.mpf(R)*mass),70),
            passed=all(r['passed'] for r in cases)))
        print('NONLINEAR VACUUM',q,max(r['maximum_relative_difference'] for r in cases),rows[-1]['passed'],flush=True)
    for magnitude in [1e-8,1e-4,.1]:
        assert exact(mu,magnitude)==exact(mu,-magnitude)
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_passed=all(r['passed'] for r in rows),symbolic_passed=True,exact_even_parity_controls=True,
        manufactured_boundary_data=True,self_consistent_nonlinear_stellar_interior=False,
        physical_EOS_certified=False,physical_atmosphere_matched=False,full_GR_evolution=False,
        observational_mass_or_background_inferred=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and r['symbolic_passed']
    print('PASS nonlinear exterior source and compactification controls; prescribed vacuum data only',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
