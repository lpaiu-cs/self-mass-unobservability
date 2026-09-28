"""Native EOS line integrals for coupled density/temperature material updates."""
import json,sys
import numpy as np
import sympy as s
import mpmath as mp
import direct_eos_gr as g

OUT=g.OUT/'gr-material-thermo-increment'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def increment(eos,lr,lt,X,eta,theta,number):
    knots,weights=np.polynomial.legendre.leggauss(number);ld=np.longdouble
    eta=ld(eta);theta=ld(theta);terms=[]
    for knot in (knots+1)/2:
        r=float(ld(lr)+ld(knot)*eta);t=float(ld(lt)+ld(knot)*theta)
        a=eos(2,r,t,X).astype(ld);rho,P=a[:2];temperature=np.exp(ld(t))
        terms.append([eta*a[9]+theta*a[10],P*(eta*a[5]+theta*a[6]),
            (eta*(a[9]-P/rho)+theta*a[10])/temperature])
    return np.asarray(weights,dtype=ld)@np.asarray(terms,dtype=ld)/2


def run():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_material_thermo_increment.py',g.OUT/'initial-state-17-4.npz',
        g.OUT/'gr-microphysics/auxiliaries.npz',g.OUT/'gr-material-rest-energy/manifest.json',
        g.OUT/'gr-caloric-increment/plan.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d8689de',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        cells=[0,1,2,64,1024,2048,2972,4096,5734],
        probes=[[1e-4,2e-4],[-1e-4,-2e-4],[1e-4,-2e-4],[1e-25,2e-25],[-1e-25,2e-25],[0.,0.]],
        quadrature_counts=[2,4],relative_quadrature_tolerance=1e-8,relative_endpoint_tolerance=1e-8,
        endpoint_energy_roundoff='Sum max(2 erg/g,32 ulp(abs(u))) at the two native endpoints. This finite acceptance budget is not a rigorous native evaluation error bound.',
        known_cut_probe='Also evaluate a constructed log-temperature path from log(1e6)-1e-6 to log(1e6)+1e-6 at the profile density/composition nearest that temperature. A known native branch boundary is not silently excluded or certified smooth.',
        scope='Extend the fixed-density caloric increment to density and temperature together. Direct native EOS values at line quadrature points supply du, dP and ds. Conditional line-integral identities require the EOS derivatives and first law along the path; native branch jumps require additional jump terms. Finite comparisons do not supply derivative, root, continuum or physical EOS certificates, and are not GR evolution.'))
    eta,theta,u0,du,rho=s.symbols('eta theta u0 du rho',real=True)
    assert s.expand(rho*s.exp(eta)*(u0+du)-rho*u0-rho*(u0*(s.exp(eta)-1)+s.exp(eta)*du))==0
    z=s.symbols('z',real=True);lr=s.Function('lr')(z);lt=s.Function('lt')(z);u=s.Function('u')(lr,lt)
    assert s.simplify(s.diff(u,z).subs({s.diff(lr,z):eta,s.diff(lt,z):theta})-
        eta*s.diff(u,lr)-theta*s.diff(u,lt))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        path='lnrho(z)=lnrho0+eta*z, lnT(z)=lnT0+theta*z, 0<=z<=1, fixed composition.',
        increments='Delta u=integral[eta*u_lnrho+theta*cvT]dz; Delta P=integral P*(eta*chi_rho+theta*chi_T)dz. If the Gibbs relation holds, Delta s=integral[eta*(u_lnrho-P/rho)+theta*cvT]/T dz.',
        stable_density_energy='Delta(rho*u)=rho0*[u0*expm1(eta)+exp(eta)*Delta u]. Rest mass and changing metric are handled by the separate material-rest-energy identity; this is the internal-energy density term.',
        conditional_error='Four-point Gauss remainder on [0,1] has bound max|f^(8)|/1778112000 for each directional integrand f. Its existence/magnitude and native node errors are not provided by a 2/4 comparison. If a path crosses a branch jump, add each one-sided primitive jump; differentiating only smooth pieces need not reproduce the endpoint.',
        no_full_GR_evolution=True))
    def manufactured(mode,r,t,X):
        rho=np.exp(r);T=np.exp(t);R=2.;cv=3.;rad=.01;P=R*rho*T+rad*T**4/3
        a=np.zeros(21);a[:4]=[rho,P,cv*T+rad*T**4/rho,cv*t-R*r+4*rad*T**3/(3*rho)]
        a[5]=R*rho*T/P;a[6]=(R*rho*T+4*rad*T**4/3)/P
        a[9]=-rad*T**4/rho;a[10]=cv*T+4*rad*T**4/rho
        return a
    mp.mp.dps=80;controls=[]
    for e,t in [(.02,-.01),(1e-25,2e-25)]:
        value=increment(manufactured,0.,0.,None,e,t,4)
        E=mp.mpf(e);T=mp.mpf(t);rad=mp.mpf(float(.01))
        truth=[3*mp.expm1(T)+rad*mp.expm1(4*T-E),
            2*mp.expm1(E+T)+rad/3*mp.expm1(4*T),
            3*T-2*E+4*rad/3*mp.expm1(3*T-E)]
        errors=[]
        for actual,expected in zip(value,truth):
            num,den=actual.as_integer_ratio();errors.append(float(abs((mp.mpf(num)/den)/expected-1)))
        assert max(errors)<1e-12,errors
        controls.append(dict(eta=e,theta=t,relative_errors=errors,passed=True))
    save('analytic-controls.json',dict(classification='Counterexample candidate',rows=controls,
        model='Ideal gas plus radiation with an analytic Gibbs-consistent energy/pressure/entropy. Large and sub-ULP changes compared at 80 decimal digits; high-precision comparison is not an interval proof.'))
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));eos=g.EOS();plan=json.loads((OUT/'plan.json').read_text())
    cases=[(i,float(state['lnd'][i]),float(state['lnT'][i]),e,t,'profile probe')
        for i in plan['cells'] for e,t in plan['probes']]
    cut=int(np.argmin(abs(state['lnT']-np.log(1e6))))
    cases.append((cut,float(state['lnd'][cut]),float(np.log(1e6)-1e-6),0.,2e-6,'known temperature cut'))
    rows=[];saved=[]
    for number,(i,lr,lt,eta,theta,label) in enumerate(cases):
        X=state['X'][i];a=eos(2,lr,lt,X).astype(np.longdouble)
        end=eos(2,float(np.longdouble(lr)+np.longdouble(eta)),float(np.longdouble(lt)+np.longdouble(theta)),X).astype(np.longdouble)
        values=[increment(eos,lr,lt,X,eta,theta,k) for k in plan['quadrature_counts']]
        scale=np.array([abs(eta*a[9])+abs(theta*a[10]),a[1]*(abs(eta*a[5])+abs(theta*a[6])),
            (abs(eta*(a[9]-a[1]/a[0]))+abs(theta*a[10]))/np.exp(np.longdouble(lt))])
        scale=np.maximum(scale,np.longdouble('1e-300'))
        energy_budget=sum(max(2.,32*np.spacing(float(abs(x[2])))) for x in [a,end])
        budgets=np.array([energy_budget,64*(np.spacing(float(abs(a[1])))+np.spacing(float(abs(end[1])))),
            energy_budget/np.exp(np.longdouble(lt))+64*(np.spacing(float(abs(a[3])))+np.spacing(float(abs(end[3]))))])
        endpoint=end[1:4][[1,0,2]]-a[1:4][[1,0,2]]
        quadrature_error=abs(values[1]-values[0])/scale
        endpoint_error=abs(values[1]-endpoint)/(budgets+plan['relative_endpoint_tolerance']*scale)
        row=dict(case=number,cell=i,eta=eta,theta=theta,kind=label,
            quadrature_scaled_errors=list(map(float,quadrature_error)),endpoint_budget_scores=list(map(float,endpoint_error)),
            native_endpoint_energy_unchanged=bool(end[2]==a[2]),nonzero_integrated_energy=bool(values[1][0]!=0))
        row['passed']=bool(np.all(quadrature_error<=plan['relative_quadrature_tolerance']) and np.all(endpoint_error<=1))
        rows.append(row);saved.append([*values[0],*values[1],*endpoint,*scale,*budgets])
        if number%6==0:print('MATERIAL THERMO',row,flush=True)
    np.savez_compressed(OUT/'increments.npz',values=np.array(saved,dtype=np.longdouble))
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_passed=all(r['passed'] for r in rows),failed_cases=[r for r in rows if not r['passed']],
        native_smoothness_certified=False,native_derivative_error_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS material thermodynamic increment bindings; inspect failed_cases for numerical gates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
