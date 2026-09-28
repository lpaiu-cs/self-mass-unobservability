"""Constructive complex-Q response and zero-free dielectric neighborhoods."""
from fractions import Fraction as F
import json,sys
import numpy as np
from mpmath import iv
import sympy as sp
import gr_response_momentum_tail as tail

cusp=tail.cusp;ROOT=tail.ROOT;OUT=tail.OUT.parent/'gr-response-complex-domain'
I=cusp.I;low=cusp.low;high=cusp.high;sha=tail.sha
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    tail.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_response_complex_domain.py',cusp.OUT/'manifest.json',cusp.OUT/'inputs.npz',cusp.THERMO/'states.npz',tail.OUT/'manifest.json']
    save('plan.json',dict(classification='Proven',checkpoint='6292a62',bits=128,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},relative_Q_radius='1/128',eta_upper_limit=24,
        parameter_order=cusp.PARAMETERS,
        target='A complex analytic response disk about every real Q>0, uniform over the declared beta>0 and eta<=24. Explicit conservative six-field majorants at all 3206 actual states; constructive zero-free dielectric neighborhoods at the 12 saved Q values.',
        no_finite_response_assumption='Use only S(Q)>=0 on the real axis and X=Q^2/(B*Sref)>0. The zero-free radius is conservative and can later be enlarged by a certified positive S(Q) lower bound.',
        boundary='Analytic derivative bounds and zero-free neighborhoods, not completed outer integration. Very small conservative radii can require impractical panel counts. Full physical EOS and actual GR/observational closure remain open.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def symbolic():
    plan=bindings();iv.prec=plan['bits'];r=I(plan['relative_Q_radius']);phase=2*(24+iv.log(4))*(1+r)*r/(1-2*r)
    assert high(phase)<1 and high(1/iv.cos(I(1)/2))<F(4,3)
    middle=I(27)/4+I(15)/4*iv.log(3)+I(5)/4*iv.log(2);assert high(middle)<12
    u=sp.symbols('u',positive=True)
    assert sp.simplify(1+(u+1/u)*u/(1-u*u)-2/(1-u*u))==0
    z=sp.symbols('z',positive=True);assert sp.expand((1+z)**2-(1-z)**2-4*z)==0
    # Exact nonlinear derivatives of the same self integrand.
    x,s,sa,sb,sab=sp.symbols('x s sa sb sab',positive=True)
    first=sp.diff(s/(x+s),s);second=sp.diff(first,s)
    assert sp.simplify(first-x/(x+s)**2)==0 and sp.simplify(second+2*x/(x+s)**3)==0
    save('symbolic.json',dict(classification='Proven',passed=True,uniform_phase_upper=str(high(phase)),middle_integral_upper=str(high(middle)),
        analytic_disk='For |Q-q0|<=r*q0, a=q0*sqrt(1-2r), kappa=(1+r)/sqrt(1-2r). With p=Q*u and g=sqrt(1+a^2*u^2), Re gamma>=g, Re t>=t_a=(g-1)/beta and |t|<=kappa^2*t_a. For t_a<=max(eta,0)+log4 the phase is <=2*(24+log4)*(1+r)*r/(1-2r)<1. Outside it, |exp(t-eta)|>=4. There are no Fermi poles or square-root branch points in the disk.',
        complex_occupation='On the disk: |q|<=(5/3)*f(t_a-eta), |1-q|<=4/3, |1-2q|<=5/3. In the inner region use |1+E|>=cos(1/2)*(1+|E|); in the outer region use |E|-1 and |E|>=4.',
        complex_kernel='For A=Q^3*u^2/gamma, B=Q*u*(1+Q^2*(u^2-1))/(2*gamma), L=log((1+u)/abs(1-u)), the absolute kernel is <=kappa^3*(K_a+2*A_a). For u>=1 it is <=kappa^3*K_a; for u<1, the extra term is a^3*u*(1-u^2)*L/gamma_a<=2*A_a, since atanh(u)<=u/(1-u^2). K_a du is the positive real response kernel at a.',
        thermal_majorant='For d0=max(eta,0), f(t-eta)*t^j <= Cj*f(t/2-eta), C0=1, Cj=[2*(d0+log2+j)]^j. The ratio is <=min(1,2*exp(d0-t/2)); split at t=2*(d0+log2+j). Thus only a value-response bound at doubled beta is needed.',
        uniform_real_kernel='At any real a>0, |K_p|<=gamma*H(p/a), H(u)=1+(u+1/u)*atanh(min(u,1/u)). Outside [1/2,2], H<=8/3. Inside, integral H(u)/u du <=27/4+(15/4)log3+(5/4)log2<12. Therefore S(a)<=8Jgamma/3+12*sup(p*gamma*f). This deliberately overcounts the middle region in Jgamma.',
        thermal_integrals='Put b=2*beta, d=d0+1, w=min(1,exp(eta)), e=exp(-1). Bounds on integral t^s*f(t-eta)dt for s=-1/2,1/2,1,3/2 are w*[2sqrt(d)+e/sqrt(d)], w*[(2/3)d^(3/2)+e*sqrt(d+1)], w*[d^2/2+e*(d+1)], w*[(2/5)d^(5/2)+e*sqrt((d+1)*(d^2+2d+2))]. Tail square-root bounds use Jensen and Cauchy-Schwarz for d+Exp(1).',
        explicit_moments='Jgamma<=sqrt(b/2)*(J_-1/2+2b*J_1/2+b^2*J_3/2); W0<=b*sqrt(2b)*J_1/2+b^2*J_1. sup(p*gamma*f)<=w*[sqrt(2b)*((d0+1/2)^(1/2)+b*(d0+3/2)^(3/2))+b*((d0+1)+b*(d0+2)^2)]. Each power bound follows by splitting t^s*f at d0+s.',
        six_majorants='Let U=(8Jgamma/3+12Pmax+2W0)/Sref and Hj=(5/3)*kappa^(3+2j)*Cj*U. Then [H0,(4/3)H0,(20/9)H0,(4/3)H1,(20/9)H1,(20/9)H2+(4/3)H1] bounds the six normalized response fields throughout the complex disk.',
        zero_free='D0_lower=X+S_lower, X=q0^2/(B*Sref), with S_lower=0 allowed. Choose dyadic sigma<=min(1/4,D0_lower/(8H0)). Cauchy gives |S(Q)-S(q0)|<=H0*sigma/(1-sigma). The X variation is <=X*(2*sigma*r+(sigma*r)^2). Their sum is <D0_lower/2. Thus the denominator has modulus >=D0_lower/2 throughout |Q-q0|<=sigma*r*q0.',
        outer_gauss='On that inner disk bound G=S/(X+S) and its six eta/tau fields by the quotient identities and the denominator lower bound. Cauchy and the certified n-point Gauss norm give a panel remainder <=2h*norm_n*(2h/(Rinner-h))^(2n)*M_G. Actual outer quadrature and endpoint tails must still be evaluated and joined.'))


def real_majorant(beta,eta):
    b=2*beta;d0=I(max(high(eta),F(0)));d=d0+1;w=I(1) if high(eta)>=0 else iv.exp(I(high(eta)));e=iv.exp(-I(1))
    jm=w*(2*iv.sqrt(d)+e/iv.sqrt(d));jh=w*(I(2)/3*d*iv.sqrt(d)+e*iv.sqrt(d+1))
    j1=w*(d*d/2+e*(d+1));j3=w*(I(2)/5*d*d*iv.sqrt(d)+e*iv.sqrt((d+1)*(d*d+2*d+2)))
    J=iv.sqrt(b/2)*(jm+2*b*jh+b*b*j3);W=b*iv.sqrt(2*b)*jh+b*b*j1
    P=w*(iv.sqrt(2*b)*(iv.sqrt(d0+I(1)/2)+b*(d0+I(3)/2)**I(F(3,2)))+b*(d0+1+b*(d0+2)**2))
    return I(8)/3*J+12*P+2*W


def run():
    plan=bindings();symbolic();iv.prec=plan['bits'];data=dict(np.load(cusp.OUT/'inputs.npz'));stored=dict(np.load(cusp.THERMO/'states.npz'));B=I(float(stored['B']))
    r=I(plan['relative_Q_radius']);kappa=(1+r)/iv.sqrt(1-2*r);rows=[];minimum=None;maximum=F(0)
    for i,cell in enumerate(data['cells']):
        beta=I(data['beta'][i]);eta=cusp.interval(str(data['root_intervals'][i]));assert high(eta)<=plan['eta_upper_limit'] and low(beta)>0
        d0=I(max(high(eta),F(0)));U=real_majorant(beta,eta)/I(data['Sref'][i]);H=[I(5)/3*kappa**(3+2*j)*([I(1),2*(d0+iv.log(2)+1),(2*(d0+iv.log(2)+2))**2][j])*U for j in range(3)]
        M=[H[0],I(4)/3*H[0],I(20)/9*H[0],I(4)/3*H[1],I(20)/9*H[1],I(20)/9*H[2]+I(4)/3*H[1]];neighborhoods=[]
        for j,Q0 in enumerate(data['Q'][i]):
            Q=I(Q0);X=Q*Q/(B*I(data['Sref'][i]));D=I(low(X));sigma=F(1,4)
            while sigma>low(D/(8*H[0])):sigma/=2
            s=I(sigma);delta=X*(2*s*r+(s*r)**2)+H[0]*s/(1-s);assert high(delta)<low(D/2)
            den=D/2;Xmax=X*(1+s*r)**2
            G=[M[0]/den,Xmax*M[1]/den**2,Xmax*(M[2]/den**2+2*M[1]**2/den**3),Xmax*M[3]/den**2,
                Xmax*(M[4]/den**2+2*M[1]*M[3]/den**3),Xmax*(M[5]/den**2+2*M[3]**2/den**3)]
            relative=sigma*F(plan['relative_Q_radius']);minimum=relative if minimum is None else min(minimum,relative);maximum=max(maximum,high(delta/D))
            neighborhoods.append(dict(z_index=j,sigma_exact=str(sigma),relative_zero_free_radius=str(relative),radius_lower=str(low(s*r*Q)),
                denominator_lower=str(low(den)),perturbation_upper=str(high(delta)),G_majorants=list(map(lambda x:str(high(x)),G))))
        rows.append(dict(cell=int(cell),response_majorants=list(map(lambda x:str(high(x)),M)),neighborhoods=neighborhoods))
    save('result.json',dict(classification='Proven',passed=True,cells=len(rows),cases=sum(len(r['neighborhoods']) for r in rows),
        minimum_relative_zero_free_radius=str(minimum),maximum_relative_denominator_perturbation=str(maximum),rows=rows,
        display_only=dict(minimum_relative_radius=float(minimum),maximum_relative_perturbation=float(maximum)),actual_outer_integral_certified=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    s=json.loads((OUT/'symbolic.json').read_text());r=json.loads((OUT/'result.json').read_text());assert s['passed'] and r['passed'] and r['cells']==3206 and r['cases']==38472
    assert F(s['uniform_phase_upper'])<1 and F(s['middle_integral_upper'])<12
    for row in r['rows']:
        assert len(row['response_majorants'])==6 and min(map(F,row['response_majorants']))>0
        for n in row['neighborhoods']:assert 0<F(n['perturbation_upper'])<F(n['denominator_lower']) and F(n['radius_lower'])>0 and min(map(F,n['G_majorants']))>0
    print('PASS uniform complex-Q response majorants and constructive zero-free dielectric disks; outer evaluation and physical closure remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
