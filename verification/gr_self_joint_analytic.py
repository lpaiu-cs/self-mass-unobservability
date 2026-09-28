"""Uniform self-integrand derivative bounds from a joint complex polydisk."""
from fractions import Fraction as F
import json,math,sys
import sympy as sp
from mpmath import iv
import gr_response_positive_sector as sector

g=sector.g;ROOT=g.ROOT;OUT=sector.OUT.parent/'gr-self-joint-analytic';I=g.I;low=g.low;high=g.high
ORDERS=[(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    sector.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_self_joint_analytic.py',sector.OUT/'manifest.json',sector.OUT/'result.json',g.cusp.OUT/'inputs.npz']
    save('plan.json',dict(classification='Proven',checkpoint='a226168',bits=128,relative_Q_radius='1/256',eta_radius='1/100',log_temperature_radius='1/1000',
        beta_reference_box=['3e-6','0.006'],eta_reference_box=[-17,24],thermal_cutoff_offset=128,
        bindings={p.relative_to(ROOT).as_posix():g.sha(p) for p in files},parameter_orders=ORDERS,
        target='Prove joint analyticity and Re S>0 for |Q-q0|<=q0/256, |eta-eta0|<=0.01 and |tau|<=0.001, beta=beta0*exp(tau), throughout the real reference box. Bound G=S/(Q^2/B+S) by its phase geometry, then use eta/tau Cauchy bounds directly instead of dividing loose response majorants by a small denominator.',
        arithmetic_scope='Reuse the previously certified common W0 lower bound and the eta=0 thermal-response upper coefficient. Increasing its thermal beta by lambda>=1 raises the bound by at most lambda^4, since all its monomials have positive coefficients and beta exponents at most four.',
        boundary='Conditional analytic and derivative theorem for the same ideal-electron self integrand. Actual finite outer integration and full EOS/GR/observational closure remain open.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.sha(ROOT/rel)==digest,rel
    return p


def run():
    p=bindings();iv.prec=p['bits'];old=json.loads((sector.OUT/'result.json').read_text())
    r=I(p['relative_Q_radius']);he=I(p['eta_radius']);ht=I(p['log_temperature_radius'])
    kappa=I(high((1+r)/iv.sqrt(1-2*r)));theta=I(high(iv.atan2(r,1-r)));psi=2*theta+ht
    cl=I(low(iv.exp(-ht)*iv.cos(psi)/kappa));cu=I(high(iv.exp(ht)*kappa*kappa));assert 0<low(cl)<=high(cl)<1<low(cu)
    T=(I(p['eta_reference_box'][1])+he+p['thermal_cutoff_offset'])/cl
    phase=cu*T*iv.sin(psi)+he;product_phase=phase+3*theta;assert high(product_phase)<low(iv.pi/2)
    L=iv.cos(product_phase)*iv.cos(2*theta)/kappa*iv.exp(-(cu-1)*T);assert low(L)>0
    W=I(old['uniform_W0_lower'])*iv.exp(-he);upper=I(old['uniform_real_response_upper_coefficient'])/cl**4
    tail=2*iv.exp((I(p['eta_reference_box'][1])+he-p['thermal_cutoff_offset'])/2)*upper
    positive=L*W-(L+2*kappa**3)*tail;assert low(positive)>0
    quotient=1/iv.cos(2*theta);majorants=[math.factorial(n)*math.factorial(m)*quotient/(he**n*ht**m) for n,m in ORDERS]
    v,s=sp.symbols('v s',real=True);assert sp.expand(1+v*v-2*v*s-(v-s)**2-(1-s*s))==0
    E,Tau,H_E,H_T,C=sp.symbols('E Tau H_E H_T C',positive=True);controls=[]
    for n,m in ORDERS:
        monomial=C*(E/H_E)**n*(Tau/H_T)**m
        observed=sp.diff(monomial,E,n,Tau,m).subs({E:0,Tau:0});bound=sp.factorial(n)*sp.factorial(m)*C/H_E**n/H_T**m
        assert sp.simplify(observed-bound)==0;controls.append(dict(order=[n,m],saturating_derivative=str(bound),passed=True))
    save('result.json',dict(classification='Proven',passed=True,product_phase_upper=str(high(product_phase)),
        c_lower=str(low(cl)),c_upper=str(high(cu)),reference_moment_lower=str(low(W)),tail_upper_coefficient=str(high(tail)),
        positive_real_response_lower_coefficient=str(low(positive)),self_integrand_modulus_upper=str(high(quotient)),
        parameter_orders=ORDERS,derivative_majorants=list(map(lambda v:str(high(v)),majorants)),monomial_controls=controls,
        display_only=dict(product_phase_upper=float(high(product_phase)),positive_coefficient=float(low(positive)),modulus_upper=float(high(quotient)),derivative_majorants=list(map(lambda v:float(high(v)),majorants))),
        joint_domain='eta=eta0+delta_eta, beta=beta0*exp(tau), with the declared complex radii. With a=q0*sqrt(1-2r), t_a=(sqrt(1+a^2*u^2)-1)/beta0: c_lower*t_a<=Re t, |t|<=c_upper*t_a, and |arg t|<=2theta+h_tau. These follow from |gamma+1|<=kappa*(gamma_a+1), |Q|>=a and |Q|<=kappa*a.',
        inner='Until t_a=(max(Re eta,0)+128)/c_lower, |Im(t-eta)|<=c_upper*Tmax*sin(2theta+h_tau)+h_eta. Add 3theta for the kernel phase. The saved product phase is below pi/2. Also |q|>=f(c_upper*t_a-Re eta)>=exp(-(c_upper-1)*Tmax)*f(t_a-Re eta). Thus the inner real part is at least L times the positive real-reference integral.',
        tail='Outside that cutoff, |q|<=2exp(Re eta-c_lower*t_a)<=4exp((24+h_eta-128)/2)*f(c_lower*t_a/2). The real-reference tail is <=2exp((24+h_eta-128)/2)*f(t_a/2). Bound both by the eta=0 response at beta=2beta0/c_lower. Its coefficient is at most the old common coefficient/c_lower^4. The complex tail has the additional factor kappa^3.',
        positive='The real-reference lower moment with Re eta>=-17-h_eta is at least exp(-h_eta) times the former W0 bound. Therefore Re S(Q,eta,beta)>=[L*Wnew-(L+2*kappa^3)*Ttail]/max(1,a^2)>0, uniformly on the joint domain. Fermi poles, square-root branch points and dielectric zeros are excluded; exponential domination proves joint analyticity.',
        modulus='Let X=Q^2/B and S have positive real part. |arg X|<=2theta and |arg S|<pi/2. Hence |X+S|^2/|S|^2 >=1+v^2-2v*sin(2theta)=(v-sin(2theta))^2+cos(2theta)^2, v=|X|/|S|. Thus |G|<=sec(2theta), uniformly for every B>0.',
        Cauchy='For each Q in the full Q disk, the six eta/tau derivatives at the real reference state are bounded by n!*m!*sec(2theta)/(h_eta^n*h_tau^m). These are also uniform complex-Q bounds, so Q Cauchy/Gauss can be applied afterward. The six holomorphic monomials saturate the factorial/radius factors and provide exact symbolic controls.',
        actual_outer_integral_certified=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():g.sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and F(r['positive_real_response_lower_coefficient'])>0 and len(r['derivative_majorants'])==6
    assert all(c['passed'] for c in r['monomial_controls'])
    print('PASS joint Q/eta/log-temperature analyticity and direct uniform self-integrand Cauchy bounds; actual outer evaluation remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
