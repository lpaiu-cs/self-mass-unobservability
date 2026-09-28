"""Remove the logarithmic cusp and certify an additive analytic Q strip."""
from fractions import Fraction as F
from functools import lru_cache
import json,math,sys
import numpy as np
import mpmath as mp
import sympy as sp
from mpmath import iv
import gr_response_positive_sector as sector
import gr_highq_moments as moments

g=sector.g;cusp=moments.cusp;ROOT=cusp.ROOT;OUT=sector.OUT.parent/'gr-response-regular-defined';I=cusp.I;low=g.low;high=g.high

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    sector.verify();assert not OUT.exists();OUT.mkdir()
    smooth=OUT.parent/'gr-smooth-response-defined'
    files=[ROOT/'verification/gr_response_regular_identity.py',sector.OUT/'manifest.json',sector.OUT/'result.json',cusp.OUT/'inputs.npz',
        smooth/'manifest.json',smooth/'complete-response.jsonl.gz',ROOT/'verification/gr_polarization_thermodynamics.py']
    save('plan.json',dict(classification='Proven',checkpoint='a87e0c5',bits=128,strip_in_sqrt_beta='1/16',thermal_split=176,
        beta_box=['3e-6','0.006'],eta_box=[-17,24],bindings={p.relative_to(ROOT).as_posix():g.sha(p) for p in files},
        control_positions=[0,1603,3205],control_z_indices=[0,5,8,11],control_digits=64,
        target='Prove an exactly equivalent positive divided-difference integral with a removable diagonal, then prove analytic Re S>0 for all real Re Q and |Im Q|<=sqrt(beta)/16. Compare its six real response jets against the already certified complete momentum integrals.',
        boundary='The additive strip applies to S at real beta and eta. It does not by itself exclude zeros of Q^2/B+S near the imaginary screening poles, certify joint eta/tau bounds, or evaluate the remaining outer self integral.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.sha(ROOT/rel)==digest,rel
    return p


def prove():
    p=bindings();iv.prec=p['bits'];d=I(p['strip_in_sqrt_beta']);bmax=I(p['beta_box'][1]);sigma=d*iv.sqrt(bmax)
    kappa=1/iv.sqrt(1-sigma*sigma);plus=iv.sqrt(1+sigma*sigma);D=d*d*kappa/2
    phase=d*kappa*iv.sqrt(2*p['thermal_split']);gamma_phase=iv.atan2(sigma*kappa*kappa,I(1));angle=phase+gamma_phase
    assert high(angle)<low(iv.pi/2)
    L=iv.cos(angle)/plus;old=json.loads((sector.OUT/'result.json').read_text());W=I(old['uniform_W0_lower']);C=I(old['uniform_real_response_upper_coefficient'])
    tail=iv.exp(I(p['eta_box'][1])-I(p['thermal_split'])/2)*C
    lower=L*W-(2*L+4*kappa*iv.exp(D))*tail;assert low(lower)>0
    P,Q,V0,V1,f,gamma,H=sp.symbols('p Q V0 V1 f gamma H',positive=True)
    V0p=-P*gamma*f;V1p=-P/gamma*f
    Wp=V0p-2*P*V1-P*P*V1p
    assert sp.cancel((Wp+P*(f/gamma+2*V1)).subs(gamma**2,1+P*P))==0
    assert sp.cancel((V0-Q*Q*V1)/(P*P-Q*Q)-V1-(V0-P*P*V1)/(P*P-Q*Q))==0
    assert sp.diff(sp.log(P+Q)-sp.log(P-Q),P)==1/(P+Q)-1/(P-Q)
    assert sp.cancel(sp.diff(sp.log(P+Q)-sp.log(P-Q),P)+2*Q/(P*P-Q*Q))==0
    x=sp.symbols('x',positive=True);assert sp.cancel(2*x-(2*x+x*x)/(1+x)**2-x*x*(2*x+3)/(x+1)**2)==0
    save('result.json',dict(classification='Proven',passed=True,phase_upper=str(high(angle)),positive_response_lower_coefficient=str(low(lower)),
        real_energy_shift_upper=str(high(D)),tail_coefficient_upper=str(high(tail)),display_only=dict(phase_upper=float(high(angle)),positive_coefficient=float(low(lower))),
        primitives='For gamma_p=sqrt(1+p^2), f=1/(1+exp((gamma_p-1)/beta-eta)), define V1(p)=integral_gamma_p^infinity f(gamma) dgamma and V0(p)=integral_gamma_p^infinity gamma^2 f(gamma) dgamma. Both are even in p. W(p)=V0(p)-p^2 V1(p) and H(p)=f(p)/gamma_p+2V1(p) satisfy W_prime(p)=-p*H(p), H>0.',
        identity='Integrate the original logarithmic kernel by parts, in the principal-value sense: S=A0-PV integral [V0(p)-Q^2 V1(p)]/(p^2-Q^2) dp. Since integral V1(p)dp=A0=integral p^2 f(p)/gamma_p dp, the constant terms cancel exactly. PV integral_0^infinity 1/(p^2-Q^2)dp=0, hence S=integral_0^infinity [W(Q)-W(p)]/(p^2-Q^2) dp. The diagonal limit is H(Q)/2. No finite-temperature approximation is introduced.',
        positive_double_integral='Using W_prime=-pH gives S(Q)=1/2 integral_0^infinity dp integral_0^1 du H(sqrt((1-u)*Q^2+u*p^2)). The integrand is positive for real p,Q. At Q=0 its limit is integral H(p)dp, agreeing with the original zero-Q kernel. The diagonal p=Q and Q=p=0 are removable. The same identity holds for all five eta/tau partials by dominated differentiation.',
        analytic_definition='For complex Q=x+iy and lambda=1-u, put a^2=lambda*x^2+u*p^2 and gamma_a=sqrt(1+a^2). In the double integral use gamma=sqrt(1+lambda*Q^2+u*p^2) directly, so no auxiliary square-root branch in p is introduced. Re gamma lies in [gamma_a/kappa,gamma_a], |gamma|<=plus*gamma_a and |arg gamma|<=atan(sigma_max*kappa^2).',
        energy='Writing t=(gamma-1)/beta and t_a=(gamma_a-1)/beta, t_a-D<=Re t<=t_a, D=delta^2*kappa/2, delta=1/16. For t_a<=Tc=176, |Im t|<=delta*kappa*sqrt(2Tc). This follows from lambda*|x|<=a, Im gamma=lambda*x*y/Re gamma and a^2/gamma_a^2<=2beta*t_a. The saved phase includes the inverse-gamma factor and is below pi/2.',
        inner='For t_a<=Tc, |arg f(t-eta)|<=|Im t| and |f(t-eta)|>=f(t_a-eta). V1=beta*integral_0^infinity f(t+s-eta)ds has the same phase bound and real lower comparison. Therefore Re H_complex >= L*H_real with L=cos(total_phase)/plus.',
        tail='For t_a>Tc, Re(t-eta)>log 2 and |f|, |V1|/beta <=2exp(eta+D-t_a). Compare with eta=0 at twice beta: f(t_a/2)>=exp(-t_a/2)/2 and log(1+exp(-t_a/2))>=exp(-t_a/2)/2. Thus |H_complex|<=4kappa exp(eta+D-Tc/2)H_half, while H_real<=2exp(eta-Tc/2)H_half. Integrating these enlarged positive tails uses the previously certified C/max(1,x^2).',
        conclusion='For every real x, beta in the declared box and eta in the declared box, Re S(x+iy)>=[L*W0-(2L+4kappa*exp(D))*exp(24-Tc/2)*C]/max(1,x^2)>0 when |y|<=sqrt(beta)/16. Exponential Fermi tails and the positive double integral supply integrable domination, proving analyticity on the open strip and continuity on its boundary. Fermi poles and gamma branch points are excluded by the same inner/tail split.',
        limitation='Re S>0 does not exclude dielectric zeros when Re(Q^2)<0. No stronger G-domain, joint parameter domain or actual outer integral is claimed.'))


def controls():
    import gzip
    p=bindings();mp.mp.dps=p['control_digits'];data=dict(np.load(cusp.OUT/'inputs.npz'));checks=[]
    def M(x):
        v=F(x);return mp.mpf(v.numerator)/v.denominator
    wanted={12*i+j for i in p['control_positions'] for j in p['control_z_indices']};source={}
    with gzip.open(OUT.parent/'gr-smooth-response-defined/complete-response.jsonl.gz','rt') as stream:
        for line in stream:
            row=json.loads(line)
            case=12*row['position']+row['z_index']
            if case in wanted:source[case]=row
    for i in p['control_positions']:
        eta=sum(map(M,cusp.endpoints(str(data['root_intervals'][i]))))/2;beta=mp.mpf(float(data['beta'][i]));Sref=mp.mpf(float(data['Sref'][i]));unit=mp.sqrt(beta)
        @lru_cache(maxsize=30000)
        def Wjet(v):
            momentum=unit*v;gamma=mp.sqrt(1+momentum*momentum);t=momentum*momentum/(beta*(gamma+1));a=eta-t;z=mp.exp(a);q=z/(1+z)
            f0=mp.log1p(z);f1=-mp.polylog(2,-z);f2=-2*mp.polylog(3,-z);terms=[(1,beta,f0,q,q*(1-q)),(2,2*beta**2*gamma,f1,f0,q),(3,beta**3,f2,2*f1,2*f0)]
            out=[mp.mpf(0)]*6
            for n,c,value,first,second in terms:
                add=[value,first,second,n*value+t*first,n*first+t*second,n*n*value+(2*n-1)*t*first+t*t*second]
                out=[old+c*x for old,x in zip(out,add)]
            return out
        for j in p['control_z_indices']:
            Q=mp.mpf(float(data['Q'][i,j]));x=Q/unit;wq=Wjet(x);points=sorted(set([mp.mpf(0),mp.mpf(1),mp.mpf(2),mp.mpf(4),mp.mpf(8),mp.mpf(16),x]));points.append(mp.inf)
            values=[]
            for k in range(6):
                def integrand(v):
                    if v==x:
                        gamma=mp.sqrt(1+Q*Q);t=Q*Q/(beta*(gamma+1));q=1/(1+mp.exp(t-eta));z=1-q;L=mp.log1p(mp.exp(eta-t))
                        occ=[q,q*z,q*z*(1-2*q),t*q*z,t*q*z*(1-2*q),(t*t*(1-2*q)-t)*q*z][k]
                        primitive=[L,q,q*z,L+t*q,q+t*q*z,L+t*q+t*t*q*z][k]
                        return unit*(occ/gamma+2*beta*primitive)/(2*Sref)
                    return (wq[k]-Wjet(v)[k])/(unit*Sref*(v*v-x*x))
                values.append(mp.quad(integrand,points))
            row=source[12*i+j];texts=row['enclosures'];inside=all(M(a)<=value<=M(b) for (a,b),value in zip(map(cusp.endpoints,texts),values))
            checks.append(dict(cell=int(data['cells'][i]),z_index=j,passed=inside,reference=list(map(str,values))));print('REGULAR RESPONSE CONTROL',int(data['cells'][i]),j,inside,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=len(checks)==12 and all(x['passed'] for x in checks),checks=checks,
        scope='Independent polylogarithmic primitive and removable divided-difference integration at the exact certified-root midpoint. Inclusion is tested in the already certified full momentum-response intervals.'))


def finalize():
    bindings()
    for name in ['result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():g.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and F(r['positive_response_lower_coefficient'])>0
    controls=json.loads((OUT/'controls.json').read_text());assert controls['passed'] and len(controls['checks'])==12
    print('PASS positive regular response identity, additive analytic strip and 72 independent controls; outer self integral remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
