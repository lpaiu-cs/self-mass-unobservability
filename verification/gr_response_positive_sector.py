"""Uniform positive-real response disks from the positive kernel coefficients."""
from fractions import Fraction as F
import json,sys
import numpy as np
import sympy as sp
from mpmath import iv
import gr_response_complex_runner as previous

g=previous.g;cusp=g.cusp;ROOT=g.ROOT;OUT=g.OUT.parent/'gr-response-positive-sector'
I=g.I;low=g.low;high=g.high
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    previous.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_response_positive_sector.py',previous.OUT/'manifest.json',previous.OUT/'symbolic.json',cusp.OUT/'inputs.npz']
    save('plan.json',dict(classification='Proven',checkpoint='84029ef',bits=128,relative_Q_radius='1/256',
        beta_box=['3e-6','0.006'],eta_box=[-17,24],thermal_cutoff_offset=128,
        bindings={p.relative_to(ROOT).as_posix():g.sha(p) for p in files},
        target='Prove Re S(Q)>0 on |Q-q0|<=q0/256 for every q0>0 throughout the complete beta/eta box, and hence a zero-free dielectric denominator for every positive B. Avoid shrinking the radius by a Cauchy perturbation bound.',
        method='Use the exact positive-coefficient kernel in u=p/Q to bound both its modulus and phase. Retain positivity up to t_a=max(eta,0)+128; compare its complement to a doubled-temperature eta=0 distribution. Uniform all-Q lower and upper response bounds both scale as 1/max(1,a^2).',
        boundary='This is a conditional ideal-electron continuum theorem and an input-domain audit. It does not evaluate the outer integral or certify a complete interacting EOS/GR/observational model.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.sha(ROOT/rel)==digest,rel
    return p


def run():
    p=bindings();iv.prec=p['bits'];r=I(p['relative_Q_radius']);bmin,bmax=map(I,p['beta_box']);emin,emax=map(I,p['eta_box'])
    kappa=(1+r)/iv.sqrt(1-2*r);theta=iv.atan2(r,1-r);T=emax+p['thermal_cutoff_offset'];phase=2*T*(1+r)*r/(1-2*r)
    assert high(3*theta)<low(phase) and high(phase)<low(iv.pi/2)
    lower_factor=iv.cos(phase)*iv.cos(2*theta)/kappa*iv.exp(-(kappa*kappa-1)*T);assert low(lower_factor)>0
    # Common lower moment on p in [sqrt(beta_min),2sqrt(beta_min)].
    W=7*bmin*iv.sqrt(bmin)/(3*iv.sqrt(1+4*bmin))/(1+iv.exp(2-emin))
    # All upper moments below are for eta=0 at b=2*beta_max.
    b=2*bmax;e=iv.exp(-I(1))
    half=[I(2)/3+e*iv.sqrt(2),I(2)/5+e*iv.sqrt(10),I(2)/7+e*iv.sqrt(80)]
    integer=[I(1)/2+2*e,I(1)/3+5*e,I(1)/4+16*e]
    J2=b*iv.sqrt(2*b)*(half[0]+2*b*half[1]+b*b*half[2])+b*b*(integer[0]+2*b*integer[1]+b*b*integer[2])
    power=lambda s:I(s)**I(s)
    P3=b*iv.sqrt(2*b)*(2*power(F(3,2))+3*b*power(F(5,2))+b*b*power(F(7,2)))+b*b*(2*power(2)+3*b*power(3)+b*b*power(4))
    upper_Q2=I(32)/3*J2+48*P3
    upper_Q0=g.real_majorant(bmax,I(0)) # Includes an extra positive 2W0; harmless enlargement.
    C=I(max(high(upper_Q0),high(upper_Q2)))
    tail_factor=2*iv.exp(emax/2-I(p['thermal_cutoff_offset'])/2)
    tail=tail_factor*C;positive=lower_factor*W-(lower_factor+I(5)/3*kappa**3)*tail
    assert low(positive)>0
    u,Q,h,j,Gamma=sp.symbols('u Q h j Gamma',positive=True)
    raw=Q/Gamma*(Q*Q*u*u+u*(1+Q*Q*(u*u-1))*u*h)
    expected=Q/Gamma*u*u*(h+Q*Q*(1-(1-u*u)*h));assert sp.expand(raw-expected)==0
    t=sp.symbols('t',positive=True);bb=sp.symbols('b',positive=True)
    assert sp.expand(t*(2+bb*t)*(1+bb*t)-(2*t+3*bb*t*t+bb*bb*t**3))==0
    data=dict(np.load(cusp.OUT/'inputs.npz'));covered=0
    for beta,text in zip(data['beta'],data['root_intervals']):
        lo,hi=cusp.endpoints(str(text));assert F(p['beta_box'][0])<=F.from_float(float(beta))<=F(p['beta_box'][1]) and p['eta_box'][0]<=lo<=hi<=p['eta_box'][1];covered+=1
    save('result.json',dict(classification='Proven',passed=True,covered_actual_states=covered,
        relative_Q_radius=p['relative_Q_radius'],phase_upper=str(high(phase)),maximum_Q_angle_upper=str(high(theta)),
        positive_inner_factor_lower=str(low(lower_factor)),uniform_W0_lower=str(low(W)),uniform_real_response_upper_coefficient=str(high(C)),
        uniform_tail_upper_coefficient=str(high(tail)),positive_real_response_lower_coefficient=str(low(positive)),
        display_only=dict(phase_upper=float(high(phase)),positive_coefficient=float(low(positive)),tail_coefficient=float(high(tail))),
        positive_kernel='For 0<u<1, K_Q du=(Q/gamma)*[a(u)+Q^2*b(u)]du with a=u^2*h, b=u^2*j, h=atanh(u)/u>0 and j=1-(1-u^2)*h>0. For u>1, a=u*atanh(1/u)>0 and b=u^2+u*(u^2-1)*atanh(1/u)>0. Thus |K_Q|<=kappa^3*K_a and |K_Q|>=cos(2theta)/kappa*K_a. For Im Q>=0 its phase lies in [0,3theta]; conjugation handles Im Q<0.',
        positive_inner='Set a=q0*sqrt(1-2r), t_a=(sqrt(1+a^2*u^2)-1)/beta. Until t_a=max(eta,0)+128 the occupation phase lies in [-Phi,0], Phi bounded by the saved phase_upper<pi/2. Hence Re(K_Q*q)>=L*K_a*f(t_a-eta), with L=cos(Phi)*cos(2theta)/kappa*exp(-(kappa^2-1)*(24+128)).',
        uniform_lower='The positive kernel gives K_p>=p^2/(gamma*max(1,a^2)) on both p<a and p>=a. On p in [sqrt(beta_min),2sqrt(beta_min)], t<=2 and gamma<=sqrt(1+4beta_min), so S(a)>=W0_lower/max(1,a^2) uniformly over the full parameter box.',
        uniform_upper='The earlier bound gives S(a)<=U0. For u<=1/2 use h<=1/(1-u^2), j<=(2/3)*u^2/(1-u^2) to obtain a^2*K_p<=(4/3)*p^2*gamma. Outside the cusp band use a^2*K_p<=(32/3)*p^2*gamma; inside it the log integral gives 48*sup(p^3*gamma*f). Thus a^2*S(a)<=U2=(32/3)*integral p^2*gamma*f dp+48*sup(p^3*gamma*f), and S(a)<=max(U0,U2)/max(1,a^2).',
        upper_moments='For eta=0,b=2beta_max, use the same split at t=1. The half-power tail moments are bounded by sqrt(E[t^(2m)]*E[t^(2m+1)]); integer moments are incomplete exponential polynomials. Expand p^2*gamma dp <= b*sqrt(2b)*(t^1/2+2b*t^3/2+b^2*t^5/2)dt+b^2*(t+2b*t^2+b^2*t^3)dt. For the supremum expand p^3*gamma <= b*sqrt(2b)*(2t^3/2+3b*t^5/2+b^2*t^7/2)+b^2*(2t^2+3b*t^3+b^2*t^4), and use sup(t^s*f(t))<=s^s.',
        tail='For t_a>=max(eta,0)+128, f(t_a-eta)<=2*exp(eta-(max(eta,0)+128)/2)*f(t_a/2)<=2exp(-52)*f(t_a/2). The real tail is <=Ttail/max(1,a^2), and the complex tail is <=(5/3)*kappa^3 times this. The same positive kernel bound now has no extra A term.',
        final_bound='Re S(Q) >= [L*W0_lower-(L+(5/3)*kappa^3)*Ttail]/max(1,a^2) >0. Also Re(Q^2/B)>=q0^2*(1-2r)/B>0 for every B>0. Therefore Q^2/B+S(Q) never vanishes on these disks, for every real q0>0 and throughout the declared continuum parameter box.',
        normalized_denominator='Divide the raw positive S coefficient by the exact positive Sref. The lower bound is q0^2*(1-2r)/(Bmax*Sref)+Cpositive/(Sref*max(1,q0^2*(1-2r))). Combine with the earlier six complex majorants and Cauchy/Gauss for the outer integrand.',
        actual_outer_integral_certified=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():g.sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['covered_actual_states']==3206 and F(r['positive_real_response_lower_coefficient'])>0
    assert F(r['relative_Q_radius'])==F(1,256)
    print('PASS positive-real response and zero-free dielectric disks for all Q>0 on the continuum beta/eta box; actual outer integral remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
