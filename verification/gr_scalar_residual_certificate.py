"""A posteriori exact-polynomial residual bound for the frozen scalar coefficient."""
from fractions import Fraction as Q
from math import comb
import json, sys
import mpmath as mp
import numpy as np
import sympy as sp
import gr_scalar_global_regular as regular

g=regular.g;OUT=g.OUT/'gr-scalar-residual-certificate'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();regular.verify()
    paths=[g.ROOT/'verification/gr_scalar_residual_certificate.py',regular.OUT/'manifest.json',
        regular.original.OUT/'field-2e-12.npz',regular.original.OUT/'result.json']
    save('plan.json',dict(classification='Proven',checkpoint='427ebf8',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        interpolant='Exact rational C1 cubic Hermite polynomial from the saved field values and derivatives. Scale by field(x1)/(1+b(0)*x1^2/(6*a(0))), set p(0)=1,p_prime(0)=0. This central estimate is only an approximation; its error is included by the global residual, never assumed exact.',
        residual='J(x)=int_0^x t^2*b(t)*p(t)dt is integrated as an exact polynomial. E_prime=p_prime-J/(x^2*a). Bernstein bounds of x^2*a*p_prime-J divided by a positive denominator bound give a rigorous L1 bound epsilon on E, with x^2 cancelled analytically in the centre piece.',
        propagation='The proven Volterra contraction gives ||phi-p|| <= epsilon/(1-kappa). Then |W(1)-J(1)| <= B2*epsilon/(1-kappa). Propagate both through the exact declared Schwarzschild normalization and R/M with outward intervals.',
        controls='An exact constant-potential sixth-order Taylor polynomial is compared with the analytic trigonometric response. Check arbitrary-degree Bernstein conversion and integral algebra symbolically. Report whether the old floating IVP coefficient is enclosed; do not change gates or interpolate away a disagreement.',
        boundary='Certificate for the exact stored coefficient model only, including its centre interval and infinite exterior. No physical EOS/GR coefficient-generation error, finite scalar backreaction, dynamical stability or observation certificate.'))


def add(a,b,sign=1):
    return [(a[i] if i<len(a) else Q(0))+sign*(b[i] if i<len(b) else Q(0)) for i in range(max(len(a),len(b)))]


def mul(a,b):
    result=[Q(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):result[i+j]+=x*y
    return result


def at(a,x):
    result=Q(0)
    for v in reversed(a):result=result*x+v
    return result


def bernstein(coeff,h):
    n=len(coeff)-1;powers=[v*h**j for j,v in enumerate(coeff)]
    return [sum((powers[j]*Q(comb(k,j),comb(n,j)) for j in range(k+1)),Q(0)) for k in range(n+1)]


def piece(left,h,a,b,p,flux):
    xsquared=[left*left,2*left,Q(1)]
    forcing=mul(xsquared,mul(b,p));J=[flux]+[v/Q(i+1) for i,v in enumerate(forcing)]
    derivative=[Q(i)*v for i,v in enumerate(p) if i]
    numerator=add(mul(xsquared,mul(a,derivative)),J,-1)
    amin=min(bernstein(a,h));assert amin>0
    if left==0:
        assert numerator[:2]==[0,0]
        error=h*max(map(abs,bernstein(numerator[2:],h)))/amin
    else:error=h*max(map(abs,bernstein(numerator,h)))/(left*left*amin)
    return at(J,h),error


def response(p1,J1,error,kappa,B2,L,scale):
    iv=regular.iv;delta=error/(1-kappa)
    field=iv(p1)+mp.iv.mpf([-1,1])*iv(delta)
    flux=iv(J1)+mp.iv.mpf([-1,1])*iv(B2*delta)
    normal=field+L*flux;assert normal.a>0
    return iv(scale)*flux/normal,delta,normal


def controls():
    q=Q(3,10);p=[Q(1),Q(0),-q*q/6,Q(0),q**4/120,Q(0),-q**6/5040]
    J,error=piece(Q(0),Q(1),[Q(1)],[-q*q],p,Q(0))
    interval,delta,normal=response(at(p,Q(1)),J,error,q*q/6,q*q/3,mp.iv.mpf(1),Q(1))
    exact=1-mp.iv.tan(regular.iv(q))/regular.iv(q)
    assert interval.a<exact.a and exact.b<interval.b
    z=sp.symbols('z');poly=[Q(2),Q(-3),Q(5),Q(0),Q(7),Q(-2),Q(3)];h=Q(2,9)
    B=bernstein(poly,h);n=len(poly)-1
    expr=sum(sp.Rational(v.numerator,v.denominator)*comb(n,j)*z**j*(1-z)**(n-j) for j,v in enumerate(B))
    target=sum(sp.Rational(v.numerator,v.denominator)*(sp.Rational(h.numerator,h.denominator)*z)**j for j,v in enumerate(poly))
    assert sp.expand(expr-target)==0
    return dict(classification='Proven',analytic_control_contained=True,error_bound=str(error),
        field_error_bound=str(delta),response_interval=regular.interval_text(interval),symbolic_Bernstein_passed=True)


def run():
    mp.mp.dps=80;mp.iv.dps=80
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    control=controls();save('controls.json',control)
    data=dict(np.load(regular.original.OUT/'field-2e-12.npz'));rat=regular.rational
    x=list(map(rat,data['x']));field=list(map(rat,data['field']));derivative=list(map(rat,data['field_derivative']))
    assert len(field)==len(x)-1==len(derivative)
    a0=rat(data['p_coefficients'][-1,0]);b0=rat(data['b_coefficients'][-1,0])
    scale=field[0]/(1+b0*x[1]**2/(6*a0));assert scale>0
    values=[Q(1)]+[v/scale for v in field];slopes=[Q(0)]+[v/scale for v in derivative]
    flux=Q(0);error=Q(0);rows=[]
    for i,(left,right) in enumerate(zip(x,x[1:])):
        h=right-left;change=values[i+1]-values[i]
        p=[values[i],slopes[i],3*change/h**2-(2*slopes[i]+slopes[i+1])/h,
            -2*change/h**3+(slopes[i]+slopes[i+1])/h**2]
        assert at(p,h)==values[i+1] and at([Q(j)*v for j,v in enumerate(p) if j],h)==slopes[i+1]
        a=list(map(rat,data['p_coefficients'][::-1,i]));b=list(map(rat,data['b_coefficients'][::-1,i]))
        flux,err=piece(left,h,a,b,p,flux);error+=err
        rows.append(dict(piece=i,exact_residual_bound=str(err),approximate_residual_bound=float(err)))
    prior=json.loads((regular.OUT/'result.json').read_text());r={k:Q(v) for k,v in prior['exact_rationals'].items()}
    L=-mp.iv.ln(1-2*regular.iv(r['mu']))/(2*regular.iv(r['mu']))
    enclosure,delta,normal=response(values[-1],flux,error,r['kappa'],r['B2'],L,r['R_over_M'])
    old=json.loads((regular.original.OUT/'result.json').read_text())['rows'][-1]['alpha_over_phi_infinity']
    value=regular.iv(rat(old));inside=bool(enclosure.a<=value.a and value.b<=enclosure.b)
    save('piece-residuals.json',dict(classification='Proven',rows=rows))
    save('result.json',dict(classification='Proven',completed=True,pieces=len(rows),
        exact_rationals=dict(interpolant_scale=str(scale),p_surface=str(values[-1]),J_surface=str(flux),
            residual_bound=str(error),field_error_bound=str(delta)),
        approximate_residual_bound=float(error),approximate_field_error_bound=float(delta),
        normalization_interval=regular.interval_text(normal),response_interval=regular.interval_text(enclosure),
        approximate_response_interval=[float(enclosure.a),float(enclosure.b)],
        maximum_piece_residual=max(rows,key=lambda r:r['approximate_residual_bound'])['piece'],
        prior_numerical_coefficient=old,prior_numerical_coefficient_inside=inside,
        analytic_control_passed=True,physical_EOS_certified=False,GR_coefficient_error_certified=False,
        finite_amplitude_backreaction=False,dynamic_orbital_observable=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    regular.verify()
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and r['analytic_control_passed']
    print('PASS exact scalar residual certificate; consult prior coefficient inclusion, physical errors remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
