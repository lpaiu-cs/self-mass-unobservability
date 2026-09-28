"""Exact-coefficient Volterra bound for the frozen linear scalar problem."""
from fractions import Fraction as Q
from math import comb
import json, sys
import mpmath as mp
import numpy as np
import sympy as sp
import gr_scalar_pressure_readout as original

g=original.g;OUT=g.OUT/'gr-scalar-global-regular'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.verify()
    paths=[g.ROOT/'verification/gr_scalar_global_regular.py',original.OUT/'manifest.json',
        original.OUT/'field-2e-12.npz',original.OUT/'result.json']
    save('plan.json',dict(classification='Proven',checkpoint='2c76707',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        model='Treat stored binary64 knots and polynomial coefficients as exact rationals defining bounded piecewise cubics on [0,1]. Interior flux and field are continuous weak solutions; tiny coefficient-rounding jumps of a,b are allowed. Use the same rounded mu=M/R for the exact declared Schwarzschild exterior and exact stored R/M for response normalization.',
        method='Exact rational Bernstein bounds certify a>=amin>0 and b<=0 on every closed piece. Exact polynomial integrals give B1=int x*(-b) and B2=int x^2*(-b). kappa=(B1-B2)/amin bounds the Volterra operator norm. Require kappa<1 and 1-kappa-L*B2>0 with outward interval L=-log(1-2mu)/(2mu).',
        conclusions='Conditional existence and uniqueness of the centre-regular unit-central static solution; global positivity including the exterior; unique normalization to prescribed phi_infinity; no nonzero centre-regular static solution with phi_infinity=0. This closes the node-free condition for the declared linear Riccati readout.',
        limits='No error enclosure for EOS, GR geometry, scalar coefficient generation or floating IVP evaluation. No finite-amplitude nonlinear scalar, physical stability/spectral gap, dynamical drive, orbital response or observational inference claim. This is a proof prerequisite, not static coefficient novelty.'))


def rational(v):return Q(float(v))


def bernstein(coeff,h):
    powers=[v*h**j for j,v in enumerate(coeff)]
    return [sum((powers[j]*Q(comb(k,j),comb(3,j)) for j in range(k+1)),Q(0)) for k in range(4)]


def moment(coeff,left,h,power):
    return sum((v*comb(power,k)*left**(power-k)*h**(j+k+1)/Q(j+k+1)
        for j,v in enumerate(coeff) for k in range(power+1)),Q(0))


def iv(q):return mp.iv.mpf(q.numerator)/mp.iv.mpf(q.denominator)


def interval_text(value):return [str(value.a),str(value.b)]


def run():
    mp.mp.dps=80;mp.iv.dps=80
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    data=dict(np.load(original.OUT/'field-2e-12.npz'));x=[rational(v) for v in data['x']]
    assert x[0]==0 and x[-1]==1 and all(a<b for a,b in zip(x,x[1:]))
    amin=None;bmax=None;B1=Q(0);B2=Q(0);bounds=[]
    for i,(left,right) in enumerate(zip(x,x[1:])):
        h=right-left;ac=[rational(v) for v in data['p_coefficients'][::-1,i]]
        bc=[rational(v) for v in data['b_coefficients'][::-1,i]]
        ab=bernstein(ac,h);bb=bernstein(bc,h)
        alo=min(ab);bhi=max(bb);assert alo>0 and bhi<=0,(i,float(alo),float(bhi))
        amin=alo if amin is None else min(amin,alo);bmax=bhi if bmax is None else max(bmax,bhi)
        B1+=moment([-v for v in bc],left,h,1);B2+=moment([-v for v in bc],left,h,2)
        bounds.append([float(alo),float(bhi)])
    assert 0<B2<B1;kappa=(B1-B2)/amin;assert 0<kappa<1
    row=json.loads((original.OUT/'result.json').read_text())['rows'][-1]
    mu=rational(row['mass_geom_m']/row['radius_m']);scale=rational(row['radius_m'])/rational(row['mass_geom_m'])
    L=-mp.iv.ln(1-2*iv(mu))/(2*iv(mu));normal=1-iv(kappa)-L*iv(B2)
    assert normal.a>0
    lower=-iv(scale*B2)/normal;upper=-iv(scale*(1-kappa)*B2)
    coefficient=iv(rational(row['alpha_over_phi_infinity']))
    assert lower.b<coefficient.a and coefficient.b<upper.a
    # Independent polynomial controls: constant a=1, b=-q^2 has phi=sin(q*x)/(q*x).
    z,q=sp.symbols('z q',positive=True);phi=sp.sin(q*z)/(q*z)
    flux=z*z*sp.diff(phi,z)
    assert sp.simplify(sp.diff(flux,z)+q*q*z*z*phi)==0
    assert sp.simplify((phi+flux).subs(z,1)-sp.cos(q))==0
    assert sp.cos(sp.pi/2)==0
    controls=[]
    for wave,expected in [(Q(3,10),True),(Q(2),False)]:
        potential=[-wave*wave,Q(0),Q(0),Q(0)]
        u=moment([-v for v in potential],Q(0),Q(1),1)
        v=moment([-v for v in potential],Q(0),Q(1),2)
        assert u==wave*wave/2 and v==wave*wave/3
        bound=1-(u-v)-v;accepted=u-v<1 and bound>0
        assert accepted==expected
        if accepted:assert iv(bound).b<mp.iv.cos(iv(wave)).a
        controls.append(dict(q=str(wave),normalization_lower_bound=str(bound),certificate_accepted=accepted,
            interpretation='The negative control refuses a positive exterior normalization; it does not assert nonuniqueness.'))
    # Manufactured cubic verifies both the exact moments and Bernstein conversion.
    polynomial=[Q(2),Q(-3),Q(5),Q(7)];left=Q(1,7);h=Q(2,9)
    for power in [1,2]:
        expr=sum(sp.Rational(v.numerator,v.denominator)*(z-sp.Rational(left.numerator,left.denominator))**j for j,v in enumerate(polynomial))*z**power
        exact=sp.integrate(expr,(z,sp.Rational(left.numerator,left.denominator),sp.Rational((left+h).numerator,(left+h).denominator)))
        assert moment(polynomial,left,h,power)==Q(int(exact.p),int(exact.q))
    coeff=bernstein(polynomial,h)
    expr=sum(sp.Rational(v.numerator,v.denominator)*comb(3,j)*z**j*(1-z)**(3-j) for j,v in enumerate(coeff))
    target=sum(sp.Rational(v.numerator,v.denominator)*(sp.Rational(h.numerator,h.denominator)*z)**j for j,v in enumerate(polynomial))
    assert sp.expand(expr-target)==0
    np.savez_compressed(OUT/'piece-controls.npz',bounds=np.array(bounds))
    save('result.json',dict(classification='Proven',completed=True,pieces=len(bounds),
        exact_rationals={k:str(v) for k,v in dict(a_min=amin,b_max=bmax,B1=B1,B2=B2,kappa=kappa,mu=mu,R_over_M=scale).items()},
        approximate_summary={k:float(v) for k,v in dict(a_min=amin,b_max=bmax,B1=B1,B2=B2,kappa=kappa).items()},
        exterior_L=interval_text(L),normalization_lower_bound=interval_text(normal),
        response_lower_bound=interval_text(lower),response_upper_bound=interval_text(upper),
        saved_numerical_coefficient_inside=True,controls=controls,symbolic_controls_passed=True,
        centre_regular_solution_unique=True,global_node_free=True,static_zero_boundary_kernel_trivial=True,
        physical_EOS_certified=False,finite_amplitude_backreaction=False,dynamic_orbital_observable=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['completed'] and r['symbolic_controls_passed'] and r['global_node_free']
    assert r['controls'][0]['certificate_accepted'] and not r['controls'][1]['certificate_accepted']
    print('PASS exact-coefficient global scalar regularity and manufactured controls; conditional frozen model only',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
