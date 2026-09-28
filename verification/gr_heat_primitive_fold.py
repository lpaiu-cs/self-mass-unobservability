"""Two certified point primitives with identical conserved data and proper heat."""
from fractions import Fraction as F
import json,sys
import sympy as sp
from mpmath import iv
import gr_heat_primitive_global as prior
import gr_outer_product_pilot as outer

g=prior.g;OUT=g.OUT/'gr-heat-primitive-fold';sha=prior.sha;cusp=outer.cusp;I=cusp.I;low=cusp.low;high=cusp.high


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    prior.verify();assert not OUT.exists();OUT.mkdir()
    files=[g.ROOT/'verification/gr_heat_primitive_fold.py',prior.OUT/'manifest.json',prior.OUT/'symbolic.json']
    save('plan.json',dict(classification='Proven',checkpoint='2f8cbd15',bits=256,
        initial_offset_halfwidth='1e-4',momentum_shift='-1e-12',root_offset_width_gate='1e-25',maximum_bisections=96,
        bindings={p.relative_to(g.ROOT).as_posix():sha(p) for p in files},
        target='Use the exact polynomial Helmholtz EOS from the previous local-jet control. Hold D,E,Q fixed at its singular state and perturb only J by the declared negative rational amount. Prove two distinct admissible roots by outward interval signs, continuity and bisection; do not infer nonuniqueness from a singular Jacobian alone.',
        source_conventions='Dimensionless natural units. The EOS is a constructed local thermodynamic counterexample, not the stellar EOS. The interval calculations use mpmath.iv at256 bits, not the native CAPD/MPFR runtime.',
        boundary='Equilibrium EOS causality and algebraic dominant energy only. No causal nonequilibrium transport model, microscopic realization, actual stellar ambiguity or GR trajectory is asserted.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(g.ROOT/rel)==digest,rel
    return p


def coefficients():
    rho,T,v=sp.symbols('rho T v',positive=True)
    free=sp.sympify(json.loads((prior.OUT/'symbolic.json').read_text())['negative_control_specific_free_energy'],locals={'rho':rho,'T':T})
    eps=sp.expand(rho*(free-T*sp.diff(free,T)));P=sp.expand(rho*rho*sp.diff(free,rho))
    A=sp.expand((eps/rho).subs(T,0));B=sp.simplify(sp.diff(eps,T,2)/(2*rho))
    assert B==sp.Rational(11,2000) and sp.expand(eps-rho*(A+B*T*T))==0
    assert sp.expand(P-rho*rho*(sp.Rational(32,125)-sp.Rational(189,1000)*rho+sp.Rational(33,1000)*T))==0
    c=sp.sqrt(46)-6;W2=1/(1-c*c);Q=sp.Rational(33,100);E=(1+c*c/10+2*c*Q)*W2;J=(sp.Rational(11,10)*c+Q*(1+c*c))*W2
    point={rho:1,T:1};er,et,err,ert,ett=[sp.diff(eps,*args).subs(point) for args in [(rho,),(T,),(rho,rho),(rho,T),(T,T)]]
    pr,pt,prr,prt,ptt=[sp.diff(P,*args).subs(point) for args in [(rho,),(T,),(rho,rho),(rho,T),(T,T)]]
    rp=-c*W2;rpp=-W2*W2;tp=sp.simplify((-(J+Q)-er*rp)/et)
    tpp=sp.simplify(-(err*rp*rp+2*ert*rp*tp+ett*tp*tp+er*rpp)/et)
    pp=pr*rp+pt*tp;ppp=prr*rp*rp+2*prt*rp*tp+ptt*tp*tp+pr*rpp+pt*tpp
    second=sp.simplify(2*pp+c*ppp);assert second<0
    parameter=sp.simplify(-1-c*c*pt/et);assert parameter<0
    save('symbolic.json',dict(classification='Proven',passed=True,free_energy=str(free),singular_velocity=str(c),
        residual_second_velocity_derivative=str(second),residual_momentum_derivative=str(parameter),
        reduction='epsilon/rho=A(rho)+(11/2000)*T^2. The positive thermal solution is the explicit positive square root. Both the second velocity derivative and the fixed-(D,E,Q) momentum derivative of F are strictly negative at the original singular state. The actual two roots are established separately by interval signs, not by this local derivative test alone.'))
    return [F(str(x)) for x in sp.Poly(A,rho).all_coeffs()]


def polynomial(coeff,x):
    value=I(0)
    for c in coeff:value=value*x+I(c)
    return value


def run():
    p=bindings();iv.prec=p['bits'];coef=coefficients();c=iv.sqrt(I(46))-6;Q=I(F(33,100));W2=1/(1-c*c)
    D=iv.sqrt(W2);E=(1+c*c/10+2*c*Q)*W2;J0=(I(F(11,10))*c+Q*(1+c*c))*W2;J=J0+I(F(p['momentum_shift']))
    def evaluate(a,b=None):
        offset=I(a) if b is None else iv.mpf([I(a).a,I(b).b]);v=c+offset;root=iv.sqrt(1-v*v);rho=D*root
        target=E-v*(J+Q);T2=(target/rho-polynomial(coef,rho))/I(F(11,2000));assert low(T2)>0
        T=iv.sqrt(T2);P=rho*rho*(I(F(32,125))-I(F(189,1000))*rho+I(F(33,1000))*T)
        epsilon=rho*(polynomial(coef,rho)+I(F(11,2000))*T*T);w=epsilon+P
        residual=v*(E+P)-(J-Q)
        ar=2*I(coef[0])*rho+I(coef[1]);er=polynomial(coef,rho)+rho*ar+I(F(11,2000))*T*T;et=I(F(11,1000))*rho*T
        pr=2*P/rho-I(F(189,1000))*rho*rho;pt=I(F(33,1000))*rho*rho
        BB=T*et/w;RR=rho*pr/w;DD=T*pt/w;EE=1-rho*er/w;sound=RR+DD*EE/BB
        Wsq=1/(1-v*v);rp=-rho*v*Wsq;tp=(-(J+Q)-er*rp)/et;Pprime=pr*rp+pt*tp;derivative=E+P+v*Pprime
        discriminant=w*w-4*Q*Q;assert low(discriminant)>0
        radical=iv.sqrt(discriminant);landau=(epsilon-P+radical)/2;radial=(-epsilon+P+radical)/2
        margins=[BB,RR,sound,1-sound,landau-abs(radial),landau-abs(P)]
        assert min(map(low,margins))>0
        return dict(residual=residual,derivative=derivative,velocity=v,rho=rho,T=T,epsilon=epsilon,P=P,thermal_square=T2,
            physical_margins=margins,forward=[rho*iv.sqrt(Wsq),(epsilon+P*v*v+2*v*Q)*Wsq,(w*v+Q*(1+v*v))*Wsq])
    def sign(value):
        if low(value)>0:return 1
        if high(value)<0:return -1
        raise AssertionError('unresolved interval sign')
    h=F(p['initial_offset_halfwidth']);whole=evaluate(-h,h);left=evaluate(-h);center=evaluate(F(0));right=evaluate(h)
    assert [sign(x['residual']) for x in [left,center,right]]==[-1,1,-1]
    rows=[]
    for index,(a,b) in enumerate([(-h,F(0)),(F(0),h)]):
        sa=sign(evaluate(a)['residual']);sb=sign(evaluate(b)['residual']);assert sa*sb==-1
        steps=0
        while b-a>F(p['root_offset_width_gate']):
            assert steps<p['maximum_bisections'];mid=(a+b)/2;sm=sign(evaluate(mid)['residual'])
            if sm==sa:a=mid;sa=sm
            else:b=mid;sb=sm
            steps+=1
        value=evaluate(a,b);slope=sign(value['derivative']);assert slope==(1 if index==0 else -1)
        for actual,target in zip(value['forward'],[D,E,J],strict=True):assert low(actual)<=low(target)<=high(target)<=high(actual)
        rows.append(dict(index=index,offset_bracket=list(map(str,[a,b])),endpoint_signs=[sa,sb],bisections=steps,
            primitive_enclosures={key:cusp.interval_text(value[key]) for key in ['velocity','rho','T','epsilon','P']},
            residual_derivative=cusp.interval_text(value['derivative']),derivative_sign=slope,
            minimum_admissibility_margin=str(min(map(low,value['physical_margins']))),conserved_forward_enclosures=list(map(cusp.interval_text,value['forward']))))
    assert F(rows[0]['offset_bracket'][1])<F(rows[1]['offset_bracket'][0])
    save('result.json',dict(classification='Proven',passed=True,distinct_primitives=2,
        identical_conserved_targets={name:cusp.interval_text(value) for name,value in zip(['D','E','J','Q'],[D,E,J,Q],strict=True)},
        initial_residual_enclosures=list(map(cusp.interval_text,[left['residual'],center['residual'],right['residual']])),
        full_initial_interval_thermal_square=cusp.interval_text(whole['thermal_square']),
        full_initial_interval_minimum_admissibility_margin=str(min(map(low,whole['physical_margins']))),rows=rows,
        proof='On the full initial velocity interval T(v)>0 is continuous by the positive square-root formula. The certified signs[-,+,-] prove at least one root in each disjoint side bracket by IVT. Bisection preserves those signs; the final boxes are disjoint and have strictly opposite derivative signs. Every root reconstructs exactly the same D,E,J,Q by the conserved identities. The explicit forward enclosures are an additional consistency check. Equilibrium heat-capacity/stiffness/sound and DEC margins are positive throughout the initial bracket.',
        clarification='This strengthens the preceding singular-Jacobian control to genuine nonuniqueness. A singular Jacobian alone would only exclude a nonsingular C1 inverse, not prove two primitives.',
        physical_stellar_EOS=False,nonequilibrium_transport_causality_certified=False,actual_stellar_nonuniqueness=False))
    save('manifest.json',dict(sha256={f.relative_to(g.ROOT).as_posix():sha(f) for f in OUT.iterdir() if f.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['distinct_primitives']==2
    assert F(r['full_initial_interval_minimum_admissibility_margin'])>0
    for row in r['rows']:
        a,b=map(F,row['offset_bracket']);assert 0<b-a<=F(p['root_offset_width_gate']) and row['endpoint_signs'][0]*row['endpoint_signs'][1]==-1
        assert F(row['minimum_admissibility_margin'])>0
    assert F(r['rows'][0]['offset_bracket'][1])<F(r['rows'][1]['offset_bracket'][0])
    assert [x['derivative_sign'] for x in r['rows']]==[1,-1]
    print('PASS two disjoint admissible primitive boxes with identical D,E,J,Q for an explicit thermodynamic counterexample',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
