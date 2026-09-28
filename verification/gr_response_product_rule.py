"""Reusable Chebyshev interpolation and analytic logarithmic product moments."""
from fractions import Fraction as F
import json,math,sys
import numpy as np
import sympy as sp
import mpmath as mp
from mpmath import iv
import gr_self_regular_domain as domain

highq=domain.highq;cusp=domain.cusp;ROOT=domain.ROOT;OUT=highq.OUT.parent/'gr-response-product-rule';I=domain.I;low=domain.low;high=domain.high;sha=domain.sha

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    domain.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_response_product_rule.py',domain.OUT/'manifest.json',domain.OUT/'result.json',
        highq.moments.OUT/'manifest.json',highq.moments.OUT/'inputs.npz',cusp.THERMO/'states.npz',highq.sector.OUT/'manifest.json']
    save('plan.json',dict(classification='Proven',checkpoint='0442d55',bits=256,interpolation_order=32,analytic_radius_ratio=4,
        far_log_series_terms=96,control_degrees=[0,1,2,8,16,32],control_z=['0','1','-1','7/4','2','-2','10','-10','1000000'],
        control_digits=320,control_allowance='1e-280',thermal_tail_beta_ratio=16,tail_response_budget='1e-15',
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        target='Factor the same response through H=f/gamma+2beta*log(1+exp(eta-t)). Certify a reusable degree-31 interpolation rule and the real logarithmic moments needed to integrate each polynomial at arbitrary Q. Also bound the entire omitted H tail, including its cusp, uniformly in Q.',
        scope='Rule, moment evaluation and real-state tail bounds only. No actual stellar H interpolation table, complete response product integral or outer self integral is claimed.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def logarithmic_moment(n,z,terms=96):
    z=F(z);x=I(z)
    if abs(z)==1:return 2*iv.log(2)-2 if n==0 else I(F(-2*(-1 if z<0 and n%2 else 1),n*(n+1)))
    if abs(z)<2:
        if n==0:return (x+1)*iv.log(abs(x+1))-(x-1)*iv.log(abs(x-1))-2
        q0=iv.log(abs((x+1)/(x-1)))/2;P=[I(1),x];W=[I(0),I(1)]
        for k in range(1,n+1):
            P.append(((2*k+1)*x*P[k]-k*P[k-1])/(k+1));W.append(((2*k+1)*x*W[k]-k*W[k-1])/(k+1))
        return 2*((P[n+1]-P[n-1])*q0-W[n+1]+W[n-1])/(2*n+1)
    k=2 if n==0 else n;moment=F(2**(n+1)*math.factorial(k)*math.factorial((k+n)//2),math.factorial((k-n)//2)*math.factorial(k+n+1))
    term=I(moment)/k/x**k;total=I(0)
    for _ in range(terms):
        total+=term;term*=I(k*(k+1))/((k-n+2)*(k+n+3)*x*x);k+=2
    remainder=2/(I(k*(k+1))*abs(x)**k*(1-1/(x*x)))
    return (2*iv.log(abs(x)) if n==0 else I(0))-total+cusp.symmetric(high(remainder))


def real_response_upper(b):
    # The existing positive Q^2 response bound evaluated at the actual enlarged beta.
    e=iv.exp(-I(1));half=[I(2)/3+e*iv.sqrt(2),I(2)/5+e*iv.sqrt(10),I(2)/7+e*iv.sqrt(80)];integer=[I(1)/2+2*e,I(1)/3+5*e,I(1)/4+16*e]
    J2=b*iv.sqrt(2*b)*(half[0]+2*b*half[1]+b*b*half[2])+b*b*(integer[0]+2*b*integer[1]+b*b*integer[2]);power=lambda s:I(s)**I(s)
    P3=b*iv.sqrt(2*b)*(2*power(F(3,2))+3*b*power(F(5,2))+b*b*power(F(7,2)))+b*b*(2*power(2)+3*b*power(3)+b*b*power(4))
    return I(max(high(highq.sector.g.real_majorant(b/2,I(0))),high(I(32)/3*J2+48*P3)))


def prove():
    p=bindings();iv.prec=p['bits'];N=p['interpolation_order'];x=sp.symbols('x');polys=[sp.Poly(sp.chebyshevt(n,x),x) for n in range(2*N-1)]
    coefficients=polys[N].monic().all_coeffs();power_sums=[sp.Integer(N)]
    for k in range(1,2*N-1):
        if k<=N:power_sums.append(-sum(coefficients[j]*power_sums[k-j] for j in range(1,k))-k*coefficients[k])
        else:power_sums.append(-sum(coefficients[j]*power_sums[k-j] for j in range(1,N+1)))
    traces=[sum(poly.nth(k)*power_sums[k] for k in range(poly.degree()+1)) for poly in polys];assert traces==[N]+[0]*(2*N-2)
    for n in range(N):
        for m in range(N):assert (traces[n+m]+traces[abs(n-m)])/2==(N if n==m==0 else N//2 if n==m else 0)
    nodes=[]
    for j in range(N):
        assert sp.cos(N*sp.Rational(2*j+1,2*N)*sp.pi)==0
        node=iv.cos(I(F(2*j+1,2*N))*iv.pi);nodes.append(cusp.interval_text(node))
    connections=[]
    for k in range(N):
        remaining=polys[k];row=[sp.Integer(0)]*N
        for n in range(k,-1,-1):
            legendre=sp.Poly(sp.legendre(n,x),x);row[n]=remaining.nth(n)/legendre.LC();remaining-=row[n]*legendre
        assert remaining.is_zero;connections.append(list(map(str,row)))
    h,Q,z=sp.symbols('h Q z',positive=True);assert sp.cancel(2*h/(4*h-h)-sp.Rational(2,3))==0
    for n in range(1,8):
        for k in range(n,n+10,2):
            want=sp.Rational(2**(n+1)*math.factorial(k)*math.factorial((k+n)//2),math.factorial((k-n)//2)*math.factorial(k+n+1))
            assert sp.integrate(x**k*sp.legendre(n,x),(x,-1,1))==want
    primitive=x/2+(x*x-Q*Q)/(4*Q)*sp.log((x+Q)/(x-Q))
    assert sp.simplify(sp.diff(primitive,x)-x/(2*Q)*sp.log((x+Q)/(x-Q)))==0
    save('rule.json',dict(classification='Proven',passed=True,order=N,nodes=nodes,trace_chebyshev=list(map(str,traces)),chebyshev_to_legendre=connections,
        interpolation_error_factor_exact=str(F(2,6**N)),
        identity='From W_prime=-pH and the regular response identity, S(Q)=integral_0^infinity p*H(p)/(2Q)*log|(p+Q)/(p-Q)| dp, with its continuous Q=0 limit. This positive logarithmic kernel permits one H table to be reused across Q.',
        discrete_rule='The roots cos((2j+1)pi/64) are distinct roots of T32. Exact Newton sums give trace(Tk)=0 for k=1..62. The Chebyshev product identity proves the discrete orthogonality matrix. Thus a0=sum H_j/32 and a_n=2sum H_j*Tn(x_j)/32 recover the exact interpolating polynomial. The saved rational connection matrix maps it to Legendre coefficients.',
        pointwise_error='On a real panel with center c and half width h, require analyticity and modulus M on |p-c|<=R=4h. Cauchy bounds the 32nd derivative by 32!M/(R-h)^32. The Chebyshev node polynomial is h^32*T32((p-c)/h)/2^31, hence |H-H_interp|<=2M/6^32. Include all directed node, value and transform rounding separately.',
        global_error='The positive kernel has A_Q(P)=P/2+(P^2-Q^2)/(4Q)*log|(P+Q)/(P-Q)|, with A_0(P)=P and A_Q(Q)=Q/2. For Q<P use log((1+x)/(1-x))<=2x/(1-x^2), x=Q/P, to get A_Q(P)<=P; for Q>=P the second term is nonpositive. Thus any uniform normalized H error epsilon on [0,P] changes each normalized response by at most P*epsilon, even when Q lies inside a panel.',
        product_integration='Write p*H_interp=sum b_n Pn((p-c)/h). Then its panel contribution is h/(2Q)*sum b_n[Jn(-(c+Q)/h)-Jn((Q-c)/h)], Jn(z)=integral_-1^1 Pn(t)log|z-t|dt. Multiplication by p uses the exact three-term Legendre recurrence. The two log(h) constants cancel.',
        logarithmic_moments='For n>=1, Jn=2(Q_(n+1)-Q_(n-1))/(2n+1), with the real Legendre second-kind boundary value. Near the cut use Qn=Pn*Q0-Wn and polynomial recurrences; at z=+/-1 take Jn=-2*(sign z)^n/[n(n+1)]. J0=(z+1)log|z+1|-(z-1)log|z-1|-2. For |z|>=2 expand log|z-t| in powers of t/z, using the exact Legendre moments. After the last retained parity power, the tail is <=2/[K(K+1)|z|^K(1-z^-2)].'))
    data=dict(np.load(highq.moments.OUT/'inputs.npz'));state=dict(np.load(cusp.THERMO/'states.npz'));rows=[];maximum=F(0);ratio=p['thermal_tail_beta_ratio'];a=I(1)-I(1)/ratio
    for i,cell in enumerate(data['cells']):
        beta=I(data['beta'][i]);P=I(data['pcut'][i]);T=I(low(P*P/(beta*(iv.sqrt(1+P*P)+1))));assert low(a*T)>2
        eta=max(cusp.endpoints(str(data['root_intervals'][i]))[1],F.from_float(float(state['eta'][i])));C=real_response_upper(ratio*beta)
        base=2*iv.exp(I(eta)-a*T)*C/I(data['Sref'][i]);bounds=[base,base,base,base*(1+T),base*(1+T),base*(1+T+T*T)]
        upper=[high(v) for v in bounds];maximum=max(maximum,max(upper));rows.append(dict(position=i,cell=int(cell),tc_lower=str(low(T)),upper_coefficients=list(map(str,upper))))
    save('tail.json',dict(classification='Proven',mathematical_bounds_valid=True,cells=len(rows),rows=rows,budget_passed=maximum<F(p['tail_response_budget']),maximum_upper_exact=str(maximum),display_only=float(maximum),
        proof='For p>=P, t>=T, both f and log(1+exp(eta-t)) are <=exp(eta-t). The six H jets are bounded by exp(eta-t)*(1/gamma+2beta)*[1,1,1,1+t,1+t,1+t+t^2]. Compare with H at eta=0 and beta enlarged by lambda=16, which is >=exp(-t/lambda)*(1/gamma+2beta)/2. Each polynomial times exp(-(1-1/lambda)t) decreases for t>=T. Integrate the positive logarithmic kernel over the enlarged full momentum axis to obtain the saved coefficient divided by max(1,Q^2). It includes the entire omitted tail cusp. The real-response coefficient uses the existing positive polynomial bound evaluated at the actual enlarged beta.',
        full_response_product_integral_certified=False))


def controls():
    p=bindings();iv.prec=p['bits'];mp.mp.dps=p['control_digits'];checks=[]
    def M(x):
        q=F(x);return mp.mpf(q.numerator)/q.denominator
    for ztext in p['control_z']:
        z=F(ztext);mz=M(z)
        for n in p['control_degrees']:
            value=logarithmic_moment(n,z,p['far_log_series_terms']);lo,hi=cusp.endpoints(cusp.interval_text(value));cuts=[mp.mpf(-1)]
            if -1<z<1:cuts.append(mz)
            cuts.append(mp.mpf(1));reference=mp.quad(lambda t:mp.legendre(n,t)*mp.log(abs(mz-t)),cuts)
            passed=M(lo)-M(p['control_allowance'])<=reference<=M(hi)+M(p['control_allowance']);assert passed,(n,z,reference,lo,hi)
            checks.append(dict(degree=n,z=ztext,enclosure=cusp.interval_text(value),reference=str(reference),passed=bool(passed)))
        print('PRODUCT LOG MOMENTS',ztext,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=True,checks=checks,
        scope='Independent 320-digit direct polynomial/logarithm integration with declared 1e-280 numerical comparison allowance. The analytic moment enclosures do not depend on these numerical controls.'))


def finalize():
    bindings();assert json.loads((OUT/'rule.json').read_text())['passed'] and json.loads((OUT/'controls.json').read_text())['passed']
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'rule.json').read_text());assert r['passed'] and r['order']==32 and len(r['nodes'])==32
    tail=json.loads((OUT/'tail.json').read_text());assert tail['mathematical_bounds_valid'] and tail['cells']==3206
    maximum=max(F(v) for row in tail['rows'] for v in row['upper_coefficients']);assert maximum==F(tail['maximum_upper_exact']) and tail['budget_passed']==(maximum<F(p['tail_response_budget']))
    controls=json.loads((OUT/'controls.json').read_text());assert controls['passed'] and len(controls['checks'])==54
    print('PASS reusable interpolation/product rule and actual tail bounds; tail budget passed:',tail['budget_passed'],'actual H coefficients/product response remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
