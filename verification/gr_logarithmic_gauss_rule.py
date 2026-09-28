"""Exact four-point Gauss rules for the RPA kernel's logarithmic singularity."""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,math,sys
import sympy as sp
from mpmath import iv
from interval_records import exact_endpoint,interval_text

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/direct-eos-gr33/gr-logarithmic-gauss-rule'
N=4


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def I(value):
    q=F(str(value));return iv.mpf(q.numerator)/q.denominator
def low(value):return exact_endpoint(value._mpi_[0])
def high(value):return exact_endpoint(value._mpi_[1])
def evaluate(poly,value):
    result=iv.mpf(0)
    for coefficient in poly.all_coeffs():result=result*value+I(coefficient)
    return result


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_logarithmic_gauss_rule.py',ROOT/'verification/interval_records.py']
    save('plan.json',dict(classification='Proven',checkpoint='e43523b',order=N,bits=128,root_isolation_bits=160,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        weights=['1','-log(x)'],domain=[0,1],
        objective='Certify the ordinary and logarithmic Gauss rules needed to integrate the explicit log|1-u| term of the positive finite-wavenumber response kernel.',
        method='Exact rational moments, monic orthogonal polynomials, rational root isolation, positive Christoffel weights, and exact companion-matrix traces prove degree-seven exactness. Hermite interpolation and Rolle give the degree-eight derivative remainder.',
        physical_application='On 1-h<=u<=1+h, -B(u)*log|1-u| becomes h*[(-log(v))-log(h)]*B(1 +/- h*v), 0<=v<=1. The ordinary and logarithmic rules integrate the two terms separately.',
        actual_integrand_derivative_bounds_certified=False,whole_EOS_certified=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def rule(power):
    x=sp.symbols('x');mu=lambda k:sp.Rational(1,(k+1)**power)
    polys=[];norms=[]
    for degree in range(N+1):
        if degree:
            hankel=sp.Matrix(degree,degree,lambda i,j:mu(i+j))
            coeff=hankel.inv()*sp.Matrix([-mu(degree+i) for i in range(degree)])
            poly=sp.Poly(x**degree+sum(coeff[j]*x**j for j in range(degree)),x)
        else:poly=sp.Poly(1,x)
        for k in range(degree):assert sum(poly.nth(j)*mu(j+k) for j in range(degree+1))==0
        square=poly*poly;norm=sum(square.nth(j)*mu(j) for j in range(2*degree+1));assert norm>0
        polys.append(poly);norms.append(norm)
    final=polys[N];kernel=sum(polys[k].as_expr()**2/norms[k] for k in range(N))
    weight_poly=sp.Poly(sp.invert(kernel,final.as_expr(),x),x)
    companion=sp.zeros(N)
    for i in range(1,N):companion[i,i-1]=1
    for i in range(N):companion[i,N-1]=-final.nth(i)
    W=sum((weight_poly.nth(k)*companion**k for k in range(N)),sp.zeros(N))
    errors=[sp.factor(mu(k)-sp.trace(W*companion**k)) for k in range(2*N+1)]
    assert all(value==0 for value in errors[:2*N]) and errors[2*N]==norms[N]
    isolation=sp.intervals(final,eps=sp.Rational(1,2**160));assert len(isolation)==N
    nodes=[];weights=[];encoded=[]
    for (a,b),multiplicity in isolation:
        assert multiplicity==1 and 0<a<=b<1
        node=iv.mpf([I(a).a,I(b).b]);weight=1/sum(evaluate(polys[k],node)**2/I(norms[k]) for k in range(N))
        assert low(weight)>0;nodes.append(node);weights.append(weight)
        encoded.append(dict(root_isolation=[str(a),str(b)],node=interval_text(node),weight=interval_text(weight)))
    for k in range(2*N):
        observed=sum(w*v**k for v,w in zip(nodes,weights));exact=F(str(mu(k)))
        assert low(observed)<=exact<=high(observed)
    # Independent non-polynomial reference: the positive exponential moment series.
    partial=sum(F(1,math.factorial(k)*(k+1)**power) for k in range(61))
    first=F(1,math.factorial(61)*62**power);tail=first/(1-F(1,62))
    reference=iv.mpf([I(partial).a,I(partial+tail).b])
    quadrature=sum(w*iv.exp(v) for v,w in zip(nodes,weights));error=reference-quadrature
    constant=norms[N]/math.factorial(2*N)
    assert low(error)>0 and low(error)>=low(I(constant)) and high(error)<=high(iv.exp(1)*I(constant))
    return dict(weight='1' if power==1 else '-log(x)',moment_formula=f'1/(k+1)^{power}',
        monic_polynomial=str(final.as_expr()),norm_exact=str(norms[N]),error_coefficient_exact=str(constant),
        weight_polynomial_mod_orthogonal=str(weight_poly.as_expr()),
        degree_0_to_8_errors_exact=list(map(str,errors)),nodes=encoded,
        exponential_control=dict(passed=True,reference=interval_text(reference),quadrature=interval_text(quadrature),error=interval_text(error)))


def run():
    p=bindings();iv.prec=p['bits'];rows=[rule(1),rule(2)]
    assert F(rows[0]['norm_exact'])==F(1,9*70**2)
    save('result.json',dict(classification='Proven',passed=True,rules=rows,
        proof='The orthogonal polynomial has exactly four isolated simple roots in (0,1). Its Christoffel weights are positive. Reducing their rational expression modulo the polynomial and taking companion-matrix traces proves all eight exact moments over Q, not just interval containment. A degree-seven Hermite interpolant at the four roots has pointwise remainder f^(8)(xi)*pi4(x)^2/8! by repeated Rolle. Integrating against the positive weight gives |error|<=norm(pi4)^2*sup|f^(8)|/8!, where norm(pi4)^2 is the stored squared norm.',
        logarithmic_window='For 0<h<=1, the two-sided log contribution is h*sum_{sign +/-}[int_0^1 (-log(v))*B(1+sign*h*v)dv-log(h)*int_0^1 B(1+sign*h*v)dv]. Its absolute quadrature error is <=h^9*(C_log+abs(log(h))*C_unit)*(M8_minus+M8_plus), with M8 bounding the eighth u derivative of B on each side.',
        constant_control='For B=1, both rules integrate exactly and the window gives 2*h*(1-log(h)), the direct primitive result.',
        scope='Exact weighted quadrature and its derivative-based remainder. Actual B derivative bounds, the other momentum ranges, outer wavenumber integration and whole EOS certification remain open.'))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and len(r['rules'])==2
    for rule_data in r['rules']:
        assert all(F(x)==0 for x in rule_data['degree_0_to_8_errors_exact'][:8])
        assert F(rule_data['degree_0_to_8_errors_exact'][8])==F(rule_data['norm_exact'])
        assert rule_data['exponential_control']['passed']
    print('PASS exact ordinary/logarithmic Gauss rules, positive weights, remainder and exponential controls',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
