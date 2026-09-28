"""Certify Q >= 4P self-integrals using positive thermal moments."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
import gzip,json,math,sys
import numpy as np
import sympy as sp
from mpmath import iv
import gr_highq_moments as moments
import gr_response_positive_sector as sector
import gr_polarization_neutral_error as neutral

cusp=moments.cusp;ROOT=cusp.ROOT;OUT=moments.OUT.parent/'gr-highq-self'
I=cusp.I;low=sector.g.low;high=sector.g.high;sha=moments.sha;FIELDS=neutral.thermo.FIELDS
PAIRS=[[(0,0)],[(0,1),(1,0)],[(0,2),(1,1),(2,0)],[(0,3),(3,0)],[(0,4),(1,3),(3,1),(4,0)],[(0,5),(3,3),(5,0)]]

def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    moments.verify();sector.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_highq_self_certificate.py',moments.OUT/'manifest.json',moments.OUT/'result.json',moments.OUT/'inputs.npz',
        sector.OUT/'manifest.json',sector.OUT/'result.json',cusp.ROOTS/'manifest.json',cusp.ROOTS/'result.json',cusp.THERMO/'states.npz',
        neutral.density.ionic.OUT/'inputs.npz',neutral.density.ionic.OUT/'constants.json',ROOT/'verification/gr_polarization_neutral_error.py',
        ROOT/'verification/gr_polarization_thermodynamics.py']
    files+=sorted(moments.OUT.glob('block-*.jsonl.gz'))
    save('plan.json',dict(classification='Proven',checkpoint='baf8a9d',bits=128,cells=3206,series_order=16,Q_over_P=4,neumann_order=4,
        field_budget='2e-7',processes=3,block_size=32,control_positions=[0,1603,3205],control_relative_vector_gate='2e-7',
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        target='Integrate G=S/(Q^2/B+S) and its five eta/tau partials over Q in [4P,infinity), then divide by fixed scale and apply the certified neutral Hessian and exact ion inventory.',
        scope='All cutoffs P,QH=4P, scale and Sref are held fixed under differentiation at each reference state. This high-Q contribution includes the entire infinite tail and replaces, rather than adds to, the old far-Q tail when assembling this partition. Intermediate Q remains separate. Same declared fully ionized ideal-electron RPA term only; no complete physical EOS/GR/observation claim.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def zeros(n):
    a=np.empty(n,dtype=object)
    for i in range(n):a[i]=I(0)
    return a


def product(a,b):
    return [sum((np.convolve(a[i],b[j]) for i,j in pairs),zeros(len(a[0])+len(b[0])-1)) for pairs in PAIRS]


def symbolic():
    y=sp.symbols('y',positive=True);L=4;r=y**(L+1)/(1+y)
    assert sp.cancel(y/(1+y)-sum((-1)**(k+1)*y**k for k in range(1,L+1))-r)==0
    assert sp.cancel(sp.diff(r,y,2)-y**(L-1)*(L*(L+1)+2*(L*L-1)*y+L*(L-1)*y*y)/(1+y)**3)==0
    assert all(c>=0 for c in sp.Poly(sp.cancel(((L+1)*y**L-sp.diff(r,y))*(1+y)**2/y**L),y).all_coeffs())
    assert all(c>=0 for c in sp.Poly(sp.cancel((L*(L+1)*y**(L-1)-sp.diff(r,y,2))*(1+y)**3/y**(L-1)),y).all_coeffs())
    e,t,v=sp.symbols('e t v');f=(1+2*e+3*t+e*t+e*e+t*t)*(v+v*v);orders=[(0,0),(1,0),(2,0),(0,1),(1,1),(0,2)]
    jet=[]
    for n,m in orders:
        poly=sp.Poly(sp.diff(f,e,n,t,m).subs({e:0,t:0})/sp.factorial(n)/sp.factorial(m),v)
        a=zeros(3)
        for k in range(3):a[k]=I(int(poly.nth(k)))
        jet.append(a)
    out=product(jet,jet)
    for a,(n,m) in zip(out,orders):
        exact=sp.Poly(sp.diff(f*f,e,n,t,m).subs({e:0,t:0})/sp.factorial(n)/sp.factorial(m),v)
        for k,x in enumerate(a):assert low(x)==high(x)==F(int(exact.nth(k)))
    save('symbolic.json',dict(classification='Proven',passed=True,
        series='Set v=(QH/Q)^2, QH=4P. R_N=sum_{m=1}^16 s_m v^m, s_m=16^-m[mu_m/(2m-1)+2P^2 mu_(m+1)/((2m-1)(2m+1))]. All true value coefficients are positive. mu_18<=mu_17 implies value remainder <= mu_17*[1/33+2P^2/(33*35)]*v^17/(15*16^16). Its six jets are bounded by multiplying [1,1,1,Tcap,Tcap,Tcap^2+Tcap], where Tcap is the upper t(P).',
        tail_inner='For P<=p<=Q/2, the positive kernel series gives Kp<=4p^2 gamma/(3Q^2). At t>=tc=64+max(eta_upper,0), p<=[sqrt(2beta/tc)+beta]t and gamma=1+beta*t. Integrate exp(eta_upper-t)*t^j with exact incomplete-gamma polynomials to obtain D_j/Q^2.',
        tail_outer='For p>=Q/2>=2P let T2 be the lower t(2P), verified >=4. For j<=2, t^j exp(-t/2) decreases for t>=T2. Hence f(t-eta)*t^j<=2exp(eta_upper-T2/2)*T2^j*f(t/2). The existing eta=0,beta<=0.012 upper-response coefficient C yields K_j/max(1,Q^2)<=K_j/Q^2, including the logarithmic cusp.',
        response_error='Combine the positive-series remainder, using v^17<=v, and both p-tail regions. Each raw response-jet error is <=delta_i*v. Then Y=R/X, X=Q^2/(B*Sref), and |Y_i-Y_N,i|<=d_i*v^2 with d_i=B*Sref*delta_i/QH^2. If A_i sums absolute upper polynomial coefficients, |Y_N,i|<=A_i*v^2.',
        quotient_error='For true Y,Y_N>=0, g=y/(1+y) has |g_prime|<=1, |g_second|<=2, |g_third|<=6. The value difference is d0*v^2; first partial a is <=da*v^2+2*A_a*d0*v^4; second ab is <=dab*v^2+2*A_ab*d0*v^4+2*(da*A_b+db*A_a+da*db)*v^4+6*A_a*A_b*d0*v^6. Integration gives divisors 3,7,11.',
        neumann='For L=4 use sum_{l=1}^L(-1)^(l+1)*Y_N^l and signed remainder Y_N^(L+1)/(1+Y_N). For positive y its value, first and second derivatives are bounded by y^(L+1),(L+1)y^L,L(L+1)y^(L-1). Chain-rule errors all have v^(2L+2), with integral divisor 4L+3. No assumption Y_N<1 is needed.',
        exact_integration='Each polynomial v^k integrates over QH..infinity to QH/(2k-1). Six jets are multiplied using normalized bivariate Taylor coefficients; second diagonal components are divided by two before multiplication and restored afterward. The product rule and remainder derivatives have exact symbolic controls.',
        physical='Apply the existing seven-by-six neutral matrix at the certified true eta root and multiply by -2B<Z^2>scale/beta. Cutoffs and coordinate normalizations are fixed in this decomposition. Directed interval arithmetic includes input intervals, moment enclosures, composition and constants.'))


def context():
    iv.prec=128
    state=dict(np.load(cusp.THERMO/'states.npz'));ions=dict(np.load(neutral.density.ionic.OUT/'inputs.npz'))
    roots=json.loads((cusp.ROOTS/'result.json').read_text())['rows'];data=dict(np.load(moments.OUT/'inputs.npz'))
    constants=json.loads((neutral.density.ionic.OUT/'constants.json').read_text())['native_binary64']
    return state,ions,roots,data,I(constants['alpha'])/iv.pi,I(json.loads((sector.OUT/'result.json').read_text())['uniform_real_response_upper_coefficient'])


def evaluate(row,ctx):
    state,ions,roots,data,B,C=ctx;i=row['position'];root=roots[i]
    assert row['cell']==root['cell']==int(state['cells'][i])==int(data['cells'][i])
    P=I(data['pcut'][i]);QH=4*P;beta=I(data['beta'][i]);Sref=I(data['Sref'][i]);scale=I(state['scale'][i]);eta=cusp.interval(root['root'])
    assert low(beta)>0 and high(2*beta)<=F('.012')
    mu=[[iv.mpf([I(a).a,I(b).b]) for a,b in values] for values in row['moments']]
    k=B*Sref/QH**2;factor=QH/scale;Y=[zeros(18) for _ in range(6)]
    for m in range(1,17):
        for j in range(6):Y[j][m+1]=k*(mu[m-1][j]/(2*m-1)+2*P*P*mu[m][j]/((2*m-1)*(2*m+1)))/16**m
    A=[I(sum(high(abs(x)) for x in values)) for values in Y]
    tc=I(64+max(F(0),high(eta)));cap=I(high(P*P/(beta*(iv.sqrt(1+P*P)+1))));T2=I(low(4*P*P/(beta*(iv.sqrt(1+4*P*P)+1))))
    assert low(cap)>=high(tc) and low(T2)>=4
    def integral(n):return iv.exp(I(high(eta))-tc)*sum(I(math.factorial(n))/math.factorial(j)*tc**j for j in range(n+1))
    D=[I(4)/3*beta*(iv.sqrt(2*beta/tc)+beta)*(integral(j+1)+2*beta*integral(j+2)+beta**2*integral(j+3))/Sref for j in range(3)]
    K=[2*iv.exp(I(high(eta))-T2/2)*T2**j*C/Sref for j in range(3)]
    tail=[D[0]+K[0],D[0]+K[0],D[0]+K[0],D[1]+K[1],D[1]+K[1],D[1]+D[2]+K[1]+K[2]]
    series=I(high(mu[16][0]))*(I(1)/33+2*P*P/(33*35))/(15*16**16);assert low(series)>0
    delta=[series*t+v/QH**2 for t,v in zip([I(1),I(1),I(1),cap,cap,cap*cap+cap],tail)]
    d=[k*x for x in delta]
    error=[d[0]/3,d[1]/3+2*A[1]*d[0]/7,I(0),d[3]/3+2*A[3]*d[0]/7,I(0),I(0)]
    for j,a,b in [(2,1,1),(4,1,3),(5,3,3)]:error[j]=d[j]/3+2*A[j]*d[0]/7+2*(d[a]*A[b]+d[b]*A[a]+d[a]*d[b])/7+6*A[a]*A[b]*d[0]/11
    L=4;rem=[A[0]**5,5*A[0]**4*A[1],I(0),5*A[0]**4*A[3],I(0),I(0)]
    for j,a,b in [(2,1,1),(4,1,3),(5,3,3)]:rem[j]=5*A[0]**4*A[j]+20*A[0]**3*A[a]*A[b]
    error=[factor*(x+y/(4*L+3)) for x,y in zip(error,rem)]
    normal=[a.copy() for a in Y]
    normal[2]=normal[2]/2;normal[5]=normal[5]/2
    power=normal;raw=[I(0)]*6
    for order in range(1,L+1):
        for j in range(6):raw[j]+=(-1)**(order+1)*factor*sum(power[j][n]/(2*n-1) for n in range(2,len(power[j])))
        if order<L:power=product(power,normal)
    raw[2]*=2;raw[5]*=2
    raw=[x+cusp.symmetric(high(err)) for x,err in zip(raw,error)]
    coeff=[cusp.interval(root['derivatives'][name]) for name in neutral.density.FIELDS]
    counts=[I(x)/I(a) for x,a in zip(ions['X'][i],ions['A'])];Z2=sum(x*I(z)**2 for x,z in zip(counts,ions['Z']))/sum(counts);pref=2*B*Z2*scale/beta
    fields=[-pref*sum(c*x for c,x in zip(line,raw)) for line in neutral.matrix(coeff)]
    texts=list(map(cusp.interval_text,fields));points=[];errors=[]
    for text in texts:
        lo,hi=cusp.endpoints(text);point=float((lo+hi)/2);points.append(point);errors.append(max(abs(F.from_float(point)-lo),abs(hi-F.from_float(point))))
    return dict(position=i,cell=root['cell'],QH_exact=str(low(QH)),raw_integrals=list(map(cusp.interval_text,raw)),
        raw_analytic_error_upper=[str(high(x)) for x in error],Y_coefficient_upper=[str(high(x)) for x in A],
        field_enclosures=texts,field_approximation=points,field_error_upper=list(map(str,errors)),passed=max(errors)<F('2e-7'))


def preflight():
    bindings();symbolic();save('preflight.json',evaluate(next(moments.records()),context()))


def block(rows):
    ctx=context();path=OUT/f"block-{rows[0]['position']:04d}.jsonl.gz";assert not path.exists()
    with gzip.open(path,'xt',encoding='utf-8',compresslevel=1) as stream:
        for row in rows:stream.write(json.dumps(evaluate(row,ctx))+'\n')
    return rows[0]['position']


def records():
    for path in sorted(OUT.glob('block-*.jsonl.gz')):
        with gzip.open(path,'rt',encoding='utf-8') as stream:
            for line in stream:yield json.loads(line)


def run():
    p=bindings();rows=list(moments.records());finished=[]
    with ProcessPoolExecutor(max_workers=p['processes']) as pool:
        work=[pool.submit(block,rows[i:i+p['block_size']]) for i in range(0,len(rows),p['block_size'])]
        for future in as_completed(work):
            finished.append(future.result());save('progress.json',dict(completed_blocks=sorted(finished)));print('HIGH Q SELF BLOCKS',len(finished),flush=True)
    maximum=[F(0)]*7;count=0
    for row in records():
        assert row['position']==count;count+=1;maximum=[max(a,F(b)) for a,b in zip(maximum,row['field_error_upper'])]
    assert count==p['cells']
    save('result.json',dict(classification='Proven',passed=max(maximum)<F(p['field_budget']),cells=count,fields=FIELDS,
        maximum_errors=dict(zip(FIELDS,map(str,maximum))),display_only=dict(zip(FIELDS,map(float,maximum))),
        high_Q_integral_certified=True,complete_outer_integral_certified=False,physical_EOS_certified=False))


def controls():
    from scipy.integrate import quad_vec
    p=bindings();state,ions,roots,data,B,C=context();checks=[]
    alpha=json.loads((neutral.density.ionic.OUT/'constants.json').read_text())['native_binary64']['alpha'];b=float(alpha)/np.pi
    for row in records():
        i=row['position']
        if i not in p['control_positions']:continue
        lo,hi=cusp.endpoints(roots[i]['root']);eta=float((lo+hi)/2);beta=float(data['beta'][i]);Sref=float(data['Sref'][i]);scale=float(state['scale'][i]);QH=4*float(data['pcut'][i]);k=b*Sref/QH**2
        def integrand(w):
            if w==0:return np.zeros(6)
            r=neutral.thermo.response(np.array([eta]),np.array([beta]),np.array([QH/w]),np.array([Sref]),2e-13).reshape(6)
            x=k*w*w;den=1+x*r[0]
            return np.array([r[0]/den,r[1]/den**2,r[2]/den**2-2*x*r[1]**2/den**3,r[3]/den**2,r[4]/den**2-2*x*r[1]*r[3]/den**3,r[5]/den**2-2*x*r[3]**2/den**3])
        value,estimate=quad_vec(integrand,0,1,epsabs=2e-13,epsrel=2e-10);value*=QH/scale*k
        target=np.array([float(sum(cusp.endpoints(x))/2) for x in row['raw_integrals']]);score=float(np.max(np.abs(target-value))/max(np.max(np.abs(value)),1e-300))
        checks.append(dict(cell=row['cell'],reference=value.tolist(),relative_vector_difference=score,adaptive_error_estimate=float(estimate),passed=score<float(p['control_relative_vector_gate'])))
        print('HIGH Q INDEPENDENT CONTROL',row['cell'],score,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=len(checks)==3 and all(x['passed'] for x in checks),checks=checks,
        scope='Independent binary64 nested adaptive quadrature at a rounded neutral-root midpoint; this finite comparison does not assert rigorous reference inclusion.'))


def finalize():
    bindings()
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    save('manifest.json',dict(sha256={path.relative_to(ROOT).as_posix():sha(path) for path in OUT.iterdir() if path.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    count=0;maximum=[F(0)]*7
    for row in records():
        assert row['position']==count and row['passed'];count+=1
        for j,(text,point,bound) in enumerate(zip(row['field_enclosures'],row['field_approximation'],row['field_error_upper'],strict=True)):
            lo,hi=cusp.endpoints(text);error=max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)));assert error==F(bound)<F(p['field_budget']);maximum[j]=max(maximum[j],error)
    r=json.loads((OUT/'result.json').read_text());assert count==p['cells']==3206 and r['maximum_errors']==dict(zip(FIELDS,map(str,maximum)))
    assert next(records())==json.loads((OUT/'preflight.json').read_text())
    print('PASS high-Q nonlinear self integral and seven neutral thermodynamic fields; intermediate Q/full EOS remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
