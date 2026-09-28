"""Certify actual response integrals around the logarithmic momentum cusp."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
import gzip,json,math,sys
import numpy as np
import mpmath as mp
from mpmath import iv
import sympy as sp
import gr_logarithmic_gauss_rule as rule
from interval_records import interval_text

ROOT=rule.ROOT;OUT=ROOT/'outputs/direct-eos-gr33/gr-logarithmic-response'
THERMO=OUT.parent/'gr-polarization-thermodynamics';ROOTS=OUT.parent/'gr-electron-density-refined-export'
PARAMETERS=['value','eta','etaeta','tau','etatau','tautau']
low=rule.low;high=rule.high


def I(x):
    if isinstance(x,(float,np.floating)):return iv.mpf(float(x))
    return rule.I(x)
def endpoints(text):return tuple(map(F,text[1:-1].split(',')))
def interval(text):
    lo,hi=endpoints(text);assert lo<=hi
    return iv.mpf([I(lo).a,I(hi).b])
def symmetric(bound):return iv.mpf([I(-bound).a,I(bound).b])
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    rule.verify();assert not OUT.exists();OUT.mkdir()
    roots=json.loads((ROOTS/'result.json').read_text());assert roots['passed'] and roots['cells']==3206
    a=dict(np.load(THERMO/'states.npz'));z=np.array([.001,.01,.03,.1,.3,1.,3.,10.,30.,100.,300.,1000.])
    assert np.array_equal(a['cells'],[r['cell'] for r in roots['rows']])
    np.savez_compressed(OUT/'inputs.npz',cells=a['cells'],beta=a['beta'],Sref=a['Sref'],scale=a['scale'],
        z=z,Q=a['scale'][:,None]*z[None,:],root_intervals=np.array([r['root'] for r in roots['rows']]))
    files=[ROOT/'verification/gr_logarithmic_response_certificate.py',rule.OUT/'manifest.json',rule.OUT/'result.json',
        ROOTS/'manifest.json',ROOTS/'result.json',THERMO/'states.npz',OUT/'inputs.npz']
    save('plan.json',dict(classification='Proven',checkpoint='fa17efb',bits=128,cells=3206,points_per_cell=12,
        bindings={p.relative_to(ROOT).as_posix():rule.sha(p) for p in files},parameter_order=PARAMETERS,
        processes=3,block_size=128,initial_radius='1/8',maximum_radius_halvings=40,phase_bound='1',window_radius_ratio=16,
        per_component_absolute_error_budget='1e-11',bound_exponent_floor=-128,
        control_positions=[0,1603,3205],control_z_indices=[0,5,8,11],control_digits=70,
        target='The complete normalized momentum contribution of u in [1-h,1+h], for each actual certified neutral root and each specified finite Q. Includes all six eta/tau response derivatives, with tau=ln(T).',
        analytic_domain='On |u-1|<=R, Re gamma >= sqrt(1+Q^2*(1-2R))=gmin. Require |Im t|<=Q^2*(1+R)*R/(beta*gmin)<=1<pi/2. Then Re exp(t-eta)>0, so q and 1-q have modulus at most one, |1-2q|<=1, and |q|<=min(1,exp(eta_max-tmin)).',
        method='Choose R only from the fixed analyticity rule, h=R/16. Bound the eighth u derivatives by Cauchy on the disk. Apply the certified ordinary/logarithmic Gauss rules, including interval nodes, weights and evaluation rounding. Midpoints and error budgets are computed from the exact encoded interval endpoints.',
        small_terms='Only the analytic upper envelope of q is enlarged to at least exp(-128). If a proved bound on the whole window is below the same budget, return a zero approximation with that full signed enclosure. No floating underflow is treated as an exact zero.',
        scope='The cusp window and six partial derivatives in the declared full-ionization ideal-electron medium. Other momentum regions, the outer wavenumber integral, full correlated EOS, GR and observations remain open.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert rule.sha(ROOT/rel)==digest,rel
    return p


def symbolic():
    Q,u,G,L=sp.symbols('Q u G L',positive=True);p=Q*u
    A=Q**3*u**2/G;B=Q*u*(1+Q**2*(u**2-1))/(2*G)
    Kp=p*p/G+p*(1+p*p-Q*Q)*L/(2*Q*G)
    assert sp.simplify(Q*Kp-A-B*L)==0
    x,y=sp.symbols('x y',real=True);den=(1+x)**2+y*y
    assert sp.expand(den-x*x-y*y-(1+2*x))==0
    assert sp.expand(den-(x-1)**2-y*y-4*x)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        decomposition='K_p dp/Sref=[A+B*(log(1+u)-log|1-u|)]du, A=Q^3*u^2/(gamma*Sref), B=Q*u*(1+Q^2*(u^2-1))/(2*gamma*Sref). Multiply by the same occupation derivative for each field.',
        complex_bounds='For Re E>=0, E=exp(t-eta), |1/(1+E)|<=min(1,1/|E|), |E/(1+E)|<=1 and |(E-1)/(E+1)|<=1. The six occupation bounds are qmax*[1,1,1,tmax,tmax,tmax^2+tmax].',
        gamma='For |u-1|<=R<=1/8, Re(1+Q^2*u^2)>=1+Q^2*(1-2R)>0. The principal gamma is analytic, Re gamma>=gmin; Im gamma=Q^2*Re(u)*Im(u)/Re(gamma). tmin=Q^2*(1-2R)/(beta*(gmin+1)); |t|<=Q^2*(1+R)^2/(beta*(gmin+1)).',
        prefactors='|A|<=Q^3*(1+R)^2/(gmin*Sref), |B|<=Q*(1+R)*(1+Q^2*R*(2+R))/(2*gmin*Sref), |log(1+u)|<=log(2)-log(1-R/2).',
        Cauchy='For each real u in [1-h,1+h], a disk of radius R-h lies in the analytic R disk. Hence M8<=8!*M/(R-h)^8. The total two-sided window remainder is <=2*h^9*[norm_unit*M_smooth+(norm_log+abs(log(h))*norm_unit)*M_B]/(R-h)^8.',
        omission='The whole window modulus is <=2*h*(M_smooth+(1-log(h))*M_B). This permits a signed enclosure around zero only when that full bound satisfies the fixed budget.',
        boundary='This proves a local momentum-window contribution at true certified eta, not a full response or EOS certificate.'))


def initialize():
    global PLAN,DATA,RULES
    PLAN=json.loads((OUT/'plan.json').read_text());iv.prec=PLAN['bits'];DATA=dict(np.load(OUT/'inputs.npz'))
    rows=json.loads((rule.OUT/'result.json').read_text())['rules']
    RULES=[dict(norm=I(row['norm_exact']),nodes=[(interval(n['node']),interval(n['weight'])) for n in row['nodes']]) for row in rows]


def kernel(u,eta,beta,Q,Sref,ctx):
    gamma=ctx.sqrt(1+Q*Q*u*u);t=Q*Q*u*u/(beta*(gamma+1));q=1/(1+ctx.exp(t-eta));v=1-q
    factors=[q,q*v,q*v*(1-2*q),t*q*v,t*q*v*(1-2*q),(t*t*(1-2*q)-t)*q*v]
    A=Q**3*u*u/(gamma*Sref);B=Q*u*(1+Q*Q*(u*u-1))/(2*gamma*Sref);smooth=A+B*ctx.log(1+u)
    return [smooth*f for f in factors],[B*f for f in factors]


def certified_window(i,j):
    eta=interval(str(DATA['root_intervals'][i]));beta=I(DATA['beta'][i]);Q=I(DATA['Q'][i,j]);Sref=I(DATA['Sref'][i])
    radius=F(PLAN['initial_radius'])
    for _ in range(PLAN['maximum_radius_halvings']+1):
        R=I(radius);gmin=iv.sqrt(1+Q*Q*(1-2*R));phase=Q*Q*(1+R)*R/(beta*gmin)
        if high(phase)<=F(PLAN['phase_bound']):break
        radius/=2
    else:raise AssertionError(('analytic disk',i,j))
    h=R/PLAN['window_radius_ratio'];tmin=Q*Q*(1-2*R)/(beta*(gmin+1));tmax=Q*Q*(1+R)**2/(beta*(gmin+1))
    exponent=max(F(PLAN['bound_exponent_floor']),high(eta-tmin))
    qmax=I(1) if exponent>=0 else iv.exp(I(exponent))
    factors=[I(1),I(1),I(1),tmax,tmax,tmax*tmax+tmax]
    A=Q**3*(1+R)**2/(gmin*Sref)*qmax;B=Q*(1+R)*(1+Q*Q*R*(2+R))/(2*gmin*Sref)*qmax
    Mb=[B*f for f in factors];Ms=[(A+B*(iv.log(2)-iv.log(1-R/2)))*f for f in factors];logh=iv.log(h)
    full=[2*h*(s+(1-logh)*b) for s,b in zip(Ms,Mb)]
    omitted=max(map(high,full))<F(PLAN['per_component_absolute_error_budget'])
    if omitted:enclosures=[symmetric(high(v)) for v in full]
    else:
        total=[iv.mpf(0) for _ in PARAMETERS]
        for sign in [-1,1]:
            for node,weight in RULES[0]['nodes']:
                smooth,b=kernel(1+sign*h*node,eta,beta,Q,Sref,iv)
                total=[old+h*weight*(s-logh*v) for old,s,v in zip(total,smooth,b)]
            for node,weight in RULES[1]['nodes']:
                _,b=kernel(1+sign*h*node,eta,beta,Q,Sref,iv)
                total=[old+h*weight*v for old,v in zip(total,b)]
        errors=[2*h**9*(RULES[0]['norm']*s+(RULES[1]['norm']-logh*RULES[0]['norm'])*b)/(R-h)**8 for s,b in zip(Ms,Mb)]
        enclosures=[value+symmetric(high(error)) for value,error in zip(total,errors)]
    texts=list(map(interval_text,enclosures));points=[];errors=[]
    for text in texts:
        lo,hi=endpoints(text);point=float((lo+hi)/2);points.append(point)
        errors.append(max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point))))
    score=max(errors)
    return dict(position=i,cell=int(DATA['cells'][i]),z_index=j,Q_hex=float(DATA['Q'][i,j]).hex(),
        cauchy_radius=str(radius),window_half_width=str(radius/PLAN['window_radius_ratio']),phase_upper=str(high(phase)),
        maximum_B_eighth_derivative_upper=str(max(high(math.factorial(8)*b/(R-h)**8) for b in Mb)),
        maximum_smooth_eighth_derivative_upper=str(max(high(math.factorial(8)*s/(R-h)**8) for s in Ms)),
        bounded_omission=omitted,enclosures=texts,approximation=points,maximum_error_upper=str(score),
        passed=score<F(PLAN['per_component_absolute_error_budget']))


def block(first,count):
    path=OUT/f'block-{first:04d}.jsonl.gz'
    with gzip.open(path,'xt',encoding='utf-8',compresslevel=1) as stream:
        for i in range(first,first+count):
            for j in range(PLAN['points_per_cell']):stream.write(json.dumps(certified_window(i,j),separators=(',',':'))+'\n')
    return path.name


def records():
    for path in sorted(OUT.glob('block-*.jsonl.gz')):
        with gzip.open(path,'rt',encoding='utf-8') as stream:
            for line in stream:yield json.loads(line)


def run():
    p=bindings();symbolic();labels=[]
    with ProcessPoolExecutor(max_workers=p['processes'],initializer=initialize) as pool:
        work=[pool.submit(block,i,min(p['block_size'],p['cells']-i)) for i in range(0,p['cells'],p['block_size'])]
        for done in as_completed(work):
            labels.append(done.result());save('progress.json',dict(completed_blocks=sorted(labels)));print('CERTIFIED CUSP BLOCKS',len(labels),flush=True)
    rows=list(records());assert [(r['position'],r['z_index']) for r in rows]==[(i,j) for i in range(p['cells']) for j in range(p['points_per_cell'])]
    maximum=max(F(row['maximum_error_upper']) for row in rows)
    save('result.json',dict(classification='Proven',passed=all(row['passed'] for row in rows),cases=len(rows),
        bounded_omissions=sum(row['bounded_omission'] for row in rows),maximum_encoded_error_exact=str(maximum),
        display_only=dict(maximum_error=float(maximum),minimum_radius=float(min(F(row['cauchy_radius']) for row in rows))),
        other_momentum_regions_certified=False,outer_wavenumber_integral_certified=False,full_physical_EOS_certified=False))
    print('CERTIFIED CUSP RESULT',float(maximum),flush=True)


def controls():
    p=bindings();initialize();mp.mp.dps=p['control_digits'];selected={(i,j) for i in p['control_positions'] for j in p['control_z_indices']}
    chosen={(r['position'],r['z_index']):r for r in records() if (r['position'],r['z_index']) in selected};checks=[]
    def M(f):return mp.mpf(f.numerator)/f.denominator
    for (i,j),row in chosen.items():
        lo,hi=endpoints(str(DATA['root_intervals'][i]));eta=(M(lo)+M(hi))/2
        beta,Q,Sref=map(lambda x:mp.mpf(float(x)),[DATA['beta'][i],DATA['Q'][i,j],DATA['Sref'][i]])
        h=M(F(row['window_half_width']));_,B0=kernel(mp.mpf(1),eta,beta,Q,Sref,mp);values=[]
        for k in range(6):
            def integrand(v):
                value=mp.mpf(0)
                for sign in [-1,1]:
                    smooth,b=kernel(1+sign*h*v,eta,beta,Q,Sref,mp)
                    value+=smooth[k]-mp.log(h)*b[k]-(mp.log(v)*(b[k]-B0[k]) if v else 0)
                return value
            values.append(h*(2*B0[k]+mp.quad(integrand,[0,mp.mpf('.5'),1])))
        inside=all(M(endpoints(text)[0])<=value<=M(endpoints(text)[1]) for text,value in zip(row['enclosures'],values))
        checks.append(dict(cell=row['cell'],z_index=j,passed=inside,reference=list(map(str,values))))
        print('CUSP INDEPENDENT CONTROL',row['cell'],j,inside,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=all(c['passed'] for c in checks),controls=checks));assert all(c['passed'] for c in checks)


def finalize():
    bindings()
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():rule.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert rule.sha(ROOT/rel)==digest,rel
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    r=json.loads((OUT/'result.json').read_text());assert r['cases']==3206*12
    print('PASS actual certified-root logarithmic momentum windows and six partial derivatives; remaining regions/outer integral open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
