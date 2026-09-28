"""Rigorous response tails outside the already certified logarithmic windows."""
from fractions import Fraction as F
import gzip,json,math,sys
import numpy as np
import mpmath as mp
from mpmath import iv
import sympy as sp
import gr_logarithmic_response_certificate as cusp

ROOT=cusp.ROOT;OUT=cusp.OUT.parent/'gr-response-momentum-tail'
I=cusp.I;low=cusp.low;high=cusp.high;sha=cusp.rule.sha
def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    cusp.verify();assert not OUT.exists();OUT.mkdir()
    files=[ROOT/'verification/gr_response_momentum_tail.py',cusp.OUT/'manifest.json',cusp.OUT/'inputs.npz']
    save('plan.json',dict(classification='Proven',checkpoint='c85063d',bits=128,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        cells=3206,points_per_cell=12,parameter_order=cusp.PARAMETERS,cutoff_offset=64,
        per_component_absolute_error_budget='1e-11',control_positions=[0,1603,3205],control_z_indices=[0,5,8,11],control_digits=70,
        target='Integral over p >= the saved upward-rounded cutoff, excluding the complete cusp window. All six normalized response partials at the certified neutral eta interval.',
        cutoff='tc=64+max(eta_upper,0), p_theory=sqrt(beta*tc*(2+beta*tc)). Save binary64 pcut >= the certified upper endpoint of p_theory. The finite integrator must use exactly that saved binary cutoff.',
        proof='Outside |p/Q-1|<h, the logarithm is at most L=log((2+h)/h). The absolute kernel is at most p+(1+p^2+Q^2)*L/(2*Q). Substitute p^2=beta*t*(2+beta*t), dp=beta*(1+beta*t)/p dt and beta/p<=sqrt(beta/(2*tc)). Bound occupations by exp(eta_upper-t)*[1,1,1,t,t,t^2+t]. Integrate the resulting positive polynomial over all t>=tc, which enlarges the excluded-window domain.',
        derivative_scope='Differentiate the original fixed-p,Q integrand first, then split the resulting integrals. No derivative of this computational cutoff is a physical term.',
        boundary='Conditional ideal-electron RPA momentum tails only. Finite smooth momentum intervals, continuous outer-Q integration and the complete physical EOS remain open.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return p


def symbolic():
    b,t,Q=sp.symbols('b t Q',positive=True)
    assert sp.expand((1+b*t)*(1+b*t*(2+b*t)+Q**2)-(1+3*b*t+3*b*b*t*t+b**3*t**3+Q**2*(1+b*t)))==0
    for n in range(6):
        polynomial=sum(sp.factorial(n)/sp.factorial(k)*t**k for k in range(n+1))
        assert sp.simplify(-sp.diff(sp.exp(-t)*polynomial,t)-t**n*sp.exp(-t))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        moments='I_n=exp(eta_upper-tc)*sum_{k=0}^n (n!/k!)*tc^k = integral_tc^infinity exp(eta_upper-t)*t^n dt.',
        envelope='E_j=beta*(I_j+beta*I_(j+1))+L/(2Q)*sqrt(beta/(2tc))*(I_j+3beta*I_(j+1)+3beta^2*I_(j+2)+beta^3*I_(j+3)+Q^2*(I_j+beta*I_(j+1))).',
        partial_bounds='[E0,E0,E0,E1,E1,E2+E1]/Sref',
        domain='The logarithm bound applies only outside the cusp. Only the positive polynomial envelope is integrated over the enlarged full tail.'))


def records():
    with gzip.open(OUT/'bounds.jsonl.gz','rt',encoding='utf-8') as stream:
        for line in stream:yield json.loads(line)


def run():
    plan=bindings();symbolic();iv.prec=plan['bits'];data=dict(np.load(cusp.OUT/'inputs.npz'));previous=None;maximum=F(0);count=0
    with gzip.open(OUT/'bounds.jsonl.gz','xt',encoding='utf-8',compresslevel=1) as stream:
        for row in cusp.records():
            i,j=row['position'],row['z_index']
            if i!=previous:
                b=I(data['beta'][i]);Sref=I(data['Sref'][i]);eta_hi=cusp.endpoints(str(data['root_intervals'][i]))[1]
                tc=plan['cutoff_offset']+max(eta_hi,F(0));T=I(tc);theory=iv.sqrt(b*T*(2+b*T));pcut=float(high(theory))
                if F.from_float(pcut)<high(theory):pcut=float(np.nextafter(pcut,np.inf))
                assert F.from_float(pcut)>=high(theory)
                moments=[iv.exp(I(eta_hi-tc))*sum(I(F(math.factorial(n),math.factorial(k)))*T**k for k in range(n+1)) for n in range(6)]
                previous=i
            Q=I(data['Q'][i,j]);h=I(row['window_half_width']);L=iv.log((2+h)/h);E=[]
            for k in range(3):
                E.append((b*(moments[k]+b*moments[k+1])+L/(2*Q)*iv.sqrt(b/(2*T))*(moments[k]+3*b*moments[k+1]+3*b*b*moments[k+2]+b**3*moments[k+3]+Q*Q*(moments[k]+b*moments[k+1])))/Sref)
            bounds=list(map(high,[E[0],E[0],E[0],E[1],E[1],E[2]+E[1]]));score=max(bounds);maximum=max(maximum,score);count+=1
            stream.write(json.dumps(dict(position=i,cell=row['cell'],z_index=j,Q_hex=row['Q_hex'],window_half_width=row['window_half_width'],
                tc_exact=str(tc),pcut=pcut,pcut_hex=pcut.hex(),pcut_theory_upper=str(high(theory)),upper_bounds=list(map(str,bounds)),
                passed=score<F(plan['per_component_absolute_error_budget'])),separators=(',',':'))+'\n')
    assert count==plan['cells']*plan['points_per_cell']
    save('result.json',dict(classification='Proven',passed=maximum<F(plan['per_component_absolute_error_budget']),cases=count,
        maximum_upper_exact=str(maximum),display_only=dict(maximum_upper=float(maximum)),finite_smooth_momentum_certified=False,outer_Q_certified=False,physical_EOS_certified=False))
    print('MOMENTUM TAIL CERTIFICATE',count,float(maximum),flush=True)


def controls():
    plan=bindings();mp.mp.dps=plan['control_digits'];data=dict(np.load(cusp.OUT/'inputs.npz'));checks=[]
    selected={(i,j) for i in plan['control_positions'] for j in plan['control_z_indices']}
    def M(x):
        q=F(x);return mp.mpf(q.numerator)/q.denominator
    for row in records():
        i,j=row['position'],row['z_index']
        if (i,j) not in selected:continue
        beta,Q,Sref=map(lambda x:mp.mpf(float(x)),[data['beta'][i],data['Q'][i,j],data['Sref'][i]])
        eta=sum(map(M,cusp.endpoints(str(data['root_intervals'][i]))))/2;h=M(row['window_half_width']);pcut=mp.mpf(row['pcut'])
        kinetic=lambda p:p*p/(beta*(mp.sqrt(1+p*p)+1))
        start=kinetic(pcut);left=kinetic(Q*(1-h));right=kinetic(Q*(1+h));domains=[]
        if start<left:domains.append((start,left))
        domains.append((max(start,right),mp.inf));values=[]
        for k in range(6):
            def integrand(t):
                p=mp.sqrt(beta*t*(2+beta*t));u=p/Q
                smooth,B=cusp.kernel(u,eta,beta,Q,Sref,mp)
                return (smooth[k]-B[k]*mp.log(abs(1-u)))*beta*(1+beta*t)/(p*Q)
            value=mp.mpf(0)
            for a,b in domains:
                cuts=[a]+[a+d for d in [1,5,20,80] if a+d<b]+[b]
                value+=mp.quad(integrand,cuts)
            values.append(value)
        passed=all(abs(value)<=M(bound) for value,bound in zip(values,row['upper_bounds']))
        checks.append(dict(cell=row['cell'],z_index=j,passed=passed,reference=list(map(str,values)),
            maximum_fraction_of_bound=str(max(abs(v)/M(b) for v,b in zip(values,row['upper_bounds'])))))
        print('TAIL INDEPENDENT CONTROL',row['cell'],j,passed,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=all(c['passed'] for c in checks),controls=checks));assert len(checks)==12 and all(c['passed'] for c in checks)


def finalize():
    bindings()
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    plan=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    for name in ['symbolic.json','result.json','controls.json']:assert json.loads((OUT/name).read_text())['passed'],name
    count=0;maximum=F(0)
    for row in records():
        assert (row['position'],row['z_index'])==divmod(count,plan['points_per_cell'])
        assert row['pcut'].hex()==row['pcut_hex'] and F.from_float(row['pcut'])>=F(row['pcut_theory_upper'])
        bounds=list(map(F,row['upper_bounds']));assert len(bounds)==6 and min(bounds)>0
        score=max(bounds);assert row['passed'] and score<F(plan['per_component_absolute_error_budget']);maximum=max(maximum,score);count+=1
    result=json.loads((OUT/'result.json').read_text());assert count==result['cases']==3206*12 and maximum==F(result['maximum_upper_exact'])
    print('PASS exact exported momentum-tail bounds and upward binary cutoffs; finite smooth intervals/outer Q/full EOS remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
