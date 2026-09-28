"""Exact damping-residue audit of the declared frozen nonzero-heat system."""
from fractions import Fraction as F
import gzip,json,sys,urllib.request
import numpy as np
import sympy as sp
import gr_heat_entropy_closure as prior

g=prior.g;OUT=g.OUT/'gr-heat-rank-one-damping'
URL='https://people.maths.ox.ac.uk/chengq/preprints/cll/cll.pdf'


def save(name,value):(OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    prior.verify();assert not OUT.exists();OUT.mkdir()
    with urllib.request.urlopen(URL,timeout=60) as response:(OUT/'chen-levermore-liu1994.pdf').write_bytes(response.read())
    assert (OUT/'chen-levermore-liu1994.pdf').read_bytes().startswith(b'%PDF')
    files=[g.ROOT/'verification/gr_heat_rank_one_damping.py',OUT/'chen-levermore-liu1994.pdf',prior.OUT/'manifest.json',prior.OUT/'root-brackets.json',
        prior.OUT/'initial-rates.npz',prior.OUT/'constant-times.json',prior.old.char.OUT/'coefficients.npz']
    save('plan.json',dict(classification='Conjectural',checkpoint='a7ab5430',bindings={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in files},
        cells=5735,models=2,maximum_bisections=256,
        declared_linear_model='M*y_t+c*N*y_x+E44*y/tau=0 using the exact frozen nonzero-heat entropy-closure M,N and original two constant proper tau values. A constant external forcing balances the background heat source. This is a declared constant-coefficient diagnostic, not the full linearization about the evolving inhomogeneous star.',
        omitted_terms='Background gradients, perturbations of their coefficients, geometry and reactions, metric feedback and boundary/atmosphere conditions are absent. Do not transfer a pass/failure to the complete GR system. Physical heat times and EOS derivatives remain uncertified.',
        target='Test the signs d_i=q3(c_i)/p4_prime(c_i), using only exact rational coefficients and certified simple real characteristic-root brackets. All d_i>0 imply strict damping for every finite nonzero real k. A negative d_i implies a growing branch for sufficiently large |k| in this frozen system. Zero/undecided cases remain unresolved.',
        method='Reuse each original causal root bracket. Exact interval Horner bounds for q3 and p4_prime, with exact bisection of p4 when required. No floating eigenvalues or selected finite k grid establish the verdict.',
        source=dict(classification='Imported from prior work',url=URL,doi='10.1002/cpa.3160470602',
            scope='Context for entropy, subcharacteristic interlacing and relaxation stability. The rank-one finite-dimensional proof and application below are derived explicitly and do not assume the paper supplies nonlinear GR closure.')))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return p


def polynomial(coeff,x):
    value=F(0)
    for a in coeff:value=value*x+a
    return value


def interval_polynomial(coeff,a,b):
    lo=hi=F(0)
    for value in coeff:
        candidates=[lo*a,lo*b,hi*a,hi*b];lo=min(candidates)+value;hi=max(candidates)+value
    return lo,hi


def symbolic():
    b,r,d,e,j,a,kr,kt,z,eps=sp.symbols('b r d e j a kr kt z eps',real=True)
    M=sp.Matrix([[1,0,0,0],[0,b,2*j,0],[0,0,1,1],[-j*kr/2,-j*kt/2,a,1]])
    N=sp.Matrix([[0,0,1,0],[0,0,e,1],[r,d,2*j,0],[0,a,0,0]]);E=sp.diag(0,0,0,1)
    p4=(z*M-N).det();q3=z*(b*z*z+2*j*(d-b)*z-b*r-d*e)
    assert sp.expand((z*M-N+eps*E).det()-p4-eps*q3)==0
    u=sp.Matrix(sp.symbols('u:4',nonzero=True,real=True));v=sp.Matrix(sp.symbols('v:4',nonzero=True,real=True));H=sp.diag(*[v[i]/u[i] for i in range(4)])
    assert H*(u*v.T)==v*v.T
    s,k=sp.symbols('s k',real=True);A=sp.diag(-1,1);good=sp.ones(2)/2;bad=sp.Matrix([[-1,2],[-1,2]])
    assert sp.expand((s*sp.eye(2)+sp.I*k*A+good).det()-(s*s+s+k*k))==0
    badpoly=sp.expand((s*sp.eye(2)+sp.I*k*A+bad).det());kk=sp.sqrt(sp.Rational(3,5));ss=sp.Rational(1,2)+sp.I*3*kk/2
    assert sp.simplify(badpoly.subs({s:ss,k:kk}))==0 and bad.det()==0
    # Small runnable check of the interval Horner primitive, with mixed signs.
    for x in [F(-2),F(-1),F(0),F(1)]:
        lo,hi=interval_polynomial([F(2),F(-3),F(1)],F(-2),F(1));assert lo<=polynomial([F(2),F(-3),F(1)],x)<=hi
    save('symbolic.json',dict(classification='Proven',passed=True,
        determinant='p4(z)=det(z*M-N), q3(z)=cofactor44(z*M-N)=z*(b*z^2+2*j*(d-b)*z-b*r-d*e). det(z*M-N+epsilon*E44)=p4(z)+epsilon*q3(z). Verified symbolically.',
        residue='A=M^-1*N, B=M^-1*E44=u*v^T. In a real eigenbasis of A with distinct eigenvalues c_i, the diagonal entries of B are d_i=q3(c_i)/p4_prime(c_i). This follows by comparing tr((zI-A)^-1 B) with the determinant lemma q3(z)/p4(z). Time is scaled by the same positive tau; the physical damping rate is d_i/tau.',
        stable_proof='If all d_i=ubar_i*vbar_i>0, every ubar_i and vbar_i is nonzero and H=diag(vbar_i/ubar_i)>0. H*diag(c_i) is symmetric and H*Bbar=vbar*vbar^T>=0. For an eigenpair of -ik*diag(c_i)-Bbar, Re(s)*(z^*H z)=-|vbar^T z|^2<=0. Equality at k!=0 would make z a single characteristic eigenvector orthogonal to vbar, impossible because all vbar_i are nonzero. Thus every mode strictly decays for every finite nonzero k.',
        unstable_proof='Simple real eigenvalues of A continue analytically under epsilon*B. Their first corrections are epsilon*d_i. Therefore the eigenvalue branch of -ik*A-B has s_i(k)=-ik*c_i-d_i+O(1/k). If d_i<0 its real part is positive for all sufficiently large |k|. This is an asymptotic existence statement, not a chosen physical wavenumber.',
        interlacing='Positive residues imply the three zeros of q3 strictly interlace the four real zeros of p4. This follows from alternating q3(c_i) signs and its degree. The direct residue/symmetrizer proof above establishes the damping criterion without assuming a physical entropy Hessian.',
        controls='A=diag(-1,1), B=ones(2)/2 gives s^2+s+k^2 and strict decay for k!=0. The rank-one B=[[-1,2],[-1,2]] has a negative first residue and the exact growing eigenvalue s=1/2+3*i*sqrt(3/5)/2 at k=sqrt(3/5). Hence real causal characteristics alone do not imply damping.',
        boundary='No claim about omitted background lower-order terms, finite nonlinear time evolution, physical transport calibration or native EOS errors.'))


def run():
    p=bindings();symbolic();fields=np.load(prior.old.char.OUT/'coefficients.npz')['coefficients'];rates=dict(np.load(prior.OUT/'initial-rates.npz'))
    times=list(map(F,json.loads((prior.OUT/'constant-times.json').read_text())['exact_proper_seconds']));brackets=json.loads((prior.OUT/'root-brackets.json').read_text())['rows']
    counts=dict(stable=0,unstable=0,unresolved=0);examples={};max_refinements=0
    with gzip.open(OUT/'certificates.jsonl.gz','xt',compresslevel=1) as stream:
        for record in brackets:
            i,model=record['cell'],record['model'];B,R,D,E,H=map(lambda x:F(float(x)),fields[i]);J,KR,KT=[F(float(rates[key][i])) for key in ['j','kr','kt']];A=H/times[model]
            p4=[B*(1-A)-J*J*KT,J*(2*(D-B)-2*A-(1-E)*KT/2+B*KR/2),-(B*R+D*E+A*(1-E-D))+J*J*KT,J*(2*A+R*KT/2-D*KR/2),R*A]
            q3=[B,2*J*(D-B),-(B*R+D*E),F(0)];derivative=[(4-k)*x for k,x in enumerate(p4[:-1])];assert p4[0]>0
            roots=[];signs=[]
            for index,(lo,hi) in enumerate(record['bounds']):
                lo,hi=F(lo),F(hi);pl=polynomial(p4,lo);ph=polynomial(p4,hi);assert -1<lo<hi<1 and pl*ph<0
                for refinement in range(p['maximum_bisections']+1):
                    ql,qh=interval_polynomial(q3,lo,hi);dl,dh=interval_polynomial(derivative,lo,hi)
                    if (ql>0 or qh<0) and (dl>0 or dh<0):
                        vals=[ql/dl,ql/dh,qh/dl,qh/dh];rl,rh=min(vals),max(vals);sign=1 if rl>0 else -1;break
                    if refinement==p['maximum_bisections']:sign=0;rl=rh=None;break
                    mid=(lo+hi)/2;pm=polynomial(p4,mid)
                    if pm==0:lo=hi=mid
                    elif pl*pm<0:hi=mid;ph=pm
                    else:lo=mid;pl=pm
                max_refinements=max(max_refinements,refinement);signs.append(sign)
                roots.append(dict(bounds=[str(lo),str(hi)],q3_interval=[str(ql),str(qh)],derivative_interval=[str(dl),str(dh)],
                    residue_interval=None if rl is None else [str(rl),str(rh)],sign=sign,bisections=refinement))
            verdict='unstable' if -1 in signs else 'stable' if all(signs) else 'unresolved';counts[verdict]+=1
            row=dict(cell=i,model=model,verdict=verdict,roots=roots);stream.write(json.dumps(row)+'\n');examples.setdefault(verdict,row)
            if sum(counts.values())%512==0:print('HEAT DAMPING',counts,flush=True)
    save('result.json',dict(classification='Proven',audit_complete=True,cells=p['cells'],models=p['models'],counts=counts,maximum_bisections_used=max_refinements,examples=examples,
        all_frozen_modes_stable=counts['stable']==p['cells']*p['models'],all_classified=counts['unresolved']==0,
        full_nonequilibrium_GR_stability_certified=False,physical_transport_calibrated=False,physical_EOS_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in OUT.iterdir() if x.is_file()}));verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed'];r=json.loads((OUT/'result.json').read_text());assert r['audit_complete'] and sum(r['counts'].values())==p['cells']*p['models']
    with gzip.open(OUT/'certificates.jsonl.gz','rt') as stream:assert sum(1 for _ in stream)==p['cells']*p['models']
    print('PASS exact frozen rank-one damping audit; inspect stable/unstable/unresolved counts:',r['counts'],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
