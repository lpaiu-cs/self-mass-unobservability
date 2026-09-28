"""Conditional mixture WK coefficient and certified BC2022 fit derivatives.

Proven: the certificate concerns the declared fit and fixed binary inputs.
Conjectural: physical many-body remainders, screening and ionization closure.
"""
from fractions import Fraction as Q
import json, shutil, sys
import numpy as np
import sympy as sp
import mpmath as mp
from mpmath import iv
from interval_records import exact_endpoint, interval_text
import gr_dense_plasma as d

g=d.g;OUT=g.OUT/'gr-quantum-mixture-certificate'
COEFF=[Q(1,24),-Q(1,2880),Q(1,181440),-Q(1,9676800),Q(1,479001600),-Q(691,15692092416000)]
WK=[Q(1),Q(4),Q(8,3),Q(3),Q(12),Q(8),Q(10)]
TAIL=[Q(1),Q(14),Q(28,3),Q(13),Q(182),Q(364,3),Q(875,9)]


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def rat(v): return iv.mpf(v.numerator)/v.denominator
def upper(v): return exact_endpoint(v._mpi_[1])
def lower(v): return exact_endpoint(v._mpi_[0])
def ends(v): return [Q(s.strip()) for s in v[1:-1].split(',')]


def prepare():
    assert not OUT.exists();OUT.mkdir()
    for name in ['wigner-kirkwood-samaj2007.pdf','wigner-kirkwood-samaj2007.txt']:
        shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    stored=dict(np.load(d.OUT/'stellar-comparison.npz'));state,_=d.state_data()
    cells=stored['cells'];w=[];rs=[];ths=[];gammas=[];plasma=[]
    for k,i in enumerate(cells):
        count=state['X'][i]/g.c.A;x=count/count.sum();z=g.c.Z
        rse,ge=stored['parameters'][k,:2]
        parameters=[]
        for j in range(len(z)):
            gamma=ge*z[j]**(5/3)
            theta=gamma/np.sqrt(rse)*np.sqrt(3/(1822.88848*g.c.A[j]))/z[j]**(7/6)
            parameters.append([3*(gamma/theta)**2,theta,gamma])
        R,theta,gamma=np.array(parameters).T
        w.append(x);rs.append(R);ths.append(theta);gammas.append(gamma)
        # Two distinct frequency combinations, with the common prefactor removed.
        lmr=float(x@z)*float(x@(z/g.c.A));collective=float(x@(z*z/g.c.A))
        plasma.append([lmr,collective,lmr/collective])
    np.savez_compressed(OUT/'inputs.npz',cells=cells,weights=w,R=rs,theta=ths,gamma=gammas,
        candidate=stored['quantum'],frequency_combinations=plasma)
    paths=[g.ROOT/'verification/gr_quantum_mixture_certificate.py',g.ROOT/'verification/interval_records.py',
        g.ROOT/'verification/gr_dense_plasma.py',d.OUT/'manifest.json',d.OUT/'stellar-comparison.npz',
        d.plasma.split.reference.OUT/'reference-state.npz',d.OUT/'potekhin-chabrier2010.pdf',
        d.OUT/'baiko-chugunov2022.pdf']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='9493039',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        arithmetic_sources={str(p):g.c.sha(p) for p in g.Path(mp.__file__).parent.rglob('*.py')},
        precision_bits=128,fields=d.FIELDS,absolute_implementation_error_gate='1e-16',
        theta_uniform_max='0.1',R_domain='all R>0',
        exact_fit='BC2022 eq34 with exact decimal constants .351 and .294; positive R and theta. Differentiation uses R~n^(-1/3), theta~n^(1/2)/T at fixed species weights.',
        inputs='R,theta,positive weights are exact stored binary64 values. Their physical uncertainty and precomputation errors are excluded. The stored binary weights need not sum to exactly one; all mixture bounds explicitly sum those weights.',
        method='Positive sinh product gives exact bounds on the WK and degree-12 remainders and both logarithmic derivatives. Evaluate the degree-12 polynomial with outward interval arithmetic, then add its analytic remainder for all seven thermodynamic fields. Audit exported endpoints as exact rationals.',
        controls='Exact-rational polynomial replay and independent 60-digit free-energy derivatives at (R,theta)=(90,.002),(500,.1),(1900,.05),(120000,.01). Full 3206-cell interval evaluation without a trace cutoff.',
        physical_assumptions='The mixture WK coefficient is conditional on the bulk Boltzmann WK expansion, repulsive positive point ions, no exchange, and a rigid uniform neutralizing background. Finite-system surface terms, mobile-electron screening and partial ionization are excluded.',
        sources=['https://arxiv.org/abs/cond-mat/0701773','https://arxiv.org/abs/1001.0690','https://arxiv.org/abs/2112.04822'],
        physical_EOS_certified=False,native_EOS_replaced=False))
    print('PREPARED quantum mixture certificate',len(cells),flush=True)


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in plan['arithmetic_sources'].items():assert g.c.sha(g.Path(path))==digest,path
    return plan


def polynomial(R,theta,convert):
    """Reuse the prior derivative algebra; squared frequencies avoid sqrt."""
    c1=convert(Q(351,1000))*R/(90+R);s1=c1*c1;s2=convert(Q(294,1000))**2
    squares=[s1,s2,1-s1-s2];d1=90/(90+R);d1p=-R/(90+R)*d1
    D=[d1,convert(Q(0)),-s1/squares[2]*d1]
    DP=[d1p,convert(Q(0)),-s1/squares[2]*(2*(d1-D[2])*d1+d1p)]
    total=[convert(Q(0)) for _ in d.FIELDS]
    for square,dd,dp in zip(squares,D,DP):
        z=square*theta*theta;f=u=cv=convert(Q(0));power=convert(Q(1))
        for k,c in enumerate(COEFF,1):
            power*=z;term=convert(c)*power
            f+=term;u+=2*k*term;cv+=2*k*(1-2*k)*term
        a=convert(Q(1,2))-dd/3;p=a*u
        row=[f,u,p,u-f,cv,a*cv,(a+a*a+dp/9)*u-a*a*cv]
        total=[old+new for old,new in zip(total,row)]
    b2=theta**4*sum(s*s for s in squares)/2880
    b7=theta**14*sum(s**7 for s in squares)/1046139494400
    return total,b2,b7


def symbolic():
    q=sp.symbols('q',positive=True)
    L=lambda v:2*q*sp.diff(v,q)
    r=q-sp.log(1+q)
    assert sp.simplify(L(r)-2*q*q/(1+q))==0
    assert sp.simplify(L(L(r))-4*q*q*(2+q)/(1+q)**2)==0
    assert sp.simplify(L(L(r))-L(r)-2*q*q*(3+q)/(1+q)**2)==0
    assert sp.simplify(sp.diff(L(r)-r,q)-q*(3+q)/(1+q)**2)==0
    remainder=sp.log(1+q)-sum((-1)**(k+1)*q**k/k for k in range(1,7))
    assert sp.simplify(sp.diff(remainder,q)-q**6/(1+q))==0
    assert sp.simplify(L(remainder)-2*q**7/(1+q))==0
    assert sp.simplify(L(L(remainder))-L(remainder)-2*q**7*(13+11*q)/(1+q)**2)==0
    assert sp.simplify(sp.diff(L(remainder)-remainder,q)-q**6*(13+11*q)/(1+q)**2)==0
    # Positive coefficient differences prove both rational-function upper bounds.
    for numerator,limit in [(3+q,3),(13+11*q,13)]:
        assert all(c>=0 for c in sp.Poly(sp.expand(limit*(1+q)**2-numerator),q).all_coeffs())
    assert sp.simplify(sp.zeta(14)/(7*(2*sp.pi)**14)-sp.Rational(1,1046139494400))==0
    kap=Q(351,1000)**2/(1-Q(351,1000)**2-Q(294,1000)**2)
    assert kap<Q(1,6) and Q(1,6)*(2*Q(7,6)+Q(1,4))<Q(1,2)
    assert Q(1,2)+Q(1,6)/3<Q(2,3)
    n,T,K=sp.symbols('n T K',positive=True);f=K*n/T**2
    nder=lambda v:n*sp.diff(v,n);tder=lambda v:T*sp.diff(v,T)
    u=-tder(f);p=nder(f)
    seven=[f,u,p,u-f,u+tder(u),p+tder(p),p+nder(p)]
    assert [sp.simplify(v/f) for v in seven]==[1,2,1,1,-2,-1,2]
    # Arbitrary three-species identity demonstrates the general bilinear sums.
    ns=sp.symbols('n0:3',positive=True);zs=sp.symbols('Z0:3',positive=True);ms=sp.symbols('m0:3',positive=True)
    ni=sum(ns);ne=sum(a*b for a,b in zip(ns,zs));Qm=sum(a*b/c for a,b,c in zip(ns,zs,ms))
    lmr=sum(ns[j]/ni*ne*zs[j]/ms[j] for j in range(3))
    assert sp.simplify(ni*lmr-ne*Qm)==0
    energy=K*ne*Qm/T
    for j in range(3):
        assert sp.simplify(sp.diff(energy,ns[j])-K*zs[j]*(Qm+ne/ms[j])/T)==0
        for k in range(3):
            assert sp.simplify(sp.diff(energy,ns[j],ns[k])-K*zs[j]*zs[k]*(1/ms[j]+1/ms[k])/T)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        WK_bulk='Delta F^(2)/V = hbar^2*4*pi*e^2*n_e*sum(n_j Z_j/m_j)/(24*k_B*T). This is the coefficient of the assumed bulk WK expansion, not a bound on its physical remainder.',
        proof='Mass-weighted WK trace contains sum_a <Laplacian_a U>/m_a. Coulomb Poisson: Laplacian_a U = -4*pi*e^2*Z_a*sum_{b!=a} Z_b delta(r_a-r_b)+4*pi*e^2*Z_a*n_e. Repulsive classical contact weights vanish. The background term gives the displayed coefficient. Bulk/surface and WK existence assumptions are explicit.',
        chemical='For K=hbar^2*4*pi*e^2/(24*k_B), Qm=sum n_j Z_j/m_j: mu_j^(2)=K Z_j(Qm+n_e/m_j)/T, Hessian_jk=K Z_j Z_k(1/m_j+1/m_k)/T. Holding n_e fixed while changing all ion numbers is a different derivative.',
        positive_product='log(sinh(y/2)/(y/2))=sum_{l>=1} log(1+y^2/(4*pi^2*l^2)). Uniform convergence on compact positive y permits the two derivatives used here.',
        WK_remainder='0<=y^2/24-f(y)<=y^4/2880; |Delta(u,s,cv)| <= (4,3,12)*y^4/2880.',
        polynomial_remainder='For degree12 P6, log(1+q)-P6(q)=integral_0^q t^6/(1+t)dt. Value,u,s,cv bounds are (1,14,13,182)*y^14/1046139494400.',
        R_uniform='C1^2/C3^2<1/6, |D3|<1/6, |D1_prime|<=1/4, |D3_prime|<1/2. Therefore 0<a_i<=2/3 and |a_i+a_i^2+D_i_prime/9|<=7/6 for all R>0.',
        field_order=d.FIELDS,WK_multipliers=list(map(str,WK)),polynomial_multipliers=list(map(str,TAIL)),
        uniform_theta_point_one_WK_bounds=[str(v*Q(1,10)**4/2880) for v in WK],
        uniform_theta_point_one_polynomial_bounds=[str(v*Q(1,10)**14/1046139494400) for v in TAIL],
        scope='A conditional coefficient theorem and exact-real fit remainder theorem. No physical higher-order quantum-mixture or screened-EOS error bound.'))


def controls():
    plan=bindings();iv.prec=plan['precision_bits'];symbolic();rows=[]
    for R,theta in [(90.,.002),(500.,.1),(1900.,.05),(120000.,.01)]:
        enclosure,b2,b7=polynomial(iv.mpf(R),iv.mpf(theta),rat)
        exact,_,_=polynomial(Q(R),Q(theta),lambda v:v)
        for v,e in zip(enclosure,exact):assert lower(v)<=e<=upper(v)
        mp.mp.dps=60;r=mp.mpf(R);th=mp.mpf(theta)
        def free(s,t):
            rr=r*mp.exp(-s/3);tt=th*mp.exp(s/2-t)
            c1=mp.mpf('.351')*rr/(90+rr);c2=mp.mpf('.294');c3=mp.sqrt(1-c1*c1-c2*c2)
            return sum(mp.log(2*mp.sinh(c*tt/2)/(c*tt)) for c in [c1,c2,c3])
        f=free(0,0);u=-mp.diff(free,(0,0),(0,1));p=mp.diff(free,(0,0),(1,0))
        reference=[f,u,p,u-f,u-mp.diff(free,(0,0),(0,2)),p+mp.diff(free,(0,0),(1,1)),p+mp.diff(free,(0,0),(2,0))]
        for value,ref,mult in zip(enclosure,reference,TAIL):
            err=rat(mult)*b7;lo,hi=ends(interval_text(value+iv.mpf([-1,1])*err))
            assert lo<=Q(str(ref))<=hi
        rows.append(dict(R=R,theta=theta,exact_rational_polynomial_passed=True,independent_free_derivatives_contained=True))
    save('controls.json',dict(classification='Counterexample candidate',passed=True,rows=rows))


def run():
    plan=bindings();assert json.loads((OUT/'controls.json').read_text())['passed']
    iv.prec=plan['precision_bits'];a=dict(np.load(OUT/'inputs.npz'));records=[]
    maxima=[Q(0)]*7;wkmax=[Q(0)]*7;ratios=[]
    for k,cell in enumerate(a['cells']):
        result=[iv.mpf(0) for _ in d.FIELDS];b2=b7=leading=iv.mpf(0);candidate=np.zeros(7)
        for w,R,theta in zip(a['weights'][k],a['R'][k],a['theta'][k]):
            if w==0:continue
            assert w>0 and R>0 and 0<Q(float(theta))<=Q(plan['theta_uniform_max'])
            wi=iv.mpf(float(w));th=iv.mpf(float(theta))
            values,e2,e7=polynomial(iv.mpf(float(R)),th,rat)
            result=[s+wi*v for s,v in zip(result,values)];b2+=wi*e2;b7+=wi*e7;leading+=wi*th*th/24
            candidate+=w*d.quantum(float(R),float(theta))
        assert np.array_equal(candidate,a['candidate'][k]),int(cell)
        intervals=[];errors=[];wkerrors=[]
        for field,(value,mult,wmult) in enumerate(zip(result,TAIL,WK)):
            # Compute each final budget FROM the serialized outward endpoints.
            text=interval_text(value+iv.mpf([-1,1])*rat(mult)*b7);lo,hi=ends(text)
            err=max(abs(lo-Q(float(candidate[field]))),abs(hi-Q(float(candidate[field]))))
            werr=upper(rat(wmult)*b2)
            intervals.append(text);errors.append(str(err));wkerrors.append(str(werr))
            maxima[field]=max(maxima[field],err);wkmax[field]=max(wkmax[field],werr)
        lo,hi=ends(intervals[0]);leadlo,leadhi=ends(interval_text(leading))
        assert lo<=leadhi and leadlo-hi<=upper(b2)
        ratio=upper(b2)/lower(leading);ratios.append(ratio)
        records.append(dict(cell=int(cell),exact_fit_intervals=intervals,implementation_absolute_bounds=errors,
            WK_absolute_bounds=wkerrors,WK_free_interval=interval_text(leading),WK_relative_free_bound=str(ratio)))
        if len(records)%128==0:
            save('progress.json',dict(classification='Counterexample candidate',completed_cells=len(records)))
            print('QUANTUM INTERVAL',len(records),'/',len(a['cells']),flush=True)
    save('records.json',dict(classification='Proven',fields=d.FIELDS,rows=records))
    save('result.json',dict(classification='Proven',cells=len(records),fields_per_cell=7,
        passed=all(v<Q(plan['absolute_implementation_error_gate']) for v in maxima),
        maximum_implementation_bounds=list(map(str,maxima)),maximum_WK_bounds=list(map(str,wkmax)),
        maximum_relative_WK_free_bound=str(max(ratios)),
        floating_diagnostics=dict(classification='Counterexample candidate',
            implementation_max=max(map(float,maxima)),WK_max=list(map(float,wkmax)),WK_relative_free_max=float(max(ratios)),
            R_range=[float(a['R'].min()),float(a['R'].max())],theta_range=[float(a['theta'].min()),float(a['theta'].max())],
            naive_collective_frequency_ratio_range=[float(a['frequency_combinations'][:,2].min()),float(a['frequency_combinations'][:,2].max())]),
        fixed_input_only=True,physical_EOS_certified=False,native_EOS_replaced=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    plan=bindings();manifest=json.loads((OUT/'manifest.json').read_text())
    for rel,digest in manifest['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    a=dict(np.load(OUT/'inputs.npz'));r=json.loads((OUT/'result.json').read_text())
    rows=json.loads((OUT/'records.json').read_text())['rows']
    assert len(rows)==r['cells']==len(a['cells'])==3206
    assert [v['cell'] for v in rows]==list(a['cells'])
    maximum=[Q(0)]*7;wkmaximum=[Q(0)]*7
    for k,row in enumerate(rows):
        for j in range(7):
            lo,hi=ends(row['exact_fit_intervals'][j]);point=Q(float(a['candidate'][k,j]))
            budget=Q(row['implementation_absolute_bounds'][j]);assert lo<=hi
            assert max(abs(lo-point),abs(hi-point))==budget and budget<Q(plan['absolute_implementation_error_gate'])
            maximum[j]=max(maximum[j],budget);wkmaximum[j]=max(wkmaximum[j],Q(row['WK_absolute_bounds'][j]))
    assert list(map(str,maximum))==r['maximum_implementation_bounds']
    assert list(map(str,wkmaximum))==r['maximum_WK_bounds'] and r['passed']
    assert json.loads((OUT/'symbolic.json').read_text())['passed'] and json.loads((OUT/'controls.json').read_text())['passed']
    print('PASS 3206 x 7 fixed-input quantum thermodynamic enclosures and exact exported budgets; physical closure remains',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
