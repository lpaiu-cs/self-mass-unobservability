"""Thermodynamic derivatives of the finite-temperature RPA self/volume term."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from fractions import Fraction as F
import json, math, sys
import numpy as np
import sympy as sp
from mpmath import iv
from scipy.integrate import quad_vec
from scipy.special import expit
import gr_electron_density_certificate as density
import gr_polarization_self_scaled as previous
from interval_records import interval_text, exact_endpoint

g=density.g;OUT=g.OUT/'gr-polarization-thermodynamics'
FIELDS=['F','U','P','S','CV','PDT','PDR']
PARAMETERS=['value','eta','etaeta','tau','etatau','tautau']


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=previous.original.finite.previous.original.previous.scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();density.bindings();previous.verify()
    assert json.loads((density.OUT/'controls.json').read_text())['passed']
    paths=[g.ROOT/'verification/gr_polarization_thermodynamics.py',density.OUT/'plan.json',density.OUT/'candidate-manifest.json',
        density.OUT/'candidates.npz',density.OUT/'controls.json',previous.OUT/'manifest.json',density.ionic.OUT/'inputs.npz',density.ionic.OUT/'constants.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='2565748',
        bindings={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in paths},
        orders=[512,1024],block_size=32,processes=4,inner_tolerance=2e-11,z_cut=2**20,
        finite_refinement_gate=2e-7,control_positions=[0,1603,3205],control_steps=[.002,.001],
        finite_difference_gate=2e-5,control_inner_tolerance=2e-13,control_outer_tolerance=2e-12,
        tail_field_gate='1e-5',fields=FIELDS,parameter_order=PARAMETERS,
        target='Differentiate the actual finite-T finite-k ideal-electron RPA polarization self term at fixed ion composition, including the neutral chemical-potential response. Compute F,U,P,S,CV,PDT,PDR per ion in the usual kBT/kB units.',
        convention='Q=hbar*k/(2*m_e*c), B=alpha/pi, G=S/(Q^2/B+S), C_pol/N_i=-(2*alpha*<Z^2>*m_e*c^2/pi)*integral G dQ. The integration rescaling a and Sref are held fixed during differentiation at each reference state; physical alpha is fixed. Numerical B and prefactors are binary approximations, while tail bounds use outward arithmetic.',
        actual_density='Use the separate neutral-root predictors for each fully ionized saved baryon/isotope inventory. Final verification requires all3206 interval root/Hessian certificates and binds them after completion. Numerical interior is evaluated at predictor centers, not declared exact roots.',
        method='Differentiate the occupation analytically at fixed p,Q. Use the already checked positive momentum kernel and occupied-region rescaling. Integrate six parameter derivatives, then apply implicit neutral-density derivatives. Independent finite differences re-solve neutral eta at shifted ln(number),ln(T) and use separate adaptive outer value integration.',
        rigorous_tail='Positive response and thermal moment Cauchy-Schwarz bound all first/second parameter-derivative tails. Propagate them with certified neutral-root derivative intervals; this does not certify the interior numerical quadrature or its predictor-to-root displacement.',
        physical_EOS_certified=False,native_EOS_replaced=False))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    density.bindings();previous.bindings();return p


def inputs():
    a=dict(np.load(density.OUT/'candidates.npz'));ions=dict(np.load(density.ionic.OUT/'inputs.npz'))
    c=json.loads((density.ionic.OUT/'constants.json').read_text())['native_binary64']
    weights=ions['X']/ions['A'];weights/=weights.sum(1)[:,None];Z2=weights@(ions['Z']**2)
    Sref=np.sqrt(2*a['beta'])*a['momentum'][1];B=c['alpha']/np.pi;scale=np.sqrt(B*Sref)
    a.update(Sref=Sref,scale=scale,B=B,Z2=Z2,prefactor=-2*B*Z2*scale/a['beta'])
    return a


def response(eta,beta,Q,Sref,tolerance):
    def values(p2,base):
        gamma=np.sqrt(1+p2);t=p2/(beta*(gamma+1));q=expit(eta-t);v=expit(t-eta)
        return base*q*np.array([np.ones_like(q),v,v*(1-2*q),t*v,t*v*(1-2*q),(t*t*(1-2*q)-t)*v])
    if np.all(Q==0):
        def f(w):
            p2=beta*w*w
            return values(p2,np.sqrt(beta)*(1+2*p2)/np.sqrt(1+p2)/Sref)
        return quad_vec(f,0,np.inf,epsabs=tolerance,epsrel=tolerance,norm='max',limit=1600)[0]
    def f(u):return values(Q*Q*u*u,previous.BASE['positive_kernel'](u,Q)/Sref)
    t=np.maximum(eta,0)+1;s=float(np.max(np.sqrt(beta*t*(2+beta*t))/Q))
    options=dict(epsabs=tolerance/3,epsrel=tolerance/3,norm='max',limit=1600)
    if s<1/16:
        return quad_vec(lambda w:s*f(s*w),0,16,**options)[0]+quad_vec(f,16*s,1,**options)[0]+quad_vec(f,1,np.inf,**options)[0]
    return quad_vec(f,0,1,**options)[0]+quad_vec(f,1,np.inf,**options)[0]


def neutral_coefficients(m):
    N,E,EE,T,ET,TT=m;en=N/E;et=-T/E
    return np.array([en,et,(N-EE*en*en)/E,-en*(EE*et+ET)/E,-(TT+2*ET*et+EE*et*et)/E])


def physical_fields(raw,coeff,pref):
    J,Je,Jee,Jt,Jet,Jtt=raw;en,et,enn,ent,ett=coeff
    n=Je*en;t=Jt+Je*et;nn=Jee*en*en+Je*enn
    nt=Jee*en*et+Jet*en+Je*ent;tt=Jtt+2*Jet*et+Jee*et*et+Je*ett
    return pref*np.array([J,J-t,n,-t,t-tt,nt,n+nn])


def symbolic():
    e,b,c=sp.symbols('e b c',positive=True);t=c/b;q=1/(1+sp.exp(t-e));v=1-q
    expected=[q,q*v,q*v*(1-2*q),t*q*v,t*q*v*(1-2*q),(t*t*(1-2*q)-t)*q*v]
    actual=[q,sp.diff(q,e),sp.diff(q,e,2),b*sp.diff(q,b),b*sp.diff(sp.diff(q,e),b),b*sp.diff(b*sp.diff(q,b),b)]
    for left,right in zip(actual,expected):assert sp.simplify(left-right)==0
    x,S,Sa,Sb,Sab=sp.symbols('x S Sa Sb Sab',positive=True)
    assert sp.simplify(sp.diff(S/(x+S),S)-x/(x+S)**2)==0
    assert sp.simplify(sp.diff(S/(x+S),S,2)+2*x/(x+S)**3)==0
    n,T=sp.symbols('n T',positive=True);free=sp.sqrt(n/T)
    dn=lambda f:n*sp.diff(f,n);dt=lambda f:T*sp.diff(f,T)
    fields=[free,free-dt(free),dn(free),-dt(free),dt(free)-dt(dt(free)),dn(dt(free)),dn(free)+dn(dn(free))]
    factors=[1,sp.Rational(3,2),sp.Rational(1,2),sp.Rational(1,2),-sp.Rational(3,4),-sp.Rational(1,4),sp.Rational(3,4)]
    for f,k in zip(fields,factors):assert sp.simplify(f-k*free)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        occupation='At fixed p, tau=ln(beta), t=(gamma-1)/beta: f_eta=f(1-f); f_etaeta=f(1-f)(1-2f); f_tau=t*f(1-f); f_etatau=t*f(1-f)(1-2f); f_tautau=[t^2*(1-2f)-t]*f(1-f). The kernel has no eta/beta dependence at fixed p,Q.',
        nonlinear='For G=S/(X+S), X=Q^2/B: G_a=X*S_a/(X+S)^2; G_ab=X*S_ab/(X+S)^2-2*X*S_a*S_b/(X+S)^3. All six derivatives are of the same self free-energy function.',
        thermodynamics='For per-ion dimensional free energy f(n_e,T) at fixed full-ion composition, let nu=ln(n_e), tau=ln(T): U=f-f_tau, P/n_i=f_nu, S/n_i=-f_tau/T, Cv/n_i=(f_tau-f_tautau)/T, P_T/n_i=f_nutau/T, (n_e/n_i)*P_ne=f_nu+f_nunu. Rescaled integration coordinates are held fixed in these derivatives.',
        positive_control='For the classical TF self term f proportional to -n_e^(1/2)*T^(-1/2), the seven relative multipliers are [1,3/2,1/2,1/2,-3/4,-1/4,3/4]. Holding eta fixed would give a different temperature derivative and fails this thermodynamic construction.',
        tail='Define M_j(Q)=integral K_p*f*t^j dp and A_j=(integral M_j dQ)/(a*Sref), c=a^2/(B*Sref). Beyond z=Z the absolute tails of [G,G_eta,G_etaeta,G_tau,G_etatau,G_tautau] are bounded by [A0,A0,3*A0,A1,3*A1,A1+3*A2]/(c*Z^2). Use |S_eta|,|S_etaeta|<=S; |S_tau|,|S_etatau|<=M1; |S_tautau|<=M1+M2; and M1^2<=S*M2.',
        moment='integral M_j(Q)dQ=(pi^2*beta/4)*integral t^j*(1+beta*t)^2*f(t)dt. Polynomial incomplete exponential moments bound this at any compact eta,beta box; together with the finite-Q positive smooth response they justify first/second differentiation under the self integral.',
        scope='Conditional exact derivative identities and tail bounds. Numerical interior, predictor displacement, full correlated/partially-ionized EOS and actual GR/observational closure remain separate.'))


def block(order,first,count):
    plan=json.loads((OUT/'plan.json').read_text());a=inputs();nodes,weights=np.polynomial.legendre.leggauss(order)
    end=np.arctan(float(plan['z_cut']));angles=(nodes+1)*end/2;weights*=end/2;total=np.zeros((6,len(a['cells'])))
    for j in range(first,first+count):
        z=np.tan(angles[j]);Q=a['scale']*z;R=response(a['eta_center'],a['beta'],Q,a['Sref'],plan['inner_tolerance'])
        assert np.all(R[0]>0);X=Q*Q/(a['B']*a['Sref']);den=X+R[0];k=X/den**2;h=2*X/den**3
        G=np.array([R[0]/den,k*R[1],k*R[2]-h*R[1]**2,k*R[3],k*R[4]-h*R[1]*R[3],k*R[5]-h*R[3]**2])
        total+=weights[j]*(1+z*z)*G
    target=OUT/f'block-{order}-{first:04d}.npz';assert not target.exists()
    np.savez_compressed(target,order=order,first=first,count=count,weighted_sum=total,angles=angles[first:first+count],weights=weights[first:first+count])
    return target.name


def run():
    plan=bindings();symbolic();labels=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        tasks=[pool.submit(block,n,j,min(plan['block_size'],n-j)) for n in plan['orders'] for j in range(0,n,plan['block_size'])]
        for done in as_completed(tasks):
            labels.append(done.result());save('progress.json',dict(completed_blocks=sorted(labels)));print('SELF THERMO BLOCKS',len(labels),flush=True)
    a=inputs();coef=neutral_coefficients(a['momentum']);orders=[];fields=[]
    for n in plan['orders']:
        raw=sum(np.load(OUT/f'block-{n}-{j:04d}.npz')['weighted_sum'] for j in range(0,n,plan['block_size']))
        orders.append(raw);fields.append(physical_fields(raw,coef,a['prefactor']))
    score=float(np.max(abs(fields[1]-fields[0])))
    np.savez_compressed(OUT/'states.npz',cells=a['cells'],eta=a['eta_center'],beta=a['beta'],Sref=a['Sref'],scale=a['scale'],prefactor=a['prefactor'],Z2=a['Z2'],neutral_coefficients=coef,
        coarse_parameters=orders[0],fine_parameters=orders[1],coarse_fields=fields[0],fine_fields=fields[1])
    save('numerical.json',dict(classification='Counterexample candidate',passed=score<plan['finite_refinement_gate'],maximum_field_refinement_difference=score,
        field_ranges={name:[float(v.min()),float(v.max())] for name,v in zip(FIELDS,fields[1])},density_certificate_attached=False))
    print('SELF THERMO REFINEMENT',score,flush=True)


def controls():
    plan=bindings();a=inputs();ix=np.array(plan['control_positions']);eta0=a['eta_center'][ix];beta0=a['beta'][ix];Sref=a['Sref'][ix];scale=a['scale'][ix]
    target=a['target'][ix];pref=a['prefactor'][ix];end=np.arctan(float(plan['z_cut']));base=dict(np.load(OUT/'states.npz'))
    def free(nu,tau):
        beta=beta0*np.exp(tau);wanted=target*np.exp(nu-1.5*tau);eta=eta0.copy()
        for _ in range(4):
            m=density.momentum(eta,beta,2e-13);eta-=(m[0]-wanted)/m[1]
        def f(angle):
            z=np.tan(angle);Q=scale*z;R=previous.response(eta,beta,Q,Sref,plan['control_inner_tolerance']);X=Q*Q/(a['B']*Sref)
            return pref*(1+z*z)*R/(X+R)
        return quad_vec(f,0,end,epsabs=plan['control_outer_tolerance'],epsrel=plan['control_outer_tolerance'],norm='max',limit=600)[0]
    zero=free(0,0);rows=[]
    for h in plan['control_steps']:
        p,m,t,b=free(h,0),free(-h,0),free(0,h),free(0,-h)
        n=(p-m)/(2*h);dt=(t-b)/(2*h);nn=(p-2*zero+m)/h**2;tt=(t-2*zero+b)/h**2
        nt=(free(h,h)-free(h,-h)-free(-h,h)+free(-h,-h))/(4*h*h)
        observed=np.array([zero,zero-dt,n,-dt,dt-tt,nt,n+nn]);expected=base['fine_fields'][:,ix]
        differences=abs(observed-expected);score=float(differences.max())
        rows.append(dict(step=h,score=score,field_maxima=dict(zip(FIELDS,map(float,differences.max(1)))),passed=score<plan['finite_difference_gate']))
        print('SELF THERMO DIFFERENCES',h,score,flush=True)
    save('controls.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),controls=rows,positions=plan['control_positions']))


def certify():
    plan=bindings();density.verify();iv.prec=128;a=inputs();root=json.loads((density.OUT/'result.json').read_text())
    roots={r['position']:r for r in root['rows']};ions=dict(np.load(density.ionic.OUT/'inputs.npz'));c=json.loads((density.ionic.OUT/'constants.json').read_text())['native_binary64']
    I=density.I;upper=density.high;records=[];maxima=[F(0)]*7
    for i in range(len(a['cells'])):
        r=roots[i];eta=density.read_interval(r['root'][1:-1].split(','));beta=I(a['beta'][i]);scale=I(a['scale'][i]);Sref=I(a['Sref'][i]);B=I(c['alpha'])/iv.pi
        coeff=[density.read_interval(r['derivatives'][name][1:-1].split(',')) for name in density.FIELDS];en,et,enn,ent,ett=coeff
        ec=eta.b;cut=I(max(F(0),upper(eta)));moments=[]
        def bound(k):return cut**(k+1)/(k+1)+iv.exp(ec-cut)*sum(I(math.factorial(k))/math.factorial(j)*cut**j for j in range(k+1))
        for j in range(3):moments.append(iv.pi**2*beta/4*(bound(j)+2*beta*bound(j+1)+beta**2*bound(j+2))/(scale*Sref))
        m0,m1,m2=moments;factor=1/((scale**2/(B*Sref))*plan['z_cut']**2)
        q0,qe,qee,qt,qet,qtt=[x*factor for x in [m0,m0,3*m0,m1,3*m1,m1+3*m2]]
        n=abs(en)*qe;t=qt+abs(et)*qe;nn=en**2*qee+abs(enn)*qe
        nt=abs(en*et)*qee+abs(en)*qet+abs(ent)*qe;tt=qtt+2*abs(et)*qet+et**2*qee+abs(ett)*qe
        counts=[I(x)/I(m) for x,m in zip(ions['X'][i],ions['A'])];Z2=sum(x*I(z)**2 for x,z in zip(counts,ions['Z']))/sum(counts)
        pref=2*B*Z2*scale/beta;bounds=[pref*x for x in [q0,q0+t,n,t,t+tt,nt,n+nn]]
        maxima=[max(old,upper(v)) for old,v in zip(maxima,bounds)]
        records.append(dict(cell=int(a['cells'][i]),field_tail_upper={name:interval_text(v) for name,v in zip(FIELDS,bounds)}))
    save('tail-certificate.json',dict(classification='Proven',passed=max(maxima)<F(plan['tail_field_gate']),records=records,
        maxima_exact=dict(zip(FIELDS,map(str,maxima))),display_only=dict(zip(FIELDS,map(float,maxima))),
        scope='Infinite-wavenumber tails at certified neutral roots, including first/second neutral derivatives and exact declared ion inventories. No certified numerical interior or predictor-to-root error is included.'))
    save('density-binding.json',dict(sha256={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in [density.OUT/'manifest.json',density.OUT/'result.json']}))
    save('manifest.json',dict(sha256={x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in OUT.iterdir() if x.is_file()}))
    verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'density-binding.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for name in ['symbolic.json','numerical.json','controls.json','tail-certificate.json']:assert json.loads((OUT/name).read_text())['passed'],name
    print('PASS neutral-density self thermodynamics finite checks and certified derivative tails; interior/full EOS remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
