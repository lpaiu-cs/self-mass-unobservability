"""Finite-temperature static, relativistic ideal-electron medium response.

Proven: occupation transform, positive kernel and thermal-positron bound.
Counterexample candidate: finite stellar quadrature comparisons, not a full EOS.
"""
import json, shutil, sys
import numpy as np
import mpmath as mp
from mpmath import iv
import sympy as sp
from scipy.integrate import quad_vec
from scipy.special import expit
from interval_records import exact_endpoint, interval_text
import gr_screened_hamiltonian_runner as previous

g=previous.g;d=previous.original.d;OUT=g.OUT/'gr-finite-wavenumber-response'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=previous.scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();previous.verify()
    for stem in ['kozhberov-potekhin2021','decarvalho2017-response']:
        for ext in ['pdf','txt']:shutil.copy2(g.ROOT/'outputs'/f'{stem}.{ext}',OUT/f'{stem}.{ext}')
    paths=[g.ROOT/'verification/gr_finite_wavenumber_response.py',g.ROOT/'verification/interval_records.py',
        previous.OUT/'manifest.json',previous.OUT/'states.npz',
        d.plasma.OUT/'stellar-fermi-plasma.npz']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='5f4ccca',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        arithmetic_sources={str(p):g.c.sha(p) for p in g.Path(mp.__file__).parent.rglob('*.py')},
        sources=['https://arxiv.org/abs/2104.09964','https://arxiv.org/abs/1704.05944'],
        equations='Kozhberov/Potekhin2021 eq8 (second term coefficient 2/3); de Carvalho2017 eq22-23,76-81 for medium response linear in electron+positron occupations. Vacuum part is separate.',
        target='Derive and evaluate the finite-k finite-T ideal-electron RPA medium response; compare direct momentum-kernel quadrature to independent high-precision thermal convolution of the cold Jancovici function.',
        k_over_kappa=[0.,.001,.01,.03,.1,.3,1.,3.,10.,30.,100.,300.,1000.],
        quadrature_tolerances=[1e-8,2e-11],agreement_absolute_gate=1e-7,
        control_positions=[0,801,1603,2404,3205],control_columns=[0,1,4,6,8,12],
        mp_decimal_digits=60,control_absolute_gate=1e-8,
        kernel_controls_x=[.00001,.01,1.,10.],kernel_controls_y=[.00001,.1,.999999,1.,1.000001,10.,1000.],
        kernel_relative_gate=1e-11,series_terms=16,series_boundary=.25,
        interval_bits=128,thermal_pair_bound_gate='1e-100',
        input_scope='Frozen binary eta,beta,kappa and density scale, 3206 selected cells. k=0 normalization is the prior independently checked ideal-electron compressibility. No density/chemical-potential retuning.',
        physics_scope='Uniform noninteracting relativistic electron medium, static linear response, no external magnetic field. Include an analytic bound on omitted equilibrium thermal positrons at the same chemical potential. Vacuum polarization, local fields, exchange/correlations, nonlinear ion response and full EOS are not certified.',
        self_energy_scope='Finite-k response grid only; no finite grid treated as a continuum or ultraviolet-tail certificate. Polarization self-energy integral and thermodynamic derivatives require a subsequent calculation.',
        physical_EOS_certified=False,native_EOS_replaced=False))
    print('PREPARED finite-k finite-T response',flush=True)


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in p['arithmetic_sources'].items():assert g.c.sha(g.Path(path))==digest,path
    return p


def symbolic():
    x,Q,t,eta=sp.symbols('x Q t eta',positive=True)
    gamma=sp.sqrt(1+x*x);s=sp.sqrt(1+Q*Q)
    A=gamma*(1+x*x-3*Q*Q)/(6*Q);B=(2*Q*Q-1)*s/(6*Q)
    dL1=-2*Q/(x*x-Q*Q);dL2=-2*Q*s/(gamma*(x*x-Q*Q))
    assert sp.simplify(sp.diff(A,x)-x*(1+x*x-Q*Q)/(2*Q*gamma))==0
    assert sp.simplify(sp.diff(sp.log((Q*gamma+x*s)/(Q*gamma-x*s)),x)-dL2)==0
    nonlog=sp.diff(sp.Rational(2,3)*(x*gamma-Q*Q*sp.asinh(x)),x)+A*dL1+B*dL2
    assert sp.simplify(nonlog-x*x/gamma)==0
    assert sp.simplify(A+B-(x*x-Q*Q)*(gamma+(1-2*Q*Q)/(gamma+s))/(6*Q))==0
    f=1/(1+sp.exp(t-eta));assert sp.simplify(-sp.diff(f,t)-f*(1-f))==0
    u=sp.symbols('u',positive=True)
    h=sp.atanh(u)/u;j=1-(1-u*u)*h
    direct=Q*Q*u*u+u*(1+Q*Q*(u*u-1))*sp.atanh(u)
    assert sp.simplify(direct-u*u*(h+Q*Q*j))==0
    assert sp.simplify(sp.limit((x*x/gamma+x*(1+x*x-Q*Q)/(2*Q*gamma)*sp.log((x+Q)/(x-Q))),Q,0)-(2*x*x+1)/gamma)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        definitions='p and x are momenta/(m_e c); Q=hbar*k/(2*m_e*c); gamma(p)=sqrt(1+p^2); beta=kBT/(m_e*c^2); eta=(mu-m_e*c^2)/(kBT). K(x,Q)=x*gamma(x)*Psi_J(Q/x,x).',
        derivative='dK/dx=x^2/gamma+x*(1+x^2-Q^2)/(2*Q*gamma)*log|(x+Q)/(x-Q)|. The cusp at x=Q is integrable; K itself is continuous.',
        transform='f(t;eta)=integral_t^infinity f(s;eta)*(1-f(s;eta)) ds. Since the vacuum-subtracted one-loop medium response is linear in occupation, S(Q)=integral_0^infinity K(sqrt(beta*t*(2+beta*t)),Q)*f(t)*(1-f(t)) dt = integral_0^infinity f(p)*dK/dp dp. This transforms susceptibility, not epsilon inverse or a nonlinear RPA-resummed energy.',
        dielectric='epsilon_medium(k)=1+(4*alpha/pi)*(m_e*c/hbar)^2*S(Q)/k^2. For Q=0, S(0)=integral f(p)*(1+2*p^2)/sqrt(1+p^2) dp = integral p*gamma*f(t)*(1-f(t)) dt, reproducing kappa_0^2.',
        positivity='For u=p/Q<1 the momentum kernel times dp/du is Q*u^2/sqrt(1+Q^2*u^2)*(h+Q^2*j), h=atanh(u)/u>0, j=sum_{n>=1}2*u^(2n)/((2n-1)*(2n+1))>0. For u>1 both terms of Q/sqrt(1+Q^2*u^2)*(Q^2*u^2+u*(1+Q^2*(u^2-1))*atanh(1/u)) are positive. Hence S(Q)>0 at every finite Q and T>0.',
        series='For u<1/4, use h=sum_{n=0}^{16}u^(2n)/(2n+1), j=sum_{n=1}^{16}2*u^(2n)/((2n-1)*(2n+1)). Each omitted tail is positive; bounds are u^34/(35*(1-u^2)) and 2*u^34/(33*35*(1-u^2)), respectively. This is a local analytic truncation bound, not an interval certificate of the numerical quadrature.',
        thermal_pairs='For eta+1/beta>0 and t>=0, f_plus/f_minus <= exp(-eta-2/beta)*(1+exp(-eta)). Positivity of the common response kernel gives 0<=S_plus(Q)/S_minus(Q)<=the same bound for ALL Q, including Q=0. This compares fixed chemical potential occupations, not a density-retuned EOS or vacuum polarization.',
        scope='Exact conditional algebra of the declared ideal-electron medium model. Full interacting physical EOS remains open.'))


def cold(x,Q):
    """High-precision source eq8, combining the two cancelling cusp logs."""
    if x==0:return mp.mpf(0)
    gamma=mp.sqrt(1+x*x)
    if Q==0:return x*gamma
    s=mp.sqrt(1+Q*Q)
    cusp=0 if x==Q else (x*x-Q*Q)*(gamma+(1-2*Q*Q)/(gamma+s))*mp.log(abs((x+Q)/(x-Q)))/(6*Q)
    return 2*(x*gamma-Q*Q*mp.asinh(x))/3+cusp+(2*Q*Q-1)*s*mp.log((Q*gamma+x*s)/(Q+x))/(3*Q)


def positive_kernel(u,Q):
    gamma=np.sqrt(1+Q*Q*u*u)
    if u==0:return np.zeros_like(Q)
    if u<.25:
        z=u*u;power=1.;h=1.;j=0.
        for n in range(1,17):
            power*=z;h+=power/(2*n+1);j+=2*power/((2*n-1)*(2*n+1))
        bracket=u*u*(h+Q*Q*j)
    elif u<1:
        h=np.arctanh(u)/u;j=1-(1-u*u)*h;bracket=u*u*(h+Q*Q*j)
    else:
        assert u>1,'Integrable logarithmic cusp must be an integration boundary'
        bracket=Q*Q*u*u+u*(1+Q*Q*(u*u-1))*np.arctanh(1/u)
    return Q/gamma*bracket


def response(eta,beta,Q,S0,tolerance):
    if np.all(Q==0):
        # Independent momentum integral, rather than returning the normalization.
        scale=np.sqrt(beta)
        def f(u):
            p=scale*u;gamma=np.sqrt(1+p*p);t=p*p/(beta*(gamma+1))
            return scale*(1+2*p*p)/gamma*expit(eta-t)/S0
        return quad_vec(f,0,np.inf,epsabs=tolerance,epsrel=tolerance,norm='max',limit=1200)[0]
    def f(u):
        p=Q*u;gamma=np.sqrt(1+p*p);t=p*p/(beta*(gamma+1))
        return positive_kernel(u,Q)*expit(eta-t)/S0
    return sum(quad_vec(f,a,b,epsabs=tolerance/2,epsrel=tolerance/2,norm='max',limit=1600)[0]
        for a,b in [(0,1),(1,np.inf)])


def convolution(eta,beta,Q,S0):
    eta,beta,Q,S0=map(lambda x:mp.mpf(float(x)),[eta,beta,Q,S0])
    tq=Q*Q/(beta*(mp.sqrt(1+Q*Q)+1));top=max(mp.mpf(2),eta)
    points=sorted(set([mp.mpf(0),mp.mpf(1),top,top+8,top+40]+([tq] if 0<tq<top+80 else [])))+[mp.inf]
    def f(t):
        x=mp.sqrt(beta*t*(2+beta*t));w=1/(1+mp.exp(t-eta))/(1+mp.exp(eta-t))
        return cold(x,Q)*w/S0
    return float(mp.quad(f,points))


def controls():
    plan=bindings();mp.mp.dps=plan['mp_decimal_digits'];rows=[]
    for xv in plan['kernel_controls_x']:
        for yv in plan['kernel_controls_y']:
            x=mp.mpf(float(xv));Q=x*mp.mpf(float(yv));gamma=mp.sqrt(1+x*x)
            points=[mp.mpf(0)]+([Q] if Q<x else [])+[x]
            def f(p):
                if p==0:return mp.mpf(0)
                return p*p/mp.sqrt(1+p*p)+p*(1+p*p-Q*Q)/(2*Q*mp.sqrt(1+p*p))*mp.log(abs((p+Q)/(p-Q)))
            integral=mp.quad(f,points);source=cold(x,Q)
            score=float(abs(integral-source)/abs(source))
            rows.append(dict(x=xv,y=yv,source=float(source),integral=float(integral),score=score,passed=score<plan['kernel_relative_gate']))
    # Explicit tiny-u controls cover cancellation and the series/direct boundary.
    stable=[]
    for u in [1e-9,.001,.249999,.25,.999999,1.000001,10.,1e6]:
        for q in [.00001,1.,100.]:
            um,qm=mp.mpf(u),mp.mpf(q);p=um*qm
            exact=qm/mp.sqrt(1+p*p)*(p*p+um*(1+p*p-qm*qm)/2*mp.log(abs((um+1)/(um-1))))
            value=float(positive_kernel(u,np.array(q)));score=float(abs(mp.mpf(value)/exact-1))
            stable.append(dict(u=u,Q=q,score=score,passed=score<plan['kernel_relative_gate']))
    save('kernel-controls.json',dict(classification='Counterexample candidate',source_integrals=rows,stable_kernel=stable,passed=all(r['passed'] for r in rows+stable)))
    assert all(r['passed'] for r in rows+stable)


def run():
    plan=bindings();symbolic();controls();mp.mp.dps=plan['mp_decimal_digits'];iv.prec=plan['interval_bits']
    prior=dict(np.load(previous.OUT/'states.npz'));fermi=dict(np.load(d.plasma.OUT/'stellar-fermi-plasma.npz'))
    cells=prior['cells'];eta=fermi['eta'][cells];beta=fermi['beta'][cells]
    scale=fermi['dimensionless_density_scale'][cells];S0=scale*prior['normalized_integrals'][1]/beta
    c=d.plasma.constants();lam=c['hbar']/(c['m_e_g']*c['c'])
    ratios=np.array(plan['k_over_kappa']);Q=lam*prior['kappa_cm_inverse'][:,None]*ratios[None,:]/2
    arrays=[]
    for tolerance in plan['quadrature_tolerances']:
        columns=[]
        for j,r in enumerate(ratios):
            values=response(eta,beta,Q[:,j],S0,tolerance);columns.append(values)
            print('RESPONSE',tolerance,float(r),float(values.min()),float(values.max()),flush=True)
        arrays.append(np.array(columns).T)
    coarse,fine=arrays;checks=[]
    for i in plan['control_positions']:
        for j in plan['control_columns']:
            value=convolution(eta[i],beta[i],Q[i,j],S0[i]);delta=abs(value-fine[i,j])
            checks.append(dict(cell=int(cells[i]),column=j,k_over_kappa=float(ratios[j]),reference=value,direct=float(fine[i,j]),difference=float(delta),passed=delta<plan['control_absolute_gate']))
        print('CONVOLUTION CONTROL',int(cells[i]),flush=True)
    pairs=[];upper=[]
    for cell,e,b in zip(cells,eta,beta):
        ie,ib=iv.mpf(float(e)),iv.mpf(float(b));assert exact_endpoint((ie+1/ib)._mpi_[0])>0
        bound=iv.exp(-ie-2/ib)*(1+iv.exp(-ie));text=interval_text(bound)
        hi=exact_endpoint(bound._mpi_[1]);upper.append(hi)
        pairs.append(dict(cell=int(cell),ratio_bound=text))
    pair_max=max(upper);gate=previous.original.qc.Q(plan['thermal_pair_bound_gate'])
    pair_pass=pair_max<gate
    save('thermal-positron-bound.json',dict(classification='Proven',scope='Exact saved eta,beta; outward 128-bit arithmetic. Uniform all-Q relative response bound at fixed chemical potential.',records=pairs,passed=pair_pass,maximum_upper=str(pair_max),maximum_upper_float=float(pair_max)))
    epsilon=np.full_like(fine,np.inf);epsilon[:,1:]=1+fine[:,1:]/ratios[None,1:]**2
    np.savez_compressed(OUT/'states.npz',cells=cells,eta=eta,beta=beta,k_over_kappa=ratios,Q=Q,S0=S0,coarse=coarse,response_over_k0=fine,epsilon_medium=epsilon)
    agreement=float(np.max(abs(coarse-fine)));zero=float(np.max(abs(fine[:,0]-1)))
    norm_delta=float(np.max(abs((4*c['alpha']/np.pi)/lam**2*S0/prior['kappa_cm_inverse']**2-1)))
    passed=agreement<plan['agreement_absolute_gate'] and zero<plan['agreement_absolute_gate'] and norm_delta<plan['agreement_absolute_gate'] and pair_pass and all(v['passed'] for v in checks) and np.all(fine>0)
    save('result.json',dict(classification='Counterexample candidate',passed=passed,cells=len(cells),points=fine.size,
        finite_quadrature_agreement=agreement,k0_compressibility_difference=zero,normalization_identity_difference=norm_delta,
        independent_controls=checks,response_ranges=[dict(k_over_kappa=float(r),min=float(fine[:,j].min()),max=float(fine[:,j].max())) for j,r in enumerate(ratios)],
        maximum_positron_relative_bound=float(pair_max),physical_EOS_certified=False,native_EOS_replaced=False,
        scope='Finite all-selected-state RPA medium-response grid. No interacting dielectric, continuum quadrature enclosure, UV self-energy integral, EOS derivative, full GR evolution or observation closure is implied.'))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    p=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());a=dict(np.load(OUT/'states.npz'))
    assert r['passed'] and r['cells']==3206 and r['points']==3206*len(p['k_over_kappa'])
    assert np.all(a['response_over_k0']>0) and np.all(a['epsilon_medium'][:,1:]>1)
    assert np.all(np.isinf(a['epsilon_medium'][:,0]))
    for name in ['symbolic.json','kernel-controls.json','thermal-positron-bound.json']:assert json.loads((OUT/name).read_text())['passed']
    print('PASS finite-k finite-T ideal-electron response; whole physical EOS and continuum self-energy remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
