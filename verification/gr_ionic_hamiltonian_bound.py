"""Nonperturbative free-energy brackets for positive ionic jellium.

Proven: conditional on the declared Boltzmann Coulomb Hamiltonian.
Conjectural: electron response, exchange, ionization and derivative closure.
"""
import json, shutil, sys
from fractions import Fraction as Q
import numpy as np
import sympy as sp
from mpmath import iv
import gr_quantum_mixture_certificate as qc

g=qc.g;OUT=g.OUT/'gr-ionic-hamiltonian-bound'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();qc.verify()
    for stem in ['kleinert-smearing1998','golden-thompson2014']:
        for ext in ['pdf','txt']:shutil.copy2(g.ROOT/'outputs'/f'{stem}.{ext}',OUT/f'{stem}.{ext}')
    a=dict(np.load(qc.OUT/'inputs.npz'));state,_=qc.d.state_data();cells=a['cells']
    np.savez_compressed(OUT/'inputs.npz',cells=cells,rho=np.exp(state['lnd'][cells]),
        temperature=np.exp(state['lnT'][cells]),X=state['X'][cells],A=g.c.A,Z=g.c.Z)
    save('constants.json',dict(classification='Imported from prior work',
        native_binary64=qc.d.plasma.constants(),provider_AUM=1822.88848,
        masses='Compare declared nominal A/N_A grams and provider A*1822.88848*m_e grams. Neither is a measured isotope mass certificate. Stored native constants are exact parameters of these two Hamiltonians.'))
    paths=[g.ROOT/'verification/gr_ionic_hamiltonian_bound.py',qc.OUT/'manifest.json',
        qc.OUT/'records.json',qc.d.plasma.split.reference.OUT/'reference-state.npz']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='29cf073',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},precision_bits=128,
        input_semantics='The saved rho,T,X and native constants are exact binary64 parameters. Number densities are recomputed with outward arithmetic: n_j=N_A*rho*X_j/A_j, n_e=sum Z_j*n_j. No rounded sum is assumed exactly neutral.',
        Hamiltonian='H=sum p_a^2/(2*m_a)+sum_{a<b}q_a*q_b*G_L(r_a-r_b)+C, q_a>0. G_L is the zero-mean periodic 3D Coulomb Green function with Laplacian G_L=-4*pi*(delta_L-1/V). C is any position-independent self/background convention common to classical and quantum systems. Boltzmann labelled trace divided by identical species factorials; no exchange.',
        proof='Retain finite-box winding factors in Golden-Thompson. Retain only positive zero-winding paths for the Jensen lower partition bound. Centroid covariance and positive periodic heat kernel give the global smearing bound. Take bulk limit only when the two free-energy densities exist.',
        controls='Symbolic Brownian bridge and Fourier centroid variance, heat-smearing derivative, three-particle pair-sum reduction, and a value-to-derivative counterexample. Exact rational audit of every exported physical-error bracket.',
        reduced_upper_diagnostic_ceiling='0.001',
        excluded='Mobile/polarizable electrons; Fermi/Bose exchange; negative point charges; finite nuclear radii; reactions/partial ionization; thermodynamic derivatives; uncertainty in the classical correlation free energy and in constants/state parameters; curved or evolving stellar cells.',
        sources=['https://journals.aps.org/pr/abstract/10.1103/PhysRev.137.B1127',
            'https://arxiv.org/abs/1408.2008','https://link.aps.org/doi/10.1103/PhysRevA.34.5080','https://arxiv.org/abs/quant-ph/9806016'],
        whole_physical_EOS_certified=False,native_EOS_replaced=False))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text());qc.bindings()
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return plan


def symbolic():
    s,t=sp.symbols('s t',nonnegative=True)
    covmean=sp.integrate(t*(1-s),(t,0,s))+sp.integrate(s*(1-t),(t,s,1))
    varmean=sp.integrate(covmean,(s,0,1))
    assert varmean==sp.Rational(1,12)
    assert sp.simplify(s*(1-s)-2*covmean+varmean)==sp.Rational(1,12)
    assert sp.simplify(2*sp.zeta(2)/(2*sp.pi)**2)==sp.Rational(1,12)
    r,t=sp.symbols('r t',positive=True)
    smeared=sp.erf(r/(2*sp.sqrt(t)))/r
    kernel=sp.exp(-r*r/(4*t))/(4*sp.pi*t)**sp.Rational(3,2)
    assert sp.simplify(sp.diff(smeared,t)+4*sp.pi*kernel)==0
    assert sp.simplify(sp.diff(smeared,r,2)+2*sp.diff(smeared,r)/r+4*sp.pi*kernel)==0
    charges=sp.symbols('q0:3',positive=True);v=sp.symbols('v0:3',positive=True)
    pair=sum(charges[a]*charges[b]*(v[a]+v[b]) for a in range(3) for b in range(a+1,3))
    one=sum(v[a]*charges[a]*(sum(charges)-charges[a]) for a in range(3))
    assert sp.expand(pair-one)==0
    h,beta,V=sp.symbols('h beta V',positive=True);m=sp.symbols('m0:3',positive=True)
    exact=2*sp.pi*pair.subs({v[a]:beta*h*h/(12*m[a]) for a in range(3)})/V
    expected=beta*h*h*4*sp.pi*sum(charges[a]*(sum(charges)-charges[a])/m[a] for a in range(3))/(24*V)
    assert sp.simplify(exact-expected)==0
    B,k,x=sp.symbols('B k x',positive=True);example=B*(1+sp.sin(k*x))/2
    assert sp.diff(example,x).subs(x,0)==B*k/2
    save('symbolic.json',dict(classification='Proven',passed=True,
        finite_bracket='-k_B*T*sum_a log W_a <= F_Q-F_cl <= C_N, W_a=sum_{nu in Z^3}exp[-m_a L^2 |nu|^2/(2 beta hbar^2)], C_N=beta*hbar^2*4*pi/(24 V)*sum_a q_a*(sum_b q_b-q_a)/m_a.',
        bulk_bracket='If the bulk free-energy densities exist at fixed positive T and composition: 0 <= f_Q-f_cl <= hbar^2*4*pi*e^2*n_e*sum_j(n_j Z_j/m_j)/(24*k_B*T). No small-theta or WK-convergence assumption is used.',
        covariance='Free zero-winding loops with their centroid removed have one-coordinate variance v_a=beta*hbar^2/(12*m_a); Brownian bridge subtraction and Fourier sums agree.',
        smearing='For t_ab=(v_a+v_b)/2, exp(t_ab Laplacian)G_L-G_L=4*pi*t_ab/V-4*pi*integral_0^t_ab K_L(s,r)ds <=4*pi*t_ab/V. K_L is the positive periodic heat kernel.',
        Jensen='At fixed centroids, E exp[-beta integral U] >= exp[-beta E integral U] >= exp[-beta(U+C_N)]. Integrating centroids and retaining only the positive zero-winding sector gives Z_Q>=Z_cl exp(-beta C_N).',
        trace='Golden-Thompson gives Z_Q<=Z_cl product W_a using the periodic free kinetic diagonal. For positive Coulomb singularities the finite-box self-adjoint form/Feynman-Kac trace is used; bounded-potential regularization followed by its semigroup limit justifies the trace inequality.',
        winding='Writing a=m L^2/(2 beta hbar^2)>0, log W<=6 exp(-a)/(1-exp(-a)), since n^2>=n and log(1+x)<=x. Hence sum log W/V vanishes in the fixed-density bulk limit.',
        derivative_boundary='A bracket on free energy cannot be differentiated as an inequality. Functions B*(1+sin(k*x))/2 remain between 0 and B while their derivative at zero is B*k/2. A physical derivative certificate requires additional control.',
        scope='A nonperturbative theorem for the declared positive-ion rigid-background Hamiltonian. The total stellar EOS and its derivatives remain open.'))


def run():
    plan=bindings();symbolic();iv.prec=plan['precision_bits']
    a=dict(np.load(OUT/'inputs.npz'));constants=json.loads((OUT/'constants.json').read_text())
    c={k:iv.mpf(v) for k,v in constants['native_binary64'].items()}
    oldrows=json.loads((qc.OUT/'records.json').read_text())['rows'];rows=[];maxima=[Q(0),Q(0)]
    massratio=1/(c['N_A']*c['m_e_g']*iv.mpf(constants['provider_AUM']))
    ceilings=[]
    for k,cell in enumerate(a['cells']):
        rho=iv.mpf(float(a['rho'][k]));T=iv.mpf(float(a['temperature'][k]));numbers=[]
        for X,A in zip(a['X'][k],a['A']):numbers.append(c['N_A']*rho*iv.mpf(float(X))/int(A))
        ni=sum(numbers);ne=sum(n*int(z) for n,z in zip(numbers,a['Z']))
        sums=[iv.mpf(0),iv.mpf(0)]
        for n,A,Z in zip(numbers,a['A'],a['Z']):
            masses=[int(A)/c['N_A'],int(A)*iv.mpf(constants['provider_AUM'])*c['m_e_g']]
            sums=[s+n*int(Z)/mass for s,mass in zip(sums,masses)]
        factor=c['hbar']**2*4*iv.pi*c['e_esu']**2*ne/(24*(c['k_B']*T)**2*ni)
        limits=[factor*s for s in sums]
        assert oldrows[k]['cell']==int(cell)
        fitlo,fithi=qc.ends(oldrows[k]['exact_fit_intervals'][0]);branches=[]
        for j,limit in enumerate(limits):
            text=qc.interval_text(limit);lo,hi=qc.ends(text);assert 0<lo<=hi
            err=max(abs(fitlo),abs(fithi),abs(hi-fitlo),abs(hi-fithi))
            maxima[j]=max(maxima[j],hi)
            branches.append(dict(upper_coefficient_interval=text,physical_quantum_free_bracket=['0',str(hi)],
                physical_minus_declared_fit=[str(-fithi),str(hi-fitlo)],absolute_fit_error_bound=str(err)))
        ceilings.append(max(qc.ends(qc.interval_text(v))[1] for v in limits)<Q(plan['reduced_upper_diagnostic_ceiling']))
        rows.append(dict(cell=int(cell),nominal_mass=branches[0],provider_mass=branches[1]))
    save('records.json',dict(classification='Proven',normalization='F/(N_i k_B T), bulk homogeneous rigid-background model at each saved state.',rows=rows))
    save('result.json',dict(classification='Proven',cells=len(rows),nonperturbative=True,
        maximum_upper_coefficients=list(map(str,maxima)),mass_normalization_ratio=qc.interval_text(massratio),
        numerical_diagnostic=dict(classification='Counterexample candidate',all_below_declared_ceiling=all(ceilings),
            ceiling=plan['reduced_upper_diagnostic_ceiling'],upper_coefficients=list(map(float,maxima))),
        free_energy_only=True,thermodynamic_derivatives_certified=False,whole_physical_EOS_certified=False,native_EOS_replaced=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    plan=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    a=dict(np.load(OUT/'inputs.npz'));rows=json.loads((OUT/'records.json').read_text())['rows']
    fits=json.loads((qc.OUT/'records.json').read_text())['rows'];maxima=[Q(0),Q(0)]
    assert [r['cell'] for r in rows]==list(a['cells']) and len(rows)==3206
    for row,fit in zip(rows,fits):
        fitlo,fithi=qc.ends(fit['exact_fit_intervals'][0]);assert row['cell']==fit['cell']
        for j,name in enumerate(['nominal_mass','provider_mass']):
            b=row[name];lo,hi=qc.ends(b['upper_coefficient_interval']);assert 0<lo<=hi
            assert b['physical_quantum_free_bracket']==['0',str(hi)]
            assert b['physical_minus_declared_fit']==[str(-fithi),str(hi-fitlo)]
            assert Q(b['absolute_fit_error_bound'])==max(abs(fitlo),abs(fithi),abs(hi-fitlo),abs(hi-fithi))
            maxima[j]=max(maxima[j],hi)
    result=json.loads((OUT/'result.json').read_text());assert list(map(str,maxima))==result['maximum_upper_coefficients']
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS positive-ion Hamiltonian bulk free-energy brackets at 3206 states; screened EOS and derivatives remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
