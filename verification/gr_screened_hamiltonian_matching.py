"""Match the full static-screening Hamiltonian, including its volume term.

Proven: source algebra and exact non-identifiability from kappa alone.
Counterexample candidate: ideal-electron screening diagnostics at saved states.
"""
import json, shutil, sys
import numpy as np
import mpmath as mp
import sympy as sp
from scipy.integrate import quad, quad_vec
from scipy.special import expit
import gr_ionic_hamiltonian_bound as ionic

g=ionic.g;qc=ionic.qc;d=qc.d;OUT=g.OUT/'gr-screened-hamiltonian-matching'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();ionic.verify()
    for stem in ['chabrier-potekhin1998','potekhin-chabrier-screening2013']:
        for ext in ['pdf','txt']:shutil.copy2(g.ROOT/'outputs'/f'{stem}.{ext}',OUT/f'{stem}.{ext}')
    plasma=d.plasma.OUT/'stellar-fermi-plasma.npz'
    paths=[g.ROOT/'verification/gr_screened_hamiltonian_matching.py',ionic.OUT/'manifest.json',
        ionic.OUT/'inputs.npz',qc.OUT/'inputs.npz',plasma,d.plasma.OUT/'manifest.json']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='154c99d',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        sources=['https://www.ioffe.ru/astro/Stars/Paper/cp98e.pdf','https://arxiv.org/abs/1310.3162'],
        source_equations='CP1998 equations 7-10; PC2013 section 4.1. Preserve the unscreened self subtraction and do not identify a long-wavelength kappa with the full epsilon(k).',
        target='Derive pair/volume decomposition for all isotope charges, state-dependent Hamiltonian thermodynamics, and a positive static dielectric family with identical kappa but different polarization self energies and real-space sign.',
        dielectric_family='1/epsilon_a(k)=x^2/(1+x^2)+a*x^4/(1+x^2)^3, x=k/kappa, 0<=a<=1. This is an analytic insufficiency witness, not a microscopic electron dielectric prediction.',
        electrons='Use frozen native eta,beta and ideal relativistic electron-only occupation; this is the k=0 response, not finite-k Lindhard/QED, local-field corrections, electron correlations or a total stellar EOS.',
        quadrature_tolerances=[1e-9,2e-12],finite_agreement_gate=1e-8,
        control_positions=[0,801,1603,2404,3205],mp_decimal_digits=50,
        kernel_controls=[.1,1.,3.,5.,10.],kernel_absolute_gate=1e-10,
        self_energy_scope='The two analytic dielectric witnesses are evaluated at the actual ideal-electron kappa. Their difference is a missing-information diagnostic, not a physical error bar or a new correction to add to native/PC EOS.',
        physical_EOS_certified=False,native_EOS_replaced=False))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return plan


def symbolic():
    x,a=sp.symbols('x a',nonnegative=True);inv=x*x/(1+x*x)+a*x**4/(1+x*x)**3
    assert sp.simplify(1-inv-((1-a)*x**4+2*x*x+1)/(1+x*x)**3)==0
    assert sp.limit(inv/x**2,x,0)==1 and sp.limit(inv,x,sp.oo)==1
    integral=sp.integrate(1/(1+x*x),(x,0,sp.oo))
    correction=sp.integrate(x**4/(1+x*x)**3,(x,0,sp.oo))
    assert integral==sp.pi/2 and correction==3*sp.pi/16
    k,r=sp.symbols('k r',positive=True)
    yuk=sp.exp(-k*r)/r
    squared=-sp.diff(yuk,k)/(2*k)
    cubed=-sp.diff(squared,k)/(4*k)
    real=sp.simplify(k*k*(squared-k*k*cubed))
    assert sp.simplify(real-k*sp.exp(-k*r)*(3-k*r)/8)==0
    total=yuk+a*real
    assert sp.limit(total-1/r,r,0)==-k+3*a*k/8
    assert sp.simplify(total.subs({k:1,r:5,a:1})+sp.exp(-5)/20)==0
    rho2,Q2,eps,v=sp.symbols('rho2 Q2 eps v',positive=True)
    full=v*(rho2/eps-Q2)
    pair=v/eps*(rho2-Q2);volume=Q2*(v/eps-v)
    assert sp.simplify(full-pair-volume)==0
    n,T,C=sp.symbols('n T C',positive=True);kap=sp.Function('kap')(n,T);alpha=sp.Function('a')(n,T)
    free=n*C*kap*(-sp.Rational(1,2)+3*alpha/16)
    pressure=sp.expand(n*sp.diff(free,n)-free);energy=sp.expand(free-T*sp.diff(free,T))
    assert sp.simplify(pressure-n*n*C*((-sp.Rational(1,2)+3*alpha/16)*sp.diff(kap,n)+3*kap*sp.diff(alpha,n)/16))==0
    assert sp.simplify(energy-n*C*((-sp.Rational(1,2)+3*alpha/16)*(kap-T*sp.diff(kap,T))-3*kap*T*sp.diff(alpha,T)/16))==0
    E=sp.Function('E')(T);H=sp.Function('H')(T)
    assert sp.diff(H-T*sp.diff(H,T),T)==-T*sp.diff(H,T,2)
    save('symbolic.json',dict(classification='Proven',passed=True,
        mixture='rho_Z(k)=sum_a Z_a exp(i k.r_a); Q2=sum_a Z_a^2. H_eff=K+(1/2V)sum_{k!=0}v_c(k)(|rho_Z|^2/epsilon-Q2). Pair part uses v_c/epsilon times (|rho_Z|^2-Q2); the remaining volume term is Q2/(2V) sum(v_c/epsilon-v_c).',
        self='C_pol/N_i=<Z^2>*e^2/pi*integral_0^infinity [1/epsilon(k)-1] dk. The Yukawa value is -<Z^2>*e^2*kappa/2. State-dependent volume terms cancel only in a fixed-state quantum-minus-classical comparison; they do not disappear from the EOS.',
        family='For 0<=a<=1 the declared inverse dielectric lies between zero and one, has the same leading k^2/kappa^2 at k=0, and tends to one at infinity. Its self energy is <Z^2>*e^2*kappa*(-1/2+3a/16). Changing a from zero to one changes its magnitude by 37.5% at the same kappa.',
        real_space='v_a(r)=e^2 exp(-kappa*r)[1/r+a*kappa*(3-kappa*r)/8]. At a=1,kappa*r=5 it is -e^2*kappa*exp(-5)/20. A positive Fourier dielectric and a known kappa do not prove a positive real-space pair potential or the previous Coulomb Poisson bound.',
        thermodynamics='For a classical canonical effective Hamiltonian H(T,V), U=<H-T H_T>, P=-<H_V>, C_V=-T<H_TT>+Var(H-T H_T)/(k_B T^2). Configuration and momentum coordinates, measure and external parameters must be specified when taking H_V. For quantum noncommuting operators the variance is replaced by the Kubo covariance.',
        self_derivatives='With fixed composition and f_self=n*C*kappa*A(a), A=-1/2+3a/16: P_self=n^2*C*(A*kappa_n+3*kappa*a_n/16), U_self=n*C*(A*(kappa-T*kappa_T)-3*kappa*T*a_T/16). Knowledge of kappa alone does not fix these derivatives.',
        compressibility='For ideal electrons without positrons and a T-independent dispersion, chi=(kBT/n_e)*dn_e/dmu=integral f(1-f)/integral f, hence 0<chi<=1; kappa_0^2=4*pi*e^2*n_e*chi/(kBT). This is only a k=0 result.',
        scope='Exact conditional algebra and a static-response insufficiency witness. No microscopic dielectric or full physical EOS certificate.'))


def integrals(eta,beta,scale,tolerance):
    def f(t):
        gamma=1+beta*t;mom=np.sqrt(beta*t*(2+beta*t));occ=expit(eta-t)
        density=beta*gamma*mom*occ/scale
        return np.array([density,density*expit(t-eta)])
    return quad_vec(f,0,np.inf,epsabs=tolerance,epsrel=tolerance,limit=600)[0]


def independent(eta,beta,scale):
    eta,beta,scale=map(lambda v:mp.mpf(float(v)),[eta,beta,scale])
    def f(p,j):
        gamma=mp.sqrt(1+p*p);t=p*p/((gamma+1)*beta);occ=1/(1+mp.exp(t-eta))
        return p*p*occ*(1-occ if j else 1)/scale
    points=[mp.mpf(0)]+[mp.sqrt(beta*t*(2+beta*t)) for t in [1,max(2,float(eta)),max(2,float(eta))+8,max(2,float(eta))+40]]+[mp.inf]
    return np.array([float(mp.quad(lambda p:f(p,j),points)) for j in range(2)])


def run():
    plan=bindings();symbolic();mp.mp.dps=plan['mp_decimal_digits']
    a=dict(np.load(ionic.OUT/'inputs.npz'));fermi=dict(np.load(d.plasma.OUT/'stellar-fermi-plasma.npz'))
    cells=a['cells'];eta=fermi['eta'][cells];beta=fermi['beta'][cells];scale=fermi['dimensionless_density_scale'][cells]
    coarse,fine=[integrals(eta,beta,scale,tol) for tol in plan['quadrature_tolerances']]
    score=float(np.max(abs(coarse-fine)/np.maximum(abs(fine),1e-12)))
    controls=[]
    for k in plan['control_positions']:
        ref=independent(eta[k],beta[k],scale[k]);err=float(np.max(abs(ref-fine[:,k])/np.maximum(abs(ref),1e-12)))
        controls.append(dict(cell=int(cells[k]),momentum_reference=ref.tolist(),energy_quadrature=fine[:,k].tolist(),score=err,passed=err<plan['finite_agreement_gate']))
    kernels=[]
    for x in plan['kernel_controls']:
        integral,error=quad(lambda u:u**3/(1+u*u)**3,0,np.inf,weight='sin',wvar=x,epsabs=1e-12,limlst=300)
        calculated=2/(np.pi*x)*integral;exact=np.exp(-x)*(3-x)/8;delta=abs(calculated-exact)
        kernels.append(dict(kappa_r=x,transform=calculated,formula=exact,difference=delta,quad_estimate=2*error/(np.pi*x),passed=delta<plan['kernel_absolute_gate']))
    chi=fine[1]/fine[0];c=d.plasma.constants();density=scale*c['number_prefactor_cm3']*fine[0]
    derivative=scale*c['number_prefactor_cm3']*fine[1];T=a['temperature'];kappa=np.sqrt(4*np.pi*c['e_esu']**2*derivative/(c['k_B']*T))
    counts=a['X']/a['A'];weights=counts/counts.sum(1)[:,None];Z2=weights@(a['Z']**2)
    unit=c['e_esu']**2*kappa*Z2/(c['k_B']*T)
    ne_full=c['N_A']*a['rho']*(counts@a['Z']);density_difference=abs(density/ne_full-1)
    np.savez_compressed(OUT/'states.npz',cells=cells,chi=chi,kappa_cm_inverse=kappa,
        density_reconstruction=density,full_ion_density=ne_full,normalized_integrals=fine,
        self_a0=-unit/2,self_a1=-5*unit/16,self_difference=3*unit/16)
    passed=score<plan['finite_agreement_gate'] and all(v['passed'] for v in controls+kernels)
    save('result.json',dict(classification='Counterexample candidate',cells=len(cells),passed=passed,
        adaptive_agreement=score,independent_controls=controls,kernel_controls=kernels,
        ideal_chi_range=[float(chi.min()),float(chi.max())],kappa_range_cm_inverse=[float(kappa.min()),float(kappa.max())],
        maximum_reconstructed_vs_full_ion_density=float(density_difference.max()),
        self_a0_range=list(map(float,[-unit.max()/2,-unit.min()/2])),
        same_kappa_self_difference_range=list(map(float,[3*unit.min()/16,3*unit.max()/16])),
        whole_physical_EOS_certified=False,native_EOS_replaced=False,
        scope='Same-kappa witness magnitude, not a measured uncertainty interval. Finite k=0 ideal response checks do not certify finite-k response or interacting electron compressibility. No correction added to the existing EOS.'))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());a=dict(np.load(OUT/'states.npz'))
    assert r['passed'] and r['cells']==len(a['cells'])==3206
    assert np.all(a['chi']>0) and np.all(a['chi']<=1) and np.all(a['kappa_cm_inverse']>0)
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS screened Hamiltonian matching and 3206-state k=0 response audit; full dielectric and physical EOS remain open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
