"""Match the classical ring limit and expose finite-branch thermodynamic boundaries."""
import json, shutil, sys
import mpmath as mp
import numpy as np
import sympy as sp
import gr_fermi_plasma as plasma

g=plasma.g;split=plasma.split;OUT=g.OUT/'gr-plasma-free-energy-matching'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();plasma.verify()
    for name in ['coulomb.f90','mod_pi_fit.f90']:
        shutil.copy2(split.model.CACHE/'source/src'/name,OUT/name)
    shutil.copy2(split.model.CACHE/'source/document_FreeEOS/eos_papers/coulomb/coulomb.tex',OUT/'coulomb.tex')
    shutil.copy2(split.OUT/'mod_free_eos.f90',OUT/'mod_free_eos.f90')
    shutil.copy2(split.model.OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    for name in ['debye-huckel-two-component2017.pdf','debye-huckel-two-component2017.txt']:
        shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    paths=[g.ROOT/'verification/gr_plasma_free_energy_matching.py',plasma.OUT/'manifest.json',split.model.OUT/'manifest.json',
           split.reference.OUT/'manifest.json',g.ROOT/'verification/direct_eos_gr.py']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='fbcfefe',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        sources=['https://arxiv.org/abs/1704.06502','https://freeeos.sourceforge.net/coulomb.pdf','https://arxiv.org/abs/hep-ph/9302213'],
        original_freeeos_pdf_download='HTTP 403 from urllib; retain native repository Coulomb paper TeX and active implementation instead. Public PDF was independently readable through web retrieval.',
        target='Prove the classical weak-coupling Coulomb ring functional matches the already-present native Debye-Huckel free-energy term. Derive the finite longitudinal-branch boundary term before any EOS substitution. Count actual saved native Coulomb regimes without treating their fit as physically certified.',
        numerical_control='80-decimal independent logarithmic ring integral against its analytic value, absolute tolerance 1e-30. Symbolic DH thermodynamic derivatives, native algebra matching, and finite-boundary integration by parts.',
        scope='Conditional matching/no-double-count boundary, not a complete RPA/free-energy replacement, a fitted dense-plasma EOS, atmosphere, transport or GR evolution.'))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    source=(OUT/'coulomb.f90').read_text();fit=(OUT/'mod_pi_fit.f90').read_text()
    assert 'gc = lambda/3._fp_kind' in source and 'fcoulomb = -sum0*boltzmann*t*gc' in source
    assert 'gamma = (lambda*lambda/3._fp_kind)**(1._fp_kind/3._fp_kind)' in source
    assert 'xdh10 = -0.4_fp_kind, xmocp10 = 0._fp_kind' in fit
    eos=(OUT/'mod_free_eos.f90').read_text();bridge=(OUT/'direct_ion_bridge.f90').read_text()
    assert 'call free_eos(0,3,1,-2,kif' in bridge
    assert '! EOS1\n          ifcoulomb = 5\n          ifpi = 3' in eos
    assert '! EOS1 without radiation pressure\n          ifcoulomb = 5\n          ifpi = 3\n          ifrad = 0' in eos
    k,K,T,e,S0,S2,n,A=sp.symbols('k K T e S0 S2 n A',positive=True)
    ring=k*k*sp.log(1+K*K/(k*k))-K*K
    assert sp.simplify(sp.diff(ring,K)+2*K**3/(k*k+K*K))==0
    integral=sp.integrate(1/(k*k+K*K),(k,0,sp.oo));assert sp.simplify(integral-sp.pi/(2*K))==0
    ring_free=-T*K**3/(12*sp.pi)
    lam=2*e**3*sp.sqrt(sp.pi)*S2**sp.Rational(3,2)/(T**sp.Rational(3,2)*S0)
    kappa=sp.sqrt(4*sp.pi*e*e*S2/T)
    assert sp.simplify(-T*S0*lam/3-ring_free.subs(K,kappa))==0
    f=-A*n**sp.Rational(3,2)/sp.sqrt(T)
    U=f-T*sp.diff(f,T);P=n*sp.diff(f,n)-f
    assert sp.simplify(U-sp.Rational(3,2)*f)==0 and sp.simplify(P-f/2)==0
    assert sp.simplify(-sp.diff(f,T)-f/(2*T))==0
    assert sp.simplify(sp.diff(U,T)+3*f/(4*T))==0
    eps=sp.Function('epsilon')(k);L=sp.log(1-sp.exp(-eps/T));occ=1/(sp.exp(eps/T)-1)
    assert sp.simplify(sp.diff(k**3*L,k)-3*k*k*L-k**3*sp.diff(eps,k)*occ/T)==0
    x=sp.symbols('x',positive=True);boundary=sp.Function('K')(x);integrand=sp.Function('a')(k,x)
    lhs=sp.diff(sp.Integral(integrand,(k,0,boundary)),x)
    rhs=sp.Integral(sp.diff(integrand,x),(k,0,boundary))+integrand.subs(k,boundary)*sp.diff(boundary,x)
    assert sp.simplify(lhs-rhs)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        ring='f_ring=(T/2)*integral d^3k/(2*pi)^3 [log(1+kappa^2/k^2)-kappa^2/k^2]=-T*kappa^3/(12*pi), with T in energy units. Differentiate the convergent integral in kappa and fix f(0)=0. The unsubtracted first-order term has a linear ultraviolet divergence.',
        native_match='kappa^2=4*pi*e^2*S2/T; lambda=2*e^3*sqrt(pi)*S2^(3/2)/(T^(3/2)*S0). Then -T*S0*lambda/3 is exactly the same ring free-energy density. S0,S2 are the declared effective Coulomb sums; the joint classical weak-coupling physical limit uses the full charged-species sums.',
        duplicate='If f_native=f_ideal+f_DH+o(f_DH) and a separately added correlation correction has f_new=f_DH+o(f_DH), then f_native+f_new=f_ideal+2*f_DH+o(f_DH). A residual matched correction must subtract the common DH term at free-energy level, including its state/composition derivatives.',
        thermodynamics='At fixed classical composition kappa^2 proportional to n/T: f_DH=-A*n^(3/2)*T^(-1/2), U_DH=3*f_DH/2, P_DH=f_DH/2, S_DH=f_DH/(2*T), C_V_DH=-3*f_DH/(4*T).',
        finite_branch='For f_L=C*T*integral_0^K k^2 log(1-exp(-epsilon(k)/T))dk, integration by parts gives f_L=C*T*K^3*log(1-exp(-epsilon(K)/T))/3-P_kin,L. Thus f_L=-P_kin,L is false at a finite thermally populated endpoint.',
        moving_boundary='For state x, d_x integral_0^K(x) a(k,x)dk = integral_0^K partial_x a dk+a(K,x)*K_x. At the timelike longitudinal endpoint epsilon(K)=K>0 and T>0, the Bose logarithm is nonzero. The endpoint cannot silently be dropped from thermodynamic derivatives.',
        limits='The time-like branch cutoff from a neutrino-decay calculation does not specify the entire equilibrium spectral determinant or damped/spacelike response. These identities do not supply the missing physical continuum or a full consistent plasma free energy.'))
    mp.mp.dps=80
    def original(y):
        if not y:return -mp.mpf(1)
        if y<=2:return y*y*mp.log1p(1/(y*y))-1
        term=1/(y*y);total=mp.mpf(0)
        for j in range(2,202):total+=(-1)**(j+1)*term/j;term/=y*y
        return total
    value=mp.quad(original,[0,1,2,4,mp.inf]);expected=-mp.pi/3;error=abs(value-expected)
    # The integrated alternating-series remainder above y=2 is explicitly bounded.
    series_tail=mp.mpf(2)**(3-2*202)/(202*(2*202-3))
    passed=bool(error+series_tail<mp.mpf('1e-30'))
    state=np.load(split.reference.OUT/'reference-state.npz')
    native=np.concatenate([np.load(split.reference.OUT/f'block-{i:04}.npz')['new'] for i in range(0,5735,128)])
    lam_native=native[:,19];assert np.all(lam_native>0)
    gamma=(lam_native*lam_native/3)**(1/3);loggamma=np.log10(gamma)
    masks={'DH':loggamma<=-.4,'transition':(loggamma>-.4)&(loggamma<0),'modified_OCP':loggamma>=0}
    census={name:dict(cells=int(mask.sum()),baryon_mass_fraction=float(np.sum(state['dm'][mask])/np.sum(state['dm']))) for name,mask in masks.items()}
    assert sum(x['cells'] for x in census.values())==5735
    np.savez_compressed(OUT/'native-coulomb-regimes.npz',lambda_native=lam_native,gamma_native=gamma,log10_gamma=loggamma)
    save('result.json',dict(classification='Counterexample candidate',completed=True,numerical_control_passed=passed,
        logarithmic_ring_integral=str(value),analytic_ring_integral=str(expected),finite_numerical_error=str(error),
        high_wave_number_series_tail_upper=str(series_tail),cells=5735,regimes=census,
        gamma_range=[float(gamma.min()),float(gamma.max())],active_native_coulomb_option=5,active_native_exchange_option=14,
        no_additive_DH_recount=True,new_physical_plasma_free_energy_completed=False,native_EOS_replaced=False,full_GR_evolution=False))
    assert passed
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['numerical_control_passed']
    print('PASS conditional Coulomb matching and finite-branch boundary; complete plasma free energy remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
