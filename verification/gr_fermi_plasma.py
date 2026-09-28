"""Match frozen native electron states to finite-T leading-order plasma moments."""
import json, shutil, sys
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.integrate import quad_vec
from scipy.special import expit
import gr_plasma_photon_runner as photon

g=photon.g;split=photon.original.split;OUT=g.OUT/'gr-fermi-plasma'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def constants():
    # Match the frozen native CGS source, including its derived electron mass.
    text=(split.OUT/'mod_free_eos_constants.f90').read_text()
    c=2.99792458e10;h=6.62607015e-27;na=6.02214076e23;k=1.380649e-16
    for token in ['109737.31568160_fp_kind','3.1415926535897932385_fp_kind',
                  '(planck/(echarge*echarge))*','(planck/echarge)*(planck/echarge)*clight/(2._fp_kind*pi*pi)']:
        assert token in text
    e=1.602176634e-19*.1*c;hb=h/(2*np.pi)
    me=109737.31568160*(h/(e*e))*(h/e)**2*c/(2*np.pi**2)
    return dict(c=c,h=h,hbar=hb,N_A=na,k_B=k,e_esu=e,m_e_g=me,
                alpha=e*e/(hb*c),mc2_erg=me*c*c,number_prefactor_cm3=(me*c/hb)**3/np.pi**2)


def prepare():
    assert not OUT.exists();OUT.mkdir();photon.verify()
    for name in ['braaten-segel1993.pdf','braaten-segel1993.txt']:
        shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    save('constants.json',dict(classification='Imported from prior work',values=constants(),
        source='Frozen native mod_free_eos_constants.f90. Its CGS charge and Rydberg-derived mass convention is retained; no experimental-constant uncertainty certificate.'))
    paths=[g.ROOT/'verification/gr_fermi_plasma.py',split.OUT/'mod_free_eos_constants.f90',
           split.OUT/'mod_free_eos.f90',split.reference.OUT/'manifest.json',photon.OUT/'manifest.json']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='8e41192',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source='https://arxiv.org/abs/hep-ph/9302213',doi='10.1103/PhysRevD.48.1478',equations=[1,4,5,6,11,17,18,19,20,21],
        model='Ideal relativistic electron/positron distributions and leading-order alpha plasma moments. Use the actual saved native kinetic eta and n_e=N_A*rmue; test their agreement before any density inversion or EOS substitution.',
        quadrature='All 5735 states, kinetic-energy coordinate t=(E-mc2)/(kT); two adaptive vector quadratures. Six actual states and three synthetic controls independently integrated in momentum coordinate at 50 decimal digits.',
        actual_controls=[0,1175,2972,3043,4352,5734],synthetic_controls=[[-20,1e-5],[0,.1],[3,3]],
        quadrature_relative_tolerance=1e-8,native_density_relative_tolerance=1e-8,
        moment_bounds='w_p^2<=m_t^2<=1.5*w_p^2; 0<=w_1^2<=w_p^2; k_max^2>=w_p^2. Approximate Braaten-Segel high-k moments are compared to direct integral moments with errors reported, not fitted or retrospectively thresholded.',
        cold_controls=[.01,.4,3.],scope='Finite state consistency and moment diagnostics; symbolic conditional identities. Not physical plasma EOS certification, a full dispersion/transport solution or a change to the running native GR calculation.'))


def symbolic():
    t,b,x=sp.symbols('t b x',positive=True)
    gamma=1+b*t;p=sp.sqrt(b*t*(2+b*t));v2=p*p/gamma**2
    assert sp.simplify(p*p*sp.diff(p,t)-b*gamma*p)==0
    assert sp.simplify(p*p/gamma*sp.diff(p,t)-b*p)==0
    # Exact cold-limit primitives, with momentum x in units of m_e*c.
    gamma=sp.sqrt(1+x*x);v2=x*x/gamma**2
    assert sp.simplify(sp.diff(x**3/(3*gamma),x)-x*x/gamma*(1-v2/3))==0
    assert sp.simplify(sp.diff(x**5/(3*gamma**3),x)-x*x/gamma*(sp.Rational(5,3)*v2-v2*v2))==0
    r=sp.symbols('r',positive=True)
    G=3/(2*r*r)*(1-(1-r*r)*sp.atanh(r)/r)
    K=3/(r*r)*(sp.atanh(r)/r-1)
    assert sp.simplify(2*G+(1-r*r)*K-3)==0
    series=sum(sp.Rational(3,(2*j+1)*(2*j+3))*r**(2*j) for j in range(6))
    assert sp.series(G-series,r,0,12).removeO()==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        cold_limit='I_p=xF^3/(3*gammaF), I_1=xF^5/(3*gammaF^3), hence v_star=vF and w_p^2/w_Drude^2=1/gammaF at fixed density.',
        kernel='G(r)=3*sum_{j>=0} r^(2j)/[(2j+1)(2j+3)], 0<=r<1. Positive coefficients imply G(0)=1, G increasing, G(1-)=3/2.',
        tail='After terms j=0..N, 0<=G-G_N<=3*r^(2N+2)/[(2N+3)(2N+5)*(1-r^2)]. This exact-real series bound is not a floating-point roundoff bound.',
        root='For q>=0 and 0<=v<1, F(u)=u^2-q^2-G(v*q/u) increases strictly on [sqrt(q^2+1),sqrt(q^2+G(v))]. Its endpoints bracket a unique transverse root. F_u=2u+G_prime(r)*r/u>0.',
        approximation_boundary='This root theorem concerns the approximate Braaten-Segel kernel, not all physical plasma modes. Exact leading-order integral moments and approximate high-k moments are reported separately.'))


def atanh_over_v(v):
    v=np.asarray(v);z=v*v
    # Direct expression loses the small correction at small v; keep it explicitly.
    return np.where(v<.02,1+z*(1/3+z*(1/5+z*(1/7+z*(1/9+z/11)))),np.arctanh(v)/v)


def integrals(eta,beta,scale,tol):
    eta,beta,scale=np.broadcast_arrays(np.atleast_1d(eta),np.atleast_1d(beta),np.atleast_1d(scale))
    def fn(t):
        gamma=1+beta*t;x=np.sqrt(beta*t*(2+beta*t));v=x/gamma;v2=v*v
        fe=expit(eta-t);fp=expit(-eta-t-2/beta);weight=beta*x*(fe+fp)/scale
        return np.array([beta*gamma*x*(fe-fp)/scale,weight*(1-v2/3),weight,
            weight*((5/3)*v2-v2*v2),weight*(2*np.arcsinh(x)/v-1),beta*gamma*x*fp/scale])
    return quad_vec(fn,0,np.inf,epsabs=tol,epsrel=tol,limit=600)


def independent(eta,beta,scale):
    mp.mp.dps=50;eta,beta,scale=map(lambda x:mp.mpf(str(x)),[eta,beta,scale])
    def fn(x,j):
        gamma=mp.sqrt(1+x*x);v=x/gamma;v2=v*v;t=x*x/((gamma+1)*beta)
        fe=1/(1+mp.exp(t-eta));fp=1/(1+mp.exp(t+eta+2/beta));weight=x*x/gamma*(fe+fp)
        values=[x*x*(fe-fp),weight*(1-v2/3),weight,weight*(mp.mpf(5)/3*v2-v2*v2),
                weight*(2*mp.asinh(x)/v-1) if x else mp.mpf(0),x*x*fp]
        return values[j]/scale
    points=[mp.mpf(0)]+[mp.sqrt(beta*t*(2+beta*t)) for t in [1,max(2,float(eta)),max(2,float(eta))+8,max(2,float(eta))+40]]+[mp.inf]
    return np.array([float(mp.quad(lambda x:fn(x,j),points)) for j in range(6)])


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    symbolic();c=constants();state=dict(np.load(split.reference.OUT/'reference-state.npz'))
    native=np.concatenate([np.load(split.reference.OUT/f'block-{i:04}.npz')['new'] for i in range(0,5735,128)])
    eta=native[:,12];T=np.exp(state['lnT']);ne=c['N_A']*native[:,13]
    beta=c['k_B']*T/c['mc2_erg'];scale=ne/c['number_prefactor_cm3']
    first,e1=integrals(eta,beta,scale,1e-9);last,e2=integrals(eta,beta,scale,2e-12)
    difference=float(np.max(abs(first-last)/np.maximum(abs(last),1e-10)))
    controls=[]
    entries=[(int(i),eta[i],beta[i],scale[i]) for i in plan['actual_controls']]
    entries += [('synthetic',a,b,b**1.5*np.exp(min(a,0))) for a,b in plan['synthetic_controls']]
    for index,a,b,s in entries:
        val,err=integrals(a,b,s,2e-12);ref=independent(a,b,s)
        score=float(np.max(abs(val[:,0]-ref)/np.maximum(abs(ref),1e-10)))
        controls.append(dict(cell=index,eta=float(a),beta=float(b),energy_coordinate=val[:,0].tolist(),
            momentum_50decimal=ref.tolist(),relative_error=score,passed=score<plan['quadrature_relative_tolerance']))
    n,ip,im,i1,ik,pairs=last;v2=i1/ip;v=np.sqrt(v2);aov=atanh_over_v(v)
    # Stable positive series for the small-velocity characteristic kernel.
    G=np.where(v<.02,1+v2/5+3*v2*v2/35+v2**3/21+v2**4/33,3/(2*v2)*(1-(1-v2)*aov))
    K=np.where(v<.02,1+3*v2/5+3*v2*v2/7+v2**3/3+3*v2**4/11,3/v2*(aov-1))
    plasma=c['mc2_erg']*np.sqrt(4*c['alpha']/np.pi*ip*scale)
    drude=c['hbar']*np.sqrt(4*np.pi*ne*c['e_esu']**2/c['m_e_g'])
    density_error=abs(n-1);finite=bool(np.all(np.isfinite(last)) and np.all(ip>0) and np.all(v2>=0) and np.all(v2<=1)
        and np.all(im>=ip) and np.all(im<=1.5*ip) and np.all(ik>=ip))
    cold=[]
    for x in plan['cold_controls']:
        def f(y):
            gamma=np.sqrt(1+y*y);v2=y*y/(1+y*y)
            return np.array([y*y/gamma*(1-v2/3),y*y/gamma*((5/3)*v2-v2*v2)])
        val,err=quad_vec(f,0,x,epsabs=1e-14,epsrel=1e-12)
        ref=np.array([x**3/(3*np.sqrt(1+x*x)),x**5/(3*(1+x*x)**1.5)])
        score=float(np.max(abs(val/ref-1)));cold.append(dict(xF=x,relative_error=score,passed=score<1e-10))
    selected=[dict(cell=int(i),eta=float(eta[i]),T_K=float(T[i]),ne_cm3=float(ne[i]),density_relative_error=float(n[i]-1),
        v_star=float(v[i]),plasma_to_Drude=float(plasma[i]/drude[i]),plasma_to_kT=float(plasma[i]/(c['k_B']*T[i])),
        exact_mt2_over_wp2=float(im[i]/ip[i]),approx_mt2_relative_error=float(G[i]*ip[i]/im[i]-1),
        exact_kmax2_over_wp2=float(ik[i]/ip[i]),approx_kmax2_relative_error=float(K[i]*ip[i]/ik[i]-1)) for i in plan['actual_controls']]
    np.savez_compressed(OUT/'stellar-fermi-plasma.npz',eta=eta,beta=beta,T_K=T,ne_native_cm3=ne,
        density_reconstructed_over_native=n,dimensionless_density_scale=scale,normalized_moments=last,
        plasma_energy_erg=plasma,Drude_energy_erg=drude,v_star=v,approx_mt2_over_wp2=G,approx_kmax2_over_wp2=K)
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(eta),
        quadrature_passed=bool(difference<plan['quadrature_relative_tolerance'] and all(r['passed'] for r in controls)),
        native_density_passed=bool(np.max(density_error)<plan['native_density_relative_tolerance']),moment_bounds_passed=finite,
        cold_controls=cold,controls=controls,quadrature_relative_difference=difference,quadrature_error_estimates=[float(e1),float(e2)],
        maximum_native_density_relative_error=float(np.max(density_error)),maximum_density_error_cell=int(np.argmax(density_error)),
        eta_range=[float(eta.min()),float(eta.max())],beta_range=[float(beta.min()),float(beta.max())],
        plasma_to_Drude_range=[float(np.min(plasma/drude)),float(np.max(plasma/drude))],
        v_star_range=[float(v.min()),float(v.max())],max_approx_mt2_relative_error=float(np.max(abs(G*ip/im-1))),
        max_approx_kmax2_relative_error=float(np.max(abs(K*ip/ik-1))),selected_states=selected,
        maximum_pair_to_net_density=float(np.max(pairs/n)),native_EOS_replaced=False,physical_plasma_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS Fermi plasma provenance; consult separate numerical and physical verdicts',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
