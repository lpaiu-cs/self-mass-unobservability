"""Conditional dispersive photon thermodynamics and actual stellar diagnostics."""
import json, shutil, sys
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.constants import hbar, epsilon_0, elementary_charge, m_e, Avogadro, Boltzmann, c
from scipy.integrate import quad_vec
import gr_radiation_eos_split as split
import lanl_tops_cutoff_reader as cutoff

g=split.g;OUT=g.OUT/'gr-plasma-photon-model'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();split.verify();cutoff.verify()
    for name in ['plasma-radiation2019.pdf','plasma-radiation2019.txt']:
        shutil.copy2(g.ROOT/'outputs'/name,OUT/name)
    constants=dict(hbar=hbar,epsilon_0=epsilon_0,e=elementary_charge,m_e=m_e,N_A=Avogadro,k_B=Boltzmann,c=c)
    save('constants.json',dict(classification='Imported from prior work',source='installed scipy.constants CODATA',values=constants))
    paths=[g.ROOT/'verification/gr_plasma_photon_model.py',split.OUT/'manifest.json',split.reference.OUT/'manifest.json',
           cutoff.OUT/'manifest.json']+list(OUT.iterdir())
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='f305247',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        source_url='https://arxiv.org/abs/1903.05616',source_DOI='10.1103/PhysRevE.100.023202',
        model='Two transverse Bose quasiparticle modes, zero chemical potential, E(p)=sqrt(c^2*p^2+b^2). No longitudinal modes, zero-point terms, nonlocal damping or physical model-error bound. Use the nonrelativistic electron plasma formula b=hbar*sqrt(n_e*e^2/(epsilon0*m_e)) only as a declared diagnostic. Actual n_e comes from the frozen new molecular FreeEOS outputs. Do not substitute the result into the running native EOS or GR structure.',
        thermodynamics='Derive all quantities from the quasiparticle free-energy density. At fixed b, F=-P_kin, U=3P_kin+D. For b(rho,T), the composed Helmholtz function gives U_thermo=U-T*b_T*F_b and P_thermo=P_kin+rho*b_rho*F_b. A cutoff-only opacity substitution does not provide these derivatives or a common matter/radiation free energy.',
        quadrature='Evaluate U,P_kin,D and T*C_V at fixed b at every one of 5735 actual states using momentum-space quadrature at two tolerances. Compare at relative/absolute normalized 1e-8. Independently use 50-decimal energy-space quadrature at a=0,0.1,1,4.3,10 and the minimum/maximum actual plasma-to-temperature ratio.',
        relative_absolute_tolerance=1e-8,controls=[0.,.1,1.,4.3,10.],
        dynamics='Differentiate the canonical Hamiltonian H=N*sqrt(b^2+p_r^2/a^2+L^2/r^2). Keep time-dependent lapse, radial metric and plasma energy. Refractive force and time-dependent quasiparticle energy exchange must be balanced by matter; they cannot be represented by just changing an absorption coefficient.',
        scope='Theorem progress for the declared model plus finite stellar diagnostics. Not physical EOS certification, dense-plasma opacity certification, a refit or actual GR/scalar evolution.'))


def symbolic():
    p,b,T,rho=sp.symbols('p b T rho',positive=True);E=sp.sqrt(p*p+b*b);occupation=1/(sp.exp(E/T)-1)
    density=p*p;f=T*density*sp.log(1-sp.exp(-E/T))
    assert sp.simplify(sp.diff(f,b)-density*b/E*occupation)==0
    assert sp.simplify(f-T*sp.diff(f,T)-density*E*occupation)==0
    ibp=sp.diff(p**3*sp.log(1-sp.exp(-E/T)),p)
    assert sp.simplify(ibp-3*p*p*sp.log(1-sp.exp(-E/T))-p**4/(E*T)*occupation)==0
    assert sp.simplify(density*E*occupation-p**4/E*occupation-b*b*density/E*occupation)==0
    Bt=sp.Function('b')(rho,T);F=sp.Function('F');composed=F(T,Bt)
    dT=sp.diff(composed,T);drho=sp.diff(composed,rho)
    assert sp.simplify(dT-sp.Subs(sp.Derivative(F(T,b),T),b,Bt)-sp.Subs(sp.Derivative(F(T,b),b),b,Bt)*sp.diff(Bt,T))==0
    assert sp.simplify(drho-sp.Subs(sp.Derivative(F(T,b),b),b,Bt)*sp.diff(Bt,rho))==0
    # The energy-space Jacobian supplies sqrt(E^2-b^2), not a hard-cut massless DOS.
    e=sp.symbols('e',positive=True);q=sp.sqrt(e*e-b*b)
    assert sp.simplify(q*q*sp.diff(q,e)-e*q)==0
    r,t,pr,L=sp.symbols('r t pr L',positive=True)
    N=sp.Function('N')(t,r);a=sp.Function('a')(t,r);B=sp.Function('B')(t,r)
    eps=sp.sqrt(B*B+pr*pr/(a*a)+L*L/(r*r));H=N*eps
    velocity=N*pr/(a*a*eps)
    force=-sp.diff(N,r)*eps+N/eps*(pr*pr*sp.diff(a,r)/a**3+L*L/r**3-B*sp.diff(B,r))
    temporal=sp.diff(N,t)*eps+N/eps*(B*sp.diff(B,t)-pr*pr*sp.diff(a,t)/a**3)
    assert sp.simplify(sp.diff(H,pr)-velocity)==0 and sp.simplify(-sp.diff(H,r)-force)==0
    assert sp.simplify(sp.diff(H,t)-temporal)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        free_energy='F(T,b)=C*T*integral_0^infinity p^2*log(1-exp(-sqrt(p^2+b^2)/T)) dp; C=1/(pi^2*hbar^3*c^3), p here has energy units and T means k_B*T.',
        fixed_b='F_b=C*b*integral p^2/E*n_B dp; F=-P_kin; U=3*P_kin+D, D=b*F_b>0 for b>0; T*C_V=C*integral p^2*E^2/T*n_B*(1+n_B) dp.',
        composed_b='For b(rho,T), U_thermo=U-T*b_T*F_b, P_thermo=P_kin+rho*b_rho*F_b. These are the derivatives of the composed quasiparticle Helmholtz density, not a full physical plasma free energy.',
        density_of_states='p^2 dp = E*sqrt(E^2-b^2) dE. At mu=0, U_dispersion < U_vacuum_hard_cut < U_vacuum for b>0. Nonnegative integrands and sqrt(E^2-b^2)<E prove the bounds.',
        ray_H='H=N*sqrt(B^2+p_r^2/a^2+L^2/r^2)',ray_velocity=str(velocity),ray_force=str(force),ray_energy_change=str(temporal),
        boundary='Quasiparticle/matter polarization energy and stress require consistent matching. Changing only grey opacity cannot supply the dispersion, refractive force, state-dependent thermodynamics or actual observation map.'))


def moments(a,tolerance):
    a=np.atleast_1d(np.array(a,float));normal=15/np.pi**4
    def integrand(y):
        E=np.sqrt(y*y+a*a);em=np.exp(-E);f=em/(-np.expm1(-E))
        return normal*np.array([y*y*E*f,y**4/(3*E)*f,a*a*y*y/E*f,y*y*E*E*f*(1+f)])
    return quad_vec(integrand,0,np.inf,epsabs=tolerance,epsrel=tolerance,limit=400)


def independent(a):
    mp.mp.dps=50;a=mp.mpf(str(a));normal=15/mp.pi**4
    def value(e,j):
        z=mp.sqrt(e*e-a*a);f=1/mp.expm1(e)
        return [e*e*z*f,z**3*f/3,a*a*z*f,e**3*z*f*(1+f)][j]
    if not a:return np.array([1,1/3,0,4.])
    return np.array([float(normal*mp.quad(lambda x:value(x,j),[a,a+1,a+4,a+16,mp.inf])) for j in range(4)])


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    symbolic();state=dict(np.load(split.reference.OUT/'reference-state.npz'))
    native=np.concatenate([np.load(split.reference.OUT/f'block-{i:04}.npz')['new'] for i in range(0,len(state['X']),128)])
    assert len(native)==5735
    rho=np.exp(state['lnd']);T=np.exp(state['lnT']);ne_cm3=Avogadro*native[:,13];assert np.all(ne_cm3>0)
    b=hbar*np.sqrt(ne_cm3*1e6*elementary_charge**2/(epsilon_0*m_e));ratio=b/(Boltzmann*T)
    values,error=moments(ratio,1e-9);tight,tighterror=moments(ratio,5e-12)
    score=float(np.max(abs(values-tight)/np.maximum(abs(tight),1)))
    trace=float(np.max(abs(tight[0]-3*tight[1]-tight[2])))
    controls=[]
    for x in plan['controls']+[float(ratio.min()),float(ratio.max())]:
        v,err=moments(x,5e-12);v=v[:,0];ref=independent(x);difference=float(np.max(abs(v-ref)/np.maximum(abs(ref),1)))
        controls.append(dict(a=x,momentum_space=v.tolist(),energy_space_50decimal=ref.tolist(),maximum_error=difference,
            passed=difference<plan['relative_absolute_tolerance']))
    arad_si=np.pi**2*Boltzmann**4/(15*hbar**3*c**3);arad_cgs=10*arad_si
    native_arad=json.loads((split.OUT/'result.json').read_text())['compiled_a_rad_cgs']
    uvac=arad_cgs*T**4;u,p,d,cv=tight*uvac
    # Finite momentum/energy reference only; actual electrons can be degenerate/relativistic.
    xF=hbar*(3*np.pi**2*ne_cm3*1e6)**(1/3)/(m_e*c)
    fermi=m_e*c*c*xF*xF/(np.sqrt(1+xF*xF)+1);theta=Boltzmann*T/fermi
    positive=bool(np.all(u>0) and np.all(p>0) and np.all(d>0) and np.all(cv>0))
    passed=bool(positive and score<plan['relative_absolute_tolerance'] and trace<plan['relative_absolute_tolerance'] and all(x['passed'] for x in controls))
    np.savez_compressed(OUT/'stellar-photon-diagnostic.npz',rho_B=rho,T_K=T,ne_cm3=ne_cm3,plasma_energy_J=b,
        plasma_to_temperature=ratio,U_vacuum_cgs=uvac,U_quasiparticle_cgs=u,P_kinetic_cgs=p,trace_deficit_cgs=d,
        CvT_fixed_gap_cgs=cv,fermi_momentum_over_mec=xF,kT_over_ideal_fermi_energy=theta,dm=state['dm'])
    selected=[]
    for i in [0,1175,2972,3043,4352,5734]:
        selected.append(dict(cell=i,T_K=float(T[i]),rho_B=float(rho[i]),plasma_to_temperature=float(ratio[i]),
            quasiparticle_U_over_vacuum=float(tight[0,i]),kinetic_P_over_vacuum_P=float(3*tight[1,i]),
            xF=float(xF[i]),kT_over_ideal_fermi_energy=float(theta[i])))
    save('result.json',dict(classification='Counterexample candidate',completed=True,all_finite_gates_passed=passed,
        cells=len(native),plasma_to_temperature_range=[float(ratio.min()),float(ratio.max())],
        quadrature_difference=score,trace_identity_maximum_residual=trace,quadrature_error_estimates=[float(error),float(tighterror)],
        controls=controls,selected_states=selected,radiation_constant_relative_to_native=float(arad_cgs/native_arad-1),
        maximum_kinetic_pressure_change_over_total_native_pressure=float(np.max(abs(p-uvac/3)/native[:,1])),
        maximum_fixed_gap_capacity_change_over_total_native_capacity=float(np.max(abs(cv-4*uvac)/(rho*native[:,10]))),
        cell_uniform_proper_energy_difference_erg=float(np.sum(state['dm']/rho*(u-uvac))),
        energy_difference_is_time_release=False,nonrelativistic_plasma_formula_physically_certified=False,
        thermodynamic_derivatives_of_actual_plasma_gap_supplied=False,matter_polarization_double_count_resolved=False,
        replaced_native_EOS_or_GR=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS dispersive photon model bindings; finite/physical/actual evolution remain separate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
