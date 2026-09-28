"""Energy-resolved electron-ion response on the saved stellar background.

Declared model: nonmagnetic, relativistic ideal electrons, classical Debye
ions, screened Born elastic collisions and zero electric current. No fitted
relaxation time. Correlations, non-Born and electron-electron terms are absent.
"""
from pathlib import Path
import argparse
import json
import time
import urllib.request
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.constants import hbar,m_e,c,k,epsilon_0,elementary_charge,Avogadro
from scipy.special import expit
from scipy.integrate import quad
import def_free_surface_thermal as thermal

h=thermal.h
OUT=thermal.OUT.parent/'def-electron-collision-response'
BENCHMARK=thermal.surface.old.exterior.s.OUT/'companion-benchmark.json'
E2=elementary_charge**2/(4*np.pi*epsilon_0)


def inputs():
    d=np.load(thermal.OUT/'coefficients.npz');rho=d['raw'][:,0];T=np.exp(d['lnT'])
    ne=Avogadro*d['raw'][:,13]*1e6;ion=rho[:,None]*Avogadro*d['X']/thermal.g.c.A*1e6
    Z=thermal.g.c.Z;full=ion@Z
    x=hbar*(3*np.pi**2*ne)**(1/3)/(m_e*c);ef=m_e*c*c*x*x/(np.sqrt(1+x*x)+1)
    theta=k*T/ef;deficit=1-ne/full
    # These cohort cuts define this experiment, not universal accuracy bounds.
    selected=(theta<.1)&(abs(deficit)<1e-4)
    return d,dict(T=T,ne=ne,ion=ion,theta=theta,ionization_deficit=deficit,selected=selected,xF=x)


def angular(u,v2):
    u=np.asarray(u,float);large=u>=8;small=np.where(large,1,u)
    logarithm=np.log1p(1/small)
    l1=.5*(logarithm-1/(1+small))
    l2=.5-small*logarithm+small/(2*(1+small))
    # The convergent large-u series avoids cancellation at low momentum.
    if np.any(large):
        ul=u[large];l1[large]=0;l2[large]=0
        for j in range(18):
            term=.5*(-1.)**j*(j+1)/ul**(j+2)
            l1[large]+=term/(j+2);l2[large]+=term/(j+3)
    return l1-v2*l2


def distribution(T,ne,order):
    t=k*T/(m_e*c*c);xf=hbar*(3*np.pi**2*ne)**(1/3)/(m_e*c)
    eta=xf*xf/(np.sqrt(1+xf*xf)+1)/t
    gx,gw=leggauss(order)
    # x=(epsilon-mc^2)/(kT)=u^2 removes the square-root endpoint.
    top=np.sqrt(eta+60);u=top[:,None]*(gx+1)/2;x=u*u
    dx=top[:,None]*gw*u;pp=np.sqrt(t[:,None]*x*(2+t[:,None]*x))
    density=(m_e*c/hbar)**3/np.pi**2*t[:,None]*(1+t[:,None]*x)*pp
    for _ in range(12):
        f=expit(eta[:,None]-x);slope=f*(1-f)
        number=np.sum(density*f*dx,axis=1);derivative=np.sum(density*slope*dx,axis=1)
        error=number/ne-1
        if np.max(abs(error))<2e-13:break
        eta-=(number-ne)/derivative
    assert np.max(abs(error))<2e-11,float(np.max(abs(error)))
    return eta,x,pp,slope,dx,derivative/(k*T),float(np.max(abs(error)))


def kinetic(indices,order):
    d,data=inputs();T=data['T'][indices];ne=data['ne'][indices];ions=data['ion'][indices]
    eta,x,pdim,fp,dx,dndmu,density_error=distribution(T,ne,order)
    Z=thermal.g.c.Z;charge2=ions@(Z*Z)
    qe2=elementary_charge**2/epsilon_0*dndmu
    qi2=elementary_charge**2/epsilon_0*charge2/(k*T);qs2=qe2+qi2
    momentum=m_e*c*pdim;energy=m_e*c*c*np.sqrt(1+pdim*pdim);velocity=momentum*c*c/energy
    u=hbar*hbar*qs2[:,None]/(4*momentum*momentum)
    coulomb=angular(u,(velocity/c)**2)
    nu=4*np.pi*E2**2*charge2[:,None]/(momentum*momentum*velocity)*coulomb
    assert np.min(nu)>0 and np.all(np.isfinite(nu))
    tau=1/nu;z=x-eta[:,None]
    W=momentum**3/(3*np.pi**2*hbar**3)*c*c/energy*fp*dx
    weights=W*tau;moments=np.array([np.sum(weights*z**j,axis=1) for j in range(3)])
    mean=moments[1]/moments[0];variance=np.sum(weights*(z-mean[:,None])**2,axis=1)
    conductivity=k*k*T*variance
    memory=np.sum(weights*tau*(z-mean[:,None])**2,axis=1)/variance
    gamma=Z[None,:]**(5/3)*E2*(4*np.pi*ne[:,None]/3)**(1/3)/(k*T[:,None])
    # Report the scattering-weighted share outside weak coupling, not a
    # probability that the microscopic approximation is physically correct.
    strong_share=np.sum(ions*(Z*Z)[None,:]*(gamma>1),axis=1)/charge2
    actual=d['K'][indices,1]*1e-5  # erg/(cm s K) -> W/(m K)
    return dict(indices=indices,order=order,eta=eta,x=x,z=z,weights=weights,tau=tau,
                moments=moments,variance=variance,conductivity_SI=conductivity,tau_eff=memory,
                native_conductivity_SI=actual,density_relative_error=density_error,
                strong_ion_scattering_share=strong_share,screening_wave_number=np.sqrt(qs2))


def response(value,omega):
    """Relative K(i omega)-K(0), without subtracting nearly equal K values."""
    omega=np.broadcast_to(omega,value['tau_eff'].shape)
    q=1j*omega[:,None]*value['tau'];factor=-q/(1+q)
    delta=np.array([np.sum(value['weights']*value['z']**j*factor,axis=1) for j in range(3)])
    a,b,_=value['moments'];da,db,dc=delta
    correction=dc-((2*b*db+db*db)*a-b*b*da)/(a*(a+da))
    return correction/value['variance']


def symbolic():
    A,B,C,Ap,Bp,Cp=sp.symbols('A B C Ap Bp Cp');s=sp.symbols('s')
    delta=sp.diff(C-s*Cp-(B-s*Bp)**2/(A-s*Ap),s).subs(s,0)
    assert sp.simplify(delta+Cp-2*B*Bp/A+B*B*Ap/A**2)==0
    u,z,beta=sp.symbols('u z beta',positive=True)
    primitive1=sp.log(z+u)+u/(z+u)
    primitive2=z-2*u*sp.log(z+u)-u*u/(z+u)
    assert sp.simplify(sp.diff(primitive1,z)-z/(z+u)**2)==0
    assert sp.simplify(sp.diff(primitive2,z)-z*z/(z+u)**2)==0
    errors=[]
    for a in [1e-8,.01,1.,7.9,8.,100.,1e6]:
        for b in [0.,.16,.9]:
            exact=quad(lambda q:q*(1-b*q)/(q+a)**2/2,0,1,epsabs=1e-25,epsrel=2e-12,limit=300)[0]
            got=angular(np.array([a]),np.array([b]))[0];errors.append(abs(got/exact-1))
    assert max(errors)<2e-10,max(errors)
    return dict(classification='Proven',passed=True,angular_quadrature_max_relative=max(errors),
        moment_operator='L_j(s)=integral W(epsilon)*tau/(1+s*tau)*[(epsilon-mu)/(kT)]^j; K(s)=k_B^2*T*(L_2-L_1^2/L_0). Zero electric current uses the frequency-dependent Schur complement.',
        first_memory='-K_prime(0)/K(0)=integral W*tau^2*(z-L1/L0)^2 / integral W*tau*(z-L1/L0)^2 > 0.',
        passivity='The two-current moment matrix is a positive sum of rank-one matrices divided by 1+s*tau. Its Hermitian part is positive for Re(s)>=0; the zero-current Schur complement retains nonnegative real part.',
        scope='Identities for this elastic diagonal collision model; no certification of its plasma approximation.')


def prepare():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    url='https://arxiv.org/html/astro-ph/9604130'
    with urllib.request.urlopen(url,timeout=20) as reply:(OUT/'potekhin-yakovlev1996.html').write_bytes(reply.read())
    d,data=inputs();indices=np.flatnonzero(data['selected']);assert len(indices)>0
    selected=indices[np.unique(np.linspace(0,len(indices)-1,32).astype(int))]
    paths=[Path(__file__),thermal.OUT/'coefficients.npz',Path(thermal.__file__),BENCHMARK,OUT/'potekhin-yakovlev1996.html']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='2ceb61f0',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        primary_source=dict(classification='Imported from prior work',url=url,
            equations='Nonmagnetic limit of sections 3.2 and 4.2-4.3, transport cross sections and statistical energy averaging. The Born screened potential is the declared approximation; full stellar validity is not inherited.'),
        model='B=0 ideal relativistic Fermi electrons matched to actual native free-electron density, full nuclear charges within selected near-full-ionization cohort, classical Debye ion screening plus ideal electron compressibility, elastic screened Born scattering summed over ions, no electron-electron collision or correlated ion structure factor.',
        selection='T/T_F<0.1 and absolute fractional native free-electron deficit relative to fully stripped nuclei <1e-4. These are declared cohort cuts, not physical error certificates.',
        decision='Construct actual energy-dependent thermal memory without an adjustable tau; compare the static conductivity to the saved native value before replacing a full-star closure. A mismatch is preserved and cannot be fitted away.',
        gates=dict(density_relative=2e-11,angular_relative=2e-10,energy_quadrature_relative=1e-6,native_conductivity_compatibility_relative=.2),
        budget=dict(pilot_cells=len(selected),production_orders=[128,256],production_hard_seconds=60,native_calls=0,new_stellar_time_steps=0,CPU_workers=1,automatic_expansion=False)))
    np.savez_compressed(OUT/'domain.npz',**data,dm=d['dm'])
    began=time.monotonic();pilot=kinetic(selected,64);elapsed=time.monotonic()-began
    h.write(OUT/'pilot.json',dict(seconds=elapsed,selected_cells=len(indices),forecast_seconds=elapsed*len(indices)/len(selected)*6,
        estimate='Linear cohort size and quadrature-node scaling; Python import and symbolic overhead excluded. Hard limit includes both.',total_seconds=time.monotonic()-start))
    h.write(OUT/'symbolic.json',symbolic());print('PILOT',elapsed,'COHORT',len(indices),flush=True)


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic()
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert h.digest(h.ROOT/rel)==sha,rel
    pilot=json.loads((OUT/'pilot.json').read_text());assert pilot['forecast_seconds']<50,pilot
    d,data=inputs();indices=np.flatnonzero(data['selected']);coarse=kinetic(indices,128);fine=kinetic(indices,256)
    differences={key:float(np.max(abs(coarse[key]/fine[key]-1))) for key in ['conductivity_SI','tau_eff']}
    benchmark=json.loads(BENCHMARK.read_text())
    period=benchmark['leading_drive']['period_seconds'];frequencies=np.arange(1,4)*2*np.pi/period
    responses=np.array([response(fine,w/(d['A'][indices]*d['N'][indices])) for w in frequencies])
    controls=[]
    for factor in [.01,1.,100.]:
        omega=factor/fine['tau_eff'];v=response(fine,omega);ref=response(coarse,factor/coarse['tau_eff'])
        assert np.min((1+v).real)>0
        controls.append(dict(scaled_frequency=factor,quadrature_max_absolute=float(np.max(abs(v-ref))),
                             non_Drude_max=float(np.max(abs(v-(-1j*factor)/(1+1j*factor))))))
    mismatch=abs(fine['conductivity_SI']/fine['native_conductivity_SI']-1)
    numerical=max(differences.values())<1e-6 and max(q['quadrature_max_absolute'] for q in controls)<1e-6
    physical_compatibility=bool(np.max(mismatch)<.2)
    np.savez_compressed(OUT/'kinetic-response.npz',**fine,orbital_relative_responses=responses,coordinate_frequencies=frequencies)
    result=dict(classification='Counterexample candidate',cohort_cells=len(indices),
        cohort_mass_fraction=float(d['dm'][indices].sum()/d['dm'].sum()),
        temperature_to_Fermi_range=[float(data['theta'][indices].min()),float(data['theta'][indices].max())],
        numerical_gates_passed=bool(numerical),quadrature_relative=differences,frequency_controls=controls,
        physical_conductivity_compatibility_gate_passed=physical_compatibility,
        native_conductivity_ratio_range=[float(np.min(fine['conductivity_SI']/fine['native_conductivity_SI'])),float(np.max(fine['conductivity_SI']/fine['native_conductivity_SI']))],
        tau_eff_seconds_range=[float(fine['tau_eff'].min()),float(fine['tau_eff'].max())],
        strong_ion_scattering_share_range=[float(fine['strong_ion_scattering_share'].min()),float(fine['strong_ion_scattering_share'].max())],
        orbital_memory_max_relative_by_harmonic=np.max(abs(responses),axis=1).tolist(),
        whole_star_closure_replaced=False,physical_model_error_certified=False,full_dynamic_charge_solved=False,
        missing='Electron-electron collisions, ion correlations and non-Born corrections, partially ionized/nondegenerate envelope and photon transfer; static compatibility is not a physical error bound.',seconds=time.monotonic()-start)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
