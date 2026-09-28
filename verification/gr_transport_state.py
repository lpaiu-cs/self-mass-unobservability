"""Counterexample candidate: fresh GR transport and buoyancy diagnostics.

The opacity includes conduction. It must not be used as a photon mean free
path, and a local convective diagnostic is not a mixing/evolution solution.
"""
import json, sys
import numpy as np
import sympy as sp
import gr_microphysics as micro

g=micro.g;OUT=g.OUT/'gr-transport'
SOURCE='https://ntrs.nasa.gov/api/citations/19760017020/downloads/19760017020.pdf'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='156009b',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/gr_transport_state.py',g.ROOT/'verification/gr_microphysics.py',
            g.ROOT/'verification/conservative_star.py',
            g.OUT/'transport-sources/thorne-1976.pdf']},
        sources=[dict(url=SOURCE,author='K. S. Thorne',
            title='The relativistic equations of stellar structure and evolution',
            inspected='NASA report printed pages 6-7, PDF pages 9-10, equations 8-10',
            imported_scope='Local proper-frame mixing-length variables and homogeneous-composition convective equations. The energy-density buoyancy diagnostic below is derived separately; turbulent closure is not inferred.')],
        state='outputs/direct-eos-gr33/initial-state-17-4.npz',
        state_existed_at_prepare=(g.OUT/'initial-state-17-4.npz').exists(),
        criterion='Report all signed buoyancy diagnostics; no physical stability pass from finite radial differences.',
        thermal_identity_relative_tolerance=1e-4,algebraic_tolerance=1e-12,
        approximation='Static pressure-equilibrated local parcels, frozen entropy and nuclear composition; radiation/conduction diffusion on the actual GR baryon faces. Report first-order and three-point radial derivatives separately.',
        missing=['Convective composition and energy flux','Radiating atmosphere','Time-dependent fluid and metric','Physical EOS and opacity errors'],
        physical_EOS_certified=False,full_GR_evolution=False))
    symbolic()


def symbolic():
    x,y,B=sp.symbols('x y B',positive=True)
    entropy=B*(y**4-x**4)*(1/x-1/y)
    factored=B*(x-y)**2*(x+y)*(x*x+y*y)/(x*y)
    assert sp.simplify(entropy-factored)==0
    rho,P,cv,cr,ct=sp.symbols('rho P cv cr ct',positive=True)
    gamma=cr+P/rho*ct**2/cv;cp=cv+P/rho*ct**2/cr
    ad=P/rho*(ct/cr)/cp
    assert sp.simplify(ad-(gamma-cr)/(gamma*ct))==0
    w,gamma1,dp,de,gravity,xi=sp.symbols('w gamma1 dp de gravity xi',positive=True)
    A=de/w-dp/(gamma1*P);density_contrast=(w/(gamma1*P)*dp-de)*xi
    acceleration=-gravity*density_contrast/w
    assert sp.simplify(acceleration-gravity*A*xi)==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        nonlinear_face_entropy_production=str(factored),
        theorem='For any positive face conductance, including a state-dependent positive one, the closed-face frozen-geometry redshifted-energy exchange conserves its sum and produces nonnegative total entropy. Zero production requires equal redshifted temperatures on every connected face.',
        buoyancy='With local proper gravity g, pressure equilibrium, frozen parcel entropy/composition, the local acceleration is g*A*xi, A=(d epsilon/dell)/(epsilon+P)-(dP/dell)/(Gamma1*P). Hence the local frozen-background buoyancy frequency squared is -g*A.',
        assumptions='A smooth first-law-consistent EOS, epsilon including physical nuclear rest energy, positive enthalpy and sound response; negligible parcel pressure lag, heat exchange and metric perturbation. No radial global-mode theorem.',
        scope='Conditional algebraic identities, not a full GR evolution or a physical transport certificate.'))
    print('PASS nonlinear face entropy and parcel buoyancy identities',flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=micro.inputs();a=aux['eos']
    op=g.OUT/'gr-opacity';assert json.loads((op/'result.json').read_text())['passed']
    opacity=np.load(op/'evaluation.npz')['values'][:,0]
    r=state['r_mid_m']*100;m=state['m_mid_geom']*100;c=g.c.gr.C*100
    rho=np.exp(state['lnd']);T=np.exp(state['lnT']);P=a[:,1];u=a[:,2]
    cx=(state['X']/g.c.A)@g.c.W;energy=rho*(cx*c*c+u);w=energy+P
    assert np.all(np.diff(r)<0) and np.all(w>0)
    cpT=a[:,10]+P/rho*a[:,6]**2/a[:,5]
    gamma=a[:,5]+P/rho*a[:,6]**2/a[:,10]
    score=float(abs(gamma/a[:,4]-1).max())
    assert np.all(cpT>0) and score<plan['thermal_identity_relative_tolerance'],score
    grad_ad=P/rho*(a[:,6]/a[:,5])/cpT
    f=1-2*m/r;assert np.all(f>0)
    # P is cgs; G/c^4 converts it to inverse-square centimetres.
    gravity=c*c*(m+4*np.pi*r**3*P*(g.c.gr.G*1000)/c**4)/(r*r*np.sqrt(f))
    dp_exact=-w*gravity/c**2
    # These derivatives describe the saved discrete profile, not an enclosure.
    buoyancy=[];pressure_errors=[]
    for order in [1,2]:
        de=np.gradient(energy,r,edge_order=order)*np.sqrt(f)
        dp=np.gradient(P,r,edge_order=order)*np.sqrt(f)
        A=de/w-dp/(a[:,4]*P)
        buoyancy.append(-gravity*A)
        pressure_errors.append(abs(dp/dp_exact-1))
    buoyancy=np.array(buoyancy)
    # Edge orders differ only at endpoints; use wider interior neighbours as
    # an additional diagnostic, without calling it a convergence certificate.
    wide=np.full(len(r),np.nan)
    for i in range(2,len(r)-2):
        de=(energy[i+2]-energy[i-2])/(r[i+2]-r[i-2])*np.sqrt(f[i])
        dp=(P[i+2]-P[i-2])/(r[i+2]-r[i-2])*np.sqrt(f[i])
        wide[i]=-gravity[i]*(de/w[i]-dp/(a[i,4]*P[i]))
    flux=g.s.baryon_face_diffusion(state,opacity)
    theta=np.exp(state['lnT']+state['nu'])
    entropy=flux*(1/theta[:-1]-1/theta[1:])
    scale=np.maximum(abs(flux/theta[:-1])+abs(flux/theta[1:]),1e-300)
    assert np.all(entropy>=-64*np.finfo(float).eps*scale)
    divergence=np.r_[flux,0.]-np.r_[0.,flux]
    conservation=float(abs(np.sum(divergence,dtype=np.longdouble))/np.sum(abs(divergence),dtype=np.longdouble))
    assert conservation<plan['algebraic_tolerance']
    dm=state['dm'];agreed=(buoyancy[1]<0)&(wide<0)
    tolman={**state,'lnT':-state['nu']};assert np.all(g.s.baryon_face_diffusion(tolman,opacity)==0)
    registered=np.load(g.s.OLD/'registered-face-flux.npz')['canonical'][1:-1]
    mismatch=abs(flux-registered)/np.maximum(abs(registered),1e-10*g.c.fresh.LSUN)
    np.savez_compressed(OUT/'diagnostics.npz',opacity=opacity,heat_capacity_P_times_T=cpT,
        gradient_adiabatic=grad_ad,proper_gravity=gravity,proper_buoyancy_squared=buoyancy,
        wide_stencil_proper_buoyancy_squared=wide,interior_Linf=flux,
        pair_entropy_production=entropy,pressure_derivative_errors=np.array(pressure_errors),
        registered_flux_mismatch=mismatch)
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(r),
        state_sha256=g.c.sha(g.ROOT/plan['state']),thermodynamic_gamma_identity_score=score,
        face_energy_conservation_residual=conservation,negative_face_entropy_cells=int((entropy<0).sum()),
        three_point_negative_buoyancy_mass_fraction=float(dm@(buoyancy[1]<0)/dm.sum()),
        wide_and_three_point_negative_buoyancy_mass_fraction=float(dm@agreed/dm.sum()),
        proper_buoyancy_squared_range=[float(buoyancy.min()),float(buoyancy.max())],
        registered_flux_mismatch_quantiles=np.quantile(mismatch,[0,.5,.9,.99,1]).tolist(),
        radiation_optical_depth_inferred_from_effective_opacity=False,
        physical_stability_certified=False,convection_solved=False,full_GR_evolution=False))
    print('NEW GR TRANSPORT STATE',score,conservation,float(dm@agreed/dm.sum()),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
