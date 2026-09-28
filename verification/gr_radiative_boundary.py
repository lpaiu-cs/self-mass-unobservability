"""Counterexample candidate: separate photon transport from effective opacity."""
import json,sys
import numpy as np
import sympy as sp
import gr_microphysics as micro
import opacity_tables as tables

g=micro.g;OUT=g.OUT/'gr-radiative-boundary'
URL='https://docs.mesastar.org/en/latest/kap/overview.html'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_radiative_boundary.py',g.ROOT/'verification/opacity_tables.py',
        g.OUT/'initial-state-17-4.npz',g.OUT/'gr-opacity/new-GR-captured.npz',
        g.OUT/'gr-opacity/evaluation.npz',g.OUT/'opacity/internal/plan.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='9cf504d',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        imported_definition=dict(classification='Imported from prior work',url=URL,
            statement='Radiative opacity is Rosseland mean; the effective opacity is the harmonic combination of radiation and conduction. Current documentation is used for this definition only, not for the version-specific numerical formula.'),
        native_radiative_value_relative_tolerance=1e-4,
        opacity_units='Declared table model evaluated at rho_B and its opacity per baryon gram, as used by the existing diffusion operator. A physical composition/mass normalization error is not certified.',
        exterior_optical_depth='Unknown; cumulative integrals start at the finite saved outer boundary, not at infinity.',
        criterion='Report Rosseland transport length / redshifted-temperature scale and energy-flux / (c*a*T^4) without claiming a physical diffusion or atmosphere certificate.',
        physical_opacity_certified=False,atmosphere_solved=False,full_GR_evolution=False))
    symbolic()


def symbolic():
    # Two positive frequency groups with equal normalized Rosseland and
    # Planck weights. This is an explicit inverse-data nonuniqueness example.
    pairs=[(sp.Rational(1),sp.Rational(1)),(sp.Rational(2,3),sp.Rational(2))]
    rosseland=[1/(sp.Rational(1,2)/a+sp.Rational(1,2)/b) for a,b in pairs]
    planck=[(a+b)/2 for a,b in pairs]
    assert rosseland==[1,1] and planck==[1,sp.Rational(4,3)]
    save('spectral-nonuniqueness.json',dict(classification='Proven',passed=True,
        theorem='A Rosseland mean alone does not determine the frequency-dependent opacity or even a second differently weighted mean. Two positive equal-weight frequency-group opacity pairs (1,1) and (2/3,2) both have Rosseland harmonic mean 1, but arithmetic absorption means 1 and 4/3. With LTE pure absorption and equal Planck group weights they imply different integrated emissivities.',
        scope='Mathematical nonuniqueness for the stated positive two-group model. It does not assert these weights/opacities describe this actual star. Physical frequency/absorption data or an explicit additional closure is necessary for atmosphere/spectral inference from a Rosseland-only input.'))
    print('PASS Rosseland-only spectral nonuniqueness control',flush=True)


def radiative(model,p):
    _,X,Z,r,t=p[:5];grid=model.zgrid
    if Z<=grid[0]:return model.hydrogen(0,X,r,t)
    if Z>=grid[-1]:return model.hydrogen(len(grid)-1,X,r,t)
    j=int(np.searchsorted(grid,Z,side='right'))-1;w=(Z-grid[j])/(grid[j+1]-grid[j])
    return (1-w)*model.hydrogen(j,X,r,t)+w*model.hydrogen(j+1,X,r,t)


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state,aux=micro.inputs();raw=dict(np.load(g.OUT/'gr-opacity/new-GR-captured.npz'))
    contract=json.loads((g.OUT/'opacity/internal/plan.json').read_text())['internal_columns']
    assert contract[8:]==['kap_rad','kap_rad_rho','kap_rad_T']
    model=tables.Opacity();parts=np.array([radiative(model,p) for p in raw['parameters']]);rad=10**parts[:,0]
    effective=dict(np.load(g.OUT/'gr-opacity/evaluation.npz'))['values'][:,0]
    error=abs(rad/raw['inner'][:,8]-1);assert error.max()<plan['native_radiative_value_relative_tolerance'],error.max()
    assert np.all(rad>0) and np.all(effective<=rad)
    dm=state['dm'];r=state['r_mid_m']*100;rho=np.exp(state['lnd']);T=np.exp(state['lnT'])
    f=1-2*state['m_mid_geom']/state['r_mid_m'];theta_log=state['lnT']+state['nu']
    transport_length=1/(rad*rho)
    gradient=np.gradient(theta_log,r,edge_order=2)*np.sqrt(f)
    knudsen=transport_length*abs(gradient)
    # Midpoint quadrature in baryon mass; no claim of an enclosed continuum
    # optical depth or of a solved omitted exterior.
    dtau=rad*dm/(4*np.pi*r*r);tau=np.cumsum(dtau)-dtau/2
    photon=g.s.baryon_face_diffusion(state,rad);total=g.s.baryon_face_diffusion(state,effective)
    assert np.all(abs(photon)<=abs(total)*(1+1e-14))
    Nface=np.exp(state['nu_faces'][1:-1]);area=4*np.pi*(state['radius_faces_m'][1:-1]*100)**2
    Tface=np.exp((state['lnT'][:-1]+state['lnT'][1:])/2)
    # c*a_rad*T^4 = 4*sigma*T^4. This compares a diffusive photon flux
    # with the LTE energy-density light-crossing scale; it is not a flux limiter.
    flux_ratio=abs(photon)/(Nface*Nface*area*4*5.670400e-5*Tface**4)
    rows=dict(classification='Counterexample candidate',completed=True,cells=len(dm),
        native_radiative_value_relative_error=float(error.max()),
        effective_to_radiative_opacity_range=[float((effective/rad).min()),float((effective/rad).max())],
        maximum_Rosseland_transport_Knudsen=float(knudsen.max()),
        mass_fraction_Knudsen_above_point_one=float(dm@(knudsen>.1)/dm.sum()),
        maximum_diffusive_photon_flux_to_LTE_c_energy_density=float(flux_ratio.max()),
        photon_faces_above_LTE_c_energy_density=int((flux_ratio>1).sum()),
        outermost_midpoint_optical_depth_from_saved_boundary=float(tau[0]),
        total_midpoint_Rosseland_optical_depth_increment=float(dtau.sum()),
        exterior_optical_depth_known=False,opacity_mass_normalization_physical_error_certified=False,
        radiative_absorption_emissivity_known=False,spectral_likelihood_connected=False,
        physical_opacity_certified=False,atmosphere_solved=False,full_GR_evolution=False)
    np.savez_compressed(OUT/'diagnostics.npz',radiative_log_opacity_and_derivatives=parts,
        radiative_opacity=rad,effective_opacity=effective,transport_length_cm=transport_length,
        Rosseland_transport_Knudsen=knudsen,optical_depth_from_saved_boundary=tau,
        interior_photon_Linf=photon,interior_total_Linf=total,photon_LTE_flux_ratio=flux_ratio)
    save('result.json',rows);print('RADIATIVE BOUNDARY',rows,flush=True);verify()


def verify():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    path=OUT/'manifest.json'
    if not path.exists():save('manifest.json',dict(classification='Counterexample candidate',
        sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    for rel,digest in json.loads(path.read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS RADIATIVE BOUNDARY SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
