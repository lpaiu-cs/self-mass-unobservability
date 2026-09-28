"""Phase276: final-charge synthesis and J0337 observational closure (arithmetic and a symbolic check; no new physics run).

Reads the accepted phase 268-275 manifests (each hash-checked against paper/revision-manifest.json), checks the free-fall closed
form symbolically, normalizes the endpoint charge by eta phi_inf^2, and computes the J0337 lag bound for a linear, passive,
scalar-driven lag of the inner white dwarf's monopole charge with the Cassini-limited weak-field charge.
Usage: python phase276-closure.py <out json>
"""
import hashlib, json, math, sys
from pathlib import Path
import sympy as sp
root = Path('E:/lab/self-mass-unobservability/.claude/worktrees/eft-massive-objects-gravity-f20672'); O = root/'outputs/direct-eos-gr33'
master = json.loads((root/'paper/revision-manifest.json').read_text(encoding='utf-8'))['sha256']
NAMES = {'268': 'native-quad-refined-primary', '269': 'native-eos-sensitivity', '270-271': 'native-charge-mechanism',
         '272-273': 'native-structure-eft-boundary', '274-275': 'native-tidal-photosphere'}
m = {}
for ph, name in NAMES.items():
    p = O/f'{name}-manifest.json'; rel = p.relative_to(root).as_posix()
    assert hashlib.sha256(p.read_bytes()).hexdigest() == master[rel], rel
    m[ph] = json.loads(p.read_text(encoding='utf-8'))
q4, q_rich, order, change = (m['268']['rows']['64'][k] for k in ['q4', 'richardson', 'observed_order', 'change_2x_to_4x'])
ff = m['270-271']['freefall']; comp = m['270-271']['composition_4x_T']; st = m['272-273']; tp = m['274-275']
assert m['268']['converged_2x_at_2pct'] and m['269']['endpoint_negative_ploff'] and st['passed'] and tp['passed'] and tp['photosphere']['sign_kept_everywhere']

# symbolic check of the free-fall closed form (phase 271): d2xi/dt2 = -c^2 (a^2/B^2) d/dr[alpha(phi0) dphi], alpha = -4 phi
t, r, c, D, eta, Rd, a, B, s = sp.symbols('t r c D eta R_d a B s', positive=True)
phi0, x = sp.Function('phi0')(r), sp.Function('x')(r)
f = (4*s*(1 - s))**4; F1 = sp.integrate(f, (s, 0, s)); F2 = sp.integrate(F1, (s, 0, s)); p = (t + x/c)/D
xi = 4*c**2*(a**2/B**2)*eta*Rd*D**2*((sp.diff(phi0, r)/r - phi0/r**2)*F2.subs(s, p) + phi0*sp.diff(x, r)*F1.subs(s, p)/(r*c*D))
force = -c**2*(a**2/B**2)*sp.diff(-4*phi0*eta*Rd*f.subs(s, p)/r, r)
residual = sp.expand(sp.diff(xi, t, 2) - force)
assert residual == 0 and F1.subs(s, 0) == 0 and F2.subs(s, 0) == 0

# normalization (declared protocol: eta = 1e-30, phi_inf = 1e-3; force and source coupling both proportional to phi_inf)
eta_, phi_inf = 1e-30, 1e-3
norm = dict(chain_continuum=q_rich, freefall_continuum=ff['richardson_freefall']['limit'], per_eta_phi2=q_rich/(eta_*phi_inf**2))
ph = tp['photosphere']
obs = dict(kaplan_central=[q_rich*v for v in ph['kaplan_central_ratio']], one_sigma=[q_rich*v for v in ph['one_sigma_ratio_range']],
           two_sigma=[q_rich*v for v in ph['two_sigma_ratio_range']], per_eta_phi2_kaplan=[q_rich*v/(eta_*phi_inf**2) for v in ph['kaplan_central_ratio']])

# J0337 closure: Cassini gamma-1 = (2.1 +- 2.3)e-5 (Bertotti et al. 2003); scalar-tensor gamma-1 = -2 alpha0^2/(1+alpha0^2)
g_low = 2.1e-5 - 2*2.3e-5
alpha_max = math.sqrt(-g_low/(2 + g_low))
phi_inf_max = alpha_max/4
drive = 3.0760784580697786e-10  # pulsar field modulation at the white dwarf per unit pulsar charge (phase 274, e_in = 6.92e-4)
limit = 1.7e-9  # Paper B lag-responding SEP limit (tau_chi = 2 d, worst phase)
kappa_decl = tp['dissipation']['kappa_struct_max']; kappa_cassini = kappa_decl*(phi_inf_max/phi_inf)**2
closure = dict(cassini_alpha0_max_2sigma=alpha_max, phi_inf_max=phi_inf_max, drive_modulation_per_unit_alpha_p=drive, paper_b_lag_limit=limit,
               required_kappa_lag_times_alpha_p=limit/(alpha_max*drive),
               freefall_delta_bound_cassini=alpha_max*kappa_cassini*drive, freefall_kappa_at_cassini=kappa_cassini,
               full_beta_relaxing_delta_bound=alpha_max*4*drive,
               freefall_delta_bound_phase274=tp['dissipation']['lagged_monopole_delta_bound'])
res = dict(accepted=dict(q4=q4, richardson=q_rich, observed_order=order, change_2x_to_4x=change, eos_max_relative=m['269']['max_abs_relative_change'],
                         state_share=comp['state_share'], baryon_share=comp['baryon_state']/comp['total'], freefall_limits_relative=ff['limits_relative'],
                         sign_margin=st['background_structure']['sign_margin_Linf'], omega0=st['static_eft']['omega0'],
                         force_ratio=st['static_eft']['long_wavelength_force_ratio'], j0337_over_omega0_sq=st['static_eft']['orbital_drives']['1.629 d (J0337 inner)']['over_omega0_squared'],
                         k2=tp['tidal']['k2_apsidal'], tau_wave=tp['tidal']['drives']['2 n_in (non-rotating tide)']['tau_wave_thick_envelope'],
                         kappa_struct=kappa_decl, photosphere_tau=tp['dissipation']['charge_layer_gray_tau']),
           symbolic_freefall_residual=str(residual), normalization=norm, observed_gravity_charge=obs, closure=closure,
           manifests={ph_: f'outputs/direct-eos-gr33/{n}-manifest.json' for ph_, n in NAMES.items()})
json.dump(res, open(sys.argv[1], 'w', encoding='utf-8'), indent=1)
print(json.dumps(dict(richardson=q_rich, per_eta_phi2=norm['per_eta_phi2'], kaplan=obs['kaplan_central'], one_sigma=obs['one_sigma'], **closure), indent=1))
