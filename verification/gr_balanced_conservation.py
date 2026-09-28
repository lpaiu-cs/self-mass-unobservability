"""Covariant conservation and a conditional TOV-balanced radial discretization."""
import json,sys
import sympy as s
import direct_eos_gr as g

OUT=g.OUT/'gr-balanced-conservation'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_balanced_conservation.py',g.OUT/'gr-spherical-mass-balance/manifest.json',
        g.OUT/'gr-cell-average-identity/manifest.json',g.OUT/'gr-spatial-preflight/manifest.json']
    save('plan.json',dict(classification='Proven',checkpoint='f161e59',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        scope='Exact smooth spherical equations and conditional discrete identities, not a numerical stellar time trajectory. c=G=1, signature -+++, fixed areal radial cells, transverse stress P and radial normal stress S.'))
    t,r,theta,phi=s.symbols('t r theta phi',real=True);coords=[t,r,theta,phi]
    N,a,E,J,S,P,D,v=[s.Function(k)(t,r) for k in ['N','a','E','J','S','P','D','v']]
    metric=s.diag(-N**2,a**2,r**2,r**2*s.sin(theta)**2);inverse=metric.inv()
    tensor=s.Matrix([[E/N**2,J/(N*a),0,0],[J/(N*a),S/a**2,0,0],
        [0,0,P/r**2,0],[0,0,0,P/(r**2*s.sin(theta)**2)]])
    mixed=tensor*metric;volume=N*a*r**2*s.sin(theta)
    # The mixed-index divergence includes the metric-derivative source.
    def divergence(nu):
        return s.simplify(sum(s.diff(volume*mixed[mu,nu],coords[mu]) for mu in range(4))-
            volume*sum(tensor[i,j]*s.diff(metric[i,j],coords[nu]) for i in range(4) for j in range(4))/2)
    energy=s.diff(a*r**2*E,t)+s.diff(N*r**2*J,r)+r**2*J*s.diff(N,r)+r**2*S*s.diff(a,t)
    momentum=s.diff(a*r**2*J,t)+s.diff(N*r**2*S,r)+r**2*E*s.diff(N,r)-2*N*r*P+r**2*J*s.diff(a,t)
    assert s.simplify(divergence(0)/(-N*s.sin(theta))-energy)==0
    assert s.simplify(divergence(1)/(a*s.sin(theta))-momentum)==0
    number_current=s.Matrix([D/N,D*v/a,0,0])
    baryon=sum(s.diff(volume*number_current[i],coords[i]) for i in range(4))/s.sin(theta)
    assert s.simplify(baryon-s.diff(a*r**2*D,t)-s.diff(N*r**2*D*v,r))==0
    # Einstein constraints turn the energy projection into a conservative
    # coordinate-volume normal energy, equivalent to shell Misner-Sharp mass.
    at=-4*s.pi*r*N*a*a*J
    Nr=N*(4*s.pi*r*a*a*(E+S)-s.diff(a,r)/a)
    Et=s.solve(energy,s.diff(E,t))[0]
    mass_conservative=s.diff(r*r*E,t)+s.diff(r*r*N*J/a,r)
    assert s.simplify(mass_conservative.subs(s.diff(E,t),Et).subs(s.diff(a,t),at).subs(s.diff(N,r),Nr))==0
    N0,P0,E0=[s.Function(k)(r) for k in ['N0','P0','E0']]
    reference_flux=r*r*N0*P0
    reference_source=-r*r*E0*s.diff(N0,r)+2*r*N0*P0
    tov={s.diff(P0,r):-(E0+P0)*s.diff(N0,r)/N0}
    assert s.simplify((s.diff(reference_flux,r)-reference_source).subs(tov))==0
    subtracted=s.diff(a*r*r*J,t)+s.diff(r*r*(N*S-N0*P0),r)+r*r*(E*s.diff(N,r)-E0*s.diff(N0,r))-2*r*(N*P-N0*P0)+r*r*J*s.diff(a,t)
    assert s.simplify((subtracted-momentum).subs(tov))==0
    # Shared faces telescope exactly; independent quadrature cancellation is
    # not presumed. Nonzero heat is a perturbation and is never subtracted.
    f=s.symbols('f0:4');rates=[f[i]-f[i+1] for i in range(3)]
    assert s.simplify(sum(rates)-f[0]+f[-1])==0
    perturbations=s.symbols('df0 df1 ds')
    balanced=perturbations[0]-perturbations[1]+perturbations[2]
    assert balanced.subs(dict.fromkeys(perturbations,s.Integer(0)))==0
    assert balanced.subs(dict(zip(perturbations,[0,0,1])))==1
    save('symbolic.json',dict(classification='Proven',passed=True,
        baryon='D=rho_B*W: partial_t(a*r^2*D)+partial_r(N*r^2*D*v)=0.',
        normal_energy='partial_t(a*r^2*E)+partial_r(N*r^2*J)=-r^2*J*N_r-r^2*S*a_t.',
        radial_momentum='partial_t(a*r^2*J)+partial_r(N*r^2*S)=-r^2*E*N_r+2*N*r*P-r^2*J*a_t.',
        constrained_energy='With a_t=-4*pi*r*N*a^2*J and N_r/N+a_r/a=4*pi*r*a^2*(E+S), partial_t(r^2*E)+partial_r(r^2*N*J/a)=0. Thus a shell integral of 4*pi*r^2*E obeys a common-face conservative mass flux, rather than differences of nearly equal total masses.',
        reference='A fixed piecewise smooth TOV reference (N0,E0,P0) must obey P0_r=-(E0+P0)*N0_r/N0 and have shared face flux r^2*N0*P0. Its momentum source integrates to that flux difference; composition/entropy jumps require the actual continuous pressure/lapse interface and cannot be independently fitted on each side.',
        subtracted_momentum='partial_t(a*r^2*J)+partial_r[r^2*(N*S-N0*P0)]=-r^2*(E*N_r-E0*N0_r)+2*r*(N*P-N0*P0)-r^2*J*a_t.',
        conditional_discrete_balance='Integrate reference source using its exact shared face flux difference; quadrature only the perturbation source. At the exact represented zero-flux equilibrium all perturbations vanish algebraically, so the momentum update is zero regardless of the large background gradient. Store the reference and perturbations separately. Shared face baryon/mass fluxes telescope. This does not bound a reference reconstruction defect or the quadrature of nonzero perturbations.',
        no_physical_force_subtraction='A nonzero heat flux has J=Q at v=0 and generates metric/time terms; keep these terms and the full momentum/energy/heat coupling. Never zero the actual initial heat acceleration to obtain a numerical equilibrium.',
        primitive_recovery='The conserved integrals are subcell moments. Their recovery as a single uniform rho,T,v,Q state requires a separately justified closure and generally cannot preserve all original entropy/energy moments. The existing point samples are not automatically these cell averages.',
        full_GR_evolution=False,continuous_space_or_time_error_certified=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('PASS covariant baryon/energy/momentum and conditional TOV-balanced conservation',flush=True);verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS balanced conservation SHA; not a finite-time GR trajectory',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
