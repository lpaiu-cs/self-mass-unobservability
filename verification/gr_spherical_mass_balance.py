"""Exact nonlinear spherical GR mass balance and its closed-surface null test."""
import json,sys
import sympy as sp
import direct_eos_gr as g

OUT=g.OUT/'gr-spherical-mass-balance'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir()
    paths=[g.ROOT/'verification/gr_spherical_mass_balance.py',
        g.OUT/'gr-heat-initial-constraints/gourgoulhon-0703035.pdf']
    save('plan.json',dict(classification='Proven',checkpoint='ed0f495',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        convention='c=G=1, signature -+++, arbitrary time-dependent spherical polar-areal metric, no assumption that the material is at rest or that the heat current is zero.',
        task='Derive the radial/time Einstein constraints directly from the four-dimensional metric; evaluate mass balance on a moving material boundary. This is a conditional theorem, not a numerical evolution of the saved stellar model.'))
    t,r,theta,phi=sp.symbols('t r theta phi',real=True);coords=[t,r,theta,phi]
    N=sp.Function('N')(t,r);a=sp.Function('a')(t,r)
    metric=sp.diag(-N*N,a*a,r*r,r*r*sp.sin(theta)**2);inverse=metric.inv()
    Gamma=[[[sp.simplify(sum(inverse[k,l]*(sp.diff(metric[l,j],coords[i])+sp.diff(metric[l,i],coords[j])-sp.diff(metric[i,j],coords[l])) for l in range(4))/2)
        for j in range(4)] for i in range(4)] for k in range(4)]
    def Ricci(i,j):
        return sp.simplify(sum(sp.diff(Gamma[k][i][j],coords[k])-sp.diff(Gamma[k][i][k],coords[j])+
            sum(Gamma[k][i][j]*Gamma[l][k][l]-Gamma[l][i][k]*Gamma[k][j][l] for l in range(4)) for k in range(4)))
    diagonal=[Ricci(i,i) for i in range(4)];scalar=sp.simplify(sum(inverse[i,i]*diagonal[i] for i in range(4)))
    Gtt=sp.simplify(inverse[0,0]*diagonal[0]-scalar/2)
    Grr=sp.simplify(inverse[1,1]*diagonal[1]-scalar/2)
    Grt=sp.simplify(inverse[1,1]*Ricci(1,0))
    mass=r*(1-a**-2)/2
    assert sp.simplify(Gtt+2*sp.diff(mass,r)/r**2)==0
    assert sp.simplify(Grt-2*sp.diff(mass,t)/r**2)==0
    assert sp.simplify(Grr+2*mass/r**3-2*sp.diff(N,r)/(r*N*a*a))==0
    v,eps,P,Q=sp.symbols('v eps P Q',real=True)
    E=(eps+P*v*v+2*v*Q)/(1-v*v)
    J=((eps+P)*v+Q*(1+v*v))/(1-v*v)
    S=(P+eps*v*v+2*v*Q)/(1-v*v)
    assert sp.simplify(v*E-J+P*v+Q)==0
    material_surface_rate=sp.simplify(-4*sp.pi*r*r*N*J/a+(N*v/a)*4*sp.pi*r*r*E)
    assert sp.simplify(material_surface_rate+4*sp.pi*r*r*N*(P*v+Q)/a)==0
    assert sp.simplify(material_surface_rate.subs({P:0,Q:0}))==0
    # Full normal-frame energy projection propagates the Hamiltonian
    # constraint, now for arbitrary normal radial stress rather than P.
    En,Sn=sp.symbols('En Sn',real=True);Jn=sp.Function('Jn')(t,r)
    A=4*sp.pi*r*a*Jn;mdot=-4*sp.pi*r*r*N*Jn/a
    Edot=N*A*(En+Sn)-N*sp.diff(r*r*Jn,r)/(a*r*r)-2*Jn*sp.diff(N,r)/a
    relation=N*(4*sp.pi*r*a*a*(En+Sn)-sp.diff(a,r)/a)
    assert sp.simplify((sp.diff(mdot,r)-4*sp.pi*r*r*Edot).subs(sp.diff(N,r),relation))==0
    # Vacuum lapse can depend on time only through a removable time factor.
    M=sp.symbols('M',real=True);vacuum_lapse=sp.sqrt(1-2*M/r)
    assert sp.simplify(sp.diff(vacuum_lapse,r)/vacuum_lapse-(1-2*M/r)**-1*M/r**2)==0
    result=dict(classification='Proven',passed=True,
        Einstein_components='For ds^2=-N(t,r)^2 dt^2+a(t,r)^2 dr^2+r^2 dOmega^2 and m=r*(1-a^-2)/2: G^t_t=-2*m_r/r^2, G^r_t=2*m_t/r^2, G^r_r=-2*m/r^3+2*N_r/(r*N*a^2). These were computed directly from the 4-metric Ricci tensor.',
        nonlinear_constraints='T^t_t=-E, T^r_t=-N*J/a, T^r_r=S give m_r=4*pi*r^2*E, m_t=-4*pi*r^2*N*J/a, N_r/N=a^2*(m/r^2+4*pi*r*S). No static-metric or zero-velocity assumption.',
        moving_material_boundary='With material speed v relative to slice normals, R_dot=N*v/a. For isotropic material pressure P and comoving radial heat density Q, v*E-J=-P*v-Q exactly, so d m(t,R(t))/dt=-4*pi*R^2*N/a*(P*v+Q).',
        closed_surface_null='If the surface is material, P_surface=0, Q_surface=0, and there is no other stress/energy channel, the enclosed gravitational mass is exactly constant, including nonlinear internal motion, heating and redistribution consistent with the total stress-energy conservation law.',
        exterior='In a connected spherical vacuum exterior, m_r=m_t=0. Integrating the radial lapse equation gives N=f(t)*sqrt(1-2*M/r); absorb f into time. Thus the exterior monopole is Schwarzschild with constant M on a regular polar-areal patch. An internal relaxation state alone gives no changing exterior mass in this closed GR sector.',
        Hamiltonian_propagation='The full normal-energy projection with S, not just comoving P, and N_r/N+a_r/a=4*pi*r*a^2*(E+S) implies partial_r m_t=4*pi*r^2 E_t. This remains a continuum identity, not a discrete time/space error guarantee.',
        necessary_escape='A varying exterior mass/readout requires boundary flux or work, external coupling, non-spherical/multipolar structure, or additional fields/stresses outside these assumptions. A nonzero scalar charge or free-fall anomaly is not generated merely by adding GR heat relaxation.',
        saved_star_boundary_not_certified='The stored outermost EOS boundary has not been shown to be a free P=0 material surface joined to a vacuum atmosphere. Closing its numerical heat faces does not certify the physical surface assumptions of this theorem.',
        no_claims='No universal SEP result for arbitrary external fields, no physical EOS/transport calibration, no complete nonlinear observational inference, and no finite-time stellar trajectory are supplied by this theorem.')
    save('symbolic.json',result)
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('PASS nonlinear spherical GR mass balance and closed-surface null',flush=True);verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS nonlinear spherical GR mass-balance SHA',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
