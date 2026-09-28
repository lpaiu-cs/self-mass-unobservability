"""Exact spherical matter/radiation balance with a time-dependent GR metric.

Proven: conditional tensor identities. No grey opacity or angular closure is
silently supplied; all physical inputs remain separate from these identities.
"""
import json, shutil, sys
import sympy as s
import gr_radiation_eos_split as split
import gr_radiative_boundary as opacity

g=split.g;OUT=g.OUT/'gr-radiation-metric-equations'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();split.verify();opacity.verify()
    shutil.copy2(g.ROOT/'outputs/radiation-weih2020.pdf',OUT/'weih2020.pdf')
    paths=[g.ROOT/'verification/gr_radiation_metric_equations.py',split.OUT/'manifest.json',
        opacity.OUT/'manifest.json',OUT/'weih2020.pdf']
    save('plan.json',dict(classification='Proven',checkpoint='fb7d9e3',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        imported_source=dict(classification='Imported from prior work',
            url='https://arxiv.org/abs/2003.13580',doi='10.1093/mnras/staa1297',
            equations='15-26: radiation stress tensor, Eulerian moments and radiation four-force. Derive spherical specialization independently from the metric and tensor divergence.'),
        convention='G=c=1, ds^2=-N(t,r)^2 dt^2+a(t,r)^2 dr^2+r^2 dOmega^2, N,a>0. E,J,S,Pperp are Eulerian energy density, radial energy flux, radial stress and tangential stress. Sources sE,sJ are orthonormal Eulerian four-force components. Opacity in the four-force has inverse proper-length units, not mass-specific opacity units.',
        checks=['Einstein tt/tr/rr from Christoffel/Ricci tensors',
            'Time-dependent energy/momentum/baryon conservation from covariant divergences',
            'Misner-Sharp mass constraint propagation and material-surface balance',
            'Exact radial Lorentz transformation of radiation and its four-force',
            'Moving LTE matter/photon recombination',
            'Tolman equilibrium and outgoing null radiation in fixed Schwarzschild background'],
        physical_boundary='No actual radiation moments, absorption, scattering, emissivity, angular closure, atmosphere or neutrino stress are inferred from an LTE EOS or Rosseland opacity. Fixed-background null radiation is a transport control, not a self-gravitating static spacetime. Actual stellar time evolution remains open.'))


def check_bindings():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():
        assert g.c.sha(g.ROOT/rel)==digest,rel


def geometry():
    t,r,theta,phi=s.symbols('t r theta phi',real=True);x=[t,r,theta,phi]
    N=s.Function('N')(t,r);a=s.Function('a')(t,r)
    metric=s.diag(-N*N,a*a,r*r,r*r*s.sin(theta)**2);inv=metric.inv()
    G=[[[s.simplify(sum(inv[i,k]*(s.diff(metric[k,j],x[l])+s.diff(metric[k,l],x[j])
        -s.diff(metric[j,l],x[k])) for k in range(4))/2) for l in range(4)] for j in range(4)] for i in range(4)]
    R=s.Matrix(4,4,lambda i,j:s.simplify(sum(s.diff(G[k][i][j],x[k])-s.diff(G[k][i][k],x[j])
        +sum(G[k][i][j]*G[l][k][l]-G[l][i][k]*G[k][j][l] for l in range(4)) for k in range(4))))
    curvature=s.simplify(s.trace(inv*R));Einstein=s.simplify(R-metric*curvature/2)
    m=r*(1-a**-2)/2
    assert s.simplify(Einstein[0,0]/N**2-2*s.diff(m,r)/r**2)==0
    assert s.simplify(Einstein[0,1]-2*s.diff(a,t)/(a*r))==0
    assert s.simplify(Einstein[1,1]/a**2-2*(s.diff(N,r)/(N*a**2*r)-m/r**3))==0
    E,J,S,P=[s.Function(k)(t,r) for k in ['E','J','S','Pperp']]
    T=s.Matrix([[E/N**2,J/(N*a),0,0],[J/(N*a),S/a**2,0,0],
        [0,0,P/r**2,0],[0,0,0,P/(r*r*s.sin(theta)**2)]])
    div=[]
    for nu in range(2):
        div.append(s.simplify(sum(s.diff(T[mu,nu],x[mu])+sum(
            G[mu][mu][lam]*T[lam,nu]+G[nu][mu][lam]*T[mu,lam] for lam in range(4)) for mu in range(4))))
    e_lhs=s.diff(a*r*r*E,t)+s.diff(N*r*r*J,r)+r*r*(s.diff(a,t)*S+J*s.diff(N,r))
    j_lhs=s.diff(a*a*r*r*J,t)+s.diff(N*a*r*r*S,r)-N*a*r*r*(S*s.diff(a,r)/a+2*P/r-E*s.diff(N,r)/N)
    assert s.simplify(e_lhs-N*a*r*r*N*div[0])==0
    assert s.simplify(j_lhs-N*a*a*r*r*a*div[1])==0
    D=s.Function('D')(t,r);v=s.Function('v')(t,r);baryon=[D/N,D*v/a,0,0]
    divergence=sum(s.diff(baryon[mu],x[mu])+sum(G[mu][mu][lam]*baryon[lam] for lam in range(4)) for mu in range(4))
    assert s.simplify(N*a*r*r*divergence-s.diff(a*r*r*D,t)-s.diff(N*r*r*D*v,r))==0
    # Independent m_rt=m_tr compatibility. Do not assume a_t=0 when J!=0.
    mass=s.symbols('m');at=-4*s.pi*r*a*a*N*J
    ar=a**3*(4*s.pi*r*E-mass/r**2);Nr=N*a*a*(mass/r**2+4*s.pi*r*S)
    commutator=-4*s.pi*r*r/a*(at*(E+S)+J*Nr+N*J*ar/a)
    assert s.simplify(commutator)==0
    save('geometry.json',dict(classification='Proven',passed=True,
        metric_constraints=['m_r=4*pi*r^2*E','a_t=-4*pi*r*a^2*N*J',
            'N_r/N=a^2*(m/r^2+4*pi*r*S)','m_t=-4*pi*N*r^2*J/a'],
        energy='d_t(a*r^2*E)+d_r(N*r^2*J)=-r^2*(a_t*S+J*N_r)+N*a*r^2*sE',
        momentum='d_t(a^2*r^2*J)+d_r(N*a*r^2*S)=N*a*r^2*(S*a_r/a+2*Pperp/r-E*N_r/N)+N*a^2*r^2*sJ',
        baryon='d_t(a*r^2*D)+d_r(N*r^2*D*v)=0; D=rho_B/sqrt(1-v^2)',
        constraint_propagation='The energy equation with zero TOTAL source and the three metric constraints imply m_rt=m_tr. Matter/radiation four-forces must cancel in the total stress tensor. An unrepresented volumetric neutrino loss does not satisfy this zero-total-source premise.',
        limitation='Smooth spherical tensor identities in this gauge, not a numerical evolution, nonlinear stability proof or EOS certificate.'))
    return dict(t=t,r=r,N=N,a=a,E=E,J=J,S=S,P=P,e_lhs=e_lhs,j_lhs=j_lhs)


def frames():
    v=s.symbols('v',real=True);W=1/s.sqrt(1-v*v)
    boost=W*s.Matrix([[1,v],[v,1]]);inverse=boost.subs(v,-v)
    E,F,P=s.symbols('E F P',real=True)
    lab=s.Matrix([[E,F],[F,P]]);comoving=s.simplify(inverse*lab*inverse.T)
    assert s.simplify(boost*comoving*boost.T-lab)==s.zeros(2)
    J=W**2*(E-2*v*F+v*v*P);H=W**2*((1+v*v)*F-v*(E+P))
    assert s.simplify(comoving[0,0]-J)==0 and s.simplify(comoving[0,1]-H)==0
    eta,ka,ks=s.symbols('eta ka ks',nonnegative=True);A=eta-ka*J;B=-(ka+ks)*H
    force=s.simplify(boost*s.Matrix([A,B]));u=s.Matrix([W,W*v]);metric=s.diag(-1,1)
    assert s.simplify((u.T*metric*force)[0]+A)==0
    assert s.simplify(force[0]-W*(A+v*B))==0
    assert s.simplify(force[1]-W*(B+v*A))==0
    # LTE reconstruction for arbitrary velocity, including momentum and stress.
    eg,pg,erad=s.symbols('eg pg erad',positive=True)
    matter=s.Matrix([[(eg+pg)*W**2-pg,(eg+pg)*W**2*v],
        [(eg+pg)*W**2*v,(eg+pg)*W**2*v*v+pg]])
    rad=s.simplify(boost*s.diag(erad,erad/3)*boost.T)
    total=matter.subs({eg:eg+erad,pg:pg+erad/3},simultaneous=True)
    assert s.simplify(matter+rad-total)==s.zeros(2)
    assert s.simplify(rad[0,0]-rad[1,1]-2*erad/3)==0
    # A material shell moves at dr/dt=N*v/a. General stress identity first.
    q=s.symbols('q',real=True);fluid=s.Matrix([[eg,q],[q,pg]])
    moving=s.simplify(boost*fluid*boost.T)
    assert s.simplify(moving[0,1]-v*moving[0,0]-q-v*pg)==0
    save('frames.json',dict(classification='Proven',passed=True,
        comoving='Jrad=W^2*(Erad-2*v*Frad+v^2*Prad); Hrad=W^2*((1+v^2)*Frad-v*(Erad+Prad))',
        source='A=eta-ka*Jrad; B=-(ka+ks)*Hrad; sE=W*(A+v*B); sJ=W*(B+v*A). Matter receives the exact negative four-force.',
        moving_scattering='Even purely elastic comoving scattering (A=0) transfers Eulerian energy W*v*B when the fluid moves. Dropping this term breaks four-force covariance.',
        LTE='Boosting diag(a_rad*T^4,a_rad*T^4/3) and adding the gas-only fluid tensor reproduces the original total LTE perfect-fluid tensor for every |v|<1, with both transverse pressures included.',
        material_mass='Along dr/dt=N*v/a, dm/dt=-4*pi*N*r^2*(J-v*E)/a. For a fluid with comoving heat flux q this is -4*pi*N*r^2*(q+v*p)/a; for separate matter/radiation use their summed E,J, not just a comoving photon luminosity.',
        physical_boundary='The algebra does not supply ka,ks,eta, an angular distribution or scalar/neutral-particle radiation.'))


def equilibria(parts):
    t,r=parts['t'],parts['r'];N=parts['N'];C,M=s.symbols('C M',positive=True)
    def residual(replacements):
        return [s.simplify(parts[k].subs(replacements,simultaneous=True).doit()) for k in ['e_lhs','j_lhs']]
    lapse=s.Function('n')(r);rad=C/lapse**4
    tolman={N:lapse,parts['a']:s.Function('a0')(r),parts['E']:rad,parts['J']:0,parts['S']:rad/3,parts['P']:rad/3}
    assert residual(tolman)==[0,0]
    lapse=s.sqrt(1-2*M/r);rad=C/(r*r*lapse*lapse)
    null={N:lapse,parts['a']:1/lapse,parts['E']:rad,parts['J']:rad,parts['S']:rad,parts['P']:0}
    assert residual(null)==[0,0]
    save('equilibria.json',dict(classification='Proven',passed=True,
        Tolman='In an arbitrary prescribed static spherical metric, isotropic radiation E=C*N^-4, J=0, S=Pperp=E/3 obeys both source-free moment equations. This is the exact Tolman temperature law.',
        null='In prescribed Schwarzschild geometry, outgoing E=J=S=C/(r^2*N^2), Pperp=0 obeys both source-free moment equations. Its redshifted luminosity 4*pi*r^2*N^2*J is constant.',
        boundary='These are analytic radiation transport controls. Adding the null stress while keeping the Schwarzschild metric exact would violate the Einstein equations; self-gravity has not been evolved.'))


def run():
    check_bindings();parts=geometry();frames();equilibria(parts)
    save('result.json',dict(classification='Proven',completed=True,
        tensor_identities_passed=True,mass_constraint_propagation_passed=True,
        moving_LTE_and_four_force_passed=True,analytic_transport_controls_passed=True,
        physical_absorption_scattering_identified=False,angular_closure_certified=False,
        stellar_radiation_evolved=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    check_bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS spherical GR radiation tensor/metric identities and controls; physical evolution remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
