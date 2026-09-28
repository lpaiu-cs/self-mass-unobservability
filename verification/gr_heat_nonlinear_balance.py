"""Exact polar-areal balances and nonlinear constraint propagation, not a trajectory."""
import json,sys
import sympy as sp
import gr_heat_entropy_closure as closure

g=closure.g;OUT=g.OUT/'gr-heat-nonlinear-balance'


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def run():
    assert not OUT.exists();closure.verify();OUT.mkdir()
    files=[g.ROOT/'verification/gr_heat_nonlinear_balance.py',closure.OUT/'manifest.json',closure.OUT/'symbolic.json',
        g.OUT/'gr-heat-initial-constraints/symbolic.json',g.OUT/'gr-heat-initial-constraints/gourgoulhon-0703035.pdf',
        g.OUT/'gr-heat-primitive-inverse/symbolic.json']
    save('plan.json',dict(classification='Proven',checkpoint='36f75f07',bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in files},
        assumptions=['G=c=1, signature -+++; N>0,a>0, r>0; polar areal coordinates and zero shift',
            'Smooth matter, metric and heat fields; no distributional shock/interface assertion',
            'Baryon-frame isotropic pressure with one radial comoving heat Q; fixed composition',
            'The earlier quadratic heat-only entropy ansatz and positive K,T,tau; no newly calibrated transport coefficient'],
        source=dict(classification='Imported from prior work',url='https://arxiv.org/abs/gr-qc/0703035',
            scope='General 3+1 equations and conventions; the explicit spherical identities below are independently differentiated from the metric and tensor.'),
        target='Replace the initial-v=0 restriction by exact nonlinear matter/metric identities and identify the finite evolution equations still to be solved.',
        boundary='No numerical time trajectory, weak solution, well-posedness theorem, spatial error bound, physical EOS/transport calibration or exterior radiation matching.'))
    t,r,theta,phi=sp.symbols('t r theta phi',real=True);coords=[t,r,theta,phi]
    N=sp.Function('N')(t,r);a=sp.Function('a')(t,r);metric=sp.diag(-N*N,a*a,r*r,r*r*sp.sin(theta)**2);inverse=metric.inv()
    Gamma=[[[sp.simplify(sum(inverse[k,l]*(sp.diff(metric[l,j],coords[i])+sp.diff(metric[l,i],coords[j])-sp.diff(metric[i,j],coords[l])) for l in range(4))/2)
        for j in range(4)] for i in range(4)] for k in range(4)]
    Ricci=sp.zeros(4)
    for i in range(4):
        for j in range(4):
            Ricci[i,j]=sp.simplify(sum(sp.diff(Gamma[k][i][j],coords[k])-sp.diff(Gamma[k][i][k],coords[j])+
                sum(Gamma[k][i][j]*Gamma[l][k][l]-Gamma[l][i][k]*Gamma[k][j][l] for l in range(4)) for k in range(4)))
    scalar=sp.simplify(sp.trace(inverse*Ricci));Einstein=sp.simplify(Ricci-metric*scalar/2);mass=r*(1-a**-2)/2
    assert sp.simplify(Einstein[0,0]/N**2-2*sp.diff(mass,r)/r**2)==0
    assert sp.simplify(Einstein[0,1]-2*sp.diff(a,t)/(a*r))==0
    assert sp.simplify(Einstein[1,1]/a**2+2*mass/r**3-2*sp.diff(N,r)/(N*a*a*r))==0
    E,S,R,P,D,v=[sp.Function(name)(t,r) for name in ['E','S','R','P','D','v']]
    tensor=sp.Matrix([[E/N**2,S/(N*a),0,0],[S/(N*a),R/a**2,0,0],[0,0,P/r**2,0],[0,0,0,P/(r*r*sp.sin(theta)**2)]])
    root=N*a*r*r*sp.sin(theta)
    def divergence_covariant(tensor,nu):
        mixed=tensor*metric
        return sp.simplify(sum(sp.diff(root*mixed[mu,nu],coords[mu])/root for mu in range(4))-
            sum(sp.diff(metric[i,j],coords[nu])*tensor[i,j] for i in range(4) for j in range(4))/2)
    BE=sp.diff(a*E,t)+sp.diff(N*r*r*S,r)/r**2+sp.diff(a,t)*R+sp.diff(N,r)*S
    BS=sp.diff(a*S,t)+sp.diff(N*r*r*R,r)/r**2+sp.diff(a,t)*S+sp.diff(N,r)*E-2*N*P/r
    BD=sp.diff(a*D,t)+sp.diff(N*r*r*D*v,r)/r**2
    assert sp.simplify(divergence_covariant(tensor,0)+N*BE/a)==0
    assert sp.simplify(divergence_covariant(tensor,1)-BS/N)==0
    current=sp.Matrix([D/N,D*v/a,0,0]);divJ=sum(sp.diff(root*current[i],coords[i])/root for i in range(4))
    assert sp.simplify(divJ-BD/(N*a))==0
    # The angular Einstein equation follows when these radial equations hold:
    # the exact radial contracted Bianchi identity supplies the missing stress.
    assert divergence_covariant(inverse*Einstein*inverse,1)==0
    print('PASS nonlinear Einstein and matter balances from the four-metric',flush=True)
    adot=-4*sp.pi*r*N*a*a*S;mdot=-4*sp.pi*r*r*N*S/a
    Edot=(-sp.diff(N*r*r*S,r)/r**2-adot*(E+R)-sp.diff(N,r)*S)/a
    C=sp.diff(mass,r)-4*sp.pi*r*r*E
    residual=sp.diff(mdot,r)-4*sp.pi*r*r*Edot
    lapse= N*a*a*(mass/r**2+4*sp.pi*r*R)
    assert sp.simplify((residual-4*sp.pi*r*N*a*S*C).subs(sp.diff(N,r),lapse))==0
    assert sp.simplify(sp.diff(mass,t).subs(sp.diff(a,t),adot)-mdot)==0
    V,eps,press,Q=sp.symbols('V eps press Q',real=True);W2=1/(1-V*V);w=eps+press
    EE=w*W2-press+2*Q*W2*V;SS=w*W2*V+Q*W2*(1+V*V);RR=w*W2*V*V+press+2*Q*W2*V
    assert sp.cancel(SS-V*EE-press*V-Q)==0
    assert sp.cancel(EE-V*(SS+Q)-eps)==0
    assert sp.cancel(RR-V*SS-press-V*Q)==0
    # A moving material boundary has dr_s/dt=N*v/a. Its baryon flux cancels.
    assert sp.expand(-N*D*V+a*D*N*V/a)==0
    moving=mdot+N*V/a*4*sp.pi*r*r*E
    assert sp.cancel(moving.subs({S:SS,E:EE})+4*sp.pi*r*r*N*(press*V+Q)/a)==0
    # Project the actual covariant acceleration and heat derivative in the
    # comoving radial orthonormal frame, without setting the velocity to zero.
    z=sp.Function('zeta')(t,r);Temp=sp.Function('T')(t,r);Heat=sp.Function('Q')(t,r)
    u=sp.Matrix([sp.cosh(z)/N,sp.sinh(z)/a,0,0]);e=sp.Matrix([sp.sinh(z)/N,sp.cosh(z)/a,0,0]);ecov=metric*e
    def covariant_along_u(vector):
        return sp.Matrix([sum(u[i]*(sp.diff(vector[k],coords[i])+sum(Gamma[k][i][j]*vector[j] for j in range(4))) for i in range(4)) for k in range(4)])
    def along(vector,field):return sum(vector[i]*sp.diff(field,coords[i]) for i in range(4))
    simplify=lambda x:sp.simplify(sp.expand_trig(x))
    accel=(ecov.T*covariant_along_u(u))[0]
    expected=along(u,z)+sp.cosh(z)*sp.diff(N,r)/(N*a)+sp.sinh(z)*sp.diff(a,t)/(N*a)
    assert simplify(accel-expected)==0
    assert simplify((ecov.T*covariant_along_u(Heat*e))[0]-along(u,Heat))==0
    expansion=sum(sp.diff(root*u[i],coords[i])/root for i in range(4))
    expected_theta=along(e,z)+sp.cosh(z)*sp.diff(a,t)/(N*a)+sp.sinh(z)*(sp.diff(N,r)/(N*a)+2/(a*r))
    assert simplify(expansion-expected_theta)==0
    kn,kt,nrate,Trate=sp.symbols('k_n k_T Dlnn DlnT',real=True)
    assert sp.expand(-nrate-kn*nrate-(kt+2)*Trate+(kn+1)*nrate+(kt+2)*Trate)==0
    save('result.json',dict(classification='Proven',passed=True,
        Eulerian_fields='W=(1-v^2)^(-1/2), D=nW, E=(epsilon+P)W^2-P+2QW^2v, S=(epsilon+P)W^2v+QW^2(1+v^2), R=(epsilon+P)W^2v^2+P+2QW^2v.',
        baryon='(aD)_t+(N*r^2*D*v)_r/r^2=0.',
        energy='(aE)_t+(N*r^2*S)_r/r^2=-a_t*R-N_r*S.',
        momentum='(aS)_t+(N*r^2*R)_r/r^2=-a_t*S-N_r*E+2*N*P/r.',
        metric='m=r*(1-a^-2)/2; m_r=4*pi*r^2*E; a_t=-4*pi*r*N*a^2*S; m_t=-4*pi*r^2*N*S/a; N_r/N=a^2*(m/r^2+4*pi*r*R).',
        angular_equation='The exact radial contracted Bianchi identity is checked from the full metric. If the tt,tr,rr Einstein equations and radial matter conservation hold, the angular Einstein residual vanishes for r>0, by its nonzero 2/r coefficient and spherical symmetry.',
        nonlinear_constraint='C=m_r-4*pi*r^2*E satisfies C_t=4*pi*r*N*a*S*C under the mass-rate law, energy balance and polar lapse. Hence C(t,r)=C(0,r)*exp(integral_0^t 4*pi*r*N*a*S dt) for a smooth solution. C=0 propagates without setting v or Q to zero; a numerical constraint defect is not automatically damped.',
        material_worldtube='For dr_s/dt=N*v/a, baryon flux cancels exactly and dm(t,r_s)/dt=-4*pi*r_s^2*N/a*(P*v+Q). The algebraic identity S-vE=P*v+Q is checked at every |v|<1.',
        comoving_derivatives='D_u=cosh(zeta)/N*partial_t+sinh(zeta)/a*partial_r; D_s=sinh(zeta)/N*partial_t+cosh(zeta)/a*partial_r; v=tanh(zeta).',
        acceleration='a_s=D_u zeta+cosh(zeta)*N_r/(N*a)+sinh(zeta)*a_t/(N*a).',
        expansion='theta=D_s zeta+cosh(zeta)*a_t/(N*a)+sinh(zeta)*(N_r/(N*a)+2/(a*r)).',
        full_heat_law='tau*D_u Q+Q=-K*(D_s T+T*a_s)-(tau*Q/2)*(theta+D_u ln[tau/(K*T^2)]). The projected covariant derivative has radial component D_u Q; no extra spurious basis derivative is retained.',
        constant_tau='For fixed composition, constant proper tau and differentiable K(n,T), baryon conservation reduces the last parenthesis to -(1+dlnK/dlnn)*D_u lnn-(2+dlnK/dlnT)*D_u lnT. This reproduces the earlier kr,kt entropy-closure coefficients at finite velocity.',
        scope='Conditional smooth nonlinear equations. They identify the coupled matter, heat and metric problem; they do not solve a finite trajectory or prove nonlinear well-posedness. The moving-surface mass flux is not by itself a measured luminosity or a radiative exterior matching condition.',
        actual_time_trajectory_computed=False,physical_EOS_certified=False,physical_transport_calibrated=False,observational_closure=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and not r['actual_time_trajectory_computed'] and not r['observational_closure']
    print('PASS full nonlinear spherical balances, homogeneous Hamiltonian-constraint propagation and complete scalar heat law; finite GR trajectory remains open',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
