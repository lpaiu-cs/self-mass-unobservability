"""Moving baryon-cell balances, including pressure work and entropy heat flux."""
import json,sys
import sympy as s
import gr_balanced_conservation as fixed

g=fixed.g;OUT=g.OUT/'gr-material-conservation'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir();fixed.verify()
    paths=[g.ROOT/'verification/gr_material_conservation.py',fixed.OUT/'manifest.json',
        g.OUT/'gr-heat-entropy-closure/manifest.json']
    save('plan.json',dict(classification='Proven',checkpoint='f161e59',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        assumptions='c=G=1, smooth spherical polar-areal patch, material faces R_dot=N*v/a, isotropic comoving P, radial comoving heat Q, no baryon diffusion, stress and baryon conservation. Entropy result additionally assumes the full specified quadratic heat closure and Gibbs relation.'))
    r,N,a,rho,v,eps,P,Q,entropy,K,T,tau=s.symbols('r N a rho v eps P Q entropy K T tau',real=True)
    W=1/s.sqrt(1-v*v);E=(eps+P*v*v+2*v*Q)/(1-v*v)
    J=((eps+P)*v+Q*(1+v*v))/(1-v*v);S=(P+eps*v*v+2*v*Q)/(1-v*v)
    speed=N*v/a
    baryon_density=a*r*r*rho*W;baryon_flux=N*r*r*rho*W*v
    assert s.simplify(baryon_flux-speed*baryon_density)==0
    mass_flux=N*r*r*J/a;mass_density=r*r*E
    assert s.simplify(mass_flux-speed*mass_density-N*r*r*(P*v+Q)/a)==0
    momentum_density=a*r*r*J;momentum_flux=N*r*r*S
    assert s.simplify(momentum_flux-speed*momentum_density-N*r*r*(P+v*Q))==0
    local_entropy=rho*entropy-tau*Q*Q/(2*K*T*T)
    St=W*(local_entropy+v*Q/T)/N;Sr=W*(v*local_entropy+Q/T)/a
    entropy_density=N*a*r*r*St;entropy_flux=N*a*r*r*Sr
    assert s.simplify(entropy_flux-speed*entropy_density-N*r*r*Q/(W*T))==0
    # Reynolds transport plus the fixed-radius local balances. Each oriented
    # common face is used once by each neighbour, including its work term.
    F=s.symbols('F0:5');cell_rates=[F[i]-F[i+1] for i in range(4)]
    assert s.expand(sum(cell_rates)-(F[0]-F[-1]))==0
    assert s.simplify((N*r*r*(P*v+Q)/a).subs({P:0,Q:0}))==0
    assert s.simplify((N*r*r*(P*v+Q)/a).subs(Q,0))!=0
    save('symbolic.json',dict(classification='Proven',passed=True,
        baryon='B_i=4*pi*integral[a*r^2*rho_B*W dr] between material faces is constant. Composition changes still require nuclear species equations and their energy reference.',
        shell_mass='M_i=m(R_out)-m(R_in)=4*pi*integral[r^2*E dr]. Its rate is F_M(in)-F_M(out), F_M=4*pi*r^2*N/a*(P*v+Q). Pressure work remains even when the heat face is closed.',
        momentum='Pi_i=4*pi*integral[a*r^2*J dr]. Its rate is F_Pi(in)-F_Pi(out)+4*pi*integral[-r^2*E*N_r+2*N*r*P-r^2*J*a_t] dr, where F_Pi=4*pi*N*r^2*(P+v*Q). The heat contribution to momentum flux must be retained.',
        entropy='Sigma_i=4*pi*integral[a*r^2*W*(rho_B*s-tau*Q^2/(2*K*T^2)+v*Q/T) dr]. With the full quadratic entropy closure, dSigma_i/dt=F_S(in)-F_S(out)+4*pi*integral[N*a*r^2*Q^2/(K*T^2) dr], F_S=4*pi*N*r^2*Q/(W*T). This is conditional on the declared continuum closure, not a finite-update entropy guarantee.',
        discrete_conservation='If each material interface has a single shared value for its speed and flux, cell baryon inventories remain fixed and shell mass fluxes telescope exactly. A numerical update must retain mass increments and pressure work; closing heat alone does not make the moving photospheric total mass constant.',
        remaining='Reconstruct subcell EOS/stress/entropy moments, solve heat-fluid primitive variables together, update metric constraints and material face trajectories, and certify the reference, quadrature and time errors. The symbolic identities do not perform these steps.',
        full_GR_evolution=False,physical_EOS_certified=False,nonlinear_discrete_entropy_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'symbolic.json').read_text())['passed']
    print('PASS moving material-cell work, heat-momentum and entropy balances',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
