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
    dummy=sp.Dummy('temperature')
    assert sp.simplify(dT-sp.Subs(sp.diff(F(dummy,Bt),dummy),dummy,T)-sp.Subs(sp.Derivative(F(T,b),b),b,Bt)*sp.diff(Bt,T))==0
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
