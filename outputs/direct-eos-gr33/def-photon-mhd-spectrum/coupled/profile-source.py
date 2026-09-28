def integrate(b,profiles,order,tail_budget=1e-4):
    """Positive cell quadrature; discarded wings have a uniform operator bound."""
    n=int((b['u']<=60).sum());edges=b['edges_u'][:n+1];T=float(b['T'])
    h=2*np.pi*a.previous.old.base.hbar;scale=h/(a.previous.old.base.k*T)
    center=profiles[:,1]*scale;width=scale/profiles[:,5]
    strength=profiles[:,4]*profiles[:,7];damping=profiles[:,6]
    weights=np.cbrt(strength*damping*width**2);weights/=weights.sum()
    factor=float(b['H'])*(1+float(b['Ci'][:n].sum()/b['Cm']))
    bounds=tail_budget/factor*weights
    radius=np.maximum(16,2*np.sqrt(strength*damping/(np.sqrt(np.pi)*bounds)))
    tail=strength*(damping/np.sqrt(np.pi)/((radius/2)**2+damping**2)
        +erfc(radius/2)/(np.sqrt(np.pi)*damping))
    total=np.zeros(n);g,w=np.polynomial.legendre.leggauss(order)
    capacity0=15*float(b['arad'])*T**3/np.pi**4
    stimulated_scale=1.
    segments=0;begin=time.monotonic()
    for u0,du,U,ag,R in zip(center,width,strength,damping,radius):
        left=max(0.,u0-R*du);right=min(edges[-1],u0+R*du)
        if right<=left:continue
        first=max(0,np.searchsorted(edges,left,side='right')-1)
        last=min(n,np.searchsorted(edges,right))
        sd=du*max(1.,ag)
        ta,tb=np.arcsinh((np.array([left,right])-u0)/sd)
        # A unit interval in asinh coordinates resolves the core and wings
        # independently of the continuum source spacing and bin boundaries.
        knots=u0+sd*np.sinh(np.arange(np.ceil(ta),tb))
        mesh=np.unique(np.r_[left,edges[first+1:last],knots[(knots>left)&(knots<right)],right])
        t=np.arcsinh((mesh-u0)/sd);mid=(t[:-1]+t[1:])/2;half=np.diff(t)/2
        tq=mid[:,None]+half[:,None]*g
        uq=u0+sd*np.sinh(tq)
        x=(uq-u0)/du
        profile=wofz(x+1j*ag).real
        density=capacity0*uq**4*np.exp(-uq)/(-np.expm1(-uq))**2
        amount=(density*U*profile*(-np.expm1(-stimulated_scale*uq))
            *sd*np.cosh(tq)*half[:,None]*w).sum(axis=1)
        ids=np.searchsorted(edges,(mesh[:-1]+mesh[1:])/2,side='right')-1
        np.add.at(total,ids,amount);segments+=len(ids)
    return total,dict(seconds=time.monotonic()-begin,segments=segments,order=order,
        omitted_wing_response_bound=float(tail.sum()*factor),tail_budget=tail_budget,
        profiles=len(profiles),bound_scope='Exact Voigt model and exact positive bin integration; floating-point quadrature is checked separately. Streaming is skew and the common collision operator is dissipative.')
