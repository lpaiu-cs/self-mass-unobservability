def infinity():
    """Integrate each saved packet on its exact causal support, same4/8 rules."""
    import def_native_global_scalar_closure as exterior
    import sympy as sp
    from scipy.optimize import brentq
    start=time.monotonic();dest=OUT/'infinity';dest.mkdir(exist_ok=True);m=exterior.Exterior();LD=np.longdouble
    positive=m.model.bulk.edges_mu[m.model.bulk.edges_mu>=0];geometries={};kernels={};cuts={};inverse_error=0.
    def packet(age,a,r):
        nonlocal inverse_error
        if age<=0:return np.zeros((2,4))
        key=(age,a,r)
        if key in kernels:return kernels[key]
        if age not in cuts:
            cuts[age]=brentq(lambda mu:float(m.delay(np.array([mu]))[0].sol(1.)[0])-age,
                float(positive[-2]),1-1e-14,xtol=5e-15)
        mu0=cuts[age];gkey=(age,a)
        if gkey not in geometries:
            edges=np.sort(np.r_[positive,mu0]);gx,gw=np.polynomial.legendre.leggauss(a)
            angles=[];weights=[]
            for lo,hi in zip(edges[:-1],edges[1:]):
                if hi==1.:
                    angles.extend(lo+(hi-lo)*(gx+1)/2);weights.extend((hi-lo)*gw/2)
                else:
                    # Exact change of variable resolves the near-radial
                    # scale 1-mu; physical bins and Gauss orders stay fixed.
                    left,right=np.log1p(-hi),np.log1p(-lo);z=left+(right-left)*(gx+1)/2
                    angles.extend(-np.expm1(z));weights.extend(np.exp(z)*(right-left)*gw/2)
            mu=np.asarray(angles);mw=np.asarray(weights)
            sol,_=m.delay(mu);delay=sol.sol(1.);end=np.ones(len(mu))
            for j in np.flatnonzero(delay>age):
                end[j]=brentq(lambda s:float(sol.sol(s)[j])-age,0.,1.,xtol=5e-15,rtol=1e-14)
                inverse_error=max(inverse_error,abs(float(sol.sol(end[j])[j])-age))
            geometries[gkey]=(mu,mw,sol,end)
        mu,mw,sol,end=geometries[gkey];gx,gw=np.polynomial.legendre.leggauss(r)
        ss=(end[:,None]*(gx+1)/2).ravel();owners=np.repeat(np.arange(len(mu)),r)
        weights=(end[:,None]*gw/2*np.arccos(mu)[:,None]*mw[:,None]*mu[:,None]).ravel()
        mm=mu[owners];th=np.arccos(mm)*(1-ss);si=np.sin(th);co=np.cos(th);impact=m.r0*np.sqrt(1-mm*mm)
        N,b,cc=m.metric(si/np.sqrt(1-mm*mm));v=np.sqrt(1-si*si*(N/m.N0)**2)
        delay=sol.sol(ss)[owners,np.arange(len(ss))];assert np.max(delay-age)<1e-12
        mass=m.K*si*co/(cc*b*impact**2);stress=m.K*si*si*co*(N/m.N0)**2/(cc*cc*impact*v)
        bins=np.searchsorted(positive,mm,side='right')-1
        value=weights*(exterior.G/C**3*mass*(age-delay)+exterior.G/(2*C**4)*stress)/m.M
        charge=np.array([np.sum(value[bins==j],dtype=LD) for j in range(4)])
        arrival=np.maximum(positive[1:]**2-np.maximum(positive[:-1],mu0)**2,0)/2
        kernels[key]=np.asarray([charge,arrival],LD);return kernels[key]
    return packet,m,kernels,lambda:inverse_error
