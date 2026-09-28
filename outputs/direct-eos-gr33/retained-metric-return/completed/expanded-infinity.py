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
    # Measure the new causal-support method before all saved packets.
    costs=[]
    for age in [(1-(1-1/np.sqrt(2)))*m.T/128,m.T/2,m.T*(128-(1-1/np.sqrt(2)))/128]:
        mark=time.monotonic()
        for a,r in [(8,8),(4,8),(8,4)]:packet(age,a,r)
        costs.append(time.monotonic()-mark)
    old=m.kernel(8,8);age=m.T/2;mask=old['delay']<=age
    integrand=old['weights']*(exterior.G/C**3*old['mass']*np.maximum(age-old['delay'],0)+exterior.G/(2*C**4)*old['stress']*mask)/m.M
    reference=np.array([np.sum(integrand[old['node_bins']==j],dtype=LD) for j in range(4)])
    owner_error=float(np.max(abs(packet(age,8,8)[0]-reference))/max(np.max(abs(reference)),1e-300));assert owner_error<.002
    upper=2*319*max(costs)+5;eligible=upper<CAPS['infinity']-(time.monotonic()-start)
    write(dest/'packet-pilot.json',dict(classification='Counterexample candidate',point_seconds=costs,
        upper_remaining_seconds=upper,eligible=eligible,original_owner_relative=owner_error,original_quad_failure_preserved=True));assert eligible
    emission={};gamma=1-1/np.sqrt(2);aw=np.arange(1,8,2,dtype=LD)/32
    for n in [64,128]:
        p=np.load(PHOTON/f'steps-{n}-reference-128.npz');h=m.T/n;t=p['accepted_angular_times'];L=p['accepted_angular_luminosity']
        assert L.shape==(2*n,4) and np.max(abs(t-h*(np.arange(n)[:,None]+[gamma,1.]).ravel()))<1e-18
        packets=LD(h)*np.tile([1-gamma,gamma],n)[:,None]*L.astype(LD)
        cumulative=np.r_[0.,np.cumsum(packets@aw)];ids=2*np.rint(p['t']/h).astype(int)
        error=float(np.max(abs(cumulative[ids]-p['radial_ports'][:,1,1]))/max(np.sum(abs(packets)@aw),1.));assert error<1e-12
        emission[n]=(packets,error)
    prior=dict(np.load(BEFORE/'direct-motion/charge-128-a8-r8.npz'));alpha=-m.K/m.M;paths=[];rows=[]
    for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]:
        mark=time.monotonic();packets,error=emission[n];inc=[];h=m.T/n
        for count in np.arange(17)*(n//16):
            values=np.zeros((2,4),LD)
            for j in range(count):
                for stage,offset in enumerate([gamma,1.]):values+=packet((count-j-offset)*h,a,r)*packets[2*j+stage]
            inc.append(np.sum(values,axis=1,dtype=LD))
        inc=np.asarray(inc,float);background=np.load(BEFORE/f'direct-motion/charge-{n}-a{a}-r{r}.npz')
        d=dict(t=m.t,normalized_exterior_parts=np.c_[background['normalized_exterior'],inc[:,0]],arrived_parts_erg=np.c_[background['arrived_energy_erg'],inc[:,1]])
        d.update(normalized_exterior=d['normalized_exterior_parts'].sum(1),arrived_energy_erg=d['arrived_parts_erg'].sum(1))
        row=dict(steps=n,angular=a,radial=r,response_port_relative=error,seconds=time.monotonic()-mark)
        wave=np.load(GR/f'wave-{n}-g8.npz')
        assert np.max(abs(wave['t']-m.t))<1e-18
        eps0=exterior.G/C**4*d['arrived_parts_erg'][:,0]/m.M
        deps=exterior.G/C**4*d['arrived_parts_erg'][:,1]/m.M;eps=eps0+deps
        previous=background['normalized']
        delta=(wave['free_scalar']+d['normalized_exterior_parts'][:,1]+(alpha+previous)*deps)/(1-eps)
        compact=background['compact'].astype(LD)+wave['free_scalar'];scalar=compact+d['normalized_exterior']
        direct=(scalar+alpha*eps)/(1-eps);normalized=previous+delta
        identity=float(np.max(abs(normalized-direct))/max(np.max(abs(direct)),1e-300));assert identity<1e-12
        d.update(compact=compact,epsilon=eps,mass_term=alpha*eps,scalar=scalar,previous=previous,
            metric_return=delta,normalized=normalized);paths.append(d);row['normalization_identity']=identity;rows.append(row)
        np.savez_compressed(dest/f'charge-{n}-a{a}-r{r}.npz',**d)
    fine=paths[0];controls={}
    for key in ['normalized','metric_return','normalized_exterior_parts','arrived_parts_erg']:
        extract=lambda d:d[key][:,1] if key.endswith('_parts') or key=='arrived_parts_erg' else d[key]
        ref=extract(fine);norm=max(float(np.max(abs(ref))),1e-300)
        controls[key]={tag:float(np.max(abs(extract(d)-ref))/norm) for tag,d in zip(['angular','radial','time'],paths[1:])}
    previous_error=float(np.max(abs(fine['previous']-prior['normalized']))/max(np.max(abs(prior['normalized'])),1e-300));assert previous_error<1e-12
    q,s0,dc,de,e0,ep,alpha_s=sp.symbols('q s dc de e ep alpha')
    q0=(s0+alpha_s*e0)/(1-e0)
    assert sp.factor((s0+dc+de+alpha_s*(e0+ep))/(1-e0-ep)-q0-(dc+de+(alpha_s+q0)*ep)/(1-e0-ep))==0
    mu,lo,hi=sp.symbols('mu lo hi');assert sp.integrate(mu,(mu,lo,hi))==(hi**2-lo**2)/2
    passed=all(v['angular']<.002 and v['radial']<.002 and v['time']<.02 for v in controls.values())
    write(dest/'result.json',dict(classification='Counterexample candidate',passed=passed,controls=controls,rows=rows,
        endpoint_normalized=float(fine['normalized'][-1]),previous_normalized=float(fine['previous'][-1]),
        metric_return_normalized=float(fine['metric_return'][-1]),endpoint_compact=float(fine['compact'][-1]),
        endpoint_exterior=float(fine['normalized_exterior'][-1]),endpoint_mass_term=float(fine['mass_term'][-1]),
        previous_same_emission_relative=previous_error,normalization_increment_symbolic=True,
        packet_fronts_integrated=True,nonradial_log_angle_coordinate=True,arrival_angular_integral_exact=True,delay_inverse_seconds=inverse_error,
        computed_packet_kernels=len(kernels),seconds=time.monotonic()-start,
        same_fine_background_for_both_response_clocks=True,original_full_cadence_compact_preserved=True,
        actual_new_angular_emission_at_null_infinity=True,actual_emission_mass_normalization=True,
        continuous_emission_error_enclosed=False,conditional_interval_reused=False,coupled_fixed_point_verified=False,final_charge_solved=False,full_goal_complete=False))
    assert passed
