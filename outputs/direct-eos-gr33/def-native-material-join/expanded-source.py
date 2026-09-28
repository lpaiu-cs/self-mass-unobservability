def source(model,steps):
    z=np.load(OUT/f'coupled-{steps}.npz');assert int(z['completed_steps'])==steps
    keys=['U','I','bulk_I','u','theta','eta','h','j','Pi','t']
    data={k:z['snapshot_'+k] for k in keys};t=np.arange(17)*run.old.END/16
    keep=np.array([np.argmin(abs(data['t']-tt)) for tt in t]);assert np.max(abs(data['t'][keep]-t))<1e-18
    data={k:v[keep] for k,v in data.items()};f=model.flow;m=model.m;b=model.bulk
    volumes=np.r_[b.volume,4*np.pi*m.RJ**2*m.vol];edges=np.r_[b.d['edges'][:-1],m.rf];radius=np.r_[b.d['r'],m.r]
    weight,delay,a,B,re=prior.green.Geometry(m)(radius-m.RJ)
    f.eos.y=f.eos.y0;ip=f.eos(f.initial[0],f.initial_temperature)[0];tau0=(f.initial[2]-(m.a-m.a0)*f.eos.cx*f.initial[0])/m.a
    bp0=b.eos.base.gas(np.zeros(b.n),np.zeros(b.n))[0]
    trace=[];stress=[];baryon=[];gasE=[];press=[];phE=[];phP=[];velocity=[];lum=[];native_u=[]
    for k,state in enumerate(data['U']):
        model.Pi=data['Pi'][k];model.h=data['h'][k];model.j=data['j'][k];dm=-np.diff(model.h);model.mass=model.mass0-model.h[1:]+model.h[:-1];model.set_material(t[k])
        rho,v,lt,y=f.primitive(state);p,*_=f.eos(rho,lt);tau=(state[2]-(m.a-m.a0)*f.eos.cx*state[0])/m.a
        pp,uu,*_=b.eos.gas(data['theta'][k],data['eta'][k]);u=data['u'][k];native_u.append(float(np.max(abs(uu-u)/np.maximum(abs(u),1.))))
        beta=model.velocity();kinetic_trace=-model.mass*(model.cx*C*C+u)*beta**2/(1+np.sqrt(1-beta**2))
        de=dm*b.u0+model.mass*(u-b.u0);deep_nonrest=de+kinetic_trace;dpb=(pp-bp0)*b.volume
        unit=f.eos.rho0*C*C*volumes[b.n:];dp=p-ip
        trace.append(np.r_[deep_nonrest-3*dpb,((tau-tau0)-state[1]*v-3*dp)*unit])
        stress.append(np.r_[deep_nonrest-dpb,((tau-tau0)-state[1]*v-dp)*unit])
        gasE.append(np.r_[de+model.kinetic(),(tau-tau0)*unit]);press.append(np.r_[dpb,dp*unit])
        baryon.append(np.r_[dm,(state[0].astype(np.longdouble)-f.initial[0])*f.eos.rho0*volumes[b.n:]])
        velocity.append(np.r_[beta,v]);II=data['I'][k].sum(0);IB=data['bulk_I'][k];en=b.d['num']*b.d['Einf']
        phE.append(np.r_[np.einsum('iqf,q,f->i',IB,b.w,en)/b.d['a']**4,np.einsum('iqf,q,f->i',II,b.w,en)/m.a**4]*volumes)
        phP.append(np.r_[np.einsum('iqf,q,f->i',IB,b.w*b.mu2,en)/b.d['a']**4,np.einsum('iqf,q,f->i',II,b.w*b.mu2,en)/m.a**4]*volumes)
        lum.append(2*np.pi*C*model.area[-1]*(II[-1,-1]@(model.number*model.E)))
    f.seed=f.initial_temperature.copy();rho,v,lt,y=f.primitive(data['U'][-1]);p,*_=f.eos(rho,lt)
    check=((tau-tau0)-data['U'][-1,1]*v-3*(p-ip))*unit
    seed=float(abs((check-np.array(trace)[-1,b.n:])@m.a)/max(np.max(abs(np.array(trace)[:,b.n:]@m.a)),1.))
    assert seed<1e-8 and max(native_u)<1e-11
    d=dict(t=t,radius=radius,edges=edges,volume=volumes,weight=weight,delay=delay,a=a,B=B,re=re,
        nonrest_trace_erg=np.array(trace),nonrest_stress_erg=np.array(stress),baryon_g=np.array(baryon),
        gas_nonrest_energy_erg=np.array(gasE),pressure_volume_erg=np.array(press),photon_energy_erg=np.array(phE)-phE[0],
        photon_radial_pressure_erg=np.array(phP)-phP[0],velocity=np.array(velocity),luminosity_per_mu=np.array(lum),
        deep_cells=b.n,cx=f.eos.cx,M_cm=m.bg.M*m.R,K_cm=m.bg.K*m.R,inner_material_face_cm=m.rf[0],RJ=m.RJ)
    np.savez_compressed(OUT/f'source-{steps}.npz',**d)
    return d,dict(steps=steps,primitive_seed_relative=seed,saved_u_EOS_relative=max(native_u),maximum_primitive_residual=f.max_recovery)
