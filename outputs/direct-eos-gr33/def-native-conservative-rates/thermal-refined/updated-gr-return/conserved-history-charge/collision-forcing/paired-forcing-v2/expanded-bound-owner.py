def point(data,native,it):
    m,ds,ats,z,by=data;b=m.bulk;e=b.eos;f=m.flow;nb=b.n;g=m.m
    history.restore(m,z,it);theta=z['snapshot_theta'][it];eta=z['snapshot_eta'][it];I=z['snapshot_bulk_I'][it];beta=m.velocity()
    oldgas=e.gas(theta,eta);oldrad=e.radiation(theta,eta);base=e.base.radiation(theta+e.theta0,eta)
    nr=[];ne=[]
    for j in range(nb):
        r=by['deep',it,j];state=dict(r,fraction=np.array(r['fraction']))
        chi,emit=chem.prior.coefficients(native,state,b.d['Einf']/b.d['a'][j]);nr.append([chi+emit,emit]);ne.append(r['raw'][13]/AMU)
    nr=np.array(nr);ne=np.array(ne)-e.inventory[2]*e.xi
    df=e.spectral['frequency'][0]*(1-e.f)+e.spectral['frequency'][1]*e.f
    for channel in [0,1]:nr[:,channel]-=base[channel]*e.rinventory[channel]*e.xi[:,None]*(1+df[:,channel]*e.frequency_shift[:,None])
    assert nr.min()>=0 and ne.min()>0
    EW=b.photon_energy_weight;NW=EW/b.d['Einf'];before=m.collision(I,theta,eta,beta)
    gas,rad=e.gas,e.radiation
    def delta_gas(t,y):
        v=list(gas(t,y));v[6]=ne-oldgas[6];return v
    def delta_rad(t,y):
        v=list(rad(t,y));v[0]=nr[:,0]-oldrad[0];v[1]=nr[:,1]-oldrad[1];return v
    e.gas=delta_gas;e.radiation=delta_rad;saved_error=m.scatter_number_error
    try:delta=m.collision(I,theta,eta,beta)
    finally:e.gas=gas;e.radiation=rad;m.scatter_number_error=saved_error
    bound=b.factor[:,None,None]*(oldrad[1][:,None,:]-(oldrad[0]-oldrad[1])[:,None,:]*I)+before[3]
    dBound=b.factor[:,None,None]*((nr[:,1]-oldrad[1])[:,None,:]-((nr[:,0]-oldrad[0])-(nr[:,1]-oldrad[1]))[:,None,:]*I)+delta[3]
    scatter=b.scfactor[:,None,None]*(ne-oldgas[6])[:,None,None]*(np.einsum('qk,ikf->iqf',b.S,I)-I)
    expected=dBound-delta[3]+delta[1]+scatter
    deep_error=float(np.sum(abs(delta[0]-expected)*EW)/max(np.sum(abs(delta[0])*EW),1.));assert deep_error<1e-10
    old_mom=[np.sum(before[0]*EW,axis=(1,2))+b.volume*before[2][:,1],np.sum(bound*NW,axis=(1,2)),
             -(np.sum(before[0]*EW*b.mu[None,:,None],axis=(1,2))+b.volume*before[2][:,2])/b.d['a']]
    def_mom=[np.sum(delta[0]*EW,axis=(1,2))+b.volume*delta[2][:,1],np.sum(dBound*NW,axis=(1,2)),
             -(np.sum(delta[0]*EW*b.mu[None,:,None],axis=(1,2))+b.volume*delta[2][:,2])/b.d['a']]
    positive=[]
    for channel,field in [(0,I),(1,np.ones_like(I)),(1,I)]:
        w=field*EW*b.factor[:,None,None];positive.append(float(np.sum(abs(nr[:,channel]-oldrad[channel])[:,None,:]*w)/max(np.sum(nr[:,channel,None,:]*w),1.)))
    # Atmospheric full native primitives at exactly the saved D,S,K,H.
    d=ats[it];U=d['U'];V=np.array([d['rho'],d['v'],d['lt'],d['y']]);ids=d['active'];field=z['snapshot_I'][it].sum(0)
    f.eos.y=V[3];kap=f.eos(V[0],V[2])[4];newV=V.copy();newkap=kap.copy()
    roots=[by['atmosphere',it,int(j)] for j in ids]
    for j,r in zip(ids,roots):
        newV[:,j]=[r['rho'],r['v'],r['lt'],r['y']];newkap[j]=r['raw'][13]/r['raw'][0]*6.6524587321e-25/AMU
    fraction=np.array([r['fraction'] for r in roots]);affinity=np.array([r['affinity'] for r in roots]);coeff=m.spectrum.coefficients;cross=m.spectrum.cross
    def native_coeff(rho,lt,y,energy):
        assert len(rho)==len(ids);ab=np.zeros_like(energy);em=np.zeros_like(energy);T=np.exp(lt)[:,None,None]
        for level in range(10):
            sigma=cross(energy,level+1);pop=y*fraction[:,level];ab+=pop[:,None,None]*sigma
            logpop=np.full_like(pop,-np.inf);good=pop>0;logpop[good]=np.log(pop[good])
            em+=np.exp(logpop[:,None,None]+affinity[:,None,None]-energy/(chem.old.K*T))*sigma
        return ab,em
    def channels(prim,fn):
        rho,v,lt,y=prim[:,ids];D=(1-v[:,None]*m.mu[None,:])/np.sqrt(1-v*v)[:,None]
        ab,em=fn(rho*f.eos.rho0,lt,y,m.E[None,None,:]*D[:,:,None]/g.a[ids,None,None])
        factor=g.a[ids,None,None]*C*(rho*f.eos.rho0*f.eos.nH)[:,None,None]*D[:,:,None]
        return factor*ab,factor*em
    ab,em=channels(V,coeff);an,en=channels(newV,native_coeff)
    db=(en-em)*(1+field[ids])-(an-ab)*field[ids]
    vol=4*np.pi*g.RJ*g.RJ*g.vol;weight=m.energy_weight;number=weight/m.E
    oldscatter,oldesc=m.scattering(field[ids],V[0,ids]*f.eos.rho0,V[1,ids],kap[ids],g.a[ids])
    newscatter,newesc=m.scattering(field[ids],newV[0,ids]*f.eos.rho0,newV[1,ids],newkap[ids],g.a[ids])
    dp=np.zeros_like(field);dbfull=dp.copy();dep=np.zeros((m.n,3));oldfull=dp.copy();oldbound=dp.copy();oldexit=dep.copy()
    dp[ids]=db+newscatter-oldscatter;dbfull[ids]=db;dep[ids]=(newesc-oldesc)*vol[ids,None]
    oldbound[ids]=em*(1+field[ids])-ab*field[ids];oldfull[ids]=oldbound[ids]+oldscatter;oldexit[ids]=oldesc*vol[ids,None]
    def moments(p,bound,esc):return np.array([np.sum(p*weight,axis=(1,2))+esc[:,1],np.sum(bound*number,axis=(1,2)),
        -(np.sum(p*weight*m.mu[None,:,None],axis=(1,2))+esc[:,2])/g.a])
    am=moments(oldfull,oldbound,oldexit);dm=moments(dp,dbfull,dep)
    hydro=f.hydro
    def evaluate(prim,opacity,fn):
        f.hydro=lambda state,t:(np.zeros_like(state),np.zeros(3),1e100,prim,opacity);m.spectrum.coefficients=fn
        return m.local_rhs(U,field,ds[it]['t'])
    try:actual=evaluate(V,kap,coeff);changed=evaluate(newV,newkap,native_coeff)
    finally:f.hydro=hydro;m.spectrum.coefficients=coeff
    unit=vol*f.eos.rho0*C*C;checks=[]
    for rhs in [actual,changed]:
        ph=float(np.sum(rhs[1]*weight,dtype=np.longdouble));ga=float(np.sum(rhs[0][2]*unit,dtype=np.longdouble));esc=float(rhs[2][-1])
        checks.append(abs(ph+ga+esc)/max(abs(ph),abs(ga),1.))
    assert max(checks)<1e-10
    owner=np.array([-(changed[0][2]-actual[0][2])*unit,(changed[0][3]-actual[0][3])*vol*f.eos.rho0*f.eos.nH,
                    (changed[0][1]-actual[0][1])*unit])
    resolution=(np.sum(abs(owner-dm),axis=1)/np.maximum(np.sum(abs(dm),axis=1),1.)).tolist();assert max(resolution)<.002
    for before_c,after_c,packet in [(ab,an,field[ids]),(em,en,np.ones_like(field[ids])),(em,en,field[ids])]:
        w=weight[ids]*packet;positive.append(float(np.sum(abs(after_c-before_c)*w)/max(np.sum(after_c*w),1.)))
    photon=np.r_[delta[0],dp];bound=np.r_[dBound,dbfull];escape=np.r_[delta[2]*b.volume[:,None],dep].T
    original=np.c_[np.array(old_mom),am];change=np.c_[np.array(def_mom),dm]
    np.savez_compressed(OUT/f'point-{it}.npz',t=ds[it]['t'],photon=photon,bound=bound,escape=escape,
        original_moments=original,defect_moments=change)
    result=dict(classification='Counterexample candidate',it=it,t=ds[it]['t'],positive_relative=positive,
        deep_difference_owner=deep_error,atmosphere_full_owner_energy=checks,atmosphere_subtractive_resolution=resolution,
        actual_deep_velocity_max=float(np.max(abs(beta))),native_cells=nb+len(ids))
    write(OUT/f'point-{it}.json',result);return result
