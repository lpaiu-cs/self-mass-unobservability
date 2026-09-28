def tangent(m,k,z,probe):
    row=m.point(k);zero=np.zeros_like(z);field=np.asarray(m.full_field,LD);nb=m.nb
    raw=m.raw(k,zero,field,0.);f=m.model.flow;b=m.model.bulk
    V=raw[3]['primitive'].astype(LD);join=f.join_state.astype(LD).copy()
    if not hasattr(m,'face_banks'):m.face_banks={}
    if k not in m.face_banks:m.face_banks[k]=dict(np.load(feedback.old.OUT/f'bank-{m.reference}/point-{k}.npz'))
    et=z.astype(LD).copy();et[2]-=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*C*C*et[0]
    p=primitive(m,k,et,field,m.face_banks[k])
    xi=m.model.mech.xi@np.r_[LD(0),-np.cumsum(z[0,:nb],dtype=LD)]
    p['dt'][:nb]+=b.eos.inventory[1]*xi/b.eos.gas(row['theta'],row['eta'])[2]
    dV=np.array([V[0]*p['dr'][nb:],p['dv'][nb:],p['dt'][nb:],V[3]*p['dy'][nb:]],dtype=LD)
    djoin=np.array([join[0]*p['dr'][nb-1],p['dv'][nb-1],p['dt'][nb-1],join[3]*p['dy'][nb-1]],dtype=LD)
    L,R,dL,dR=reconstruction(f,V,dV,join,djoin)
    baseline=face_flux(f,L,R)
    factor=4*np.pi*m.model.m.RJ**2*f.eos.rho0*C*np.array([1,C*C,C*C,f.eos.nH],dtype=LD)
    owner=float(np.max(np.sum(abs(baseline[:3]*factor[:3,None]-row['flux'][:3,nb:]),axis=1)/np.maximum(np.sum(abs(row['flux'][:3,nb:]),axis=1),1.)))
    assert owner<1e-8,('Reconstructed face owner',owner)
    m.face_owner_error=max(m.face_owner_error,owner)
    flux=flux_direction(m,k,L,R,dL,dR,probe)*factor[:,None]
    uf,ef=[np.interp(m.rEf,m.rE,np.asarray(v,float))[nb:] for v in field[:2]]
    flux[:3]+=(2*uf+ef)*baseline[:3]*factor[:3,None]
    # Differentiate the true donor. The amplified probe cannot select it.
    base=row['flux'][0,nb:];donor=np.where(base==0,flux[0]>=0,base>=0)
    yy=np.where(donor,np.r_[join[3],V[3]],np.r_[V[3],f.eos.y0])
    dy=np.where(donor,np.r_[djoin[3],dV[3]],np.r_[dV[3],LD(0)])
    flux[3]=f.eos.nH*(flux[0]*yy+base*dy)
    df,dg=m.deep_tangent(k,z,field)
    F=np.c_[df,flux]
    g=m.model.m;u,ell,lam,aden,ap=field[:,nb:];vol=3*u+lam
    E0=(row['Q'][2,nb:]+LD(m.rest)*row['Q'][0,nb:])/(m.a[nb:]*m.V[nb:])
    dE=(z[2,nb:].astype(LD)+LD(m.rest)*z[0,nb:])/(m.a[nb:]*m.V[nb:])-E0*vol
    P0=row['Pg'][nb:]/m.V[nb:]
    base=-g.r*g.r*g.ap*E0+2*g.a*g.r*P0
    variation=-g.r*g.r*(g.ap*dE+ap*E0+2*u*g.ap*E0)+2*g.a*g.r*(p['dp'][nb:]+(ell+u)*P0)
    gravity=4*np.pi*C*(np.diff(g.rf)*variation+np.diff(g.rf*uf)*base)
    G=np.r_[dg,gravity]
    nonzero=row['flux'][0]!=0
    ratio=float(AMP*np.max(abs(F[0,nonzero]/row['flux'][0,nonzero]),initial=0))
    m.physical_branch_ratio=max(m.physical_branch_ratio,ratio);assert ratio<.01
    return F,G,row['dt']
