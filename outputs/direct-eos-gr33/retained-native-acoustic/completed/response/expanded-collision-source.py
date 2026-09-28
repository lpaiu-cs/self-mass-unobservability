def point(self,k,half=False):
    m=self.m;p=self.p;mark=time.monotonic();t=m.t[k];assert p['t'][k]==t
    self.precision(False);m.model.flow.seed=np.asarray(m.model.flow.seed,float)
    target=m.point(k);point=m.point;m.point=lambda index:target
    try:c=m.local(t)
    finally:m.point=point
    x=p['photon_history_scaled_occupation'][k]/(m.scale*AMP);g=p['material_history'][k]/AMP
    linear,_,escape,bound=m.collision(c,x,g,False)
    linear=linear*m.scale*AMP;bound=bound*m.scale*AMP;escape=escape*AMP
    # Full actual free-material state, including its noncollisional E/H.
    z=m.motion[k];self.precision(True);a=self.coefficients(k,np.zeros_like(z),0.)
    sa,ea=m.scattering_matrix(a);bank=np.load(BEFORE/f'photons/bank-128/point-{k}.npz')
    owner=max(float(np.max(abs(a[q]-bank[q]))/max(np.max(abs(bank[q])),1e-300)) for q in ['emit','loss','sc','rho','beta','u','p'])
    I=m.I[k].astype(LD);delta=p['photon_history_scaled_occupation'][k].astype(LD);X=I/m.scale;dx=delta/m.scale
    weight=m.Eweight/m.scale;N=m.Nweight/m.scale;rows=[];full=None
    for factor in ([1.,.5] if half else [1.]):
        b=self.coefficients(k,z,factor);sb,eb=m.scattering_matrix(b)
        db=(b['emit']-a['emit'])-(b['loss']-a['loss'])*I-b['loss']*delta*factor
        scatter=((sb-sa)@X.ravel()).reshape(X.shape)*m.scale+(sb@(dx*factor).ravel()).reshape(X.shape)*m.scale
        de=np.einsum('knqf,nqf->kn',eb-ea,X)+np.einsum('knqf,nqf->kn',eb,dx*factor)
        residual=db+scatter-factor*linear;br=db-factor*bound;er=de-factor*escape
        norm=max(np.sum(abs(residual)*weight),1.)
        rounding=16*np.finfo(LD).eps*np.sum((abs(a['emit'])+abs(b['emit'])+(abs(a['loss'])+abs(b['loss']))*abs(I))*weight)
        if not np.any(z) and not np.any(delta):rounding=0.
        original=residual.copy();elo,ehi=LD(m.E[0]),LD(m.E[-1])
        for _ in range(2):
            number=np.sum((residual-br)*N,axis=(1,2))+er[0]
            residual[:,0,0]-=number*ehi/(ehi-elo)/N[:,0,0]
            residual[:,0,-1]+=number*elo/(ehi-elo)/N[:,0,-1]
        projection=float(np.sum(abs(residual-original)*weight)/norm)
        number=np.sum((residual-br)*N,axis=(1,2))+er[0]
        ns=np.maximum(np.sum((abs(residual)+abs(br))*N,axis=(1,2))+abs(er[0]),1.)
        negative=float(np.sum(np.maximum(-(I+factor*delta),0)*weight)/max(np.sum(abs(I)*weight),1.))
        rows.append(dict(factor=factor,residual_over_full_source=float(np.sum(abs(residual)*weight)/norm),
            rounding_over_full_source=float(rounding/norm),number_relative=float(np.max(abs(number)/ns)),
            projection_over_full_source=projection,negative_photon_energy_relative=negative))
        if factor==1.:full=dict(t=t,photon=residual,bound=br,escape=er)
    row=dict(classification='Counterexample candidate',steps=self.n,k=k,t=float(t),base_owner_relative=owner,rows=rows,
        seconds=time.monotonic()-mark,passed=bool(owner<1e-9 and all(r['rounding_over_full_source']<.002 and r['number_relative']<1e-10 and r['projection_over_full_source']<1e-8 and r['negative_photon_energy_relative']<1e-12 for r in rows)))
    np.savez_compressed(self.folder/f'point-{k}.npz',**full);write(self.folder/f'point-{k}.json',row)
    self.precision(False);print(json.dumps(row),flush=True);return row
