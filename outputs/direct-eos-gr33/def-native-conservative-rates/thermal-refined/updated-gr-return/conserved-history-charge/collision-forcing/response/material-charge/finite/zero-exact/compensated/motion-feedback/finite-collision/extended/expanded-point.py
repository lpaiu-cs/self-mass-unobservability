def point(self,k):
    m=self.m;p=self.p;start=time.monotonic();t=m.t[k];assert p['t'][k]==t
    # Evaluate the tangent at exactly the state used by the accepted solve.
    target=m.point(k);point=m.point;m.point=lambda index:target
    try:c=run.base.Response.local(m,t)
    finally:m.point=point
    x=p['photon_history_scaled_occupation'][k]/(m.scale*AMP);g=p['material_history'][k]/AMP
    linear,_,escape,bound=m.collision(c,x,g,True)
    linear=linear*m.scale*AMP;bound=bound*m.scale*AMP;escape=escape*AMP
    z=m.motion[k].copy();z[[2,3]]=p['moments'][k,[1,2]]/AMP
    zero=np.zeros_like(z);a=self.coefficients(k,zero,0.);sa,ea=m.scattering_matrix(a)
    bank=dict(np.load(run.prior.OUT/f'bank-128/point-{k}.npz'))
    owner=max(float(np.max(abs(a[q]-bank[q]))/max(np.max(abs(bank[q])),1e-300)) for q in ['emit','loss','sc','rho','beta','u','p'])
    I=m.I[k];delta=p['photon_history_scaled_occupation'][k];X=I/m.scale;dx=delta/m.scale
    rows=[];full=[];floors=[]
    for factor in [1.,.5]:
        b=self.coefficients(k,z,factor);sb,eb=m.scattering_matrix(b)
        db=(b['emit'].astype(LD)-a['emit'])-(b['loss'].astype(LD)-a['loss'])*I-b['loss'].astype(LD)*delta*factor
        scatter=((sb-sa)@X.ravel()).reshape(X.shape)*m.scale+(sb@(dx*factor).ravel()).reshape(X.shape)*m.scale
        de=np.einsum('knqf,nqf->kn',eb-ea,X)+np.einsum('knqf,nqf->kn',eb,dx*factor)
        dp=db+scatter
        remainder=dp-factor*linear;br=db-factor*bound;er=de-factor*escape
        weight=m.Eweight/m.scale;N=m.Nweight/m.scale
        norm=max(np.sum(abs(factor*linear)*weight),1.)
        absolute=float(np.sum(abs(remainder)*weight));relative=absolute/norm
        rounding=16*np.finfo(LD).eps*np.sum((abs(a['emit'])+abs(b['emit'])+(abs(a['loss'])+abs(b['loss']))*abs(I))*weight)
        if not np.any(z) and not np.any(delta):rounding=0.  # Identical arguments have exactly zero difference.
        before=remainder.copy();elo,ehi=LD(m.E[0]),LD(m.E[-1])
        for _ in range(2):
            number=np.sum((remainder-br)*N,axis=(1,2))+er[0]
            remainder[:,0,0]-=number*ehi/(ehi-elo)/N[:,0,0]
            remainder[:,0,-1]+=number*elo/(ehi-elo)/N[:,0,-1]
        projection=float(np.sum(abs(remainder-before)*weight)/norm)
        number=np.sum((remainder-br)*N,axis=(1,2))+er[0]
        ns=np.maximum(np.sum((abs(remainder)+abs(br))*N,axis=(1,2))+abs(er[0]),1.)
        negative=float(np.sum(np.maximum(-(I+factor*delta),0)*weight)/max(np.sum(abs(I)*weight),1.))
        rows.append(dict(factor=factor,remainder_over_linear=float(relative),energy_weighted_L1=absolute,rounding_over_linear=float(rounding/norm),
            number_relative=float(np.max(abs(number)/ns)),projection_over_linear=projection,negative_photon_energy_relative=negative,negative_packet_count=int(np.sum(I+factor*delta<0))))
        full.append((remainder,br,er));floors.append(float(rounding))
    ratio=rows[1]['energy_weighted_L1']/max(rows[0]['energy_weighted_L1'],1.)
    row=dict(classification='Counterexample candidate',k=k,t=float(t),base_owner_relative=owner,rows=rows,half_over_full=ratio,seconds=time.monotonic()-start,
             passed=bool(owner<1e-9 and max(r['negative_photon_energy_relative'] for r in rows)<1e-12 and max(r['rounding_over_linear'] for r in rows)<.002 and max(r['number_relative'] for r in rows)<1e-10 and max(r['projection_over_linear'] for r in rows)<1e-8))
    np.savez_compressed(OUT/f'point-{k}.npz',t=t,photon=full[0][0],bound=full[0][1],escape=full[0][2],half_photon=full[1][0],
                        linear=linear,linear_bound=bound,linear_escape=escape)
    write(OUT/f'point-{k}.json',row);print(json.dumps(row),flush=True);return row
