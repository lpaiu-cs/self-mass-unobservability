def stages(m,t,h,x,g,lus):
    cs=[m.local(t+c*h) for c in C];ss=[m.source(t+c*h) for c in C];v=m.pack(x,g);dim=len(v)
    guides=[m.guide(t+c*h) for c in C];guess=np.array([m.pack(x,gg) for gg in guides]).ravel();audit=[]
    seed=getattr(m,'resume_seed',None)
    continuing=seed is not None and t==seed['time']
    limit=8
    if continuing:
        guess=seed['solution'].copy();guides=[m.unpack(row)[1].copy() for row in guess.reshape(2,dim)]
        audit=list(seed['equations'])
    for newton in range(limit):
        maps=[m.jacobian(t+c*h,gg) for c,gg in zip(C,guides)];Js=[r[0] for r in maps]
        affine=[r[1]-(J@gg.ravel()).reshape(m.n,4) for r,J,gg in zip(maps,Js,guides)]
        inv=[m.inverse_pair(c,h*A[j,j],Js[j]) for j,c in enumerate(cs)]
        def L(j,value):
            xx,gg=m.unpack(value);p,q,*_=m.collision(cs[j],xx,gg)
            stream=(m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)
            return m.pack(stream+p,q+(Js[j]@gg.ravel()).reshape(m.n,4))
        def mat(value):
            vv=value.reshape(2,dim);return (vv-h*(A@np.array([L(j,row) for j,row in enumerate(vv)]))).ravel()
        def pre(value):
            rows=[]
            for j,row in enumerate(value.reshape(2,dim)):
                xx,gg=m.unpack(row);xx=lus[j].solve(xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape);rows.append(inv[j](xx,gg))
            return np.array(rows).ravel()
        src=np.array([m.pack(s[0]/(m.scale*AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+a) for s,c,a in zip(ss,cs,affine)])
        rhs=(np.tile(v,(2,1))+h*(A@src)).ravel();rhs=precise_rhs(m,t,h,v,maps,guides,rhs);op=LinearOperator((2*dim,)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
        sol=solve(m,op,P,rhs,guess);pairs=[m.unpack(row) for row in sol.reshape(2,dim)];rates=[];details=[]
        for j,(xx,gg) in enumerate(pairs):
            p,q,e,b=m.collision(cs[j],xx,gg,True);native=m.native(t+C[j]*h,gg,details=True);native=precise_native(m,t+C[j]*h,gg,native)
            rates.append(m.pack((m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)+p+ss[j][0]/(m.scale*AMP),q+native[0]))
            details.append((p,q,e,native))
        defect=(sol.reshape(2,dim)-v-h*(A@np.array(rates))).ravel();defect=precise_defect(m,t,h,v,sol,defect);relative=float(np.linalg.norm(defect)/max(np.linalg.norm(rhs),1e-290))
        moments=physical_norm(m,defect)/scales(m,rhs,sol);audit.append(dict(relative=relative,moments=moments.astype(float).tolist()))
        if relative<1e-12 and max(moments)<1e-13:break
        if newton==limit-1:
            np.savez_compressed(OUT/'rejected-joint-stage.npz',time=t,step=h,initial=v,solution=sol,
                guides=guides,defect=defect,physical_scales=scales(m,rhs,sol),native_rates=[d[3][0] for d in details])
            write(OUT/'rejected-joint-stage.json',dict(classification='Counterexample candidate',time=float(t),step=float(h),equations=audit))
            raise AssertionError(('True native joint Radau equation',audit))
        guess=sol;guides=[gg.copy() for _,gg in pairs]
    m.newton_iterations.append(audit);result=[];transport=[]
    for j,((xx,gg),(p,q,e,native)) in enumerate(zip(pairs,details)):
        now=t+C[j]*h;nr,raw,discard,F,gravity=native
        probes=[m.native(now,gg,v) for v in [.5,2.]]
        den=np.maximum(np.sum(abs(nr)*m.units,axis=0,dtype=LD),LD('1e-290'))
        probe=[(np.sum(abs(v-nr)*m.units,axis=0,dtype=LD)/den).astype(float).tolist() for v in probes]
        assert np.max(probe)<.002,('Native constitutive control',now,probe)
        native_balance=np.sum(raw,axis=1,dtype=LD)-(F[:,0]-F[:,-1]+np.sum(gravity,axis=1,dtype=LD))
        balance=(abs(native_balance)/np.maximum(np.sum(abs(raw),axis=1,dtype=LD),LD('1e-290'))).astype(float).tolist()
        assert max(balance)<1e-8,('Native shared-face balance',balance)
        m.stage_log.append(dict(time=float(now),probe=probe,native_balance=balance))
        m.stage_t.append(now);m.stage_h.append(h*B[j]);m.stage_states.append(m.conserved(gg));m.stage_native.append(nr*m.units);m.stage_collision.append(q*m.units);m.stage_discard.append(discard)
        result.append((xx,gg,p,q+nr,e,ss[j][1]/AMP,ss[j][2]));transport.append(nr[:,:2])
    m.guide_g=result[-1][1].copy()
    return result,sum(b*r for b,r in zip(B,transport))
