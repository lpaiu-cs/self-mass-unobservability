def stages(m,t,h,x,g,lus):
    """Solve both collocation stages as one linear physical system."""
    cs=[];sources=[];ls=[];errors=[]
    for fraction in RK_C:
        now=t+fraction*h;c=m.local(now);s,l,e=m.source(now)
        cs.append(c);sources.append(s/(m.scale*AMP));ls.append(l/AMP);errors.append(e)
    # Canonical coefficient knots are step boundaries. The closing stage
    # takes the mechanical derivative from this same interval, as before.
    cs,sources,ls=fit(m,t,h,cs,sources,ls)
    mechanical=cs[0]['mechanical'].copy()
    for c in cs:c['mechanical']=mechanical
    def stream(xx):return (m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)
    def L(j,v):
        xx,gg=m.unpack(v);p,q,*_=m.collision(cs[j],xx,gg)
        return m.pack(stream(xx)+p,q)
    v=m.pack(x,g);dim=len(v)
    q=np.array([m.pack(s+c['q'],m.gas(c['q'],c['qb'],c['qe'])+mechanical) for s,c in zip(sources,cs)])
    rhs=np.tile(v,(2,1))+h*(RK_A@q)
    def mat(value):
        u=value.reshape(2,dim);rates=np.array([L(j,u[j]) for j in range(2)])
        return (u-h*(RK_A@rates)).ravel()
    inverses=[m.inverse(c,h*RK_A[j,j]) for j,c in enumerate(cs)]
    def pre(value):
        result=[]
        for j,vv in enumerate(value.reshape(2,dim)):
            xx,gg=m.unpack(vv);xx=lus[j].solve(xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)
            result.append(inverses[j](m.pack(xx,gg)))
        return np.array(result).ravel()
    def unpack(value):
        pairs=[m.unpack(vv) for vv in value.reshape(2,dim)]
        return tuple(np.concatenate([p[j] for p in pairs]) for j in range(2))
    moment_owner=SimpleNamespace(unpack=unpack,Nweight=np.tile(m.Nweight,(2,1,1)),
        Eweight=np.tile(m.Eweight,(2,1,1)),eu=np.tile(m.eu,2),nu=np.tile(m.nu,2))
    op=LinearOperator((2*dim,)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
    iterations=[]
    sol,info=prior.owner.moment_gmres(op,rhs.ravel(),owner=moment_owner,x0=np.tile(v,2),M=P,
        rtol=1e-14,atol=0.,restart=20,maxiter=5,callback=iterations.append,callback_type='pr_norm')
    residual=float(np.linalg.norm(mat(sol)-rhs.ravel())/max(np.linalg.norm(rhs),1e-290))
    m.max_residual=max(m.max_residual,residual);m.max_iterations=max(m.max_iterations,len(iterations))
    assert info==0 and residual<1e-12,('Coupled Radau stage',info,residual,len(iterations))
    result=[]
    for j,vv in enumerate(sol.reshape(2,dim)):
        xx,gg=m.unpack(vv);p,q,e,_=m.collision(cs[j],xx,gg,True)
        result.append((xx,gg,p,q,e,ls[j],errors[j]))
    return result,mechanical
