def stages(m,t,h,x,g,lus):
    """Solve both collocation stages as one linear physical system."""
    cs=[];sources=[];ls=[];errors=[]
    for fraction in RK_C:
        now=t+fraction*h;c=m.local(now);s,l,e=m.source(now)
        cs.append(c);sources.append(s/(m.scale*AMP));ls.append(l/AMP);errors.append(e)
    # Canonical coefficient knots are step boundaries. The closing stage
    # takes the mechanical derivative from this same interval, as before.
    for c in cs:c['mechanical'][:,0]=cs[0]['mechanical'][:,0]
    def stream(xx):return (m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)
    def L(j,v):
        xx,gg=m.unpack(v);p,q,*_=m.collision(cs[j],xx,gg)
        return m.pack(stream(xx)+p,q)
    v=m.pack(x,g);dim=len(v)
    q=np.array([m.pack(s+c['q'],m.gas(c['q'],c['qb'],c['qe'])+c['mechanical']) for s,c in zip(sources,cs)])
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
    native_equation=[]
    rx,rg=unpack(rhs.ravel())
    denom=np.maximum([np.sum(abs(rx)*moment_owner.Nweight)+np.sum(abs(rg[:,1])*moment_owner.nu),
                      np.sum(abs(rx)*moment_owner.Eweight)+np.sum(abs(rg[:,0])*moment_owner.eu)],1.)
    for correction in range(3):
        additions=[]
        for j,vv in enumerate(sol.reshape(2,dim)):
            xx,gg=m.unpack(vv);native=m.native_neutral(cs[j],gg)
            approximate=cs[j]['mechanical'][:,1]+m.neutral_linear(cs[j],gg)
            delta=np.zeros_like(gg);delta[:,1]=native-approximate
            additions.append(m.pack(np.zeros_like(xx),delta))
        defect=rhs.ravel()-mat(sol)+(h*(RK_A@np.array(additions))).ravel()
        relative=float(np.linalg.norm(defect)/max(np.linalg.norm(rhs),1e-290))
        dx,dg=unpack(defect)
        moments=np.array([np.sum(abs(dx)*moment_owner.Nweight)+np.sum(abs(dg[:,1])*moment_owner.nu),
                          np.sum(abs(dx)*moment_owner.Eweight)+np.sum(abs(dg[:,0])*moment_owner.eu)])/denom
        native_equation.append(dict(relative=relative,moments=moments.astype(float).tolist()))
        if relative<1e-12 and max(moments)<1e-13:break
        assert correction<2,('Native nonlinear Radau equation',native_equation)
        delta,info=prior.owner.moment_gmres(op,defect,owner=moment_owner,x0=np.zeros_like(sol),M=P,
            rtol=1e-14,atol=0.,restart=20,maxiter=5,callback=iterations.append,callback_type='pr_norm')
        assert info==0;sol+=delta
    m.native_equations.append(native_equation)
    result=[];transport=[]
    for j,vv in enumerate(sol.reshape(2,dim)):

        xx,gg=m.unpack(vv);p,q,e,_=m.collision(cs[j],xx,gg,True)
        mech=np.asarray(cs[j]['mechanical'],np.longdouble).copy();mech[:,1]=m.native_neutral(cs[j],gg)
        approximate=cs[j]['mechanical'][:,1]+m.neutral_linear(cs[j],gg)
        q[:,1]+=mech[:,1]-approximate
        m.check_neutral(cs[j],gg,mech[:,1],h*RK_B[j]);transport.append(mech)
        result.append((xx,gg,p,q,e,ls[j],errors[j]))
    return result,sum(b*v for b,v in zip(RK_B,transport))
