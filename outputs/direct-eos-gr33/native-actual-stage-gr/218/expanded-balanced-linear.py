def solve(m,op,P,rhs,guess):
    iterations=[];options=dict(M=P,rtol=1e-14,atol=0.,restart=20,maxiter=5,callback=iterations.append,callback_type='pr_norm')
    sol=guess.copy();info=0
    for k in range(12):
        residual=rhs-op.matvec(sol);relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290))
        moments=physical_norm(m,residual)/scales(m,rhs,sol)
        if relative>=1e-14 and relative<1e-8 and max(moments)<1e-14:
            sol=polish(m,op,rhs,sol);residual=rhs-op.matvec(sol)
            relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290));moments=physical_norm(m,residual)/scales(m,rhs,sol)
        if relative<1e-14 and max(moments)<1e-13:
            radau.prior.owner.reuse.LINEAR.append(dict(initial_info=int(info),corrections=k,extended_residual=relative))
            m.max_residual=max(m.max_residual,relative);m.max_iterations=max(m.max_iterations,len(iterations));return sol
        assert k<11,('Four-moment linear residual',relative,moments.tolist())
        delta,code=gmres(op,np.asarray(residual,float),x0=None,**dict(options,rtol=1e-12));sol+=delta
        if k==0:info=code
