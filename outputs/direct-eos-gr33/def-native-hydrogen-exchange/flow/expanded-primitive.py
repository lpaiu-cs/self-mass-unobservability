def primitive(self,U):
    m=self.base;D=np.maximum(U[0],0);active=D>=self.eos.floor
    tau=(U[2]-(m.a-m.a0)*self.eos.cx*D)/m.a
    all_y=np.divide(U[3],D,out=np.full_like(D,self.eos.y0),where=active);self.eos.y=all_y
    sigma=self.seed.copy();p=np.zeros_like(D);lower=np.full_like(D,np.log(100.));upper=np.full_like(D,np.log(1e7))
    # Start inside the existing native support; extend it only when the
    # energy equation, rather than the vacuum seed, requires a new state.
    lo,hi=temperature.parent.EOS.limits(self.eos,D);sigma=np.clip(sigma,lo,hi)
    for iteration in range(48):
        revision=self.eos.fan.calls
        # Settle momentum at this fixed temperature before recording an
        # energy-residual sign. A bracket built from the previous
        # temperature's pressure is not a bracket for one function.
        for _ in range(3):
            previous_pressure=p.copy()
            v=np.divide(U[1],self.eos.cx*D+tau+p,out=np.zeros_like(D),where=active)
            assert max(abs(v))<1,'Subluminal conservative inverse'
            root=np.sqrt(1-v*v);rho=D*root;W=1/root;wm=v*v/(root*(1+root))
            lo,hi=self.eos.limits(rho);sigma=np.clip(sigma,lo,hi)
            p,u,gamma,T,kap,cvT,entropy=self.eos.evaluate(rho,sigma)
            if np.max(np.where(active,abs(p-previous_pressure)/np.maximum(abs(p),1e-100),0))<1e-12:break
        if self.eos.fan.calls!=revision:lower[:]=np.log(100.);upper[:]=np.log(1e7)
        recovered=self.eos.cx*D*wm+(rho*u+p)*W*W-p
        error=recovered-tau
        errors=np.where(active,abs(error)/np.maximum(abs(tau),1e-100),0);relative=np.max(errors)
        if relative<2e-11:break
        derivative=np.maximum(rho*W*W*cvT,1e-100)
        correction=np.where(active,error/derivative,0)
        # sigma is the reconstructed log temperature in this producer.
        # Domain-bracketed iterates do not clamp the final physical state:
        # a root that cannot satisfy energy within the native bank fails.
        moving=active&(errors>=2e-11)
        lower=np.where(moving&(error<0),sigma,lower);upper=np.where(moving&(error>0),sigma,upper)
        proposed=sigma-np.clip(correction,-.5,.5)
        proposed=np.where((proposed<=lower)|(proposed>=upper),(lower+upper)/2,proposed)
        sigma=np.where(moving,proposed,sigma)
    else:
        # The rare unresolved cells use the independent scalar solve that
        # succeeded on the preserved failing state. No energy tolerance
        # is relaxed, and the same EOS/unknown is solved.
        for i in np.flatnonzero(active&(errors>=2e-11)):
            def residual(theta):
                self.eos.y=np.array([all_y[i]])
                pressure=0.
                for _ in range(4):
                    vv=U[1,i]/(self.eos.cx*D[i]+tau[i]+pressure);root=np.sqrt(1-vv*vv);rr=D[i]*root;ww=1/root;wm=vv*vv/(root*(1+root))
                    pp,uu,*_=self.eos(np.array([rr]),np.array([theta]));pressure=float(pp[0])
                return float((self.eos.cx*D[i]*wm+(rr*uu[0]+pressure)*ww*ww-pressure-tau[i])/tau[i])
            lowerT,upperT=temperature.parent.EOS.limits(self.eos,np.array([rho[i]]))
            sigma[i]=brentq(residual,float(lowerT[0])+1e-10,float(upperT[0])-1e-10,xtol=5e-15,rtol=1e-15);self.scalar_roots+=1
        self.eos.y=all_y
        for _ in range(4):
            v=np.divide(U[1],self.eos.cx*D+tau+p,out=np.zeros_like(D),where=active);root=np.sqrt(1-v*v);rho=D*root;W=1/root;wm=v*v/(root*(1+root));p,u,*_=self.eos(rho,sigma)
        error=self.eos.cx*D*wm+(rho*u+p)*W*W-p-tau
        relative=np.max(np.where(active,abs(error)/np.maximum(abs(tau),1e-100),0))
        if relative>=2e-11:
            np.savez_compressed(OUT/'primitive-failure.npz',U=U,rho=rho,sigma=sigma,error=error,tau=tau,relative=relative)
            raise AssertionError(('Conservative primitive root',relative))
    self.max_recovery=max(self.max_recovery,float(relative));self.seed=sigma.copy()
    self.minimum_sigma=min(self.minimum_sigma,float(min(sigma[active])));self.maximum_sigma=max(self.maximum_sigma,float(max(sigma[active])))
    self.eos.y=all_y
    return np.array([rho,v,sigma,all_y])
