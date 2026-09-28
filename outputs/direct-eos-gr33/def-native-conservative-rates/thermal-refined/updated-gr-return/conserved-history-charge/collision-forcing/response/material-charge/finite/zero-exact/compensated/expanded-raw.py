def raw(self,k,delta,field,eps,row=None):
    row=self.point(k) if row is None else row;m=self.model;b=m.bulk;f=m.flow;nb=self.nb;self.raw_calls+=1
    V,a=self.geometry(field,eps);Q=row['Q']+eps*delta*np.asarray(row['active'])[None]
    E=Q[2]+eps*field[1]*(Q[2]+self.rest*Q[0])
    h=row['h'].copy();h[1:]-=eps*np.cumsum(delta[0,:nb])
    state=[h,Q[1,:nb]/C,E[:nb],Q[3,:nb]]
    u,theta,eta=m.recover_material(state,self.t[k],row['theta'])
    U=Q[:,nb:]/(V[nb:]*f.eos.rho0*np.array([1,C*C,C*C,f.eos.nH])[:,None]);U[2]=E[nb:]/(V[nb:]*f.eos.rho0*C*C)
    f.seed=row['seed'].copy();flux,gravity,dt,prim=self.atmosphere(U,self.t[k])
    deep,dg,ddt=self.deep(u,theta,eta)
    factor=4*np.pi*m.m.RJ**2*f.eos.rho0*C*np.array([1,C*C,C*C,f.eos.nH])
    aflux=flux*factor[:,None];ag=gravity*V[nb:]*f.eos.rho0*C*C
    original_shared=abs((deep[1,-1]+C*m.base_momentum_flux[-1])-aflux[1,0])/max(abs(aflux[1,0]),1.)
    assert original_shared<1e-12
    aflux[1,0]=deep[1,-1]
    ag[0]+=C*m.base_momentum_flux[-1]
    shared=float(np.max(abs(deep[:,-1]-aflux[:,0])/np.maximum(abs(aflux[:,0]),1.)))
    assert shared<1e-12,('Shared face',shared)
    return np.c_[deep[:,:-1],aflux],np.r_[dg,ag],min(dt,ddt),dict(theta=theta,eta=eta,primitive=prim)
