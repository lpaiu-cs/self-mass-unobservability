def conserved(self,rho,v,sigma,a):
    p,u,gamma,T,kap=self.eos(rho,sigma);rho,v,a,p,u,gamma=[np.asarray(x,np.longdouble) for x in (rho,v,a,p,u,gamma)];root=np.sqrt(1-v*v);W=1/root;wm=v*v/(root*(1+root));D=rho*W
    h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor);S=rho*h*W*W*v
    K=a*(self.eos.cx*D*wm+(rho*u+p)*W*W-p)+(a-self.base.a0)*self.eos.cx*D
    U=np.array([D,S,K]);F=np.array([D*v,S*v+p,(K+a*p)*v])
    cs=np.sqrt(gamma*p/np.maximum(rho*h,1e-100))
    assert max(cs)<1,'Causal sound speed'
    return U,F,(p,u,gamma,T,kap,cs)
