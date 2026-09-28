def scattering(self,I,rho,v,kap,a):
    """Elastic Thomson in moving gas; positive photon packet remapping.

    The declared frequency box has absorbing exits. Record their number,
    energy and momentum separately; never deposit escaped photons as heat.
    """
    n=len(rho);q=self.q;m=self.freq;W=1/np.sqrt(1-v*v)
    edges=(self.bulk.edges_mu[None,:]-v[:,None])/(1-v[:,None]*self.bulk.edges_mu[None,:]);dw=np.diff(edges)/2
    p2=(edges[:,:-1]**2+edges[:,:-1]*edges[:,1:]+edges[:,1:]**2-1)/2
    num=I*self.w[None,:,None]*self.number[None,None,:]/a[:,None,None]**3
    rate=a[:,None]*C*rho[:,None]*kap[:,None]*W[:,None]*(1-v[:,None]*self.mu[None,:])
    loss=num*rate[:,:,None];net=-loss.copy();escape=np.zeros((n,3),dtype=I.dtype);rows=np.arange(n)[:,None]
    # Complete the positive interpolation at both ends. Dumping a whole
    # endpoint packet for an infinitesimal Doppler shift has a nonzero
    # escape rate as v tends to zero. Ghost packets leave this finite
    # frequency representation with their own number/energy/momentum.
    # ponytail: one extrapolated node per edge; spectral closure still needs a domain comparison.
    extended=np.r_[self.E[0]**2/self.E[1],self.E,self.E[-1]**2/self.E[-2]]
    for j,muin in enumerate(self.mu):
        for k,muout in enumerate(self.mu):
            probability=dw[:,k]*(1+.5*p2[:,j]*p2[:,k]);amount=loss[:,j,:]*probability[:,None]
            dest=self.E[None,:]*(1-v[:,None]*muin)/(1-v[:,None]*muout)
            high=np.searchsorted(extended,dest)
            assert high.min()>0 and high.max()<len(extended),'Scattering ghost support'
            low=high-1;upper=high;fraction=(dest-extended[low])/(extended[upper]-extended[low])
            for index,weight in [(low,1-fraction),(upper,fraction)]:
                outside=(index==0)|(index==m+1);values=amount*weight
                np.add.at(net[:,k,:],(rows,np.clip(index-1,0,m-1)),np.where(outside,0.,values))
                escaped=np.where(outside,values,0.);escape[:,0]+=escaped.sum(1);escape[:,1]+=(escaped*extended[index]).sum(1);escape[:,2]+=(escaped*extended[index]*muout).sum(1)
    error=np.max(abs(net.sum((1,2))+escape[:,0])/np.maximum(abs(loss).sum((1,2)),1e-250));self.scatter_number_error=max(self.scatter_number_error,float(error))
    return net*a[:,None,None]**3/(self.w[None,:,None]*self.number[None,None,:]),escape
