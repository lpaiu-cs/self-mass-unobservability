def rhs(self,U,t):
    m=self.base;V=self.primitive(U);rho,v,sigma,y=V
    _,_,thermo=self.conserved(rho,v,sigma,y,m.a);p,u,gamma,T,kap,cs=thermo
    L,R=self.reconstruct(V,t)
    UL,FL,tl=self.conserved(*L,m.af);UR,FR,tr=self.conserved(*R,m.af)
    sl=np.minimum(0,np.minimum((L[1]-tl[-1])/(1-L[1]*tl[-1]),(R[1]-tr[-1])/(1-R[1]*tr[-1])))
    sr=np.maximum(0,np.maximum((L[1]+tl[-1])/(1+L[1]*tl[-1]),(R[1]+tr[-1])/(1+R[1]*tr[-1])))
    flux=np.divide(sr*FL-sl*FR+sl*sr*(UR-UL),sr-sl,out=np.zeros_like(FL),where=sr>sl)*m.af*m.area
    # A passive species must follow the same mass flux. Reconstructing rho
    # and y independently does not preserve positivity of their product.
    # ponytail: first-order donor y; retain only if the registered grid comparison passes.
    donor=np.where(flux[0]>=0,np.r_[self.incoming_y(y),y],np.r_[y,self.eos.y0])
    flux[3]=flux[0]*donor
    flux[:,0]=self.join_flux(R[:,0],y[0])
    rate=-C*np.diff(flux)/m.vol
    rho,v=[np.asarray(x,np.longdouble) for x in (rho,v)]
    W=1/np.sqrt(1-v*v);h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor);E=rho*h*W*W-p
    rate[1]+=C*np.diff(m.rf)/m.vol*(-(m.r/m.RJ)**2*E*m.ap+2*m.a*m.r/m.RJ**2*p)
    self.eos.y=y
    ledger=np.array([C*(flux[0,0]-flux[0,-1]),C*(flux[2,0]-flux[2,-1]),C*(flux[3,0]-flux[3,-1])])
    dt=.35*np.min(np.diff(m.rf)/(C*m.a/m.B*(abs(v)+cs+1e-100)))
    self.max_speed=max(self.max_speed,float(max(abs(v))))
    return rate,ledger,dt,V,kap
