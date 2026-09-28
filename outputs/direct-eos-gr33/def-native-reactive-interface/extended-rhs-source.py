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
    rate=-C*np.diff(flux)/m.vol
    W=1/np.sqrt(1-v*v);h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor);E=rho*h*W*W-p
    rate[1]+=C*m.dx/m.vol*(-(m.r/m.RJ)**2*E*m.ap+2*m.a*m.r/m.RJ**2*p)
    area=(m.r/m.RJ)**2;Fbase=m.F0/(area*(m.a/m.a0)**2)
    mu=np.sqrt(np.maximum(0,1-(m.RJ/m.r)**2*(m.a/m.a0)**2));Eg=2*Fbase/(1+mu);Pg=2*Fbase*(1+mu+mu*mu)/(3*(1+mu))
    work=np.zeros(self.n)
    for _ in range(2):
        F=Fbase+np.r_[0,np.cumsum(work[:-1])]/(m.a*m.a*area*C)
        fcom=self.eos.rho0*rho*kap*W*W*((1+v*v)*F-v*(Eg+Pg))
        force=m.a*C*W*fcom;work=-m.vol*m.a*m.a*C*W*v*fcom
    rate[1]+=force;rate[2]-=work/m.vol
    self.max_optical=max(self.max_optical,float(np.sum(m.vol/area*self.eos.rho0*rho*kap)))
    ledger=np.array([C*(flux[0,0]-flux[0,-1]),C*(flux[2,0]-flux[2,-1])-sum(work),-sum(work)])
    dt=.35*m.dx/np.max(C*m.a/m.B*(abs(v)+cs+1e-100))
    self.eos.y=y
    reaction=self.eos.reactions(rho,sigma)
    dilution=(1-mu)/2
    absorbed=reaction[:,0]*dilution[:,None]
    emitted=reaction[:,1]+reaction[:,2]*dilution[:,None]
    net=absorbed-emitted
    species=-m.a*rho*net[:,0]
    Q=rho*self.eos.nH*net[:,1]/C**2
    # Declared bath approximation: evaluate photons in the local static
    # frame, use its hemisphere mean for absorption; isotropic spontaneous
    # emission. A velocity-frame error estimate is recorded, not hidden.
    force_energy=rho*self.eos.nH*(absorbed[:,1]-dilution*reaction[:,2,1])*(1+mu)/2/C**2
    momentum=m.a*W*(force_energy+v*Q)
    energy=m.a*m.a*W*(Q+v*force_energy)
    rate[1]+=momentum;rate[2]+=energy;rate[3]+=species
    transfer=float(np.sum(energy*m.vol));ledger[1]+=transfer;ledger[2]+=transfer
    ledger=np.r_[ledger,C*(flux[3,0]-flux[3,-1]),np.sum(species*m.vol)]
    stiffness=dilution*reaction[:,0,0]/y+(reaction[:,1,0]+dilution*reaction[:,2,0])/(1-y)
    dt=min(dt,.2/max(float(max(stiffness)),1e-100))
    for amount,change in [(U[3],rate[3]),(U[0]-U[3],rate[0]-rate[3])]:
        losses=(change<0)&(U[0]>=self.eos.floor)
        if np.any(losses):dt=min(dt,.8*float(min(amount[losses]/(-change[losses]))))
    optical=float(np.sum(m.B*m.dx*rho*self.eos.nH*absorbed[:,1]/(C**3*Fbase)))
    self.max_absorption=max(self.max_absorption,optical)
    self.max_speed=max(self.max_speed,float(max(abs(v))))
    self.min_y=min(self.min_y,float(min(y[rho>=self.eos.floor])))
    j=self.join
    local_photon=float(np.sum(energy[j:]*m.vol[j:])-np.sum(work[j:]))
    local=np.array([C*(flux[0,j]-flux[0,-1]),C*(flux[2,j]-flux[2,-1])+local_photon,
        local_photon,C*(flux[3,j]-flux[3,-1]),np.sum(species[j:]*m.vol[j:])])
    ledger=np.r_[ledger,local]
    self.maximum_inner_characteristic=max(self.maximum_inner_characteristic,float(max(C*m.a[:j]/m.B[:j]*(abs(v[:j])+cs[:j]))))
    return rate,ledger,dt
