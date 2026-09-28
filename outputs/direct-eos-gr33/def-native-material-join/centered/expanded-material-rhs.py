def material_rhs(self,u,theta,eta):
    b=self.bulk;d=b.d;v=self.velocity()*C;p=b.eos.gas(theta,eta)[0];dp=p-self.f0['p0']
    vl=np.r_[0.,v];vr=np.r_[v,0.];pl=np.r_[dp[0],dp];pr=np.r_[dp,dp[-1]]
    sound=np.sqrt(self.face_K/(self.cx*self.face_rho));Z=self.cx*self.face_rho*sound
    self.deep_dt=float(.35*np.min(np.diff(self.edge)*self.f0['B']/d['a']/(np.maximum(sound[:-1],sound[1:])+abs(v))))
    vf=(vl+vr)/2;ps=(pl+pr)/2
    vf[0]=0.;ps[0]=dp[0]
    mass=self.area_gas*self.face_rho*vf
    momentum=self.area_gas*(self.face_p+ps+self.cx*self.face_rho*vf*vf)
    donor_u=np.where(mass>=0,np.r_[u[0],u],np.r_[u,u[-1]])
    nh=d['thermo'][:,4];y=d['y0']*(1+eta)
    donor_y=np.where(mass>=0,np.r_[nh[0]*y[0],nh*y],np.r_[nh*y,nh[-1]*y[-1]])
    energy=mass*((self.f0['af']-self.m.a0)*self.cx*C*C+self.f0['af']*(donor_u+(self.face_p+ps)/self.face_rho))
    neutral=mass*donor_y
    mass[-1],momentum[-1],energy[-1],neutral[-1]=self.mflux
    de=(self.mass-self.mass0)/b.volume*self.cx*C*C+(self.mass*u-self.mass0*b.u0)/b.volume
    geometry=b.volume/self.f0['B']*(-de*self.f0['ap']+2*d['a']*dp/d['r'])
    # Background quadrature retains the saved nonzero gas/gravity force;
    # it does not remove the radiation-driven physical initial imbalance.
    force=-np.diff(momentum-self.base_momentum_flux)+b.volume*self.f0['initial_support']+geometry
    return [mass,force,-np.diff(energy),-np.diff(neutral)]
