    def rhs(self,U,t):
        m=self.base;V=self.primitive(U);rho,v,sigma=V
        _,_,thermo=self.conserved(rho,v,sigma,m.a);p,u,gamma,T,kap,cs=thermo
        L,R=self.reconstruct(V,t)
        UL,FL,tl=self.conserved(*L,m.af);UR,FR,tr=self.conserved(*R,m.af)
        sl=np.minimum(0,np.minimum((L[1]-tl[-1])/(1-L[1]*tl[-1]),(R[1]-tr[-1])/(1-R[1]*tr[-1])))
        sr=np.maximum(0,np.maximum((L[1]+tl[-1])/(1+L[1]*tl[-1]),(R[1]+tr[-1])/(1+R[1]*tr[-1])))
        flux=np.divide(sr*FL-sl*FR+sl*sr*(UR-UL),sr-sl,out=np.zeros_like(FL),where=sr>sl)*m.af*m.area
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
        return rate,ledger,dt
