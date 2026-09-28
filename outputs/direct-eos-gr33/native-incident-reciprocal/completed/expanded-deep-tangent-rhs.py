def rhs(self,t,z,probe=1.):
    j,f,field,rates=self.fields(t);points=[(1-f,self.point(j)),(f,self.point(j+1))]
    size=max(np.max(abs(field[:4]))/1e-5,np.max(abs(field[4])/np.maximum(abs(self.ap),self.a/self.R))/1e-3,1e-200)
    for weight,p in points:
        if weight==0:continue
        q=p['Q'];active=p['active'];mass=np.maximum(q[0],1.)
        size=max(size,float(np.max(abs(z[0,active])/mass[active])/1e-5),float(np.max(abs(z[1,active])/(mass[active]*C*C)))/1e-7)
        thermal=z[2]-(self.a-self.model.m.a0)*self.model.cx*C*C*z[0]
        units=np.maximum(abs(q[2]-(self.a-self.model.m.a0)*self.model.cx*C*C*q[0]),1.)
        size=max(size,float(np.max(abs(thermal[active])/units[active]))/1e-5,float(np.max(abs(z[3,active])/np.maximum(abs(q[3,active]),1.)))/1e-5)
    eps=probe/size;self.min_probe=min(self.min_probe,eps);F=np.zeros((4,self.n+1));G=np.zeros((4,self.n));dt=np.inf;error=np.zeros((4,self.n))
    for k,(weight,p) in enumerate(points,j):
        if weight==0:continue
        a=self.raw(k,z,field,eps);b=self.raw(k,z,field,eps/2)
        fa=(a[0].astype(LD)-p['flux'].astype(LD))/eps;fb=(b[0].astype(LD)-p['flux'].astype(LD))/(eps/2)
        ga=(a[1].astype(LD)-p['gravity'].astype(LD))/eps;gb=(b[1].astype(LD)-p['gravity'].astype(LD))/(eps/2)
        df,dg=self.deep_tangent(k,z,field)
        fa[:,:self.nb]=df;fb[:,:self.nb]=df;ga[:self.nb]=dg;gb[:self.nb]=dg
        F+=weight*np.asarray(2*fb-fa,float);G[1]+=weight*np.asarray(2*gb-ga,float);dt=min(dt,a[2],b[2])
        comparison=np.asarray(np.sum(abs(fb-fa),axis=1)/np.maximum(np.sum(abs(fb),axis=1),1.),float)
        error-=weight*np.asarray(np.diff(fb-fa,axis=1),float);error[1]+=weight*np.asarray(gb-ga,float)
        G[2]-=weight*field[1]*(p['base_rate'][2]+self.rest*p['base_rate'][0])
        G[1]-=weight*(rates[0]+rates[1])*p['Q'][1]
        G[2]-=weight*self.a*(p['Pr']*(rates[0]+rates[1])+2*p['Pg']*rates[0])
    drive=(self.transfer[j+1]-self.transfer[j])/(self.t[j+1]-self.t[j]);G+=drive
    rate=-np.diff(F,axis=1)+G
    comparison=np.sum(abs(error),axis=1)/np.maximum(np.sum(abs(rate),axis=1),1.)
    self.probe_error=max(self.probe_error,float(max(comparison)));self.probe_rows.append(dict(t=float(t),comparison=comparison.tolist()))
    return rate,F[:,0]-F[:,-1]+np.sum(G,axis=1),dt
