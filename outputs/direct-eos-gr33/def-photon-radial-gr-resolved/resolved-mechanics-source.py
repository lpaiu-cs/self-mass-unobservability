def __init__(self,degree=4):
    patch.base.Background=Background
    self.problem=p=go.Problem(degree);self.m=m=p.model
    b,d,R,self.edges=geometry();self.R=R;self.tc=m.original.radiation.geometry.tc
    self.geo=go.task.fem.base.task.h.gr.G*.1*R**2/go.task.fem.base.task.h.gr.C**4
    self.gr=go.task.fem.base.task.h.gr
    inside=(m.points['r']>=self.edges[0])&(m.points['r']<self.edges[-1])
    r=m.points['r'][inside];dr=m.weights[inside];ids=(r>=self.edges[1]).astype(int)
    point=m.bg.sample(r);N,a=m.original.radiation.geometry.metric(r);A4=np.exp(-8*point['phi']**2)
    weight=4*np.pi*(100*R)**3*dr*N*a*A4*r*r
    self.volumes=np.bincount(ids,weights=weight,minlength=2)
    avg=sparse.csr_matrix((weight/self.volumes[ids],(ids,np.arange(len(r)))),shape=(2,len(r)))
    V,D=m.evaluation(r);mass,pressure,phi,v=[point[k] for k in ['m','p','phi','v']]
    bb=1-2*mass/r;alpha=-4*phi
    dlz=r*r*v*v/2-4*np.pi*r*r*A4*pressure/bb-mass/(r*bb)
    Dq=-sparse.diags(r)@D[0]-sparse.diags(3+dlz+3*alpha*r*v)@V[0]-sparse.diags(r*v+3*alpha)@V[1]
    self.Dq=(avg@Dq).tocsr()
    self.DS=(avg@(-sparse.diags(1/(r*bb))@self.sources(point)[2][:,:2])).toarray()
    self.gammaP=np.asarray(avg@(point['gamma']*point['p']/self.geo)).ravel()
    # Assemble all pressure, energy and enclosed-mass contributions together.
    data=m.data;B=data[:,4:8].reshape(-1,2,2);gs=data[:,16:22].reshape(-1,2,3);hs=data[:,22:28].reshape(-1,2,3)
    src=self.sources(m.points)
    g=[sum(sparse.diags(gs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
    h=[sum(sparse.diags(hs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
    F=sum(m.cov[i].T@sparse.diags(m.weights*B[:,i,j])@g[j] for i in range(2) for j in range(2))
    F-=sum(m.V[i].T@sparse.diags(m.weights)@h[i] for i in range(2))
    self.F=F.toarray().astype(np.longdouble)
    rf=np.array([self.edges[1]]);pf=m.bg.sample(rf);Nf,af=m.original.radiation.geometry.metric(rf)
    Af=np.exp(-2*pf['phi'][0]**2);self.luminosity_factor=float(4*np.pi*(100*R*rf[0])**2*Nf[0]**2*Af**4)
    self.face_volume=float(4*np.pi*(100*R*rf[0])**2*Nf[0]*af[0]*Af**4*(self.edges[2]-self.edges[0])*100*R/2)
    self.rate_clock=float(Af*Nf[0])*self.tc
    nf=m.surface_index+1;r=m.grid[:nf];point=m.bg.sample(r);N,a=m.original.radiation.geometry.metric(r)
    enclosed=self.sources(point)[2][:,:2].toarray()*(self.gr.C**4*R*(N*a)[:,None])/(self.gr.G*1e-7)
    shape=enclosed[:,0]/self.volumes[0]-enclosed[:,1]/self.volumes[1]
    shape[(r<=self.edges[0])|(r>=self.edges[-1])]=0.
    coeff=np.zeros(nf);inside=(shape!=0)&(r>0)&(point['p']>0)
    coeff[inside]=self.gr.G*1e-7*self.tc*self.luminosity_factor*shape[inside]/(4*np.pi*self.gr.C**4*R*N[inside]*a[inside]*np.exp(-8*point['phi'][inside]**2)*(point['e'][inside]+point['p'][inside])*r[inside]**3)
    HQ=np.zeros(m.size);HQ[m.indices[:nf,0]]=coeff
    self.HQ=np.asarray(m.to_hierarchical@HQ,np.longdouble)
    value,_=m.evaluation(rf)
    self.Cv=(value[0]*(af[0]/Nf[0]*rf[0]*self.gr.C)).tocsr()
    self.error=0.
