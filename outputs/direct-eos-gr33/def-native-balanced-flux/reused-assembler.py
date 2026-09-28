def __init__(self,degree,resolved=False):
    started=time.monotonic();self.bg=bg=Background();d=np.load(OUT/'inputs.npz');th=np.load(OUT/'coefficients.npz')
    self.native=d['radius']/bg.R;self.edges=d['edges']/bg.R;self.dm=d['dm'];self.n=len(self.native)
    self.raw=th['raw'];self.thermo=th['thermo'];self.temperature=d['temperature'];self.core_count=int(d['core_count'])
    self.degree=degree;self.outer=1.1
    # ponytail: retain every thermal cell; mechanical p2/p4 controls share
    # the registered mesh, not an automatically expanded acoustic mesh.
    k=self.core_count
    ci=np.unique(np.r_[np.arange(0,k,8),np.arange(max(0,k-40),k+1)])
    self.cells=np.unique(np.r_[self.edges[ci],self.edges[k::4],1.,np.linspace(1,1.1,257)])
    if resolved:
        ids=np.arange(k-40,k+40)
        extra=(self.edges[ids,None]+np.diff(self.edges)[ids,None]*np.arange(1,64)[None,:]/64).ravel()
        self.cells=np.unique(np.r_[self.edges,extra,np.linspace(1,1.1,257)])
    self.cells=np.unique(np.r_[self.edges,self.native,np.linspace(1,1.1,257)])
    self.local,self.polynomials=hierarchy.task.basis(degree);dx=np.diff(self.cells)
    self.grid=np.r_[(self.cells[:-1,None]+dx[:,None]*self.local[:-1]).ravel(),self.outer]
    self.surface_index=int(np.where(self.cells==1)[0][0])*degree
    nf=self.surface_index+1;ng=len(self.grid);self.size=ng+nf-1
    self.indices=np.full((ng,2),-1,int);self.indices[:nf,0]=2*np.arange(nf);self.indices[:nf,1]=2*np.arange(nf)+1
    self.indices[nf:,1]=2*nf+np.arange(ng-nf);self.indices[-1,1]=-1
    cuts=np.unique(np.r_[self.cells,bg.x,self.edges,self.native]);gx,gw=np.polynomial.legendre.leggauss(6)
    points=(cuts[:-1,None]+np.diff(cuts)[:,None]*(gx+1)/2).ravel()
    self.weights=(np.diff(cuts)[:,None]*gw/2).ravel();self.points=bg.sample(points)
    import sympy as sp
    saved=json.loads((OUT.parent/'def-gr-canonical/regular/symbolic.json').read_text())
    names='r m p e phi v Gamma A4 C e_prime Gamma_prime';symbols=sp.symbols(names)
    expr=[sp.sympify(x,locals=dict(zip(names.split(),symbols))) for x in saved['expressions']]
    fn=sp.lambdify(symbols,expr,'numpy',cse=True)
    p=self.points;r=points;b=1-2*p['m']/r;c=p['N']*np.sqrt(b);inside=r<1
    ids=np.clip(np.searchsorted(bg.x,r[inside],side='right')-1,0,len(bg.x)-2)
    ep=np.diff(bg.saved['e'])[ids]/np.diff(bg.x)[ids];gp=np.diff(bg.saved['gamma'])[ids]/np.diff(bg.x)[ids]
    data=np.zeros((len(r),36));vals=fn(*[p[x][inside] for x in ['r','m','p','e','phi','v','gamma']],np.exp(-8*p['phi'][inside]**2),c[inside],ep,gp)
    data[inside]=np.array([np.broadcast_to(v,ep.shape) for v in vals]).T
    data[~inside,7]=r[~inside]**2*c[~inside];data[~inside,11]=-2*r[~inside]**2*c[~inside]*p['v'][~inside]**2/b[~inside]
    data[~inside,15]=r[~inside]**2/c[~inside];data[~inside,27]=-2*c[~inside]*p['v'][~inside]/b[~inside]**2
    assert np.all(np.isfinite(data));self.V,self.D=self.evaluation(points)
    A=data[:,:4].reshape(-1,2,2);B=data[:,4:8].reshape(-1,2,2);CC=data[:,8:12].reshape(-1,2,2);W=data[:,12:16].reshape(-1,2,2)
    cov=[self.D[i]-sum(diags(A[:,i,j])@self.V[j] for j in range(2)) for i in range(2)]
    def form(left,coeff,right):return sum(left[i].T@diags(self.weights*coeff[:,i,j])@right[j] for i in range(2) for j in range(2)).tocsc()
    self.K=form(cov,B,cov)+form(self.V,CC,self.V);self.M=form(self.V,W,self.V)
    # Proper redshift volume, evaluated with the same input geometry.
    vx=(self.edges[:-1,None]+np.diff(self.edges)[:,None]*(gx+1)/2)
    vp=bg.sample(vx.ravel());vb=1-2*vp['m']/np.maximum(vp['r'],1e-100)
    self.volumes=4*np.pi*bg.R**3*np.diff(self.edges)/2*np.sum((vp['N']/np.sqrt(vb)*vp['r']**2).reshape(vx.shape)*gw,axis=1)
    gs=data[:,16:22].reshape(-1,2,3);hs=data[:,22:28].reshape(-1,2,3);src=self.sources(p)
    g=[sum(diags(gs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
    h=[sum(diags(hs[:,i,j])@src[j] for j in range(3)) for i in range(2)]
    if resolved:
        self.F=sum(cov[i].astype(np.longdouble).T@diags(self.weights.astype(np.longdouble)*B[:,i,j])@g[j].astype(np.longdouble) for i in range(2) for j in range(2))
        self.F-=sum(self.V[i].astype(np.longdouble).T@diags(self.weights.astype(np.longdouble))@h[i].astype(np.longdouble) for i in range(2))
    else:
        self.F=sum(cov[i].T@diags(self.weights*B[:,i,j])@g[j] for i in range(2) for j in range(2))-sum(self.V[i].T@diags(self.weights)@h[i] for i in range(2))
    pr=bg.sample(self.grid[:nf]);rn=pr['r'];bn=1-2*pr['m']/np.maximum(rn,1e-100)
    factor=np.zeros(nf);valid=rn>0
    factor[valid]=G/(4*np.pi*C**4*bg.R*pr['N'][valid]/np.sqrt(bn[valid])*np.exp(-8*pr['phi'][valid]**2)*(pr['e'][valid]+pr['p'][valid])*rn[valid]**3)
    hmap=(diags(factor)@self.face_map(rn)).tolil()
    hmap[0,0]=G/(4*np.pi*C**4*bg.R*pr['N'][0]*np.exp(-8*pr['phi'][0]**2)*(pr['e'][0]+pr['p'][0])*self.edges[1]**3)
    hz=hmap.tocoo()
    H=coo_matrix((hz.data,(self.indices[hz.row,0],hz.col)),shape=(self.size,self.n)).tocsc()
    # Convert only the heat lift; evaluation already uses endpoint/bubbles.
    rows=[];cols=[];values=[]
    for j in range(1,degree):
        nodes=np.arange(len(self.cells)-1)*degree+j
        for shift,value in [(-j,-(1-self.local[j])),(degree-j,-self.local[j])]:
            a=self.indices[nodes,0];z=self.indices[nodes+shift,0];valid=(a>=0)&(z>=0)
            rows.extend(a[valid]);cols.extend(z[valid]);values.extend(np.full(valid.sum(),value))
    self.H=(eye(self.size)+coo_matrix((values,(rows,cols)),shape=(self.size,self.size)))@H
    self.nativeV,self.nativeD=self.evaluation(self.native)
    self.temperature_maps()
    # Finite material pressure requires a scalar traction jump, even when
    # the emitted rays are evaluated on the background geometry.
    ps=bg.sample(np.array([1.]));bs=1-2*ps['m'][0];cs=ps['N'][0]*np.sqrt(bs)
    aa=np.exp(-8*ps['phi'][0]**2);alpha=-4*ps['phi'][0];v=ps['v'][0];pressure=ps['p'][0]
    boundary,_=self.evaluation(np.array([1.]));z,f=boundary
    jump=8*np.pi*aa*cs*pressure*(2*alpha+v)/bs
    env=np.load(prior.OUT/'final-envelope.npz');Pr=float(env['Prad'][-1])*G*bg.R**2/C**4
    traction=-16*np.pi*aa*cs*Pr/bs
    self.K=(self.K-f.T@(jump*z)-z.T@(traction*self.TLq[-1])).tocsc()
    self.F=(self.F+z.T@(traction*self.TE[-1])).tocsc()
    # Actual nonzero grey constitutive currents and their temperature loop.
    pp=bg.sample(self.native);AN=np.exp(-2*pp['phi']**2)*pp['N'];self.theta=AN*self.temperature
    face=bg.sample(self.edges[1:-1]);af=np.exp(-2*face['phi']**2);bf=1-2*face['m']/face['r']
    arad=prior.envelope.Envelope().arad
    conductivity=4*arad*C*self.temperature**3/(3*self.raw[:,0]*th['opacity'])
    pref=-4*np.pi*(face['r']*bg.R)**2*face['N']*af**2*np.sqrt(bf)*np.sqrt(conductivity[:-1]*conductivity[1:])/np.diff(d['radius'])
    diff=coo_matrix((np.tile([-1.,1.],self.n-1),(np.repeat(np.arange(self.n-1),2),np.c_[np.arange(self.n-1),np.arange(1,self.n)].ravel())),shape=(self.n-1,self.n)).tocsr()
    op=diags(bg.tc*pref)@diff@diags(self.theta)
    self.coefficient_data=data;self.cov=cov;self.thermal_operator=op
    from scipy.sparse import vstack
    self.L0=np.r_[pref*np.diff(self.theta),float(env['Linfinity'])]
    self.L0[k-1:]=float(env['Linfinity'])
    self.Lq=vstack([op@self.Tq,4*bg.tc*self.L0[-1]*self.TLq[-1]]).tocsc()
    self.LE=vstack([op@self.TE,4*bg.tc*self.L0[-1]*self.TE[-1]]).tocsc()
    self.lam=bg.tc*C*np.sqrt((self.raw[:-1,0]*th['opacity'][:-1])*(self.raw[1:,0]*th['opacity'][1:]))*np.sqrt(AN[:-1]*AN[1:])
    self.active=np.arange(self.n);self.energy_scale=1/np.maximum(np.asarray(abs(self.LE).max(axis=0).toarray()).ravel(),1e-100)
    positions=np.zeros(self.size)
    for field in range(2):
        valid=self.indices[:,field]>=0;positions[self.indices[valid,field]]=self.grid[valid]
    self.permutation=np.argsort(np.r_[positions,self.edges[1:]],kind='stable')
    outside=points>1;self.outside=outside;self.rays=bg.rays(points[outside],48)
    self.photon_test=self.V[1][outside].T@diags(self.weights[outside])
    self.history_t=[0.];self.history_E=[0.];self.history_f=[self.L0[-1]*bg.tc]
    self.max_error=self.max_heat_error=0.;self.setup_seconds=time.monotonic()-started
