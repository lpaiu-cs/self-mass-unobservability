def __init__(self,n,flat=False):
    self.eos=EOS();self.env=dict(np.load(old.prior.OUT/'final-envelope.npz'));env=self.env;self.R=float(env['r'][-1]);self.As=float(env['A'][-1]);self.Ns=float(env['N'][-1]);self.bs=float(env['b'][-1]);self.flat=flat
    self.xf=np.linspace(-20000,120000,n+1);self.x=(self.xf[:-1]+self.xf[1:])/2;self.dx=self.xf[1]-self.xf[0]
    self.RJ=self.As*self.R;self.r=self.RJ+self.x;self.rf=self.RJ+self.xf
    self.area=(self.rf/self.RJ)**2 if not flat else np.ones(n+1)
    self.a0=self.As*self.Ns
    # Jordan areal metric: rJ=A*rE, bJ^1/2=1/[sqrt(bE)*(1+rE*alpha*Phi)].
    self.bg=old.Background()
    def metric(x):
        re=1+x/(self.As*self.R);p=self.bg.sample(re);A=np.exp(-2*p['phi']**2);b=1-2*p['m']/re
        a=A*p['N'];B=1/(np.sqrt(b)*(1-4*re*p['phi']*p['v']))
        nr=p['m']/(re*re*b)+4*np.pi*re*A**4*p['p']/b+re*p['v']**2/2
        ap=a*(nr-4*p['phi']*p['v'])/(self.R*A*(1-4*re*p['phi']*p['v']))
        if flat:a[:]=self.a0;B[:]=1.;ap[:]=0
        return a,B,ap
    self.a,self.B,self.ap=metric(self.x);self.af,self.Bf,_=metric(self.xf)
    self.vol=self.dx*self.B*(self.r/self.RJ)**2 if not flat else np.full(n,self.dx)
    d=self.eos.d;re=self.R+self.x/self.As
    rho=np.exp(np.interp(np.minimum(re,self.R),env['r'],np.log(env['rho'])))/self.eos.rho0
    sigma=np.interp(np.minimum(re,self.R),d['initial_r'],d['initial_sigma']);sigma[self.x>0]=0;rho[self.x>0]=0
    if flat:rho[self.x<0]=1.;sigma[:]=0.
    self.U=self.conserved(rho,np.zeros(n),sigma)[0]
    self.initial=self.U.copy();self.left=(float(rho[0]),0.,float(sigma[0]))
    previous=np.load(prior.prior.OUT/'p4-64.npz');self.hist_t=previous['emission_times']*self.bg.tc
    ri=(self.R+self.x[0]/self.As)/self.R
    self.hist_v=np.array([np.interp(ri,previous['radius'],np.asarray(v,float)) for v in previous['velocity']])*100/C
    self.F0=float(env['Linfinity'])/(4*np.pi*self.RJ**2*self.a0**2*C)/(self.eos.rho0*C*C)
    self.max_scattering_work=0.;self.max_optical=0.
