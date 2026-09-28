"""Analytic deep central-flux tangent; retain the actual shared HLL owner.

Counterexample candidate. Differentiate the existing retained constitutive
table and geometry directly, without subtracting large deep background forces.
"""
import inspect,resource,sys,time
import numpy as np
import sympy as sp

LD=np.longdouble


def deep_tangent(self,k,z,field):
    m=self.model;b=m.bulk;n=self.nb;C=LD(29979245800.);q=self.point(k)['Q'].astype(LD)
    if not hasattr(self,'deep_tangent_cache'):self.deep_tangent_cache={}
    if k not in self.deep_tangent_cache:
        raw=self.raw(k,np.zeros_like(q),np.zeros_like(field),0.)
        theta,eta=raw[3]['theta'],raw[3]['eta'];p,u,ut,uy,pt,py,*_=b.eos.gas(theta,eta)
        self.deep_tangent_cache[k]=dict(p=p.copy(),u=u.copy(),ut=ut.copy(),uy=uy.copy(),pt=pt.copy(),py=py.copy(),eta=eta.copy(),
            drp=b.eos.coeff('pr').copy()*(1+b.eos.x),dru=b.eos.coeff('ur').copy()*(1+b.eos.x))
    bg=self.deep_tangent_cache[k];p,u=[bg[key].astype(LD) for key in ['p','u']]
    a,V,B,R=[np.asarray(v[:n],LD) for v in [self.a,self.V,self.B,self.R]]
    gu,ell,lam,aden,ap=field.astype(LD);vol=3*gu+lam
    db=z[0,:n]/q[0,:n];beta=q[1,:n]/(m.cx*q[0,:n]*C*C)
    dv=z[1,:n]/(m.cx*q[0,:n]*C)-beta*C*db
    du=(z[2,:n]-(a-m.m.a0)*m.cx*C*C*z[0,:n])/(a*q[0,:n])
    du-=(u+LD('.5')*m.cx*C*C*beta*beta)*db+m.cx*C*beta*dv
    dy=z[3,:n]/q[3,:n]-db;dr=db-vol[:n]
    xi=m.mech.xi@np.r_[LD(0),-np.cumsum(z[0,:n],dtype=LD)]
    deta=(1+bg['eta'])*dy
    dt=(du-bg['dru']*dr-bg['uy']*deta+b.eos.inventory[1]*xi)/bg['ut']
    dp=bg['drp']*dr+bg['pt']*dt+bg['py']*deta-b.eos.inventory[0]*xi
    def face(v):return np.asarray(np.interp(self.rEf,self.rE,np.asarray(v,float))[:n+1],LD)
    uf,ef,lf=face(gu),face(ell),face(lam)
    af=np.asarray(self.af[:n+1],LD);rho=np.asarray(self.frho,LD)
    A=4*LD(np.pi)*np.asarray(self.Rf[:n+1],LD)**2*af
    dA=A*(2*uf+ef);vf=(np.r_[LD(0),beta*C]+np.r_[beta*C,LD(0)])/2
    dvf=(np.r_[LD(0),dv]+np.r_[dv,LD(0)])/2;vf[0]=0;dvf[0]=0
    pg=p-self.f0['p0'];ps=(np.r_[pg[0],pg]+np.r_[pg,pg[-1]])/2
    dps=(np.r_[dp[0],dp]+np.r_[dp,dp[-1]])/2
    mass=A*rho*vf;dmass=A*rho*(dvf+(ef-uf-lf)*vf)
    momentum=A*(m.face_p+ps+m.cx*rho*vf*vf)
    dmomentum=dA*(m.face_p+ps+m.cx*rho*vf*vf)+A*(dps+m.cx*rho*(2*vf*dvf-(3*uf+lf)*vf*vf))
    donor=np.where(mass==0,dmass>=0,mass>=0)
    pick=lambda v:np.where(donor,np.r_[v[0],v],np.r_[v,v[-1]])
    donor_u,delta_u=pick(u),pick(du)
    nh=b.d['thermo'][:,4];y=b.d['y0']*(1+bg['eta']);donor_h=pick(nh*y);delta_h=pick(nh*y*dy)
    enthalpy=(af-m.m.a0)*m.cx*C*C+af*(donor_u+(m.face_p+ps)/rho)
    denthalpy=af*ef*(m.cx*C*C+donor_u+(m.face_p+ps)/rho)+af*(delta_u+dps/rho+(m.face_p+ps)*(3*uf+lf)/rho)
    flux=np.array([dmass,C*dmomentum,dmass*enthalpy+mass*denthalpy,dmass*donor_h+mass*delta_h])
    e0=np.asarray(m.mass0,LD)/V*(m.cx*C*C+self.f0['u0']);p0=np.asarray(self.f0['p0'],LD)
    ap0=np.asarray(self.ap[:n],LD);u0=np.asarray(b.u0,LD)
    support=V*(self.f0['initial_support']*(2*gu[:n]+ell[:n])+ap0/B*(e0*(2*gu[:n]+ell[:n]+lam[:n]-aden[:n])+p0*(ell[:n]-gu[:n]-aden[:n]))-(e0+p0)/B*ap[:n])
    N=(q[0,:n]-m.mass0)*m.cx*C*C+q[0,:n]*u-m.mass0*u0
    dN=z[0,:n]*(m.cx*C*C+u)+q[0,:n]*du
    geom=-ap0/B*dN-N/B*(ap[:n]-ap0*(lam[:n]-aden[:n]))+2*V*a/(B*R)*(dp+pg*(2*gu[:n]+ell[:n]+aden[:n]))
    gravity=C*(np.diff(dA*m.face_p)+support+geom)
    # The shared face n is deliberately left to the common atmospheric HLL
    # owner. Only faces0..n-1 and the deep volume source are analytic here.
    return flux[:,:n],gravity


def rhs(owner):
    source=owner.rhs_source
    anchor="ga=(a[1].astype(LD)-p['gravity'].astype(LD))/eps;gb=(b[1].astype(LD)-p['gravity'].astype(LD))/(eps/2)"
    assert source.count(anchor)==1
    source=source.replace(anchor,anchor+"\n        df,dg=self.deep_tangent(k,z,field)\n        fa[:,:self.nb]=df;fb[:,:self.nb]=df;ga[:self.nb]=dg;gb[:self.nb]=dg")
    namespace=dict(owner.ns);exec(compile(source,__file__,'exec'),namespace)
    return namespace['rhs'],source


def symbolic():
    h,u,e,l,d,a,B,V,p,pg,E,ap,dap,S=sp.symbols('h u e l d a B V p pg E ap dap S')
    vr=V*(1+h*(3*u+l));aa=a*(1+h*e);bb=B*(1+h*(l-d));apr=ap+h*dap
    residual=S+a/B*(pg+(E+p)*ap/a)
    expression=vr*(-aa/bb*(pg/(1+h*(u+d))+(E*V/vr+p)*apr/aa)+residual*(aa/bb)/(a/B)/(1+h*(u+d)))
    tangent=V*(S*(2*u+e)+ap/B*(E*(2*u+e+l-d)+p*(e-u-d))-(E+p)/B*dap)
    assert sp.simplify(sp.diff(expression,h).subs(h,0)-tangent)==0
    return dict(classification='Proven',passed=True,scope='Exact derivative of the existing deep geometric support, including retained background residual. Constitutive jets differentiate the retained table; no full native EOS or continuum certificate.')


if __name__=='__main__':
    import solve_native_incident_reciprocal as solve
    sweep=2;out=solve.OUT;start=time.monotonic();cpu=time.process_time();error=None
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));solve.base.drive.native.deadline(45)
    try:
        solve.initialize(sweep);m=solve.Material(128,128);rows=[]
        samples=[(out/'sweep-2/material-probed/pilot-128.npz',-1),
                 (out/'sweep-1/material-probed/steps-128-reference-128.npz',64),
                 (out/'sweep-1/material-probed/steps-128-reference-128.npz',128)]
        for path,i in samples:
            data=np.load(path);t=float(data['t'][i]);z=data['history_scaled'][i]
            exact=m.rhs(t,z)[0];original=m.directional_rhs
            m.directional_rhs=solve.base.old.aligned.ns['rhs'].__get__(m,type(m))
            try:reference=[m.rhs(t,z,p)[0] for p in [1.,2.,4.]]
            finally:m.directional_rhs=original
            norm=np.maximum(np.sum(abs(exact),axis=1),1.)
            errors=[(np.sum(abs(v-exact),axis=1)/norm).astype(float).tolist() for v in reference]
            rows.append(dict(time=t,finite_probe_512_1024_2048=errors))
            assert max(errors[-1])<.002,rows[-1]
        result=dict(classification='Counterexample candidate',passed=True,rows=rows,symbolic=symbolic(),source_sha256=solve.sha(__file__))
        solve.write(out/'deep-tangent-check.json',result);print(result,flush=True)
    except Exception as exc:error=repr(exc);raise
    finally:solve.write(out/'deep-tangent-check-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=solve.sha(__file__)))
