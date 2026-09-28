"""Counterexample candidate: actual shared material response to GR/photons.

Directional differences retain the existing nonlinear native recovery, shared
HLL face and limiter branches. Saved radiation forces drive this material
sweep; returning its motion to photons and GR is still required.
"""
from pathlib import Path
from types import MethodType
import inspect
import textwrap
import json
import signal
import sys
import time
import numpy as np
import sympy as sp
import def_native_monolithic_response as photons

flow=photons.old.flow;C=flow.C;write=flow.write;sha=flow.sha;LD=np.longdouble
OUT=photons.OUT.parent.parent/'def-native-material-response';AMP=1e-26


def replace(s,a,b):assert s.count(a)==1,(a,s.count(a));return s.replace(a,b)


class Material:
    def __init__(self,reference):
        self.reference=reference;self.model=m=flow.Coupled();b=m.bulk;f=m.flow;g=m.m;self.nb=b.n;self.n=b.n+m.n
        z=np.load(flow.OUT/f'coupled-{reference}.npz');self.d={k:z[k] for k in z.files if k.startswith('snapshot_') and k not in ['snapshot_I','snapshot_bulk_I']}
        self.metric=dict(np.load(photons.old.prior.metric.OUT/'corrected'/f'metric-{reference}-g8.npz'));self.t=self.metric['t'];self.rE=self.metric['radius_E']
        self.ids=[int(np.argmin(abs(self.d['snapshot_t']-t))) for t in self.t];assert np.max(abs(self.d['snapshot_t'][self.ids]-self.t))<1e-18
        self.a=np.r_[b.d['a'],g.a];self.V=np.r_[b.volume,4*np.pi*g.RJ*g.RJ*g.vol];self.R=np.r_[b.d['r'],g.r]
        self.Rf=np.r_[b.d['edges'][:-1],g.rf];self.af=np.r_[b.d['face_a'],g.af[1:]];self.B=np.r_[b.d['B'],g.B];self.Bf=np.r_[b.d['face_B'],g.Bf[1:]];self.ap=np.r_[m.f0['ap'],g.ap]
        bg=g.bg.fields(self.rE);self.den=1-4*bg['phi']*self.rE*bg['Phi'];self.Aden=np.exp(-2*bg['phi']**2)*self.den
        self.rEf=g.bg.edges[-(self.n+1):];self.rho0=b.d['rho'].copy();self.frho=m.face_rho.copy();self.f0=dict(m.f0)
        self.restore_set=m.set_material;self.volume_shift=np.zeros(self.n)
        def set_material(t):
            self.restore_set(t);b.eos.x=(1+b.eos.x)/(1+self.volume_shift[:self.nb])-1
        m.set_material=set_material
        # Extract fluxes before divergence; the shared face is retained once.
        src=flow.hydro
        anchor='    rate[1]+=C*np.diff(m.rf)/m.vol*(-(m.r/m.RJ)**2*E*m.ap+2*m.a*m.r/m.RJ**2*p)'
        src=replace(src,anchor,'    gravity=C*np.diff(m.rf)/m.vol*(-(m.r/m.RJ)**2*E*m.ap+2*m.a*m.r/m.RJ**2*p)\n    rate[1]+=gravity')
        src=replace(src,'    return rate,ledger,dt,V,kap','    return flux,gravity,dt,V')
        ns=dict(flow.hydro_ns);exec(compile(src,__file__,'exec'),ns);self.atmosphere=MethodType(ns['rhs'],f)
        method=inspect.getsource(m.material_rhs)
        if 'def material_rhs' not in method:raise AssertionError('Material owner unavailable')
        method=textwrap.dedent(method)
        method=replace(method,'return [mass,force,-np.diff(energy),-np.diff(neutral)]','return np.array([mass,C*momentum,energy,neutral]),C*(b.volume*self.f0["initial_support"]+geometry),self.deep_dt')
        ns=dict(m.material_rhs.__func__.__globals__);exec(compile(method,__file__,'exec'),ns);self.deep=MethodType(ns['material_rhs'],m)
        self.expanded=(src,method);self.cache={};self.raw_calls=0;self.probe_error=0.;self.min_probe=np.inf
        x=np.load(photons.OUT/f'steps-128-reference-{reference}.npz');idx=[int(np.argmin(abs(x['t']-t))) for t in self.t];rows=x['moments'][idx]
        self.transfer=np.stack([np.zeros_like(rows[:,1]),rows[:,3]/self.a,rows[:,1],rows[:,2]],axis=1)/AMP
        self.rest=m.m.a0*m.cx*C*C

    def geometry(self,field,eps):
        m=self.model;b=m.bulk;g=m.m;nb=self.nb
        u,ell,lam,aden,ap=field;interp=lambda v:np.interp(self.rEf,self.rE,v)
        uf=interp(u);lf=interp(lam);df=interp(aden)
        R=self.R*(1+eps*u);Rf=self.Rf*(1+eps*uf);a=self.a*(1+eps*ell);af=self.af*(1+eps*interp(ell))
        B=self.B*(1+eps*(lam-aden));Bf=self.Bf*(1+eps*(lf-df));V=self.V*(1+eps*(3*u+lam));apnew=self.ap+eps*ap
        assert np.min(V)>0 and np.min(a)>0
        self.volume_shift=eps*(3*u+lam)
        b.volume=V[:nb];b.d.update(r=R[:nb],edges=Rf[:nb+1],a=a[:nb],B=B[:nb],rho=self.rho0/(1+self.volume_shift[:nb]))
        for key,value in [('r',R[nb:]),('rf',Rf[nb:]),('a',a[nb:]),('af',af[nb:]),('B',B[nb:]),('Bf',Bf[nb:]),('ap',apnew[nb:]),('vol',V[nb:]/(4*np.pi*g.RJ*g.RJ))]:setattr(g,key,value)
        g.area=(g.rf/g.RJ)**2
        m.edge=Rf[:nb+1];m.area_gas=4*np.pi*m.edge**2*af[:nb+1];m.face_rho=self.frho/(1+eps*(3*uf[:nb+1]+lf[:nb+1]));m.base_momentum_flux=m.area_gas*m.face_p
        e0=m.mass0/self.V[:nb]*(m.cx*C*C+self.f0['u0']);p0=self.f0['p0'];pg=self.f0['pgprime']
        residual=self.f0['initial_support']+self.a[:nb]/self.B[:nb]*(pg+(e0+p0)*self.ap[:nb]/self.a[:nb])
        grad=pg/(1+eps*(u[:nb]+aden[:nb]));en=m.mass0/V[:nb]*(m.cx*C*C+self.f0['u0'])
        support=-a[:nb]/B[:nb]*(grad+(en+p0)*apnew[:nb]/a[:nb])
        support+=residual*(a[:nb]/B[:nb])/(self.a[:nb]/self.B[:nb])/(1+eps*(u[:nb]+aden[:nb]))
        m.f0=dict(self.f0,r=R[:nb],edges=Rf[:nb+1],a=a[:nb],B=B[:nb],ap=apnew[:nb],af=af[:nb+1],Bf=Bf[:nb+1],volume=V[:nb],initial_support=support,pgprime=grad)
        return V,a

    def point(self,k):
        if k in self.cache:return self.cache[k]
        m=self.model;b=m.bulk;f=m.flow;nb=self.nb;j=self.ids[k];d=self.d
        self.geometry(np.zeros((5,self.n)),0.);m.h=d['snapshot_h'][j].copy();m.Pi=d['snapshot_Pi'][j].copy();m.mass=d['snapshot_mass'][j].copy();m.set_material(self.t[k])
        theta=d['snapshot_theta'][j];eta=d['snapshot_eta'][j];u=d['snapshot_u'][j];state=m.material_state(u,eta);U=d['snapshot_U'][j]
        rho,v,lt,y=f.primitive(U);f.eos.y=y;p,uu,*_=f.eos(rho,lt);active=rho>=f.eos.floor
        Q=np.zeros((4,self.n));Q[:,:nb]=np.array([m.mass,m.Pi*C,state[2],state[3]])
        Q[:,nb:]=U*self.V[nb:]*f.eos.rho0*np.array([1,C*C,C*C,f.eos.nH])[:,None]
        pp=b.eos.gas(theta,eta)[0];Pg=np.r_[pp*self.V[:nb],p*f.eos.rho0*C*C*self.V[nb:]]
        Pr=Pg.copy();gamma=1/(1-v*v);enthalpy=(rho*(f.eos.cx+uu)+p)*f.eos.rho0*C*C
        Pr[nb:]+=self.V[nb:]*enthalpy*gamma*v*v
        row=dict(Q=Q,h=state[0],theta=theta,eta=eta,seed=lt,active=np.r_[np.ones(nb,bool),active],Pg=Pg,Pr=Pr)
        self.cache[k]=row
        flux,gravity,dt,_=self.raw(k,np.zeros_like(Q),np.zeros((5,self.n)),0.,row)
        row.update(flux=flux,gravity=gravity,dt=dt,base_rate=-np.diff(flux,axis=1));row['base_rate'][1]+=gravity
        # Match the original complete rate, independently of the flux export.
        original=f.hydro(U,self.t[k])[0]*self.V[nb:]*f.eos.rho0*np.array([1,C*C,C*C,f.eos.nH])[:,None]
        r=m.material_rhs(u,theta,eta);actual=np.c_[np.array([-np.diff(r[0]),r[1]*C,r[2],r[3]]),original]
        row['owner_error']=float(np.max(np.sum(abs(actual-row['base_rate']),axis=1)/np.maximum(np.sum(abs(actual),axis=1),1.)))
        return row

    def raw(self,k,delta,field,eps,row=None):
        row=self.point(k) if row is None else row;m=self.model;b=m.bulk;f=m.flow;nb=self.nb;self.raw_calls+=1
        V,a=self.geometry(field,eps);Q=row['Q']+eps*delta*np.asarray(row['active'])[None]
        E=(Q[2]+self.rest*Q[0])*(a/self.a)-self.rest*Q[0]
        h=row['h'].copy();h[1:]-=eps*np.cumsum(delta[0,:nb])
        state=[h,Q[1,:nb]/C,E[:nb],Q[3,:nb]]
        u,theta,eta=m.recover_material(state,self.t[k],row['theta'])
        U=Q[:,nb:]/(V[nb:]*f.eos.rho0*np.array([1,C*C,C*C,f.eos.nH])[:,None]);U[2]=E[nb:]/(V[nb:]*f.eos.rho0*C*C)
        f.seed=row['seed'].copy();flux,gravity,dt,prim=self.atmosphere(U,self.t[k])
        deep,dg,ddt=self.deep(u,theta,eta)
        factor=4*np.pi*m.m.RJ**2*f.eos.rho0*C*np.array([1,C*C,C*C,f.eos.nH])
        aflux=flux*factor[:,None];ag=gravity*V[nb:]*f.eos.rho0*C*C
        shared=float(np.max(abs(deep[:,-1]-aflux[:,0])/np.maximum(abs(aflux[:,0]),1.)))
        assert shared<1e-12,('Shared face',shared)
        return np.c_[deep[:,:-1],aflux],np.r_[dg,ag],min(dt,ddt),dict(theta=theta,eta=eta,primitive=prim)

    def fields(self,t):
        j=max(0,min(np.searchsorted(self.t,t,side='right')-1,15));f=(t-self.t[j])/(self.t[j+1]-self.t[j]);blend=lambda name:((1-f)*self.metric[name][j]+f*self.metric[name][j+1])/AMP
        u=blend('delta_u');ell=blend('delta_log_lapse');lam=blend('delta_lambda');aden=self.rE*blend('delta_u_prime')/self.den
        ap=self.ap*(ell-u-aden)+self.a/self.Aden*(blend('delta_nu_prime')+blend('delta_u_prime'))
        rates=(self.metric['delta_u'][j+1]-self.metric['delta_u'][j],self.metric['delta_lambda'][j+1]-self.metric['delta_lambda'][j])
        rates=np.array(rates)/(AMP*(self.t[j+1]-self.t[j]))
        return j,f,np.array([u,ell,lam,aden,ap]),rates

    def rhs(self,t,z,probe=1.):
        j,f,field,rates=self.fields(t);points=[(1-f,self.point(j)),(f,self.point(j+1))]
        size=max(np.max(abs(field[:4]))/1e-5,np.max(abs(field[4])/np.maximum(abs(self.ap),1e-50))/1e-3,1.)
        for _,p in points:
            q=p['Q'];active=p['active'];mass=np.maximum(q[0],1.)
            size=max(size,float(np.max(abs(z[0,active])/mass[active])/1e-5),float(np.max(abs(z[1,active])/(mass[active]*C*C)))/1e-7)
            thermal=z[2]-(self.a-self.model.m.a0)*self.model.cx*C*C*z[0]
            units=np.maximum(abs(q[2]-(self.a-self.model.m.a0)*self.model.cx*C*C*q[0]),1.)
            size=max(size,float(np.max(abs(thermal[active])/units[active]))/1e-5,float(np.max(abs(z[3,active])/np.maximum(abs(q[3,active]),1.)))/1e-5)
        eps=probe/size;self.min_probe=min(self.min_probe,eps);F=np.zeros((4,self.n+1));G=np.zeros((4,self.n));dt=np.inf;diff=[]
        for k,(weight,p) in enumerate(points,j):
            if weight==0:continue
            a=self.raw(k,z,field,eps);b=self.raw(k,z,field,eps/2)
            fa=(a[0].astype(LD)-p['flux'].astype(LD))/eps;fb=(b[0].astype(LD)-p['flux'].astype(LD))/(eps/2)
            ga=(a[1].astype(LD)-p['gravity'].astype(LD))/eps;gb=(b[1].astype(LD)-p['gravity'].astype(LD))/(eps/2)
            F+=weight*np.asarray(2*fb-fa,float);G[1]+=weight*np.asarray(2*gb-ga,float);dt=min(dt,a[2],b[2])
            diff.append(float(np.sum(abs(fb-fa))/max(np.sum(abs(fb)),1.)))
            G[2]-=weight*field[1]*(p['base_rate'][2]+self.rest*p['base_rate'][0])
            G[1]-=weight*(rates[0]+rates[1])*p['Q'][1]
            G[2]-=weight*self.a*(p['Pr']*(rates[0]+rates[1])+2*p['Pg']*rates[0])
        self.probe_error=max(self.probe_error,max(diff,default=0.))
        drive=(self.transfer[j+1]-self.transfer[j])/(self.t[j+1]-self.t[j]);G+=drive
        return -np.diff(F,axis=1)+G,F[:,0]-F[:,-1]+np.sum(G,axis=1),dt

    def active(self,t):
        j,f,_,_=self.fields(t);D=(1-f)*self.d['snapshot_U'][self.ids[j],0]+f*self.d['snapshot_U'][self.ids[j+1],0]
        return np.r_[np.ones(self.nb,bool),D>=self.model.flow.eos.floor]


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='5acae2934',
        claim='Apply the actual GR geometry and saved paired photon energy/H/momentum transfer to the shared native material mass,momentum,energy and neutral fluxes; evolve the additional material response across the original interior/atmosphere mesh.',
        decision='Determine the missing mass redistribution, velocity and pressure response before closing it back into radiation and GR. Do not replace material displacement by a small-force argument.',
        reuse='All531 material cells, corrected64/128 nonlinear backgrounds, Phase124 collision integrals, Phase123 metric increments and original3.434ms. No new EOS bank or full nonlinear path.',
        method='Forward Richardson directional flux differences through the actual native recovery, HLL shared face and limiter owners; retain direction-dependent branches. Apply SSP2 with the existing hydrodynamic CFL. Keep increments scaled separately by1e-26.',
        geometry='Fixed Einstein mesh with perturbed Jordan area radius, proper volumes, lapse, radial length and lapse gradient. Hold the fixed material-label inventory gradient. Retain the original pressure/gravity residual as a reference covector instead of imposing equilibrium.',
        time_metric='Physical radial momentum receives-(u_t+lambda_t)P. Reference material energy receives-a0*V*[Pr*(u_t+lambda_t)+2*Pt*u_t]. Convert the dynamic lapse energy of the existing owner back to a0 energy, including the background weighted-rate term.',
        boundary='Same shared HLL face and native vacuum floor. The derivative of the floor projection is recorded as a signed discard; no unrecorded deletion or artificial heat source.',
        limits='Radiation/metric histories are prescribed for this sweep. Additional fluid motion has not yet been returned into their evolution. Native derivative sampling and finite reconstruction are not uniform physical certification.',
        budgets=dict(check_seconds=35,pilot_seconds=35,production_seconds=160,CPU_threads=1,memory_GB=2,new_native_bank_calls=0),
        gates=dict(owner=1e-8,directional=.002,conservation=1e-8,time=.02,background_time=.02,small_relative_state=.000001),
        stop='Stop on support, original owner guards, directional/conservation error or wall cap. Measure actual prefixes before production with2x cost margin. No automatic extra clocks,domain,EOS support or weakened gates.',
        reference='https://arxiv.org/abs/gr-qc/0201064',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(flow.__file__),Path(photons.__file__),photons.OUT/'result.json',photons.OUT/'steps-128-reference-128.npz',photons.OUT/'steps-128-reference-64.npz']}))


def check():
    assert not (OUT/'check.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(35)
    m=Material(128);rows=[]
    for k in [0,8,16]:
        p=m.point(k);r,l,dt=m.rhs(m.t[k],np.zeros((4,m.n)))
        rows.append(dict(time=float(m.t[k]),owner_relative=p['owner_error'],flux_sum_relative=float(np.max(abs(np.sum(r,axis=1,dtype=LD)-l)/np.maximum(np.sum(abs(r),axis=1),1.))),hydro_CFL_seconds=float(dt),maximum_rates=np.max(abs(r),axis=1).tolist()))
    # Covariant conservation in a homogeneous anisotropically expanding cell.
    E,Pr,Pt,hr,ht,a,at,V=sp.symbols('E Pr Pt hr ht a at V')
    Vdot=V*(hr+2*ht);Edot=-(E+Pr)*hr-2*(E+Pt)*ht
    assert sp.expand(Vdot*E+V*Edot+V*(Pr*hr+2*Pt*ht))==0
    assert sp.expand(at*V*E+a*(Vdot*E+V*Edot)-at*V*E+a*V*(Pr*hr+2*Pt*ht))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Homogeneous zero-shift anisotropic metric work and conversion to fixed-reference lapse energy. A shared flux telescopes; this is not a complete nonlinear GR theorem.'))
    for name,source in zip(['atmosphere','deep'],m.expanded):(OUT/f'flux-owner-{name}.py').write_text(source)
    result=dict(classification='Counterexample candidate',passed=max(r['owner_relative'] for r in rows)<1e-8 and m.probe_error<.002,rows=rows,maximum_probe_difference=m.probe_error,minimum_probe=m.min_probe,raw_owner_calls=m.raw_calls,seconds=time.monotonic()-start)
    write(OUT/'check.json',result);print(json.dumps(result),flush=True);signal.alarm(0);assert result['passed']


if __name__=='__main__':globals()[sys.argv[1]]()
