"""Counterexample candidate: compensated radiation/heat/H response.

Reuse the actual moving collision kernel and corrected native EOS. Density,
velocity and advected inventories are prescribed by the accepted trajectory;
their *additional* mechanical response remains open, with its impulse saved.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu, LinearOperator, gmres
import def_native_compensated_transport as prior

flow=prior.flow; C=prior.C; write=prior.write; sha=prior.sha
OUT=prior.metric.OUT.parent/'def-native-collision-response'
AMPLITUDE=1e-26


def prepare():
    assert not OUT.exists(); OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='e71ed7eba',
        claim='Evolve actual moving absorption/emission/Thomson, native thermal energy and neutral H response together with the GR-driven photon transport, using separately scaled increments on all531 cells.',
        decision='Determine whether thermochemical feedback amplifies or screens the geometric photon response; save the paired material energy/species and momentum impulse for the remaining mechanical and GR closure.',
        reuse='Reuse the corrected64/128 backgrounds,17 lapse/source times,531 cells,8 angles,152 frequencies and original3.434ms horizon. No new native bank or nonlinear trajectory.',
        boundary='Same fixed deep incoming packet perturbation and external vacuum. Density, velocity and other advected inventories follow the saved background. The additional force is measured, not silently treated as a free fluid evolution.',
        method='Strang composition of shared implicit transport and a coupled local linear radiation/heat/H solve. Actual moving scattering is retained. A static Thomson plus absorption/heat/H Schur inverse preconditions the local solve. Time comparison decides this partition.',
        metric='Use delta_F=delta_N/measure-delta_log_volume*F0. Include the volume/time-rate change and canonical density/isentropic temperature response in collisions. Do not put frequency-box losses or gravitational work into gas heat.',
        derivatives='Centered relative temperature/neutral/density coefficient probes in the existing native interpolation model; compare1e-5 and5e-6 at representative times. This is sampled derivative evidence, not a uniform EOS derivative bound.',
        gates=dict(derivative=.0001,owner=1e-10,linear=1e-10,energy=1e-8,species=1e-8,time=.02,background_time=.02),
        budget=dict(inspect_seconds=25,bank_seconds=90,pilot_seconds=45,production_seconds=240,CPU_threads=1,memory_GB=3,new_native_bank_calls=0),
        production='Only64/128 on background128 and128 on background64, after a measured pilot with2x projected cost fits240s. Halt on budget or gate failure, preserve outputs; no automatic extra mesh, paths, horizon or relaxed tolerance.',
        limitations='Prescribed material motion and finite EOS/angle/frequency model; no full hydrodynamic/GR feedback or final physical charge.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(prior.__file__),Path(flow.__file__),prior.OUT/'result.json',flow.OUT/'coupled-128.npz',flow.OUT/'coupled-64.npz',prior.metric.OUT/'corrected/metric-128-g8.npz',prior.metric.OUT/'corrected/metric-64-g8.npz']}))


class Response(prior.Transport):
    def __init__(self,reference):
        super().__init__(reference); self.data=dict(np.load(flow.OUT/f'coupled-{reference}.npz'))
        self.keep=np.array([np.argmin(abs(self.data['snapshot_t']-t)) for t in self.t])
        b=self.model.bulk; self.scale=b.scale.copy(); self.a=np.r_[b.d['a'],self.model.m.a]
        self.volume=4*np.pi*self.W*self.a**3; self.nb=b.n
        self.Eweight=self.weights*self.E*self.scale
        self.Nweight=self.weights*self.scale
        self.cache={}

    def restore(self,k):
        m=self.model; b=m.bulk; f=m.flow; j=self.keep[k]; d=self.data
        m.h=d['snapshot_h'][j].copy(); m.Pi=d['snapshot_Pi'][j].copy(); m.mass=d['snapshot_mass'][j].copy(); m.set_material(self.t[k])
        self.theta=d['snapshot_theta'][j].copy(); self.eta=d['snapshot_eta'][j].copy()
        self.rho,self.beta,self.lt,self.y=f.primitive(d['snapshot_U'][j]); self.active=self.rho>=f.eos.floor
        self.bulk_beta=m.velocity(); self.bulk_x=b.eos.x.copy()

    def coefficients(self,dr=0.,dt=0.,dy=0.):
        """Existing bound-free owners, with relative density/T/neutral probes."""
        m=self.model; b=m.bulk; f=m.flow; nb=self.nb
        b.eos.x=self.bulk_x+(1+self.bulk_x)*np.expm1(dr)
        theta=self.theta+dt; eta=(1+self.eta)*np.exp(dy)-1
        p,u,ut,uy,pt,py,ne,*_=b.eos.gas(theta,eta)
        ab,em,*_=b.eos.radiation(theta,eta)
        df=b.eos.spectral['frequency'][0]*(1-b.eos.f)+b.eos.spectral['frequency'][1]*b.eos.f
        boost=-self.bulk_beta[:,None,None]*self.mu[None,:,None]
        absorb=b.factor[:,None,None]*np.exp(dr)*ab[:,None,:]*(1+boost*(1+df[:,0,None,:]))
        emit=b.factor[:,None,None]*np.exp(dr)*em[:,None,:]*(1+boost*(1+df[:,1,None,:]))
        sc=b.scfactor*ne
        energy=m.mass*b.d['a']*u
        neutral=m.mass*b.d['thermo'][:,4]*b.d['y0']*(1+eta)
        rho=m.mass/b.volume*np.exp(dr); kap=ne*6.6524587321e-25/rho
        n=self.n; E=np.zeros((n,self.q,self.nf)); loss=np.zeros_like(E)
        E[:nb]=emit; loss[:nb]=absorb-emit
        ss=np.zeros(n); ss[:nb]=sc
        gg=np.zeros(n); gg[:nb]=energy; hh=np.zeros(n); hh[:nb]=neutral
        rr=np.zeros(n); rr[:nb]=rho; vv=np.r_[self.bulk_beta,self.beta]; kk=np.zeros(n); kk[:nb]=kap
        specific=np.zeros(n); specific[:nb]=u; pressure=np.zeros(n); pressure[:nb]=p
        idx=np.flatnonzero(self.active); rho=self.rho[idx]*np.exp(dr); lt=self.lt[idx]+dt; y=self.y[idx]*np.exp(dy)
        f.eos.y=y
        if len(idx):
            pp,uu,_,_,kap,_,_=f.eos.evaluate(rho,lt); rr0=rho*f.eos.rho0; a=self.a[nb+idx]
            gamma=1/np.sqrt(1-self.beta[idx]**2); D=gamma[:,None]*(1-self.beta[idx,None]*self.mu)
            ab,em=m.spectrum.coefficients(rr0,lt,y,self.E[None,None,:]*D[:,:,None]/a[:,None,None])
            factor=a[:,None,None]*C*rr0[:,None,None]*f.eos.nH*D[:,:,None]
            E[nb+idx]=factor*em; loss[nb+idx]=factor*(ab-em)
            ss[nb+idx]=a*C*rr0*kap
            # Exact lab energy derivative at prescribed density/velocity;
            # its opposite photon transfer is evolved, momentum is retained.
            pp=pp*f.eos.rho0*C*C; eps=rr0*uu*C*C
            # Remove the unchanged rest energy before taking thermal jets.
            gg[nb+idx]=a*self.volume[nb+idx]*((eps+pp)*gamma**2-pp)
            hh[nb+idx]=self.volume[nb+idx]*gamma*rr0*f.eos.nH*y
            rr[nb+idx]=rr0; kk[nb+idx]=kap
            specific[nb+idx]=uu*C*C; pressure[nb+idx]=pp
        b.eos.x=self.bulk_x
        return dict(emit=E,loss=loss,sc=ss,energy=gg,neutral=hh,rho=rr,beta=vv,kap=kk,u=specific,p=pressure)

    def scattering_matrix(self,c):
        """Exact sparse form of the existing positive two-node packet map."""
        n,q,nf=self.n,self.q,self.nf; v=c['beta']; gamma=1/np.sqrt(1-v*v)
        edges=(self.edges_mu[None,:]-v[:,None])/(1-v[:,None]*self.edges_mu[None,:]); dw=np.diff(edges)/2
        p2=(edges[:,:-1]**2+edges[:,:-1]*edges[:,1:]+edges[:,1:]**2-1)/2
        rate=c['sc'][:,None]*gamma[:,None]*(1-v[:,None]*self.mu)
        ids=np.arange(n*q*nf).reshape(n,q,nf); rows=[ids.ravel()]; cols=[ids.ravel()]; vals=[np.broadcast_to(-rate[:,:,None],ids.shape).ravel()]
        escaped=np.zeros((3,n,q,nf)); packet=self.w[None,:,None]*self.num*self.scale
        for j,muin in enumerate(self.mu):
            for k,muout in enumerate(self.mu):
                probability=dw[:,k]*(1+.5*p2[:,j]*p2[:,k])
                dest=self.E[None,:]*(1-v[:,None]*muin)/(1-v[:,None]*muout)
                high=np.searchsorted(self.extended,dest); low=high-1
                assert high.min()>0 and high.max()<len(self.extended)
                frac=(dest-self.extended[low])/(self.extended[high]-self.extended[low])
                for ind,w in [(low,1-frac),(high,frac)]:
                    inside=(ind>0)&(ind<=nf); destid=np.clip(ind-1,0,nf-1)
                    gain=rate[:,j,None]*probability[:,None]*w
                    value=gain*self.w[j]/self.w[k]*(self.num*self.scale)[None,:]/(self.num*self.scale)[destid]
                    rows.append(ids[:,k][np.arange(n)[:,None],destid][inside]); cols.append(ids[:,j][inside]); vals.append(value[inside])
                    lost=gain*(~inside)*packet[:,j]
                    escaped[0,:,j]+=lost; escaped[1,:,j]+=lost*self.extended[ind]; escaped[2,:,j]+=lost*self.extended[ind]*muout
        S=sparse.coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(ids.size,ids.size)).tocsr(); S.eliminate_zeros()
        return S,escaped* (4*np.pi*self.W)[None,:,None,None]

    def inspect(self):
        started=time.monotonic(); rows=[]
        for k in [0,8,16]:
            self.restore(k); c=self.coefficients(); then=time.monotonic(); S,esc=self.scattering_matrix(c)
            I=self.I[k]; X=I/self.scale; actual=(S@X.ravel()).reshape(X.shape)*self.scale
            expected,escape=self.model.scattering(I,c['rho'],c['beta'],c['kap'],self.a)
            weight=self.weights*self.E
            error=float(np.sum(abs(actual-expected)*weight)/max(np.sum(abs(expected)*weight),1.))
            actualesc=np.einsum('knqf,nqf->kn',esc,X)
            target=escape.T*self.volume
            # Owner escape is per proper volume, with Eref units.
            err=float(np.max(abs(actualesc-target))/max(np.max(abs(target)),1.))
            rows.append(dict(time=float(self.t[k]),active=int(self.active.sum()),matrix_nnz=int(S.nnz),matrix_seconds=time.monotonic()-then,owner_relative=error,escape_relative=err,max_loss=float(c['loss'].max()),max_sc=float(c['sc'].max())))
        row=dict(classification='Counterexample candidate',seconds=time.monotonic()-started,rows=rows,passed=max(max(r['owner_relative'],r['escape_relative']) for r in rows)<1e-10)
        write(OUT/'inspect.json',row); print(json.dumps(row),flush=True); assert row['passed']


def inspect():
    signal.signal(signal.SIGALRM,flow.old.optical.timeout); signal.alarm(25)
    Response(128).inspect(); signal.alarm(0)


def bank():
    assert not (OUT/'bank-result.json').exists(); start=time.monotonic()
    signal.signal(signal.SIGALRM,flow.old.optical.timeout); signal.alarm(90)
    errors=[]; files=[]
    for reference in [128,64]:
        m=Response(reference); folder=OUT/f'bank-{reference}'; folder.mkdir()
        Eunit=[]; Nunit=[]
        for k in range(17):
            m.restore(k); c=m.coefficients(); jets={}
            for param in ['dr','dt','dy']:
                plus=m.coefficients(**{param:1e-5}); minus=m.coefficients(**{param:-1e-5})
                jets[param]={key:(plus[key]-minus[key])/2e-5 for key in ['emit','loss','sc','energy','neutral','u','p']}
                if reference==128 and k in [0,8,16]:
                    plus=m.coefficients(**{param:5e-6}); minus=m.coefficients(**{param:-5e-6})
                    for key in ['emit','loss','sc','energy','neutral','u','p']:
                        fine=(plus[key]-minus[key])/1e-5; coarse=jets[param][key]
                        if key in ['emit','loss']:
                            w=m.weights*m.E*(1+m.I[k])
                            err=np.sum(abs(fine-coarse)*w)/max(np.sum(abs(fine)*w),1.)
                        else: err=np.max(abs(fine-coarse))/max(np.max(abs(fine)),1.)
                        errors.append(dict(reference=reference,k=k,parameter=param,quantity=key,relative=float(err)))
            active=c['rho']>0; assert np.all(jets['dt']['energy'][active]>0)
            gammaT=np.divide(c['p']/np.maximum(c['rho'],1e-300)-jets['dr']['u'],jets['dt']['u'],out=np.zeros(m.n),where=active)
            record=dict(c,adiabatic_logT=gammaT)
            record.update({param+'_'+key:value for param,v in jets.items() for key,value in v.items()})
            p=folder/f'point-{k}.npz'; np.savez_compressed(p,**record); files.append(str(p))
            Eunit.append(jets['dt']['energy']); Nunit.append(c['neutral'])
        np.savez_compressed(folder/'units.npz',energy=np.maximum(np.max(Eunit,axis=0),1.),neutral=np.maximum(np.max(Nunit,axis=0),1.))
    result=dict(classification='Counterexample candidate',passed=max(e['relative'] for e in errors)<1e-4,derivative_checks=errors,seconds=time.monotonic()-start,
        new_native_bank_calls=0,coefficient_model='existing corrected native interpolation, centered relative jets',full_uniform_derivative_bound=False)
    write(OUT/'bank-result.json',result); signal.alarm(0); print(json.dumps({k:v for k,v in result.items() if k!='derivative_checks'}),flush=True); assert result['passed']


class CoupledResponse(Response):
    def __init__(self,reference):
        super().__init__(reference); units=np.load(OUT/f'bank-{reference}/units.npz')
        self.eu=units['energy']; self.nu=units['neutral']; self.size=self.n*self.q*self.nf
        self.max_iterations=0; self.max_residual=0.; self.last_key=None

    def point(self,k):
        if k in self.cache:return self.cache[k]
        c=dict(np.load(OUT/f'bank-{self.reference}/point-{k}.npz')); S,esc=self.scattering_matrix(c)
        X=self.I[k]/self.scale; scat=(S@X.ravel()).reshape(X.shape); escaped=np.einsum('knqf,nqf->kn',esc,X)
        em=c['emit']/self.scale; bound=em-c['loss']*X; coll=bound+scat
        derivative={}; dbound={}; desc={}
        for param in ['dr','dt','dy']:
            ratio=np.divide(c[param+'_sc'],c['sc'],out=np.zeros(self.n),where=c['sc']>0)
            dbound[param]=c[param+'_emit']/self.scale-c[param+'_loss']*X
            derivative[param]=dbound[param]+ratio[:,None,None]*scat
            desc[param]=ratio[None,:]*escaped
        et=c['dt_energy']; ey=c['dy_energy']; neutral=c['neutral']; active=et>0
        mapT=np.zeros((self.n,2)); mapY=np.zeros_like(mapT)
        mapT[:,0]=np.divide(self.eu,et,out=np.zeros(self.n),where=active)
        mapY[:,1]=np.divide(self.nu,neutral,out=np.zeros(self.n),where=neutral>0)
        mapT[:,1]=-np.divide(ey,et,out=np.zeros(self.n),where=active)*mapY[:,1]
        B=derivative['dt'][...,None]*mapT[:,None,None,:]+derivative['dy'][...,None]*mapY[:,None,None,:]
        Bb=dbound['dt'][...,None]*mapT[:,None,None,:]+dbound['dy'][...,None]*mapY[:,None,None,:]
        Be=desc['dt'][...,None]*mapT[None,:,:]+desc['dy'][...,None]*mapY[None,:,:]
        # These fields multiply the actual interpolated volume/lapse changes.
        vol=em-derivative['dr']-c['adiabatic_logT'][:,None,None]*derivative['dt']
        volb=em-dbound['dr']-c['adiabatic_logT'][:,None,None]*dbound['dt']
        vole=-desc['dr']-c['adiabatic_logT'][None,:]*desc['dt']
        row=dict(S=S,esc=esc,loss=c['loss'],sc=c['sc'],B=B,Bb=Bb,Be=Be,vol=vol,volb=volb,vole=vole,coll=coll,bound=bound,escaped=escaped,
            pressure_map=c['dt_p'][:,None]*mapT+c['dy_p'][:,None]*mapY)
        self.cache[k]=row
        for old in list(self.cache):
            if old not in [k,k-1,k+1]:del self.cache[old]
        return row

    def local(self,t):
        k=max(0,min(np.searchsorted(self.t,t,side='right')-1,15)); f=(t-self.t[k])/(self.t[k+1]-self.t[k]); a=self.point(k); b=self.point(k+1)
        c={key:(1-f)*a[key]+f*b[key] for key in ['loss','sc','B','Bb','Be','pressure_map']}
        c['S']=(1-f)*a['S']+f*b['S']; c['esc']=(1-f)*a['esc']+f*b['esc']
        volume=((1-f)*(3*self.g['delta_u'][k]+self.g['delta_lambda'][k])+f*(3*self.g['delta_u'][k+1]+self.g['delta_lambda'][k+1]))/AMPLITUDE
        lapse=((1-f)*self.g['delta_log_lapse'][k]+f*self.g['delta_log_lapse'][k+1])/AMPLITUDE
        c['q']=sum(w*(volume[:,None,None]*r['vol']+lapse[:,None,None]*r['coll']) for w,r in [(1-f,a),(f,b)])
        c['qb']=sum(w*(volume[:,None,None]*r['volb']+lapse[:,None,None]*r['bound']) for w,r in [(1-f,a),(f,b)])
        c['qe']=sum(w*(volume[None,:]*r['vole']+lapse[None,:]*r['escaped']) for w,r in [(1-f,a),(f,b)])
        return c

    def gas(self,p,b,e):
        return np.stack([-(np.sum(p*self.Eweight,axis=(1,2))+e[1])/self.eu,np.sum(b*self.Nweight,axis=(1,2))/self.nu],axis=1)

    def pack(self,x,g):return np.r_[x.ravel(),g.ravel()]
    def unpack(self,x):return x[:self.size].reshape(self.n,self.q,self.nf),x[self.size:].reshape(self.n,2)

    def collision(self,c,x,g,source=False):
        bound=-c['loss']*x+np.einsum('nqfj,nj->nqf',c['Bb'],g)
        p=-c['loss']*x+(c['S']@x.ravel()).reshape(x.shape)+np.einsum('nqfj,nj->nqf',c['B'],g)
        e=np.einsum('knqf,nqf->kn',c['esc'],x)+np.einsum('knj,nj->kn',c['Be'],g)
        if source:p+=c['q'];bound+=c['qb'];e+=c['qe']
        return p,self.gas(p,bound,e),e,bound

    def inverse(self,c,h):
        # Static Thomson has rank two in angle; retain exact directional
        # absorption and the two material Schur variables in this inverse.
        D=1+h*(c['loss']+c['sc'][:,None,None]); U=h*c['sc'][:,None,None]*np.stack([np.ones(self.q),.5*self.model.bulk.P2],axis=1)[None]
        V=np.stack([self.w,self.w*self.model.bulk.P2]); DU=U[:,:,None,:]/D[:,:,:,None]
        M=np.eye(2)[None,None]-np.einsum('jq,nqfk->nfjk',V,DU)
        def photon(rhs):
            z=rhs/D[:,:,:,None]; v=np.einsum('jq,nqfl->nfjl',V,z)
            return z+np.einsum('nqfj,nfjl->nqfl',DU,np.linalg.solve(M,v))
        PB=photon(c['B'])
        def G(x):
            bound=-c['loss'][:,:,:,None]*x
            return np.stack([-np.sum(bound*self.Eweight[:,:,:,None],axis=(1,2))/self.eu[:,None],np.sum(bound*self.Nweight[:,:,:,None],axis=(1,2))/self.nu[:,None]],axis=1)
        DB=np.stack([self.gas(c['B'][...,j],c['Bb'][...,j],c['Be'][...,j]) for j in range(2)],axis=2)
        Mgas=np.eye(2)[None]-h*DB-h*h*G(PB)
        def apply(rhs):
            x,g=self.unpack(rhs); z=photon(x[...,None])[...,0]
            y=np.linalg.solve(Mgas,(g+h*G(z[...,None])[...,0])[...,None])[...,0]
            return self.pack(z+h*np.einsum('nqfj,nj->nqf',PB,y),y)
        return apply

    def step_collision(self,c,x,g,h):
        inverse=self.inverse(c,h); iterations=[]
        def mat(v):
            xx,gg=self.unpack(v); p,q,*_=self.collision(c,xx,gg); return v-h*self.pack(p,q)
        rhs=self.pack(x+h*c['q'],g+h*self.gas(c['q'],c['qb'],c['qe']))
        op=LinearOperator((len(rhs),len(rhs)),mat,dtype=float); pre=LinearOperator(op.shape,inverse,dtype=float)
        sol,info=gmres(op,rhs,M=pre,rtol=1e-10,atol=0.,restart=15,maxiter=4,callback=iterations.append,callback_type='pr_norm')
        relative=float(np.linalg.norm(mat(sol)-rhs)/max(np.linalg.norm(rhs),1e-300)); self.max_residual=max(self.max_residual,relative); self.max_iterations=max(self.max_iterations,len(iterations))
        assert info==0 and relative<1e-10,('Collision solve',info,relative,len(iterations))
        return self.unpack(sol)

    def run(self,steps,label,limit=None):
        assert not (OUT/f'{label}.npz').exists(); start=time.monotonic(); h=self.t[-1]/steps; gamma=prior.GAMMA
        lu=splu(sparse.eye(self.n*self.q,format='csc')-gamma*h/2*self.A)
        x=np.zeros_like(self.I[0]); g=np.zeros((self.n,2)); ledger=np.zeros(2); escape=np.zeros(3); impulse=np.zeros(self.n); error=0.; records=[]; times=[]; moment=0.
        def record(t):
            times.append(t); records.append(np.stack([np.sum(x*self.Eweight,axis=(1,2)),g[:,0]*self.eu,g[:,1]*self.nu,impulse,np.sum(abs(x)*self.Eweight,axis=(1,2))])*AMPLITUDE)
        record(0.)
        def streaming(x,t,dt):
            s1,l1,e1=self.source(t+gamma*dt); s1=s1/(self.scale*AMPLITUDE); l1=l1/AMPLITUDE
            z=lu.solve((x+gamma*dt*s1).reshape(self.n*self.q,self.nf)).reshape(x.shape)
            first=(self.A@z.reshape(self.n*self.q,self.nf)).reshape(x.shape)+s1
            s2,l2,e2=self.source(t+dt); s2=s2/(self.scale*AMPLITUDE); l2=l2/AMPLITUDE
            y=lu.solve((x+(1-gamma)*dt*first+gamma*dt*s2).reshape(self.n*self.q,self.nf)).reshape(x.shape)
            port=dt*((1-gamma)*(self.port(z*self.scale)+l1[3:])+gamma*(self.port(y*self.scale)+l2[3:]))
            return y,port,dt*((1-gamma)*l1[:3]+gamma*l2[:3]),max(e1,e2)
        count=steps if limit is None else limit
        for k in range(count):
            t=k*h; x,port,ghost,mm=streaming(x,t,h/2); ledger+=port; escape+=ghost; moment=max(moment,mm)
            c=self.local(t+h/2)
            z,gg=self.step_collision(c,x,g,gamma*h); p1,g1,e1,b1=self.collision(c,z,gg,True)
            y,gy=self.step_collision(c,x+(1-gamma)*h*p1,g+(1-gamma)*h*g1,gamma*h); p2,g2,e2,b2=self.collision(c,y,gy,True)
            e=h*((1-gamma)*e1+gamma*e2); esc=e.sum(1); escape+=esc
            neutral=h*np.sum(((1-gamma)*b1+gamma*b2)*self.Nweight)
            ledger+=np.array([neutral-esc[0],-esc[1]])
            impulse-=h*np.sum(((1-gamma)*p1+gamma*p2)*self.Eweight*self.mu[None,:,None],axis=(1,2))+e[2]
            x,g=y,gy; x,port,ghost,mm=streaming(x,t+h/2,h/2); ledger+=port; escape+=ghost; moment=max(moment,mm)
            photon=self.moments(x*self.scale); total=np.array([photon[0],photon[1]+np.sum(g[:,0]*self.eu)])
            denominator=max(abs(total[1]),abs(ledger[1]),np.sum(abs(x)*self.Eweight),np.sum(abs(g[:,0])*self.eu),1.)
            error=max(error,float(abs(total[1]-ledger[1])/denominator))
            if (k+1)%max(1,steps//16)==0 or k+1==count:record((k+1)*h)
        np.savez_compressed(OUT/f'{label}.npz',t=times,moments=records,delta_packet_scaled_occupation=x*self.scale*AMPLITUDE,delta_material=g*AMPLITUDE,material_energy_units=self.eu,material_neutral_units=self.nu,ledger=ledger*AMPLITUDE,escape=escape*AMPLITUDE,radius_E=self.r)
        row=dict(classification='Counterexample candidate',reference=self.reference,steps=steps,completed_steps=count,seconds=time.monotonic()-start,energy_balance_relative=error,linear_relative=self.max_residual,max_Krylov_iterations=self.max_iterations,frequency_moment_relative=moment,endpoint_photon_reference_energy_erg=float(np.sum(x*self.Eweight)*AMPLITUDE),endpoint_material_reference_energy_erg=float(np.sum(g[:,0]*self.eu)*AMPLITUDE),endpoint_material_energy_L1_erg=float(np.sum(abs(g[:,0])*self.eu)*AMPLITUDE),material_momentum_impulse_saved=True,additional_material_motion_evolved=False,full_GR_feedback=False,final_charge_solved=False)
        row['passed']=error<1e-8 and self.max_residual<1e-10
        write(OUT/f'{label}.json',row); print(json.dumps(row),flush=True); return row


def pilot():
    assert json.loads((OUT/'bank-result.json').read_text())['passed']; assert not (OUT/'pilot.json').exists()
    start=time.monotonic(); signal.signal(signal.SIGALRM,flow.old.optical.timeout); signal.alarm(45)
    rows=[CoupledResponse(128).run(n,f'pilot-{n}',2) for n in [64,128]]
    estimate=rows[0]['seconds']/2*64+rows[1]['seconds']/2*256
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',rows=rows,forecast_seconds=estimate,upper_seconds=2*estimate,eligible=all(r['passed'] for r in rows) and 2*estimate<240,seconds=time.monotonic()-start));signal.alarm(0)


if __name__=='__main__': globals()[sys.argv[1]]()
