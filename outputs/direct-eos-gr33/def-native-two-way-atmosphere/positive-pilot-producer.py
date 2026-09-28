"""Counterexample candidate: shared photons and an actually moving atmosphere.

The deep material/photon equations remain unsplit. Atmospheric transport is
eliminated in angular direction order; local moving collisions and fluid flow
are subcycled with their shared conservative energy and species transfers.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
import def_native_unsplit_photons as deep
import def_native_reactive_interface as interface
import def_native_atmosphere_spectrum as optical

OUT=optical.OUT;write=optical.write;sha=optical.sha;C=interface.C;END=deep.prior.prior.END

# Keep the validated physical fluid reconstruction, gravity, species donor
# flux and primitive recovery. Remove the prescribed photon bath at its owner.
hydro_source=interface.rhs_source.split('    area=(m.r/m.RJ)**2;Fbase=')[0]+'''    self.eos.y=y
    ledger=np.array([C*(flux[0,0]-flux[0,-1]),C*(flux[2,0]-flux[2,-1]),C*(flux[3,0]-flux[3,-1])])
    dt=.35*m.dx/np.max(C*m.a/m.B*(abs(v)+cs+1e-100))
    self.max_speed=max(self.max_speed,float(max(abs(v))))
    return rate,ledger,dt,V,kap
'''
ns=dict(interface.ns,OUT=OUT);exec(compile(hydro_source,__file__,'exec'),ns)


class Flow(interface.PhysicalFlow):
    hydro=ns['rhs']


def transport(W,area,mu,edges_mu):
    n=len(W);q=len(mu);dm=np.diff(edges_mu);ids=np.arange(n*q).reshape(n,q);rr=[];cc=[];vv=[]
    def add(i,j,v):rr.append(int(i));cc.append(int(j));vv.append(float(v))
    for k,u in enumerate(mu):
        for face in range(n+1):
            f=C*area[face]*u
            if face==0:
                if u<0:add(ids[0,k],ids[0,k],f/W[0])
            elif face==n:
                if u>0:add(ids[-1,k],ids[-1,k],-f/W[-1])
            else:
                j=face-1 if u>0 else face
                add(ids[face-1,k],ids[j,k],-f/W[face-1]);add(ids[face,k],ids[j,k],f/W[face])
    for i in range(n):
        for k in range(1,q):
            f=C*(area[i+1]-area[i])/2*(1-edges_mu[k]**2)
            assert f>=0
            add(ids[i,k-1],ids[i,k-1],-f/(W[i]*dm[k-1]));add(ids[i,k],ids[i,k-1],f/(W[i]*dm[k]))
    return sparse.coo_matrix((vv,(rr,cc)),shape=(n*q,n*q)).tocsc()


class Coupled:
    def __init__(self,cells=448,angles=8):
        self.flow=f=Flow(cells);self.m=m=f.base;self.bulk=b=deep.Model(angles);self.spectrum=optical.Spectrum()
        self.n=f.n;self.q=angles;self.mu=b.mu;self.w=b.w;self.freq=b.m
        assert abs(m.RJ-b.d['edges'][-1])<1e-4
        old_width=np.diff(b.d['edges'])[-1];new_edge=m.rf[0];ratio=(new_edge-b.d['edges'][-2])/old_width;assert 0<ratio<1
        b.d=dict(b.d);b.d['edges']=b.d['edges'].copy();b.d['edges'][-1]=new_edge
        b.W[-1]*=ratio;b.volume[-1]*=ratio;b.gas_weight[-1]*=ratio;b.photon_energy_weight[-1]*=ratio
        b.area[-1]=new_edge**2/m.af[0]**2;b.A=transport(b.W,b.area,b.mu,b.edges_mu);b.stream=b.A.tocoo()
        self.W=m.RJ**2*m.vol/m.a**3;self.area=m.rf*m.rf/(m.af*m.af);self.A=transport(self.W,self.area,self.mu,b.edges_mu)
        self.ids=np.arange(self.n*self.q).reshape(self.n,self.q);self.neg=self.ids[:,self.mu<0].ravel();self.pos=self.ids[:,self.mu>0].ravel()
        assert self.A[self.neg][:,self.pos].nnz==0,'Angular elimination order'
        self.pm=self.A[self.pos][:,self.neg];self.mm=self.A[self.neg][:,self.neg];self.pp=self.A[self.pos][:,self.pos];self.lu={}
        self.energy_weight=4*np.pi*self.W[:,None,None]*self.w[None,:,None]*b.d['num']*b.d['Einf']
        self.gas_scale=4*np.pi*m.RJ*m.RJ*f.eos.rho0*C*C
        self.number=b.d['num'];self.E=b.d['Einf'];self.S=np.ones((angles,1))*self.w[None,:]+.5*b.P2[:,None]*(b.P2*self.w)[None,:]
        bg=deep.prior.prior.Background();inside=m.r<=m.RJ;I=np.zeros((self.n,self.q,self.freq))
        z=bg.sample(m.r[inside]);J=1/np.expm1(self.E[None,:]/(z['a']*optical.ex.K*z['T'])[:,None])
        Erad=(J*self.number*self.E).sum(1)/z['a']**4;frac=bg.luminosity(z)/(4*np.pi*z['r']**2*z['a']**2*C*Erad)
        I[inside]=J[:,None,:]*b.shape(frac)[:,:,None]
        cut=np.sqrt(np.maximum(0,1-(m.RJ/m.r[~inside])**2*(m.a[~inside]/m.a0)**2))
        angular=np.maximum(0,b.edges_mu[1:][None,:]-np.maximum(b.edges_mu[:-1][None,:],cut[:,None]))/(2/self.q)
        js=1/np.expm1(self.E/(m.a0*self.spectrum.native.fan.T*optical.ex.K))
        I[~inside]=angular[:,:,None]*js
        self.initial_I=I;self.maximum_frame_speed=0.;self.maximum_source_optical_step=0.;self.scatter_number_error=0.

    def stage(self,xb,xa,u,y,h,theta,eta):
        b=self.bulk;q=self.q;n=self.n
        if h not in self.lu:self.lu[h]=(splu(sparse.eye(len(self.neg),format='csc')-h*self.mm),splu(sparse.eye(len(self.pos),format='csc')-h*self.pp))
        lm,lp=self.lu[h];flat=xa.reshape(2,n*q,self.freq)
        def solve(lu,rhs):return lu.solve(rhs.transpose(1,0,2).reshape(rhs.shape[1],2*self.freq)).reshape(rhs.shape[1],2,self.freq).transpose(1,0,2)
        minus=solve(lm,flat[:,self.neg]);at=np.zeros_like(flat);at[:,self.neg]=minus
        inward=minus[:,:q//2].sum(0)*b.scale
        b.boundary[-1,self.mu<0]=-C*b.area[-1]*self.mu[self.mu<0,None]/b.W[-1]*inward
        xb,uu,tt,yy,it,err=b.implicit(xb,u,y,h,theta,eta)
        rhs=flat[:,self.pos]+h*np.array([self.pm@v for v in minus])
        rhs[0,:q//2]+=h*C*self.area[0]*self.mu[self.mu>0,None]/self.W[0]*xb[-1,self.mu>0]
        plus=solve(lp,rhs);at[:,self.pos]=plus;assert at.sum(0).min()>=0,('Atmospheric photon positivity',float(at.sum(0).min()))
        return xb,at.reshape(xa.shape),uu,tt,yy,it,err

    def port(self,xb,xa):
        b=self.bulk;inner=np.where((self.mu>0)[:,None],b.incoming,xb[0]*b.scale)
        full=xa.sum(0);outer=np.where((self.mu>0)[:,None],full[-1]*b.scale,0.)
        Hi=self.w*self.mu@inner;Ho=self.w*self.mu@outer
        delta=xa[1]*b.scale;di=np.where((self.mu<0)[:,None],delta[0],0.);do=np.where((self.mu>0)[:,None],delta[-1],0.)
        change=float(4*np.pi*C*((self.area[0]*(self.w*self.mu@di)-self.area[-1]*(self.w*self.mu@do))*(self.number*self.E)).sum())
        return np.array([float(4*np.pi*C*((b.area[0]*Hi-self.area[-1]*Ho)*(self.number*self.E)).sum()),change])

    def radiate(self,xb,I,u,theta,eta,h):
        xa=I/self.bulk.scale
        bb,aa,uu,tt,yy,it,err=self.stage(xb,xa,u,eta,h,theta,eta)
        return bb,aa*self.bulk.scale,uu,tt,yy,h*self.port(bb,aa),it,err

    def scattering(self,I,rho,v,kap,a):
        """Elastic Thomson in moving gas; positive photon packet remapping.

        The declared frequency box has absorbing exits. Record their number,
        energy and momentum separately; never deposit escaped photons as heat.
        """
        n=len(rho);q=self.q;m=self.freq;W=1/np.sqrt(1-v*v)
        edges=(self.bulk.edges_mu[None,:]-v[:,None])/(1-v[:,None]*self.bulk.edges_mu[None,:]);dw=np.diff(edges)/2
        p2=(edges[:,:-1]**2+edges[:,:-1]*edges[:,1:]+edges[:,1:]**2-1)/2
        num=I*self.w[None,:,None]*self.number[None,None,:]/a[:,None,None]**3
        rate=a[:,None]*C*rho[:,None]*kap[:,None]*W[:,None]*(1-v[:,None]*self.mu[None,:])
        loss=num*rate[:,:,None];net=-loss.copy();escape=np.zeros((n,3));rows=np.arange(n)[:,None]
        for j,muin in enumerate(self.mu):
            for k,muout in enumerate(self.mu):
                probability=dw[:,k]*(1+.5*p2[:,j]*p2[:,k]);amount=loss[:,j,:]*probability[:,None]
                dest=self.E[None,:]*(1-v[:,None]*muin)/(1-v[:,None]*muout)
                high=np.searchsorted(self.E,dest);outside=(high==0)|(high==m)
                # An unchanged endpoint belongs to its represented bin.
                exact=(dest==self.E[None,:]);outside&=~exact
                low=np.clip(high-1,0,m-2);upper=low+1;fraction=(dest-self.E[low])/(self.E[upper]-self.E[low])
                for index,weight in [(low,1-fraction),(upper,fraction)]:
                    values=np.where(outside,0.,amount*weight);np.add.at(net[:,k,:],(rows,index),values)
                escaped=np.where(outside,amount,0.);escape[:,0]+=escaped.sum(1);escape[:,1]+=(escaped*dest).sum(1);escape[:,2]+=(escaped*dest*muout).sum(1)
        error=np.max(abs(net.sum((1,2))+escape[:,0])/np.maximum(loss.sum((1,2)),1e-250));self.scatter_number_error=max(self.scatter_number_error,float(error))
        return net*a[:,None,None]**3/(self.w[None,:,None]*self.number[None,None,:]),escape

    def local_rhs(self,U,I,t):
        f=self.flow;m=self.m;gas,ledger,dt,V,kap=f.hydro(U,t);rho,v,lt,y=V;active=rho>=f.eos.floor;idx=np.flatnonzero(active)
        change=np.zeros_like(I);g=np.zeros_like(U);escape=np.zeros((self.n,3));stiffness=0.
        if len(idx):
            a=m.a[idx];rr=rho[idx]*f.eos.rho0;vv=v[idx];yy=y[idx];tt=lt[idx];W=1/np.sqrt(1-vv*vv);D=W[:,None]*(1-vv[:,None]*self.mu[None,:])
            energy=self.E[None,None,:]*D[:,:,None]/a[:,None,None]
            ab,em=self.spectrum.coefficients(rr,tt,yy,energy);field=I[idx]
            factor=a[:,None,None]*C*rr[:,None,None]*f.eos.nH*D[:,:,None]
            absorbed=factor*ab*field;emitted=factor*em*(1+field);bound=emitted-absorbed
            scatter,esc=self.scattering(field,rr,vv,kap[idx],a);total=bound+scatter;change[idx]=total;escape[idx]=esc
            number=np.einsum('iqf,q,f->i',bound,self.w,self.number)/a**3
            force=np.einsum('iqf,q,f->i',total,self.w*self.mu,self.number*self.E)/a**4+esc[:,2]/a
            power=np.einsum('iqf,q,f->i',total,self.w,self.number*self.E)/a**3+esc[:,1]
            g[1,idx]=-force/(f.eos.rho0*C*C);g[2,idx]=-power/(f.eos.rho0*C*C);g[3,idx]=number/(f.eos.rho0*f.eos.nH)
            destroy=np.einsum('iqf,q,f->i',absorbed,self.w,self.number)/a**3/(f.eos.rho0*f.eos.nH)
            create=np.einsum('iqf,q,f->i',emitted,self.w,self.number)/a**3/(f.eos.rho0*f.eos.nH)
            stiffness=float(max(destroy/U[3,idx]+create/(U[0,idx]-U[3,idx])))
            dt=min(dt,.2/max(stiffness,1e-100));self.maximum_frame_speed=max(self.maximum_frame_speed,float(max(abs(vv))))
        gas+=g
        for value,der in [(U[0],gas[0]),(U[3],gas[3]),(U[0]-U[3],gas[0]-gas[3])]:
            losing=(der<0)&active
            if losing.any():dt=min(dt,.8*float(min(value[losing]/(-der[losing]))))
        losing=change<0
        if losing.any():dt=min(dt,.8*float(np.min(I[losing]/(-change[losing]))))
        transfer=float(np.sum(g[2]*m.vol));species=float(np.sum(g[3]*m.vol));escaped=float(4*np.pi*m.RJ**2*(escape[:,1]@m.vol))
        return gas,change,np.r_[ledger,transfer,species,escaped],dt

    def local(self,U,I,t,duration):
        elapsed=0.;ledger=np.zeros(6);discard=np.zeros(4);steps=0;m=self.m
        while elapsed<duration:
            k,p,l,dt=self.local_rhs(U,I.sum(0),t+elapsed);dt=min(dt,duration-elapsed);assert dt>0
            for _ in range(12):
                trial=U+dt*k;phot=I[1]+dt*p;assert (I[0]+phot).min()>=0
                k2,p2,l2,limit=self.local_rhs(trial,I[0]+phot,t+elapsed+dt)
                if dt<=limit*(1+1e-12):break
                dt=min(dt/2,limit)
            else:raise AssertionError('Local SSP stage time cap')
            U=(U+trial+dt*k2)/2;I[1]=(I[1]+phot+dt*p2)/2;assert I.sum(0).min()>=0
            tiny=U[0]<self.flow.eos.floor;discard+=np.sum(U[:,tiny]*m.vol[tiny],axis=1);U[:,tiny]=0
            ledger+=dt*(l+l2)/2;elapsed+=dt;steps+=1;assert steps<20000
        return U,I,ledger,discard,steps

    def run(self,steps,label,stop_after=None):
        assert not (OUT/(label+'.npz')).exists();start=time.monotonic();b=self.bulk;f=self.flow;m=self.m;h=END/steps;count=steps if stop_after is None else stop_after
        U=f.initial.copy();I=np.stack([self.initial_I,np.zeros_like(self.initial_I)]);xb=b.initial/b.scale;u=b.u0.copy();theta=np.zeros(b.n);eta=np.zeros(b.n)
        ledger=np.zeros(6);discard=np.zeros(4);boundary=np.zeros(2);balance=0.;gas_balance=0.;local_balance=0.;iterations=0;residual=0.;substeps=0;history=[];snapshots=[];failure=None
        p0=b.eos.gas(theta,eta)[0];f.eos.y=f.eos.y0;iu=f.eos(U[0],f.initial_temperature)[1];ip=f.eos(U[0],f.initial_temperature)[0]
        def record(t):
            rho,v,lt,y=f.primitive(U);p,uu,*_=f.eos(rho,lt)
            trace=-f.eos.cx*U[0]*v*v/(1+np.sqrt(1-v*v))+rho*uu-3*p-(f.initial[0]*iu-3*ip)
            bulk_trace=b.d['rho']*(u-b.u0)-3*(b.eos.gas(theta,eta)[0]-p0)
            history.append(dict(t=t,bulk_trace=float(bulk_trace@(b.volume*b.d['a'])),atmosphere_trace=float(trace@m.vol*self.gas_scale),
                outside_mass=float(np.sum(U[0,m.x>=0]*m.vol[m.x>=0])*self.gas_scale/(C*C)),maximum_speed=float(max(abs(v)))))
        record(0.)
        try:
            for j in range(count):
                U,I,ll,dd,ss=self.local(U,I,j*h,h/2);ledger+=ll;discard+=dd;substeps+=ss
                xb,I,u,theta,eta,port,it,err=self.radiate(xb,I,u,theta,eta,h);boundary+=port;iterations=max(iterations,it);residual=max(residual,err)
                U,I,ll,dd,ss=self.local(U,I,j*h+h/2,h/2);ledger+=ll;discard+=dd;substeps+=ss
                record((j+1)*h)
                photon=float(np.sum((xb*b.scale-b.initial)*b.photon_energy_weight)+np.sum((I.sum(0)-self.initial_I)*self.energy_weight))
                material=float(b.gas_weight@(u-b.u0)+(np.sum((U[2]-f.initial[2])*m.vol)+discard[2])*self.gas_scale)
                expected=boundary[0]+ledger[1]*self.gas_scale-ledger[5];balance=max(balance,abs(photon+material-expected))
                gas_res=(np.sum((U[2]-f.initial[2])*m.vol)+discard[2]-ledger[1]-ledger[3])*self.gas_scale;gas_balance=max(gas_balance,abs(gas_res))
                atmosphere_res=float(np.sum(I[1]*self.energy_weight))+(np.sum((U[2]-f.initial[2])*m.vol)+discard[2]-ledger[1])*self.gas_scale-boundary[1]+ledger[5]
                local_balance=max(local_balance,abs(atmosphere_res))
                if (j+1)%max(1,steps//16)==0 or j+1==count:snapshots.append(dict(U=U.copy(),I=I.copy(),bulk_I=xb*b.scale,theta=theta.copy(),eta=eta.copy(),t=(j+1)*h))
        except Exception as exc:failure=repr(exc)
        arrays={k:np.array([row[k] for row in history]) for k in history[0]};response=max(abs(arrays['bulk_trace']).max(),abs(arrays['atmosphere_trace']).max(),abs(boundary[0]),1.)
        gas_response=max(abs(arrays['atmosphere_trace']).max(),abs(ledger[3]*self.gas_scale),1.)
        mass=abs(np.sum((U[0]-f.initial[0])*m.vol)+discard[0]-ledger[0])/np.sum(f.initial[0]*m.vol)
        species=abs(np.sum((U[3]-f.initial[3])*m.vol)+discard[3]-ledger[2]-ledger[4])/np.sum(f.initial[3]*m.vol)
        np.savez_compressed(OUT/(label+'.npz'),**arrays,U=U,I=I,bulk_I=xb*b.scale,theta=theta,eta=eta,initial_U=f.initial,initial_I=self.initial_I,
            ledger=ledger,discard=discard,boundary_energy=boundary,**{f'snapshot_{k}':np.array([row[k] for row in snapshots]) for k in ['U','I','bulk_I','theta','eta','t']})
        result=dict(classification='Counterexample candidate',passed=bool(failure is None and balance/response<1e-8 and gas_balance/gas_response<1e-8 and local_balance/gas_response<1e-8 and mass<1e-10 and species<1e-9 and residual<1e-9),
            failure=failure,atmosphere_cells=f.atmosphere_cells,actual_fluid_cells=f.n,angles=self.q,steps=steps,completed_steps=len(history)-1,local_SSP_steps=substeps,
            seconds=time.monotonic()-start,total_energy_relative=balance/response,gas_energy_relative=gas_balance/gas_response,atmosphere_photon_energy_relative=local_balance/gas_response,baryon_relative=float(mass),species_relative=float(species),
            maximum_deep_residual=residual,maximum_deep_Newton=iterations,maximum_velocity_over_c=self.maximum_frame_speed,scatter_number_relative=self.scatter_number_error,
            spectral_escape_energy_erg=float(ledger[5]),atmosphere_trace_erg=float(arrays['atmosphere_trace'][-1]),bulk_trace_erg=float(arrays['bulk_trace'][-1]),
            coupled_photons_and_moving_atmosphere=failure is None,space_time_comparison_passed=False,full_GR_scalar_feedback=False,final_charge_solved=False,source_sha256=sha(__file__))
        write(OUT/(label+'.json'),result);print(json.dumps(result),flush=True);return result


def pilot():
    assert json.loads((OUT/'repaired-spectrum-controls.json').read_text())['passed'];signal.signal(signal.SIGALRM,optical.timeout);signal.alarm(45)
    (OUT/'hydrodynamic-source.py').write_text(hydro_source)
    row=Coupled(448,8).run(128,'pilot-448',1)
    write(OUT/'coupled-pilot-budget.json',dict(pilot=row,forecast_seconds=row['seconds']*(64+128+2*128),
        assumption='A fine fluid/photon path is assumed2x measured coarse per macrostep; later recovery/collision cost remains unmeasured.',eligible=bool(row['passed'] and row['seconds']*(64+128+2*128)*1.5<600)))
    signal.alarm(0)


def positive_pilot():
    assert not json.loads((OUT/'pilot-448.json').read_text())['passed']
    write(OUT/'positive-time-plan.json',dict(classification='Counterexample candidate',
        failure='The second SDIRK stage of the first actual joint photon/fluid step produced a negative atmospheric scaled occupation-4.1742e-9. The atmosphere light-crossing time is much shorter than the deep128step clock. Preserve this failure and all source code.',
        repair='Use backward Euler for the simultaneous deep reaction and whole photon transport. Its frozen positive transfer generator has a nonnegative resolvent. This is a first-order temporal method; keep the same64/128 time acceptance gate rather than claiming inherited second-order accuracy.',
        unchanged='Same actual deep/fluid domains,EOS,frequency grid,angular bins,initial occupations,paired Lorentz collision transfers and original positivity/error gates. No clipping or extra time paths.',
        source_sha256=sha(__file__)))
    signal.signal(signal.SIGALRM,optical.timeout);signal.alarm(45)
    row=Coupled(448,8).run(128,'positive-pilot-448',1)
    forecast=row['seconds']*(64+128+2*128)
    write(OUT/'positive-pilot-budget.json',dict(pilot=row,forecast_seconds=forecast,upper_seconds=forecast*1.5,
        eligible=bool(row['passed'] and forecast*1.5<600),assumption='Fine fluid/photon per-macrostep cost assumed2x coarse, then1.5x overall margin; actual later and fine behavior remains unmeasured.'))
    signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
