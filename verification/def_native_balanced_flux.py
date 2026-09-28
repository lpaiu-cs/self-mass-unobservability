"""Resolve the constant-current cancellation in the native coupled GR loop.

Counterexample candidate: the same frozen C1 tangent equations, in deviation
coordinates. The physical surface history is retained for moving-ray coupling.
"""
from pathlib import Path
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import bmat, diags
from scipy.sparse.linalg import splu
from scipy.interpolate import CubicHermiteSpline
import def_native_conservative_source as prior
from def_native_source_order import condense_factor

old=prior.old
OUT=old.OUT.parent/'def-native-balanced-flux'
write=old.write
source=prior.prior.source
anchor='op=diags(bg.tc*pref)@diff@diags(self.theta)'
assert source.count(anchor)==1
source=source.replace(anchor,anchor+'\n    self.coefficient_data=data;self.cov=cov;self.thermal_operator=op')
namespace=dict(prior.assemble.__globals__)
exec(compile(source,__file__,'exec'),namespace)
assemble=namespace['__init__']


class Model(prior.Model):
    def __init__(self):
        started=time.monotonic();previous=prior.assemble;prior.assemble=assemble
        try:super().__init__(4,lumped=False)
        finally:prior.assemble=previous
        self.f0=(self.bg.tc*self.L0).astype(np.longdouble)
        rr,loss,J=self.source_values(self.points,self.f0)
        src=np.array([rr,loss,J]).T;data=self.coefficient_data.astype(np.longdouble)
        B=data[:,4:8].reshape(-1,2,2)
        g=np.einsum('nij,nj->ni',data[:,16:22].reshape(-1,2,3),src)
        h=np.einsum('nij,nj->ni',data[:,22:28].reshape(-1,2,3),src)
        weights=self.weights.astype(np.longdouble)
        self.F0=sum(self.cov[i].astype(np.longdouble).T@(weights*B[:,i,j]*g[:,j]) for i in range(2) for j in range(2))
        self.F0-=sum(self.V[i].astype(np.longdouble).T@(weights*h[:,i]) for i in range(2))
        p=self.bg.sample(self.native);_,loss,self.J0=self.source_values(p,self.f0)
        b=1-2*p['m']/self.native;geo=old.G*self.bg.R**2/old.C**4
        self.TE0=-self.thermo[:,4]/(self.native*b)*self.J0-loss/(self.raw[:,0]*geo*self.thermo[:,3])
        self.LE0=np.r_[self.thermal_operator@self.TE0,4*self.f0[-1]*self.TE0[-1]]
        self.surfaceV=self.evaluation(np.array([1.]))[0]
        ps=self.bg.sample(np.array([1.]));bs=1-2*ps['m'][0]
        env=np.load(old.prior.OUT/'final-envelope.npz')
        Pr=float(env['Prad'][-1])*geo
        traction=-16*np.pi*np.exp(-8*ps['phi'][0]**2)*ps['N'][0]*np.sqrt(bs)*Pr/bs
        self.F0+=np.asarray(self.surfaceV[0].T@(traction*self.TE0[-1:])).ravel()
        self.Jmap=self.sources(p)[2]
        del self.coefficient_data,self.cov,self.V,self.D
        self.setup_seconds=time.monotonic()-started

    def source_values(self,p,face):
        # Difference the face inventory BEFORE interpolation. Constant through
        # current then has exactly zero local debit outside its central source.
        r=p['r'];ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,self.n-1)
        face=np.asarray(face,np.longdouble);volume=self.volumes.astype(np.longdouble)
        inc=np.diff(np.r_[np.longdouble(0),face]);sec=inc/volume
        slope=np.r_[sec[0],(volume[1:]*sec[:-1]+volume[:-1]*sec[1:])/(volume[1:]+volume[:-1]),sec[-1]]
        base=old.Model.face_map(self,r);x=np.asarray(base[np.arange(len(r)),ids]).ravel().astype(np.longdouble)
        left=np.r_[np.longdouble(0),face[:-1]][ids]
        dl=slope[ids]-sec[ids];dr=slope[ids+1]-sec[ids]
        enclosed=left+x*inc[ids]+volume[ids]*(x*(1-x)**2*dl+x*x*(x-1)*dr)
        density=sec[ids]+(3*x*x-4*x+1)*dl+(3*x*x-2*x)*dr
        inside=r<=1;geo=old.G*self.bg.R**2/old.C**4;A4=np.exp(-8*p['phi']**2)
        loss=geo/A4*density*inside;ratio=np.zeros_like(r)
        np.divide(np.interp(r,self.native,self.thermo[:,4]),p['gamma']*p['p'],out=ratio,where=inside)
        b=1-2*p['m']/np.maximum(r,1e-100)
        J=-old.G/(old.C**4*self.bg.R)*np.sqrt(b)/p['N']*enclosed*inside
        return ratio*loss,loss,J

    def photon_force(self,t):
        mu,weights,delay=self.rays;at=t-delay;u=np.maximum(at,0)
        de=np.zeros_like(u);df=np.zeros_like(u)
        if len(self.history_t)>1:
            history=CubicHermiteSpline(self.history_t,self.history_e,self.history_d,extrapolate=False)
            clipped=np.minimum(u,self.history_t[-1]);de=history(clipped);df=history(clipped,1)
            future=u>self.history_t[-1];dt=u[future]-self.history_t[-1]
            slope=(self.history_d[-1]-self.history_d[-2])/(self.history_t[-1]-self.history_t[-2])
            de[future]=self.history_e[-1]+self.history_d[-1]*dt+slope*dt**2/2
            df[future]=self.history_d[-1]+slope*dt
        energy=self.f0[-1]*u+de;flux=self.f0[-1]+df
        energy[at<=0]=0;flux[at<=0]=0
        p={k:v[self.outside] for k,v in self.points.items()};r=p['r'];b=1-2*p['m']/r;c=p['N']*np.sqrt(b)
        factor=old.G/(old.C**4*self.bg.R)
        j=-factor*(energy@weights)*np.sqrt(b)/p['N']
        ep=factor*((flux*(1/mu-mu))@weights)/(4*np.pi*r*r*p['N']**2)
        return self.photon_test@(2*c*p['v']*j/b**2-4*np.pi*r**3*c*p['v']*ep/b)

    def stage(self,h):
        h=np.longdouble(h);K=self.K.astype(np.longdouble);M=self.M.astype(np.longdouble)
        F=self.F.astype(np.longdouble);H=self.H.astype(np.longdouble)
        Lq=self.Lq.astype(np.longdouble);LE=self.LE.astype(np.longdouble)
        gain=np.r_[h*h*self.lam/(1+h*self.lam),h]
        GG=K+M/h**2;BB=F-M@H/h**2
        block=bmat([[GG,-BB@diags(self.energy_scale)],[-Lq,(diags(1/gain)-LE)@diags(self.energy_scale)]],format='csc')
        perm=self.permutation;AA=block[perm,:][:,perm].astype(float)
        row=np.asarray(abs(AA).max(axis=1).toarray()).ravel();aa=diags(1/row)@AA
        col=np.asarray(abs(aa).max(axis=0).toarray()).ravel()
        factor=condense_factor((aa@diags(1/col)).tocsc(),self)
        scale=np.sqrt(abs(GG.diagonal()))
        gr=splu((diags(1/scale)@GG@diags(1/scale)).astype(float).tocsc(),permc_spec='NATURAL')
        def invert(rhs):
            ans=np.empty(len(rhs));ans[perm]=factor.solve(np.asarray(rhs[perm]/row,float))/col
            return ans.astype(np.longdouble)
        def step(state,t):
            qr,vr,er,dr=state;pred=er+h*dr
            rhs=M@(qr/h**2+vr/h)+t*self.F0+F@pred+self.photon_force(float(t))
            vec=np.r_[rhs,t*self.LE0+LE@pred-dr];answer=invert(vec)
            for _ in range(3):answer+=invert(vec-block@answer)
            delta=answer[self.size:]*self.energy_scale;target=rhs+BB@delta
            q=(gr.solve(np.asarray(target/scale,float))/scale).astype(np.longdouble)
            for _ in range(3):q+=(gr.solve(np.asarray((target-K@q-M@q/h**2)/scale,float))/scale).astype(np.longdouble)
            answer[:self.size]=q;defect=vec-block@answer
            error=float(max(abs(defect)/(abs(vec)+abs(block)@abs(answer)+1e-100)))
            self.max_error=max(self.max_error,error);assert error<1e-9,error
            e=pred+delta;targetflux=Lq@q+t*self.LE0+LE@e
            d=np.r_[(dr[:-1]+h*self.lam*targetflux[:-1])/(1+h*self.lam),targetflux[-1]]
            # Local check in deviation coordinates: the large background
            # current cannot mask a failed small energy update.
            heat_error=float(max(abs(e-er-h*d)/(abs(e)+abs(er)+abs(h*d)+1e-100)))
            self.max_heat_error=max(self.max_heat_error,heat_error);assert heat_error<1e-9,heat_error
            return q,(q-qr)/h,e,d
        return step

    def evolve(self,seconds,steps,label):
        assert not (OUT/f'{label}.npz').exists()
        start=time.monotonic();gamma=1-1/np.sqrt(np.longdouble(2))
        times=np.longdouble(seconds/self.bg.tc)*(np.arange(steps+1,dtype=np.longdouble)/steps)**2
        state=tuple(np.zeros(n,np.longdouble) for n in [self.size,self.size,self.n,self.n])
        self.history_t=[0.];self.history_e=[0.];self.history_d=[0.]
        self.max_error=self.max_heat_error=0.
        rows=[];temp=[];vel=[];scalar=[];surface=[];weights=self.dm/self.dm.sum()
        p=self.bg.sample(self.native);b=1-2*p['m']/self.native
        speed=old.C/100*self.native/(p['N']*np.sqrt(b))
        for j,t in enumerate(times):
            if j:
                dt=t-times[j-1];stage=self.stage(gamma*dt)
                first=stage(state,times[j-1]+gamma*dt)
                base=tuple(a+(1-gamma)/gamma*(c-a) for a,c in zip(state,first));state=stage(base,t)
                self.history_t.append(float(t));self.history_e.append(float(state[2][-1]));self.history_d.append(float(state[3][-1]))
            q,v,e,d=state;T=self.Tq@q+t*self.TE0+self.TE@e
            velocity=speed*(self.nativeV[0]@v);s=self.nativeV[1]@q
            surf=[float((a@z)[0]) for z in [q,v] for a in self.surfaceV]
            row=dict(t=float(t*self.bg.tc),maximum_delta_lnT=float(max(abs(T))),base_cell_delta_lnT=float(T[self.core_count-1]),
                velocity_RMS_m_s=float(np.sqrt(weights@velocity**2)),maximum_velocity_m_s=float(max(abs(velocity))),
                scalar_RMS=float(np.sqrt(weights@s**2)),surface_displacement=surf[0],surface_scalar=surf[1],
                surface_velocity_coordinate=surf[2],outgoing_luminosity_relative=float(d[-1]/self.f0[-1]))
            rows.append(row);temp.append(T);vel.append(velocity);scalar.append(s);surface.append(surf)
            assert row['maximum_delta_lnT']<.05,('Tangent window exceeded',row)
        dmgeom=self.native**2*b*p['v']*s-4*np.pi*self.native**3*np.exp(-8*p['phi']**2)*(p['e']+p['p'])*(self.nativeV[0]@q)+t*self.J0+self.Jmap@e
        np.savez_compressed(OUT/f'{label}.npz',temperature=temp,velocity=vel,scalar=scalar,surface=surface,
            q=q,v=v,e=e,d=d,f0=self.f0,radius=self.native,edges=self.edges,grid=self.grid,indices=self.indices,
            Eulerian_mass_geom_increment_cm=dmgeom*self.bg.R,emission_times=self.history_t,
            emission_energy_deviation=self.history_e,emission_flux_deviation=self.history_d)
        result=dict(classification='Counterexample candidate',steps=steps,degree=4,graded=True,history=rows,
            seconds=time.monotonic()-start,setup_seconds=self.setup_seconds,max_linear_residual=self.max_error,
            max_local_deviation_heat_identity=self.max_heat_error,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
            moving_surface_solved=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/f'{label}.json',result);print('BALANCED',label,result['seconds'],rows[-1],flush=True)
        return result


def controls():
    import sympy as s
    q,qr,v,vr,e,er,d,dr,h,t,f0,M,H,K,F,Lq,LE,lam=s.symbols('q qr v vr e er d dr h t f0 M H K F Lq LE lam',nonzero=True)
    delta=s.symbols('delta');pred=er+h*dr
    mechanical=M*((q-qr)/h-vr)/h+M*H*delta/h**2+K*q-F*(f0*t+pred+delta)
    block=(K+M/h**2)*q-(F-M*H/h**2)*delta-M*(qr/h**2+vr/h)-F*(f0*t+pred)
    assert s.expand(mechanical-block)==0
    thermal=(dr+delta/h-dr)/h-lam*(Lq*q+LE*(f0*t+pred+delta)-dr-delta/h)
    reduced=-Lq*q+(1/(h*h*lam/(1+h*lam))-LE)*delta-LE*(f0*t+pred)+dr
    assert s.simplify(thermal/lam-reduced)==0
    assert s.expand(er+h*(dr+delta/h)-(pred+delta))==0
    return dict(classification='Proven',passed=True,
        identity='E=f0*t+e, f=f0+d, v=w-H*f. The transformed stage is exactly the original semidiscrete stage in exact arithmetic. Nonzero initial physical currents remain.',
        boundary='This is a numerical-coordinate identity, not physical convergence or a moving-atmosphere theorem.')


def prepare():
    assert not (OUT/'plan.json').exists()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='9418a4d9b',
        claim='Remove the measured constant-current source and velocity cancellation from the actual native coupled evolution before using surface motion for photon pressure work.',
        change='Analytic reference current plus deviations; factor source action after face differences. Identical p4 consistent mass, C1 source, fixed background and surface physics.',
        paths=['p4-8 pilot','p4-32','p4-64'],
        gates=dict(temperature=.02,velocity=.03,scalar=.03,linear_residual=1e-9,local_heat_identity=1e-9,uniform_debit_absolute=0),
        comparison='Time refinement and saved Phase95 same-order trajectories. No new spatial-convergence verdict. Preserve the failed p4/p6 order control.',
        budget=dict(total_seconds=180,CPU_threads=1,memory_GB=5,native_EOS_calls=0,new_background_roots=0),
        stop='Stop on local residual, tangent or measured budget failure. No degree, mesh, horizon, threshold or source-profile expansion.',
        symbolic=controls(),bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),Path(prior.__file__),Path(old.__file__),
            Path(__file__).with_name('def_native_source_order.py'),old.OUT/'inputs.npz',old.OUT/'coefficients.npz',prior.OUT/'consistent-p4-64.npz',OUT/'cancellation-control.json']}))
    (OUT/'reused-assembler.py').write_text(source,encoding='utf-8')


def run():
    assert not (OUT/'result.json').exists();spec=json.loads((OUT/'plan.json').read_text())
    for p,h in spec['bindings'].items():assert old.photons.digest(old.ROOT/p)==h,p
    signal.alarm(175);resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)));start=time.monotonic()
    m=Model();p=m.bg.sample(m.native)
    _,loss,_=m.source_values(p,np.full(m.n,m.f0[-1],np.longdouble))
    assert np.max(abs(loss[2:]))==0
    horizon=json.loads((old.OUT/'plan.json').read_text())['horizon_seconds']
    pilot=m.evolve(horizon,8,'p4-8')
    forecast=time.monotonic()-start+1.5*(pilot['seconds']*96/8+5)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',forecast_seconds=forecast,
        setup_seconds=m.setup_seconds,evolution8_seconds=pilot['seconds'],uniform_debit_absolute=float(np.max(abs(loss[2:]))),
        assumption='Measured same p4 eight-step cost,96 further steps with50percent margin; fixed physical model.'))
    assert forecast<175,'Measured budget exceeded; no automatic expansion'
    for n in [32,64]:m.evolve(horizon,n,f'p4-{n}')
    a=np.load(OUT/'p4-32.npz');b=np.load(OUT/'p4-64.npz');c=np.load(prior.OUT/'consistent-p4-64.npz')
    errors={};oldnew={}
    for f in ['temperature','velocity','scalar']:
        scale=max(float(np.max(abs(b[f]))),1e-100)
        errors[f]=float(np.max(abs(a[f]-b[f][::2]))/scale)
        oldnew[f]=float(np.max(abs(c[f]-b[f]))/scale)
    write(OUT/'result.json',dict(classification='Counterexample candidate',time_passed=all(errors[f]<spec['gates'][f] for f in errors),
        time_relative=errors,old_new_same_p4_relative=oldnew,seconds=time.monotonic()-start,
        prior_spatial_failure_resolved=False,moving_surface_solved=False,final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run']);globals()[parser.parse_args().action]()
