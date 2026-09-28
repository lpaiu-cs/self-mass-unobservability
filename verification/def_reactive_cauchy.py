"""Linear reactive fluid/scalar/metric evolution with causal neutrino stress.

Frozen initial reaction directions, initially empty collisionless neutrinos,
and the saved mechanical free surface. This is not a nonlinear thermal star.
"""
from pathlib import Path
import argparse
import inspect
import json
import time
import numpy as np
import sympy as sp
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu
import def_causal_neutrinos as transport
import def_reactive_structure as reactive

surface=transport.surface
thermal=transport.thermal
h=surface.h
OUT=thermal.OUT.parent/'def-reactive-cauchy'


def emitters(n,flat=False):
    assert not flat
    d=np.load(thermal.OUT/'coefficients.npz');s=np.load(thermal.OUT/'sources.npz')
    geometry=transport.Geometry();support=float(d['faces_cm'][0]/(100*geometry.R))
    edges=np.linspace(0,support,n+1);r=d['radius_cm']/(100*geometry.R)
    power=d['dm']*d['A']**2*d['N']**2*(s['neutrino']+s['thermal_neutrino'])
    ids=np.searchsorted(edges,r,side='right')-1
    assert np.min(ids)>=0 and np.max(ids)<n and support<1
    sums=np.bincount(ids,weights=power,minlength=n)
    centres=np.divide(np.bincount(ids,weights=power*r,minlength=n),sums,
                      out=(edges[:-1]+edges[1:])/2,where=sums>0)
    assert abs(sums.sum()/power.sum()-1)<1e-14
    return edges,centres,sums


# Reuse the frozen, tested null geodesics; only align the mesh with the actual
# native source support. Smearing a material debit into vacuum is not allowed.
source=inspect.getsource(transport.build)
old='_,radius,power=emitters(n,flat);edges=np.linspace(0,2,2*n+1)'
new='inside,radius,power=emitters(n,flat);edges=np.r_[inside,1.,np.linspace(1,2,n+1)[1:]]'
assert source.count(old)==1
namespace=dict(vars(transport),emitters=emitters)
exec(compile(source.replace(old,new),__file__,'exec'),namespace)
build_rays=namespace['build']


def symbolic():
    r,m,p,e,phi,v,beta=sp.symbols('r m p e phi v beta',real=True)
    xi,f,V,J,eta,ga,rr,loss,Er,Pr=sp.symbols('xi f V J eta ga rr loss Er Pr')
    b=1-2*m/r;A4=sp.exp(2*beta*phi**2);alpha=beta*phi;w=e+p
    M=4*sp.pi*r*r*A4*e+r*r*b*v*v/2
    g=m/(r*r*b)+4*sp.pi*r*A4*p/b+r*v*v/2+alpha*v
    F=4*sp.pi*A4/b*(alpha*(e-3*p)+r*v*(e-p))-2*(r-m)*v/(r*r*b)
    H=r*r*b*v*f-(4*sp.pi*r*r*A4*p+r*r*b*v*v/2)*xi
    Dm=H+J;dl=(Dm/r-m*xi/r**2)/b
    xp=-eta/ga-2*xi/r-dl-3*alpha*f
    inc=[xi,Dm,p*(eta-ga*rr),w*eta/ga-loss,f,V]
    DM=sum(sp.diff(M,t)*u for t,u in zip([r,m,p,e,phi,v],inc))+4*sp.pi*r*r*Er
    Hp=sum(sp.diff(H,t)*u for t,u in zip([r,m,p,phi,v,xi,f],[1,M,-w*g,v,F,xp,V+xp*v]))
    target=-(4*sp.pi*r*A4*w/b+r*v*v)*J+4*sp.pi*r*r*(Er-A4*loss)
    assert sp.simplify(DM+xp*M-Hp-target)==0
    # The lifted pressure eta_ad=Delta p/p+Gamma*rho_ref removes cancellation
    # in density, while its exact node jump stays in the pressure equation.
    dg=sp.symbols('dg')
    lhs=-g*(w*eta/ga-loss+p*(eta-ga*rr)+w*xp)-w*dg
    bracket=g*(2*xi/r+dl+3*alpha*f+eta-ga*rr)+g*loss/w-dg
    assert sp.simplify(lhs/p+w*g*(eta-ga*rr)/p-(w*bracket/p-g*(eta-ga*rr)))==0
    return dict(classification='Proven',passed=True,
        mass_identity='(N*a*J)_prime=4*pi*r^2*N*a*(E_nu-A^4*loss_J)',
        conservative_mass='J=-(G/c^4)*1e-7*integrated_redshifted_face_energy/(Rstar*N*a)',
        lifted_pressure='eta_ad=Delta p/p+Gamma*rho_ref; Delta ln rho=eta_ad/Gamma',
        source_order='With initial neutrino intensity zero, metric and displacement corrections to the O(epsilon) rays are O(epsilon^2). This is a first-order frozen-source expansion, not nonlinear metric transport.',
        radiation_terms='Delta g includes 4*pi*r*P_nu/b; Delta F includes 4*pi*r*Phi*(E_nu-P_nu)/b.',
        scope='Linear identities about the declared static mechanical background; photon/conductive heat transport and physical opacity certification are excluded.')


class Radiation:
    def __init__(self,ray):
        self.ray=ray;self.geometry=transport.Geometry();self.edges=ray['edges']
        self.power=np.r_[ray['emitter_power'],np.zeros(len(self.edges)-1-len(ray['emitter_power']))]
        self.volumes=4*np.pi*(100*self.geometry.R)**3*self.integral(self.edges[:-1],self.edges[1:])
        self.events=[];seg=ray['segments'];ids=seg[:,0].astype(int)
        weight=ray['emitter_power'][seg[:,1].astype(int)]*seg[:,2]
        for j in range(len(self.power)):
            selected=ids==j;mu=seg[selected,5]
            val=weight[selected,None]*np.array([np.ones_like(mu),mu,mu*mu]).T
            times=np.r_[seg[selected,3],seg[selected,4]];vals=np.r_[val,-val]
            order=np.argsort(times);times=times[order];vals=vals[order]
            self.events.append((times,np.vstack([np.zeros(3),np.cumsum(vals,axis=0)]),
                               np.vstack([np.zeros(3),np.cumsum(times[:,None]*vals,axis=0)])))

    def integral(self,lo,hi):
        gx,gw=leggauss(16);x=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*gx/2
        N,a=self.geometry.metric(x)
        return (hi-lo)/2*np.sum(gw*N*a*x*x,axis=1)

    def inventory(self,tau):
        out=[]
        for times,slope,offset in self.events:
            i=np.searchsorted(times,tau,side='right')
            out.append((tau*slope[i]-offset[i])*self.geometry.tc)
        return np.array(out)

    def projection(self,r):
        ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,len(self.power)-1)
        fraction=self.integral(self.edges[ids],np.minimum(r,self.edges[ids+1]))
        fraction*=4*np.pi*(100*self.geometry.R)**3/self.volumes[ids]
        fraction=np.clip(fraction,0,1)
        return ids,fraction

    def source(self,tau,bg,projection):
        ids,fraction=projection;inv=self.inventory(tau);debit=self.power*(tau*self.geometry.tc)
        # Each ray crossing telescopes exactly into the same matter debit.
        net=inv[:,0]-debit;prefix=np.r_[0,np.cumsum(net)]
        enclosed=prefix[ids]+fraction*net[ids]
        r=bg['r'];N,a=self.geometry.metric(r);A4=np.exp(-8*bg['phi']**2)
        geo=h.gr.G*.1*self.geometry.R**2/h.gr.C**4
        E=inv[ids,0]/self.volumes[ids]*geo;P=inv[ids,2]/self.volumes[ids]*geo
        loss=debit[ids]/self.volumes[ids]*geo/A4
        J=enclosed*h.gr.G*1e-7/(h.gr.C**4*self.geometry.R*N*a)
        return np.array([bg['rho_ref_rate']*tau*self.geometry.tc,loss,E,P,J])


class Background:
    def __init__(self,radiation,outer):
        self.radiation=radiation;self.saved=np.load(surface.OUT/'background.npz')
        self.Rs=float(self.saved['r'][-1]);self.grid=np.r_[self.saved['grid']/self.Rs,np.linspace(1,outer,128*(outer-1)+1)[1:]]
        self.surface_index=len(self.saved['grid'])-1
        geom=radiation.geometry;phi0=float(self.saved['phi'][-1]);flux=float(self.saved['v'][-1]*self.Rs)
        N,a=geom.metric(np.array([1.]));self.flux=float(N[0]/a[0]*flux)
        def rhs(r,y):
            N,a=geom.metric(np.array([r]));return [self.flux*a[0]/(N[0]*r*r)]
        self.phi=solve_ivp(rhs,[1,outer],[phi0],rtol=2e-12,atol=1e-16,dense_output=True).sol
        self.native=np.load(thermal.OUT/'coefficients.npz');self.reactions=np.load(thermal.OUT/'sources.npz')
        self.forcing=np.load(reactive.chemical.OUT/'forcing.npz')['rows'][::-1,1]
        self.nodes=self.sample(self.grid);self.mid=self.sample((self.grid[:-1]+self.grid[1:])/2)
        self.node_projection=radiation.projection(self.grid);self.mid_projection=radiation.projection(self.mid['r'])

    def sample(self,r):
        s=self.saved;rx=s['r']/self.Rs
        out={k:np.interp(r,rx,s[k]) for k in ['m','p','e','phi','v','N','gamma']};out['r']=r
        out['m']/=self.Rs;out['p']*=self.Rs**2;out['e']*=self.Rs**2;out['v']*=self.Rs
        central=r<rx[1];out['m'][central]=(s['m'][1]/self.Rs)*(r[central]/rx[1])**3
        out['v'][central]=s['v'][1]*self.Rs*r[central]/rx[1]
        outside=r>1
        if np.any(outside):
            q=r[outside];N,a=self.radiation.geometry.metric(q)
            out['m'][outside]=q*(1-1/a**2)/2;out['p'][outside]=0;out['e'][outside]=0
            out['N'][outside]=N;out['phi'][outside]=self.phi(q)[0];out['v'][outside]=self.flux*a/(N*q*q)
        d=self.native;rx=d['radius_cm'][::-1]/(100*self.radiation.geometry.R)
        interp=lambda values:np.interp(r,rx,values[::-1])
        active=r<self.radiation.edges[len(self.radiation.ray['emitter_power'])]
        rho=interp(d['raw'][:,0]);A4=np.exp(-8*out['phi']**2)
        ids,_=self.radiation.projection(r)
        remapped=self.radiation.power[ids]/self.radiation.volumes[ids]/(A4*rho)
        native_loss=interp(d['A']*d['N']*(self.reactions['neutrino']+self.reactions['thermal_neutrino']))
        delta=native_loss-remapped
        rr=np.interp(r,rx,self.forcing[:,0])+interp(d['raw'][:,8]/d['thermo'][:,5])*delta
        theta=np.interp(r,rx,self.forcing[:,2])-interp(d['thermo'][:,4])*np.interp(r,rx,self.forcing[:,0])+delta/interp(d['thermo'][:,3])
        out['rho_ref_rate']=np.where(active,rr,0);out['theta_rate']=np.where(active,theta,0)
        out['adiabatic_T_rho']=interp(d['thermo'][:,4])
        return out


def operators(bg,fn):
    r,m,p,e,phi,v,N,ga=[bg[k] for k in ['r','m','p','e','phi','v','N','gamma']]
    b=1-2*m/r;A4=np.exp(-8*phi*phi);alpha=-4*phi;w=e+p
    vals=fn(r,m,p,e,phi,v,-4.);M,g,F=vals[:3];dg,dF=vals[9:15],vals[15:21]
    basis=np.eye(9)[:,:,None]+np.zeros((9,9,len(r)))
    z,eta,f,V,rr,loss,E,P,J=basis;xi=r*z
    Dm=r*r*b*v*f-(4*np.pi*r*r*A4*p+r*r*b*v*v/2)*xi+J
    dl=(Dm/r-m*xi/r**2)/b;xp=-eta/ga-2*z-dl-3*alpha*f
    outside=r>1;xp[:,outside]=z[:,outside]
    inc=[xi,Dm,p*(eta-ga*rr),w*eta/ga-loss,f,V]
    Dg=sum(q*u for q,u in zip(dg,inc))+4*np.pi*r*P/b
    DF=sum(q*u for q,u in zip(dF,inc))+4*np.pi*r*v*(E-P)/b
    bracket=g*(2*z+dl+3*alpha*f+eta-ga*rr)-Dg
    bracket+=g*loss/np.where(w>0,w,1)
    with np.errstate(divide='ignore',invalid='ignore'):
        result=np.array([(xp-z)/r,w/p*bracket-g*(eta-ga*rr),V+xp*v,DF+xp*F])
    result[:2,:,outside]=0
    inertia=np.zeros((len(r),4,4));inside=~outside
    inertia[inside,1,0]=-w[inside]/p[inside]*r[inside]/(N[inside]**2*b[inside])
    inertia[:,3,0]=-r*v/(N*N*b);inertia[:,3,2]=1/(N*N*b)
    return np.moveaxis(result,-1,0),inertia,bracket,g


def assemble(bg,fn):
    grid=bg.grid;n=len(grid)-1;dx=np.diff(grid);ops,B,_,_=operators(bg.mid,fn)
    A=ops[:,:,:4];S=ops[:,:,4:];rows=[];cols=[];kv=[];dv=[]
    for shift,sign in [(0,-1),(4,1)]:
        rows.extend((np.arange(4*n).reshape(n,4,1)+np.zeros((1,1,4),int)).ravel())
        cols.extend((4*np.arange(n)[:,None,None]+np.arange(4)[None,None,:]+shift+np.zeros((1,4,1),int)).ravel())
        kv.extend((sign*np.eye(4)[None]-dx[:,None,None]*A/2).ravel())
        dv.extend((dx[:,None,None]*B/2).ravel())
    si=bg.surface_index;end={k:np.array([v[si]]) for k,v in bg.nodes.items()}
    with np.errstate(divide='ignore',invalid='ignore'):_,_,bracket,g=operators(end,fn)
    ga=bg.nodes['gamma'][0];alpha=-4*bg.nodes['phi'][0]
    N=end['N'][0];b=1-2*end['m'][0];outer=grid[-1]
    boundary=[(0,[3*ga,1,3*ga*alpha,0],[0]*4),(0,[0,0,0,1],[0]*4),
              (4*si,bracket[:4,0]/g[0],[-1/(N*N*b*g[0]),0,0,0]),
              (4*n,[-outer*bg.nodes['v'][-1],0,1,0],[0]*4)]
    for j,(offset,kr,dr) in enumerate(boundary):
        rows.extend([4*n+j]*4);cols.extend(offset+np.arange(4));kv.extend(kr);dv.extend(dr)
    shape=(4*(n+1),)*2
    K=coo_matrix((kv,(rows,cols)),shape=shape).tocsc();D=coo_matrix((dv,(rows,cols)),shape=shape).tocsc()
    def forcing(tau):
        src=bg.radiation.source(tau,bg.mid,bg.mid_projection)
        out=np.zeros(shape[0]);out[:4*n]=(dx[:,None]*np.einsum('nij,jn->ni',S,src)).ravel()
        node=bg.radiation.source(tau,bg.nodes,bg.node_projection)
        out[1:4*n:4]+=np.diff(bg.nodes['gamma']*node[0])
        out[4*n+2]=-bracket[4:,0]@node[:,si]/g[0]
        return out
    return K,D,forcing


def solve(ray,steps,outer=2,label=None):
    began=time.monotonic();radiation=Radiation(ray);bg=Background(radiation,outer);fn,_=reactive.symbolic()
    K,D,forcing=assemble(bg,fn);dt=1/steps;matrix=K-4*D/dt**2
    scale=np.asarray(abs(matrix).sum(1)).ravel();lu=splu(matrix.multiply((1/scale)[:,None]).tocsc());extended=matrix.astype(np.longdouble)
    y=np.zeros(matrix.shape[0]);velocity=y.copy();acceleration=y.copy();max_residual=0.;history=[]
    native_r=bg.native['radius_cm'][::-1]/(100*radiation.geometry.R);weights=bg.native['dm'][::-1];weights/=weights.sum()
    def readout(tau):
        yy=y.reshape(-1,4);vv=velocity.reshape(-1,4);r=bg.grid
        N,a=radiation.geometry.metric(r);physical_v=a/N*r*vv[:,0]*h.gr.C
        scalar=yy[:,2]-r*yy[:,0]*bg.nodes['v']
        cv=np.interp(native_r,r,physical_v);cf=np.interp(native_r,r,scalar)
        return dict(tau=tau,seconds=tau*radiation.geometry.tc,
                    velocity_mass_RMS_m_s=float(np.sqrt(weights@(cv*cv))),
                    scalar_mass_RMS=float(np.sqrt(weights@(cf*cf))),
                    surface_displacement_m=float(yy[bg.surface_index,0]*radiation.geometry.R),
                    surface_Euler_scalar=float(scalar[bg.surface_index]))
    history.append(readout(0.))
    for j in range(1,steps+1):
        rhs=forcing(j*dt)-D@(4*y/dt**2+4*velocity/dt+acceleration)
        answer=lu.solve(rhs/scale)
        for _ in range(2):
            defect=rhs.astype(np.longdouble)-extended@answer.astype(np.longdouble)
            answer+=lu.solve(np.asarray(defect/scale,float))
        residual=float(np.max(abs(matrix@answer-rhs)/(abs(matrix)@abs(answer)+abs(rhs)+1e-100)))
        max_residual=max(max_residual,residual)
        new_acceleration=4*(answer-y-dt*velocity)/dt**2-acceleration
        velocity+=dt*(acceleration+new_acceleration)/2;acceleration=new_acceleration;y=answer
        assert np.all(np.isfinite(y)) and max_residual<1e-9,(j,max_residual)
        history.append(readout(j*dt))
    row=dict(classification='Counterexample candidate',steps=steps,outer=outer,cells=len(bg.grid)-1,
             seconds=time.monotonic()-began,linear_residual=max_residual,history=history)
    if label:
        np.savez_compressed(OUT/(label+'.npz'),grid=bg.grid,response=y.reshape(-1,4),velocity=velocity.reshape(-1,4),acceleration=acceleration.reshape(-1,4))
        h.write(OUT/(label+'.json'),row)
    print(json.dumps({k:v for k,v in row.items() if k!='history'}),flush=True)
    return row


def prepare():
    assert not OUT.exists();OUT.mkdir();began=time.monotonic()
    prior=thermal.OUT.parent/'gr-causal-neutrino-milestone-manifest.json'
    for rel,sha in json.loads(prior.read_text())['sha256'].items():assert h.digest(h.ROOT/rel)==sha,rel
    paths=[Path(__file__),Path(transport.__file__),Path(reactive.__file__),surface.OUT/'background.npz',
           thermal.OUT/'coefficients.npz',thermal.OUT/'sources.npz',reactive.chemical.OUT/'forcing.npz',prior]
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='5a72dea4',symbolic=symbolic(),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        claim='Evolve actual native initial reaction directions, material inertia, scalar waves and metric constraints with causal neutrino energy and radial pressure on the same free-surface background.',
        approximation='First order in the frozen native thermochemical source; no photons/conductive heat flux, no nonlinear reaction update, no external orbital drive. Exterior displacement is a coordinate label, not vacuum matter.',
        sources='Support-aligned ray cells; identical redshifted energy removed from matter and credited to neutrinos. Replace the native neutrino loss in the fixed-pressure thermochemical tangent by that exact shared-bin debit; preserve its chemical contribution.',
        grid='Original 5863 interior intervals; 128 exterior intervals per Rstar. Radiation 48x24 and 96x48, no adaptive refinement.',
        boundary='Regular centre, dynamic pressure regularity at the true free surface, zero Eulerian scalar at 2Rstar. Control at 3Rstar with identical dx. Horizon Rstar/c.',
        time_steps=[16,32,64],horizon_crossings=1.,
        gates=dict(linear_residual=1e-9,event_replay=2e-12,energy_balance=2e-13,
                   time_readout_relative=.02,time_order_min=1.5,radiation_readout_relative=.03,outer_readout_relative=.002),
        readouts=['velocity_mass_RMS_m_s','scalar_mass_RMS'],
        budget=dict(pilot_hard_seconds=60,production_hard_seconds=180,workers=1,BLAS_threads=1,native_calls=0,
                    GPU=False,estimated_peak_memory_MiB=600,maximum_production_paths=5,automatic_expansion=False)))
    ray,measurement=build_rays(12,12);row=solve(ray,8,label='pilot')
    forecast=measurement['seconds']*((96/12)**2*4+(48/12)**2*2)*1.5+row['seconds']*(16+32+64+64+64)/8*1.5
    h.write(OUT/'pilot-budget.json',dict(ray=measurement,solve=row['seconds'],production_forecast_seconds=forecast,
            assumptions='Linear solve time proportional to steps; fine ray build cubic work scaling. Fine source/angle and extended domain speed not yet measured.',seconds=time.monotonic()-began))
    print('FORECAST_SECONDS',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();began=time.monotonic();plan=json.loads((OUT/'plan.json').read_text())
    for rel,sha in plan['bindings'].items():assert h.digest(h.ROOT/rel)==sha,rel
    assert json.loads((OUT/'pilot-budget.json').read_text())['production_forecast_seconds']<165
    cases={};controls={}
    for name,n,angles in [('coarse',48,24),('fine',96,48)]:
        ray,measurement=build_rays(n,angles);np.savez_compressed(OUT/(name+'-rays.npz'),**ray)
        rad=Radiation(ray);times=np.array([.125,.5,1.]);rows,_,stress=transport.moments(ray,times)
        replay=max(float(np.max(abs(rad.inventory(t)-ref.T)))/(ray['emitter_power'].sum()*ray['tc']) for t,ref in zip(times,stress))
        balance=max(q['energy_balance_relative'] for q in rows)
        assert replay<2e-12 and balance<2e-13,(replay,balance)
        controls[name]=dict(measurement=measurement,event_replay=replay,energy_balance=balance,history=rows)
        if name=='fine':
            for steps in [16,32,64]:cases[str(steps)]=solve(ray,steps,label=f'fine-{steps}')
            cases['outer']=solve(ray,64,3,label='outer-64')
        else:cases['coarse']=solve(ray,64,label='coarse-64')
        assert time.monotonic()-began<175
    fields=plan['readouts'];comparisons={}
    for field in fields:
        a=np.array([r[field] for r in cases['16']['history']]);b=np.array([r[field] for r in cases['32']['history']]);c=np.array([r[field] for r in cases['64']['history']])
        norm=max(abs(c).max(),1e-100);d1=np.max(abs(a-b[::2]))/norm;d2=np.max(abs(b-c[::2]))/norm
        space=np.max(abs(c-np.array([r[field] for r in cases['coarse']['history']])))/norm
        outer=np.max(abs(c-np.array([r[field] for r in cases['outer']['history']])))/norm
        comparisons[field]=dict(time_last_relative=float(d2),time_previous_relative=float(d1),
             order=float(np.log2(max(d1,1e-100)/max(d2,1e-100))),radiation_relative=float(space),outer_relative=float(outer))
    passed=all(q['time_last_relative']<.02 and q['order']>1.5 and q['radiation_relative']<.03 and q['outer_relative']<.002 for q in comparisons.values())
    result=dict(classification='Counterexample candidate',coupled_linear_time_paths_completed=True,all_gates_passed=passed,
                comparisons=comparisons,controls=controls,seconds=time.monotonic()-began,
                nonlinear_reactions_evolved=False,photon_heat_transport=False,physical_neutrino_opacity_certified=False,
                mechanical_spatial_error_certified=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
