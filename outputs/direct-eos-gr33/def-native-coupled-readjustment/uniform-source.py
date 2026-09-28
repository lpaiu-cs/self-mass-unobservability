"""Native core/envelope thermal readjustment in the regular GR weak form.

Counterexample candidate: tangent LTE grey transport with fixed coefficients,
finite-pressure outgoing boundary and causal photons. This is not nonlinear
radiation hydrodynamics or an asymptotic charge measurement.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse
import json
import resource
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix, diags, bmat, eye
from scipy.sparse.linalg import splu
from scipy.interpolate import CubicHermiteSpline
import def_native_whole_star_match as prior
import def_native_photon_exterior as photons
import def_gr_hierarchical as hierarchy
import def_native_temperature_closure as closure

ROOT=prior.ROOT
OUT=prior.OUT.parent/'def-native-coupled-readjustment'
write=prior.write
G,C=photons.G,photons.C


def symbolic():
    import sympy as s
    r,m,p,e,A,c,v,alpha=s.symbols('r m p e A c v alpha', nonzero=True)
    b=1-2*m/r
    T10=4*s.pi*A*c*r**3*(e+p)*(alpha+r*v)/b
    jump=-4*s.pi*A/b*(alpha*(e-3*p)+r*v*(e-p))
    result=8*s.pi*A*c*r**3*p*(2*alpha+r*v)/b
    assert s.simplify(r**3*c*jump+T10-result)==0
    return dict(classification='Proven',passed=True,temperature=closure.symbolic(),
        finite_surface='Continuity of the Lagrangian scalar gradient gives Pcanonical_inside-Pcanonical_outside=8*pi*A4*C*r^3*p*(2*alpha+r*Phi)*zeta. It vanishes only at zero pressure.',
        pressure='Pcanonical_fluid=-4*pi*A4*C*r^3*DeltaP/b. The frozen emitting geometry boundary has DeltaP=4*Prad*Delta lnT, with fixed external gas pressure.',
        scope='Algebraic boundary and thermodynamic identities in the declared tangent model. Photon corrections to background coefficients, moving-ray geometry, nuclear evolution and non-LTE opacity are outside this model.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    inputs=[Path(__file__),Path(prior.__file__),Path(photons.__file__),Path(hierarchy.__file__),
        Path(hierarchy.task.__file__),Path(closure.__file__),prior.OUT/'background.npz',
        prior.OUT/'final-core.npz',prior.OUT/'final-envelope.npz',prior.OUT/'audit.json',
        OUT.parent/'def-gr-canonical/regular/symbolic.json']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='4c6f99a28',
        claim='Apply the actual Phase92 native core/envelope background to simultaneous thermal, material, scalar and metric evolution. Retain the base cooling, nonzero outgoing flux and finite-pressure scalar junction.',
        equations='qdot=w-H*f; M*wdot=-K*q+F*E+causal_photon_force; Edot=f; internal fdot=lambda*(tc*L0+Lq*q+LE*E-f); outgoing f=tc*L0+4*tc*L0*surface_Lagrangian_delta_lnT.',
        initial='q=E=0, physical qdot=0, nonzero stored constitutive currents f=tc*L0 and w=H*f. No artificial zero-current startup. Internal core gradients set L0; shared base and radiative envelope use the matched L0.',
        scope='Frozen LTE grey Cattaneo coefficients and composition, temperature feedback, fixed emitting geometry. Finite gas pressure is held by an explicit bath. Rays use the same emitted-energy history. No claimed physical non-LTE atmosphere, complete metric-dependent transport, continuum proof or final charge.',
        horizon_seconds=json.loads((prior.OUT/'audit.json').read_text())['frozen_density_one_percent_timescale_s'],
        paths=['p4-16','p4-32','p4-64','p2-64'],
        gates=dict(temperature_time=.02,temperature_space=.02,GR_time=.03,GR_space=.03,
            linear_residual=1e-9,heat_identity=1e-9,inventory=1e-12,tangent_temperature=.05),
        decision='Accept only the fields whose independent time/space controls pass. If any gate fails preserve the failure and diagnose before changing resolution, horizon or physical equations. No automatic extension.',
        budget=dict(total_seconds=600,pilot_seconds=60,CPU_threads=1,memory_GB=5,native_calls=6003,
            maximum_production_paths=4,automatic_expansion=False),
        bindings={str(p.relative_to(ROOT)):photons.digest(p) for p in inputs},symbolic=symbolic()))
    start=time.monotonic();signal.alarm(60)
    bg=np.load(prior.OUT/'background.npz');core=np.load(prior.OUT/'final-core.npz');env=np.load(prior.OUT/'final-envelope.npz')
    n=len(core['dm']);rc=bg['radius_cm'][1:n+1]
    re=(env['r'][:-1]+env['r'][1:])/2
    radius=np.r_[rc,re]
    pressure=np.r_[bg['pressure_cgs'][1:n+1],np.interp(re,env['r'],env['Ptotal'])]
    temperature=np.r_[bg['temperature_K'][1:n+1],np.exp(np.interp(re,env['r'],np.log(env['T'])))]
    X=np.r_[core['X'][::-1],np.broadcast_to(core['X'][0],(len(re),26))]
    envmodel=prior.envelope.Envelope();total=envmodel.total_baryon
    dm=np.r_[core['dm'][::-1],np.diff(env['baryon_outside_base'])*total]
    edges=np.r_[bg['core_faces'][::-1,0]*float(core['R_scale_m'])*100,env['r'][1:]]
    assert len(edges)==len(radius)+1 and np.all(np.diff(edges)>0) and np.all(dm>0)
    ids=np.unique(np.r_[np.linspace(0,len(radius)-1,31).astype(int),n-1])
    before=time.monotonic();rows=[]
    for i in ids:
        raw=envmodel.eos(1,float(np.log(pressure[i])),float(np.log(temperature[i])),X[i])
        opacity=prior.envelope.prior.two.opacity_parts(envmodel.opacity,(float(np.log(raw[0])),float(np.log(temperature[i])),X[i]))[0]
        rows.append(np.r_[raw,opacity])
    elapsed=time.monotonic()-before
    np.savez_compressed(OUT/'inputs.npz',radius=radius,pressure=pressure,temperature=temperature,X=X,dm=dm,edges=edges,
        pilot_ids=ids,pilot_rows=rows,core_count=n)
    forecast=elapsed*len(radius)/len(ids)*1.5
    row=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,
        pilot_calls=len(ids),pilot_call_seconds=elapsed,coefficient_forecast_seconds=forecast,
        cells=len(radius),whole_baryon_g=float(dm.sum()),envelope_baryon_g=float(dm[n:].sum()),
        measured_core_and_envelope=True,model_setup_and_stage_speed_unmeasured=True)
    write(OUT/'native-pilot.json',row);signal.alarm(0);print('NATIVE PILOT',json.dumps(row),flush=True)


def coefficients():
    assert not (OUT/'coefficients.npz').exists()
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():
        actual=OUT/'registered-source.py' if p==str(Path(__file__).relative_to(ROOT)) else ROOT/p
        assert photons.digest(actual)==sha,p
    review=json.loads((OUT/'budget-review.json').read_text());cap=review['native_coefficient_cap_seconds']
    pilot=json.loads((OUT/'native-pilot.json').read_text());assert pilot['coefficient_forecast_seconds']<cap
    signal.alarm(cap);start=time.monotonic();d=np.load(OUT/'inputs.npz');model=prior.envelope.Envelope()
    rows=np.empty((len(d['radius']),22));rows[d['pilot_ids']]=d['pilot_rows'];known=set(d['pilot_ids'])
    for i,(P,T,X) in enumerate(zip(d['pressure'],d['temperature'],d['X'])):
        if i in known:continue
        raw=model.eos(1,float(np.log(P)),float(np.log(T)),X)
        op=prior.envelope.prior.two.opacity_parts(model.opacity,(float(np.log(raw[0])),float(np.log(T)),X))[0]
        rows[i]=np.r_[raw,op]
    thermo=prior.envelope.prior.transformed(rows[:,:21]);assert np.all(thermo[:,[3,5]]>0)
    error=float(max(abs(thermo[:,0]+thermo[:,1]*thermo[:,4]-rows[:,4])/rows[:,4]))
    assert error<1e-7,error
    np.savez_compressed(OUT/'coefficients.npz',raw=rows[:,:21],opacity=rows[:,21],thermo=thermo)
    row=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,native_calls=len(rows)-len(known),
        gamma_chain_relative=error,minimum_cvT=float(thermo[:,3].min()))
    write(OUT/'coefficients.json',row);signal.alarm(0);print('COEFFICIENTS',json.dumps(row),flush=True)


class Background(photons.Exterior):
    def __init__(self):
        self.d=dict(np.load(prior.OUT/'background.npz'));self.R=float(self.d['radius_cm'][-1]);self.tc=self.R/C
        self.vb=float(self.d['phi_prime_cm'][-1])*self.R;self.mu=float(self.d['mass_geom_cm'][-1])/self.R
        self.flat=False;self.M,self.K,self.solution=photons.vacuum.background(self.mu,self.vb,2e-13)
        self.Nb=float(self.metric(np.array([1.]))[1][0]);self.x=self.d['radius_cm']/self.R
        geo=G*self.R**2/C**4
        self.saved=dict(r=self.x,e=self.d['energy_cgs']*geo,p=self.d['pressure_cgs']*geo,gamma=self.d['gamma1'])
        self.Rs=1.

    def sample(self,r):
        r=np.asarray(r);d=self.d
        keys=['mass_geom_cm','pressure_cgs','energy_cgs','phi','phi_prime_cm','lapse','gamma1']
        scales=[1/self.R,G*self.R**2/C**4,G*self.R**2/C**4,1,self.R,1,1]
        z={k:np.interp(r,self.x,d[name])*scale for k,name,scale in zip(['m','p','e','phi','v','N','gamma'],keys,scales)}
        inside=r<self.x[1]
        z['m'][inside]=d['mass_geom_cm'][1]/self.R*(r[inside]/self.x[1])**3
        z['v'][inside]=d['phi_prime_cm'][1]*self.R*r[inside]/self.x[1]
        outside=r>1
        if np.any(outside):
            mass,n,b,c,v=self.metric(r[outside]);z['m'][outside]=mass;z['N'][outside]=n;z['v'][outside]=v
            z['p'][outside]=z['e'][outside]=0;z['gamma'][outside]=1
        z['r']=r
        return z


class Model:
    evaluation=hierarchy.Model.evaluation

    def __init__(self,degree):
        started=time.monotonic();self.bg=bg=Background();d=np.load(OUT/'inputs.npz');th=np.load(OUT/'coefficients.npz')
        self.native=d['radius']/bg.R;self.edges=d['edges']/bg.R;self.dm=d['dm'];self.n=len(self.native)
        self.raw=th['raw'];self.thermo=th['thermo'];self.temperature=d['temperature'];self.core_count=int(d['core_count'])
        self.degree=degree;self.outer=1.1
        # ponytail: retain every thermal cell; mechanical p2/p4 controls share
        # the registered mesh, not an automatically expanded acoustic mesh.
        k=self.core_count
        ci=np.unique(np.r_[np.arange(0,k,8),np.arange(max(0,k-40),k+1)])
        self.cells=np.unique(np.r_[self.edges[ci],self.edges[k::4],1.,np.linspace(1,1.1,257)])
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

    def face_map(self,r):
        ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,self.n-1)
        left=self.edges[ids];right=np.minimum(r,self.edges[ids+1]);width=np.maximum(right-left,0)
        gx,gw=np.polynomial.legendre.leggauss(6)
        x=left[:,None]+width[:,None]*(gx+1)/2;p=self.bg.sample(x.ravel())
        b=1-2*p['m']/np.maximum(p['r'],1e-100)
        volume=4*np.pi*self.bg.R**3*width/2*np.sum((p['N']/np.sqrt(b)*p['r']**2).reshape(x.shape)*gw,axis=1)
        frac=np.clip(volume/self.volumes[ids],0,1)
        rows=np.repeat(np.arange(len(r)),2);cols=np.c_[ids-1,ids].ravel();data=np.c_[1-frac,frac].ravel();valid=cols>=0
        return coo_matrix((data[valid],(rows[valid],cols[valid])),shape=(len(r),self.n)).tocsr()

    def sources(self,p):
        r=p['r'];ids=np.clip(np.searchsorted(self.edges,r,side='right')-1,0,self.n-1);inside=r<=1
        rows=np.repeat(np.arange(len(r)),2);cols=np.c_[ids-1,ids].ravel();valid=(cols>=0)&np.repeat(inside,2)
        inc=coo_matrix((np.tile([1.,-1.],len(r))[valid],(rows[valid],cols[valid])),shape=(len(r),self.n)).tocsr()
        geo=G*self.bg.R**2/C**4;A4=np.exp(-8*p['phi']**2)
        loss=diags(-geo/(self.volumes[ids]*A4))@inc
        ratio=self.raw[ids,8]/self.thermo[ids,5]/self.raw[ids,0]
        rr=diags(-ratio/geo)@loss
        b=1-2*p['m']/np.maximum(r,1e-100)
        J=diags(-G/(C**4*self.bg.R)*np.sqrt(b)/p['N']*inside)@self.face_map(r)
        return rr.tocsc(),loss.tocsc(),J.tocsc()

    def temperature_maps(self):
        r=self.native;p=self.bg.sample(r);b=1-2*p['m']/r;ad=self.thermo[:,4];alpha=-4*p['phi']
        V,D=self.nativeV,self.nativeD
        slope=np.gradient(np.log(self.temperature),r)
        dlz=r*r*p['v']**2/2-4*np.pi*r*r*np.exp(-8*p['phi']**2)*p['p']/b-p['m']/(r*b)
        Q=-diags(ad*r)@D[0]-diags(ad*(3+dlz+3*alpha*r*p['v']))@V[0]-diags(ad*(r*p['v']+3*alpha))@V[1]
        self.TLq=Q.tocsc();self.Tq=(Q-diags(r*slope)@V[0]).tocsc()
        _,loss,J=self.sources(p);geo=G*self.bg.R**2/C**4
        self.TE=(-diags(ad/(r*b))@J-diags(1/(self.raw[:,0]*geo*self.thermo[:,3]))@loss).tocsc()

    def photon_force(self,t):
        mu,weights,delay=self.rays;at=t-delay;clipped=np.maximum(at,0)
        if len(self.history_t)>1:
            history=CubicHermiteSpline(self.history_t,self.history_E,self.history_f,extrapolate=False)
            where=np.minimum(clipped,self.history_t[-1]);energy=history(where);flux=history(where,1)
            future=clipped>self.history_t[-1];delta=clipped[future]-self.history_t[-1]
            slope=(self.history_f[-1]-self.history_f[-2])/(self.history_t[-1]-self.history_t[-2])
            energy[future]=self.history_E[-1]+self.history_f[-1]*delta+slope*delta**2/2
            flux[future]=self.history_f[-1]+slope*delta
        else:
            energy=clipped*self.history_f[0];flux=np.full_like(clipped,self.history_f[0])
        flux[at<=0]=0;energy[at<=0]=0
        p={k:v[self.outside] for k,v in self.points.items()};r=p['r'];b=1-2*p['m']/r;c=p['N']*np.sqrt(b)
        factor=G/(C**4*self.bg.R)
        j=-factor*(energy@weights)*np.sqrt(b)/p['N']
        ep=factor*((flux*(1/mu-mu))@weights)/(4*np.pi*r*r*p['N']**2)
        force=2*c*p['v']*j/b**2-4*np.pi*r**3*c*p['v']*ep/b
        return self.photon_test@force

    def stage(self,h):
        # Same SDIRK elimination as Phase83/89, including an algebraic outgoing
        # face. Scaling and residuals use the actual new thermal/GR blocks.
        h=np.longdouble(h);K=self.K.astype(np.longdouble);M=self.M.astype(np.longdouble);H=self.H.astype(np.longdouble)
        F=self.F.astype(np.longdouble);Lq=self.Lq.astype(np.longdouble);LE=self.LE.astype(np.longdouble)
        gain=np.r_[h*h*self.lam/(1+h*self.lam),h]
        GG=K+M/(h*h);BB=F-M@H/(h*h)
        block=bmat([[GG,-BB@diags(self.energy_scale)],[-Lq,(diags(1/gain)-LE)@diags(self.energy_scale)]],format='csc')
        perm=self.permutation;AA=block[perm,:][:,perm].astype(float)
        row=np.asarray(abs(AA).max(axis=1).toarray()).ravel();aa=diags(1/row)@AA
        col=np.asarray(abs(aa).max(axis=0).toarray()).ravel();factor=splu((aa@diags(1/col)).tocsc(),permc_spec='NATURAL')
        scale=np.sqrt(abs(GG.diagonal()));grfactor=splu((diags(1/scale)@GG@diags(1/scale)).astype(float).tocsc(),permc_spec='NATURAL')
        def invert(rhs):
            ans=np.empty(len(rhs));ans[perm]=factor.solve(np.asarray(rhs[perm]/row,float))/col
            return ans.astype(np.longdouble)
        def step(state,t):
            qr,wr,Er,fr=state
            rhs=M@(qr/(h*h)+wr/h+H@Er/(h*h))+self.photon_force(float(t))
            Eb=Er.copy();Eb[:-1]+=h*(fr[:-1]+h*self.lam*self.bg.tc*self.L0[:-1])/(1+h*self.lam)
            Eb[-1]+=h*self.bg.tc*self.L0[-1]
            vec=np.r_[rhs,Eb/gain];answer=invert(vec)
            for _ in range(3):answer+=invert(vec-block@answer)
            E=answer[self.size:]*self.energy_scale;target=rhs+BB@E
            q=(grfactor.solve(np.asarray(target/scale,float))/scale).astype(np.longdouble)
            for _ in range(3):q+=(grfactor.solve(np.asarray((target-K@q-M@q/(h*h))/scale,float))/scale).astype(np.longdouble)
            answer[:self.size]=q;defect=vec-block@answer
            error=float(max(abs(defect)/(abs(vec)+abs(block)@abs(answer)+1e-100)));self.max_error=max(self.max_error,error)
            assert error<1e-9,error
            targetflux=self.bg.tc*self.L0+Lq@q+LE@E
            f=np.r_[(fr[:-1]+h*self.lam*targetflux[:-1])/(1+h*self.lam),targetflux[-1]]
            heat_error=float(max(abs(E-Er-h*f))/max(max(abs(E)),max(abs(Er)),1e-100));self.max_heat_error=max(self.max_heat_error,heat_error)
            assert heat_error<1e-9,heat_error
            w=(q-qr+H@(E-Er))/h
            return q,w,E,f
        return step

    def evolve(self,seconds,steps,label):
        start=time.monotonic();dt=np.longdouble(seconds/self.bg.tc/steps);gamma=1-1/np.sqrt(np.longdouble(2));stage=self.stage(gamma*dt)
        f=(self.bg.tc*self.L0).astype(np.longdouble);state=(np.zeros(self.size,np.longdouble),self.H@f,np.zeros(self.n,np.longdouble),f)
        self.history_t=[0.];self.history_E=[0.];self.history_f=[float(f[-1])]
        temp=[];vel=[];scalar=[];rows=[];weights=self.dm/self.dm.sum();p=self.bg.sample(self.native)
        speed=C/100*self.native/(p['N']*np.sqrt(1-2*p['m']/self.native))
        for j in range(steps+1):
            if j:
                first=stage(state,(j-1+gamma)*dt);base=tuple(a+(1-gamma)/gamma*(b-a) for a,b in zip(state,first));state=stage(base,j*dt)
                self.history_t.append(float(j*dt));self.history_E.append(float(state[2][-1]));self.history_f.append(float(state[3][-1]))
            q,w,E,f=state
            T=self.Tq@q+self.TE@E;v=speed*(self.nativeV[0]@(w-self.H@f));s=self.nativeV[1]@q
            dmgeom=self.native**2*(1-2*p['m']/self.native)*p['v']*s-4*np.pi*self.native**3*np.exp(-8*p['phi']**2)*(p['e']+p['p'])*(self.nativeV[0]@q)+self.sources(p)[2]@E
            balance=float(abs(np.sum(-np.diff(np.r_[np.longdouble(0),E]),dtype=np.longdouble)+E[-1])/max(abs(E).max(),1e-100))
            row=dict(t=float(j*dt*self.bg.tc),maximum_delta_lnT=float(max(abs(T))),base_cell_delta_lnT=float(T[self.core_count-1]),
                velocity_RMS_m_s=float(np.sqrt(weights@v**2)),scalar_RMS=float(np.sqrt(weights@s**2)),
                surface_scalar=float((self.evaluation(np.array([1.]))[0][1]@q)[0]),
                outgoing_luminosity_relative=float(f[-1]/(self.bg.tc*self.L0[-1])-1),energy_balance=balance)
            rows.append(row);temp.append(T);vel.append(v);scalar.append(s)
            if row['maximum_delta_lnT']>=.05:
                write(OUT/f'{label}-stopped.json',dict(classification='Counterexample candidate',reason='Registered tangent window exceeded',history=rows))
                raise AssertionError(('Tangent window exceeded',row))
        np.savez_compressed(OUT/f'{label}.npz',temperature=temp,velocity=vel,scalar=scalar,q=state[0],w=state[1],E=state[2],f=state[3],
            radius=self.native,edges=self.edges,grid=self.grid,indices=self.indices,Eulerian_mass_geom_increment_cm=dmgeom*self.bg.R,
            emission_times=self.history_t,emission_energy=self.history_E,emission_flux=self.history_f)
        result=dict(classification='Counterexample candidate',steps=steps,degree=self.degree,history=rows,
            seconds=time.monotonic()-start,setup_seconds=self.setup_seconds,max_linear_residual=self.max_error,max_heat_identity=self.max_heat_error,
            memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,full_dynamic_charge_solved=False)
        write(OUT/f'{label}.json',result);print('EVOLVED',label,result['seconds'],rows[-1],flush=True);return result


def pilot():
    assert not (OUT/'pilot.json').exists();signal.alarm(90);start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    p=Model(4);h=json.loads((OUT/'plan.json').read_text())['horizon_seconds'];row=p.evolve(h,8,'pilot')
    prior_cost=sum(json.loads((OUT/f).read_text())['seconds'] for f in ['native-pilot.json','coefficients.json'])
    forecast=prior_cost+1.5*(3*p.setup_seconds+row['seconds']*176/8+20)
    result=dict(classification='Counterexample candidate',seconds=time.monotonic()-start,prior_seconds=prior_cost,forecast_total_seconds=forecast,
        assumption='Measured p4 setup/eight-step cost scaled to176 production steps and three setups plus50percent. The p2 cost is not measured.',
        dofs=p.size,thermal_cells=p.n,quadrature_points=len(p.weights),row=row)
    write(OUT/'pilot.json',result);signal.alarm(0);print('BUDGET',forecast,flush=True)


def run():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'plan.json').read_text());pilot=json.loads((OUT/'pilot.json').read_text())
    assert pilot['forecast_total_seconds']<600,'Reassess measured budget before proceeding'
    previous=pilot['prior_seconds']+pilot['seconds'];signal.alarm(int(600-previous));start=time.monotonic()
    resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    p=Model(4);rows=[p.evolve(plan['horizon_seconds'],n,f'p4-{n}') for n in [16,32,64]]
    fields=['temperature','velocity','scalar'];comparisons={};accepted=True
    arrays=[dict(np.load(OUT/f'p4-{n}.npz')) for n in [16,32,64]]
    for key in fields:
        scale=float(max(abs(arrays[-1][key]).ravel()));err=float(max(abs(arrays[-1][key][::2]-arrays[-2][key]).ravel())/max(scale,1e-100))
        before=float(max(abs(arrays[-2][key][::2]-arrays[0][key]).ravel())/max(scale,1e-100))
        limit=plan['gates']['temperature_time' if key=='temperature' else 'GR_time']
        comparisons[key]=dict(time_relative=err,previous_relative=before,time_pass=err<limit,space_pass=False)
        accepted &= err<limit
    if accepted:
        del p
        other=Model(2);rows.append(other.evolve(plan['horizon_seconds'],64,'p2-64'));data=np.load(OUT/'p2-64.npz')
        for key in fields:
            scale=float(max(abs(arrays[-1][key]).ravel()));err=float(max(abs(arrays[-1][key]-data[key]).ravel())/max(scale,1e-100))
            limit=plan['gates']['temperature_space' if key=='temperature' else 'GR_space']
            comparisons[key].update(space_relative=err,space_pass=err<limit);accepted &= err<limit
    result=dict(classification='Counterexample candidate',passed=bool(accepted),actual_new_background_coupled_evolution=True,
        comparisons=comparisons,endpoint=rows[2]['history'][-1],seconds=time.monotonic()-start,total_compute_seconds=previous+time.monotonic()-start,
        memory_GB=max(x['memory_GB'] for x in rows),full_nonlinear_radiation_hydrodynamics=False,
        physical_zero_pressure_atmosphere=False,final_asymptotic_charge=False,full_goal_complete=False)
    write(OUT/'result.json',result);signal.alarm(0);print('RESULT',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','coefficients','pilot','run'])
    globals()[parser.parse_args().action]()
