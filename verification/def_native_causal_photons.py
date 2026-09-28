"""Counterexample candidate: native H, heat and spectral P1 on the causal interior.

Fixed gas density/metric and a tangent material response are explicit limits.
The photon Cauchy spectrum is specified; bolometric flux alone cannot select it.
"""
from pathlib import Path
import ctypes
import json
import signal
import sys
import time
import numpy as np
from scipy import sparse
from scipy.interpolate import PchipInterpolator
from scipy.sparse.linalg import splu
import def_native_hydrogen_exchange as old
import def_native_whole_star_match as star

OUT=old.OUT.parent/'def-native-causal-photons'
C,K,H=old.C,old.K,old.H
write=old.write
sha=old.old.cold.sha
END=.0034344311179287023
DEPTH=1.1e8


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='23694eddf',
        claim='Evolve spectral photon moments, native H inventory and material heat together throughout the scalar causal interior, including the retained core/envelope junction.',
        decision='Measure whether the deeper volumetric thermal/chemical source can be omitted from the later charge. This calculation does not itself complete the moving atmosphere or GR feedback.',
        model='Frozen physical density/metric; same-native finite H levels and thermodynamic reciprocal bound-free rates, elastic Thomson scattering, P1 angular moment closure. All non-H ionic inventories frozen separately at their local unbiased native populations. Tangent heat/H response is accepted only within its explicit small-change and native-defect gates.',
        radiation_initial='Local Planck mean occupation at saved gas T, spectral flux proportional to that Planck shape and normalized to the saved grey luminosity (native core diffusion below the junction). E and F are specified independently. At the free surface F=cE/2. This replaces the earlier diluted hotter Planck bath as an explicitly different Cauchy premise; neither spectrum is inferred uniquely from bolometric luminosity.',
        transport='Killing frequency a*h*nu; Jordan areal radius. J_t+(a*c/B)[H_r+2(1/r-a_r/a)H]=collision; H_t+(a*c/3B)J_r=-a*c*chi_transport*H. P1 is not full angle transport.',
        support=dict(depth_cm=DEPTH,horizon_seconds=END,radial_paths=[64,128],pilot_cells=16,frequency_orders=[8,12],time_steps=[64,128]),
        gates=dict(equilibrium=1e-10,native_derivative=.002,frequency_integral=.002,energy_balance=1e-8,space=.02,time=.02,temperature_change=.02,neutral_relative_change=.1,native_collision_defect=.02),
        budget=dict(native_calls=3500,native_seconds=120,transport_seconds=240,memory_GB=3,CPU_threads=1),
        forecast='Recent native constrained calls about0.01-0.04s each; actual deep states and sparse factors unmeasured. Measure16cells first; no production if the observed factor/run forecast exceeds240s.',
        stop='Preserve any failed path, no automatic extra resolution/time/domain, no final charge from a failed tangent response, no Rosseland opacity used as a spectral absorption law.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.__file__),star.OUT/'background.npz',star.OUT/'final-core.npz',star.OUT/'final-envelope.npz',old.ATOMIC,old.LEVELS]}))
    import sympy as s
    r,a,B,c=s.symbols('r a B c',positive=True)
    jr,hh,hr,ar=s.symbols('jr hh hr ar')
    w=r*r*B/a**3;area=r*r/a**2
    assert s.simplify(w*a*c/B*(hr+2*(1/r-ar/a)*hh)-c*(area*hr+(2*r/a**2-2*r*r*ar/a**3)*hh))==0
    q,energy,rho,nh,uy,cv=s.symbols('q energy rho nh uy cv',nonzero=True)
    yd=q/(rho*nh);td=(-energy*q/rho-uy*yd)/cv
    assert s.simplify(rho*(cv*td+uy*yd)+energy*q)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        photon_measure='W=r^2*B/a^3 gives d(W*J)/dt+c*d(r^2*H/a^2)/dr=W*collision at fixed Killing frequency.',
        exchange='Every retained bound-free net photon creation changes neutral H by one. rho*(cvT*theta_dot+u_y*y_dot) exactly cancels the photon energy source. These are identities of the declared moment/reaction model.'))


class Background:
    def __init__(self):
        d=self.d=dict(np.load(star.OUT/'background.npz'));self.env=np.load(star.OUT/'final-envelope.npz')
        self.core=np.load(star.OUT/'final-core.npz')
        self.re=d['radius_cm'];self.A=np.exp(-2*d['phi']**2);self.r=self.A*self.re
        self.R=self.r[-1];self.depth=self.R-self.r
        core_r=self.core['states'][:,0]*float(self.core['R_scale_m'])*100
        relevant=core_r>self.re[-1]-DEPTH*1.01
        assert np.max(abs(self.core['X'][relevant]-self.core['X'][0]))==0,'Causal core composition changes need explicit cells'
        self.X=self.core['X'][0]
        self.junction=float(self.env['r'][0]*self.env['A'][0])
        self.opacity=star.envelope.prior.two.radiative.tables.Opacity()
        self.theta=PchipInterpolator(self.r,self.A*d['lapse']*d['temperature_K'])

    def sample(self,r):
        assert min(r)>=self.R-DEPTH-1e-5 and max(r)<=self.R+1e-5
        d=self.d
        def get(key,log=False):
            v=d[key];return np.exp(np.interp(r,self.r,np.log(v))) if log else np.interp(r,self.r,v)
        re=np.interp(r,self.r,self.re);phi=get('phi');a=np.exp(-2*phi**2)*get('lapse')
        b=1-2*get('mass_geom_cm')/re;den=1-4*phi*re*get('phi_prime_cm')
        return dict(r=r,a=a,B=1/(np.sqrt(b)*den),rho=get('density_cgs',True),T=get('temperature_K',True),phi=phi)

    def luminosity(self,z):
        r,a,B,rho,T=[z[k] for k in ['r','a','B','rho','T']]
        L=np.full(len(r),float(self.env['Linfinity']));core=r<self.junction
        for j in np.flatnonzero(core):
            kap=star.envelope.prior.two.opacity_parts(self.opacity,(np.log(rho[j]),np.log(T[j]),self.X))[0]
            arad=float(self.env['Prad'][0]*3/self.env['T'][0]**4)
            F=-4*arad*C*T[j]**3/(3*rho[j]*kap*a[j]*B[j])*self.theta(r[j],1)
            L[j]=4*np.pi*r[j]**2*a[j]**2*F
        return L


def frequencies(native,order):
    # Resolve every retained threshold. Local redshift is kept in evaluation;
    # the tiny shift of thresholds within a bin is tested by the second order.
    a=Background().sample(np.array([Background().R]))['a'][0]
    edges=np.unique(np.r_[0.,native.binding*a,K*np.array([2e4,4e4,8e4])*a*np.array([1,1,1]),K*8e4*a*np.array([2,4,8,16,32,60])])
    x,w=np.polynomial.legendre.leggauss(order)
    E=((edges[:-1]+edges[1:])[:,None]/2+np.diff(edges)[:,None]*x/2).ravel()
    dE=(np.diff(edges)[:,None]*w/2).ravel()
    return E,8*np.pi*E**2*dE/(H**3*C**3)


def coefficients(native,state,energy):
    T=np.exp(state['lt']);sigma=np.zeros_like(energy)
    for n,f in enumerate(state['fraction'],1):
        good=energy>=native.binding[n-1]
        if not good.any() or f==0:continue
        nu=np.ascontiguousarray((native.oldground/n**2+energy[good]-native.binding[n-1])/float(np.float32(6.6256e-27)))
        bf=np.zeros_like(nu);ff=np.zeros_like(nu);native.cross(len(nu),1,n,float(T),nu,bf,ff)
        sigma[good]+=f*bf
    absorb=state['y']*sigma
    emit=absorb*np.exp(np.clip(state['affinity']-energy/(K*T),-745,700))
    return np.array([absorb-emit,emit])


def bank(n=128,order=8,label=None):
    label=label or f'bank-{n}-{order}';assert not (OUT/(label+'.npz')).exists()
    start=time.monotonic();signal.alarm(120);bg=Background();native=old.Native(cap=3500)
    assert np.array_equal(bg.X,native.fan.X)
    edges=np.linspace(bg.R-DEPTH,bg.R,n+1);r=(edges[:-1]+edges[1:])/2;z=bg.sample(r);zf=bg.sample(edges)
    Einf,num=frequencies(native,order);raws=[];const=[];thermo=[];ys=[];snapshots=[];derivative=[]
    for j,(rho,T,a) in enumerate(zip(z['rho'],z['T'],z['a'])):
        lt=np.log(T);eq=native.ion.snapshot(np.log(rho),lt,np.zeros(318));target=eq['number_fractions']
        native.lr=np.log(rho);native.base=eq;native.target=target;native.epsH=float(target[0,:2].sum());native.y0=float(target[0,0]/native.epsH)
        native.nH=native.epsH*native.fan.cx*old.NA
        native.active=np.concatenate([target[e,1:v+1]>target[e].sum()*1e-18 for e,v in enumerate(native.ion.Z)])
        native.prefix=dict(T=np.array([T]),fields=np.zeros((1,318)),log_density_ratio=np.zeros(1))
        y=native.y0;h=1e-4
        states=[native.state(0.,lt,y),native.state(0.,lt+h,y),native.state(0.,lt-h,y),native.state(0.,lt,y*(1+h)),native.state(0.,lt,y*(1-h))]
        rows=np.array([coefficients(native,s,Einf/a) for s in states]);rs=np.array([s['raw'] for s in states])
        bb=1/np.expm1(Einf/(a*K*T));c0=rows[0];ct=(rows[1]-rows[2])/(2*h);cy=(rows[3]-rows[4])/(2*h)
        equilibrium=np.max(abs(c0[1]-c0[0]*bb)/(abs(c0[1])+abs(c0[0]*bb)+1e-300));assert equilibrium<1e-10
        const.append([c0[0],ct[1]-ct[0]*bb,cy[1]-cy[0]*bb,ct[0],cy[0]])
        dt=(rs[1]-rs[2])/(2*h);dy=(rs[3]-rs[4])/(2*h)
        # Derivatives dy are with respect to relative neutral change eta=dy/y0.
        thermo.append([dt[2],dy[2],dt[1],dy[1],native.nH,rs[0,13]/1.66053906660e-24,dt[13]/1.66053906660e-24,dy[13]/1.66053906660e-24])
        derivative.append(abs(dt[2]/(T*dt[3])-1));ys.append(y);raws.append(rs[0]);snapshots.append(target)
        assert dt[2]>0
    # Thomson electron count uses the declared atomic mass constant, as before.
    th=np.array(thermo);th[:,5]=np.array(raws)[:,13]/1.66053906660e-24
    np.savez_compressed(OUT/(label+'.npz'),edges=edges,**z,face_a=zf['a'],face_B=zf['B'],face_T=zf['T'],face_L=bg.luminosity(zf),
        Einf=Einf,num=num,coeff=np.array(const),thermo=th,raw=np.array(raws),y0=ys,target=np.array(snapshots),cx=native.fan.cx,native_calls=native.ion.calls,derivative_errors=derivative)
    elapsed=time.monotonic()-start
    result=dict(classification='Counterexample candidate',passed=bool(max(derivative)<.002),native_derivative_relative=float(max(derivative)),native_calls=native.ion.calls,seconds=elapsed,n=n,order=order,
        same_causal_composition=True,source_sha256=sha(__file__))
    write(OUT/(label+'.json'),result);signal.alarm(0);print(json.dumps(result),flush=True);assert result['passed']


class Model:
    def __init__(self,label):
        self.d=d=np.load(OUT/(label+'.npz'));n=len(d['r']);m=len(d['Einf']);self.n=n;self.m=m;self.nj=n*m;self.nh=(n+1)*m;self.size=self.nj+self.nh+2*n
        self.it=self.nj+self.nh+np.arange(n)*2;self.iy=self.it+1
        r,a,B,rho,T=[d[k] for k in ['r','a','B','rho','T']];edges=d['edges'];dx=np.diff(edges)
        self.J=1/np.expm1(d['Einf'][None,:]/(a*K*T)[:,None]);self.scale=self.J.max(0)
        jf=1/np.expm1(d['Einf'][None,:]/(d['face_a']*K*d['face_T'])[:,None])
        ef=(jf*d['num']*d['Einf']).sum(1)/d['face_a']**4
        ff=d['face_L']/(4*np.pi*edges**2*d['face_a']**2*C*ef)
        # With K=J/3, positive moments require H^2<=J*K. The half-isotropic
        # surface value1/2 is not the maximum allowed interior moment ratio.
        assert max(abs(ff))<1/np.sqrt(3),('Initial angular moment realizability',min(ff),max(ff))
        self.H=jf*ff[:,None];self.volume=4*np.pi*r*r*B*dx;W=r*r*B/a**3*dx;area=edges**2/d['face_a']**2
        self.energyJ=(4*np.pi*W[:,None]*d['num']*d['Einf']*self.scale).ravel()
        self.gasT=self.volume*rho*d['thermo'][:,0]*a;self.gasY=self.volume*rho*d['thermo'][:,1]*a
        I=np.arange(self.nj).reshape(n,m);F=self.nj+np.arange(self.nh).reshape(n+1,m)
        rows=[];cols=[];vals=[];b=np.zeros(self.size)
        def add(row,col,val):
            rr,cc,vv=np.broadcast_arrays(row,col,val);rows.extend(rr.ravel());cols.extend(cc.ravel());vals.extend(vv.ravel())
        cl,ct,cy,clt,cly=np.moveaxis(d['coeff'],1,0)
        factor=a*C*rho*d['thermo'][:,4]
        cj=-factor[:,None]*cl;ct=factor[:,None]*ct/self.scale;cy=factor[:,None]*cy/self.scale
        add(I,I,cj);add(I,self.it[:,None],ct);add(I,self.iy[:,None],cy)
        nj=d['num'][None,:]*self.scale/(a**3*rho*d['thermo'][:,4]*d['y0'])[:,None]
        et=-d['num'][None,:]*d['Einf']*self.scale/(a**4*rho*d['thermo'][:,0])[:,None]-nj*(d['thermo'][:,1]/d['thermo'][:,0])[:,None]
        for idx,w in [(self.it,et),(self.iy,nj)]:
            add(idx[:,None],I,w*cj);add(idx, self.it,(w*ct).sum(1));add(idx,self.iy,(w*cy).sum(1))
        add(I,F[:-1],C*area[:-1,None]/W[:,None]);add(I,F[1:],-C*area[1:,None]/W[:,None])
        h0=self.H/self.scale;b[I]=C/W[:,None]*(area[:-1,None]*h0[:-1]-area[1:,None]*h0[1:])
        # Face moment equation uses a centered gradient and native extinction.
        dist=np.diff(r);speed=d['face_a'][1:-1]*C/(3*d['face_B'][1:-1]*dist)
        add(F[1:-1],I[:-1],speed[:,None]);add(F[1:-1],I[1:],-speed[:,None])
        chi=rho[:,None]*d['thermo'][:,4,None]*cl+d['thermo'][:,5,None]*6.6524587321e-25
        rate=d['face_a'][1:-1,None]*C*(chi[:-1]+chi[1:])/2
        add(F[1:-1],F[1:-1],-rate)
        b[F[1:-1]]=speed[:,None]*(self.J[:-1]-self.J[1:])/self.scale-rate*h0[1:-1]
        for idx,der,necol in [(self.it,clt,6),(self.iy,cly,7)]:
            for shift in [0,1]:
                sel=slice(shift,n-1+shift)
                extinction=(rho*d['thermo'][:,4])[sel,None]*der[sel]+d['thermo'][sel,necol,None]*6.6524587321e-25
                val=-d['face_a'][1:-1,None]*C*extinction/2*h0[1:-1]
                add(F[1:-1],idx[sel,None],val)
        # Fixed incoming flux at the causally remote inner face. The outer
        # Marshak perturbation H=J/2 is an algebraic eliminated face flux.
        self.outer=F[-1];self.inner=F[0]
        A=sparse.coo_matrix((vals,(rows,cols)),shape=(self.size,self.size)).tocsc()
        # Outer H is copied from its algebraic value at every step by replacing
        # that row in the time-discrete system; its evolution row stays zero.
        self.A=A;self.b=b;self.I=I;self.F=F
        self.initial_flux_fraction=ff

    def run(self,steps,label):
        start=time.monotonic();h=END/steps;A=self.A;n=self.n;d=self.d
        left=(sparse.eye(self.size,format='csc')-h*A/2).tolil();right=(sparse.eye(self.size,format='csc')+h*A/2).tocsr()
        for f,j in zip(self.outer,self.I[-1]):left.rows[f]=[int(j),int(f)];left.data[f]=[-.5,1.]
        lu=splu(left.tocsc());factor_seconds=time.monotonic()-start;state=np.zeros(self.size);times=[0.];history=[state.copy()];flux_energy=0.;balance=0.
        port=4*np.pi*C*d['edges']**2/d['face_a']**2
        q=d['num']*d['Einf']*self.scale
        def energy(s):return float(self.energyJ@s[:self.nj]+self.gasT@s[self.it]+self.gasY@s[self.iy])
        def power(s):return float(port[0]*(q@(s[self.inner]+self.H[0]/self.scale))-port[-1]*(q@(s[self.outer]+self.H[-1]/self.scale)))
        # The initial photons already contain finite energy. Only changes and
        # net face power enter this ledger, avoiding rest/background subtraction.
        for k in range(steps):
            rhs=right@state+h*self.b;rhs[self.outer]=0.;new=lu.solve(rhs)
            flux_energy+=h*(power(state)+power(new))/2
            balance=max(balance,abs(energy(new)-flux_energy))
            state=new;times.append((k+1)*h);history.append(state.copy())
        states=np.array(history);theta=states[:,self.it];eta=states[:,self.iy]
        trace=(d['rho']*d['thermo'][:,0]-3*d['thermo'][:,2])*theta+(d['rho']*d['thermo'][:,1]-3*d['thermo'][:,3])*eta
        trace_energy=trace@(self.volume*d['a'])
        scale=max(abs(flux_energy),np.max(abs(states[:,:self.nj]@self.energyJ)),np.max(abs(trace_energy)),1.)
        J=states[:,:self.nj].reshape(-1,n,self.m)*self.scale+self.J
        Hf=states[:,self.nj:self.nj+self.nh].reshape(-1,n+1,self.m)*self.scale+self.H
        np.savez_compressed(OUT/(label+'.npz'),times=times,theta=theta,eta=eta,J=J,H=Hf,trace_density=trace,trace_energy=trace_energy,r=d['r'],volume=self.volume)
        result=dict(classification='Counterexample candidate',passed=bool(balance/scale<1e-8 and np.max(abs(theta))<.02 and np.max(abs(eta))<.1 and J.min()>=-1e-10*self.J.max()),
            cells=n,steps=steps,frequencies=self.m,seconds=time.monotonic()-start,factor_seconds=factor_seconds,
            energy_balance_relative=float(balance/scale),boundary_energy_erg=flux_energy,trace_endpoint_erg=float(trace_energy[-1]),
            maximum_logT_change=float(np.max(abs(theta))),maximum_relative_neutral_change=float(np.max(abs(eta))),minimum_photon_occupation=float(J.min()),
            tangent_material=True,full_angular_transport=False,gas_motion_evolved=False,full_GR_feedback=False,final_charge_solved=False,source_sha256=sha(__file__))
        write(OUT/(label+'.json'),result);print(json.dumps(result),flush=True);return result


def pilot():
    signal.alarm(60)
    if not (OUT/'bank-16-8.npz').exists():bank(16,8)
    m=Model('bank-16-8');result=m.run(32,'pilot')
    forecast=result['seconds']*64*3
    write(OUT/'measured-budget.json',dict(pilot_seconds=result['seconds'],forecast_seconds=forecast,assumption='Eightfold radial size, conservative quadratic factor/run scaling, three paths. Deep native preparation is accounted separately.',remaining_transport_seconds=240-result['seconds'],eligible=forecast<240-result['seconds']))
    signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]](*map(int,sys.argv[2:]))
