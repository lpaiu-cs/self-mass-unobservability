"""Counterexample candidate: conservative, forced native interior mechanics.

The saved thermochemical/radiation history drives baryon and momentum motion.
The density/temperature response is linear and feeds mechanical pressure back;
its radiation feedback and spatial continuum error remain separate questions.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
import def_native_coupled_charge as prior
import def_native_causal_nonlinear as chemistry

OUT=prior.OUT.parent/'def-native-interior-motion'
C=prior.C;G=prior.G;write=prior.write;sha=prior.sha


def prepare():
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True)
    paths=[Path(__file__),Path(prior.__file__),Path(chemistry.__file__),
           prior.OUT/'source-896-128.npz',prior.OUT/'exterior-energy.npz',
           prior.cold.OUT/'resumed-896-128.npz',prior.prior.OUT/'cells-896-steps-128.npz']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d7708d9a2',
        claim='Release the frozen deep baryons using the actual saved pressure and photon collision force, including native fixed-inventory adiabatic pressure feedback; apply the resulting mass redistribution to the same outgoing scalar component.',
        decision='Determine whether neglecting interior motion changes the sign or size of the Phase114 component. A controlled time solve is not a spatial or full-GR certificate.',
        reuse='Same16 deep volumes, actual448/128 and896/128 histories, initial metric, local inventories,152 frequencies,8 angles and3.434ms. No repeated heat/photon/fluid production trajectory.',
        model='Staggered conservative displacement-mass flux on the existing16 volumes. Delta mass=-difference(face transported mass). Linear Lagrangian adiabatic response with native constrained dP,dU. Evolve interior face momentum and compression pressure together. Actual cumulative atmospheric baryon influx imposes the outer face displacement; inner displacement is zero outside the observed light cone.',
        forcing='Actual H bound-free plus Thomson collision first angular moment is given to gas. Initial gas/gravity imbalance comes from the saved total hydrostatic pressure and initial LTE photon gradient; it is retained, not subtracted. The actual evolving thermal pressure adds its radial force.',
        limitations='One-way saved radiation drive, linear mechanics,16 original volumes, prescribed mechanical outer face and fixed metric. Neither a tiny displacement nor time convergence proves a final charge. Non-H inventories advect with material only in this first-order constitutive closure.',
        native=dict(points='all16 initial and final states',log_step=.0001,half_step_cells=[0,10,15],call_cap=2000,seconds=60),
        paths=dict(time_steps=[64,128],controls=['initial-only stiffness','PCHIP versus centered background gradients','original448/128 forcing']),
        gates=dict(native_derivative=.002,native_anchor=.002,time_charge=.02,time_motion=.002,linear_density=.001,linear_speed_over_c=.00001,baryon_ledger=1e-12),
        budget=dict(source_seconds=30,native_seconds=60,motion_seconds=45,readout_seconds=45,CPU_threads=1,memory_GB=2,new_whole_star_steps=0),
        forecast='Previously measured native calls0.01-0.04s, at most2000 including internal constraint iterations.32 mechanical degrees of freedom; time solve expected below1s but unmeasured. Hard action alarms; no automatic grid/horizon/path expansion.',
        stop='Preserve any gate failure. No gate relaxation; a material component change requires revising the physical conclusion. Do not declare the full goal complete.',
        bindings={str(p):sha(p) for p in paths}))


def model():
    return prior.prior.Coupled(896,8)


def forcing():
    assert not (OUT/'forcing.json').exists();start=time.monotonic();signal.alarm(30)
    m=model();b=m.bulk;d=b.d;geo=prior.green.Geometry(m.m)
    bg=chemistry.prior.Background();r=d['r'];edges=d['edges'];_,af,Bf,ref,_=geo.metric(edges-m.m.RJ)
    # The saved whole-star temperature defines the initial radiation support.
    lt=PchipInterpolator(bg.r,np.log(bg.d['temperature_K']))
    Prad=float(bg.env['Prad'][0]*3/bg.env['T'][0]**4)*d['T']**4/3
    _,a,B,re,phi=geo.metric(r-m.m.RJ);pt=m.m.bg.sample(re)
    A=np.exp(-2*phi**2);bb=1-2*pt['m']/re;alpha=-4*phi;den=1+alpha*re*pt['v']
    ap=a*(pt['m']/(re*re*bb)+4*np.pi*re*A**4*pt['p']/bb+re*pt['v']**2/2+alpha*pt['v'])/(m.m.R*A*den)
    p0,u0,*_=b.eos.gas(np.zeros(b.n),np.zeros(b.n));cx=float(d['cx'])
    pgprime=-(d['rho']*(cx*C*C+u0)+p0+4*Prad)*ap/a-4*Prad*lt(r,1)
    rho_prime=d['rho']*PchipInterpolator(bg.r,np.log(bg.d['density_cgs']))(r,1)
    # u and inventory background gradients use the actual retained native nodes.
    uprime=PchipInterpolator(r,u0)(r,1)
    rho_face=bg.sample(edges)['rho']
    initial_support=a/B*(4*Prad*lt(r,1)+4*Prad*ap/a)
    rows=[]
    for cells in [896,448]:
        z=prior.load_path(cells,128);times=np.r_[0,z['t']]
        theta=np.vstack([np.zeros(b.n),z['theta']]);eta=np.vstack([np.zeros(b.n),z['eta']])
        photons=np.concatenate([b.initial[None],z['bulk_I']]);press=[];energy=[];force=[]
        for th,et,I in zip(theta,eta,photons):
            p,u,*_=b.eos.gas(th,et);ne=b.eos.gas(th,et)[6]
            ab,em,*_=b.eos.radiation(th,et)
            collision=b.factor[:,None,None]*(em[:,None,:]-(ab-em)[:,None,:]*I)
            mean=b.mean(I);p2=b.mean(I*b.P2[None,:,None])
            collision+=(a*C*ne*6.6524587321e-25)[:,None,None]*(mean[:,None,:]+.5*b.P2[None,:,None]*p2[:,None,:]-I)
            f=-np.einsum('iqf,q,f->i',collision,b.w*b.mu,d['num']*d['Einf'])/(a**4*C)
            press.append(p-p0);energy.append(d['rho']*(u-u0));force.append(f)
        source=np.load(prior.OUT/f'source-{cells}-128.npz');outer=-np.asarray(source['baryon_g'],float).sum(1)
        np.savez_compressed(OUT/f'forcing-{cells}.npz',t=times,r=r,edges=edges,volume=b.volume,a=a,B=B,ap=ap,af=af,Bf=Bf,
            rho=d['rho'],rho_face=rho_face,p0=p0,u0=u0,cx=cx,pgprime=pgprime,rho_prime=rho_prime,uprime=uprime,
            initial_support=initial_support,dp=press,denergy=energy,radiation_force=force,outer_mass=outer,theta=theta,eta=eta,
            RJ=m.m.RJ,M_cm=float(source['M_cm']))
        rows.append(dict(cells=cells,maximum_thermal_pressure_relative=float(np.max(abs(np.array(press)/p0))),
                         maximum_force_cgs=float(np.max(abs(force))),initial_support_force_cgs=initial_support.tolist()))
    write(OUT/'forcing.json',dict(classification='Counterexample candidate',passed=True,paths=rows,seconds=time.monotonic()-start))
    signal.alarm(0);print(json.dumps(rows),flush=True)


def native():
    assert not (OUT/'native.json').exists();start=time.monotonic();signal.alarm(60)
    b=model().bulk;d=b.d;z=np.load(OUT/'forcing-896.npz');n=chemistry.old.Native(cap=2000)
    h=1e-4;raw=np.zeros((2,b.n,5,21));checks=[];anchors=[]
    for j in range(b.n):
        chemistry.setup(n,d,j)
        for ti,it in enumerate([0,-1]):
            lt=np.log(d['T'][j])+z['theta'][it,j];y=d['y0'][j]*(1+z['eta'][it,j])
            states=[n.state(x,lt+dt,y)['raw'] for x,dt in [(0,0),(h,0),(-h,0),(0,h),(0,-h)]]
            raw[ti,j]=states;rr=(states[1]-states[2])/(2*h);tt=(states[3]-states[4])/(2*h)
            p,u=states[0][[1,2]];ut=tt[2];pr=rr[1];pt=tt[1];ur=rr[2]
            K=pr+pt*(p/d['rho'][j]-ur)/ut
            assert ut>0 and K>0
            exact=b.eos.gas(z['theta'][it],z['eta'][it]);anchors += [abs(p/exact[0][j]-1),abs(u/exact[1][j]-1)]
            thermo=max(abs(ut/(np.exp(lt)*tt[3])-1),abs((ur-p/d['rho'][j]-np.exp(lt)*rr[3])/(abs(ur)+p/d['rho'][j])))
            row=dict(cell=j,time_index=it,K=K,gamma_fixed_inventory=K/p,thermodynamic_relative=float(thermo))
            if j in [0,10,15]:
                half=[n.state(x,lt+dt,y)['raw'] for x,dt in [(h/2,0),(-h/2,0),(0,h/2),(0,-h/2)]]
                rh=(half[0]-half[1])/h;th=(half[2]-half[3])/h;Kh=rh[1]+th[1]*(p/d['rho'][j]-rh[2])/th[2]
                row['half_step_K_relative']=abs(Kh/K-1)
            checks.append(row)
    pr=(raw[:,:,1,1]-raw[:,:,2,1])/(2*h);ur=(raw[:,:,1,2]-raw[:,:,2,2])/(2*h)
    pt=(raw[:,:,3,1]-raw[:,:,4,1])/(2*h);ut=(raw[:,:,3,2]-raw[:,:,4,2])/(2*h)
    K=pr+pt*(raw[:,:,0,1]/d['rho']-ur)/ut
    np.savez_compressed(OUT/'native.npz',raw=raw,K=K,pr=pr,ur=ur,pt=pt,ut=ut,native_calls=n.ion.calls)
    passed=max(anchors)<.002 and max(x['thermodynamic_relative'] for x in checks)<.002 and max(x.get('half_step_K_relative',0) for x in checks)<.002
    row=dict(classification='Counterexample candidate',passed=bool(passed),checks=checks,anchor_relative=max(anchors),native_calls=n.ion.calls,seconds=time.monotonic()-start)
    write(OUT/'native.json',row);signal.alarm(0);print(json.dumps(row),flush=True);assert passed


class Mechanics:
    def __init__(self,cells=896,initial_K=False,gradient=False):
        self.d=d=dict(np.load(OUT/f'forcing-{cells}.npz'));self.n=len(d['r']);self.K=np.load(OUT/'native.npz')['K'];self.initial_K=initial_K
        r=d['r'];V=d['volume'];rho=d['rho'];edges=d['edges'];n=self.n
        self.P=prior.green.polynomial(d['t'],np.c_[d['dp'],d['denergy'],d['radiation_force']], 'linear')
        self.boundary=prior.green.polynomial(d['t'],d['outer_mass'],'linear')
        self.mass=-np.diff(np.eye(n+1),axis=0)
        # Maps integrated face baryon flux to Eulerian density and displacement.
        self.drho=self.mass/V[:,None]
        average=np.zeros((n,n+1));average[np.arange(n),np.arange(n)]=.5;average[np.arange(n),np.arange(1,n+1)]=.5
        self.xi=average/(4*np.pi*r*r*d['B']*rho)[:,None]
        if gradient:
            d['rho_prime']=np.gradient(rho,r);d['uprime']=np.gradient(d['u0'],r)
        self.compression=self.drho/rho[:,None]+(d['rho_prime']/rho)[:,None]*self.xi
        self.de=(d['cx']*C*C+d['u0'])[:,None]*self.drho+d['p0'][:,None]*self.compression-(rho*d['uprime'])[:,None]*self.xi
        self.du=d['u0'][:,None]*self.drho+d['p0'][:,None]*self.compression-(rho*d['uprime'])[:,None]*self.xi
        self.face_average=np.zeros((n-1,n))
        frac=(edges[1:-1]-r[:-1])/np.diff(r)
        self.face_average[np.arange(n-1),np.arange(n-1)]=1-frac
        self.face_average[np.arange(n-1),np.arange(1,n)]=frac
        self.grad=np.diff(np.eye(n),axis=0)/np.diff(r)[:,None]
        self.inertia=d['rho']*d['cx']+(d['rho']*d['u0']+d['p0'])/C**2
        inertia_face=self.face_average@self.inertia
        self.pref=4*np.pi*edges[1:-1]**2*d['af'][1:-1]*d['rho_face'][1:-1]/inertia_face
        self.gravity=self.face_average@(d['ap']/d['B'])
        self.force0=self.face_average@d['initial_support']

    def pressure(self,t):
        K=self.K[0] if self.initial_K else self.K[0]+t/self.d['t'][-1]*(self.K[1]-self.K[0])
        return K[:,None]*self.compression-self.d['pgprime'][:,None]*self.xi

    def operator(self,t):
        d=self.d;n=self.n;pressure=self.pressure(t);p,e,fr=np.split(self.P(t),3)
        L=-d['af'][1:-1,None]/d['Bf'][1:-1,None]*(self.grad@pressure)-self.gravity[:,None]*(self.face_average@(self.de+pressure))
        f=-d['af'][1:-1]/d['Bf'][1:-1]*(self.grad@p)-self.gravity*(self.face_average@(e+p))+self.face_average@fr+self.force0
        return self.pref[:,None]*L,self.pref*f

    def solve(self,steps,label):
        start=time.monotonic();d=self.d;n=self.n;h=d['t'][-1]/steps
        y=np.zeros(2*(n-1));states=[];ts=np.linspace(0,d['t'][-1],steps+1);eye=np.eye(len(y));residual=0.
        for i,t in enumerate(ts):
            mass=np.r_[0,y[:n-1],float(self.boundary(t))];current=np.r_[0,y[n-1:],float(self.boundary(t,1))]
            states.append((mass,current))
            if i==steps:break
            tm=t+h/2;L,f=self.operator(tm);force=f+L[:,-1]*float(self.boundary(tm))
            A=np.block([[np.zeros((n-1,n-1)),np.eye(n-1)],[L[:,1:-1],np.zeros((n-1,n-1))]])
            rhs=(eye+h*A/2)@y+h*np.r_[np.zeros(n-1),force];mat=eye-h*A/2
            yy=np.linalg.solve(mat,rhs);residual=max(residual,float(np.max(abs(mat@yy-rhs))/(np.max(abs(rhs))+1e-100)));y=yy
        transported=np.array([x[0] for x in states]);current=np.array([x[1] for x in states])
        mass=transported@self.mass.T;xi=transported@self.xi.T;rho=transported@self.drho.T
        compression=transported@self.compression.T;pressure=np.array([self.pressure(t)@z for t,z in zip(ts,transported)])
        internal=transported@self.du.T;velocity=current/(4*np.pi*d['edges']**2*d['af']*d['rho_face'])
        net=-transported[:,-1];ledger=float(np.max(abs(mass.sum(1)-net))/max(np.max(abs(mass).sum(1)),1.))
        trace=(internal-3*pressure)*d['volume'];energy=internal*d['volume']
        np.savez_compressed(OUT/(label+'.npz'),t=ts,transported_g=transported,current_g_s=current,baryon_g=mass,
            displacement_cm=xi,density_relative=rho/d['rho'],lagrangian_density_relative=compression,
            velocity_cm_s=velocity,pressure=pressure,nonrest_trace_erg=trace,internal_energy_erg=energy)
        row=dict(classification='Counterexample candidate',steps=steps,seconds=time.monotonic()-start,linear_solve_relative=residual,
            baryon_ledger_relative=ledger,maximum_density_relative=float(np.max(abs(rho/d['rho']))),maximum_lagrangian_compression=float(np.max(abs(compression))),
            maximum_displacement_cm=float(np.max(abs(xi))),maximum_velocity_over_c=float(np.max(abs(velocity))/C),
            maximum_mechanical_pressure_over_initial=float(np.max(abs(pressure)/d['p0'])),maximum_internal_baryon_redistribution_g=float(np.max(abs(mass))),
            endpoint_deep_net_mass_g=float(net[-1]),full_two_way_radiation=False,spatial_continuum_certified=False,final_charge_solved=False)
        write(OUT/(label+'.json'),row);print(label,json.dumps(row),flush=True);return row


def motion():
    assert json.loads((OUT/'native.json').read_text())['passed'];assert not (OUT/'motion.json').exists();start=time.monotonic();signal.alarm(45)
    rows=[]
    for label,cells,steps,ki,grad in [('motion-64',896,64,False,False),('motion-128',896,128,False,False),
            ('initial-stiffness',896,128,True,False),('background-gradient',896,128,False,True),('coarse-forcing',448,128,False,False)]:
        rows.append(Mechanics(cells,ki,grad).solve(steps,label))
    z=np.load(OUT/'motion-128.npz');v=np.load(OUT/'motion-64.npz');err=float(np.max(abs(z['transported_g'][::2]-v['transported_g']))/np.max(abs(z['transported_g'])))
    passed=err<.002 and all(r['maximum_density_relative']<.001 and r['maximum_velocity_over_c']<1e-5 and r['baryon_ledger_relative']<1e-12 for r in rows)
    write(OUT/'motion.json',dict(classification='Counterexample candidate',passed=bool(passed),time_motion_relative=err,paths=rows,seconds=time.monotonic()-start));signal.alarm(0);assert passed


def component(z,m,order=8):
    d=np.load(OUT/'forcing-896.npz');xg,wg=np.polynomial.legendre.leggauss(order);edges=d['edges'];r=(edges[:-1,None]+np.diff(edges)[:,None]*(xg+1)/2).ravel()
    weights=np.tile(wg/2,len(edges)-1);ids=np.repeat(np.arange(len(edges)-1),order);geo=prior.green.Geometry(m.m)
    w,delay,*_=geo(r-m.m.RJ);observer=np.load(prior.OUT/'wave-896-128.npz')['t'];answer=[]
    for key,fac in [('baryon_g',G*float(d['cx'])/(2*C*float(d['M_cm']))),('nonrest_trace_erg',G/(2*C**3*float(d['M_cm'])))]:
        H=prior.green.polynomial(z['t'],z[key]).antiderivative();rows=[]
        for u in observer:
            at=u+delay;cut=np.clip(at,0,H.x[-1]);j=np.clip(np.searchsorted(H.x,cut,side='right')-1,0,len(H.x)-2);dt=cut-H.x[j];v=np.zeros(len(dt))
            for coeff in H.c:v=v*dt+coeff[j,ids]
            v[at<=0]=0;rows.append(float(np.sum(v*w*weights,dtype=np.longdouble))*fac)
        answer.append(rows)
    return observer,np.array(answer)


def readout():
    assert json.loads((OUT/'motion.json').read_text())['passed'];assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(45)
    m=model();curves={};parts={}
    for label in ['motion-64','motion-128','initial-stiffness','background-gradient','coarse-forcing']:
        t,p=component(np.load(OUT/(label+'.npz')),m);parts[label]=p;curves[label]=p.sum(0)
    baseline=np.load(prior.OUT/'wave-896-128.npz');original=baseline['normalized_direct'];debit=baseline['components'][-1]
    q=original-debit+curves['motion-128'];scale=max(abs(q));controls={label:float(max(abs(curves[label]-curves['motion-128']))/scale) for label in curves if label!='motion-128'}
    _,p4=component(np.load(OUT/'motion-128.npz'),m,4);controls['radial_quadrature']=float(max(abs(p4.sum(0)-curves['motion-128']))/scale)
    ext=np.load(prior.OUT/'exterior-energy.npz');print('EXTERIOR_KEYS',list(ext),flush=True)
    np.savez_compressed(OUT/'wave.npz',t=t,direct_original=original,direct_with_interior=q,interior_rest=parts['motion-128'][0],interior_nonrest=parts['motion-128'][1],removed_point_debit=debit)
    passed=controls['motion-64']<.02
    row=dict(classification='Counterexample candidate',passed=bool(passed),endpoint_direct_original=float(original[-1]),endpoint_direct_with_interior=float(q[-1]),
        endpoint_interior_rest=float(parts['motion-128'][0,-1]),endpoint_interior_nonrest=float(parts['motion-128'][1,-1]),removed_point_debit=float(debit[-1]),
        controls=controls,seconds=time.monotonic()-start,actual_internal_momentum_and_baryons_evolved=True,
        actual_adiabatic_pressure_feedback=True,one_way_saved_photon_drive=True,spatial_continuum_certified=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',row);signal.alarm(0);print(json.dumps(row),flush=True);assert passed


if __name__=='__main__':
    globals()[sys.argv[1]]()
