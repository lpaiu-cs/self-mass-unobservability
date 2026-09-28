"""Replace the supported outer gas by conservative radial vacuum expansion.

Counterexample candidate: native equilibrium chemistry, advected gas entropy,
the saved Jordan metric, and elastic electron scattering. Absorption, chemistry
kinetics and the full dynamical scalar/metric response are separate closures.
"""
from pathlib import Path
import argparse
import json
import signal
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
import def_native_vacuum_release as prior

old=prior.old
OUT=prior.OUT.parent/'def-native-radial-release'
write=old.write
C=old.C


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='c71b31a34',
        claim='Replace, rather than superpose upon, the gas in the outer200m by nonlinear baryon/momentum/entropy conservation through the vacuum boundary. Couple actual electron scattering momentum and its equal photon work debit.',
        initial='Native gas EOS at the saved nonuniform envelope rho,T,X. Subtract LTE photons; the unchanged stellar photon field is independent of the gas temperature. Inner boundary reuses the saved coupled velocity history; it is outside the inward gas acoustic cone over the same3.434ms.',
        model='Spherical Valencia gas on the saved Jordan metric, exact Lorentz stresses, HLL MUSCL SSPRK2, advected gas entropy. Electron elastic Thomson force uses native rmue=rho/mu_e. Current photon work updates the outward luminosity deviation through the layer. This optically thin instantaneous radiation correction is not absorption or full retarded scattering transport.',
        EOS=dict(log_density_ratio=[-18,.02],entropy_planes=3,reuse='All145 surface isentrope states from Phase98',native_calls=1200,seconds=60),
        paths=['224-cell radial pilot to full registered horizon','224-cell flat native-fan control','896 cells','1792 cells'],
        domain_m=[-200,1200],CFL=.35,
        gates=dict(baryon_ledger_relative=1e-10,energy_response_relative=.02,outside_mass_refinement=.02,integrated_stress_refinement=.02,EOS_control_relative=.002,flat_native_fan_mass=.03),
        budget=dict(production_seconds=240,CPU_threads=1,memory_GB=2,new_whole_star_evolutions=0),
        stop='No unregistered grid, horizon, native domain or threshold expansion. Stop on positivity, EOS domain, budget or conservation failure. Preserve every failed candidate.',
        boundaries='This replaces only the gas layer. Bulk thermal/GR feedback and a final outgoing normalized scalar charge are not inferred from gas convergence.',
        references=['https://www.uv.es/astrorela/simulacionnumerica/node13.html','https://academic.oup.com/mnras/article/417/4/2899/1098574'],
        bindings={str(p.relative_to(old.ROOT)):old.photons.digest(p) for p in [Path(__file__),prior.OUT/'fine.npz',old.prior.OUT/'final-envelope.npz',old.prior.OUT/'final-core.npz',prior.prior.OUT/'p4-64.npz']}))


def bank():
    assert not (OUT/'eos.npz').exists();start=time.monotonic();signal.alarm(60)
    fan=prior.Fan(call_cap=1200,reuse=True);env=fan.env;R=float(env['r'][-1]);A=float(env['A'][-1]);b=float(env['b'][-1])
    sample_r=R+np.linspace(-20000,0,17)*np.sqrt(b)/A
    rho=np.exp(np.interp(sample_r,env['r'],np.log(env['rho'])));T=np.exp(np.interp(sample_r,env['r'],np.log(env['T'])))
    initial=np.array([fan.call(np.log(r),np.log(t)) for r,t in zip(rho,T)])
    s0=fan.entropy;unit=fan.base[10]/fan.T;sig=(initial[:,3]-s0)/unit
    assert np.max(sig)<1e-10
    entropy=np.linspace(float(min(sig))*1.02,0,3)
    x=np.r_[np.arange(-18,0.00001,.125),.01,.02]
    known=np.load(prior.OUT/'fine.npz');cache={float(xx):(tt,rr) for xx,tt,rr in zip(known['log_density_ratio'],known['T'],known['raw'])}
    rows=np.zeros((3,len(x),len(fan.base)));temps=np.zeros((3,len(x)));done=np.zeros((3,len(x)),bool)
    for j,sigma in enumerate(entropy):
        for i in range(len(x)-1,-1,-1):
            xx=float(x[i])
            if j==2 and xx in cache:tt,row=cache[xx]
            else:
                guess=np.interp(xx,known['log_density_ratio'][::-1],np.log(known['T'][::-1]))+sigma
                if xx>0:guess=np.log(fan.T)+.6*xx+sigma
                target=s0+sigma*unit
                for _ in range(8):
                    row=fan.call(np.log(fan.rho)+xx,guess);error=(row[3]-target)*np.exp(guess)/row[10]
                    if abs(error)<2e-12:break
                    assert abs(error)<.15;guess-=error
                else:raise AssertionError(('Native entropy root',j,i))
                tt=np.exp(guess)
            rows[j,i]=row;temps[j,i]=tt;done[j,i]=True
            np.savez_compressed(OUT/'eos-progress.npz',x=x,sigma=entropy,raw=rows,T=temps,done=done,native_calls=fan.calls)
    np.savez_compressed(OUT/'eos.npz',x=x,sigma=entropy,raw=rows,T=temps,rho0=fan.rho,cx=fan.cx,s0=s0,sunit=unit,
        initial_r=sample_r,initial_raw=initial,initial_T=T,initial_sigma=sig)
    lookup=EOS();controls=[]
    for xx in [-.037,-.7,-3.3,-10.3]:
        for fraction in [.25,.75]:
            sigma=entropy[0]*fraction;target=s0+unit*sigma
            guess=float(np.interp(xx,x,np.log(temps[-1])))+sigma
            for _ in range(8):
                row=fan.call(np.log(fan.rho)+xx,guess);error=(row[3]-target)*np.exp(guess)/row[10]
                if abs(error)<2e-12:break
                guess-=error
            else:raise AssertionError('EOS control root')
            p,u,gm,Tc,kap=lookup(np.array([np.exp(xx)]),np.array([sigma]))
            errors=[float(abs(p[0]*fan.rho*C*C/row[1]-1)),float(abs(u[0]*C*C/row[2]-1)),float(abs(gm[0]/row[4]-1)),float(abs(Tc[0]/np.exp(guess)-1))]
            controls.append(dict(log_density=xx,entropy_parameter=sigma,relative=errors))
    passed=max(max(r['relative']) for r in controls)<.002
    write(OUT/'eos.json',dict(classification='Counterexample candidate',native_calls=fan.calls,seconds=time.monotonic()-start,
        states=rows.shape[:2],reused_states=145,entropy_parameter_interval=entropy.tolist(),minimum_temperature_K=float(temps.min()),
        maximum_initial_pressure_relative=float(max(abs(initial[:,1]/np.interp(sample_r,env['r'],env['Pgas'])-1))),
        density_interval=x[[0,-1]].tolist(),controls=controls,passed=passed))
    assert passed,'Native EOS interpolation control failed'
    signal.alarm(0);print((OUT/'eos.json').read_text(),flush=True)


class EOS:
    def __init__(self):
        self.d=np.load(OUT/'eos.npz');self.rho0=float(self.d['rho0']);self.cx=float(self.d['cx']);self.sigma=self.d['sigma']
        raw=self.d['raw'];self.x=self.d['x'];self.floor=np.exp(self.x[0]);self.top=np.exp(self.x[-1])
        # rmue is a density, not inverse mu_e. The atomic mass constant
        # converts it directly to the native free-electron number density.
        values=np.stack([np.log(raw[:,:,1]/(self.rho0*C*C)),raw[:,:,2]/C**2,raw[:,:,4],np.log(self.d['T']),raw[:,:,13]/raw[:,:,0]],axis=-1)
        self.fun=PchipInterpolator(self.x,values,axis=1,extrapolate=False)

    def __call__(self,rho,sigma):
        assert np.max(rho)<self.top*(1+1e-9),('EOS high-density domain',np.max(rho),self.top)
        active=rho>=self.floor
        assert not np.any(active & ((sigma<self.sigma[0]-1e-8)|(sigma>1e-8))),('EOS entropy domain',float(np.min(sigma[active])),float(np.max(sigma[active])))
        xx=np.log(np.maximum(rho,self.floor));data=self.fun(xx)
        z=np.clip(sigma,self.sigma[0],0);frac=(z-self.sigma[0])/(self.sigma[1]-self.sigma[0]);i=np.minimum(frac.astype(int),1);f=frac-i
        ids=np.arange(len(rho));v=(1-f[:,None])*data[i,ids]+f[:,None]*data[i+1,ids]
        p=np.exp(v[:,0])*active;u=v[:,1]*active;gamma=v[:,2];T=np.exp(v[:,3]);kap=6.6524587321e-25/1.66053906660e-24*v[:,4]*active
        return p,u,gamma,T,kap


def minmod(a,b):return np.where(a*b>0,np.sign(a)*np.minimum(abs(a),abs(b)),0.)


class Flow:
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
        self.hist_v=np.array([np.interp(ri,previous['radius'],v) for v in previous['velocity']])*100/C
        self.F0=float(env['Linfinity'])/(4*np.pi*self.RJ**2*self.a0**2*C)/(self.eos.rho0*C*C)
        self.max_scattering_work=0.;self.max_optical=0.

    def conserved(self,rho,v,sigma):
        p,u,gamma,T,kap=self.eos(rho,sigma);W=1/np.sqrt(1-v*v);D=rho*W;h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor)
        S=rho*h*W*W*v;U=np.array([D,S,D*sigma]);flux=np.array([D*v,S*v+p,D*sigma*v])
        cs=np.sqrt(gamma*p/np.maximum(rho*h,1e-100));return U,flux,(p,u,gamma,T,kap,cs)

    def primitive(self,U):
        D=np.maximum(U[0],0.);sigma=np.divide(U[2],D,out=np.zeros_like(D),where=D>=self.eos.floor);v=np.zeros_like(D)
        for _ in range(3):
            rho=D*np.sqrt(1-v*v);p,u,_,_,_=self.eos(rho,sigma);h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor)
            q=np.divide(U[1],D*h,out=np.zeros_like(D),where=D>=self.eos.floor);v=q/np.sqrt(1+q*q)
        return np.array([D*np.sqrt(1-v*v),v,sigma])

    def rhs(self,U,t):
        V=self.primitive(U);rho,v,sigma=V;_,_,thermo=self.conserved(rho,v,sigma);p,u,gamma,T,kap,cs=thermo
        ghost=np.array(self.left);ghost[1]=np.interp(t,self.hist_t,self.hist_v) if not self.flat else 0.
        ext=np.column_stack([ghost,V,np.zeros(3)]);slope=np.zeros_like(ext);slope[:,1:-1]=minmod(ext[:,1:-1]-ext[:,:-2],ext[:,2:]-ext[:,1:-1])
        L=ext[:,:-1]+slope[:,:-1]/2;R=ext[:,1:]-slope[:,1:]/2
        UL,FL,tl=self.conserved(*L);UR,FR,tr=self.conserved(*R)
        sl=np.minimum(0,np.minimum((L[1]-tl[-1])/(1-L[1]*tl[-1]),(R[1]-tr[-1])/(1-R[1]*tr[-1])))
        sr=np.maximum(0,np.maximum((L[1]+tl[-1])/(1+L[1]*tl[-1]),(R[1]+tr[-1])/(1+R[1]*tr[-1])))
        flux=np.divide(sr*FL-sl*FR+sl*sr*(UR-UL),sr-sl,out=np.zeros_like(FL),where=sr>sl)
        flux*=self.af*self.area;rate=-C*np.diff(flux)/self.vol
        W=1/np.sqrt(1-v*v);h=self.eos.cx+u+p/np.maximum(rho,self.eos.floor);E=rho*h*W*W-p
        geom=-(self.r/self.RJ)**2*E*self.ap
        if not self.flat:geom+=2*self.a*self.r/self.RJ**2*p
        rate[1]+=C*self.dx/self.vol*geom
        # Elastic electron scattering: in the gas frame G0=0, Gr=rho*kappa*F0/c.
        # Current mechanical work is removed from the independent photon flux.
        area_c=(self.r/self.RJ)**2
        localF=self.F0/(area_c*(self.a/self.a0)**2);mu=np.sqrt(np.maximum(0,1-(self.RJ/self.r)**2*(self.a/self.a0)**2))
        Eg=2*localF/(1+mu);Pg=2*localF*(1+mu+mu*mu)/(3*(1+mu))
        work=np.zeros_like(v)
        for _ in range(2):
            F=localF+np.r_[0,np.cumsum(work[:-1])]/(self.a*self.a*area_c*C)
            fcom=self.eos.rho0*rho*kap*W*W*((1+v*v)*F-v*(Eg+Pg))
            momentum=self.a*C*W*fcom
            work=-self.vol*self.a*self.a*C*W*v*fcom
        if self.flat:momentum[:]=0;work[:]=0
        rate[1]+=momentum
        self.max_scattering_work=max(self.max_scattering_work,float(abs(work.sum())))
        self.max_optical=max(self.max_optical,float(np.sum(self.vol/area_c*self.eos.rho0*rho*kap)))
        # Independent isentropic Killing-energy flux for the conservation audit.
        def energy(vv,tt):
            rr,vv,ss=vv;pp,uu=tt[:2];ww=1/np.sqrt(1-vv*vv);wm=vv*vv/(np.sqrt(1-vv*vv)*(1+np.sqrt(1-vv*vv)))
            DD=rr*ww
            ee=self.af*(rr*self.eos.cx*ww*wm+(rr*uu+pp)*ww*ww-pp)+(self.af-self.a0)*self.eos.cx*DD
            ff=self.af*((rr*self.eos.cx*ww*wm+(rr*uu+pp)*ww*ww)*vv)+(self.af-self.a0)*self.eos.cx*DD*vv
            return ee,ff
        eL,fL=energy(L,tl);eR,fR=energy(R,tr)
        ef=np.divide(sr*fL-sl*fR+sl*sr*(eR-eL),sr-sl,out=np.zeros_like(sr),where=sr>sl)*self.af*self.area
        ledgers=np.array([C*(flux[0,0]-flux[0,-1]),C*(ef[0]-ef[-1])-work.sum(),-work.sum()])
        dt=.35*self.dx/np.max(C*self.a/self.B*(abs(v)+cs+1e-100))
        return rate,ledgers,dt

    def energy(self,U):
        rho,v,sigma=self.primitive(U);p,u,_,_,_=self.eos(rho,sigma);W=1/np.sqrt(1-v*v);wm=v*v/(np.sqrt(1-v*v)*(1+np.sqrt(1-v*v)))
        return self.a*(rho*self.eos.cx*W*wm+(rho*u+p)*W*W-p)+(self.a-self.a0)*self.eos.cx*U[0]

    def run(self,label):
        assert not (OUT/f'{label}.json').exists();start=time.monotonic();end=float(self.hist_t[-1]);t=0.;U=self.U.copy();ledger=np.zeros(3);history=[];steps=0;discard=np.zeros(3)
        e0=self.energy(U);snapshots=[];next_dump=0.
        while t<end:
            k,l,dt=self.rhs(U,t);dt=min(dt,end-t);trial=U+dt*k
            assert np.min(trial[0])>-1e-13,('Positive density stage',np.min(trial[0]))
            k2,l2,_=self.rhs(trial,t+dt);nxt=(U+trial+dt*k2)/2
            tiny=nxt[0]<self.eos.floor;discard+=np.sum(nxt[:,tiny]*self.vol[tiny],axis=1);nxt[:,tiny]=0
            ledger+=dt*(l+l2)/2;U=nxt;t+=dt;steps+=1
            assert steps<20000,'Step budget'
            if t>=next_dump or t>=end:
                rho,v,sigma=self.primitive(U);outside=float(np.sum(U[0,self.x>=0]*self.vol[self.x>=0]))
                history.append([t,outside,float(np.sum((self.energy(U)-e0)*self.vol)),*ledger])
                snapshots.append(U.copy());next_dump+=end/32
        rho,v,sigma=self.primitive(U);p,u,_,T,_=self.eos(rho,sigma);response=max(abs(history[-1][2]),float(np.sum(self.vol*(rho*v*v+p))),1e-100)
        baryon=abs(np.sum((U[0]-self.initial[0])*self.vol)-ledger[0]+discard[0])/np.sum(self.initial[0]*self.vol)
        energy_error=abs(history[-1][2]-ledger[1])/response
        scale=4*np.pi*self.RJ**2*self.eos.rho0
        trace=self.eos.cx*(rho-U[0])+rho*u-3*p
        initialrho,iv,iss=self.primitive(self.initial);ip,iu,_,_,_=self.eos(initialrho,iss)
        trace-=initialrho*iu-3*ip
        trace_moment=float(np.sum(self.vol*trace))
        np.savez_compressed(OUT/f'{label}.npz',U=U,initial=self.initial,x_cm=self.x,volume=self.vol,rho=rho*self.eos.rho0,velocity_cm_s=v*C,T=T,pressure=p*self.eos.rho0*C*C,
            history=history,snapshots=snapshots,conserved_trace_moment=trace_moment,baryon_discard=discard)
        row=dict(classification='Counterexample candidate',cells=len(self.x),steps=steps,seconds=time.monotonic()-start,
            baryon_ledger_relative=float(baryon),isentropic_Killing_energy_response_relative=float(energy_error),
            gas_outside_original_radius_g=history[-1][1]*scale,discarded_dilute_mass_g=float(discard[0]*scale),
            maximum_velocity_m_s=float(max(v)*C/100),integrated_trace_energy_erg=trace_moment*scale*C*C,
            matter_received_scattering_work_erg=float(ledger[2]*scale*C*C),maximum_electron_scattering_optical_depth=self.max_optical,
            horizon_seconds=end,full_GR_scalar_feedback=False,absorption_and_chemistry_kinetics=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/f'{label}.json',row);print(label,json.dumps(row),flush=True);return row


def run():
    assert not (OUT/'result.json').exists();start=time.monotonic();signal.alarm(240)
    spec=json.loads((OUT/'plan.json').read_text())
    for p,h in spec['bindings'].items():assert old.photons.digest(old.ROOT/p)==h,p
    assert json.loads((OUT/'eos.json').read_text())['passed']
    pilot=Flow(224).run('pilot-224');forecast=1.3*pilot['seconds']*(2+16+64)
    write(OUT/'pilot-budget.json',dict(classification='Counterexample candidate',measured_pilot_seconds=pilot['seconds'],forecast_seconds=forecast,
        assumption='Explicit CFL cost proportional to cells squared for the same fixed horizon;30percent margin. Finer-grid wave speeds and costs are not yet measured.'))
    assert forecast<235,'Production forecast over budget'
    flat=Flow(224,flat=True).run('flat-224');reference=json.loads((prior.OUT/'result.json').read_text())['gas_outside_original_cut_g']
    flat_error=abs(flat['gas_outside_original_radius_g']/reference-1)
    write(OUT/'flat-control.json',dict(classification='Counterexample candidate',native_exact_fan_mass_relative=flat_error,passed=flat_error<.03))
    assert flat_error<.03,'Independent native flat fan control failed'
    rows=[Flow(n).run(f'cells-{n}') for n in [896,1792]]
    errors={k:abs(rows[0][k]/rows[1][k]-1) for k in ['gas_outside_original_radius_g','integrated_trace_energy_erg']}
    passed=all(r['baryon_ledger_relative']<1e-10 and r['isentropic_Killing_energy_response_relative']<.02 for r in rows) and max(errors.values())<.02
    write(OUT/'result.json',dict(classification='Counterexample candidate',passed=bool(passed),refinement=errors,seconds=time.monotonic()-start,
        gas_replaced_not_added=True,actual_radial_nonlinear_gas_evolved=True,current_elastic_scattering_work_paired=True,
        complete_photon_transport=False,full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False))
    signal.alarm(0)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','bank','run']);globals()[parser.parse_args().action]()
