"""Counterexample candidate: whole-star momentary scalar equilibrium.

Coordinate gas baryons/entropy/ions and photon canonical momenta are fixed
in the metric substep. Free mechanical and thermal equilibrium are NOT imposed.
The two core photon debits are sensitivity examples, not evolved completions
or an uncertainty enclosure. Reuse the actual current source; no fluid replay.
"""
from pathlib import Path
import hashlib, json, resource, signal, sys, time
import numpy as np
import mpmath as mp
import sympy as sp
import def_retained_native_acoustic as native
import def_native_anisotropic_gr as gr
import gr_scalar_nonlinear_exterior as exterior

OUT=Path('retained-static-response154-work')
BASE=Path('retained-metric-return152-work')
ACOUSTIC=native.OUT/'response-centered'
LD=np.longdouble; C=gr.C; G=gr.G
read=native.read; write=native.write; sha=native.sha
KEYS=['baryon_g','gas_nonrest_energy_erg','photon_energy_erg',
      'nonrest_trace_erg','metric_stress_erg','inner_cumulative_energy_erg',
      'outer_cumulative_energy_erg']


def symbolic():
    r,a,ap,f,fp,fpp,V,S=sp.symbols('r a ap f fp fpp V S',nonzero=True)
    wave=a*a*(r*fpp+2*fp)+a*ap*(r*fp+f)-V*r*f+S
    flux=a*r*r*fpp+(ap*r*r+2*a*r)*fp-(V-a*ap/r)*r*r*f/a+r*S/a
    assert sp.simplify(wave*r/a-flux)==0
    x,alpha,k,N,A,b=sp.symbols('x alpha k N A b',nonzero=True)
    dDphi=-12*sp.pi*r*N*N*A**4*(x*alpha+3*alpha**2)*k
    dDl=-4*sp.pi*r*N*N*A**4*(x+3*alpha)*k
    assert sp.simplify(-(dDphi+dDl*x)/r-4*sp.pi*N*N*A**4*k*(x+3*alpha)**2)==0
    return dict(classification='Proven',passed=True,
        flux='B=a*r^2*f_prime; B_prime=(V-Vgeom)*r^2*f/a-r*S/a; Vgeom=a*a_prime/r; a=N*sqrt(b).',
        scope='Algebraic rewriting of the existing anisotropic first variation only; no new equilibrium or physical closure theorem.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(gr.__file__),Path(exterior.__file__),
           gr.flow.INPUT/'balanced-20.npz',BASE/'fields/source-128.npz',
           ACOUSTIC/'infinity/result.json',ACOUSTIC/'infinity/charge-128-a8-r8.npz']
    paths += [folder/f'source-{n}-reference-128.npz' for folder in [BASE/'gr',ACOUSTIC/'gr-precision'] for n in [64,128]]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='13cadadfdac7802cc47f1b2083c1047921f380fe',
        claim='Construct the missing full-radius, regular-center, exact-static-vacuum scalar/metric operator on the CURRENT momentarily balanced background. Apply current saved source endpoints and expose dependence on the unrepresented core energy-debit profile.',
        comparator='Momentary scalar and polar metric balance at fixed coordinate baryons, entropy, ionic inventory and photon canonical momenta, with declared additional frozen source increments. No free radial mechanical equilibrium or thermal stationarity is imposed. This is not the complete static EFT comparator.',
        core='Restore the saved core Gamma profile and isotropic fourth photon moment below the19+512 causal shell. The previous shell coefficient lookup is not licensed as a core EOS. Deep explicit source is unknown; compare two predeclared isotropic photon depletions proportional to the saved photon energy in0.25-0.5 and0.5-0.75 times the inner represented radius. Each has exactly the same supplied Killing-energy debit and zero baryon/ionic change. These are static source sensitivity examples, NOT alternative validated3.434ms trajectories, matched photon inventories or rigorous extrema.',
        current_source='Phase152 fields/source-128 contains the full retained background, actual EOS correction, Phase150 and direct Phase151 sources. Add Phase152 GR return and Phase153 native-acoustic return once, interpolating their existing17 source knots. Evaluate time0,midpoint,endpoint only.',
        exterior='Use the differentiated exact Just vacuum map beyond the initial Cauchy edge. Escaped photons outside that edge are intentionally NOT replaced by vacuum in a claim about the whole physical star: their static source is unclosed here, so these are compact-source components only.',
        decision='Determine the whole-star momentary susceptibility and whether the available shell history/inner port uniquely fixes a static readout. Do not label a mismatch with the autonomous retarded transient as external-drive nonabsorption.',
        reuse='No new fluid steps, EOS evaluations, stellar mesh, time horizon, angular/frequency bank or orbital fit. Existing full-star source panels at orders4/8; quadrature is not physical mesh refinement.',
        budget=dict(pilot_seconds=45,production_seconds=120,audit_seconds=45,CPU_threads=1,virtual_GiB=3,max_picard=12),
        gates=dict(quadrature=.002,response_clock=.02,picard=1e-12,contraction=.05,energy_inventory=1e-12,flat_control=1e-10,positive_remaining_radiation=True),
        stop='Preserve failed controls and timings; do not automatically add source profiles, quadrature orders, clocks or relaxed gates. Core sensitivity examples cannot certify global error or exclusion.',
        bindings={str(p):sha(p) for p in paths}))
    write(OUT/'symbolic.json',symbolic())


def source(n):
    d={k:v.copy() for k,v in np.load(BASE/'fields/source-128.npz').items()}
    for folder in [BASE/'gr',ACOUSTIC/'gr-precision']:
        inc=np.load(folder/f'source-{n}-reference-128.npz')
        for k in KEYS:
            value=inc[k];i=np.clip(np.searchsorted(inc['t'],d['t'],side='right')-1,0,len(inc['t'])-2)
            w=(d['t']-inc['t'][i])/(inc['t'][i+1]-inc['t'][i]);w=w.reshape((-1,)+(1,)*(value.ndim-1))
            d[k]=d[k]+(1-w)*value[i]+w*value[i+1]
    return d


class Whole(gr.Response):
    def __init__(self):
        super().__init__()
        # Reuse the actual imported core table, not a shell gamma extrapolation.
        self.core=gr.flow.initial.Data().bg.d

    def coeff(self,r):
        z=super().coeff(r);inside=np.asarray(r)<self.faces[0]
        gamma=np.interp(np.asarray(r)[inside],self.core['radius_cm'],self.core['gamma1'])
        dK=gamma*z['Pg'][inside]-z['Kg'][inside]
        dR=z['Er'][inside]/5-z['R4'][inside]
        x=np.asarray(r)[inside]*z['Phi'][inside];a=z['alpha'][inside]
        k=4*np.pi*z['lapse'][inside]**2*z['A4'][inside]
        z['V'][inside]+=k*(dK*(x+3*a)**2-dR*x*x)
        z['K'][inside]+=k*(-dK*(x+3*a)+dR*x)/z['b'][inside]
        z['Kg'][inside]+=dK;z['R4'][inside]+=dR
        return z


class Static:
    def __init__(self,model,order):
        self.model=model;self.q=q=gr.flow.initial.Quadrature(model.bg.edges,order)
        self.r=r=q.r.astype(LD);self.z=z={k:np.asarray(v,LD).reshape(r.shape) for k,v in model.coeff(q.r.ravel()).items()}
        self.a=a=z['lapse']*np.sqrt(z['b']);self.R=LD(model.bg.rend)
        self.geom=z['lapse']**2*(2*z['mass']/r**3+4*LD(np.pi)*z['A4']*(z['Pg']+z['Pr']-z['Eg']-z['Er']))
        self.W=(z['V']-self.geom)*r*r/a
        self.physical_volume=4*LD(np.pi)*r*r*np.exp(-6*z['phi']**2)/np.sqrt(z['b'])
        self.outer_start=len(q.h)-len(model.gamma)
        self.volumes=q.h[self.outer_start:]*(self.physical_volume[self.outer_start:]@q.w)
        bg=model.bg.z;mu=LD(bg['mass_faces'][-1])/self.R
        self.qs=LD(bg['panel_Phi'][-1,-1]) # replaced by exact face flux below
        self.asurf=np.exp(LD(bg['nu_faces'][-1]))*np.sqrt(1-2*mu)
        self.qs=LD(bg['scalar_flux_faces'][-1])/(self.asurf*self.R)
        self.bsurf=1-2*mu
        mp.mp.dps=65;u,v=mp.mpf(str(mu)),mp.mpf(str(self.qs))
        fun=[lambda u,v:u+v*v*exterior.exact(u,v)[0],
             lambda u,v:v*exterior.exact(u,v)[1],lambda u,v:v*exterior.exact(u,v)[2]]
        self.der=np.array([[LD(str(mp.diff(lambda x:fn(x,v),u))),LD(str(mp.diff(lambda x:fn(u,x),v)))] for fn in fun])
        self.M=LD(bg['ADM_mass']);self.K=LD(bg['scalar_K'])
        self.den=1+self.der[2,0]*self.bsurf*self.qs
        absB,absfaces=q.integrate(abs(self.W));_,norm=q.integrate(absB/(a*r*r))
        self.eta=float(abs(self.der[2,1]/(self.asurf*self.R*self.den))*absfaces[-1]+norm[-1])
        assert self.eta<.05,self.eta

    def solve(self,forcing,J=LD(0),incident=LD(0)):
        q=self.q;r=self.r;f=np.zeros_like(r);history=[]
        for _ in range(12):
            B,bfaces=q.integrate(self.W*f+forcing)
            integ,faces=q.integrate(B/(self.a*r*r))
            fs=(incident-self.der[2,1]*bfaces[-1]/(self.asurf*self.R)-self.der[2,0]*J/self.R)/self.den
            new=fs-(faces[-1]-integ);scale=max(np.max(abs(new)),LD('1e-290'))
            change=np.max(abs(new-f));err=float(change/scale);history.append(err);f=new
            if err<1e-12:break
        else:raise AssertionError(('static Picard',history))
        du=self.bsurf*self.qs*fs+J/self.R;dv=bfaces[-1]/(self.asurf*self.R)
        dM,dK,far=self.R*(self.der@np.array([du,dv]));far=fs+far/self.R
        return dict(phi=f,flux=B,history=history,eta=self.eta,
            scalar_charge=float(-dK/self.M),mass_term=float(self.K*dM/self.M**2),
            normalized_charge=float(-dK/self.M+self.K*dM/self.M**2),
            delta_ADM_cm=float(dM),delta_K_cm=float(dK),boundary_error=float(abs(far-incident)/scale),
            finite_operator_remainder=float(self.eta/(1-self.eta)*change),maximum_phi=float(scale))

    def forcing(self,d,it,profile):
        z=self.z;r=self.r;q=self.q;s=self.outer_start;factor=LD(G)/LD(C)**4
        e=np.zeros_like(r);trace=e.copy();stress=e.copy()
        rest=np.asarray(d['baryon_g'][it],LD)*LD(d['cx'])*LD(C)**2
        for arr,values in [(e,rest+d['gas_nonrest_energy_erg'][it]+d['photon_energy_erg'][it]),
                           (trace,rest+d['nonrest_trace_erg'][it]),(stress,d['metric_stress_erg'][it])]:
            arr[s:]=np.asarray(values,LD)[:,None]/self.volumes[:,None]*factor
        debit=LD(0);fraction=LD(0)
        if profile!='shell':
            lo,hi={'inner':(.25,.5),'outer':(.5,.75)}[profile]
            x=r/LD(self.model.faces[0]);shape=np.where((x>lo)&(x<hi),np.sin(np.pi*(x-lo)/(hi-lo))**2,0)
            weight=4*LD(np.pi)*r*r*z['A4']*z['lapse']/np.sqrt(z['b'])
            _,capacity=q.integrate(weight*z['Er']*shape)
            target=LD(d['inner_cumulative_energy_erg'][it])*factor
            fraction=target/capacity[-1];assert 0<=fraction<1,(profile,float(fraction))
            de=-fraction*z['Er']*shape;e+=de;stress+=2*de/3
            _,check=q.integrate(weight*de);debit=check[-1]/factor
        density=4*LD(np.pi)*r*r*z['A4']*z['lapse']/np.sqrt(z['b'])*e
        integral,faces=q.integrate(density);J=np.sqrt(z['b'])/z['lapse']*integral
        force=-r*z['K']*J/self.a+4*LD(np.pi)*r*r*z['lapse']**2*z['A4']/self.a*(z['alpha']*trace+r*z['Phi']*stress)
        Jout=np.sqrt(self.bsurf)/np.exp(LD(self.model.bg.z['nu_faces'][-1]))*faces[-1]
        return force,Jout,dict(core_debit_erg=float(debit),maximum_radiation_depletion=float(fraction),
            target_debit_erg=float(d['inner_cumulative_energy_erg'][it]),core_baryon_change=0,
            omitted_escaped_energy_erg=float(d['outer_cumulative_energy_erg'][it]))


def construct():
    native.prior.initialize()
    return Whole()


def pilot():
    start=time.monotonic();model=construct();setup=time.monotonic()-start
    tick=time.monotonic();s=Static(model,4);d=source(128);f,J,meta=s.forcing(d,-1,'shell');v=s.solve(f,J)
    cost=time.monotonic()-tick
    np.savez_compressed(OUT/'pilot.npz',r=s.r,phi=v.pop('phi'),flux=v.pop('flux'))
    write(OUT/'pilot.json',dict(classification='Counterexample candidate',setup_seconds=setup,
        solve_seconds=cost,forecast_seconds=setup+12*cost*2+10,eligible=setup+12*cost*2+10<120,
        result=v,metadata=meta,seconds=time.monotonic()-start))


def run():
    assert read(OUT/'pilot.json')['eligible'];start=time.monotonic();model=construct();rows=[]
    for order in [4,8]:
        s=Static(model,order)
        zero=s.solve(np.zeros_like(s.r));assert zero['normalized_charge']==0
        susceptibility=s.solve(np.zeros_like(s.r),incident=LD(1))
        np.savez_compressed(OUT/f'susceptibility-g{order}.npz',r=s.r,phi=susceptibility.pop('phi'),flux=susceptibility.pop('flux'))
        rows.append(dict(case='unit_static_incident',order=order,**susceptibility))
        for n in [64,128]:
            d=source(n)
            for it in [0,16,32]:
                for profile in ['shell','inner','outer']:
                    force,J,meta=s.forcing(d,it,profile);v=s.solve(force,J)
                    np.savez_compressed(OUT/f'static-{n}-g{order}-t{it}-{profile}.npz',r=s.r,phi=v.pop('phi'),flux=v.pop('flux'),forcing=force,Jout=J)
                    rows.append(dict(case=profile,order=order,steps=n,index=it,time=float(d['t'][it]),**v,**meta))
    def pick(n,o,i,c):return next(r for r in rows if r.get('steps')==n and r['order']==o and r.get('index')==i and r['case']==c)
    comparisons=[]
    for i in [0,16,32]:
        for c in ['shell','inner','outer']:
            a,b,h=pick(128,8,i,c),pick(128,4,i,c),pick(64,8,i,c)
            scale=max(abs(a['scalar_charge']),abs(a['mass_term']),1e-290)
            comparisons.append(dict(index=i,case=c,quadrature=abs(a['normalized_charge']-b['normalized_charge'])/scale,
                response_clock=abs(a['normalized_charge']-h['normalized_charge'])/scale))
    inv=max(abs(r['core_debit_erg']/r['target_debit_erg']+1) for r in rows if r.get('case') in ['inner','outer'] and r['target_debit_erg'])
    gates=dict(quadrature=max(r['quadrature'] for r in comparisons),response_clock=max(r['response_clock'] for r in comparisons),
               energy_inventory=inv,boundary=max(r['boundary_error'] for r in rows))
    result=dict(classification='Counterexample candidate',passed=gates['quadrature']<.002 and gates['response_clock']<.02 and inv<1e-12,
        controls=gates,comparisons=comparisons,rows=rows,seconds=time.monotonic()-start,
        current_background_global_scalar_metric_operator=True,regular_center_and_exact_static_vacuum=True,
        complete_thermal_static_comparator=False,free_mechanical_static_comparator=False,matched_core_evolution=False,
        escaped_photons_in_static_comparator=False,static_EFT_nonabsorption=False,external_drive_identified=False,
        full_physical_error_enclosed=False,full_goal_complete=False)
    write(OUT/'result.json',result);assert result['passed'],gates


def audit():
    result=read(OUT/'result.json');assert result['passed'];assert symbolic()['passed']
    for path,h in read(OUT/'plan.json')['bindings'].items():assert sha(path)==h,path
    # Independent analytic flat-space uniform source: (r^2 f')'=r^2.
    q=gr.flow.initial.Quadrature(np.linspace(0,1,41),8);r=q.r.astype(LD)
    B,Bf=q.integrate(r*r);integ,end=q.integrate(B/(r*r))
    f=-Bf[-1]-(end[-1]-integ);exact=(r*r-3)/6
    flat=float(np.max(abs(f-exact)));assert flat<1e-10
    # Recompute every charge from independently saved exact-map mass/charge outputs.
    z=np.load(gr.flow.INPUT/'balanced-20.npz');M=float(z['ADM_mass']);K=float(z['scalar_K'])
    identities=[]
    for row in result['rows']:
        q=-row['delta_K_cm']/M+K*row['delta_ADM_cm']/M**2
        identities.append(abs(q-row['normalized_charge'])/max(abs(q),1e-290))
    assert max(identities)<1e-12
    end={r['case']:r for r in result['rows'] if r.get('steps')==128 and r['order']==8 and r.get('index')==32}
    write(OUT/'audit.json',dict(classification='Counterexample candidate',passed=True,flat_poisson_error=flat,
        normalized_charge_identity=max(identities),symbolic=symbolic(),
        endpoint_core_profile_charge_difference=end['outer']['normalized_charge']-end['inner']['normalized_charge'],
        limits='These examples do not enclose core-source uncertainty and are not causal3.434ms solutions. No comparison to the whole retarded Bondi charge is certified.'))


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','pilot','run','audit'];start=time.monotonic();cpu=time.process_time();error=None
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    native.deadline({'prepare':30,'pilot':45,'run':120,'audit':45}[action])
    try:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():
            path=OUT/f'{action}-receipt.json';assert not path.exists()
            write(path,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
                peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
