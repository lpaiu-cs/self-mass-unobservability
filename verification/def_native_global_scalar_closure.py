"""Counterexample candidate: causal exterior return and null-infinity charge.

Use the actual projected-background photon history. No new matter/EOS run.
The global bound is for the declared linear scalar operator and source model;
it does not certify the unclosed constitutive or nonlinear source errors.
"""
from pathlib import Path
import json,signal,sys,time
import numpy as np
from numpy.polynomial import legendre as leg
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import sympy as sp
import def_native_characteristic_gr as compact
import verify_native_anisotropic_gr as centers

flow=compact.base.flow;C=flow.C;G=flow.G;LD=np.longdouble;write=flow.write;sha=flow.sha
OUT=compact.OUT.parent.parent/'def-native-global-scalar-closure'


class Exterior:
    def __init__(self):
        self.response=compact.Response();self.model=self.response.model;self.bg=self.response.bg;self.r0=self.bg.rend
        self.M=float(self.bg.z['ADM_mass']);self.K=self.bg.vk*self.r0
        self.N0=float(self.bg.metric(np.array([self.r0/self.model.m.R]))[1][0]);self.T=flow.old.END
        self.t=np.linspace(0,self.T,17);self.kernels={}

    def metric(self,q):
        mass,cc=self.bg.solution.sol(np.asarray(q));b=1-2*mass*q;N=cc/np.sqrt(b)
        return N,b,cc

    def delay(self,mu,flat=False):
        theta0=np.arccos(mu);impact=self.r0*np.sqrt(1-mu*mu)
        def rhs(s,value):
            th=theta0*(1-s);si=np.sin(th);co=np.cos(th);q=si/np.sqrt(1-mu*mu)
            N,b,cc=(np.ones_like(mu),)*3 if flat else self.metric(q)
            N0=1. if flat else self.N0;v=np.sqrt(np.maximum(1-si*si*(N/N0)**2,1e-300))
            return theta0*impact*co*(N/N0)**2/(C*cc*v*(1+v))
        sol=solve_ivp(rhs,[0,1],np.zeros(len(mu)),rtol=2e-12,atol=2e-15,dense_output=True)
        assert sol.success
        return sol,rhs

    def kernel(self,angular,radial):
        key=(angular,radial)
        if key in self.kernels:return self.kernels[key]
        start=time.monotonic();edges=self.model.bulk.edges_mu;positive=edges[edges>=0];gx,gw=leg.leggauss(angular)
        # Resolve arrival-at-infinity cuts inside the existing physical bins.
        # The last bin spans a narrow causal front; uncut4/8 rules miss it.
        if not hasattr(self,'angular_cuts'):
            lower=float(positive[-2]);limit=float(self.delay(np.array([lower]))[0].sol(1.)[0])
            assert self.T<limit
            roots=[brentq(lambda v:float(self.delay(np.array([v]))[0].sol(1.)[0]-t),lower,1-1e-14,xtol=5e-15) for t in self.t[1:]]
            self.angular_cuts=np.unique(np.r_[positive,roots])
        cuts=self.angular_cuts
        mu=(cuts[:-1,None]+np.diff(cuts)[:,None]*(gx+1)/2).ravel();mw=(np.diff(cuts)[:,None]*gw/2).ravel()
        bins=np.searchsorted(positive,mu,side='right')-1
        solution,_=self.delay(mu);infinity=solution.sol(1.);rx,rw=leg.leggauss(radial)
        nodes=[];owners=[];weights=[];inversion=[]
        for j in range(len(mu)):
            targets=self.t[(self.t>0)&(self.t<infinity[j])]
            roots=[brentq(lambda s:float(solution.sol(s)[j]-v),0.,1.,xtol=5e-15,rtol=1e-14) for v in targets]
            last=roots[-1] if len(targets) and targets[-1]==self.T else 1.
            cuts=np.unique(np.r_[0.,roots,last]);cuts=cuts[cuts<=last]
            ss=(cuts[:-1,None]+np.diff(cuts)[:,None]*(rx+1)/2).ravel()
            ww=(np.diff(cuts)[:,None]*rw/2).ravel()*np.arccos(mu[j])*mw[j]*mu[j]
            nodes.extend(ss);owners.extend([j]*len(ss));weights.extend(ww)
            inversion.extend([abs(float(solution.sol(s)[j])-v) for s,v in zip(roots,targets)])
        nodes=np.asarray(nodes);owners=np.asarray(owners);weights=np.asarray(weights)
        th=np.arccos(mu[owners])*(1-nodes);si=np.sin(th);co=np.cos(th);impact=self.r0*np.sqrt(1-mu[owners]**2)
        q=si/np.sqrt(1-mu[owners]**2);N,b,cc=self.metric(q);v=np.sqrt(1-si*si*(N/self.N0)**2)
        delay=solution.sol(nodes)[owners,np.arange(len(nodes))]
        mass=self.K*si*co/(cc*b*impact**2)
        stress=self.K*si*si*co*(N/self.N0)**2/(cc*cc*impact*v)
        row=dict(mu=mu,mw=mw,bins=bins,infinity=infinity,weights=weights,mass=mass,stress=stress,delay=delay,node_bins=bins[owners],
            inversion_error=max(inversion,default=0.),seconds=time.monotonic()-start)
        self.kernels[key]=row;return row

    def luminosity(self,reference):
        d=np.load(flow.OUT/f'coupled-{reference}.npz');ids=[np.argmin(abs(d['snapshot_t']-t)) for t in self.t]
        assert np.max(abs(d['snapshot_t'][ids]-self.t))<1e-18
        I=d['snapshot_I'][ids].sum(1)[:,-1,self.model.mu>0,:]
        return 2*np.pi*C*self.model.area[-1]*(I@(self.model.number*self.model.E))

    def evaluate(self,reference,angular=8,radial=8):
        start=time.monotonic();k=self.kernel(angular,radial);lum=self.luminosity(reference);assert np.min(lum)>=0
        poly=flow.green.polynomial(self.t,lum);H=poly.antiderivative();HH=H.antiderivative();mass=[];stress=[];arrived=[]
        for t in self.t:
            at=np.maximum(t-k['delay'],0.);h=H(at)[np.arange(len(at)),k['node_bins']];hh=HH(at)[np.arange(len(at)),k['node_bins']]
            mass.append(float(-G/C**3*np.sum(k['weights']*k['mass']*hh,dtype=LD)))
            stress.append(float(-G/(2*C**4)*np.sum(k['weights']*k['stress']*h,dtype=LD)))
            arrived.append(float(np.sum(k['mw']*k['mu']*H(np.maximum(t-k['infinity'],0.))[np.arange(len(k['mu'])),k['bins']],dtype=LD)))
        positive=self.model.bulk.edges_mu[self.model.bulk.edges_mu>=0]
        emitted=float(H(self.t[-1])@(np.diff(positive**2)/2))
        d=np.load(compact.base.OUT/f'source-{reference}.npz');port=abs(emitted/d['outer_cumulative_energy_erg'][-1]-1)
        data=dict(t=self.t,mass_U=np.array(mass),stress_U=np.array(stress),normalized_exterior=-(np.array(mass)+stress)/self.M,
            arrived_energy_erg=np.array(arrived),emitted_energy_erg=emitted)
        row=dict(classification='Counterexample candidate',reference=reference,angular=angular,radial=radial,seconds=time.monotonic()-start,
            kernel_seconds=k['seconds'],delay_inverse_error_seconds=k['inversion_error'],emitted_port_relative=float(port),
            endpoint_mass_charge=-mass[-1]/self.M,endpoint_stress_charge=-stress[-1]/self.M,
            endpoint_exterior_charge=float(data['normalized_exterior'][-1]),endpoint_arrived_energy_erg=arrived[-1])
        return data,row


def symbolic():
    r,N,b,Phi,e,mu=sp.symbols('r N b Phi e mu',positive=True)
    # Einstein-frame thin shell of Killing energy e=G*Einf/c^4.
    shell=e*sp.sqrt(b)/(4*sp.pi*N*r*r)
    stress=sp.simplify(-4*sp.pi*r*r*N*N*Phi*(1-mu*mu)*shell/(N*sp.sqrt(b)))
    assert sp.simplify(stress+e*Phi*(1-mu*mu))==0
    mass=sp.simplify(2*N*N*Phi/(r*b)*(-e*sp.sqrt(b)/N)/(N*sp.sqrt(b)))
    assert sp.simplify(mass+2*e*Phi/(r*b))==0
    return dict(classification='Proven',passed=True,
        shell='A packet with Killing energy e has integral S_stress dx=-e*Phi*(1-mu^2). Its compensating body debit gives S_mass dx=-2*e*Phi/(r*b)dr behind that packet.',
        delay='delta(r,mu0)=integral(r0,r)[1/mu-1]dr/(c*N*sqrt(b)) is the delay to null infinity relative to radial propagation. It stays finite; do not subtract two divergent flight times.',
        null_infinity='U_mass(u)=-c*G/c^4*integral dmu dr Phi/(r*b)*H2(u-delta,mu); U_stress(u)=-G/(2*c^4)*integral dmu dr Phi*(1-mu^2)/(N*sqrt(b)*mu)*H1(u-delta,mu). H1/H2 are the first/second primitives of actual outgoing Killing luminosity.',
        cone='For zero initial increments and no incoming wave, sources with optical x<x_outer-c*T cannot affect null-infinity retarded time corresponding to outer time T. For positive initial m and lapse<=1, inner travel time is at least (r_outer-r_inner)/c.',
        limits='Declared first variation and paired conserved exterior photons. Source embedding and nonlinear/constitutive errors are not removed by this algebra.')


def prepare():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(30)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='71ed2659e',
        claim='Compute actual angular exterior photon mass-debit/stress scalar forcing at null infinity on the inventory-matched background; close its declared global potential and deeper causal domain before combining the charge.',
        decision='Can exterior and unrepresented deep scalar forcing erase the existing compact residual? Preserve source/constitutive/nonlinear limits and do not inherit the old-background positive interval.',
        method='Use the same four outward angular bins and their constant sub-bin occupation, actual64/128 emission histories and exact scalar vacuum. Integrate the finite nonradial-minus-radial delay in a compactified angle coordinate all the way to infinity. Split radial and angular quadrature at saved-history causal delay knots; these are integration cuts, not new physical bins.',
        reuse='Stored Phase121 trajectories, Phase122 compact fields and Phase126 additional response. No new EOS, fluid step, frequency/physical-angle bin, time path or duration.',
        conditions='No initially exterior photons and zero initial scalar increments, as in the current initial data. Exterior packets debit the body by the same Killing energy. The compact source port mismatch must remain explicit; this is not a full physical error certificate.',
        gates=dict(quadrature=.002,time=.02,port=1e-11,manufactured=1e-10,delay_inverse_seconds=1e-12,contraction=1.),
        budget=dict(prepare_seconds=30,pilot_seconds=30,production_seconds=100,bound_seconds=40,CPU_threads=1,memory_GB=2,new_native_calls=0),
        stop='Measure the real exterior kernel and one history before production. Keep2x measured forecast within100s. Stop on gate/cap; no automatic quadrature, physical mesh or horizon enlargement.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(compact.__file__),Path(compact.base.__file__),Path(flow.__file__),flow.INPUT/'balanced-20.npz',compact.base.OUT/'audit.json',flow.OUT/'result.json']}))
    write(OUT/'symbolic.json',symbolic());m=Exterior();mu=np.array([.1,.3,.8,.99]);sol,_=m.delay(mu,True);s=np.linspace(0,1,101)
    th=np.arccos(mu)[:,None]*(1-s);impact=m.r0*np.sqrt(1-mu*mu)
    exact=m.r0/C*(1-mu[:,None])-impact[:,None]/C*np.sin(th)/(1+np.cos(th))
    error=float(np.max(abs(sol.sol(s)-exact))/(m.r0/C));assert error<1e-10
    # Independent full-support flat impulse integrals at u=2*r0/c.
    x,w=leg.leggauss(64);errors=[];K=-1.;u=2*m.r0/C
    for v in mu:
        th0=np.arccos(v);bb=m.r0*np.sqrt(1-v*v);th=th0*(x+1)/2;ww=th0*w/2
        delay=m.r0/C*(1-v)-bb/C*np.sin(th)/(1+np.cos(th))
        mass=-C*np.sum(ww*K*np.sin(th)*np.cos(th)/bb**2*(u-delay));stress=-.5*np.sum(ww*K*np.sin(th)**2/bb)
        exact_mass=-C*K/bb**2*((u-m.r0/C*(1-v))*bb**2/(2*m.r0**2)+bb/C*(np.sin(th0)-th0/2-np.sin(2*th0)/4))
        exact_stress=-K/(4*bb)*(th0-v*np.sqrt(1-v*v));errors.extend([abs(mass/exact_mass-1),abs(stress/exact_stress-1)])
    assert max(errors)<1e-10
    d=np.load(compact.base.OUT/'source-128.npz');rin=float(m.bg.edges[-(m.model.bulk.n+m.model.n+1)])
    write(OUT/'check.json',dict(classification='Counterexample candidate',passed=True,flat_delay_relative=error,flat_impulse_relative=max(errors),
        initial_inner_radius_cm=rin,outer_radius_cm=m.r0,flat_minimum_inner_travel_seconds=(m.r0-rin)/C,end_seconds=m.T,
        positive_causal_margin_seconds=(m.r0-rin)/C-m.T,actual_curved_travel_seconds=float(np.diff(m.response.geo(d['edges'][[0,-1]]-m.model.m.RJ)[1])[0]),seconds=time.monotonic()-start))
    signal.alarm(0)


def pilot():
    assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(30)
    m=Exterior();data,row=m.evaluate(128,4,4);np.savez_compressed(OUT/'pilot.npz',**data)
    # Four settings, including a second background with the same cached8/8 kernel.
    forecast=7*row['seconds']+8;upper=2*forecast
    result=dict(classification='Counterexample candidate',row=row,forecast_seconds=forecast,upper_seconds=upper,eligible=upper<100 and row['emitted_port_relative']<1e-11,seconds=time.monotonic()-start)
    write(OUT/'pilot.json',result);print(json.dumps(result),flush=True);signal.alarm(0)
    if result['eligible']:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',budget_seconds=100,forecast_seconds=forecast,upper_seconds=upper,
        settings=[[128,8,8],[128,4,8],[128,8,4],[64,8,8]],bindings={str(p):sha(p) for p in [Path(__file__),OUT/'plan.json',OUT/'pilot.json',OUT/'check.json']}))


def production():
    p=json.loads((OUT/'execution-plan.json').read_text());assert not (OUT/'result.json').exists()
    for f,h in p['bindings'].items():assert sha(f)==h,f
    start=time.monotonic();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(p['budget_seconds']);m=Exterior();data=[];rows=[]
    try:
        for ref,ang,rad in p['settings']:
            d,row=m.evaluate(ref,ang,rad);data.append(d);rows.append(row);np.savez_compressed(OUT/f'exterior-{ref}-a{ang}-r{rad}.npz',**d);write(OUT/f'exterior-{ref}-a{ang}-r{rad}.json',row);print(json.dumps(row),flush=True)
        norm=max(np.max(abs(data[0]['normalized_exterior'])),1e-300)
        controls={key:float(np.max(abs(d['normalized_exterior']-data[0]['normalized_exterior']))/norm) for key,d in zip(['angular_quadrature','radial_quadrature','time'],data[1:])}
        energy_controls={key:float(np.max(abs(d['arrived_energy_erg']-data[0]['arrived_energy_erg']))/max(data[0]['arrived_energy_erg'][-1],1.)) for key,d in zip(['angular_quadrature','radial_quadrature','time'],data[1:])}
        passed=max(controls['angular_quadrature'],controls['radial_quadrature'])<.002 and controls['time']<.02 and energy_controls['angular_quadrature']<.002 and energy_controls['time']<.02 and max(r['emitted_port_relative'] for r in rows)<1e-11
        result=dict(classification='Counterexample candidate',passed=bool(passed),controls=controls,arrival_controls=energy_controls,paths=rows,seconds=time.monotonic()-start,
            actual_exterior_scalar_mass_stress_computed=True,global_potential_enclosed=False,full_source_errors_enclosed=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/'result.json',result);print(json.dumps(result),flush=True)
    finally:signal.alarm(0)

if __name__=='__main__':globals()[sys.argv[1]]()
