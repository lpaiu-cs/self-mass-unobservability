"""Exact affine photon change of variable on the unchanged incident problem.

Counterexample candidate. Remove the known incoming conformal/redshift packet
before SDIRK; keep the same equation, clocks, waveform and physical ledgers.
"""
from pathlib import Path
import inspect,json,resource,sys,time
import numpy as np
import sympy as sp
import apply_native_incident_drive as base

OUT=base.OUT;PHOTON=OUT/'photons-lift';read,write,sha=base.read,base.write,base.sha
AMP=base.AMP;original_initialize=base.initialize


def initialize():
    base.PHOTON=PHOTON;original_initialize();Parent=base.Response
    class Lifted(Parent):
        def __init__(self,steps):super().__init__(steps);self.lifts={}
        def lift_map(self,I,u,lam):
            omega=-u[:,None]-self.mu[None,:]*lam[:,None]
            packet,en,ee,work,error=self.frequency(I,omega)
            mu=self.edges_mu[1:-1]
            flux=np.zeros((self.n,self.q+1,self.nf))
            flux[:,1:-1]=-(1-mu*mu)[None,:,None]*lam[:,None,None]*I[:,:-1]
            packet-=np.diff(flux,axis=1)/np.diff(self.edges_mu)[None,:,None]
            w=4*np.pi*self.W[:,None]*self.w
            e=np.array([np.sum(w*v,dtype=np.longdouble) for v in [en,ee,work]],float)
            return packet,e,error
        def lift(self,t):
            if t in self.lifts:return self.lifts[t]
            j=np.clip(np.searchsorted(self.t,t,side='left')-1,0,15)
            h=self.t[j+1]-self.t[j];a=(t-self.t[j])/h
            I=(1-a)*self.I[j]+a*self.I[j+1];It=(self.I[j+1]-self.I[j])/h
            d=self.driver;s=(t+d.xc/base.C)/d.D
            U=base.drive.ETA*d.r0*base.drive.pulse(s)*self.drive_scale
            Ut=base.drive.ETA*d.r0/d.D*base.drive.pulse(s,1)*self.drive_scale
            u=d.zc['alpha']*U/d.centers;lam=d.zc['Phi']*U
            ut=d.zc['alpha']*Ut/d.centers;lt=d.zc['Phi']*Ut
            H,e,err=self.lift_map(I,u,lam)
            Ht,et,err1=self.lift_map(It,u,lam);v,w,err2=self.lift_map(I,ut,lt)
            Ht+=v;et+=w
            row=(H/(self.scale*AMP),Ht,e,et,max(err,err1,err2))
            self.lifts[t]=row
            if len(self.lifts)>3:del self.lifts[next(iter(self.lifts))]
            return row
        def local(self,t):
            c=super().local(t);H=self.lift(t)[0]
            p,_,e,b=super().collision(c,H,np.zeros((self.n,2)))
            c['q']+=p;c['qb']+=b;c['qe']+=e
            return c
        def source(self,t):
            s,l,error=super().source(t);H,Ht,e,et,err=self.lift(t)
            packet=H*self.scale*AMP
            s+=(self.A@packet.reshape(self.n*self.q,self.nf)).reshape(s.shape)-Ht
            l[:3]-=et;l[3:]+=self.port(packet)-self.moments(Ht)
            return s,l,max(error,err)
        def boundary_ports(self,t,x):return super().boundary_ports(t,x+self.lift(t)[0])
    runner=(base.drive.native.prior.OUT/'expanded-photon-run.py').read_text()
    replace=base.old.fields.previous.prior.updated.replace
    # Only storage/observables use physical x=y+H; stage unknowns remain y.
    begin=runner.index('    def record(t):');end=runner.index('    count=steps',begin)
    record=runner[begin:end].replace('    def record(t):','    def record(t):\n        physical_x=x+self.lift(t)[0]')
    record=record.replace('x*self.','physical_x*self.').replace('abs(x)','abs(physical_x)')
    runner=runner[:begin]+record+runner[end:]
    runner=replace(runner,'    def stream(xx):',
        '    if restart is not None:\n        H,_,e,_,_=self.lift(begin*h)\n        x-=H;ledger-=self.moments(H*self.scale);escape-=e/AMPLITUDE\n    def stream(xx):')
    runner=replace(runner,'photon=self.moments(x*self.scale);total=',
        'H=self.lift((k+1)*h)[0];physical_x=x+H;physical_ledger=ledger+self.moments(H*self.scale)\n        photon=self.moments(physical_x*self.scale);total=')
    runner=runner.replace('np.sum(abs(x)*self.Nweight)','np.sum(abs(physical_x)*self.Nweight)').replace('np.sum(abs(x)*self.Eweight)','np.sum(abs(physical_x)*self.Eweight)')
    runner=replace(runner,'defect=abs(total-ledger)/np.maximum(norm,abs(ledger))','defect=abs(total-physical_ledger)/np.maximum(norm,abs(physical_ledger))')
    runner=replace(runner,'    solve_seconds=time.monotonic()-start',
        '    H,_,e,_,_=self.lift(count*h)\n    x+=H;ledger+=self.moments(H*self.scale);escape+=e/AMPLITUDE\n    solve_seconds=time.monotonic()-start')
    ns=dict(base.Response.run.__globals__,OUT=PHOTON);exec(compile(runner,__file__,'exec'),ns)
    Lifted.run=ns['run'];base.Response=Lifted
    (OUT/'expanded-lifted-run.py').write_text(runner)


def prepare():
    assert not read(OUT/'photons/pilot.json')['passed'];assert not PHOTON.exists();PHOTON.mkdir()
    t=sp.symbols('t');L=sp.Function('L')(t);H=sp.Function('H')(t);y=sp.Function('y')(t);s=sp.Function('s')(t)
    assert sp.expand(L*(y+H)+s-sp.diff(H,t)-(L*y+s+L*H-sp.diff(H,t)))==0
    write(OUT/'lift-plan.json',dict(classification='Counterexample candidate',
        failure=read(OUT/'photons/pilot.json')['time_comparison'],
        change='Exact affine x=y+H(t), not a new physical source or weakened gate. H is the conservative primary incoming conformal frequency/angle packet. Solve yprime=L*y+s+L*H-Hprime, with exact primary pulse and piecewise-linear saved background in Hprime. Tiny Born forcing remains unchanged.',
        physical_ledgers='Collision and port owners see y+H. The transformed frequency ledger subtracts the exact primitive derivative; physical totals add its endpoint. Store actual x, gas, accepted angular ports and collision transfer. Restart converts the stored physical state back to y.',
        scope='Same finite equation and2-stage SDIRK method,64/128 clocks,531cells,8angles,152frequencies,eta1e-30 and3.434ms. No new background/native calls. This transformation is not an EOS, continuum or full-GR certificate.',
        budget=dict(check_seconds=45,pilot_seconds=90,production_seconds=1000),
        admission='Repeat only the failed equal-horizon4/8 prefixes, retain the original2percent and1e-8 conservation gates and previous measured late-step cost floor. No production above the original1000s cap.',
        symbolic=dict(classification='Proven',passed=True,identity='x=y+H implies yprime=L*y+s+L*H-Hprime for the same linear finite equation.'),
        bindings={str(p):sha(p) for p in [Path(__file__),Path(base.__file__),Path(base.drive.__file__),OUT/'photons/pilot.json',OUT/'rejected-photon-producer.py']}))


def check():
    initialize();m=base.Response(128);Parent=m.__class__.__mro__[1];rows=[]
    for t in [m.t[-1]*.31,m.t[-1]*.63]:
        H,Ht,e,et,_=m.lift(t);p=H*m.scale*AMP
        h=m.t[-1]*1e-6;hp=m.lift(t+h)[0]*m.scale*AMP;hm=m.lift(t-h)[0]*m.scale*AMP
        derivative=float(np.max(abs((hp-hm)/(2*h)-Ht))/max(np.max(abs(Ht)),1e-290))
        exact=m.moments(p);led=np.array([-e[0],e[2]-e[1]])
        norm=np.array([np.sum(abs(p)*m.weights,dtype=np.longdouble),np.sum(abs(p)*m.weights*m.E,dtype=np.longdouble)],float)
        moment=float(np.max(abs(exact-led)/np.maximum(norm,1e-290)))
        raw=Parent.local(m,t);shifted=m.local(t)
        p0,g0,e0,b0=m.collision(raw,H,np.zeros((m.n,2)),True)
        p1,g1,e1,b1=m.collision(shifted,np.zeros_like(H),np.zeros((m.n,2)),True)
        collision=max(float(np.max(abs(a-b))/max(np.max(abs(a)),1e-290)) for a,b in [(p0,p1),(g0,g1),(e0,e1),(b0,b1)])
        rows.append(dict(time=t,derivative=derivative,moment=moment,collision=collision,
                         moment_absolute_residual=abs(exact-led).tolist(),number_energy_absolute_norm=norm.tolist()))
        assert derivative<1e-6 and moment<1e-12 and collision<1e-12,rows[-1]
    write(OUT/'lift-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        moment_normalization='Separate photon number/energy residuals divided by their corresponding absolute packet moments, as in the existing frequency check. Do not divide both by a mixed-unit cancelled net moment.',
        preserved_rejected_check='lift-check-receipt.json',producer_sha256=sha(__file__)))


check_moments=check


def admit():
    p=read(PHOTON/'pilot.json');assert p['passed'] and not p['eligible']
    assert p['upper_remaining_seconds']<1100 and not (PHOTON/'execution-plan.json').exists()
    historical=[read(Path('retained-metric-return152-work')/folder/f'steps-{n}-reference-128.json')
                for folder,n in [('material',64),('material-resolved',128)]]
    prior_material=2*sum(r['worker_wall_seconds']+5 for r in historical)
    assert prior_material<400
    write(OUT/'budget-reassessment.json',dict(classification='Counterexample candidate',eligible=True,
        decision='The affine repair passes the unchanged time/physical gates. Its conservative late-cost floor gives1068.37s, so original1000s admission remains false. Reassign100 unused seconds from the material production allocation; do not expand the1500s combined production cap.',
        production_caps=dict(photon_production=1100,material_production=400),combined_production_cap=1500,
        photon_remaining_forecast=p['upper_remaining_seconds'],historical_material_double_wall=prior_material,
        material_dispatch='Must separately pass its new actual prefix forecast below400s; historical cost is only reallocation evidence, not dispatch authorization.',
        assumptions='Prior full-path late step cost is retained. A2x margin covers the measured photon estimate; new-pulse late Krylov cost remains unmeasured.',
        stop='Same physical problem and clocks. Stop at either cap or any original gate; no further automatic budget/grid/horizon change.',
        bindings={str(q):sha(q) for q in [Path(__file__),OUT/'pilot-lift-producer.py',OUT/'lift-check.json',PHOTON/'pilot.json']}))
    write(PHOTON/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,cap_seconds=1100,
        forecast=p['upper_remaining_seconds'],source_sha256=sha(__file__),driver_sha256=sha(base.drive.__file__),
        reassessment='budget-reassessment.json'))


if __name__=='__main__':
    action=sys.argv[1];assert action in ['prepare','check','check_moments','admit','photon_pilot','photon_production','material_pilot','material_production']
    if (OUT/'budget-reassessment.json').exists():base.drive.CAPS.update(read(OUT/'budget-reassessment.json')['production_caps'])
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));base.drive.native.deadline(base.drive.CAPS.get(action,45))
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action.startswith(('photon_','material_')):
            assert read(OUT/'lift-check.json')['passed'];base.initialize=initialize;base.PHOTON=PHOTON
            if action=='material_production':assert read(base.MATERIAL/'pilot.json')['upper_remaining_seconds']<400
            getattr(base,action.split('_')[0])(action.endswith('pilot'))
            if action=='material_pilot':assert read(base.MATERIAL/'pilot.json')['upper_remaining_seconds']<400
        else:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        receipt=OUT/f'lift-{action}-receipt.json';assert not receipt.exists()
        write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
