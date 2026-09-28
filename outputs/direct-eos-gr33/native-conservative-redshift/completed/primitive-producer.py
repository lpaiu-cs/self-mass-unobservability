"""Counterexample candidate: actual coupled return of the photon work repair.

The saved legacy trajectories remain frozen. This is their separate linear
source increment, using the existing photon/material and GR owners.
"""
from pathlib import Path
from types import SimpleNamespace
import json,resource,sys,time
import numpy as np
import return_native_mixed_gr as owner
import reconcile_native_mass_energy as previous

OUT=Path('native-conservative-redshift167-work');METRIC=OUT/'metric'
reuse=owner.reuse;c=reuse.coupled;inf=owner.inf
LD=np.longdouble;C=owner.C;read,write,sha=owner.read,owner.write,owner.sha
CAPS=dict(prepare=30,check=100,photon_pilot=120,photon_production=1250,
          material_pilot=70,material_production=450,block=300,compact=180,infinity=30)
TOTAL=3600;MAX_SWEEPS=2
CAPS.update(lift_plan=30,lift_check=100,photon_retry=120)
CAPS.update(primitive_plan=30,primitive_check=100,photon_exact=120)
parent_initialize=owner.initialize


def paths(sweep):return OUT/f'sweep-{sweep}/photons',OUT/f'sweep-{sweep}/material'


class ZeroMetric:
    def __init__(self,n):
        self.factor=1.
        z=np.load(inf.SELF/'metric/metric-128-g8.npz')
        self.g={k:(np.zeros_like(v) if k.startswith('delta_') else v.copy()) for k,v in z.items()}
    def view(self,t,clock,side='left'):
        assert np.array_equal(clock,self.g['t']);return self.g


def frequency_source(model,I,field,scale):
    """Incoming-face commutator; frequency ghosts stay in their own ledger."""
    ell=field['delta_log_lapse']*scale
    old=-C*model.cc[:,None]*model.mu*(field['delta_nu_prime']+field['delta_u_prime'])[:,None]*scale
    pos=model.mu>0;donor=np.zeros_like(I);rate=np.zeros((model.n,model.q))
    donor[1:,pos]=I[:-1,pos];donor[:-1,~pos]=I[1:,~pos]
    jump=ell[:-1]-ell[1:]
    rate[1:,pos]=C*model.area[1:-1,None]/model.W[1:,None]*model.mu[pos]*jump[:,None]
    rate[:-1,~pos]=C*model.area[1:-1,None]/model.W[:-1,None]*model.mu[~pos]*jump[:,None]
    new=model.frequency(donor,rate);legacy=model.frequency(I,old)
    variance=model.model.bulk.mu2-model.mu**2
    angular=model.frequency(I,-variance[None,:]*field['delta_lambda_rate'][:,None]*scale)
    packet=new[0]-legacy[0]+angular[0]
    weights=LD(4*np.pi)*model.W[:,None]*model.w
    moments=np.array([np.sum(weights*(new[k]-legacy[k]+angular[k]),dtype=LD) for k in [1,2,3]],float)
    en,ee,work=moments
    ledger=np.array([en,ee,work,-en,work-ee])
    return packet,ledger,max(new[4],legacy[4],angular[4])


def correction(model,t):
    j=np.clip(np.searchsorted(model.t,t,side='left')-1,0,len(model.t)-2)
    f=(t-model.t[j])/(model.t[j+1]-model.t[j]);I=(1-f)*model.I[j]+f*model.I[j+1]
    return frequency_source(model,I,model.redshift_driver.at(t),model.drive_scale)


def correction_rate(model,t):
    d=model.redshift_driver;wave=d.wave;j=np.clip(np.searchsorted(d.times,t,side='left')-1,0,31)
    def derivative(now,x):
        s=(now+np.asarray(x)/C)/d.D
        ut=inf.incident.ETA*d.r0/d.D*inf.incident.pulse(s,1)
        utt=inf.incident.ETA*d.r0/d.D**2*inf.incident.pulse(s,2)
        return [base+np.interp(x,d.tx,(d.born[k][j+1]-d.born[k][j])/(d.times[j+1]-d.times[j]))
                for base,k in zip([ut,utt,utt/C],['U','Ut','Ux'])]
    field=d.at(t)
    try:d.wave=derivative;d.cache={};rate=d.at(t)
    finally:d.wave=wave;d.cache={}
    j=np.clip(np.searchsorted(model.t,t,side='left')-1,0,len(model.t)-2);h=model.t[j+1]-model.t[j]
    f=(t-model.t[j])/h;I=(1-f)*model.I[j]+f*model.I[j+1];It=(model.I[j+1]-model.I[j])/h
    a=frequency_source(model,It,field,model.drive_scale);b=frequency_source(model,I,rate,model.drive_scale)
    return a[0]+b[0],a[1]+b[1],max(a[2],b[2])


def integrated_wave(d,lo,hi,x,moment):
    """Exact polynomial primitive and declared linear Born primitive."""
    x=np.asarray(x,LD);D=LD(d.D);dt=LD(hi-lo)
    a=(LD(lo)+x/C)/D;b=(LD(hi)+x/C)/D
    coeff=np.r_[np.zeros(4,dtype=LD),LD(256)*np.array([1,-4,6,-4,1],LD)]
    def P(z,power=0):
        co=np.polynomial.polynomial.polyint(np.r_[np.zeros(power,dtype=LD),coeff])
        return np.polynomial.polynomial.polyval(np.clip(z,0,1),co)
    A=P(b)-P(a);B=P(b,1)-P(a,1)
    pulse=lambda z:np.where((z>0)&(z<1),np.polynomial.polynomial.polyval(np.clip(z,0,1),coeff),LD(0))
    if moment==0:U=D*A;Ut=pulse(b)-pulse(a)
    else:U=D*D*(B-a*A);Ut=dt*pulse(b)-D*A
    result=[LD(inf.incident.ETA)*LD(d.r0)*v for v in [U,Ut,Ut/C]]
    for j,key in enumerate(['U','Ut','Ux']):
        integ=np.zeros_like(d.born[key][0],dtype=LD)
        for k in range(len(d.times)-1):
            left=max(lo,d.times[k]);right=min(hi,d.times[k+1])
            if right<=left:continue
            slope=(d.born[key][k+1].astype(LD)-d.born[key][k])/(d.times[k+1]-d.times[k])
            v=d.born[key][k].astype(LD)+(left-d.times[k])*slope;h=LD(right-left)
            base=h*v+h*h/2*slope
            integ+=base if moment==0 else LD(left-lo)*base+h*h/2*v+h*h*h/3*slope
        result[j]+=np.interp(np.asarray(x,float),d.tx,np.asarray(integ,float))
    return result


def integrated_source(model,j,hi):
    d=model.redshift_driver;wave=d.wave;lo=model.t[j];fields=[]
    try:
        for moment in [0,1]:
            d.wave=lambda now,x:integrated_wave(d,lo,hi,x,moment)
            d.cache={};fields.append(d.at(hi))
    finally:d.wave=wave;d.cache={}
    It=(model.I[j+1]-model.I[j])/(model.t[j+1]-model.t[j])
    a=frequency_source(model,model.I[j],fields[0],1.)
    b=frequency_source(model,It,fields[1],1.)
    return a[0].astype(LD)+b[0],a[1].astype(LD)+b[1],max(a[2],b[2])


def initialize(sweep):
    parent_initialize(sweep)
    Parent=c.Response
    class Response(Parent):
        def __init__(self,n):
            super().__init__(n);self.redshift_driver=inf.incident.Driver(8);self.source_lifts={}
            self.source_primitives=[(np.zeros_like(self.I[0],dtype=LD),np.zeros(5,dtype=LD),0.)]
        def lift(self,t):
            if not (OUT/'lift-plan.json').exists():return super().lift(t)
            key=(t,self.drive_scale)
            if key in self.source_lifts:return self.source_lifts[key]
            if (OUT/'primitive-plan.json').exists():
                j=np.clip(np.searchsorted(self.t,t,side='left')-1,0,len(self.t)-2)
                while len(self.source_primitives)<=j:
                    k=len(self.source_primitives)-1;a=integrated_source(self,k,self.t[k+1]);b=self.source_primitives[-1]
                    self.source_primitives.append((a[0]+b[0],a[1]+b[1],max(a[2],b[2])))
                a=integrated_source(self,j,t);b=self.source_primitives[j]
                H=(a[0]+b[0])*self.drive_scale;e=(a[1]+b[1])[:3]*self.drive_scale
                Ht,l,err=correction(self,t)
                row=(H/(self.scale*reuse.AMP),Ht,e,l[:3],max(a[2],b[2],err))
                self.source_lifts={key:row};return row
            s,l,e=correction(self,t);ds,dl,de=correction_rate(self,t)
            H=t/4*s;Ht=(s+t*ds)/4
            row=(H/(self.scale*reuse.AMP),Ht,t/4*l[:3],(l[:3]+t*dl[:3])/4,max(e,de))
            self.source_lifts={key:row};return row
        def source(self,t):
            s,l,e=super().source(t);extra,ports,error=correction(self,t)
            return s+extra,l+ports,max(e,error)
    c.Response=Response


def install():
    owner.OUT=OUT;owner.METRIC=METRIC;owner.paths=paths;owner.StageMetric=ZeroMetric
    owner.CAPS=CAPS;owner.TOTAL=TOTAL;owner.initialize=initialize
    owner.install()


def prepare():
    assert not OUT.exists();OUT.mkdir();METRIC.mkdir()
    assert read(previous.OUT/'result.json')['passed']
    for s in range(MAX_SWEEPS+1):
        for p in paths(s):p.mkdir(parents=True)
    files=[Path(__file__),Path(owner.__file__),Path(reuse.__file__),Path(c.__file__),
           Path(previous.__file__),previous.OUT/'result.json',previous.OUT/'commutator.json',
           previous.OUT/'physical-port-128.npz',inf.prior.EV/'coupled-128.npz']
    for n in [64,128]:
        pp=inf.SELF/f'sweep-2/photons-common-gr/steps-{n}-reference-128.npz'
        mm=inf.SELF/f'sweep-2/material-common-gr/steps-{n}-reference-128.npz'
        files += [pp,mm,pp.with_suffix('.json'),mm.with_suffix('.json')]
        with np.load(pp) as p:np.savez_compressed(paths(0)[0]/pp.name,t=p['t'],moments=np.zeros_like(p['moments']),collision_transfer=np.zeros_like(p['collision_transfer']))
        with np.load(mm) as d:np.savez_compressed(paths(0)[1]/mm.name,t=d['t'],history_scaled=np.zeros_like(d['history_scaled']))
    write(OUT/'normalization.json',dict(factor=1.,power_two_exponent=0))
    write(OUT/'photon-conservation-plan.json',dict(scope='Reuse the already established same-equation residual refinement; energy/number weighted stage residual1e-13 and Euclidean1e-14.'))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply the identified spatial photon redshift discrepancy and bin-average angular work to the actual simultaneous photon/thermal/H and free material equations, then measure their GR charge.',
        reason='The measured0.38538erg work discrepancy is larger than the existing0.30736erg mass residual. It cannot be bounded away by the earlier small exterior sectors.',
        equation='Replace the point-gradient spatial frequency forcing by [Lrad,ell*D_E]I, represented by each incoming upwind donor and ell_donor-ell_receiver. Use the same odd conservative frequency owner and retain its exact spectral ghost and work ledgers. Replace mu_center^2 by its declared bin mean in the time-dependent lambda work.',
        reuse='Separate zero-initial linear source response on the same frozen background. No original input replay: additional metric is zero, additive redshift drive is the same actual primary+finite Born. Sum the new response once with the saved response only after controls pass.',
        boundary='The discrete photon state uses donor-center reference energy at its boundary. Conversion to physical outgoing energy and its donor-to-face reference shift remain distinct. No work correction is manually added to the mass or charge.',
        material='Existing conservative free-material response receives only new actual collision transfers. Reciprocally return its B,S,xi and noncollisional E/H; maximum two sweeps, only while actual block defect decreases.',
        gates=dict(stage=1e-12,stage_physical_moment=1e-13,energy_species=1e-8,source_identity=1e-12,time=.02,material_directional=.002,block=.002,GR_quadrature=.002),
        budget=dict(actions=CAPS,total_action_seconds=TOTAL,max_sweeps=MAX_SWEEPS,CPU_threads=1,virtual_GiB=3,new_background_steps=0,new_roots=0,new_rays=0),
        measured_basis='Reuse completed Phase157 per-point, per-step and late material call costs. Measure4/8 photon and two-step material prefixes at the new actual source; admit each full pair only if twice measured remaining cost fits1250s/450s. Total3600s includes all failed actions; later Krylov cost remains an extrapolation.',
        stop='Stop on gate, cost, nondecreasing block defect or two sweeps. No extra resolution, horizon, EOS bank, relaxed threshold, or residual subtraction.',
        limits='Finite retained-table response and selected GR readout. Inner physical energy matching, remaining scalar/pressure work, continuum/native derivative errors, full nonlinear ADM, static comparison and observations remain separate.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def check(sweep):
    install();c.check(sweep);m=c.Response(128);rows=[]
    for t in [m.t[-1]*.31,m.t[-1]*.63]:
        samples={}
        for scale in [0.,1.,-1.,2.]:
            m.drive_scale=scale;s,l,e=correction(m,t);samples[scale]=(s,l)
            number=np.sum(s*m.weights,dtype=LD);energy=np.sum(s*m.weights*m.E,dtype=LD)
            normN=max(np.sum(abs(s)*m.weights,dtype=LD),LD('1e-290'))
            normE=max(np.sum(abs(s)*m.weights*m.E,dtype=LD),LD('1e-290'))
            err=max(float(abs(number+l[0])/normN),float(abs(energy+l[1]-l[2])/normE),e)
            assert err<1e-12,(t,scale,err)
        assert all(np.count_nonzero(v)==0 for v in samples[0.])
        odd=max(float(np.max(abs(samples[s][k]-s*samples[1.][k]))/max(np.max(abs(samples[1.][k])),1e-290)) for s in [-1.,2.] for k in [0,1])
        assert odd<1e-12
        rows.append(dict(time=float(t),source_identity=err,odd_scaling=odd,exact_zero=True))
    write(OUT/f'sweep-{sweep}/source-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows))


def block(sweep):
    owner.block(sweep);p=OUT/f'sweep-{sweep}/block-result.json';r=read(p)
    r['next_sweep_allowed']=r['next_sweep_allowed'] and sweep<MAX_SWEEPS;write(p,r)


def lift_plan(sweep):
    failed=paths(sweep)[0]/'pilot.json';p=read(failed);assert not p['passed']
    write(OUT/'lift-plan.json',dict(classification='Counterexample candidate',
        failure='The actual new source conserves energy/H but early64/128 response moments differ up to23percent. Preserve the rejected prefixes; no production admitted.',
        change='Use the existing affine photon-variable owner with H(t)=t*S_correction(t)/4, Hprime=(S+t*Sprime)/4. The compact C3 pulse starts with a fourth-order zero. This lift removes its leading initial source curvature while leaving the exact continuous forced equation unchanged.',
        derivative='Differentiate the actual linear reference-history interpolation, exact primary polynomial and each stored linear Born field separately. The derivative of stored U interpolation is not silently replaced by independently stored Ut. Use the left interval at closing knots.',
        ledgers='Use the same lifted run owner and actual x=y+H in collisions, energy/H, force, snapshots and outgoing packets. Ghost/work primitives accompany H. No physical source or acceptance threshold changes.',
        budget_seconds=dict(lift_check=100,photon_retry=120),total_unchanged=TOTAL,
        bindings={str(v):sha(v) for v in [Path(__file__),OUT/'registered-producer.py',failed,OUT/'photon_pilot-1-receipt.json']}))
    source=paths(sweep)[0];rejected=source.with_name('photons-rejected')
    assert source.resolve().is_relative_to(OUT.resolve()) and not rejected.exists()
    source.rename(rejected);source.mkdir()


def lift_check(sweep):
    install();initialize(sweep);m=c.Response(128);rows=[]
    for t in [m.t[-1]*.031,m.t[-1]*.31,m.t[-1]*.63]:
        H,Ht,e,et,err=m.lift(t);h=m.t[-1]*1e-6
        plus=m.lift(t+h)[0]*m.scale*reuse.AMP;minus=m.lift(t-h)[0]*m.scale*reuse.AMP
        derivative=float(np.max(abs((plus-minus)/(2*h)-Ht))/max(np.max(abs(Ht)),1e-290))
        packet=H*m.scale*reuse.AMP
        value=m.moments(packet);expected=np.array([-e[0],e[2]-e[1]])
        norms=np.array([np.sum(abs(packet)*m.weights,dtype=LD),np.sum(abs(packet)*m.weights*m.E,dtype=LD)],float)
        moment=float(np.max(abs(value-expected)/np.maximum(norms,1e-290)))
        assert derivative<1e-5 and moment<1e-12 and err<1e-12,(derivative,moment,err)
        rows.append(dict(t=float(t),derivative_relative=derivative,moment_relative=moment))
    write(OUT/('primitive-check.json' if (OUT/'primitive-plan.json').exists() else 'lift-check.json'),dict(classification='Counterexample candidate',passed=True,rows=rows,
        scope='Affine source lift and its derivative/physical moment identities only. Original photon time gate must pass again.'))


def primitive_plan(sweep):
    source=paths(sweep)[0];assert not read(source/'pilot.json')['passed']
    write(OUT/'primitive-plan.json',dict(classification='Counterexample candidate',
        failure='The leading t*S/4 lift reduces the worst early clock discrepancy from23percent to5.23percent but still fails2percent. Preserve it; no production admitted.',
        repair='Use H(t)=integral_0^t S(tau)d tau for the SAME source. Integrate the exact compact polynomial and the existing piecewise-linear Born fields analytically, including the product with the existing17-knot linear photon history. Hprime is the actual source. Reuse the existing affine stage, collision, ghost and packet owner.',
        scope='No new evolution clock, background knot, quadrature approximation, EOS evaluation or physical source. Only the known source primitive changes; both prior prefix failures remain failed.',
        gates=dict(original_time=.02,source_derivative=1e-5,moments=1e-12,energy_H=1e-8),
        budget_seconds=dict(check=100,pilot=120),original_total_seconds=TOTAL,
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'affine-producer.py',OUT/'lift-plan.json',source/'pilot.json',OUT/'photon_retry-1-receipt.json']}))
    rejected=source.with_name('photons-leading-lift-rejected')
    assert source.resolve().is_relative_to(OUT.resolve()) and not rejected.exists()
    source.rename(rejected);source.mkdir()


if __name__=='__main__':
    action=sys.argv[1];sweep=int(sys.argv[2]) if len(sys.argv)>2 else 0
    receipt=OUT/f'{action}-{sweep}-receipt.json';assert action in CAPS and not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                target=OUT/'registered-producer.py' if p==str(Path(__file__)) and (OUT/'registered-producer.py').exists() else p
                assert sha(target)==h,p
            spent=sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'));assert spent+CAPS[action]<=TOTAL
            if sweep>1:assert read(OUT/f'sweep-{sweep-1}/block-result.json')['next_sweep_allowed']
            install()
        if action=='prepare':prepare()
        elif action=='photon_retry':
            assert read(OUT/'lift-check.json')['passed'];c.photon(sweep,True)
        elif action=='photon_exact':
            assert read(OUT/'primitive-check.json')['passed'];c.photon(sweep,True)
        elif action=='primitive_check':lift_check(sweep)
        elif action.startswith(('photon_','material_')):getattr(c,action.split('_')[0])(sweep,action.endswith('pilot'))
        elif action in ['compact','infinity']:getattr(owner,action)(sweep)
        else:globals()[action](sweep)
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,sweep=sweep,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
