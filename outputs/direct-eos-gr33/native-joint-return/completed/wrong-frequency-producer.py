"""Counterexample candidate: exact-center/packet lapse for the SAME joint GR.

Do not infer a final charge from a computed metric. Apply the returned input
to the actual current stage owners before admitting any new coupled evolution.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import gc,json,os,resource,sys,time
import numpy as np
import read_full_incident_joint_gr as prior
import def_native_dynamic_lapse as metric
import verify_native_anisotropic_gr as constraints
import solve_native_incident_self_gr as old
import def_native_feedback_metric_return as compensated

OUT=Path('native-joint-return187-work');GR=prior.OUT/'gr'
read,write,sha=prior.read,prior.write,prior.sha
LD,AMP,C=prior.LD,prior.AMP,prior.C
CAPS=dict(prepare=20,metric=40,apply=45)
CAPS['compensated']=90
CAPS['compensated_finish']=65
CAPS['compensated_identity']=50
KEYS=old.METRIC_KEYS


def initialize():
    prior.run.prior.OUT=OUT;prior.run.prior.initialize()


def prepare():
    assert not OUT.exists();OUT.mkdir();assert read(prior.OUT/'fields.json')['passed']
    (OUT/'metric').mkdir();reused={}
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((prior.OUT/'sweep-0').rglob('*.npz'))+[prior.OUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for src in files:
        dst=OUT/src.relative_to(prior.OUT);os.link(src,dst);reused[str(dst)]=sha(src)
    files += [prior.OUT/'sources.json',prior.OUT/'fields.json']+list(GR.glob('*.npz'))+[prior.saved(n) for n in [64,128]]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None)
        and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='1f628e407',
        claim='Construct the exact cell-integral mass constraint and asymptotically normalized lapse from186same-solution sources and actual Radau angular packets, then apply that returned input to current photon and native material stage owners.',
        decision='Only if the true returned field survives arithmetic and all original mapping/packet/metric gates may a new jointly evolved response be admitted. A field diagnostic alone is not GR closure.',
        scope='The accepted0..T/8 prefix only. Same531cells,64/128 clocks and4/8 quadrature. Incoming driver plus its own returned GR, not addition of independently evolved response states.',
        boundary='Use actual saved stage times and weights and edge lookup, never uniform SDIRK or a clipped late packet. Leading outgoing scalar-vacuum lapse and explicit ADM residual are retained; missing exterior scalar return remains open.',
        arithmetic='Keep return increments separate from background and measure cancellation when combining with the incoming input. Require the applied return to reproduce the separate input component to0.2percent; do not hide complete loss behind a small total-input relative error.',
        gates=dict(time=.02,quadrature=.002,port=1e-12,ray=1e-10,returned_component=.002),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,new_physical_steps=0,
        forecast='186three-time GR fields6.22s; current three lapse reconstructions expected5..25s with40s cap. One actual model constructor and saved-state RHS application expected15..35s with45s cap. No long integration or automatic retries/refinement.',
        stop='Any packet, metric, arithmetic, source or wall cap failure blocks new evolution. Keep failures. No background/producer edit, horizon/cadence/amplitude/gate change or separated-charge addition.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused))
    write(OUT/'symbolic.json',metric.symbolic())


class Lapse(metric.Lapse):
    def boundary(self,d,steps,order):
        t=d['t'];end=t[-1];mu,mw,ids,solution,invariant=self.rays(order,end)
        p=np.load(prior.saved(steps));stage=p['accepted_angular_times'];edges=p['actual_step_edges']
        packets=p['accepted_angular_quadrature_weights'][:,None].astype(LD)*p['accepted_angular_luminosity']
        photon=[];energy=[];norm=[]
        for now in t:
            j=int(np.argmin(abs(edges-now)));assert abs(edges[j]-now)<1e-18;n=2*j
            if not n:photon.append(0.);energy.append(0.);norm.append(0.);continue
            assert np.max(stage[:n])<=now+1e-18
            weight=packets[:n,ids]*mu*mw;ages=(now-stage[:n])/end
            assert np.min(ages)>=-1e-14
            state=solution.sol(np.maximum(ages,0)).reshape(3,len(mu),n).transpose(0,2,1)
            rp=self.r0*state[0];mp=state[1];cc=self.bg.metric(rp.ravel()/self.model.m.R)[3].reshape(rp.shape)
            kernel=state[2]/self.r0-mp*mp/(rp*cc)
            photon.append(float(metric.G/C**4*np.sum(weight*kernel,dtype=LD)))
            energy.append(float(np.sum(weight,dtype=LD)));norm.append(float(np.sum(abs(weight),dtype=LD)))
        energy=np.array(energy);error=float(np.max(abs(energy-d['outer_cumulative_energy_erg']))/max(max(norm),1e-290))
        assert error<1e-12,error
        gx,gw=np.polynomial.legendre.leggauss(order);z=(gx+1)/2
        _,N,b,_,_=self.bg.metric(self.r0/(z*self.model.m.R));kernel=float(np.sum(gw/(2*self.r0*N*b**1.5)))
        return np.array(photon),energy,dict(ray_invariant=invariant,emitted_energy_relative=error,
            unit_ADM_lapse_kernel_per_cm=kernel,only_actual_signed_incremental_packets=True)


def metric_run():
    initialize();m=Lapse()
    scope=dict(metric.Lapse.run.__globals__,OUT=OUT/'metric',
        prior=SimpleNamespace(OUT=GR,new=SimpleNamespace(OUT=GR),centers=constraints.centers))
    fn=FunctionType(metric.Lapse.run.__code__,scope)
    rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]
    fine=dict(np.load(OUT/'metric/metric-128-g8.npz'));controls={}
    for label,n,q in [('time',64,8),('quadrature',128,4)]:
        other=np.load(OUT/'metric'/f'metric-{n}-g{q}.npz')
        controls[label]={key:float(np.max(abs(fine[key]-other[key]))/max(np.max(abs(fine[key])),1e-290)) for key in KEYS}
    result=dict(classification='Counterexample candidate',passed=max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002 and max(r['ray_invariant'] for r in rows)<1e-10,
        controls=controls,rows=rows,exact_center_constraints=True,actual_stage_packet_lapse=True,
        GR_input_applied=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'metric-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


class ReturnedDriver:
    """Same fixed input for both clocks; no separate response superposition."""
    def __init__(self,primary):
        self.primary=primary;self.n=primary.n
        self.g=dict(np.load(OUT/'metric/metric-128-g8.npz'));self.clock=self.g['t']
    def __getattr__(self,name):return getattr(self.primary,name)
    def returned(self,t):
        k=int(np.argmin(abs(self.clock-t)))
        if abs(self.clock[k]-t)<=16*np.finfo(float).eps*self.clock[-1]:t=float(self.clock[k])
        assert self.clock[0]<=t<=self.clock[-1]
        j=int(np.clip(np.searchsorted(self.clock,t,side='left')-1,0,len(self.clock)-2));h=self.clock[j+1]-self.clock[j];w=(t-self.clock[j])/h
        row={k:(1-LD(w))*v[j].astype(LD)+LD(w)*v[j+1] for k,v in self.g.items()
            if k.startswith('delta_') and v.shape==(len(self.clock),self.n)}
        row['delta_u_t']=(self.g['delta_u'][j+1].astype(LD)-self.g['delta_u'][j])/h
        row['delta_lambda_rate']=self.g['delta_lambda_interval_rate'][j].astype(LD)
        row['delta_log_volume']=3*row['delta_u']+row['delta_lambda']
        row['delta_log_areal_radius']=row['delta_u']
        return row
    def at(self,t):
        base=self.primary.at(t);extra=self.returned(t)
        return {k:np.asarray(v,LD)+extra[k] for k,v in base.items()}
    view=prior.run.owner.drive.Driver.view


class ReturnOnly(ReturnedDriver):
    def __init__(self,primary):super().__init__(primary);self.factor=LD(1)
    def at(self,t):return {k:np.asarray(self.factor*self.returned(t)[k],float) for k in self.primary.at(t)}


class ProbeDriver(ReturnedDriver):
    def at(self,t):
        base=self.primary.at(t);extra=self.returned(t)
        return {k:np.asarray(np.asarray(v,LD)+self.factor*extra[k],float) for k,v in base.items()}


def compensated_apply():
    assert not read(OUT/'application-admission.json')['passed'];initialize()
    m=prior.run.owner.Model(128);primary=m.driver;returned=ReturnOnly(primary);probe=ProbeDriver(primary)
    p=np.load(prior.saved(128));rows=[]
    def select(driver):m.driver=driver;m.redshift_driver=driver;m.material.driver=driver
    def relative(delta,reference,units):
        return (np.sum(abs(delta-reference)*units,axis=0)/np.maximum(np.sum(abs(reference)*units,axis=0),LD('1e-290'))).astype(float)
    frequency=m.frequency
    for index in [1,2]:
        returned.factor=LD(1)
        t=float(returned.clock[index]);j=int(np.argmin(abs(p['t']-t)));assert abs(p['t'][j]-t)<1e-18
        g=p['material_history'][j]/AMP;select(primary);directions=[]
        def record(I,omega):directions.append(omega.copy());return frequency(I,omega)
        m.frequency=record;base_source=m.source(t);m.frequency=frequency
        tangent,reset=prior.run.owner.joint.selected_tangent()
        base_native=m.native(t,g);selected=m.native(t,g,tangent=tangent)
        assert np.array_equal(base_native,selected),'Selected original material branches must match the actual state'
        select(returned);reset(False)
        delta_native=m.native(t,np.zeros_like(g),tangent=tangent)
        cursor=0;increments=[];m.branch_crossings=0;m.zero_primary_drifts=0;m.max_drift_ratio=0.
        def increment(I,omega):
            nonlocal cursor
            m.primary_omega=directions[cursor];cursor+=1;increments.append(omega.copy())
            return compensated.Response.frequency(m,I,omega)
        m.frequency=increment;delta_source=m.source(t);m.frequency=frequency
        assert cursor==len(directions)
        physical_crossings=m.branch_crossings
        # Independent high-precision scalar check on every actual donor pair.
        import mpmath as mp
        mp.mp.dps=80;positive_error=0.
        for base,delta in zip(directions,increments):
            actual=compensated.positive_increment(base,delta)
            for a,b,v in zip(base.ravel(),delta.ravel(),actual.ravel()):
                aa,bb=mp.mpf(float(a)),mp.mpf(float(b));exact=max(aa+bb,0)-max(aa,0)
                positive_error=max(positive_error,float(abs(mp.mpf(float(v))-exact)/max(abs(bb),mp.mpf('1e-290'))))
        assert positive_error<1e-12
        c=m.local(t);delta_collision=m.collision(c,np.zeros_like(m.I[0]),np.zeros_like(g),source=True)
        # Select amplification from the current actual rate, never a target charge.
        ratios=np.sum(abs(delta_native)*m.units,axis=0)/np.maximum(np.sum(abs(base_native)*m.units,axis=0),LD('1e-290'))
        magnitude=float(max(ratios));assert 0<magnitude<1e-4
        exponent=int(np.floor(np.log2(1e-4/magnitude)));factor=LD(2)**exponent
        checks=[]
        for scale in [factor/2,factor]:
            probe.factor=scale;select(probe)
            actual=(m.native(t,g)-base_native)/scale
            native_error=relative(actual,delta_native,m.units)
            actual_source=m.source(t)
            # Amplified probes may cross upwind branches. Compare the EXACT
            # finite increment at that amplitude, not an assumed linear scaling.
            select(returned);returned.factor=scale;cursor=0;before=m.branch_crossings
            m.frequency=increment;resolved=m.source(t);m.frequency=frequency
            source_error=float(np.sum(abs((actual_source[0]-base_source[0])-resolved[0])*m.weights*m.E)/max(np.sum(abs(resolved[0])*m.weights*m.E),LD('1e-290')))
            checks.append(dict(amplification=float(scale),native=native_error.tolist(),photon=source_error,
                amplified_frequency_crossings=m.branch_crossings-before))
        row=dict(time=t,checks=checks,frequency_calls=cursor,frequency_branch_crossings=physical_crossings,
            independent_80digit_positive_increment_error=positive_error,
            photon_moment_error=float(delta_source[2]),native_return_over_primary_rate=ratios.astype(float).tolist(),
            actual_native_return_L1=np.sum(abs(delta_native)*m.units,axis=0).astype(float).tolist(),
            actual_photon_return_energy_L1=float(np.sum(abs(delta_source[0])*m.weights*m.E)),
            actual_collision_return_L1=np.sum(abs(delta_collision[1])*m.units,axis=0).astype(float).tolist())
        rows.append(row)
        np.savez_compressed(OUT/f'applied-return-{index}.npz',t=t,material=delta_native,
            photon=delta_source[0],photon_ledger=delta_source[1],collision_photon=c['q'],collision_bound=c['qb'],collision_escape=c['qe'],
            collision_gas=delta_collision[1],current_material=g)
        write(OUT/f'compensated-{index}.json',dict(classification='Counterexample candidate',**row))
        assert max(max(v['native']+[v['photon']]) for v in checks)<.002,row
        assert delta_source[2]<1e-12,row
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,
        actual_current_state_stage_return_applied=True,direct_addition_failure_preserved=True,
        geometry_components_kept_separate=True,new_physical_steps=0,
        limitation='Branch-conditioned stage forcing on two actual states, checked with two resolved probes. This is not a uniform branch-stability theorem, an evolved return, or permission to add independent nonlinear response solutions.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'compensated-result.json',result);print(json.dumps(result),flush=True)


def apply():
    assert read(OUT/'metric-result.json')['passed'];initialize()
    m=prior.run.owner.Model(128);driver=ReturnedDriver(m.driver)
    p=np.load(prior.saved(128));rows=[]
    for j in [1,2]:
        t=float(driver.clock[j]);base=m.driver.at(t);extra=driver.returned(t);total=driver.at(t)
        metrics={key:dict(returned_max=float(np.max(abs(extra[key]))),incident_max=float(np.max(abs(base[key]))),
            lost_relative=float(np.sum(abs(total[key]-base[key]-extra[key]))/max(np.sum(abs(extra[key])),LD('1e-290')))) for key in KEYS}
        rows.append(dict(time=t,components=metrics))
    result=dict(classification='Counterexample candidate',passed=max(v['lost_relative'] for r in rows for v in r['components'].values())<.002,
        rows=rows,scope='Actual input-combination arithmetic before dispatch; no new physical integration.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'application-admission.json',result);print(json.dumps(result),flush=True)
    assert result['passed'],('Returned component lost before true stage application',result)
    # The branch-sensitive native owner must see the TOTAL field and state.
    j=int(np.argmin(abs(p['t']-driver.clock[-1])));g=p['material_history'][j]/AMP;t=float(p['t'][j])
    old_rate=m.native(t,g);old_source=m.source(t)[0]
    m.driver=driver;m.redshift_driver=driver;m.material.driver=driver
    new_rate=m.native(t,g);new_source=m.source(t)[0]
    write(OUT/'applied-stage.json',dict(classification='Counterexample candidate',passed=True,
        native_rate_change_L1=np.sum(abs(new_rate-old_rate)*m.units,axis=0).astype(float).tolist(),
        photon_energy_source_change=float(np.sum(abs(new_source-old_source)*m.Eweight/(m.scale*AMP))),
        actual_full_input_stage_owners_called=True,new_physical_steps=0,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));prior.run.owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            plan=read(OUT/'plan.json')
            for p,h in dict(plan['bindings'],**plan['reused']).items():assert sha(p)==h,p
        (compensated_apply if action.startswith('compensated') else metric_run if action=='metric' else globals()[action])()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
