"""Counterexample candidate: actual GR return with saved-step source history.

Reuse only the strict-passing prefix needed by the ORIGINAL 188 experiment.
The failed longer 204/205 recovery is retained, not accepted by another gate.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import gc,json,os,resource,sys,time
import numpy as np
import recover_exact_joint_stages as recovery
import read_completed_joint_gr as precision
import evolve_same_solution_gr_return as evolution

base=precision.base;metric=evolution.prior;run=base.run
OUT=Path('native-resolved-return206-work');REC=recovery.OUT
read,write,sha=base.read,base.write,base.sha
LD,AMP,C=base.LD,base.AMP,base.C
CAPS=dict(prepare=180,source=600,fields=600,metric=600,coarse=900,fine=1500,audit=180)


def initialize():
    FunctionType(metric.initialize.__code__,dict(metric.initialize.__globals__,OUT=OUT))()


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(REC/'recovered-64.json')['passed']
    assert 'collision_relative' in read(REC/'fine-receipt.json')['error']
    old=read(evolution.OUT/'result.json');assert not old['passed']
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','gr','metric']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    files=list((REC/'sweep-0').rglob('*.npz'))+[REC/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for p in files:
        dst=OUT/p.relative_to(REC);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [REC/n for n in ['recovered-64.npz','recovered-64.json','accepted-128.npz','fine-receipt.json','proposal-audit.json']]
    files += [base.saved(n) for n in [64,128]]+[evolution.OUT/'result.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='76786f4ab',
        claim='Apply step-resolved same-solution GR input to the original failed188return experiment and test whether its returned B time error falls below the unchanged2percent gate.',
        scope='Exactly188firstT/64,2coarse/4fineexistingsteps,531cells. This is not a shortened204/205recovery pass: theirT/32fine failure remains. Use only earlier strictly accepted reconstruction stages. No new native-fluid high-component evolution, grid, period or parameter path.',
        method='At original accepted step edges use saved gas closing stages with original active-floor rule, recovered photon E/Pr/N moments, and cumulative radial ports from their own original Radau weights. Blend native pressure/background readout at actual time using the same canonical owners. Reuse186/187characteristic GR and exact-center packet lapse; both return clocks then receive the same fine-history metric. Evolve188same compensated equations with actual high/low branch decisions.',
        decision='Preserve source/metric and actual returned-equation, ledger and2percent time failures. Compare to188without inheriting any charge. A pass admits further same-solution closure, not infinity charge or whole-period completion.',
        limits='Piecewise-linear source/metric between2/4existingedges; no uniform temporal interpolation certificate, EOS/derivative certificate, spatial or exterior-scalar closure. The post-result205original-equation audit does not override its stricter failed archive identity.',
        gates=dict(recovery_collision=1e-12,recovery_packet=1e-12,pressure=.002,identity=1e-12,conservation=1e-8,time=.02,quadrature=.002,independent_GR=1e-9,ray=1e-10,stage=1e-12,physical_stage=1e-13),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='205coarse4conditional photon steps took about90s. All required photon stages now saved.188actualreturn2/4steps took57/97s; allow15/25minutes rather than repeated short-cap interruption. Source/field/metric costs expected minutes with10minutes each. Main202run remains untouched.',
        stop='Any original gate or wall cap. No wider gate, frozen physical branch or automatic extra grid. Preserve all results before deciding the next physical lever.',
        original_return_horizon_seconds=old['same_horizon_seconds'],bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',dict(classification='Proven',passed=bool(base.gr.check()['passed']),scope='Reused retarded characteristic polynomial self-check; no physical closure theorem.'))


def recovered(n):
    count=n//32
    if n==64:
        d=dict(np.load(REC/'recovered-64.npz'));logs=read(REC/'recovered-64.json')['rows']
        mom,ports=d['photon_moments'],d['radial_ports']
    else:
        d=dict(np.load(REC/'accepted-128.npz'));assert int(d['step'])>=count
        logs=json.loads(str(d['logs']));mom,ports=d['moments'],d['ports']
    for row in logs[:count]:
        assert row['linear_relative']<1e-14 and max(row['physical_relative'])<1e-13
        assert np.max(row['collision_relative'])<1e-12 and row['packet_relative']<1e-12
    assert len(logs)>=count and len(mom)>=2*count and len(ports)>=2*count
    return count,mom[:2*count],ports[:2*count]


def aligned(coarse,fine,key):
    ids=np.array([int(np.argmin(abs(fine['t']-t))) for t in coarse['t']]);assert np.max(abs(coarse['t']-fine['t'][ids]))<1e-18
    a,b=coarse[key],fine[key][ids]
    return float(np.max(abs(a-b))/max(np.max(abs(b)),LD('1e-290')))


def source():
    initialize();pressure=FunctionType(precision.pressure_function.__code__,dict(precision.pressure_function.__globals__,OUT=OUT))()
    model=base.gr.Response();rows=[];outputs=[]
    for n in [64,128]:
        m=run.owner.Model(n);p=dict(np.load(base.saved(n)));count,mom,ports=recovered(n)
        times=p['actual_step_edges'][:count+1];assert abs(times[-1]-read(OUT/'plan.json')['original_return_horizon_seconds'])<1e-18
        d=base.retained.template()
        for key,v in list(d.items()):
            if v.ndim and len(v)==17:d[key]=np.zeros((len(times),*v.shape[1:]),dtype=v.dtype)
        d['t']=times.copy();material=m.material;coeff=model.coeff(material.rE)
        Eg0,Pg0,Kg0,Er0,Pr0=[coeff[key]*C**4/base.gr.base.G*material.V for key in ['Eg','Pg','Kg','Er','Pr']]
        b=model.model.bulk;weights=4*np.pi*np.r_[b.W,model.model.W][:,None,None]*b.w[None,:,None]*b.d['num']*b.d['Einf']
        gas=[];photons=[];baryons=[];probes=[];mapping=[];references=[];ledgers=[];discard=np.zeros((m.n,4),LD)
        cumulative=np.concatenate([np.zeros((1,2,2),LD),np.cumsum(np.sum((ports*p['joint_stage_weights'][:2*count,None,None]).reshape(count,2,2,2),axis=1),axis=0)],axis=0)
        for i,t in enumerate(times):
            q=np.zeros((4,m.n),LD) if i==0 else p['joint_stage_conserved_scaled'][2*i-1].copy()
            g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su]);off=~material.active(t)
            discard[off]+=g[off]*m.units[off];g[off]=0
            z=np.array([g[:,2]*m.bu,g[:,3]*m.su,g[:,0]*m.eu,g[:,1]*m.nu],LD)
            field=m.geometry(float(t))[0];vol=(3*field[0]+field[2])*AMP
            k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));w=(t-m.t[k])/(m.t[k+1]-m.t[k])
            pp=[];qs=[]
            for j,v in [(k,1-w),(k+1,w)]:
                bank=dict(np.load(base.feedback.old.OUT/f'bank-{m.reference}/point-{j}.npz'))
                pp.append(v*np.array(pressure(material,j,z,field,bank)));qs.append(v*material.point(j)['Q'])
            p0,delta,error=sum(pp);Q=sum(qs);backgroundE=(Q[2]+material.rest*Q[0])/material.a
            nonrest=z[2]*AMP/material.a+(Eg0+Pg0-backgroundE)*vol
            pg,pr=delta*AMP+(Kg0[None]-p0)*vol
            B=z[0]*AMP;gas.append([nonrest,pg,pr]);baryons.append(B);probes.append(error*AMP)
            local=m.local(float(t));reference=(m.pressure(float(t),g)+local['pressure_source'])*m.volume*AMP
            references.append(reference);mapping.append(delta[0]*AMP-reference)
            I=(1-w)*m.I[k]+w*m.I[k+1];Ebg=np.einsum('nqf,nqf->n',I,weights)/material.a;Pbg=np.einsum('nqf,nqf,q->n',I,weights,b.mu2)/material.a
            em,pm=np.zeros((2,m.n),LD) if i==0 else mom[2*i-1,:2]
            phi=field[0]*AMP;lam=field[2]*AMP
            photons.append([em/material.a-Ebg*vol+4*Er0*phi+(Er0+Pr0)*lam,pm/material.a-Pbg*vol+4*Pr0*phi+(3*Pr0-model.ratio4*Er0)*lam])
            if i:
                expected=np.sum(p['joint_stage_weights'][:2*i,None,None]*(p['joint_native_rates_scaled'][:2*i]+p['joint_collision_rates_scaled'][:2*i]),axis=0,dtype=LD)
                actual=g*m.units+discard;ledgers.append((np.sum(abs(actual-expected),axis=0)/np.maximum(np.sum(abs(actual)+abs(expected),axis=0),LD('1e-290'))).astype(float).tolist())
        gas=np.array(gas);photons=np.array(photons);B=np.array(baryons);rest=B*LD(m.model.cx)*LD(C)**2
        d.update(baryon_g=B,gas_nonrest_energy_erg=gas[:,0],nonrest_trace_erg=gas[:,0]-gas[:,2]-2*gas[:,1],nonrest_stress_erg=gas[:,0]-gas[:,2],pressure_volume_erg=gas[:,1],photon_energy_erg=photons[:,0],photon_radial_pressure_erg=photons[:,1],metric_stress_erg=rest+gas[:,0]-gas[:,2]+photons[:,0]-photons[:,1],inner_cumulative_energy_erg=cumulative[:,0,1],outer_cumulative_energy_erg=cumulative[:,1,1])
        norm=max(np.max(np.sum(abs(np.array(references)),axis=-1)),LD('1e-290'))
        row=dict(clock=n,accepted_recovery_steps=count,pressure_probe=float(np.max(np.sum(abs(np.array(probes)),axis=-1))/norm),pressure_mapping=float(np.max(np.sum(abs(np.array(mapping)),axis=-1))/norm),local_material_ledger=ledgers)
        write(OUT/f'source-{n}-check.json',row);np.savez_compressed(OUT/'gr'/f'source-{n}.npz',**d)
        assert row['pressure_probe']<.002 and row['pressure_mapping']<1e-12 and np.max(ledgers)<1e-8,row
        rows.append(row);outputs.append(d);del m,p;gc.collect()
    keys=['baryon_g','gas_nonrest_energy_erg','nonrest_trace_erg','pressure_volume_erg','photon_energy_erg','photon_radial_pressure_erg','metric_stress_erg']
    errors={k:aligned(*outputs,k) for k in keys}
    result=dict(classification='Counterexample candidate',passed=max(errors.values())<.02,source_time=errors,rows=rows,strict_prefix_only=True,old_longer_recovery_failure_preserved=True,final_charge_conclusion='unadjudicated')
    write(OUT/'sources.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def fields():
    assert read(OUT/'sources.json')['passed'];initialize();m=base.gr.Response()
    fn=FunctionType(base.gr.Response.run.__code__,dict(base.gr.Response.run.__globals__,OUT=OUT/'gr'))
    rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]
    fine=dict(np.load(OUT/'gr/fields-128-g8.npz'));coarse=dict(np.load(OUT/'gr/fields-64-g8.npz'));low=dict(np.load(OUT/'gr/fields-128-g4.npz'))
    d=dict(np.load(OUT/'gr/source-128.npz'));direct,coordinate=base.charge.independent.direct(m,d,8)
    errors=dict(time=aligned(coarse,fine,'U'),quadrature=aligned(low,fine,'U'),independent_GR=abs(direct-rows[0]['endpoint_direct'])/max(abs(direct),1e-290))
    result=dict(classification='Counterexample candidate',passed=errors['time']<.02 and errors['quadrature']<.002 and errors['independent_GR']<1e-9,controls=errors,rows=rows,final_charge_conclusion='unadjudicated')
    write(OUT/'fields.json',result);assert result['passed'],result


def metric_run():
    assert read(OUT/'fields.json')['passed'];initialize();m=metric.Lapse()
    ns=dict(metric.metric.Lapse.run.__globals__,OUT=OUT/'metric',prior=SimpleNamespace(OUT=OUT/'gr',new=SimpleNamespace(OUT=OUT/'gr'),centers=metric.constraints.centers))
    fn=FunctionType(metric.metric.Lapse.run.__code__,ns);rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]
    fine=dict(np.load(OUT/'metric/metric-128-g8.npz'));coarse=dict(np.load(OUT/'metric/metric-64-g8.npz'));low=dict(np.load(OUT/'metric/metric-128-g4.npz'))
    controls={name:{k:aligned(other,fine,k) for k in metric.KEYS if k!='delta_lambda_rate'} for name,other in [('time',coarse),('quadrature',low)]}
    # Interval rates live on different clocks; compare their integrals on shared edges.
    for name,other in [('time',coarse),('quadrature',low)]:controls[name]['delta_lambda_rate_integral']=aligned(other,fine,'delta_lambda')
    result=dict(classification='Counterexample candidate',passed=max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002 and max(r['ray_invariant'] for r in rows)<1e-10,controls=controls,rows=rows,interval_derivative_pointwise_convergence_unadjudicated=True,final_charge_conclusion='unadjudicated')
    write(OUT/'metric-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def evolve(n):
    assert read(OUT/'metric-result.json')['passed']
    # Existing exact same equations and stage owners; only metric history changes.
    metric.OUT=OUT
    scope=dict(evolution.initialize.__globals__,OUT=OUT);FunctionType(evolution.initialize.__code__,scope)()
    FunctionType(evolution.evolve.__code__,dict(evolution.evolve.__globals__,OUT=OUT,initialize=lambda:None,Model=scope['Model']))(n)


def audit():
    FunctionType(evolution.audit.__code__,dict(evolution.audit.__globals__,OUT=OUT))()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        elif action=='metric':metric_run()
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
