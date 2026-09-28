"""Finish the actual coupled return, replaying its formerly terminal interval."""
from pathlib import Path
from types import SimpleNamespace
import gc,json,os,resource,sys,time
import numpy as np
import extend_retarded_history as field
import return_complete_history_gr as previous

ROOT=Path('native-full-return249-work');CHECK=len(sys.argv)>1 and sys.argv[1]=='check'
OUT=ROOT/('check' if CHECK else 'full');OLD=previous.OUT;INPUT=field.OUT
base=field.base;complete=field.source;actual=previous.prior;geometry=actual.geometry
read,write,sha,bind,LD=field.read,field.write,field.sha,field.bind,field.LD
saved=complete.saved
CAPS=dict(check=600,prepare=600,metric=3600,coarse=7200,fine=14400,audit=600)


def seed(folder,numbers,metric=False):
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr','metric']:(folder/part).mkdir(parents=True,exist_ok=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    if metric:files.append(OLD/'metric/metric-128-g8.npz')
    files += [OLD/f'sweep-1/photons/interval-14-{n}{ext}' for n in numbers for ext in ['.npz','.json']]
    for p in files:
        dst=folder/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    return files


def initialize():
    geometry.ReturnOnly=actual.StageDriver
    return bind(actual.prior.initialize,OUT=OUT,saved=saved)()


def check():
    assert read(ROOT/'anchor/result.json')['passed'];assert not OUT.exists()
    inputs=seed(OUT,[64],True);Model=initialize();m=Model(64)
    row=m.run(64,'zero-64',14*64//16,'interval-14-64');assert row['passed']
    a=np.load(OUT/'sweep-1/photons/interval-14-64.npz');b=np.load(OUT/'sweep-1/photons/zero-64.npz')
    assert set(a.files)==set(b.files)
    for k in a.files:assert np.array_equal(a[k],b[k]),('Restart altered saved array',k)
    write(ROOT/'restart-regression.json',dict(classification='Counterexample candidate',passed=True,
        every_saved_array_exact=True,arrays=len(a.files),reused_actual_steps=len(a['actual_step_edges'])-1,
        new_physical_steps=0,boundary_replay_required=True,
        bindings={str(p):sha(p) for p in [Path(__file__),*inputs,ROOT/'boundary-extension.json',ROOT/'anchor/result.json']},
        scope='Zero-step restart on the accepted coarse interval14. Actual full-period stages and metric remain required.',final_charge_conclusion='unadjudicated'))


class Lapse(actual.Lapse):
    boundary=bind(actual.Lapse.boundary,prior=SimpleNamespace(saved=saved))


representation=bind(actual.representation,INPUT=complete.OUT)


def prepare():
    assert read(ROOT/'restart-regression.json')['passed']
    assert read(OLD/'controller-status.json')['state']=='completed' and read(OLD/'result.json')['passed']
    assert read(INPUT/'controller-status.json')['state']=='completed' and read(INPUT/'result.json')['GR_return_admitted']
    assert not OUT.exists();files=seed(OUT,[64,128]);counts=[]
    for p in (INPUT/'gr').iterdir():
        if p.suffix in ['.npz','.json']:os.link(p,OUT/'gr'/p.name);files.append(p)
    times=[]
    for n in [64,128]:
        prefix=np.load(OUT/f'sweep-1/photons/interval-14-{n}.npz');high=np.load(saved(n))
        count=len(prefix['actual_step_edges'])-1;counts.append(count);times.append(float(prefix['actual_step_edges'][-1]))
        assert np.array_equal(prefix['actual_step_edges'],high['actual_step_edges'][:count+1])
        assert np.array_equal(prefix['joint_stage_times'],high['joint_stage_times'][:2*count])
        files += [saved(n),OLD/f'recovered-{n}.npz',OLD/f'run-{n}.json']
    assert times[0]==times[1]
    files += [OLD/f'metric/metric-{n}-g{q}.npz' for n,q in field.FIELDS.values()]
    files += [OLD/n for n in ['result.json','metric-result.json','coarse-receipt.json','fine-receipt.json']]
    files += [INPUT/n for n in ['result.json','plan.json','audit-receipt.json']]
    files += [ROOT/n for n in ['restart-regression.json','boundary-extension.json','anchor/result.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the full-period249metric of248same primary history to the actual saved coupled high/low solution and complete its original final period.',
        method='Restart at canonical14/16 with every low state/history/ledger retained. Re-evolve15/16because its formerly terminal source derivative changes on extension, then evolve the final16/16. Capture the actual new photons/ports and append them to the exact stored prefix.',
        metric='Compute the extended actual-stage metric using the same primary full-period angular emission and source derivatives. Retain applied past values only through the restart, after comparison below1e-12; keep the freshly computed right-sided derivative at the formerly terminal15/16point. Record the splice and require all original metric time/quadrature/ray gates after it.',
        gates=read(previous.prior.prior.OUT/'plan.json')['gates'],budgets=CAPS,CPU_affinity=3,virtual_GiB=16,
        reused_actual_steps=counts,restart_time_seconds=times[0],replayed_canonical_intervals=[15],new_canonical_intervals=[16],
        forecast='236coarse111steps measured2258s; only16coarse and32fine stages are needed from14/16. Roughly6..12/12..25minutes if marginal cost holds, late branches unmeasured. Allow2/4hours plus1hour metric; no automatic longer path, finer clock or tolerance change.',
        stop='Any provenance, exact restart, metric splice/original metric, high anchor, branch, whole/physical stage, constitutive, energy/port or paired-time gate. Preserve rejected proposals and actual checkpoints.',
        scope='One full-period compensated feedback iterate of the same retained operator. Uniform EOS/derivative/time/spatial/boundary errors, selfGR fixed point, nonlinear/static/observational/infinity closure remain.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},full_declared_period=True,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(OLD/'symbolic.json'))


def metric():
    geometry.prior.saved=saved
    bind(actual.metric,OUT=OUT,Lapse=Lapse,representation=representation)()
    result=read(OUT/'metric-result.json');cut=read(OUT/'plan.json')['restart_time_seconds'];splices=[]
    for n,q in field.FIELDS.values():
        path=OUT/f'metric/metric-{n}-g{q}.npz';a=dict(np.load(path));old=dict(np.load(OLD/f'metric/metric-{n}-g{q}.npz'))
        count=int(np.count_nonzero(a['t']<=cut+1e-18));assert np.array_equal(a['t'][:count],old['t'][:count]);errors={}
        for k,v in a.items():
            if v.ndim and v.shape[0]==len(a['t']) and k in old:
                norm=max(np.max(np.sum(abs(old[k][:count]),axis=-1)) if v.ndim>1 else np.max(abs(old[k][:count])),1e-290)
                errors[k]=float(np.max(np.sum(abs(v[:count]-old[k][:count]),axis=-1)) if v.ndim>1 else np.max(abs(v[:count]-old[k][:count])))/norm
                v[:count]=old[k][:count]
            elif k=='delta_lambda_interval_rate':v[:count-1]=old[k][:count-1]
        row=dict(clock=n,order=q,reused_metric_points=count,relative=errors);splices.append(row)
        assert max(errors.values())<1e-12,('Past metric differs materially',row)
        np.savez_compressed(path,**a)
    fine=dict(np.load(OUT/'metric/metric-128-g8.npz'));keys=[k for k in geometry.KEYS if k!='delta_lambda_rate']+['delta_u_t','actual_delta_lambda_rate']
    controls={name:{k:base.endpoint.aligned(dict(np.load(OUT/f'metric/metric-{n}-g{q}.npz')),fine,k) for k in keys} for name,n,q in [('time',64,8),('quadrature',128,4)]}
    assert max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002,controls
    result.update(causal_splice=splices,after_splice_controls=controls,restart_time_seconds=cut,
        formerly_terminal_derivative_recomputed=True,actual_full_period_return_executed=False)
    write(OUT/'metric-result.json',result)


def evolve(n):
    assert read(OUT/'metric-result.json')['passed'];Model=initialize()
    prefix=np.load(OUT/f'sweep-1/photons/interval-14-{n}.npz');count=len(prefix['joint_stage_times'])
    recovered=np.load(OLD/f'recovered-{n}.npz');moments=list(recovered['photon_moments'][:count]);ports=list(recovered['radial_ports'][:count])
    oldrow=read(OLD/f'run-{n}.json');cut=float(prefix['actual_step_edges'][-1])
    native=[v for v in oldrow['anchor_checks'] if v['time']<=cut+1e-18]
    branches=[v for v in oldrow['branch_checks'] if v['time']<=cut+1e-18]
    intervals=oldrow['intervals'][:14];run=Model.run;stage=run.__globals__['stages'];boundary=Model.boundary_ports
    def observed(m,t,h,x,g,lus):
        pair,mechanical=stage(m,t,h,x,g,lus)
        for v in pair:
            xx=v[0];moments.append(np.array([np.sum(xx*m.Eweight,axis=(1,2)),np.sum(xx*m.Eweight*m.model.bulk.mu2[None,:,None],axis=(1,2)),np.sum(xx*m.Nweight,axis=(1,2))])*actual.AMP)
        np.savez_compressed(OUT/f'last-pair-{n}.npz',time=t,step=h,photons=[v[0] for v in pair],gas=[v[1] for v in pair],x_initial=x,g_initial=g)
        write(OUT/f'capture-{n}.json',dict(actual_steps=len(moments)//2,new_steps=(len(moments)-count)//2))
        return pair,mechanical
    def port(m,t,x):
        value=boundary(m,t,x);ports.append(value.copy()*actual.AMP);return value
    Model.boundary_ports=port;Model.run=bind(run,stages=observed);restart=f'interval-14-{n}'
    for j in [15,16]:
        m=Model(n);label=f'interval-{j}-{n}';row=m.run(n,label,j*n//16,restart);assert row['passed']
        native+=m.anchor_checks;branches+=m.branch_checks;intervals.append(dict(row))
        path=OUT/f'sweep-1/photons/{label}.npz';z=np.load(path)
        np.savez_compressed(OUT/f'recovered-{n}.npz',times=z['joint_stage_times'],weights=z['joint_stage_weights'],photon_moments=moments,radial_ports=ports,collision_rates=z['joint_collision_rates_scaled'],angular=z['accepted_angular_luminosity'],endpoint_occupation=z['restart_x']*m.scale*actual.AMP)
        checks=dict(newton=m.newton_iterations,stages=m.stage_log);restart=label;del m;gc.collect()
    dst=OUT/f'sweep-1/photons/return-{n}.npz';os.link(path,dst);os.link(path.with_suffix('.json'),dst.with_suffix('.json'))
    for k in ['actual_step_edges','joint_stage_times','joint_stage_weights','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled','joint_discard_rates_scaled','conserved_material_history','material_floor_discard_history_scaled','material_history','photon_history_scaled_occupation','radial_ports','collision_transfer','accepted_angular_times','accepted_angular_luminosity','accepted_angular_quadrature_weights','moments','t']:
        assert np.array_equal(z[k][:len(prefix[k])],prefix[k]),('Restart prefix changed',k)
    assert np.array_equal(np.array(moments)[:count],recovered['photon_moments'][:count]) and np.array_equal(np.array(ports)[:count],recovered['radial_ports'][:count])
    anchor=np.load(saved(n));assert np.array_equal(z['actual_step_edges'],anchor['actual_step_edges']) and np.array_equal(z['joint_stage_times'],anchor['joint_stage_times'])
    audit,_,_,_=geometry.prior.run.verify(dst,dst);expected=np.sum(z['joint_stage_weights'][:,None,None]*ports,axis=0,dtype=LD)
    error=float(np.max(abs(expected-z['radial_ports'][-1])/np.maximum(abs(z['radial_ports'][-1]),LD('1e-290'))));assert error<1e-12,error
    row.update(audit=audit,intervals=intervals,anchor_checks=native,branch_checks=branches,
        actual_stage_photon_moments_captured=True,captured_radial_port_relative=error,
        maximum_true_stage=max(v[-1]['relative'] for v in checks['newton']),maximum_true_physical_stage=max(max(v[-1]['moments']) for v in checks['newton']),
        same_saved_stage_equation=True,actual_return_time_evolved=True,full_declared_period=True,
        self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/f'run-{n}.json',row);write(OUT/f'checks-{n}.json',dict(classification='Counterexample candidate',**checks))


def audit():
    bind(actual.audit,OUT=OUT)();r=read(OUT/'result.json')
    r.update(full_declared_period=True,causal_restart_with_last_interval_replay=True,final_charge_conclusion='unadjudicated',full_goal_complete=False);write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;ROOT.mkdir(exist_ok=True);receipt=ROOT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));actual.prior.joint.previous.original.inf.incident.native.deadline(CAPS[action]);start=time.monotonic();error=None
    try:
        if action not in ['check','prepare']:
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        evolve(64 if action=='coarse' else 128) if action in ['coarse','fine'] else globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
