"""Recover missing photon stages of the SAME accepted long coupled history.

No fluid replay. Reuse205/208 stage histories and solve only the missing
conditional photon blocks. Every new pair must pass the original whole equation.
"""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import complete_joint_photon_recovery as prior

base=prior.base;LD,AMP=base.LD,base.AMP
OUT=Path('native-long-photon213-work');INPUT=Path('native-flux-precision202-work')
read,write,sha=base.read,base.write,base.sha
CAPS=dict(prepare=180,prepare_retry=180,check=300,coarse=7200,fine=10800,audit=300)


def saved(n):return INPUT/f'sweep-1/photons/interval-15-{n}.npz'
def seed_folder(n):return prior.OLD if n==64 else prior.OUT


def prepare():
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True);files=[];reused={}
    assert read(prior.OUT/'result.json')['original_endpoint_and_ledger_passed']
    for name in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material','clock-64','clock-128']:(OUT/name).mkdir(parents=True,exist_ok=True)
    inputs=list((INPUT/'sweep-0').rglob('*.npz'))+[INPUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for p in inputs:
        dst=OUT/p.relative_to(INPUT)
        if dst.exists():assert sha(dst)==sha(p)
        else:os.link(p,dst)
        reused[str(dst)]=sha(p)
    files+=inputs;paths=[]
    for n in [64,128]:
        z=dict(np.load(saved(n)));folder=seed_folder(n);a=dict(np.load(folder/f'recovered-{n}.npz'));row=read(folder/f'recovered-{n}.json')
        assert row['passed'];begin=row['recovered_steps'];count=len(z['actual_step_edges'])-1
        assert np.array_equal(a['times'],z['joint_stage_times'][:2*begin]) and np.array_equal(a['weights'],z['joint_stage_weights'][:2*begin])
        original=dict(np.load(base.saved(n)))
        for k in ['joint_stage_times','joint_stage_weights','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled']:
            assert np.array_equal(z[k][:2*begin],original[k][:2*begin]),k
        assert abs(z['t'][-1]-z['actual_step_edges'][-1])<1e-18 and z['t'][1]==row['horizon_seconds']
        paths.append(dict(clock=n,reused_steps=begin,target_existing_steps=count,missing_steps=count-begin,horizon_seconds=float(z['t'][-1])))
        files += [saved(n),base.saved(n),folder/f'recovered-{n}.npz',folder/f'recovered-{n}.json']
    files += [prior.OLD/'expanded-recovery.py',prior.OUT/'result.json',INPUT/'prefix-result.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='e0cab4459',
        claim='Supply missing time-resolved photon E/Pr/N and original collision/port history on the already accepted commonT*15/16trajectory, so actual full-period GR consumers can reuse it when the running final segment is accepted.',
        decision='Each new conditional photon pair must satisfy the original whole coupled stage equation and physical gates. At every original snapshot require endpoint, reconstructed collision/material ledger and radial/angular port gates. Only the resulting same-solution history may feed GR; no final charge or full-period completion is inferred.',
        reuse='Reuse205coarse4/208fine8 complete recovered steps, their actual moments and ports. Do not replay any native-fluid trajectory or recovered photon prefix. Resume before each remaining step from an atomic photon/history checkpoint.',
        scope='Existing64/128paths through accepted interval15 only; running210final interval is unchanged. Fineprefix includes205strict first5and208whole-equation remaining3. All new stages get original whole-equation checks. Preserve205extra collision-rate identity failure; do not fit or force moments to the archive.',
        gates=dict(conditional=1e-14,conditional_physical=1e-13,original_stage=1e-12,original_physical=1e-13,native_archival_identity=0,endpoint=1e-12,angular=1e-12,radial_port=1e-12,material_ledger=1e-8),
        budgets=CAPS,paths=paths,CPU_threads_per_path=1,virtual_GiB_per_path=6,maximum_parallel_paths=2,
        forecast='205coarse4steps87.4s;208two new photon solves plus three whole-state audits88.0s. Later conditional costs are unmeasured. At20..45s per missing step, coarse roughly36..81minutes and fine69..155minutes; allow2/3hours. These are extrapolations, not completion promises. AvailableWSLmemory was50.6GiB; both6GiB paths plus existing6GiB evolution fit. Stop rather than enlarge caps if inaccurate.',
        stop='Any unchanged whole-equation,endpoint,ledger,port,12linear-refinement or wall gate. Preserve accepted checkpoint and failed pair. No grid, period, new physical path, extra fluid Newton or automatic restart.',
        precision='Original stored constitutive stages and their original operator define this conditional recovery.202high-precision prefix constitutive/ledger controls remain separate; this is not uniform EOS/derivative or every earlier new-precision local-equation certification.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(prior.OLD/'symbolic.json'))


def checkpoint(n,**values):
    p=OUT/f'accepted-{n}.npz';temp=p.with_name(p.stem+'-next.npz');np.savez_compressed(temp,**values);os.replace(temp,p)


def seed(m,z,n):
    folder=seed_folder(n);p=dict(np.load(folder/f'recovered-{n}.npz'));r=read(folder/f'recovered-{n}.json');begin=r['recovered_steps']
    x=np.asarray(p['endpoint_occupation']/(m.scale*AMP),LD)
    expected=z['photon_history_scaled_occupation'][1]
    err=float(np.sum(abs(x*m.scale*AMP-expected))/max(np.sum(abs(expected)),LD('1e-290')))
    assert err<1e-12 and np.array_equal(p['times'],z['joint_stage_times'][:2*begin]) and np.array_equal(p['weights'],z['joint_stage_weights'][:2*begin])
    return begin,x,r['rows'],list(p['photon_moments']),list(p['collision_rates']),list(p['angular']),list(p['radial_ports']),err


def check():
    import gc
    FunctionType(base.base.prior.initialize.__code__,dict(base.base.prior.initialize.__globals__,OUT=OUT))()
    rows=[]
    for n in [64,128]:
        m=base.base.owner.Model(n);z=dict(np.load(saved(n)));begin,x,logs,moments,collisions,packets,ports,error=seed(m,z,n)
        checkpoint(n,x=x,step=begin,t=z['actual_step_edges'][begin],moments=moments,collisions=collisions,packets=packets,ports=ports,logs=np.array(json.dumps(logs)))
        q=dict(np.load(OUT/f'accepted-{n}.npz'))
        assert np.array_equal(q['x'],x) and np.array_equal(q['moments'],moments) and np.array_equal(q['collisions'],collisions) and np.array_equal(q['ports'],ports)
        rows.append(dict(clock=n,reused_steps=begin,checkpoint_value_identity=True,endpoint_relative=error,new_physical_steps=0,new_photon_solves=0));del m,z,q;gc.collect()
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows))


def snapshot(n,m,z,step,photons,collisions,ports,packets):
    t=z['actual_step_edges'][step+1];idx=int(np.argmin(abs(z['t']-t)))
    if abs(z['t'][idx]-t)>1e-18:return None
    count=2*(step+1);actual=photons*m.scale*AMP;expected=z['photon_history_scaled_occupation'][idx]
    endpoint=float(np.sum(abs(actual-expected))/max(np.sum(abs(expected)),LD('1e-290')))
    q=z['conserved_material_history'][idx]/AMP;k=z['energy_offset_reference'][idx]/AMP
    matter=np.column_stack([q[2]-k,q[3],q[0],q[1]])+z['material_floor_discard_history_scaled'][idx]
    weights=z['joint_stage_weights'][:count];target=np.sum(weights[:,None,None]*(z['joint_native_rates_scaled'][:count]+np.array(collisions)),axis=0,dtype=LD)
    ledger=(np.sum(abs(matter-target),axis=0)/np.maximum(np.sum(abs(matter)+abs(target),axis=0),LD('1e-290'))).astype(float).tolist()
    p=np.sum(weights[:,None,None]*np.array(ports),axis=0,dtype=LD);old=z['radial_ports'][idx]
    radial=float(np.max(abs(p-old)/np.maximum(abs(old),LD('1e-290'))))
    angular=float(np.sum(abs(np.array(packets)-z['accepted_angular_luminosity'][:count]))/max(np.sum(abs(z['accepted_angular_luminosity'][:count])),LD('1e-290')))
    row=dict(time=float(t),snapshot_index=idx,steps=step+1,endpoint_relative=endpoint,material_ledger=ledger,radial_port_relative=radial,angular_relative=angular)
    row['passed']=endpoint<1e-12 and max(ledger)<1e-8 and radial<1e-12 and angular<1e-12
    write(OUT/f'clock-{n}'/f'snapshot-{idx:02d}.json',row);assert row['passed'],row
    return row


def run(n):
    assert read(OUT/'restart-check.json')['passed']
    source=(prior.OLD/'expanded-recovery.py').read_text()
    changes=[("if n==128:assert read(OLD/'recovered-64.json')['passed']", "assert read(seed_folder(n)/f'recovered-{n}.json')['passed']"),
        ("target=float(z['t'][1])","target=float(z['t'][-1])"),
        ("zero=np.zeros((m.n,4),LD);zeroJ=sparse.csr_matrix((4*m.n,4*m.n));logs=[];moments=[];collisions=[];packets=[];timings=[];ports=[]", "zero=np.zeros((m.n,4),LD);zeroJ=sparse.csr_matrix((4*m.n,4*m.n));timings=[]\n    begin,x,logs,moments,collisions,packets,ports,seed_error=seed(m,z,n)"),
        ('for step in range(stop):','for step in range(begin,stop):'),
        ("np.savez_compressed(OUT/f'accepted-{n}.npz',","checkpoint(n,"),
        ('logs.append(row);write(', "row['original_joint_equation']=original_equation(m,z,step,t,h,x,gas,pairs,times,cs,ss)\n        row['strict_archival_rate_identity_passed']=bool(np.max(errors)<1e-12)\n        row['snapshot']=snapshot(n,m,z,step,pairs[-1],collisions,ports,packets)\n        logs.append(row);write("),
        ('if np.max(errors)>=1e-12 or packet_error>=1e-12:','if packet_error>=1e-12:'),
        ("assert attempt<11,('Conditional photon residual',err,physical)","if attempt>=11:\n                np.savez_compressed(OUT/f'rejected-linear-{n}.npz',step=step,time=t,step_size=h,x_initial=x,solution=sol,rhs=rhs,residual=residual)\n                raise AssertionError(('Conditional photon residual',err,physical))"),
        ("expected=z['photon_history_scaled_occupation'][1]","expected=z['photon_history_scaled_occupation'][-1]"),
        ('same_material_history_unchanged=True,new_physical_steps=0', "same_material_history_unchanged=True,new_physical_steps=0,reused_steps=begin,new_conditional_steps=stop-begin,seed_endpoint_relative=seed_error,strict_archival_identity_admitted=False,original_equation_audited_steps=list(range(begin,stop))")]
    for a,b in changes:assert source.count(a)==1,(a,source.count(a));source=source.replace(a,b)
    check_original=FunctionType(prior.original_equation.__code__,dict(prior.original_equation.__globals__,OUT=OUT/f'clock-{n}'))
    def original_equation(m,z,step,t,h,x,gas,photons,times,cs,ss):
        try:return check_original(m,z,step,t,h,x,gas,photons,times,cs,ss)
        except BaseException:
            np.savez_compressed(OUT/f'rejected-original-{n}.npz',step=step,time=t,step_size=h,x_initial=x,photon_stage_solution=photons,gas=gas,stage_times=times);raise
    ns=dict(base.run.__globals__,OUT=OUT,saved=saved,seed_folder=seed_folder,seed=seed,checkpoint=checkpoint,original_equation=original_equation,snapshot=snapshot)
    exec(compile(source,__file__,'exec'),ns);(OUT/f'expanded-recovery-{n}.py').write_text(source);ns['run'](n)


def audit():
    rows=[read(OUT/f'recovered-{n}.json') for n in [64,128]];assert all(r['passed'] for r in rows)
    assert rows[0]['horizon_seconds']==rows[1]['horizon_seconds']
    for n in [64,128]:
        records=[read(p) for p in (OUT/f'clock-{n}').glob('snapshot-*.json')];assert records and all(r['passed'] for r in records)
    write(OUT/'result.json',dict(classification='Counterexample candidate',same_solution_accepted_history_recovered=True,rows=rows,original205failure_preserved=True,
        full_declared_period_recovered=False,missing_final_interval=True,GR_source_time_admission_passed=False,physical_final_charge_solved=False,final_charge_conclusion='unadjudicated',full_goal_complete=False))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));base.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action=='prepare_retry':
            assert read(OUT/'prepare-receipt.json')['error']=='AssertionError()'
            write(OUT/'preparation-fix.json',dict(change='Use the existing1e-18clock-coincidence gate for canonical-output versus accumulated-step edges. Their signed discrepancies are plus/minus4.336808689942018e-19. Exact stored stage times/weights/state-prefix comparisons remain bit-identical.',failed_producer_sha256=sha(OUT/'failed-prepare-producer.py')))
        if action not in ['prepare','prepare_retry']:
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:run(64 if action=='coarse' else 128)
        else:globals()['prepare' if action=='prepare_retry' else action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
