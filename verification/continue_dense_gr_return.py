"""Continue225 with the original canonical model lifetime and unchanged gates."""
from pathlib import Path
from types import FunctionType
import gc,json,os,resource,sys,time
import numpy as np
import evolve_dense_same_history_return as prior

OUT=Path('native-interval-return226-work');OLD=prior.OUT
read,write,sha,LD=prior.read,prior.write,prior.sha,prior.LD
CAPS=dict(prepare=180,coarse=1800,fine=3600,audit=180)


def prepare():
    assert read(OLD/'controller-status.json')['state']=='failed'
    assert 'Saved-stage anchor precision' in read(OLD/'coarse-receipt.json')['error']
    assert read(OLD/'metric-result.json')['passed'];assert not OUT.exists();OUT.mkdir()
    files=[]
    for src in list((OLD/'sweep-0').rglob('*.npz'))+list((OLD/'metric').glob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','metric-result.json']]:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    for name in ['photons','material']:(OUT/'sweep-1'/name).mkdir(parents=True)
    files += [OLD/n for n in ['coarse-receipt.json','coarse.stderr.log','plan.json','sweep-1/photons/return-64-checkpoint.npz']]
    files += [prior.saved(n) for n in [64,128]]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the same accepted224dense GR to actual same-equation stages across the canonical background transition. Restore the original producer model lifetime, not a relaxed native gate.',
        cause='225accepted8coarse actual steps throughT/16; first next-interval native anchor failed3.8747131e-12>1e-12.214already demonstrated that original saved material rates require a fresh model at each canonical continuation.',
        method='Use225unchanged equations/inverse/metric and184exact continuation: finish each canonical interval, instantiate the next model and restore every physical state/history/ledger. Require first coarse checkpoint arrays to exactly reproduce225. Reuse the accepted metric without recomputation.',
        replay='225raw checkpoint lacks the running linear/constitutive/balance maxima needed by the original complete restart serializer. Replay its8steps once (~172s measured) to recover them and require exact physical checkpoint identity. Do not invent missing diagnostics or rerun the long material/photon/GR histories.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,virtual_GiB=8,
        forecast='225eight steps172.40s including model setup.15/29steps and2models roughly6..12/12..25min; allow30/60min. Later branches unmeasured. No new grid,clock,period or scientific tolerance.',
        decision='Original10channel2percent time comparison and same-solution ledgers decide thisT/8single feedback iterate. Also retain originalT/64comparison even if global audit fails.225failure and188local failure remain.',
        stop='Any existing stage/physical/native/constitutive/branch/port/time gate or cap. Preserve all interval endpoints and rejected proposals.218/223/224remain immutable.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as sp
    x,d,a,b=sp.symbols('x d a b');assert sp.expand(a*(x+d)+b-(a*x+b)-a*d)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Same-branch affine difference identity only; canonical lifetime and restart are tested on actual arrays, not a physical closure theorem.'))


def evolve(n):
    Model=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    native=[];branches=[];intervals=[];restart=None
    for j in [1,2]:
        m=Model(n);label=f'interval-{j}-{n}';row=m.run(n,label,j*n//16,restart)
        path=OUT/f'sweep-1/photons/{label}.npz'
        checks=dict(newton=m.newton_iterations,stages=m.stage_log)
        native+=m.anchor_checks;branches+=m.branch_checks
        write(OUT/f'interval-{j}-{n}-checks.json',dict(**checks,anchor=m.anchor_checks,branches=m.branch_checks))
        assert row['passed'];intervals.append(dict(row))
        if n==64 and j==1:
            old=np.load(OLD/'sweep-1/photons/return-64-checkpoint.npz');new=np.load(path.with_name(label+'-checkpoint.npz'))
            for key in old.files:assert np.array_equal(old[key],new[key]),('225checkpoint changed',key)
            write(OUT/'prefix-replay-identity.json',dict(classification='Counterexample candidate',passed=True,all_original_checkpoint_arrays_exact=True,actual_steps=8,old_sha256=sha(OLD/'sweep-1/photons/return-64-checkpoint.npz')))
        if restart:
            previous=np.load(OUT/f'sweep-1/photons/{restart}.npz');current=np.load(path)
            for key in ['joint_stage_times','joint_stage_weights','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled','joint_discard_rates_scaled','conserved_material_history','material_floor_discard_history_scaled','material_history','photon_history_scaled_occupation','radial_ports','collision_transfer','accepted_angular_times','accepted_angular_luminosity','accepted_angular_quadrature_weights','actual_step_edges','moments','t']:
                assert np.array_equal(previous[key],current[key][:len(previous[key])]),('Restart changed prefix',key)
        restart=label;del m;gc.collect()
    dst=OUT/f'sweep-1/photons/return-{n}.npz';os.link(path,dst)
    os.link(path.with_suffix('.json'),dst.with_suffix('.json'))
    z=np.load(dst);anchor=np.load(prior.saved(n))
    assert np.array_equal(z['actual_step_edges'],anchor['actual_step_edges']) and np.array_equal(z['joint_stage_times'],anchor['joint_stage_times'])
    audit,_,_,_=prior.geometry.prior.run.verify(dst,dst)
    row.update(audit=audit,intervals=intervals,anchor_checks=native,branch_checks=branches,
        maximum_true_stage=max(v[-1]['relative'] for v in checks['newton']),maximum_true_physical_stage=max(max(v[-1]['moments']) for v in checks['newton']),
        same_saved_stage_equation=True,actual_return_time_evolved=True,self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/f'run-{n}.json',row);write(OUT/f'checks-{n}.json',dict(classification='Counterexample candidate',**checks))


def audit():
    # Keep the old local failure visible even when the global comparison fails.
    try:FunctionType(prior.original.audit.__code__,dict(prior.original.audit.__globals__,OUT=OUT))()
    finally:
        if (OUT/'result.json').exists():
            r=read(OUT/'result.json');first=[]
            for n,count in [(64,2),(128,4)]:
                z=np.load(OUT/f'sweep-1/photons/return-{n}.npz');first.append(z['joint_stage_conserved_scaled'][2*count-1])
            errors=(np.sum(abs(first[0]-first[1]),axis=-1)/np.maximum(np.sum(abs(first[1]),axis=-1),LD('1e-290'))).astype(float).tolist()
            r.update(original188failure_preserved=True,original225failure_preserved=True,original_T64_endpoint_conserved_relative=errors,limited_to_one_feedback_iterate=True)
            write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
