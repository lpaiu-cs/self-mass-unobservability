"""Capture missing stages from the original bounded coupled interval, unchanged."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import continue_full_material_front as prior
import complete_joint_photon_recovery as recovery
import recover_remaining_joint_photons as snapshots

OUT=Path('native-joint-capture222-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
LD,AMP,owner=prior.LD,prior.AMP,prior.owner
CAPS=dict(prepare=180,fine=900)


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[]
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    (OUT/'clock-128').mkdir();(OUT/'stage-pairs').mkdir()
    sources=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    sources += [OLD/f'sweep-1/photons/restored-128{ext}' for ext in ['.npz','.json']]
    for src in sources:os.link(src,OUT/src.relative_to(OLD))
    files+=sources+[OLD/n for n in ['fine-receipt.json','pilot-128-result.json','restart-check.json','expanded-full-stages.py']]
    files += [OLD/'sweep-1/photons/pilot-128.npz',recovery.OUT/'recovered-128.npz',recovery.OUT/'recovered-128.json',Path('native-archival-photon220-work/refinement-audit/result.json'),Path('native-original-state221-work/endpoint-8.json')]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Resolve the missing same-solution photon history by capturing its actual original coupled stages in the one short interval that conditional recovery cannot reproduce at the original port gate.',
        reason='Exact original photon input and11extra final-block corrections failed to resolve4.55e-12port difference. A further endpoint control could not reproduce archived collision bits even from the selected restart guide. Original184actual8coupled steps took187.67s, less than220conditional interval250.96s. Replaying this small accepted interval has a concrete missing-data purpose; do not replay the long accepted trajectory.',
        method='Use184same initialize, exact restored-128 input, stage rule, equations, original solver and n128limit8. Wrap only accepted-stage returns to save actual photon/gas pairs and E/Pr/N moments. Wrap the existing boundary call to copy its actual returned port once, without any extra call or equation/source mutation.',
        controls='Require all archived physical material/photon/native/collision/energy/port histories and exact internal endpoint arrays to reproduce184pilot-128 value bits. Combine existing208first8stage records with the16new actual stage samples, then apply the unchanged endpoint/material/angular/radial gates. Mismatch stops; no approximate old-state relabelling or port fitting.',
        gates=dict(stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,original_physical_array_identity=0),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,max_Newton=3,
        forecast='Measured original8actualsteps187.67s; capture8compressed photon pairs adds unmeasured IO. Allow900s (15minutes),8replayed original steps, no new horizon or path. Main218continues independently and its source/plan stays fixed.',
        stop='Any original integration gate, physical array mismatch, recovered-history gate or cap. No automatic full-horizon replay.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    from fractions import Fraction as F
    for k in range(3):assert F(3,4)*F(1,3)**k+F(1,4)==F(1,k+1)
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Original Radau moments0..2. Observation copies do not assert convergence or physical certification.'))


def fine():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    run=owner.Model.run;stage=run.__globals__['stages'];moments=[];ports=[];times=[];weights=[];gas_records=[]
    assert inspect.getsource(stage)==(OLD/'expanded-full-stages.py').read_text(),'Original stage producer changed'
    def observed(m,t,h,x,g,lus):
        pair,mechanical=stage(m,t,h,x,g,lus)
        index=len(moments)//2
        np.savez_compressed(OUT/f'stage-pairs/step-{index+8:03d}.npz',time=t,step=h,x_initial=x,g_initial=g,
            photons=np.array([v[0] for v in pair]),gas=np.array([v[1] for v in pair]))
        for j,v in enumerate(pair):
            xx=v[0]
            moments.append(np.array([np.sum(xx*m.Eweight,axis=(1,2)),np.sum(xx*m.Eweight*m.model.bulk.mu2[None,:,None],axis=(1,2)),np.sum(xx*m.Nweight,axis=(1,2))])*AMP)
            gas_records.append(v[1].copy());times.append(m.stage_t[-2+j]);weights.append(m.stage_h[-2+j])
        write(OUT/'capture-progress.json',dict(replayed_actual_steps=index+1,target=8,stage_time=float(m.stage_t[-1])))
        return pair,mechanical
    boundary=owner.Model.boundary_ports
    def observed_boundary(m,t,x):
        value=boundary(m,t,x);ports.append(value.copy()*AMP);return value
    owner.Model.boundary_ports=observed_boundary
    owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,stages=observed),argdefs=run.__defaults__)
    m=owner.Model(128);row=m.run(128,'capture-128',8,'restored-128');assert row['passed'] and row['actual_new_steps']==8
    original=dict(np.load(OLD/'sweep-1/photons/pilot-128.npz'));got=dict(np.load(OUT/'sweep-1/photons/capture-128.npz'))
    keys=list(prior.HISTORY.values())+['material_floor_discard_scaled','restart_guide','t','moments','radial_ports','photon_history_scaled_occupation','material_history','collision_transfer',
        'accepted_angular_times','accepted_angular_luminosity','accepted_angular_quadrature_weights','actual_step_edges','energy_offset_reference','energy_offset_t']
    keys += ['restart_'+key for key in ['x','g','ledger','escape','impulse','ports','transfer']]
    differences={k:float(np.max(abs(original[k]-got[k]))) for k in keys if not np.array_equal(original[k],got[k])}
    write(OUT/'physical-replay.json',dict(classification='Counterexample candidate',passed=not differences,checked_keys=keys,differences=differences,replayed_actual_steps=8,original_stage_producer_exact=True))
    assert not differences,('Original physical replay mismatch',differences)
    seed=dict(np.load(recovery.OUT/'recovered-128.npz'));assert len(seed['times'])==len(times)==len(ports)==16
    assert np.array_equal(times,original['joint_stage_times'][16:]) and np.array_equal(weights,original['joint_stage_weights'][16:])
    mm=np.r_[seed['photon_moments'],moments];pp=np.r_[seed['radial_ports'],ports]
    cc=np.r_[seed['collision_rates'],got['joint_collision_rates_scaled'][16:]];aa=np.r_[seed['angular'],got['accepted_angular_luminosity'][16:]]
    snapshot=FunctionType(snapshots.snapshot.__code__,dict(snapshots.snapshot.__globals__,OUT=OUT))
    control=snapshot(128,m,original,15,got['restart_x'],cc,pp,aa)
    np.savez_compressed(OUT/'recovered-128.npz',times=original['joint_stage_times'],weights=original['joint_stage_weights'],photon_moments=mm,
        collision_rates=cc,angular=aa,radial_ports=pp,endpoint_occupation=got['restart_x']*m.scale*AMP)
    rows=read(recovery.OUT/'recovered-128.json')['rows']+[dict(step=k,method='original_coupled_stage_capture',original_physical_arrays_exact=True) for k in range(8,16)]
    np.savez_compressed(OUT/'accepted-128.npz',step=16,t=original['actual_step_edges'][16],x=got['restart_x'],moments=mm,collisions=cc,packets=aa,ports=pp,logs=np.array(json.dumps(rows)))
    result=dict(classification='Counterexample candidate',passed=True,actual_original_coupled_interval_reproduced=True,all_physical_array_values_exact=True,
        reused_photon_steps=8,captured_actual_steps=8,replayed_physical_steps=8,new_physical_horizon=False,snapshot=control,
        previous_conditional_failures_preserved=True,full_history_completed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
