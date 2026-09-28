"""Separate the saved native energy discrepancy before modifying evolution."""
from pathlib import Path
import gc,json,os,resource,sys,time
import numpy as np
import finish_returned_material_accuracy as producer
import resolve_returned_native_time as previous

OUT=Path('native-mass-time258-work');OLD=producer.OUT
read,write,sha,bind=producer.read,producer.write,producer.sha,producer.bind
LD=np.longdouble;AMP=producer.base.actual.AMP
resume=len(sys.argv)>1 and sys.argv[1]=='resume'
if resume:
    assert read(OUT/'receipt.json')['error']=='AssertionError(1037.9218287335589)'
    assert sha(OUT/'initial-budget-producer.py')==read(OUT/'receipt.json')['source_sha256']
else:
    assert not OUT.exists();OUT.mkdir()
start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,)*2)
producer.joint.previous.original.inf.incident.native.deadline(1800 if resume else 900)
write(OUT/('resume-start.json' if resume else 'start.json'),dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],
    boot_id=Path('/proc/sys/kernel/random/boot_id').read_text().strip(),started_unix=time.time(),source_sha256=sha(__file__)))
try:
    files=list((OLD/'sweep-0').rglob('*.npz'))
    files += [p for part in ['metric','gr'] for p in (OLD/part).iterdir() if p.is_file()]
    files += [OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True)
        if resume:assert sha(p)==sha(dst),p
        else:os.link(p,dst)
    for part in ['sweep-1/photons','sweep-1/material']:(OUT/part).mkdir(parents=True)
    if not resume:write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Locate the native thermal-rate contribution to the failed same-solution homogeneous mass clock comparison.',
        decision='Separate coarse quadrature/base-history from evolved-state differences using existing189dense projection. Repair the identified owner in the actual joint equation; do not add a diagnostic correction to charge.',
        evidence='257energy ledger:2.3066239671e-19erg native difference, minus5.535682e-21erg collision difference; floor negligible. Cells1..3 account for99.2percent of the net thermal difference.',
        gates=dict(saved_native=1e-12,dense_stage=1e-12,decomposition=1e-8),
        physical_steps=0,new_clock_paths=0,seconds_cap=900,virtual_GiB=8,CPU_affinity=2,
        forecast='Prior native replay32calls plus initialization completed in tens of seconds. Measure three current calls, then allow the462saved fine stages only if a4x maximum-rate forecast is below15minutes. No physical re-integration.',
        original_mass_time_gate=.02,original_failed_verdict_preserved=True,
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in [Path(__file__),Path(previous.__file__),OLD/'result.json',Path('native-material-exterior257-work/full/audit.json')]}))
    if resume:
        write(OUT/'resume-plan.json',dict(classification='Conjectural',
            reason='Measured calls0.484..0.518s imply about4minutes for462calls. The conservative4x forecast1037.92s exceeded the initial15minute planning gate before projection started. Preserve that stop and continue the same finite462calls with30minutes; no accuracy or physical/grid expansion.',
            max_seconds=1800,virtual_GiB=8,new_physical_steps=0,new_clock_paths=0,
            bindings={str(p):sha(p) for p in [Path(__file__),OUT/'initial-budget-producer.py',OUT/'plan.json',OUT/'receipt.json',OUT/'measured-forecast.json']},
            old_execution_source_alias={'source':'.phase258-native-projection.py','frozen':str(OUT/'initial-budget-producer.py')}))
    Model=bind(producer.base.initialize,OUT=OUT)();m=Model(128)
    keys=['actual_step_edges','joint_stage_conserved_scaled','joint_native_rates_scaled','joint_collision_rates_scaled','joint_stage_weights','joint_stage_times','material_floor_discard_scaled','material_history','material_energy_units']
    p=[]
    for n in [64,128]:
        with np.load(OLD/f'sweep-1/photons/return-{n}.npz') as z:p.append({k:z[k] for k in keys})
    a,b=p;controls=[]
    for z in p:
        norm=[]
        for i,t in enumerate(z['joint_stage_times']):
            expected=previous.gas(m,z['joint_stage_conserved_scaled'][i]);dense=previous.dense(m,z,float(t))
            norm.append(np.sum(abs(dense-expected)*m.units,axis=0)/np.maximum(np.sum(abs(expected)*m.units,axis=0),LD('1e-290')))
        value=float(np.max(norm));assert value<1e-12,value;controls.append(value)
    forecast=[]
    for i in [0,len(b['joint_stage_times'])//2,len(b['joint_stage_times'])-1]:
        t=float(b['joint_stage_times'][i]);g=previous.gas(m,b['joint_stage_conserved_scaled'][i]);begin=time.monotonic()
        value=m.native(t,g)*m.units;seconds=time.monotonic()-begin
        rel=(np.sum(abs(value-b['joint_native_rates_scaled'][i]),axis=0)/np.maximum(np.sum(abs(b['joint_native_rates_scaled'][i]),axis=0),LD('1e-290'))).astype(float).tolist()
        assert max(rel)<1e-12,rel;forecast.append(dict(stage=i,seconds=seconds,relative=rel))
    upper=4*max(v['seconds'] for v in forecast)*len(b['joint_stage_times'])+time.monotonic()-start+60
    cap=1800 if resume else 900
    write(OUT/('resume-forecast.json' if resume else 'measured-forecast.json'),dict(rows=forecast,upper_seconds=upper,passed=upper<cap));assert upper<cap,upper
    values=[]
    for i,t in enumerate(b['joint_stage_times']):
        values.append(m.native(float(t),previous.dense(m,a,float(t)))*m.units)
        if i%64==0:write(OUT/'progress.json',dict(completed=i+1,total=len(b['joint_stage_times']),seconds=time.monotonic()-start))
    weighted=[np.sum(z['joint_stage_weights'][:,None,None].astype(LD)*z['joint_native_rates_scaled'],axis=0,dtype=LD) for z in p]
    projected=np.sum(b['joint_stage_weights'][:,None,None].astype(LD)*np.array(values),axis=0,dtype=LD)
    quadrature=(weighted[0]-projected)*AMP;path=(projected-weighted[1])*AMP
    native=(weighted[0]-weighted[1])*AMP;den=max(np.sum(abs(native[:,0])),LD('1e-290'))
    defect=float(np.sum(abs(native[:,0]-quadrature[:,0]-path[:,0]))/den);assert defect<1e-8,defect
    result=dict(classification='Counterexample candidate',passed=True,dense_controls=controls,
        native_thermal_difference=float(np.sum(native[:,0])),quadrature_and_base_difference=float(np.sum(quadrature[:,0])),
        evolved_state_difference=float(np.sum(path[:,0])),
        quadrature_L1_fraction=float(np.sum(abs(quadrature[:,0]))/den),state_L1_fraction=float(np.sum(abs(path[:,0]))/den),
        decomposition_relative=defect,physical_steps=0,new_clock_paths=0,
        scope='Exact saved-history split using coarse dense gas on fine high-anchor history; it does not establish a unique root cause or correct the physical solution.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'projection.npz',quadrature_and_base=quadrature,evolved_state=path,native=native,
        coarse_dense_native_on_fine=np.array(values),fine_times=b['joint_stage_times'],fine_weights=b['joint_stage_weights'])
    write(OUT/'result.json',result);print(json.dumps(result),flush=True)
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/('resume-receipt.json' if resume else 'receipt.json'),dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
