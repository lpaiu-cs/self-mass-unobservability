"""Resolve the failed material derivative on saved states before any replay."""
from pathlib import Path
from types import FunctionType
import json,resource,shutil,signal,sys,time
import numpy as np
import def_retained_metric_return as run

OUT=run.OUT;read,write,sha=run.read,run.write,run.sha
ORIGINAL=run.MATERIAL;MATERIAL=OUT/'material-resolved'


def prepare():
    failure=read(run.MATERIAL/'steps-128-reference-128.json');assert not failure['passed']
    assert not (OUT/'material-repair-plan.json').exists()
    write(OUT/'material-repair-plan.json',dict(classification='Counterexample candidate',
        failure='Fine material path completed, but endpoint momentum directional comparison1.0798629percent exceeds original0.2percent. Preserve this path as failed.',
        claim='Determine whether amplified native flux/recovery arithmetic, rather than physical donor changes, causes the failed direction; do not accept a failed derivative or merely loosen its tolerance.',
        method='Reuse saved fine material states at first, middle and last canonical times. Compare the unchanged Richardson owner at0.5,1,2,4,8,16 probe scales. Locate largest derivative changes. No new trajectory in diagnosis.',
        allocation='Move250s from unspent photon production reserve:50s for diagnosis and200s to the material original350s allowance. Total2050s unchanged. Any corrected replay requires measured pilot and2x forecast within the remaining material550s. No new clock, background, EOS bank or resolution.',
        diagnosis_cap_seconds=50,material_total_cap_seconds=550,total_original_seconds=2050,
        stop='Stop if no common probe regime meets original0.2percent across saved states, or corrected actual-amplitude donor/support gates fail. Do not reduce physical amplitude. Full production is not authorized by this diagnostic alone.',
        source_sha256=sha(__file__),bindings={str(p):sha(p) for p in [Path(run.__file__),OUT/'response-execution-plan.json',
            run.MATERIAL/'steps-128-reference-128.npz',run.MATERIAL/'steps-128-reference-128.json',OUT/'material_production-receipt.json']}))


def diagnose():
    run.initialize();m=run.Material(128,128);d=np.load(run.MATERIAL/'steps-128-reference-128.npz');rows=[];saved={}
    for k in [1,8,16]:
        t=m.t[k];j=int(np.argmin(abs(d['t']-t)));z=d['history_scaled'][j];rates={}
        for probe in [.5,1.,2.,4.,8.,16.]:
            mark=time.monotonic();r,ledger,dt=m.rhs(t,z,probe);rates[probe]=r;saved[f'rate-{k}-{probe}']=r
            rows.append(dict(k=k,probe=probe,seconds=time.monotonic()-mark,
                L1=np.sum(abs(r),axis=1).tolist(),physical_branch=m.physical_branch_ratio))
        errors={str(p):(np.sum(abs(rates[p]-rates[p/2]),axis=1)/np.maximum(np.sum(abs(rates[p/2]),axis=1),1.)).tolist() for p in [1.,2.,4.,8.,16.]}
        diff=abs(rates[1]-rates[.5]);ids=np.argsort(diff[1])[-8:][::-1]
        rows.append(dict(k=k,adjacent_probe_comparisons=errors,largest_momentum_cells=ids.tolist(),
            momentum_difference=diff[1,ids].tolist(),directional_momentum=rates[1][1,ids].tolist()))
    np.savez_compressed(OUT/'material-probe-rates.npz',**saved)
    write(OUT/'material-probe-diagnosis.json',dict(classification='Counterexample candidate',rows=rows,
        actual_trajectory_replayed=False,physical_branch_ratio=m.physical_branch_ratio,final_charge_solved=False))
    print(json.dumps(rows),flush=True)


def initialize():
    run.initialize();Base=run.Material
    class Material(Base):
        def rhs(self,t,z,probe=1.):return super().rhs(t,z,probe*(8 if self.steps==128 else 1))
        run=FunctionType(Base.run.__code__,dict(Base.run.__globals__,OUT=MATERIAL),argdefs=Base.run.__defaults__)
    run.Material=Material;run.MATERIAL=MATERIAL


def register():
    assert not MATERIAL.exists();d=read(OUT/'material-probe-diagnosis.json')
    rows=[r for r in d['rows'] if 'adjacent_probe_comparisons' in r]
    assert len(rows)==3 and all(max(r['adjacent_probe_comparisons'][p])<.002 for r in rows for p in ['8.0','16.0'])
    assert d['physical_branch_ratio']<.01
    p=read(OUT/'material-repair-plan.json');assert sha(OUT/'diagnostic-material-producer.py')==p['source_sha256']
    MATERIAL.mkdir();coarse=read(ORIGINAL/'steps-64-reference-128.json');assert coarse['passed']
    for suffix in ['json','npz']:shutil.copyfile(ORIGINAL/f'steps-64-reference-128.{suffix}',MATERIAL/f'steps-64-reference-128.{suffix}')
    remaining=550-read(OUT/'material_production-receipt.json')['seconds']
    write(OUT/'material-repair-execution.json',dict(classification='Counterexample candidate',
        diagnosis='Probe enlargement gives decreasing deep-cell momentum differences: fine endpoint1/0.5=1.07986percent,8/4=0.0545688percent,16/8=0.0204734percent. Physical donor change stays negligible. Evidence supports underresolved amplified flux/recovery differences, not a changed physical branch.',
        repair='Use the existing8x amplified Richardson probe on fine128 only, with4/8 and8/16 checks. The actual physical amplitude, geometry, equation, clocks, support, conservation and all0.2percent gates stay fixed. Keep accepted coarse64 unchanged.',
        decision='Recompute only the failed fine material history from a fresh corrected two-step prefix; do not merely replace its endpoint verdict or replay any photons. Then compare its actual pressure/GR result against the retained coarse path.',
        inherited_coarse_sha256=sha(ORIGINAL/'steps-64-reference-128.npz'),
        production_remaining_seconds=remaining,pilot_remaining_seconds=50-read(OUT/'material_pilot-receipt.json')['seconds'],
        original_total_seconds=2050,material_total_seconds=550,source_sha256=sha(__file__),
        bindings={str(f):sha(f) for f in [Path(run.__file__),OUT/'material-repair-plan.json',OUT/'material-probe-diagnosis.json',
            OUT/'diagnostic-material-producer.py',ORIGINAL/'steps-128-reference-128.json',ORIGINAL/'steps-64-reference-128.json']}))


def worker(pilot):
    initialize();start=time.monotonic();m=run.Material(128,128);label='pilot-128' if pilot else 'steps-128-reference-128'
    row=m.run(128,label,2 if pilot else None,restart=None if pilot else 'pilot-128')
    row.update(worker_wall_seconds=time.monotonic()-start,physical_branch_ratio=m.physical_branch_ratio,
        maximum_owner_error=max(p['owner_error'] for p in m.cache.values()),probe_multiplier=8)
    if not pilot:
        d=np.load(MATERIAL/f'{label}.npz');z=d['delta_scaled'];t=float(d['time'])
        a=m.rhs(t,z,1.)[0];b=m.rhs(t,z,2.)[0]
        row['reverse_probe_8_16']=(np.sum(abs(a-b),axis=1)/np.maximum(np.sum(abs(b),axis=1),1.)).tolist()
        row['passed']=row['passed'] and max(row['reverse_probe_8_16'])<.002
    row['passed']=row['passed'] and row['physical_branch_ratio']<.01 and row['maximum_owner_error']<1e-8
    write(MATERIAL/f'{label}.json',row);assert row['passed'],row
    return row


def pilot():
    row=worker(True);old=read(ORIGINAL/'steps-128-reference-128.json');plan=read(OUT/'material-repair-execution.json')
    forecast=old['raw_owner_calls']*row['seconds']/row['raw_owner_calls']+row['worker_wall_seconds']-row['seconds']+5
    result=dict(classification='Counterexample candidate',row=row,upper_remaining_seconds=2*forecast,
        eligible=row['passed'] and 2*forecast<plan['production_remaining_seconds'])
    write(MATERIAL/'pilot.json',result);print(json.dumps(result),flush=True);assert result['eligible']


def production():
    assert read(MATERIAL/'pilot.json')['eligible'];start=time.monotonic();row=worker(False)
    original=read(ORIGINAL/'steps-64-reference-128.json');assert original['passed']
    write(MATERIAL/'production.json',dict(classification='Counterexample candidate',passed=True,rows=[original,row],
        seconds=time.monotonic()-start,accepted_coarse_reused=True,failed_fine_preserved=True,corrected_fine_recomputed=True))


if __name__=='__main__':
    action=sys.argv[1]
    if action in ['prepare','register']:globals()[action]()
    else:
        plan=read(OUT/('material-repair-plan.json' if action=='diagnose' else 'material-repair-execution.json'))
        name='material-diagnosis' if action=='diagnose' else 'material-repair-'+action
        assert not (OUT/f'{name}-receipt.json').exists();assert sha(__file__)==plan['source_sha256']
        for p,h in plan['bindings'].items():assert sha(p)==h,p
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
        cap=50 if action=='diagnose' else plan[('pilot' if action=='pilot' else 'production')+'_remaining_seconds']
        def timeout(*_):raise TimeoutError('Registered material repair budget')
        signal.signal(signal.SIGALRM,timeout);signal.setitimer(signal.ITIMER_REAL,cap);start=time.monotonic();cpu=time.process_time();error=None
        try:globals()[action]()
        except Exception as exc:error=repr(exc);raise
        finally:
            signal.setitimer(signal.ITIMER_REAL,0);write(OUT/f'{name}-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,error=error,source_sha256=sha(__file__)))
