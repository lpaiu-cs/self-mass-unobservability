"""Resolve the failed material derivative on saved states before any replay."""
from pathlib import Path
import json,resource,signal,sys,time
import numpy as np
import def_retained_metric_return as run

OUT=run.OUT;read,write,sha=run.read,run.write,run.sha


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


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    else:
        assert not (OUT/'material-diagnosis-receipt.json').exists();plan=read(OUT/'material-repair-plan.json')
        assert sha(__file__)==plan['source_sha256']
        for p,h in plan['bindings'].items():assert sha(p)==h,p
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
        def timeout(*_):raise TimeoutError('50s material diagnosis budget')
        signal.signal(signal.SIGALRM,timeout);signal.alarm(50);start=time.monotonic();cpu=time.process_time();error=None
        try:diagnose()
        except Exception as exc:error=repr(exc);raise
        finally:
            signal.alarm(0);write(OUT/'material-diagnosis-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,error=error,source_sha256=sha(__file__)))
