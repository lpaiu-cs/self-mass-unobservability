"""Apply the whole-history GR consumer to an accepted cross-background prefix."""
from pathlib import Path
import json,os,resource,shutil,sys,time
import numpy as np
import read_complete_radau_history as r

OUT=Path('native-complete-radau224-check-work');ORIGINAL=r.INPUT
r.OUT=OUT;r.INPUT=OUT/'input';old_saved=r.saved;r.saved=lambda n:OUT/'input'/f'material-{n}.npz'
CAPS=dict(prepare=180,endpoints=900,source=1800,fields=1800)
action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3));r.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action])
try:
    if action=='prepare':
        assert not OUT.exists();OUT.mkdir();files=[Path(__file__),Path(r.__file__)];horizons=[]
        for folder in ['sweep-0','sweep-1/photons','sweep-1/material','gr','input']:(OUT/folder).mkdir(parents=True)
        for src in list((ORIGINAL/'sweep-0').rglob('*.npz'))+[ORIGINAL/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
            dst=OUT/src.relative_to(ORIGINAL);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
        for n,count in [(64,15),(128,29)]:
            gate=ORIGINAL/f'clock-{n}/snapshot-03.json';row=r.read(gate);assert row['passed'] and row['steps']==count
            src=ORIGINAL/f'accepted-{n}.npz';dst=OUT/'input'/f'accepted-{n}.npz'
            before=r.sha(src);shutil.copyfile(src,dst);assert before==r.sha(src)==r.sha(dst)
            q=dict(np.load(dst));p=dict(np.load(old_saved(n)));assert int(q['step'])>=count
            p['actual_step_edges']=p['actual_step_edges'][:count+1]
            for k in list(p):
                if k.startswith('joint_') and p[k].ndim:p[k]=p[k][:2*count]
            np.savez_compressed(r.saved(n),**p);horizons.append(float(p['actual_step_edges'][-1]))
            np.savez_compressed(OUT/'input'/f'recovered-{n}.npz',times=p['joint_stage_times'],weights=p['joint_stage_weights'],
                photon_moments=q['moments'][:2*count],radial_ports=q['ports'][:2*count])
            shutil.copyfile(gate,OUT/'input'/f'snapshot-{n}.json');files += [old_saved(n),r.saved(n),dst,OUT/'input'/f'recovered-{n}.npz',OUT/'input'/f'snapshot-{n}.json']
        assert abs(horizons[0]-horizons[1])<1e-18
        files += [r.prior.OUT/n for n in ['expanded-source.py','expanded-fields.py']]
        files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
        r.write(OUT/'plan.json',dict(classification='Conjectural',claim='Apply the generalized GR consumer to the accepted same-solutionT/8prefix crossing an actual background interval, then use the same consumer on the completed common history.',
            decision='Require source identities, native-pressure map, dense polynomial/live driver and independent GR controls before dispatching the full history. Preserve original source/field2percent time failures; this regression cannot admit self-GR or final charge.',
            original_return_horizon_seconds=horizons[0],new_physical_steps=0,reused_actual_steps=[15,29],budgets=CAPS,CPU_threads=1,virtual_GiB=12,
            forecast='212source6steps45.81s.44storedsteps plus setup estimated2..6minutes; allow15minutes endpoints,30minutes source,30minutes fields. No background/fluid/photon reintegration.223and218continue unchanged.',
            stop='Any original numerical gate or cap. No physical-grid extension or tolerance change.',bindings={str(p):r.sha(p) for p in dict.fromkeys(files)},full_goal_complete=False,final_charge_conclusion='unadjudicated'))
        r.write(OUT/'regression.json',r.check());r.write(OUT/'symbolic.json',r.prior.symbolic())
    else:
        for p,h in r.read(OUT/'plan.json')['bindings'].items():assert r.sha(p)==h,p
        getattr(r,action)()
except BaseException as exc:error=repr(exc);raise
finally:
    if OUT.exists():r.write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=r.sha(__file__),consumer_sha256=r.sha(r.__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
