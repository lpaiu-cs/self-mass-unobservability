"""Measured-budget execution of the frozen pulse model; reuse the pilot prefix."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import inspect
import json
import shutil
import time
import numpy as np
import def_resolved_scalar_pulse as s

OUT=s.OUT/'production'
source=inspect.getsource(s.run_path)
source=source.replace("folder=OUT/f'{label}-{steps}'", "folder=PRODUCTION/f'{label}-{steps}'")
source=source.replace('    history=[];budgets=[];records=[]', '''    history=[];budgets=[];records=[];start_step=0
    if label=='driven' and steps==48:
        prefix=OUT/'driven-48';saved=json.loads((prefix/'result.json').read_text())
        start_step=saved['completed_steps'];history=saved['history'];budgets=saved['budgets']
        cp=np.load(prefix/f'step-{start_step:03d}.npz')
        delta=cp['delta'].copy();z={k:cp[k].copy() for k in cp.files if k not in ['delta','residual']}
        star.previous=z;star.field_guess=z['psi'].copy(),z['Pi'].copy(),z['Phi'].copy()
        star.amplitude=ld(str(history[-1]['boundary']))
        # Reconstruct the exact original boundary expression, not JSON rounding.
        star.amplitude=peak*np.sin(ld(str(np.pi))*ld(start_step)*h/(2*tc))**8
        for q in prefix.glob('step-*.npz'):shutil.copy2(q,folder/q.name)
        shutil.copy2(prefix/'iterations.jsonl',folder/'iterations.jsonl')
        y=star.base+delta;logA=z['logA']
        star.material_cache.update({e.material_key(row):raw for row,raw in zip(zip(y[:,0]-3*logA,y[:,1]-logA,y[:,5:]),z['raw'])})''')
source=source.replace('range(1,(limit or steps)+1)', 'range(start_step+1,(limit or steps)+1)')
source=source.replace("            print('PULSE',label,steps,j,history[-1]['native_norm'],flush=True)",
    "            e.write(PRODUCTION/'progress.json',dict(path=label,steps=steps,completed=j,native_norm=history[-1]['native_norm'],path_seconds=time.monotonic()-start))")
namespace=dict(vars(s),PRODUCTION=OUT,shutil=shutil)
exec(compile(source,__file__,'exec'),namespace);run_path=namespace['run_path']


def prepare():
    p=s.bindings();pilot=json.loads((s.OUT/'pilot.json').read_text());assert pilot['passed'] and not OUT.exists()
    files=[Path(__file__),s.OUT/'plan.json',s.OUT/'pilot.json',s.OUT/'initial-manifest.json',s.OUT/'initial.npz']+list((s.OUT/'driven-48').glob('*'))
    OUT.mkdir()
    s.e.write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',
        bindings={q.relative_to(s.ROOT).as_posix():s.e.digest(q) for q in files if q.is_file()},
        expected_wall_seconds=[900,1800],hard_timeout_seconds=1800,
        budget_revision='Before production launch: 4 native pilot steps took 5.34s. The original provisional 900-second cap cannot cover 1008 stages at that measured rate. This separate execution plan budgets one unchanged nine-path experiment up to 1800s, with early coarse-pair stopping and no extra paths or tolerance changes.',
        pilot_seconds_per_step=pilot['seconds']/4,workers=4,blas_threads=1,
        early_stop='Stop on any native/energy failure; after the three 48-step paths, stop if the preregistered response floor is not met. Before refinement, stop if measured remaining time would exceed 1800 seconds.',
        comparison='At common final time: R=Jordan logT(driven)-Jordan logT(undriven); C=Jordan logT(driven)-Jordan logT(beta=0). RMS weights are the original conserved coarse baryons. Fixed-grid orders use norms of 48-96 and 96-192 vector differences. Empirical time error indicator is twice the latter difference, not a rigorous bound.',
        unchanged_gates=p['gates'],same_initial_state=True,reused_pilot_steps=4,physical_orbital_match=False))


def vectors(steps):
    data={label:np.load(OUT/f'{label}-{steps}'/f'step-{steps:03d}.npz') for label in ['driven','undriven','decoupled']}
    theta={label:cp['delta'][:,1]-cp['logA'] for label,cp in data.items()}
    return theta['driven']-theta['undriven'],theta['driven']-theta['decoupled']


def summarize(results):
    p=s.bindings();B=np.load(s.OUT/'initial.npz')['baryons'];weight=B/B.sum()
    norm=lambda v:float(np.sqrt(np.sum(weight*v*v)))
    vec=[vectors(k) for k in p['time_steps']];reads=[]
    for index,name in enumerate(['driven_minus_undriven','material_coupling_contrast']):
        a,b,c=[v[index] for v in vec];d1,d2=norm(a-b),norm(b-c);signal=norm(c)
        order=float(np.log2(d1/d2)) if min(d1,d2)>0 else None
        reads.append(dict(name=name,amplitude=signal,differences=[d1,d2],order=order,
            relative_time_difference=d2/signal,signal_over_empirical_time_error=signal/(2*d2)))
    gates=p['gates']
    passed=all(r['order'] is not None and r['order']>=gates['minimum_time_order']
        and r['relative_time_difference']<=gates['maximum_relative_time_difference']
        and r['signal_over_empirical_time_error']>=gates['minimum_control_separation_over_time_error']
        and r['amplitude']>=gates['minimum_Jordan_temperature_rms'] for r in reads)
    return dict(classification='Counterexample candidate',passed=bool(passed),readouts=reads,paths=results,
        comparison_scope=p['comparator'],temperature_is_not_scalar_charge=True,
        fixed_grid_only=True,physical_EOS_certified=False,observational_closure=False)


def run():
    plan=json.loads((OUT/'execution-plan.json').read_text())
    for rel,h in plan['bindings'].items():assert s.e.digest(s.ROOT/rel)==h,rel
    assert not (OUT/'result.json').exists() and not (OUT/'failure.json').exists()
    started=time.monotonic();results=[]
    try:
        with ProcessPoolExecutor(max_workers=4,initializer=s.old.imported.initial.original.worker_init) as pool:
            for steps in [48,96,192]:
                for label in ['driven','decoupled','undriven']:
                    result=run_path(pool,label,steps)
                    results.append(dict(label=label,steps=steps,seconds=result['seconds'],native_calls=result['native_calls']))
                    print('COMPLETED',label,steps,result['seconds'],flush=True)
                if steps==48:
                    B=np.load(s.OUT/'initial.npz')['baryons'];v=vectors(48)
                    signals=[float(np.sqrt(np.sum(B*x*x)/B.sum())) for x in v]
                    elapsed=time.monotonic()-started
                    s.e.write(OUT/'coarse-decision.json',dict(classification='Counterexample candidate',signals=signals,seconds=elapsed,
                        projected_total_seconds=elapsed*7,scope='Coarse signal and measured cost only; no time-converged response yet.'))
                    assert min(signals)>=s.bindings()['gates']['minimum_Jordan_temperature_rms'],('Response floor',signals)
                    assert elapsed*7<plan['hard_timeout_seconds'],('Remaining cost exceeds cap',elapsed*7)
        result=summarize(results);result['seconds']=time.monotonic()-started;s.e.write(OUT/'result.json',result)
    except Exception as error:
        s.e.write(OUT/'failure.json',dict(classification='Counterexample candidate',error=repr(error),seconds=time.monotonic()-started,completed_paths=results));raise
    s.e.write(OUT/'manifest.json',dict(sha256={q.relative_to(s.ROOT).as_posix():s.e.digest(q) for q in OUT.rglob('*') if q.is_file() and q.name not in ['manifest.json','run.log']}))
    print('FINAL',json.dumps(result),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','run'])
    globals()[parser.parse_args().action]()
