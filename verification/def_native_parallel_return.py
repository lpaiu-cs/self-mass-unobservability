"""Bounded execution of the three independent, unchanged response paths."""
from pathlib import Path
import json,os,resource,signal,subprocess,sys,time
import numpy as np
import verify_native_updated_gr_return as v

OUT=v.run.RESPONSE;write=v.write;sha=v.sha
PATHS=[(64,128),(128,128),(128,64)]


def worker(steps,reference,label,limit,restart):
    # Three processes share a6GiB virtual-address cap, one BLAS thread each.
    cap=2*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap))
    v.configure();start=time.monotonic();cpu=time.process_time()
    row=v.Response(reference).run(steps,label,limit=limit,restart=restart)
    row.update(worker_wall_seconds=time.monotonic()-start,worker_CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    write(OUT/f'{label}.json',row)


def dispatch(specs,cap):
    start=time.monotonic();children=[]
    try:
        for n,ref,label,limit,restart in specs:
            log=(OUT/f'{label}.log').open('w')
            args=[sys.executable,__file__,'worker',str(n),str(ref),label,str(limit),str(restart)]
            child=subprocess.Popen(args,stdout=log,stderr=subprocess.STDOUT,start_new_session=True);log.close();children.append((child,label))
        while any(p.poll() is None for p,_ in children):
            if time.monotonic()-start>cap:raise TimeoutError('Registered parallel wall cap')
            for p,label in children:
                if p.poll() not in [None,0]:raise RuntimeError(('Response worker',label,p.returncode))
            time.sleep(.25)
        assert all(p.returncode==0 for p,_ in children)
        return [json.loads((OUT/f'{label}.json').read_text()) for _,label in children],time.monotonic()-start
    finally:
        for p,_ in children:
            if p.poll() is None:os.killpg(p.pid,signal.SIGTERM)
        for p,_ in children:p.wait()


def pilot():
    assert not (OUT/'parallel-pilot.json').exists();previous=json.loads((OUT/'pilot.json').read_text());assert not previous['eligible']
    write(OUT/'parallel-plan.json',dict(classification='Counterexample candidate',
        reason='Serial measured forecast937s exceeds650s. Profiling confirms actual Krylov/preconditioner arithmetic dominates; repeated table reads are no longer the main cost. These three prescribed response paths are independent.',
        decision='Measure concurrent four-step continuations before dispatch. Keep650s production wall cap and all equations/gates; explicitly change resource allocation from one to three CPU processes, one thread each,2GiB virtual memory each. No GPU port, new path, finer grid or relaxed residual.',
        reuse='Reuse all bank/GR/background data and saved physical response prefixes. The profile continuation is reused for the fine path. No completed response segment is repeated.',
        budget=dict(pilot_wall_seconds=65,production_wall_seconds=650,CPU_processes=3,threads_each=1,total_virtual_memory_GiB=6),
        forecast='For each concurrently measured path separately:17 operator points, remaining step cost and15s setup/export;2x the maximum must fit650s. This tests actual contention rather than dividing serial time by3. Later iterations remain extrapolated.',
        stop='No dispatch on failed concurrent pilot, memory limit, forecast or physics gate. Stop all workers if one fails or650s expires. No automatic increase.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(v.__file__),Path(v.run.__file__),OUT/'pilot.json',OUT/'profile.json']}))
    specs=[(64,128,'parallel-pilot-64-128',8,'pilot-64'),(128,128,'parallel-pilot-128-128',16,'profile-128'),(128,64,'parallel-pilot-128-64',4,None)]
    rows,elapsed=dispatch(specs,65);estimates=[]
    for r in rows:
        point=r['operator_point_seconds']/r['operator_points']
        estimates.append(point*17+r['stepping_seconds']/r['new_steps']*(r['steps']-r['completed_steps'])+15)
    result=dict(classification='Counterexample candidate',rows=rows,forecast_each_seconds=estimates,upper_wall_seconds=2*max(estimates),
        eligible=all(r['passed'] for r in rows) and 2*max(estimates)<650,seconds=elapsed)
    write(OUT/'parallel-pilot.json',result);print(json.dumps(result),flush=True)
    if result['eligible']:write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=650,
        CPU_processes=3,total_virtual_memory_GiB=6,paths=[[n,ref,f'parallel-pilot-{n}-{ref}'] for n,ref in PATHS],
        bindings={str(p):sha(p) for p in [Path(__file__),Path(v.__file__),Path(v.run.__file__),OUT/'parallel-plan.json',OUT/'parallel-pilot.json',v.LAPSE/'result.json',v.run.COLL/'bank-result.json']}))


def production():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'execution-plan.json').read_text());assert plan['eligible']
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    specs=[(n,ref,f'steps-{n}-reference-{ref}',None,restart) for n,ref,restart in plan['paths']]
    try:rows,elapsed=dispatch(specs,650)
    except Exception as exc:write(OUT/'parallel-failure.json',dict(error=repr(exc)));raise
    def history(n,ref):
        d=np.load(OUT/f'steps-{n}-reference-{ref}.npz');t=np.linspace(0,v.run.flow.old.END,17)
        ids=[int(np.argmin(abs(d['t']-s))) for s in t];assert np.max(abs(d['t'][ids]-t))<1e-18
        return d['moments'][ids][:,[0,1,2,3,5,6]]
    fine=history(128,128);norm=np.maximum(np.max(np.sum(abs(fine),axis=2),axis=0),1.)
    comparisons={key:(np.max(np.sum(abs(other-fine),axis=2),axis=0)/norm).tolist()
        for key,other in [('time',history(64,128)),('background_time',history(128,64))]}
    result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(max(x) for x in comparisons.values())<.02,
        paths=rows,comparisons=comparisons,comparison_order=['photon_energy','material_energy','neutral_count','momentum_impulse','photon_radial_pressure','material_pressure'],
        seconds=elapsed,CPU_seconds=sum(r['worker_CPU_seconds'] for r in rows),sum_peak_RSS_bytes=sum(r['peak_RSS_bytes'] for r in rows),
        actual_monolithic_radiation_thermal_H_response_evolved=True,additional_material_motion_evolved=False,full_GR_feedback=False,final_charge_solved=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
