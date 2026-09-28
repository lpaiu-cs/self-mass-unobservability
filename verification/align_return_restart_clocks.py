"""Retain each accepted restart clock under the existing stage-time contract."""
from pathlib import Path
from types import SimpleNamespace
import inspect,os,resource,sys,time
import numpy as np
import complete_returned_period as prior

ROOT,OUT=prior.ROOT,prior.OUT
read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind


def aligned(times):
    # Same1e-18time matching used by the existing StageDriver and source reader.
    assert abs(times[0]-times[1])<1e-18,times


def check():
    rows=[];files=[Path(__file__),Path(prior.__file__)]
    for n in [64,128]:
        path=prior.OLD/f'sweep-1/photons/interval-14-{n}.npz'
        with np.load(path) as p,np.load(prior.saved(n)) as high:
            count=len(p['actual_step_edges'])-1
            assert np.array_equal(p['actual_step_edges'],high['actual_step_edges'][:count+1])
            assert np.array_equal(p['joint_stage_times'],high['joint_stage_times'][:2*count])
            rows.append(dict(clock=n,steps=count,restart_seconds=float(p['actual_step_edges'][-1]),own_raw_clock_exact=True))
        files += [path,prior.saved(n)]
    times=[r['restart_seconds'] for r in rows];aligned(times)
    try:aligned([times[0],times[0]+1e-17])
    except AssertionError:pass
    else:raise AssertionError('Distinct stage times accepted')
    write(ROOT/'clock-alignment.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        original_bitwise_equality=times[0]==times[1],difference_seconds=times[1]-times[0],existing_stage_alignment_threshold=1e-18,
        shifted_clock_rejected=True,raw_times_changed=False,
        bindings={str(p):sha(p) for p in files},final_charge_conclusion='unadjudicated'))


def link(src,dst):
    if Path(dst).exists():assert sha(src)==sha(dst),(src,dst)
    else:os.link(src,dst)


def prepare():
    proof=read(ROOT/'clock-alignment.json');assert proof['passed']
    for p,h in proof['bindings'].items():assert sha(p)==h,p
    s=inspect.getsource(prior.prepare)
    s=prior.base.replace(s,'assert times[0]==times[1]','aligned(times)')
    safe_os=SimpleNamespace(**dict(vars(os),link=link))
    ns=dict(prior.prepare.__globals__,aligned=aligned,os=safe_os,seed=bind(prior.seed,os=safe_os))
    exec(compile(s,__file__,'exec'),ns);ns['prepare']()
    p=read(OUT/'plan.json');p.update(restart_clocks=proof['rows'],restart_raw_times_changed=False,
        clock_alignment='Original StageDriver1e-18absolute matching; retain each actual restart clock. Physical/time/overlap gates unchanged.')
    p['bindings'].update({str(Path(__file__)):sha(__file__),str(ROOT/'clock-alignment.json'):sha(ROOT/'clock-alignment.json')})
    write(OUT/'plan.json',p)


if __name__=='__main__':
    action=sys.argv[1];assert action in prior.CAPS;receipt=ROOT/('clock-check-receipt.json' if action=='check' else f'{action}-receipt.json');assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3));prior.actual.prior.joint.previous.original.inf.incident.native.deadline(prior.CAPS[action]);start=time.monotonic();error=None
    try:
        if action not in ['check','prepare']:
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in ['check','prepare']:globals()[action]()
        elif action in ['coarse','fine']:prior.evolve(64 if action=='coarse' else 128)
        else:getattr(prior,action)()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
