"""Counterexample candidate: reuse the accepted15/16 full-input history.

Capture the actual failed linear solve before choosing a bounded repair. Never
restart the accepted prefix or change the physical equation and its gates.
"""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import complete_full_incident_horizon as prior

OUT=Path('native-full-finish190-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=25,probe=135)


def prepare():
    assert not OUT.exists();OUT.mkdir();failure=read(OLD/'failure.json')
    assert failure['last_comparison']['shared_interval']==15 and failure['last_comparison']['passed']
    assert 'Four-moment linear residual' in failure['error']
    assert read(OLD/'run-receipt.json')['error'] is not None
    files=[];reused={}
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    sources=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    sources += [OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    for src in sources:
        dst=OUT/src.relative_to(OLD);os.link(src,dst);reused[str(dst)]=sha(src);files.append(src)
    files += [OLD/name for name in ['plan.json','failure.json','run-receipt.json','comparison-15.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Resolve the terminal linear-solver failure of the original full-input coupled evolution without rerunning its accepted15/16history.',
        decision='Observe actual restarted-GMRES convergence, true residual and physical channels at the failed step before selecting a solver repair. Original equations, stage/Newton limits and acceptance gates remain unchanged.',
        scope='Replay only the first coarse macro after accepted60/64, stopping on the original error. No fine or full-period rerun is dispatched by this probe.',
        forecast='Original terminal attempt took91.63s after interval15. One reconstruction plus linear trace capture expected75..125s, hard cap135s; added I/O is unmeasured.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        stop='Any mismatch of prefix, owner, original failure, bindings or wall cap. Preserve the failed185producer and original2percent/1e-12/1e-13gates.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused))


def probe():
    prior.prior.OUT=OUT;prior.prior.initialize();owner=prior.owner;joint=owner.joint
    calls=[];solves=[];last={}
    def traced_gmres(op,rhs,**kw):
        history=[];callback=kw['callback'];start=time.monotonic()
        def record(v):history.append(float(v));callback(v)
        answer,info=joint.gmres(op,rhs,**dict(kw,callback=record))
        residual=np.asarray(rhs,np.longdouble)-op.matvec(answer.astype(np.longdouble))
        calls.append(dict(info=int(info),seconds=time.monotonic()-start,history=history,
            RHS_norm=float(np.linalg.norm(rhs)),solution_norm=float(np.linalg.norm(answer)),
            true_relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290))))
        last.update(solution=answer,residual=residual)
        write(OUT/'linear-calls.json',dict(classification='Counterexample candidate',calls=calls))
        return answer,info
    original=FunctionType(joint.solve.__code__,dict(joint.solve.__globals__,gmres=traced_gmres))
    def solve(m,op,P,rhs,guess):
        begin=len(calls)
        try:return original(m,op,P,rhs,guess)
        except BaseException as exc:
            np.savez_compressed(OUT/'failed-linear-system.npz',rhs=rhs,guess=guess,**last)
            write(OUT/'linear-failure.json',dict(classification='Counterexample candidate',error=repr(exc),
                call_begin=begin,call_end=len(calls),dimension=len(rhs),original_failure_preserved=True))
            raise
        finally:
            solves.append(dict(call_begin=begin,call_end=len(calls)))
            write(OUT/'linear-solves.json',dict(classification='Counterexample candidate',solves=solves))
    stage=owner.Model.run.__globals__['stages']
    modified=FunctionType(stage.__code__,dict(stage.__globals__,solve=solve))
    run=owner.Model.run
    owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,stages=modified),argdefs=run.__defaults__)
    m=owner.Model(64);m.run(64,'probe-64',61,restart='interval-15-64')
    raise AssertionError('Original terminal failure did not reproduce; reconsider before repair')


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));prior.owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            plan=read(OUT/'plan.json')
            for p,h in dict(plan['bindings'],**plan['reused']).items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
