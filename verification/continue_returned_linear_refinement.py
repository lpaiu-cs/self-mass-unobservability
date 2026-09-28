"""Continue the same full return with a larger inner refinement allowance.

All physical, branch, Newton and residual gates are unchanged. Preserve the
original four-solve failure and every completed source/metric input.
"""
from pathlib import Path
import inspect,os,resource,sys,time
import numpy as np
import complete_returned_period as prior

ROOT=Path('native-returned-refinement253-work');OUT=ROOT
OLD=prior.ROOT/'full';read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind
joint=prior.actual.prior.joint
CAPS=dict(prepare=600,coarse=7200,fine=14400,audit=600)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    failure=read(prior.ROOT/'coarse-receipt.json')
    assert 'Four-moment linear residual' in failure['error']
    assert read(prior.ROOT/'controller-status.json')['state']=='failed'
    assert read(OLD/'metric-result.json')['passed']
    paths=list((OLD/'sweep-0').rglob('*.npz'))
    paths += [p for folder in ['metric','gr'] for p in (OLD/folder).iterdir() if p.is_file()]
    paths += [OLD/name for name in ['normalization.json','photon-conservation-plan.json','check-result.json','symbolic.json','metric-result.json','metric-receipt.json']]
    paths += [OLD/f'sweep-1/photons/interval-14-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    for p in paths:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    (OUT/'sweep-1/material').mkdir(parents=True,exist_ok=True)
    files=paths+[prior.ROOT/'coarse-receipt.json',prior.ROOT/'coarse.stderr.log',OLD/'plan.json',Path(__file__),
        OLD/'last-pair-64.npz',OLD/'recovered-64.npz',OLD/'sweep-1/photons/interval-15-64.npz']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    plan=read(OLD/'plan.json');plan.update(
        claim='Finish the actual same full-period returned solution after the late linear solver exhausted its inner refinement budget.',
        failure=failure,repair='Keep the exact operator, preconditioner, GMRES restart20/maxiter5 and long-double residual. Increase total inner solves from4to12 only; preserve the first four-solve failure with trace and rejected vector. Each accepted solve must still pass vector1e-14 and all physical1e-13; actual stage1e-12 and three Newton proposals unchanged.',
        evidence='249full metric passed. Coarse canonical15/16 completed111actualsteps;112pair captured but113linear solve failed1.881796669776e-6 with photon number8.591209781709e-6. No resource timeout. The primary high solution used a prior larger-budget repair; this low-return route still used the original four-solve function.',
        decision='Apply the repaired solve directly to the actual remaining coupled stages. If true residual stops decreasing or12inner solves/physical gates fail, preserve the actual rejected system. Do not lower accuracy or refine physical clocks.',
        replay='Reuse all metric/GR and103/199earlier physical steps. Replay only coarse canonical15/16(~161s) to restore transient native-branch/anchor diagnostics omitted by249failure. Require every stored interval15array and original112pair exactly unchanged. Fine starts at original199checkpoint; no long primary trajectory rerun.',
        forecast='249metric1285s is reused. Coarse original16steps estimated6..12min plus161s replay; fine32steps12..25min before extra Krylov costs. Allow2/4hours and16GiB. Actual later refinement convergence is unmeasured.',
        budgets=CAPS,max_inner_solves=12,max_Newton=3,CPU_affinity=3,virtual_GiB=16,
        scientific_gates_changed=False,full_declared_period=True,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    plan['bindings'].update({str(p):sha(p) for p in dict.fromkeys(files)});write(OUT/'plan.json',plan)


def initialize():
    Model=bind(prior.initialize,OUT=OUT)()
    original=joint.solve
    source=inspect.getsource(original)
    source=prior.base.replace(source,'range(4)','range(12)')
    source=prior.base.replace(source,"        assert k<3,('Four-moment linear residual',relative,moments.tolist())",
        "        if k==3: record_failure(m,rhs,guess,sol,residual,relative,moments)\n        assert k<11,('Four-moment linear residual',relative,moments.tolist())")
    calls=[];solves=[];current=[0]
    def gmres(op,rhs,**options):
        mark=time.monotonic();history=[];callback=options['callback']
        def record(value):history.append(float(value));callback(value)
        answer,info=joint.gmres(op,rhs,**dict(options,callback=record))
        defect=np.asarray(rhs,joint.LD)-op.matvec(answer.astype(joint.LD))
        calls.append(dict(info=int(info),seconds=time.monotonic()-mark,iterations=len(history),history=history,
            true_relative=float(np.linalg.norm(defect)/max(np.linalg.norm(rhs),1e-290))))
        write(OUT/f'linear-{current[0]}.json',dict(calls=calls,solves=solves));return answer,info
    def record_failure(m,rhs,guess,sol,residual,relative,moments):
        label=f'original-limit-{current[0]}-{len(solves)}'
        np.savez_compressed(OUT/(label+'.npz'),rhs=rhs,guess=guess,solution=sol,residual=residual)
        write(OUT/(label+'.json'),dict(classification='Counterexample candidate',original_four_solve_failed=True,
            relative=relative,moments=[float(v) for v in moments],call_begin=solves[-1]['begin'],call_end=len(calls)))
    ns=dict(original.__globals__,gmres=gmres,record_failure=record_failure);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-refinement.py').write_text(source)
    def solve(m,op,P,rhs,guess):
        current[0]=int(m.anchor['base_steps']) if 'base_steps' in m.anchor else len(m.anchor['actual_step_edges'])
        solves.append(dict(begin=len(calls)))
        value=ns['solve'](m,op,P,rhs,guess)
        residual=rhs-op.matvec(value)
        solves[-1].update(end=len(calls),relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290)),
            moments=[float(v) for v in joint.physical_norm(m,residual)/joint.scales(m,rhs,value)])
        write(OUT/f'linear-{current[0]}.json',dict(calls=calls,solves=solves));return value
    run=Model.run;stage=run.__globals__['stages'];stage=bind(stage,solve=solve)
    Model.run=bind(run,stages=stage);return Model


def evolve(n):
    source=inspect.getsource(prior.evolve)
    marker="        checks=dict(newton=m.newton_iterations,stages=m.stage_log);restart=label;del m;gc.collect()"
    source=prior.base.replace(source,marker,"        if n==64 and j==15:\n            replay_check(z)\n"+marker)
    marker="        np.savez_compressed(OUT/f'last-pair-{n}.npz',time=t,step=h,photons=[v[0] for v in pair],gas=[v[1] for v in pair],x_initial=x,g_initial=g)"
    source=prior.base.replace(source,marker,marker+"\n        if n==64 and len(moments)//2==112: pair_check(t,h,x,g,pair)")
    def replay_check(z):
        with np.load(OLD/'sweep-1/photons/interval-15-64.npz') as old:
            assert set(old.files)==set(z.files)
            for k in old.files:assert np.array_equal(old[k],z[k]),('Accepted interval replay changed',k)
        write(OUT/'replay-111.json',dict(classification='Counterexample candidate',passed=True,every_saved_array_exact=True,new_physical_model=False))
    def pair_check(t,h,x,g,pair):
        with np.load(OLD/'last-pair-64.npz') as old:
            for k,v in dict(time=t,step=h,x_initial=x,g_initial=g,photons=[v[0] for v in pair],gas=[v[1] for v in pair]).items():
                assert np.array_equal(old[k],v),('Accepted112pair changed',k)
        write(OUT/'replay-112.json',dict(classification='Counterexample candidate',passed=True,actual_pair_exact=True))
    ns=dict(prior.evolve.__globals__,OUT=OUT,initialize=initialize,replay_check=replay_check,pair_check=pair_check)
    exec(compile(source,__file__,'exec'),ns);ns['evolve'](n)


def audit():
    bind(prior.audit,OUT=OUT)()
    r=read(OUT/'result.json');r.update(original_linear_failure_preserved=True,original_metric_reused=True,
        scientific_gates_changed=False,max_inner_solves=12)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        evolve(64 if action=='coarse' else 128) if action in ['coarse','fine'] else globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
