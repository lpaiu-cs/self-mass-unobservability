"""Continue the actual return after short restarted Krylov solves stagnated.

Only failed inner solves gain a larger Krylov space. The operator, physical
equations, clocks, nonlinear proposals and acceptance tolerances stay fixed.
"""
from pathlib import Path
import resource,sys,time
import numpy as np
import continue_returned_linear_refinement as prior

OUT=Path('native-returned-krylov254-work')
read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind
CAPS=prior.CAPS


def prepare():
    failed=read(prior.OUT/'coarse-receipt.json')
    assert 'Four-moment linear residual' in failed['error']
    assert read(prior.OUT/'controller-status.json')['state']=='failed'
    bind(prior.prepare,OUT=OUT)()
    p=read(OUT/'plan.json')
    p.update(previous_twelve_solve_failure=failed,
        repair='Retain the original GMRES(restart20,maxiter5). On info>0 continue its same solution using restart80,maxiter20 and the identical rtol,atol,M,operator,RHS. Keep at most12 extended-residual solves, three Newton proposals and all original physical gates.',
        evidence='253actual113 first Newton solve passed extended vector2.20249e-19 and all moments2.71950e-19 after12solves. The next Newton system still failed vector1.26740e-5; its last corrections retained0.92..0.99 of their own true residual. More identical small-space restarts alone are not established as an adequate repair.',
        decision='Execute the actual remaining coupled interval directly; accept only original linear, true nonlinear native, branch, conservation and paired-time gates. Preserve failed20-space and expanded80-space traces; any12solve or wall cap stops without new physical clocks.',
        forecast='Reuse the completed1285s metric and all primary histories. Replay only the161s coarse canonical15 interval to restore unstored transient state and require exact111/112reproduction. Ordinary successful inner calls stay unchanged. Measured100short-space iterations take10..12s;1600expanded iterations assumed3..5minutes per hard system. Allow2hours coarse/4hours fine,16GiB,CPU3. Later convergence is unmeasured.',
        fallback_restart=80,fallback_maxiter=20,scientific_gates_changed=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    p['bindings'].update({str(v):sha(v) for v in [Path(__file__),prior.OUT/'coarse-receipt.json',prior.OUT/'linear-120.json',prior.OUT/'replay-111.json',prior.OUT/'replay-112.json']})
    write(OUT/'plan.json',p)


def initialize():
    original=prior.joint.gmres;events=[]
    def gmres(op,rhs,**options):
        answer,info=original(op,rhs,**options)
        if info>0:
            residual=rhs-op.matvec(answer.astype(prior.joint.LD))
            row=dict(initial_info=int(info),initial_true_relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290)),
                state='running',restart=80,maxiter=20,dimension=len(rhs))
            events.append(row);write(OUT/'krylov-expansion.json',dict(events=events))
            start=time.monotonic();iterations=[];callback=options['callback']
            def record(value):iterations.append(float(value));callback(value)
            answer,info=original(op,rhs,**dict(options,x0=answer,restart=80,maxiter=20,callback=record))
            residual=rhs-op.matvec(answer.astype(prior.joint.LD))
            row.update(state='completed',info=int(info),seconds=time.monotonic()-start,iterations=len(iterations),history=iterations,
                true_relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290)))
            write(OUT/'krylov-expansion.json',dict(events=events))
        return answer,info
    prior.joint.gmres=gmres
    return bind(prior.initialize,OUT=OUT)()


def evolve(n):
    bind(prior.evolve,OUT=OUT,initialize=initialize)(n)


def audit():
    bind(prior.audit,OUT=OUT)()
    r=read(OUT/'result.json');r.update(failed_short_krylov_expanded=True,original_twelve_solve_failure_preserved=True)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        evolve(64 if action=='coarse' else 128) if action in ['coarse','fine'] else globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
