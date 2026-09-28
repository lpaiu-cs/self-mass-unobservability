"""Compare inner restart budgets on one stored actual coupled stage equation.

This does not modify or restart the live258 calculation. Both proposals face
the existing extended-precision vector, physical and material residual gates.
"""
from pathlib import Path
import hashlib,inspect,os,resource,sys,time
import numpy as np
import repair_returned_metric_endpoint as live

OUT=Path('native-short-krylov259-work');OLD=live.OUT
read,write,sha,bind=live.read,live.write,live.sha,live.bind
joint=live.joint;LD=joint.LD
CAPS=dict(prepare=120,check=1200)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=list((OLD/'sweep-0').rglob('*.npz'))
    files += [p for part in ['metric','gr'] for p in (OLD/part).iterdir() if p.is_file()]
    files += [OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    for part in ['sweep-1/photons','sweep-1/material']:(OUT/part).mkdir(parents=True,exist_ok=True)
    pair=OLD/'last-pair-64.npz';before=sha(pair);raw=pair.read_bytes();after=sha(pair)
    assert before==after==hashlib.sha256(raw).hexdigest(),'Live pair changed during copy; do not use it'
    (OUT/'stored-pair.npz').write_bytes(raw)
    with np.load(OUT/'stored-pair.npz') as z:
        actual_time=float(z['time']);step=float(z['step'])
    files += [OUT/'stored-pair.npz',Path(__file__)]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Avoid double-precision inner stagnation by returning after one80-vector restart to the existing long-double refinement, without relaxing a single residual gate.',
        evidence='Live258calls used1360/1223inner iterations and194.7/170.9seconds before another outer correction. Accepted actual steps remain valid; this is a cost test, not a physical repair.',
        reconstruction='Use a byte-stable copy of an accepted actual258pair and its initial state. Linearize about its accepted material stages. This is the same saved stage equation, but is not claimed to reproduce its original first-Newton matrix. Require the stored pair to satisfy the original true nonlinear and material gates before comparison.',
        decision='Compare one versus twenty inner restarts on identical operators,RHS and guess; both retain at most12outer refinements, vector1e-14, physical1e-13 and material1e-13. A passing cost comparison alone cannot authorize a physical readout. Leave the live258source,plan and process unchanged.',
        budgets=CAPS,CPU_affinity=2,virtual_GiB=16,
        forecast='Constructor/native maps are assumed15..60s; live258inner solves measured3..195s. Two saved-system proposals are assumed1..8minutes; allow20minutes. Stop on a mismatch,original residual failure or cap; no extra grid,clock,path or retry ladder.',
        stored_time=actual_time,stored_step=step,new_physical_steps=0,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},
        final_charge_conclusion='unadjudicated',full_goal_complete=False))


def check():
    Model=bind(live.initialize,OUT=OUT)();m=Model(64)
    saved=dict(np.load(OUT/'stored-pair.npz'));t=saved['time'].item();h=saved['step'].item()
    x=saved['x_initial'];g=saved['g_initial'];times=np.array([t+c*h for c in joint.C])
    def guide(now):
        i=int(np.argmin(abs(times-now)));assert abs(times[i]-now)<1e-18
        return saved['gas'][i].copy()
    m.guide=guide
    stage=Model.run.__globals__['stages'];wrapped=stage.__globals__['solve']
    linear=inspect.getclosurevars(wrapped).nonlocals['ns']['solve']
    inner=linear.__globals__['gmres'];calls=inspect.getclosurevars(inner).nonlocals['log']
    lus=[joint.splu(joint.sparse.eye(m.n*m.q,format='csc')-a*h*m.A) for a in [5/12,1/4]]
    captured={}
    class Captured(Exception):pass
    def intercept(model,op,P,rhs,guess):
        captured.update(op=op,P=P,rhs=rhs,guess=guess);raise Captured()
    try:bind(stage,solve=intercept)(m,t,h,x,g,lus)
    except Captured:pass
    assert captured
    stored=np.array([m.pack(xx,gg) for xx,gg in zip(saved['photons'],saved['gas'])]).ravel()
    rates=[]
    for now,xx,gg in zip(times,saved['photons'],saved['gas']):
        cs=m.local(now);ss=m.source(now);p,q,*_=m.collision(cs,xx,gg,True)
        rates.append(m.pack((m.A@xx.reshape(m.n*m.q,m.nf)).reshape(xx.shape)+p+ss[0]/(m.scale*joint.AMP),q+m.native(now,gg)))
    defect=(stored.reshape(2,-1)-m.pack(x,g)-h*(joint.A@np.array(rates))).ravel()
    def gates(answer,residual,vector_gate):
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(captured['rhs']),LD('1e-290')))
        moments=(joint.physical_norm(m,residual)/joint.scales(m,captured['rhs'],answer)).astype(float).tolist()
        gas=live.prior.gas_relative(m,residual,answer)
        return dict(passed=relative<vector_gate and max(moments)<1e-13 and max(gas)<1e-13,vector_relative=relative,physical_relative=moments,material_relative=gas)
    reconstruction=gates(stored,defect,1e-12)
    write(OUT/'reconstruction.json',dict(classification='Counterexample candidate',**reconstruction,
        same_actual_stage_equation=True,original_first_Newton_matrix_reproduced=False,new_physical_steps=0))
    assert reconstruction['passed'],reconstruction
    def limited(*args,**kwargs):return inner(*args,**dict(kwargs,maxiter=1))
    rows=[];answers=[]
    for name,fn in [('one_restart',bind(linear,gmres=limited)),('twenty_restarts',linear)]:
        calls.clear();started=time.monotonic()
        try:answer=fn(m,**captured)
        finally:write(OUT/f'{name}-calls.json',dict(classification='Counterexample candidate',calls=calls))
        residual=captured['rhs']-captured['op'].matvec(answer)
        row=dict(name=name,seconds=time.monotonic()-started,iterations=sum(v['iterations'] for v in calls),inner_calls=len(calls),**gates(answer,residual,1e-14))
        rows.append(row);write(OUT/'comparison-progress.json',dict(rows=rows));assert row['passed'],row
        answers.append(answer);np.savez_compressed(OUT/f'{name}.npz',solution=answer,residual=residual)
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,
        speedup=rows[1]['seconds']/rows[0]['seconds'],solution_relative_difference=float(np.linalg.norm(answers[0]-answers[1])/np.linalg.norm(answers[1])),
        identical_matrix_rhs_guess=True,physical_step_accepted=False,live258_changed=False,new_physical_steps=0,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'result.json',result);print(result,flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            assert not (OUT/f'{action}-receipt.json').exists()
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
