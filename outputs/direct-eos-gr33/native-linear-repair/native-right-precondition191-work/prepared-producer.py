"""Counterexample candidate: minimize the physical linear residual directly.

Reconstruct the saved failed Radau system without replaying accepted steps.
Right preconditioning changes the Krylov norm, not the stage equation or gates.
"""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
from scipy.sparse.linalg import LinearOperator
import finish_full_incident_horizon as prior

OUT=Path('native-right-precondition191-work');OLD=prior.OUT
owner=prior.prior.owner;joint=owner.joint;LD=joint.LD
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=25,check=125)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(OLD/'repair-coarse-receipt.json')['error'] is not None
    sources=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    reused={}
    for src in sources:
        dst=OUT/src.relative_to(OLD);os.link(src,dst);reused[str(dst)]=sha(src)
    sources += [OLD/'failed-linear-system.npz',OLD/'linear-failure.json',OLD/'linear-calls.json',OLD/'repair-coarse-receipt.json',OLD/'sweep-1/photons/interval-15-64.npz']
    sources += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Solve the stored failing physical linear system by minimizing its unpreconditioned residual in a right-preconditioned Krylov basis, while retaining the original four-refinement and all physical accuracy gates.',
        evidence='Left-preconditioned GMRES reported tiny internal residual while its true residual remained0.147/0.133/0.00933/0.000394 per correction. Merely raising the refinement count to7 failed later at2.455638e-14 vector residual. Preserve both failures; no further iteration increase.',
        reconstruction='The saved failed guess has identical photon and gas stage rows, so it is the first Newton proposal after one accepted substep. Reuse its initial photon and guide gas, apply the known floor map, and reconstruct the second substep after60/64. Require the actual RHS and guess to match the saved failing system exactly.',
        repair='Solve A P y=b-A x0 and return x0+P y, retaining this addition and residual in long double. Same P, original restart20/maxiter5/rtol1e-14, at most four extended refinements, original1e-14 vector/1e-13 four-moment gates. No physical step is accepted by this system check.',
        decision='Only a matching reconstructed system and all original residual gates admit a new, separately budgeted actual continuation. A small preconditioned norm alone is not acceptance.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,new_physical_steps=0,
        forecast='Constructor15s plus two native Jacobians and local coefficients about15s; original failing inner work40s. Allow125s for reconstruction and four right-preconditioned corrections. No blind full-run replay or new grid.',
        stop='Mismatch, original gate, four-correction limit or125s cap. No further iteration or solver ladder under this plan.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(sources)},reused=reused))


def initialize():
    prior.prior.prior.OUT=OUT;prior.prior.prior.initialize()


def right_solver(log):
    def right(op,rhs,**kwargs):
        P=kwargs.pop('M');initial=kwargs.pop('x0');callback=kwargs.pop('callback');history=[]
        x0=np.zeros_like(rhs,dtype=LD) if initial is None else np.asarray(initial,LD)
        residual=np.asarray(rhs,LD)-op.matvec(x0)
        system=LinearOperator(op.shape,lambda v:op.matvec(P.matvec(v)),dtype=float)
        def record(v):history.append(float(v));callback(v)
        start=time.monotonic()
        y,info=joint.gmres(system,np.asarray(residual,float),x0=None,callback=record,**kwargs)
        answer=x0+np.asarray(P.matvec(y),LD)
        r=np.asarray(rhs,LD)-op.matvec(answer)
        log.append(dict(info=int(info),seconds=time.monotonic()-start,iterations=len(history),history=history,
            true_relative=float(np.linalg.norm(r)/max(np.linalg.norm(rhs),1e-290))))
        return answer,info
    return FunctionType(joint.solve.__code__,dict(joint.solve.__globals__,gmres=right))


def check():
    initialize();m=owner.Model(64);saved=dict(np.load(OLD/'failed-linear-system.npz'))
    guess=saved['guess'].reshape(2,-1);assert np.array_equal(guess[0],guess[1])
    x,guide=m.unpack(guess[0]);x=x.copy();guide=guide.copy()
    with np.load(OLD/'sweep-1/photons/interval-15-64.npz') as p:parts=1+int(p['split_macro_steps'][60])
    assert parts==2
    h=m.t[-1]/64/parts;t=60*m.t[-1]/64+h
    g=guide.copy();g[~m.material.active(t)]=0.;m.guide_g=guide
    lus=[joint.splu(joint.sparse.eye(m.n*m.q,format='csc')-a*h*m.A) for a in [5/12,1/4]]
    captured={}
    class Captured(Exception):pass
    def intercept(m,op,P,rhs,initial):
        captured.update(op=op,P=P,rhs=rhs,guess=initial);raise Captured()
    stage=owner.Model.run.__globals__['stages'];stage=FunctionType(stage.__code__,dict(stage.__globals__,solve=intercept))
    try:stage(m,t,h,x,g,lus)
    except Captured:pass
    assert captured
    identical=np.array_equal(captured['rhs'],saved['rhs']) and np.array_equal(captured['guess'],saved['guess'])
    mapping=dict(classification='Counterexample candidate',passed=identical,time=float(t),step=float(h),
        RHS_relative=float(np.linalg.norm(captured['rhs']-saved['rhs'])/np.linalg.norm(saved['rhs'])),
        guess_relative=float(np.linalg.norm(captured['guess']-saved['guess'])/np.linalg.norm(saved['guess'])))
    write(OUT/'system-reconstruction.json',mapping);assert identical,mapping
    log=[];solve=right_solver(log)
    try:answer=solve(m,**captured)
    finally:write(OUT/'right-calls.json',dict(classification='Counterexample candidate',calls=log))
    residual=captured['rhs']-captured['op'].matvec(answer)
    relative=float(np.linalg.norm(residual)/np.linalg.norm(captured['rhs']))
    moments=joint.physical_norm(m,residual)/joint.scales(m,captured['rhs'],answer)
    np.savez_compressed(OUT/'accepted-linear-solution.npz',solution=answer,residual=residual)
    result=dict(classification='Counterexample candidate',passed=relative<1e-14 and max(moments)<1e-13,
        actual_saved_system_reconstructed=True,vector_relative=relative,physical_relative=moments.astype(float).tolist(),
        calls=log,new_physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'check-result.json',result);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
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
