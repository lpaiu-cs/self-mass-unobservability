"""Correct the saved118th system using joint integer representable coordinates."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
from scipy.optimize import Bounds,LinearConstraint,milp
import continue_rounding_aware_native as evolution

rounding=evolution.rounding;system_owner=rounding.prior
OUT=Path('native-integer228-work');OLD=evolution.OUT
read,write,sha,LD,joint=evolution.read,evolution.write,evolution.sha,evolution.LD,evolution.joint
CAPS=dict(prepare=180,solve=1200)


def prepare():
    assert read(OLD/'controller-status.json')['state']=='failed';assert read(OLD/'failure-64.json')['actual_accepted_steps']==117
    assert not OUT.exists();OUT.mkdir();files=[]
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    for name in ['photons','material']:(OUT/'sweep-1'/name).mkdir(parents=True)
    files += [OLD/n for n in ['failed-linear-64.npz','last-accepted-64.npz','failure-64.json','coarse-receipt.json','rounding-64.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Resolve the saved actual118linear system at unchanged precision, then use a passing proposal for actual continuation. No physical state is accepted by this linear repair alone.',
        evidence='218passed actual117and timed out after7202s before118acceptance. Paired/triple rounding searches repeatedly leave roughly2.2e-14against1e-14; all physical moment norms can pass without the original full-vector gate.',
        method='Reconstruct the SAME saved RHS and residual exactly using the existing216system and210stable B/S operator. Reuse its last accumulated solution. Project fine representable directions as217does; jointly choose ALL hard integer-ULP coordinates using installed SciPy/HiGHS, instead of pair/triple enumeration. Bounded1024ULPs, up to90s/5000nodes per search, at most4improving passes; test selected point and its one-ULP neighbours with the unchanged full operator and physical norms.',
        gates=dict(linear=1e-14,physical_linear=1e-13,actual_stage=1e-12,physical_stage=1e-13),budgets=CAPS,CPU_threads=1,virtual_GiB=8,
        forecast='Existing saved-system reconstruction20..45s;4integer searches capped90s each plus at most96full residual evaluations. Allow20minutes. No GMRES, new physical step, grid, period or accepted-prefix replay.',
        decision='A full original linear pass admits a Newton proposal only. Original nonlinear stage is evaluated separately; actual continuation must still pass it. Stop on nonimprovement or cap; no automatic repeated search.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(OLD/'symbolic.json'))


def integer_guesses(PH,pr,goal,logs):
    # Numerical solver coordinates only: all candidates face the original RHS.
    H=PH/float(goal);r=pr/float(goal);n=H.shape[1]
    assert n<=16,'Reassess an unexpectedly larger integer problem'
    matrix=np.r_[np.c_[H,-np.ones(len(r))],np.c_[-H,-np.ones(len(r))]]
    start=time.monotonic();result=milp(np.r_[np.zeros(n),1.],integrality=np.r_[np.ones(n),0],
        bounds=Bounds(np.r_[np.full(n,-1024.),0],np.r_[np.full(n,1024.),np.inf]),
        constraints=LinearConstraint(matrix,np.full(2*len(r),-np.inf),np.r_[r,-r]),
        options=dict(time_limit=90.,node_limit=5000,mip_rel_gap=.001))
    logs.append(dict(status=int(result.status),message=result.message,seconds=time.monotonic()-start,
        variables=n,objective=None if result.fun is None else float(result.fun),node_count=getattr(result,'mip_node_count',None)))
    if result.x is None:return {}
    point=np.rint(result.x[:n]).astype(int);assert np.max(abs(point))<=1024
    points=[point]
    for j in range(n):
        for shift in [-1,1]:
            p=point.copy();p[j]+=shift
            if abs(p[j])<=1024:points.append(p)
    return {tuple(p):float(np.linalg.norm(pr-PH@p)) for p in points}


def solve():
    m,op,rhs,z,p,cs,ss,t,h,v=FunctionType(system_owner.system.__code__,dict(system_owner.system.__globals__,OUT=OUT,OLD=OLD))()
    source=inspect.getsource(rounding.neighbours);a=source.index('        guesses={};count=');b=source.index('        best=sol;',a)
    source=source[:a]+"        order=list(range(len(hard)));guesses=integer_guesses(PH,pr,goal,integer_log)\n"+source[b:]
    integer_log=[];ns=dict(rounding.neighbours.__globals__,integer_guesses=integer_guesses,integer_log=integer_log)
    exec(compile(source,__file__,'exec'),ns);(OUT/'expanded-integer-corrector.py').write_text(source)
    log=[];sol=ns['neighbours'](m,op,rhs,z['solution'],log);residual=rhs-op.matvec(sol);norm=np.linalg.norm(rhs)
    physical=joint.physical_norm(m,residual)/joint.scales(m,rhs,sol);rates=[];m.precise_values={};owner=system_owner.prior.prior
    for j,row in enumerate(sol.reshape(2,-1)):
        x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True)
        native=owner.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
    defect=(sol.reshape(2,-1)-v-h*(joint.A@np.array(rates))).ravel();defect=owner.precise_defect(m,t,h,v,sol,defect)
    actual=joint.physical_norm(m,defect)/joint.scales(m,rhs,sol)
    result=dict(classification='Counterexample candidate',exact_original_RHS_and_residual=True,linear_relative=float(np.linalg.norm(residual)/norm),physical=physical.astype(float).tolist(),
        actual_stage_relative=float(np.linalg.norm(defect)/norm),actual_stage_physical=actual.astype(float).tolist(),
        linear_passed=bool(np.linalg.norm(residual)/norm<1e-14 and max(physical)<1e-13),actual_stage_passed=bool(np.linalg.norm(defect)/norm<1e-12 and max(actual)<1e-13),
        physical_state_accepted=False,new_physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'integer-search.json',integer_log);write(OUT/'full-operator-checks.json',log)
    np.savez_compressed(OUT/'proposal.npz',solution=sol,residual=residual,actual_defect=defect);write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
