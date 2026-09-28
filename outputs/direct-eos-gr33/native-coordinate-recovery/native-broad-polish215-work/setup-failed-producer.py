"""Test the saved117th system with a broader correction trigger, same gates."""
from pathlib import Path
from types import FunctionType
import inspect,json,resource,time
import numpy as np
import continue_precise_momentum as current

OUT=Path('native-broad-polish215-work');OLD=current.OUT;OUT.mkdir(exist_ok=True)
LD=current.LD;read,write,sha=current.read,current.write,current.sha
assert not (OUT/'result.json').exists();start=time.monotonic();error=None
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));current.joint.previous.original.inf.incident.native.deadline(600)
write(OUT/'plan.json',dict(classification='Conjectural',claim='Apply the existing restricted residual correction to the retained117th proposal, whose4.094e-10vector residual was outside the old1e-10correction trigger. Test actual improvement at the unchanged1e-14linear and1e-12nonlinear gates before any new evolution.',
    method='Keep24columns,4passes and all full-vector/physical improvement requirements. Widen only the correction eligibility ceiling from1e-10to1e-8; it is not an acceptance gate. Reconstruct the exact210RHS and residual. No new GMRES or physical step.',
    cap_seconds=600,virtual_GiB=6,CPU_threads=1,forecast='Constructor and2Jacobians plus bounded24columns/4passes; earlier such reconstruction took20..30s. Allow10minutes. If it does not solve the actual saved defect, preserve that failure and do not rerun the long system merely by adding iterations.',
    bindings={str(p):sha(p) for p in [Path(__file__),Path(current.__file__),OLD/'failed-linear-64.npz',OLD/'last-accepted-64.npz',OLD/'failure-64.json',OLD/'coarse-receipt.json']},final_charge_conclusion='unadjudicated',full_goal_complete=False))
try:
    FunctionType(current.initialize.__code__,dict(current.initialize.__globals__,OUT=OUT))(False)
    m=current.owner.Model(64);z=dict(np.load(OLD/'failed-linear-64.npz'));p=dict(np.load(OLD/'last-accepted-64.npz'));t,h=p['next_time'][()],p['next_step'][()]
    v=m.pack(p['x'],p['g']);dim=len(v);joint=current.joint;prior=current.prior
    guides=[m.unpack(row)[1] for row in z['guess'].reshape(2,-1)];maps=[m.jacobian(t+c*h,g) for c,g in zip(joint.C,guides)];Js=[r[0] for r in maps]
    cs=[m.local(t+c*h) for c in joint.C];ss=[m.source(t+c*h) for c in joint.C]
    def L(j,value):
        x,g=m.unpack(value);ph,q,*_=m.collision(cs[j],x,g)
        return m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph,q+(Js[j]@g.ravel()).reshape(m.n,4))
    def mat(value):
        x=value.reshape(2,dim);return (x-h*(joint.A@np.array([L(j,row) for j,row in enumerate(x)]))).ravel()
    raw=joint.LinearOperator((2*dim,)*2,mat,dtype=float);op=current.stable_operator(m,raw,[],True)
    affine=[b-(J@g.ravel()).reshape(m.n,4) for (J,b),g in zip(maps,guides)]
    src=np.array([m.pack(s[0]/(m.scale*joint.AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+a) for s,c,a in zip(ss,cs,affine)])
    rhs=(np.tile(v,(2,1))+h*(joint.A@src)).ravel();rhs=prior.precise_rhs(m,t,h,v,maps,guides,rhs)
    assert np.array_equal(rhs,z['rhs']);original=rhs-op.matvec(z['solution']);assert np.array_equal(original,z['residual'])
    source=inspect.getsource(current.polish_function).replace('source=inspect.getsource(legacy.polish)',"source=inspect.getsource(legacy.polish).replace('if norm>goal*10000:break','if norm>goal*1000000:break')")
    ns=dict(current.polish_function.__globals__,OUT=OUT);exec(compile(source,__file__,'exec'),ns);polish=ns['polish_function']()
    sol=polish(m,op,rhs,z['solution']);res=rhs-op.matvec(sol);norm=np.linalg.norm(rhs);physical=joint.physical_norm(m,res)/joint.scales(m,rhs,sol)
    rates=[];m.precise_values={}
    for j,row in enumerate(sol.reshape(2,dim)):
        x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True);native=prior.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
    defect=(sol.reshape(2,dim)-v-h*(joint.A@np.array(rates))).ravel();defect=prior.precise_defect(m,t,h,v,sol,defect)
    actual=joint.physical_norm(m,defect)/joint.scales(m,rhs,sol)
    result=dict(classification='Counterexample candidate',exact_RHS_and_residual=True,old_linear_relative=float(np.linalg.norm(original)/norm),linear_relative=float(np.linalg.norm(res)/norm),physical=physical.astype(float).tolist(),actual_stage_relative=float(np.linalg.norm(defect)/norm),actual_physical=actual.astype(float).tolist(),
        linear_passed=bool(np.linalg.norm(res)/norm<1e-14 and max(physical)<1e-13),actual_stage_passed=bool(np.linalg.norm(defect)/norm<1e-12 and max(actual)<1e-13),new_Krylov_iterations=0,new_physical_steps=0,final_charge_conclusion='unadjudicated')
    np.savez_compressed(OUT/'proposal.npz',solution=sol,residual=res,actual_defect=defect);write(OUT/'result.json',result);print(json.dumps(result))
except BaseException as exc:error=repr(exc);raise
finally:write(OUT/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
