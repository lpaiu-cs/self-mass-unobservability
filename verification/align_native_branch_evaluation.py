"""Locate the actual nonlinear defect in the saved118th Newton equation."""
from pathlib import Path
from types import FunctionType
import json,os,resource,sys,time
import numpy as np
import continue_precise_conserved_native as prior
import continue_precise_momentum as dynamics

OUT=Path('native-branch-consistency233-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
joint=prior.joint;LD=joint.LD;flux=prior.flux;precision=prior.repair.precision
CAPS=dict(prepare=180,check=900)


def prepare():
    assert read(OLD/'failure-64.json')['actual_accepted_steps']==117
    assert not OUT.exists();OUT.mkdir();files=[]
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    for part in ['photons','material']:(OUT/'sweep-1'/part).mkdir(parents=True)
    files += [OLD/n for n in ['rejected-joint-stage.npz','rejected-joint-stage.json','last-accepted-64.npz','coarse-receipt.json','linear-64.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Identify the missing actual native increment in the failed118Newton equation before another coarse integration.',
        evidence='231passed8linear solves after precise conserved input but failed actual1.21893e-6. A further unchanged Newton budget increase lacks evidence of contraction.',
        method='Reconstruct the last Newton guide, matrix, RHS and actual defect. Require the saved full actual defect bit identity and original linear gates. At40/80digits compare the independently evaluated native B increment to J times the exact gas increment, and split its Radau defect from the linear residual. Test half/double increments around the same guide to distinguish a wrong local affine prediction from branch crossing or hidden quantization. Same stored input, no physical step or accepted state.',
        decision='If an affine mismatch accounts for the defect, repair the shared native branch/derivative evaluation and apply it to actual continuation. If the true increment is nonlinear, retain that fact and use an actual nonlinear correction. Do not claim the saved proposal is accepted.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=8,forecast='231eight actual proposals263.89s; this reconstructs one pair and a few true native evaluations. Allow15minutes, no Krylov and no accepted prefix replay.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,precision=1e-25),bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as s
    x,g,J,F,G,h=s.symbols('x g J F G h')
    assert s.expand((x-h*F)-(x-h*(G+J*(x-g)))+h*(F-G-J*(x-g)))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Exact nonlinear-minus-affine residual identity only; no physical convergence theorem.'))


def check():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    m=prior.prior.owner.Model(64);z=np.load(OLD/'rejected-joint-stage.npz');t,h=z['time'][()],z['step'][()]
    sol=z['solution'];v=z['initial'];guides=z['guides'];dim=len(v)
    maps=[m.jacobian(t+c*h,g) for c,g in zip(joint.C,guides)];Js=[p[0] for p in maps]
    cs=[m.local(t+c*h) for c in joint.C];ss=[m.source(t+c*h) for c in joint.C]
    def L(j,value):
        x,g=m.unpack(value);ph,q,*_=m.collision(cs[j],x,g)
        return m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph,q+(Js[j]@g.ravel()).reshape(m.n,4))
    def mat(value):
        vv=value.reshape(2,dim);return (vv-h*(joint.A@np.array([L(j,row) for j,row in enumerate(vv)]))).ravel()
    raw=joint.LinearOperator((2*dim,)*2,mat,dtype=float);op=dynamics.stable_operator(m,raw,[],True)
    affine=[base-(J@g.ravel()).reshape(m.n,4) for (J,base),g in zip(maps,guides)]
    src=np.array([m.pack(s[0]/(m.scale*joint.AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+a) for s,c,a in zip(ss,cs,affine)])
    rhs=(np.tile(v,(2,1))+h*(joint.A@src)).ravel();rhs=flux.precise_rhs(m,t,h,v,maps,guides,rhs)
    linear=rhs-op.matvec(sol);norm=np.linalg.norm(rhs);rates=[];m.precise_values={}
    for j,row in enumerate(sol.reshape(2,dim)):
        x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True)
        native=flux.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
    actual=(sol.reshape(2,dim)-v-h*(joint.A@np.array(rates))).ravel();actual=flux.precise_defect(m,t,h,v,sol,actual)
    assert np.array_equal(actual,z['defect']),'Exact full actual defect reconstruction'
    physical=joint.physical_norm(m,linear)/z['physical_scales'];assert np.linalg.norm(linear)/norm<1e-14 and max(physical)<1e-13
    values=[];nonlinear=[]
    for digits in [40,80]:
        with precision.mp.workdps(digits):
            errors=[]
            for j,(c,guide,row) in enumerate(zip(joint.C,guides,sol.reshape(2,dim))):
                _,g=m.unpack(row);now=t+c*h;delta=precision.hp(g)-precision.hp(guide)
                f0=precision.native_B(m,now,guide,m.precise_tangent);f1=precision.native_B(m,now,g,m.precise_tangent)
                errors.append(f1-f0-precision.baryon_product(Js[j],delta))
                if digits==80:
                    for scale in [precision.mp.mpf('.5'),precision.mp.mpf(2)]:
                        fp=precision.native_B(m,now,precision.hp(guide)+scale*delta,m.precise_tangent)
                        curvature=precision.cast(fp-f0-scale*(f1-f0))
                        nonlinear.append(dict(stage=j,scale=float(scale),maximum_stage_effect=float(abs(h)*np.max(abs(curvature))/norm),norm_stage_effect=float(abs(h)*np.linalg.norm(curvature)/norm)))
            values.append(precision.cast(-precision.hp(h)*(precision.hp(joint.A)@np.array(errors))))
    baryon=np.array([m.unpack(row)[1][:,2] for row in (actual+linear).reshape(2,dim)])
    mismatch=values[-1];change=float(np.linalg.norm(values[0]-mismatch)/norm)
    result=dict(classification='Counterexample candidate',exact_saved_actual_defect=True,
        linear_relative=float(np.linalg.norm(linear)/norm),actual_relative=float(np.linalg.norm(actual)/norm),
        affine_native_B_mismatch=float(np.linalg.norm(mismatch)/norm),decomposition_remainder=float(np.linalg.norm(baryon-mismatch)/norm),precision_change=change,
        increment_linearity=nonlinear,new_physical_steps=0,physical_state_accepted=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'comparison.npz',linear=linear,actual=actual,affine_native_B_mismatch=mismatch)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert change<1e-25


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
