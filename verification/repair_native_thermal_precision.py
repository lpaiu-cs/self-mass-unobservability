"""Preserve precision through the inverse conserved thermal coordinate."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import align_native_branch_evaluation as prior

OUT=Path('native-thermal-precision234-work');OLD=prior.OLD
read,write,sha=prior.read,prior.write,prior.sha
precision,joint,LD=prior.precision,prior.joint,prior.LD
CAPS=dict(prepare=180,check=900)


def build():
    engine=joint.previous.engine
    old='et=z.astype(LD).copy();et[2]-=(m.a.astype(LD)-m.model.m.a0)*m.model.cx*C*C*et[0]'
    new='et=hp(z).copy();et[2]-=hp((m.a.astype(LD)-m.model.m.a0)*m.model.cx*C*C)*et[0]'
    assert engine.source.count(old)==1
    copied=SimpleNamespace(**dict(vars(engine),source=engine.source.replace(old,new)))
    # Keep the original rounded background coefficient, but do the changing
    # conserved-state multiplication/subtraction before rounding its tiny result.
    return precision.build(SimpleNamespace(previous=SimpleNamespace(engine=copied)),prior.prior.prior.owner)


def prepare():
    result=read(prior.OUT/'result.json');assert result['exact_saved_actual_defect'] and result['decomposition_remainder']<1e-18
    assert not OUT.exists();OUT.mkdir();files=[]
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    for part in ['photons','material']:(OUT/'sweep-1'/part).mkdir(parents=True)
    files += [OLD/n for n in ['rejected-joint-stage.npz','rejected-joint-stage.json','last-accepted-64.npz']]
    files += [prior.OUT/n for n in ['plan.json','result.json','comparison.npz','check-receipt.json','symbolic.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Remove the identified internal thermal-coordinate quantizer from the true native baryon owner, verify its actual increment, then apply it to actual continuation.',
        evidence='233exactly reproduces231actual defect1.21893e-6and linear9.97529e-16. Independent native-minus-affine B increment explains the difference with3.89e-20remainder. Half/double tiny increments around the guide change the secant by up to1.36e-5in stage units despite40/80digit agreement. The high-precision tangent still executes et=z.astype(longdouble) and its cancelling inverse-energy subtraction before entering the high-precision primitive.',
        change='Keep z and inverse thermal multiplication/subtraction in arbitrary precision. Preserve the original rounded fixed background coefficient exactly; no new EOS coefficient, state, probe, branch rule, gate, grid or clock. The running232/224sources remain frozen.',
        controls='Reuse233actual guide and saved proposal, original Jacobian and full equation. Measure both modified linear and nonlinear residual without accepting the proposal. Compare40/80digits and half/double increments about the guide; repeat the original0.5/2constitutive probes. Keep233original reconstruction evidence unchanged.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=8,forecast='233reconstruction and independent increment checks complete in its bound receipt. This repeats that one pair with one precision change;15minute cap, no Krylov or accepted physical step.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(prior.OUT/'symbolic.json'))


def check():
    source=inspect.getsource(prior.check)
    changes=[("m=prior.prior.owner.Model(64);z=", "m=prior.prior.owner.Model(64);m.precise_tangent=build();z="),
        ("assert np.array_equal(actual,z['defect']),'Exact full actual defect reconstruction'", "old_actual=z['defect'];assert read(previous/'result.json')['exact_saved_actual_defect']"),
        (';assert np.linalg.norm(linear)/norm<1e-14 and max(physical)<1e-13',''),
        ('exact_saved_actual_defect=True,','original233_reconstruction_preserved=True,same_saved_proposal=True,'),
        ("np.savez_compressed(OUT/'comparison.npz'", "result.update(old_actual_relative=float(np.linalg.norm(old_actual)/norm),arithmetic_effect=float(np.linalg.norm(actual-old_actual)/norm),constitutive_B=constitutive(m,t,h,sol),linear_passed=bool(np.linalg.norm(linear)/norm<1e-14 and max(physical)<1e-13),actual_stage_passed=bool(np.linalg.norm(actual)/norm<1e-12 and max(joint.physical_norm(m,actual)/z['physical_scales'])<1e-13))\n    assert max(result['constitutive_B'])<.002\n    np.savez_compressed(OUT/'comparison.npz'")]
    for a,b in changes:assert source.count(a)==1,(a,source.count(a));source=source.replace(a,b)
    ns=dict(prior.check.__globals__,OUT=OUT,build=build,previous=prior.OUT,constitutive=constitutive)
    exec(compile(source,__file__,'exec'),ns);(OUT/'expanded-check.py').write_text(source);ns['check']()


def constitutive(m,t,h,sol):
    result=[]
    with precision.mp.workdps(80):
        for c,row in zip(joint.C,sol.reshape(2,-1)):
            _,g=m.unpack(row);now=t+c*h;base=precision.native_B(m,now,g,m.precise_tangent)
            for probe in [.5,2.]:
                value=precision.native_B(m,now,g,m.precise_tangent,probe)
                result.append(float(np.sum(abs(value-base))/max(np.sum(abs(base)),precision.mp.mpf('1e-290'))))
    return result


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
