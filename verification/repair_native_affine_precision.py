"""Counterexample candidate: make the native affine RHS and matrix consistent."""
from pathlib import Path
from types import FunctionType
from decimal import Decimal,localcontext
import inspect,json,os,resource,sys,time
import numpy as np
import polish_full_native_roundoff as prior

OUT=Path('native-affine-precision198-work');OLD=prior.OUT
owner,joint,LD=prior.owner,prior.joint,prior.LD
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=60,check=180,primitive_check=180)


def exact(v):
    p,q=v.as_integer_ratio();return Decimal(p)/Decimal(q)


def baryon_rows(J,gas):
    values=[exact(v) for v in gas.ravel()]
    return [sum(exact(J.data[p])*values[int(J.indices[p])] for p in range(J.indptr[k],J.indptr[k+1])) for k in range(2,J.shape[0],4)]


def coherent_rhs(m,h,v,maps,guides):
    with localcontext() as context:
        context.prec=80
        affine=[[exact(base[k,2])-a for k,a in enumerate(baryon_rows(J,g))] for (J,base),g in zip(maps,guides)]
        _,initial=m.unpack(v)
        return np.array([[LD(str(exact(initial[k,2])+exact(h)*sum(exact(joint.A[i,j])*affine[j][k] for j in range(2)))) for k in range(m.n)] for i in range(2)])


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(OLD/'controller-status.json')['state']=='failed'
    assert 'True native joint Radau equation' in read(OLD/'coarse-receipt.json')['error']
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for p in files:dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['rejected-joint-stage.npz','rejected-joint-stage.json','last-accepted-64.npz','coarse-receipt.json','linear-64.json','polish.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='cba0e6374',
        claim='Locate and repair the inconsistency between the actual native equation and the accurately evaluated linear matrix, reusing the last197Newton proposal.',
        evidence='All8linear solves passed; the actual native stage still failed at5.186855e-11. The0.0271568defect norm is almost entirely B, with0.0251465at cell261. Iteration did not contract monotonically. The B matrix is evaluated as80digit preassembled rows but the affine RHS still forms base-J*guide through cancelling long-double sparse sums.',
        method='Reconstruct only the stored last system. Compare the unchanged original RHS,80digit consistent affine RHS,80digit selected-branch prediction and the separately replayed true native RHS. No Krylov work or accepted physical step in this check.',
        decision='If affine arithmetic accounts for the discrepancy, apply the consistent RHS to actual continuation using the stored last Newton proposal. If the selected branch/native mismatch dominates instead, repair that owner before another long solve. No tolerance change.',
        forecast='Previous same-model reconstruction cost about20..30s. This check evaluates2Jacobians,2native RHS and sparse80digit B rows only. Allow180s with6GiB and1CPU; no replay of113accepted steps.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,new_physical_steps=0,new_Krylov_iterations=0,
        gates=dict(linear=1e-14,stage=1e-12,physical_stage=1e-13),
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))


def initialize():
    init=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT));init()


def check():
    initialize();m=owner.Model(64);z=dict(np.load(OLD/'rejected-joint-stage.npz'))
    t,h=z['time'][()],z['step'][()];guides=z['guides'];v=z['initial'];sol=z['solution'].reshape(2,-1)
    maps=[m.jacobian(t+c*h,g) for c,g in zip(joint.C,guides)]
    gas=[m.unpack(row)[1] for row in sol];native=np.array([m.native(t+c*h,g) for c,g in zip(joint.C,gas)])
    assert np.array_equal(native,z['native_rates']),'The unchanged native RHS must replay exactly'
    _,initial=m.unpack(v)
    affine=np.array([base-(J@g.ravel()).reshape(m.n,4) for (J,base),g in zip(maps,guides)])
    oldrhs=(np.tile(initial,(2,1,1))+h*np.einsum('ij,jnk->ink',joint.A,affine))[:,:,2]
    newrhs=coherent_rhs(m,h,v,maps,guides)
    with localcontext() as context:
        context.prec=80
        products=[baryon_rows(J,g) for (J,_),g in zip(maps,gas)]
        oldproducts=[baryon_rows(J,g) for (J,_),g in zip(maps,guides)]
        predicted=np.array([[LD(str(exact(base[k,2])+products[j][k]-oldproducts[j][k])) for k in range(m.n)] for j,(_,base) in enumerate(maps)])
        lhs=np.array([[LD(str(exact(gas[i][k,2])-exact(h)*sum(exact(joint.A[i,j])*products[j][k] for j in range(2)))) for k in range(m.n)] for i in range(2)])
        native_defect=np.array([[LD(str(exact(gas[i][k,2])-exact(initial[k,2])-exact(h)*sum(exact(joint.A[i,j])*exact(native[j,k,2]) for j in range(2)))) for k in range(m.n)] for i in range(2)])
    norm=np.linalg.norm(z['defect'])/LD(read(OLD/'rejected-joint-stage.json')['equations'][-1]['relative'])
    stored=np.array([m.unpack(row)[1][:,2] for row in z['defect'].reshape(2,-1)])
    mismatch=h*joint.A@(native[:,:,2]-predicted)
    measure=lambda q:float(np.linalg.norm(q)/norm)
    result=dict(classification='Counterexample candidate',native_replay_exact=True,
        stored_actual_B_relative=measure(stored),original_affine_B_linear_relative=measure(oldrhs-lhs),
        consistent_affine_B_linear_relative=measure(newrhs-lhs),affine_RHS_arithmetic_change=measure(newrhs-oldrhs),
        native_minus_selected_affine_stage_effect=measure(mismatch),
        exact_outer_stage_sum_B_relative=measure(native_defect),outer_stage_sum_error=measure(stored-native_defect),
        original_gate=1e-12,new_physical_steps=0,new_Krylov_iterations=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'affine-comparison.npz',old_rhs=oldrhs,consistent_rhs=newrhs,lhs=lhs,
        predicted_native=predicted,actual_native=native[:,:,2],consistent_native_defect=native_defect,stored_defect=stored)
    write(OUT/'check-result-new.json',result);print(json.dumps(result),flush=True)


def install_primitive_precision():
    engine=joint.previous.engine;face=engine.face;source=face.source
    assert source.count('np.zeros(m.n)')==5 and source.count(',float)')==5
    source=source.replace('np.zeros(m.n)','np.zeros(m.n,dtype=LD)').replace(',float)',',LD)')
    ns=dict(face.primitive.__globals__);exec(compile(source,__file__,'exec'),ns)
    precise=ns['primitive']
    engine.tangent=FunctionType(engine.tangent.__code__,dict(engine.tangent.__globals__,primitive=precise))
    assert engine.tangent.__globals__['primitive'] is precise
    (OUT/'expanded-precise-native-primitive.py').write_text(source)


def primitive_check():
    old=read(OUT/'check-result-new.json')
    assert old['native_replay_exact'] and old['native_minus_selected_affine_stage_effect']>100*old['affine_RHS_arithmetic_change']
    write(OUT/'primitive-plan.json',dict(classification='Conjectural',
        rejected_hypothesis='Consistent affine RHS changes the stage by1.468e-13; it does not explain the5.184e-11true native mismatch. Preserve that result and do not deploy the unhelpful affine-only repair.',
        repair='The actual native primitive variation allocates five binary64 output arrays and explicitly casts five results back to binary64, despite long-double state/flux arithmetic. Promote only these intermediate/output conversions to long double in the actual and branch-selected native owner. Preserve stored EOS coefficients, photon maps, equations, branch rules and gates.',
        check='Re-evaluate the saved last failed pair, its selected Jacobian, and original half/nominal/double constitutive probes. Separate native arithmetic changes from true branch variation. No new Krylov or physical step. Success admits actual continuation; past accepted-history consistency must also be checked.',
        budget_seconds=180,virtual_GiB=6,CPU_threads=1,forecast='Same reconstruction just took20.79s; promotion changes only small native arrays. Allow180s including independent branch/constitutive comparisons.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'affine-producer.py',OUT/'check-result-new.json',OLD/'rejected-joint-stage.npz']},final_charge_conclusion='unadjudicated'))
    initialize();m=owner.Model(64);z=dict(np.load(OLD/'rejected-joint-stage.npz'))
    t,h=z['time'][()],z['step'][()];v=z['initial'];sol=z['solution'].reshape(2,-1);guides=z['guides']
    gas=[m.unpack(row)[1] for row in sol];_,initial=m.unpack(v)
    original=np.array([m.native(t+c*h,g) for c,g in zip(joint.C,gas)])
    assert np.array_equal(original,z['native_rates'])
    install_primitive_precision()
    current=np.array([m.native(t+c*h,g) for c,g in zip(joint.C,gas)])
    maps=[m.jacobian(t+c*h,g) for c,g in zip(joint.C,guides)]
    probes=[]
    for j,(c,g) in enumerate(zip(joint.C,gas)):
        base=current[j];den=np.maximum(np.sum(abs(base)*m.units,axis=0),LD('1e-290'))
        probes.append([(np.sum(abs(m.native(t+c*h,g,p)-base)*m.units,axis=0)/den).astype(float).tolist() for p in [.5,2.]])
    norm=np.linalg.norm(z['defect'])/LD(read(OLD/'rejected-joint-stage.json')['equations'][-1]['relative'])
    with localcontext() as context:
        context.prec=80
        products=[baryon_rows(J,g) for (J,_),g in zip(maps,gas)]
        oldproducts=[baryon_rows(J,g) for (J,_),g in zip(maps,guides)]
        predicted=np.array([[LD(str(exact(base[k,2])+products[j][k]-oldproducts[j][k])) for k in range(m.n)] for j,(_,base) in enumerate(maps)])
        defect=np.array([[LD(str(exact(gas[i][k,2])-exact(initial[k,2])-exact(h)*sum(exact(joint.A[i,j])*exact(current[j,k,2]) for j in range(2)))) for k in range(m.n)] for i in range(2)])
    measure=lambda a:float(np.linalg.norm(a)/norm)
    result=dict(classification='Counterexample candidate',old_native_replay_exact=True,
        original_actual_B_relative=old['exact_outer_stage_sum_B_relative'],promoted_actual_B_relative=measure(defect),
        native_arithmetic_stage_change=measure(h*joint.A@(current[:,:,2]-original[:,:,2])),
        promoted_native_minus_selected_affine_stage_effect=measure(h*joint.A@(current[:,:,2]-predicted)),
        constitutive_relative=probes,constitutive_passed=bool(np.max(probes)<.002),
        new_physical_steps=0,actual_stage_accepted=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'primitive-comparison.npz',original_native=original,promoted_native=current,
        promoted_B_defect=defect,selected_B_prediction=predicted)
    write(OUT/'primitive-result.json',result);print(json.dumps(result),flush=True)
    assert result['constitutive_passed']


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():
                actual=OUT/'affine-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(actual)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
