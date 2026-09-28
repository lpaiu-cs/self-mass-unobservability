"""Counterexample candidate: preserve native precision in actual continuation."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import repair_native_affine_precision as precision

prior=precision.prior;before,stable,owner,joint=prior.before,prior.stable,prior.owner,prior.joint
OUT=Path('native-precise-continuation199-work');OLD=prior.OUT
LD=joint.LD;read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=60,prefix=180,coarse=1800,fine=2700,audit=60)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    p=read(precision.OUT/'primitive-result.json')
    assert p['constitutive_passed'] and p['promoted_native_minus_selected_affine_stage_effect']<1e-12
    assert not p['actual_stage_accepted']
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    files += [OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    reused={}
    for src in files:dst=OUT/src.relative_to(OLD);os.link(src,dst);reused[str(dst)]=sha(src)
    stored=dict(np.load(OLD/'last-accepted-64.npz'));original=dict(np.load(prior.prior.OLD/'last-accepted-64.npz'))
    fingerprint=before.prior.prior.prior.fingerprint
    assert stored.keys()==original.keys()
    for k in stored:assert fingerprint(stored[k])==fingerprint(original[k]),k
    files += [OLD/n for n in ['last-accepted-64.npz','rejected-joint-stage.npz','rejected-joint-stage.json','coarse-receipt.json']]
    files += [precision.OUT/n for n in ['primitive-result.json','primitive-plan.json','check-result-new.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='5cec2d937',
        claim='Apply the identified primitive precision repair to the actual remaining photon/material coupled evolution, then its own energy/boundary and time readouts.',
        evidence='197passed every linear solve but failed the true native stage after8Newton proposals.198replayed it exactly and rejected affine RHS rounding as dominant. Five default-float outputs and five explicit float casts in the native primitive variation caused state-dependent quantization. Keeping them in long double reduced native/selected-affine mismatch5.184e-11to7.567e-13and passed original0.2percent constitutive controls.',
        repair='Promote only native primitive intermediate/output precision in BOTH actual and selected-branch owners. Keep EOS banks and photon coefficient maps unchanged, and preserve the original native/linear/constitutive/ledger/time gates. Do not deploy the affine-only change that failed to explain the error.',
        reuse='Restore the exact113accepted substeps. Before production, re-evaluate every saved gas stage with the promoted native owner and check its integrated four-material ledger and original0.2percent RHS scale control. This is an arithmetic consistency check, not a new integration or a uniform error certificate.',
        proposal='Use the last197rejected Radau pair only as the initial Newton proposal. Allow8new proposals with the repaired native arithmetic; the8failed old-arithmetic proposals remain rejected and separately recorded. No old failed linear system is reused as though its RHS remained identical.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='197cost816.53s, including679.59s in two difficult inner solves. Its accepted linear proposals and final failed native pair are saved. Allow30/45minutes again; no speed guarantee for later steps. The prefix check adds native RHS evaluations, not Jacobian/Krylov work or physical integration.',
        stop='Any prefix ledger/constitutive, actual native,8newNewton/12linear,packet or2percent time failure;180s prefix or30/45minute path caps. No new grid,period or extra parameter path.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,
        exact_all_keys_of_last_accepted_state=True,new_physical_steps=0,prior_original_state=113))


def initialize(seed=True):
    init=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT));init()
    install=FunctionType(precision.install_primitive_precision.__code__,dict(precision.install_primitive_precision.__globals__,OUT=OUT));install()
    stage=owner.Model.run.__globals__['stages'];source=(OUT/'expanded-eight-newton-stage.py').read_text()
    a='limit=5 if continuing else 8';assert source.count(a)==1;source=source.replace(a,'limit=8')
    ns=dict(stage.__globals__);exec(compile(source,__file__,'exec'),ns)
    run=owner.Model.run;owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,stages=ns['stages']),argdefs=run.__defaults__)
    (OUT/'expanded-precise-native-stage.py').write_text(source)
    if seed:
        constructor=owner.Model.__init__
        def construct(m,n):
            constructor(m,n)
            if n==64:
                z=dict(np.load(OLD/'rejected-joint-stage.npz'));z['equations']=[];m.resume_seed=z
        owner.Model.__init__=construct


def factory(log,n):
    source=inspect.getsource(prior.factory).replace('first_system=[True]','first_system=[False]')
    assert source!=inspect.getsource(prior.factory)
    polish=FunctionType(prior.polish.__code__,dict(prior.polish.__globals__,OUT=OUT))
    ns=dict(vars(prior),OUT=OUT,OLD=OLD,polish=polish);exec(compile(source,__file__,'exec'),ns)
    return ns['factory'](log,n)


def prefix():
    initialize(False);rows=[]
    for n in [64,128]:
        m=owner.Model(n)
        file=OLD/'last-accepted-64.npz' if n==64 else OLD/'sweep-1/photons/interval-15-128.npz'
        z=dict(np.load(file));rates=[];changes=[]
        for j,t in enumerate(z['joint_stage_times']):
            q=z['joint_stage_conserved_scaled'][j]
            g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])
            value=m.native(t,g)*m.units;old=z['joint_native_rates_scaled'][j]
            rates.append(value)
            changes.append(np.sum(abs(value-old),axis=0)/np.maximum(np.sum(abs(old),axis=0),LD('1e-290')))
        expected=np.sum(z['joint_stage_weights'][:,None,None].astype(LD)*(np.array(rates)+z['joint_collision_rates_scaled']),axis=0,dtype=LD)
        if n==64:actual=z['g']*m.units+z['material_floor_discard_scaled']
        else:actual=z['restart_g']*m.units+z['material_floor_discard_scaled']
        error=np.sum(abs(actual-expected),axis=0)/np.maximum(np.sum(abs(actual)+abs(expected),axis=0),LD('1e-290'))
        row=dict(clock=n,actual_stage_count=len(rates),native_rate_change=np.max(changes,axis=0).astype(float).tolist(),
            same_prefix_material_balance=error.astype(float).tolist(),passed=bool(max(error)<1e-8 and np.max(changes)<.002))
        rows.append(row);write(OUT/f'prefix-{n}.json',dict(classification='Counterexample candidate',**row));assert row['passed'],row
    write(OUT/'prefix-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        original_local_stage_audits_preserved=True,new_physical_steps=0,uniform_error_certificate=False,final_charge_conclusion='unadjudicated'))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'prefix-result.json')['passed'];stable.initialize=initialize;n=64 if action=='coarse' else 128
            old=SimpleNamespace(right_solver=lambda log:factory(log,n))
            evolve=FunctionType(before.evolve.__code__,dict(before.evolve.__globals__,OUT=OUT,old=old));evolve(n)
        elif action=='audit':
            audit=FunctionType(before.audit.__code__,dict(before.audit.__globals__,OUT=OUT));audit()
            r=read(OUT/'result.json');r.update(maximum_inner_solves=12,maximum_new_Newton_proposals=8,previous_failed_proposals_preserved=True,native_primitive_precision_preserved=True);write(OUT/'result.json',r)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
