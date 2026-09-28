"""Apply true momentum arithmetic to the remaining original coarse stage."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import repair_true_momentum_stage as repair

prior=repair.prior;OUT=Path('native-true-momentum238-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
joint,precision,flux=repair.joint,repair.precision,repair.flux
owner=prior.prior.prior.owner;legacy=prior.prior.prior.old
B_NATIVE,B_RHS,B_DEFECT=flux.precise_native,flux.precise_rhs,flux.precise_defect
momentum_product=precision.compile_function(precision.baryon_product,[('J[2::4]','J[3::4]')],{})
CAPS=dict(prepare=180,check=300,coarse=10800,prefix=3600)


def precise_native(m,t,g,native):
    native=B_NATIVE(m,t,g,native)
    with precision.mp.workdps(80):
        rate,raw,face,gravity=repair.native_S(m,t,g,details=True)
        native[0][:,3]=precision.cast(rate);native[1][1]=precision.cast(raw)
        native[3][1]=precision.cast(face);native[4][1]=precision.cast(gravity)
        m.true_S_values[t]=(g.copy(),rate)
    return native


def precise_rhs(m,t,h,v,maps,guides,rhs):
    rhs=B_RHS(m,t,h,v,maps,guides,rhs);m.true_S_values={}
    with precision.mp.workdps(80):
        affine=[]
        for c,g,(J,_) in zip(joint.C,guides,maps):
            now=t+c*h;local=m.local(now);source=m.gas(local['q'],local['qb'],local['qe'])[:,3]
            affine.append(repair.native_S(m,now,g)-momentum_product(J,g)+precision.hp(source))
        initial=m.unpack(v)[1][:,3]
        value=precision.hp(initial)+precision.hp(h)*(precision.hp(joint.A)@np.array(affine))
        for row,new in zip(rhs.reshape(2,-1),value):m.unpack(row)[1][:,3]=precision.cast(new)
    return rhs


def precise_defect(m,t,h,v,sol,defect):
    defect=B_DEFECT(m,t,h,v,sol,defect);pairs=[m.unpack(row) for row in sol.reshape(2,-1)]
    with precision.mp.workdps(80):
        rates=[]
        for c,(x,g) in zip(joint.C,pairs):
            now=t+c*h;saved,rate=m.true_S_values[now];assert np.array_equal(saved,g)
            collision=m.collision(m.local(now),x,g,True)[1][:,3];rates.append(rate+precision.hp(collision))
        initial=m.unpack(v)[1][:,3];end=np.array([g[:,3] for _,g in pairs])
        value=precision.hp(end)-precision.hp(initial)-precision.hp(h)*(precision.hp(joint.A)@np.array(rates))
        for row,new in zip(defect.reshape(2,-1),value):m.unpack(row)[1][:,3]=precision.cast(new)
    return defect


def prepare():
    r=read(repair.OUT/'result.json');assert r['controls_passed'] and r['native_S_arithmetic_effect']>1e-12
    assert not r['saved_proposal_actual_stage_passed'];assert not OUT.exists();OUT.mkdir();files=[]
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/interval-15-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    for src in inputs:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    (OUT/'sweep-1/material').mkdir(exist_ok=True)
    files += [OLD/n for n in ['last-accepted-64.npz','rejected-joint-stage.npz','rejected-joint-stage.json','failure-64.json','coarse-receipt.json','linear-64.json']]
    files += [repair.OUT/n for n in ['plan.json','result.json','comparison.npz','check-receipt.json','symbolic.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Finish actual119with the demonstrated true-native momentum arithmetic correction at all original gates.',
        evidence='235accepted118but119failed1.15246e-11.237reproduced the full defect exactly. Promoting accumulation alone gives1.22608e-11; consistent high-precision S face/gravity evaluation changes it1.19840e-11.40/80digits agree and original0.5/2constitutive probes are below5.02e-23. The saved proposal remains rejected at3.06712e-12 and requires a new RHS/solve.',
        method='Restore118state and every history exactly, not117. Use235last proposal as a guide only. Apply the same true S evaluation to Newton affine RHS and independent nonlinear residual, with consistent face/gravity and material momentum ledgers. Keep235thermal/B repair and210exact linear S operator. Recompute every RHS; no inherited linear acceptance.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,virtual_GiB=12,max_Newton=8,max_linear=12,
        forecast='235two stages1041.49s with final8failed proposals.237one-pair controls are measured in its receipt. Allow3hours for ONE remaining original stage and1hour for saved-prefix ledger verification; no additional path, grid, period or accepted-history replay.',
        stop='Exact restart or original native/linear/nonlinear/physical/constitutive/balance gate,8Newton/12linear or generous wall cap. Save the accepted state and failed proposals.232and236stay frozen. Full paired-time admission still needs explicit common-arithmetic cross-evaluation.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(repair.OUT/'symbolic.json'))


def initialize(seed=True):
    flux.precise_native,flux.precise_rhs,flux.precise_defect=precise_native,precise_rhs,precise_defect
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(seed)
    run=owner.Model.run;source=inspect.getsource(legacy.restore_substep)
    for a,b in [("int(q['macro_index'])==61 and int(q['sub_index'])==1","int(q['macro_index'])==63 and int(q['sub_index'])==1"),
                ("trace['solves'][:2]","trace['solves'][:len(checks['newton'][-1])]"),
                ("previous['resume_macro']=61;previous['resume_substeps']=1","previous['resume_macro']=63;previous['resume_substeps']=1")]:
        assert source.count(a)==1;source=source.replace(a,b)
    ns=dict(legacy.restore_substep.__globals__,OUT=OUT,OLD=OLD,restore_previous=run.__globals__['restore_substep'])
    exec(compile(source,__file__,'exec'),ns)
    owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,restore_substep=ns['restore_substep']),argdefs=run.__defaults__)
    constructor=owner.Model.__init__
    def construct(m,n):
        constructor(m,n);m.true_S_values={}
        if seed and n==64:
            z=np.load(OLD/'rejected-joint-stage.npz');m.resume_seed=dict(time=z['time'][()],solution=z['solution'].copy(),equations=[])
    owner.Model.__init__=construct


def check():
    source=inspect.getsource(legacy.check)
    for a,b in [("m.run(64,'resume-check-64',61,","m.run(64,'resume-check-64',63,"),('accepted_actual_steps=114','accepted_actual_steps=118'),('remaining_coarse_substeps=5','remaining_coarse_substeps=1')]:
        assert source.count(a)==1;source=source.replace(a,b)
    ns=dict(legacy.check.__globals__,OUT=OUT,OLD=OLD,initialize=lambda:initialize(False))
    exec(compile(source,__file__,'exec'),ns);ns['check']()


def factory(log,n):return FunctionType(prior.factory.__code__,dict(prior.factory.__globals__,OUT=OUT))(log,n)
def coarse():FunctionType(prior.coarse.__code__,dict(prior.coarse.__globals__,OUT=OUT,initialize=initialize,factory=factory))()


def prefix():
    source=inspect.getsource(flux.prefix)
    old="g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])"
    assert source.count(old)==1;source=source.replace(old,'g=restored_gas(m,q)')
    ns=dict(flux.prefix.__globals__,OUT=OUT,OLD=OLD,initialize=lambda seed=False:initialize(False),
            precise_native=precise_native,restored_gas=prior.prior.bridge.recovery.prior.restored_gas)
    exec(compile(source,__file__,'exec'),ns);ns['prefix']()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action=='coarse':assert read(OUT/'restart-check.json')['passed']
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
