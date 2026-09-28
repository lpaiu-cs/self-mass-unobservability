"""Counterexample candidate: apply consistent B/S arithmetic to actual evolution."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import repair_native_momentum_residual as repair

prior=repair.prior;old=prior.prior;legacy=old.prior.prior
OUT=Path('native-momentum-continuation210-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
owner,joint,LD=prior.owner,prior.joint,prior.LD
B_OPERATOR=prior.stable.stable_operator
CAPS=dict(prepare=180,check=300,coarse=7200,fine=10800,audit=180)


def prepare():
    assert not OUT.exists();OUT.mkdir();r=read(repair.OUT/'result.json')
    assert r['exact_saved_rhs_and_residual'] and r['precision_difference']<1e-25
    assert r['momentum_arithmetic_effect']>1e-14 and not r['linear_gates_passed']
    assert read(OLD/'failure-64.json')['actual_accepted_steps']==115
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    files += [OLD/f'sweep-1/photons/interval-15-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    reused={}
    for p in files:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['failed-linear-64.npz','last-accepted-64.npz','failure-64.json','coarse-receipt.json','linear-64.json','stage-progress-64.json','prefix-result.json']]
    files += [repair.OUT/n for n in ['result.json','residual-comparison.npz','check-receipt.json','symbolic.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='9a67c4653',
        claim='Remove the demonstrated momentum linear-operator arithmetic defect in the ACTUAL remaining same coupled solution and test every original physical gate.',
        evidence='202accepted115steps, then failed1.79440e-14linear vector despite physical moments below2.04e-19.209exactly replays RHS/residual.40/80digit S evaluations agree but differ from old arithmetic by5.58723e-14. Saved solution still fails linear4.87414e-14 and actual stage9.95053e-12; it is only a new proposal.',
        repair='Keep202true60digit native B, its reference-energy conversion and full original native acceptance. Evaluate S native preassembled rows plus their original collision contribution at80digits on EVERY Krylov and residual input. Extend existing restricted full-residual polishing to dominant B OR S rows; retain24columns/4passes and require actual full vector and physical improvement.',
        reuse='Restore115accepted states and every history/ledger exactly. Reuse failed202branch guess and its accumulated linear solution only after exact RHS/guess identity. Do not replay the accepted prefix or re-evaluate unchanged physical coefficients.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,max_Newton=8,max_linear=12,
        forecast='202coarse spent4181.74s, mostly in repeated attempts near the wrong S evaluation floor.209setup/replay is measured separately. Correcting that floor may avoid repetition but speed is not promised. Allow2/3hours for4coarse/16fine remainingsteps; no grid,period or path expansion.',
        stop='Any exact restart/seed,original accuracy/physics,8Newton/12linear or generous wall-cap failure. Preserve the complete accepted state before each actual step and any failed proposal.',
        decision='Only complete paired-time and same-solution ledger acceptance admits full-period GR reading. Existing GR-source interpolation, self-GR and final-charge failures remain.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(repair.OUT/'symbolic.json'))


def stable_operator(m,op,audit,consistent=False):
    cv=inspect.getclosurevars(op._CustomLinearOperator__matvec_impl).nonlocals
    local=inspect.getclosurevars(cv['L']).nonlocals
    baryon=B_OPERATOR(m,op,audit,True)
    momentum=repair.momentum_operator(m,baryon,local['cs'],cv['h'],local['Js'])
    get=lambda x:inspect.getclosurevars(x._CustomLinearOperator__matvec_impl).nonlocals['blocks']
    blocks_by_component={2:get(baryon),3:get(momentum)}
    def apply(value):
        assert len(blocks_by_component)==2
        return momentum.matvec(value)
    return joint.LinearOperator(op.shape,apply,dtype=float)


def polish_function():
    source=inspect.getsource(legacy.polish)
    changes=[("blocks=cv['blocks']","blocks=cv['blocks_by_component']"),
        ("target += [(stage,int(k)) for k in np.flatnonzero(abs(gas[:,2])>goal/8)]","target += [(stage,int(k),component) for component in [2,3] for k in np.flatnonzero(abs(gas[:,component])>goal/8)]"),
        ('for stage,cell in target for col,value in blocks[stage][cell]','for stage,cell,component in target for col,value in blocks[component][stage][cell]'),
        ('dict(blocks[stage][cell]).get(col,0)','dict(blocks[component][stage][cell]).get(col,0)'),
        ('for stage,cell in target),default=LD(0))','for stage,cell,component in target),default=LD(0))'),
        ('records.append(dict(before=',"records.append(dict(components=sorted({v[2] for v in target}),before=")]
    for a,b in changes:assert source.count(a)==1,a;source=source.replace(a,b)
    ns=dict(legacy.polish.__globals__,OUT=OUT);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-BS-polish.py').write_text(source);return ns['polish']


def factory(log,n):
    source=inspect.getsource(legacy.factory).replace('first_system=[True]','first_system=[False]')
    ns=dict(vars(legacy),OUT=OUT,OLD=OLD,polish=polish_function());exec(compile(source,__file__,'exec'),ns)
    solve=ns['factory'](log,n);first=[True]
    def seeded(m,op,P,rhs,guess):
        if n==64 and first[0]:
            first[0]=False;z=dict(np.load(OLD/'failed-linear-64.npz'))
            assert np.array_equal(rhs,z['rhs']) and np.array_equal(guess,z['guess']),'Stored linear branch and RHS changed'
            guess=z['solution'].copy();write(OUT/'linear-seed-identity.json',dict(classification='Counterexample candidate',passed=True,exact_RHS_and_branch_guess=True,accepted_as_physical_state=False))
        return solve(m,op,P,rhs,guess)
    return seeded


def initialize(seed=True):
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    run=owner.Model.run;source=inspect.getsource(old.restore_substep)
    changes=[("int(q['macro_index'])==61 and int(q['sub_index'])==1","int(q['macro_index'])==62 and int(q['sub_index'])==0"),
        ("trace['solves'][:2]","trace['solves'][:1]"),("previous['resume_macro']=61;previous['resume_substeps']=1","previous['resume_macro']=62;previous['resume_substeps']=0")]
    for a,b in changes:assert source.count(a)==1,a;source=source.replace(a,b)
    ns=dict(old.restore_substep.__globals__,OUT=OUT,OLD=OLD,restore_previous=run.__globals__['restore_substep']);exec(compile(source,__file__,'exec'),ns)
    owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,restore_substep=ns['restore_substep']),argdefs=run.__defaults__)
    if seed:
        constructor=owner.Model.__init__
        def construct(m,n):
            constructor(m,n)
            if n==64:
                z=np.load(OLD/'failed-linear-64.npz');p=np.load(OLD/'last-accepted-64.npz')
                m.resume_seed=dict(time=p['next_time'][()],solution=z['guess'].copy(),equations=[])
        owner.Model.__init__=construct


def check():
    source=inspect.getsource(old.check).replace("m.run(64,'resume-check-64',61,","m.run(64,'resume-check-64',62,").replace('accepted_actual_steps=114','accepted_actual_steps=115').replace('remaining_coarse_substeps=5','remaining_coarse_substeps=4')
    ns=dict(old.check.__globals__,OUT=OUT,OLD=OLD,initialize=lambda:initialize(False));exec(compile(source,__file__,'exec'),ns);ns['check']()


def evolve(n):
    prior.stable.stable_operator=stable_operator
    source=inspect.getsource(old.evolve);a='factory=FunctionType(prior.factory.__code__,dict(prior.factory.__globals__,OUT=OUT))';assert source.count(a)==1
    source=source.replace(a,'factory=current_factory')
    ns=dict(old.evolve.__globals__,OUT=OUT,initialize=initialize,current_factory=factory);exec(compile(source,__file__,'exec'),ns);ns['evolve'](n)


def audit():
    FunctionType(old.audit.__code__,dict(old.audit.__globals__,OUT=OUT))()
    r=read(OUT/'result.json');r.update(consistent_B_and_S_linear_arithmetic=True,high_precision_actual_native_flux=True)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
