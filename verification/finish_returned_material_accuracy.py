"""Meet the unchanged dense-source gate in the actual coupled solution.

Right preconditioning proposes states; original physical gates plus a stricter
separate material residual decide their acceptance. No readout is patched.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import gc,inspect,os,resource,shutil,sys,time
import numpy as np
import finish_returned_right_residual as prior

OUT=Path('native-material-accuracy257-work');OLD=prior.OUT;base=prior.base;joint=prior.joint
read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind
CAPS=dict(prepare=600,check=600,fine=21600,audit=600)
FAILED=Path('native-right-charge256-work/full/dense')


def prepare():
    assert read(OLD/'result.json')['passed']
    bad=read(FAILED/'source-128-check.json');assert not bad['passed'] and bad['dense_stage_max']>1e-12
    assert read(FAILED/'source-64-check.json')['passed']
    for n in [64,128]:assert max(max(v) for v in read(FAILED/f'source-{n}-check.json')['dense_stage'][:2*215 if n==128 else 2*119])<1e-12
    assert not OUT.exists();OUT.mkdir();files=list((OLD/'sweep-0').rglob('*.npz'))
    files += [p for d in ['metric','gr'] for p in (OLD/d).iterdir() if p.is_file()]
    files += [p for p in (OLD/'sweep-1/photons').iterdir() if p.name.endswith(('-64.npz','-64.json')) or p.name in ['interval-14-128.npz','interval-14-128.json','interval-15-128.npz','interval-15-128.json']]
    files += [OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','symbolic.json','metric-result.json','metric-receipt.json','coarse-receipt.json','run-64.json','checks-64.json','recovered-64.npz','recovered-128.npz']]
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True)
        if p.name=='recovered-128.npz':shutil.copy2(p,dst)
        else:os.link(p,dst)
    (OUT/'sweep-1/material').mkdir(parents=True,exist_ok=True)
    files += [OLD/'result.json',FAILED/'source-64-check.json',FAILED/'source-128-check.json',FAILED/'source-receipt.json',Path('phase257-neutral-precision.json')]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    p=read(OLD/'plan.json');p.update(classification='Conjectural',
        claim='Resolve the actual neutral-material stage defect that blocks the unchanged same-solution dense charge reader.',
        evidence='256original actual119/231 and paired-time gates passed. Its dense reader failed fine220second-stage H residual1.719880872203483e-12 versus1e-12. Direct80digit conserved-coordinate summation gives1.719880873310006e-12; arithmetic effect2.19e-20, so changing summation or the readout is not an adequate repair. All coarse and fine1..215dense checks passed.',
        repair='Use existing right-preconditioned GMRES80/20 with at most12 extended refinements. Keep original vector1e-14, combined physical1e-13, three Newton proposals, stage1e-12, conservation/constitutive/port/time gates. Add separate per-stage gas Etilde/H/B/S residual below1e-13 relative to its own stored state, both to linear and true nonlinear acceptance. The original dense1e-12gate is unchanged.',
        reuse='Exact zero-step restart at the accepted215fine checkpoint; retain full coarse119, the primary, applied metric and earlier fine states. Re-evolve only original216..231 because a correction at220must propagate in the same coupled solution. Preserve256accepted full path and its rejected dense reader; do not add a correction to its final charge.',
        forecast='256remaining16fine steps measured2985.50s with the left-first solver; one saved difficult right solve took1.524s/10iterations. Other systems with the extra gas gate are unmeasured: assume5..60minutes and retain a generous6hour cap,16GiB,CPU3. No new physical grid, clock, period or independent path.',
        budgets=CAPS,internal_material_gate=1e-13,original_dense_source_gate=1e-12,
        scientific_gates_changed=False,internal_acceptance_strengthened=True,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    p['bindings'].update({str(p):sha(p) for p in dict.fromkeys(files)});write(OUT/'plan.json',p)
    write(OUT/'right-symbolic.json',read(OLD/'right-symbolic.json'))


def gas_relative(m,residual,solution):
    rr=np.array([m.unpack(v)[1] for v in residual.reshape(2,-1)])
    ss=np.array([m.unpack(v)[1] for v in solution.reshape(2,-1)])
    return (np.sum(abs(rr)*m.units,axis=1)/np.maximum(np.sum(abs(ss)*m.units,axis=1),joint.LD('1e-290'))).astype(float).ravel().tolist()


def initialize():
    Model=bind(base.initialize,OUT=OUT)();calls=[]
    template=prior.right.right_solver(calls);original=template.__globals__['gmres']
    inner=FunctionType(original.__code__,dict(original.__globals__,joint=SimpleNamespace(gmres=prior.scipy_gmres)),argdefs=original.__defaults__,closure=original.__closure__)
    source=inspect.getsource(joint.solve)
    for a,b in [('restart=20,maxiter=5','restart=80,maxiter=20'),('range(4)','range(12)'),('assert k<3,','assert k<11,'),
                ('if relative<1e-14 and max(moments)<1e-13:','if relative<1e-14 and max(moments)<1e-13 and max(gas_relative(m,residual,sol))<1e-13:')]:
        assert source.count(a)==1;source=source.replace(a,b)
    ns=dict(joint.solve.__globals__,gmres=inner,gas_relative=gas_relative);exec(compile(source,__file__,'exec'),ns)
    def solve(m,op,P,rhs,guess):
        start=len(calls)
        try:return ns['solve'](m,op,P,rhs,guess)
        except BaseException as exc:
            tb=exc.__traceback__
            while tb and not {'sol','residual','rhs'}.issubset(tb.tb_frame.f_locals):tb=tb.tb_next
            if tb:
                v=tb.tb_frame.f_locals;np.savez_compressed(OUT/'failed-linear.npz',rhs=rhs,guess=guess,solution=v['sol'],residual=v['residual'])
                write(OUT/'failed-linear.json',dict(error=repr(exc),gas_relative=gas_relative(m,v['residual'],v['sol'])))
            raise
        finally:write(OUT/'right-calls.json',dict(classification='Counterexample candidate',last_solve_begin=start,calls=calls))
    stage=Model.run.__globals__['stages'];source=(OUT/'expanded-full-stages.py').read_text()
    old='if relative<1e-12 and max(moments)<1e-13:break'
    source=base.base.replace(source,old,'if relative<1e-12 and max(moments)<1e-13 and max(gas_relative(m,defect,sol))<1e-13:break')
    source=base.base.replace(source,'audit.append(dict(relative=relative,moments=moments.astype(float).tolist()))','audit.append(dict(relative=relative,moments=moments.astype(float).tolist(),material_relative=gas_relative(m,defect,sol)))')
    ns_stage=dict(stage.__globals__,solve=solve,gas_relative=gas_relative);exec(compile(source,__file__,'exec'),ns_stage)
    Model.run=bind(Model.run,stages=ns_stage['stages']);return Model


def check():
    Model=initialize();m=Model(128);row=m.run(128,'identity-128',120,'interval-15-128');assert row['passed']
    with np.load(OUT/'sweep-1/photons/interval-15-128.npz') as a,np.load(OUT/'sweep-1/photons/identity-128.npz') as b:
        assert set(a.files)==set(b.files)
        for k in a.files:assert np.array_equal(a[k],b[k]),k
    write(OUT/'restart-regression.json',dict(classification='Counterexample candidate',passed=True,every_saved_array_exact=True,physical_steps=0,reused_fine_steps=215))


def fine():
    assert read(OUT/'restart-regression.json')['passed'];Model=initialize()
    seed=bind(prior.seed_row,OUT=OUT,OLD=OLD)(Model);source=inspect.getsource(base.evolve)
    for a,b in [("interval-14-{n}.npz","interval-15-{n}.npz"),("oldrow=read(OLD/f'run-{n}.json')","oldrow=seed"),("oldrow['intervals'][:14]","oldrow['intervals'][:15]"),("restart=f'interval-14-{n}'","restart=f'interval-15-{n}'"),('for j in [15,16]:','for j in [16]:')]:
        assert source.count(a)==1;source=source.replace(a,b)
    ns=dict(base.evolve.__globals__,OUT=OUT,OLD=OLD,initialize=lambda:Model,seed=seed);exec(compile(source,__file__,'exec'),ns);ns['evolve'](128)
    row=read(OUT/'run-128.json');checks=read(OUT/'checks-128.json')
    row['maximum_new_true_material_stage']=max(max(v[-1]['material_relative']) for v in checks['newton'][215:])
    assert row['maximum_new_true_material_stage']<1e-13;write(OUT/'run-128.json',row)


def audit():
    bind(prior.audit,OUT=OUT)();r=read(OUT/'result.json')
    r.update(original256_dense_failure_preserved=True,internal_material_gate=1e-13,original_dense_source_gate=1e-12,right_preconditioning_primary=True)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2);joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
