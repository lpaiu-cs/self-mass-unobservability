"""Continue the corrected actual return with frequent true-residual updates.

Reuse canonical states verbatim. Short inner proposals retain every physical
gate; if twelve short refinements exhaust, the original solver starts from
their last proposal. No applied source, equation, clock or reader is changed.
"""
from pathlib import Path
import gc,inspect,os,resource,shutil,sys,time
import numpy as np
import repair_returned_metric_endpoint as live

OUT=Path('native-short-return259-work');OLD=live.OUT;TEST=Path('native-short-krylov259-work')
read,write,sha,bind=live.read,live.write,live.sha,live.bind
base=live.base;joint=live.joint;LD=joint.LD
CAPS=dict(prepare=600,check=900,coarse=21600,fine=28800,audit=900)


def prepare():
    test=read(TEST/'result.json');assert test['passed'] and test['speedup']>2
    assert read(TEST/'check-receipt.json')['error'] is None
    assert not OUT.exists();OUT.mkdir()
    # Copy the overwritable capture first; only immutable canonical files are
    # linked, and the matching canonical row must already have been accepted.
    recovered=OLD/'recovered-64.npz';before=sha(recovered);raw=recovered.read_bytes()
    assert before==sha(recovered)
    (OUT/'prefix-recovered-64.npz').write_bytes(raw);assert sha(OUT/'prefix-recovered-64.npz')==before
    with np.load(OUT/'prefix-recovered-64.npz') as z:count=len(z['times']);times=z['times'].copy()
    found=[]
    for j in range(1,17):
        p=OLD/f'sweep-1/photons/interval-{j}-64.npz'
        if p.exists() and p.with_suffix('.json').exists():
            with np.load(p) as z:
                if np.array_equal(z['joint_stage_times'],times):found.append(j)
    assert len(found)==1;begin=found[0];assert begin<16
    files=list((OLD/'sweep-0').rglob('*.npz'))
    files += [p for part in ['metric','gr'] for p in (OLD/part).iterdir() if p.is_file()]
    files += [OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','symbolic.json','metric-result.json','endpoint-check.json']]
    files += [OLD/f'sweep-1/photons/interval-{j}-64{ext}' for j in range(1,begin+1) for ext in ['.npz','.json']]
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    (OUT/'sweep-1/material').mkdir(parents=True,exist_ok=True)
    for name in ['charge-reader.py','exterior-reader.py']:
        src=OLD/name;text=src.read_text()
        for a,b in [('native-left-metric258-work',str(OUT)),('native-left-metric-charge258-work','native-short-return-charge259-work'),('native-left-metric-exterior258-work','native-short-return-exterior259-work'),('actual258full-period','actual259full-period'),('Only258actual','Only259actual')]:text=text.replace(a,b)
        dst=OUT/name;dst.write_text(text);files += [src,dst]
    files += [OUT/'prefix-recovered-64.npz',TEST/'result.json',TEST/'reconstruction.json',TEST/'check-receipt.json',OLD/'plan.json',OLD/'controller-start.json',Path(__file__)]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    plan=read(OLD/'plan.json');plan.update(classification='Conjectural',
        claim='Complete the same corrected actual pair and its own final mass/charge audits, removing measured inner-solver waste while retaining all physical accuracy gates.',
        repair='At most12short right-preconditioned GMRES80/1 proposals with true long-double residual updates. If exhausted, the original80/20,12-refinement solver starts from the last proposal. No physical acceptance precedes original vector1e-14,physical1e-13,material1e-13,three-Newton,constitutive,conservation,dense or2percenttime gates.',
        evidence=test,
        reuse='Retain every array of the accepted coarse canonical prefix, including histories, ledgers, guides and residual checks. Verify a zero-step restart and replay only its native rates for omitted diagnostic lists, without integrating the prefix. Fine has not started. All sources,metrics,primary histories and original119/231physical clocks are reused.',
        restart_interval=begin,reused_coarse_actual_steps=count//2,
        forecast='Saved identical stage linearization measured '+str(round(test['rows'][0]['seconds'],3))+'s short versus '+str(round(test['rows'][1]['seconds'],3))+'s original. Full-step and later-path speedup is unmeasured: assume30..120min remaining coarse,60..240min fine,50..90min reader. Keep6h/8h physical caps and16GiB perworker.',
        decision='Keep258unchanged until a new accepted canonical prefix catches its live accepted-step count, then terminate only the bound old coarse worker and preserve the intentional supersession. Any new physical failure stops; no automatic finer clock/grid/path or tolerance change.',
        solver='Frequent true-residual updates with unchanged257acceptance; short proposal followed by original solver only when needed.',
        budgets=CAPS,CPU_affinity=4,virtual_GiB=16,scientific_gates_changed=False,
        discrete_equations_changed=False,source_knot_side_corrected=True,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    plan['bindings']={str(p):sha(p) for p in dict.fromkeys(files)};write(OUT/'plan.json',plan)


def initialize():
    Model=bind(live.initialize,OUT=OUT)();stage=Model.run.__globals__['stages']
    full=inspect.getclosurevars(stage.__globals__['solve']).nonlocals['ns']['solve']
    inner=full.__globals__['gmres'];calls=inspect.getclosurevars(inner).nonlocals['log'];solves=[]
    def limited(*args,**kwargs):return inner(*args,**dict(kwargs,maxiter=1))
    short=bind(full,gmres=limited)
    def solve(m,op,P,rhs,guess):
        begin=len(calls);started=time.monotonic();fallback=False
        try:
            try:answer=short(m,op,P,rhs,guess)
            except AssertionError as exc:
                if not exc.args or not isinstance(exc.args[0],tuple) or exc.args[0][0]!='Four-moment linear residual':raise
                tb=exc.__traceback__
                while tb and tb.tb_frame.f_code!=short.__code__:tb=tb.tb_next
                assert tb is not None
                last=tb.tb_frame.f_locals['sol'];fallback=True
                write(OUT/'last-short-exhaustion.json',dict(classification='Counterexample candidate',error=repr(exc),physical_step_accepted=False,original_solver_continues_from_last_proposal=True))
                answer=full(m,op,P,rhs,last)
            return answer
        finally:
            solves.append(dict(begin=begin,end=len(calls),seconds=time.monotonic()-started,original_fallback=fallback))
            write(OUT/'right-calls.json',dict(classification='Counterexample candidate',calls=calls,solves=solves))
    Model.run=bind(Model.run,stages=bind(stage,solve=solve));return Model


def check():
    begin=read(OUT/'plan.json')['restart_interval'];Model=initialize();m=Model(64)
    label=f'interval-{begin}-64';row=m.run(64,'identity-64',begin*4,label);assert row['passed']
    with np.load(OUT/f'sweep-1/photons/{label}.npz') as a,np.load(OUT/'sweep-1/photons/identity-64.npz') as b:
        assert set(a.files)==set(b.files)
        for k in a.files:assert np.array_equal(a[k],b[k]),k
    del m;gc.collect();native=[];branches=[];errors=[];last_count=0;intervals=[]
    for j in range(1,begin+1):
        m=Model(64);p=OUT/f'sweep-1/photons/interval-{j}-64.npz';z=np.load(p)
        for t,q,rate in zip(z['joint_stage_times'][last_count:],z['joint_stage_conserved_scaled'][last_count:],z['joint_native_rates_scaled'][last_count:]):
            g=base.actual.bridge.recovery.prior.restored_gas(m,q);actual=m.native(t,g)*m.units
            error=float(np.max(np.sum(abs(actual-rate),axis=0)/np.maximum(np.sum(abs(rate),axis=0),LD('1e-290'))))
            errors.append(error);assert error<1e-12,(j,t,error)
        last_count=len(z['joint_stage_times']);native+=m.anchor_checks;branches+=m.branch_checks
        intervals.append(read(p.with_suffix('.json')));del m;gc.collect()
    write(OUT/'prefix-seed.json',dict(anchor_checks=native,branch_checks=branches,intervals=intervals))
    write(OUT/'restart-regression.json',dict(classification='Counterexample candidate',passed=True,
        every_saved_array_exact=True,physical_steps=0,reused_actual_steps=last_count//2,max_saved_native_relative=max(errors)))


def coarse():
    assert read(OUT/'restart-regression.json')['passed'];begin=read(OUT/'plan.json')['restart_interval']
    source=inspect.getsource(base.evolve)
    changes=[('interval-14-{n}',f'interval-{begin}-{{n}}'),("recovered=np.load(OLD/f'recovered-{n}.npz')","recovered=np.load(OUT/f'prefix-recovered-{n}.npz')"),
             ("oldrow=read(OLD/f'run-{n}.json')","oldrow=read(OUT/'prefix-seed.json')"),
             ("oldrow['intervals'][:14]",f"oldrow['intervals'][:{begin}]"),('for j in [15,16]:',f'for j in range({begin+1},17):')]
    for a,b in changes:
        assert a in source,a;source=source.replace(a,b)
    ns=dict(base.evolve.__globals__,OUT=OUT,OLD=OUT,initialize=initialize)
    exec(compile(source,__file__,'exec'),ns);ns['evolve'](64)


def fine():bind(live.evolve,OUT=OUT,initialize=initialize)(128)
def audit():bind(live.audit,OUT=OUT)()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            assert not (OUT/f'{action}-receipt.json').exists()
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
