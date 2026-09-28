"""Finish the photon-boundary pair with a long-double linear fallback, reusing its accepted prefix.

Counterexample candidate. The H-gate run (native-photon-boundary-hgate265-work)
accepted 118/119 coarse steps and was at fine step 219 when the last coarse
step failed its linear solve: after twelve short right-preconditioned proposals
and twelve refinements of the original double-precision solver, the true
long-double residual stayed at 1.25e-12 with photon-number moment 5.6e-9
(gates 1e-14 and 1e-13). The per-refinement true-residual ratio of the double
corrections ranged from 0.01 to 450: the double operator inside GMRES no longer
represents the long-double acceptance operator well enough.

A scratch probe rebuilt that exact system (its first GMRES call reproduced
bitwise) and solved it with flexible GMRES in long double using the same
double-precision preconditioner; the stage then passed its nonlinear gates.

This module adds that solver only as a third fallback, when both existing
double solvers have failed; systems they accept are computed exactly as before.
Hence every step accepted by the H-gate run is also what this solver accepts,
and its canonical 15/16 checkpoints are reused after an exact zero-step restart
and a replay of the saved native rates (the phase-259 restart pattern). Only
the last canonical interval of each clock is evolved again. No acceptance gate
changes: vector 1e-14, physical moments 1e-13, material Etilde/B/S 1e-13 and
H 2e-13 (user-approved), nonlinear stage, native constitutive, balance, port,
paired-time and the reader's dense 1e-12 gate.
"""
from pathlib import Path
import gc,inspect,os,resource,shutil,sys,time
import numpy as np
import apply_photon_geometric_boundary_hgate as hgate
import solve_long_double_fgmres as ld

first=hgate.first;joint=first.joint;material=hgate.material;base=first.base;LD=np.longdouble
HGATE=hgate.OUT
# First long-double attempt: passed its restart checks, then stopped at fine step231 on a JSON write of a
# numpy.bool_ inside the fallback log (preserved; code at commit0ba369426). The bool came from a correction row
# whose vector residual was already below 1e-14 while a physical moment was not, so at least one long-double
# correction had run; the write raised before the refine returned or raised, so that outcome was not recorded.
LD1=Path('native-photon-boundary-ld265-work')
OUT=hgate.OUT=first.OUT=Path('native-photon-boundary-ld2-265-work')
CHARGE=hgate.CHARGE=first.CHARGE=Path('native-photon-boundary-ld2-charge265-work')
EXTERIOR=hgate.EXTERIOR=first.EXTERIOR=Path('native-photon-boundary-ld2-exterior265-work')
read,write,sha,bind=first.read,first.write,first.sha,first.bind
BEGIN=15
CAPS=dict(prepare=1800,boundary=900,check64=3600,check128=3600,coarse=21600,fine=28800,audit=900)
LD_CORRECTIONS=12;LD_RESTART=80;LD_CYCLES=2;LD_INNER=LD('1e-10')
PROBE=Path('.phase265-ld-probe.json')
EVENTS=[]
base_initialize=first.initialize


def failed(exc):
    return bool(exc.args) and isinstance(exc.args[0],tuple) and exc.args[0][0]=='Four-moment linear residual'


def initialize():
    Model=base_initialize();stage=Model.run.__globals__['stages'];original=stage.__globals__['solve']
    full=inspect.getclosurevars(original).nonlocals['full'];linear=full.__globals__['radau'].prior.owner.reuse.LINEAR
    def refine(m,op,P,rhs,start):
        sol=np.asarray(start,LD).copy();rows=[];began=time.monotonic();step=len(m.newton_iterations)+1
        for k in range(LD_CORRECTIONS+1):
            residual=rhs-op.matvec(sol);relative=float(ld.norm(residual)/max(ld.norm(rhs),LD('1e-290')))
            moments=joint.physical_norm(m,residual)/joint.scales(m,rhs,sol);gas=material.gas_relative(m,residual,sol)
            # bool(): a numpy comparison in the chain would otherwise leave a numpy.bool_, which json cannot write
            # (this stopped the first long-double run at the log write of its fine step231 fallback).
            passed=bool(relative<1e-14 and max(moments)<1e-13 and hgate.gas_gate(gas))
            rows.append(dict(correction=k,relative=relative,moments=[float(v) for v in moments],material_relative=[float(v) for v in gas],passed=passed))
            if passed or k==LD_CORRECTIONS:break
            inner=[];sol=sol+ld.fgmres(op.matvec,P.matvec,residual,restart=LD_RESTART,cycles=LD_CYCLES,rtol=LD_INNER,log=inner)
            rows[-1].update(fgmres_iterations=len(inner),fgmres_seconds=inner[-1]['seconds'] if inner else 0.)
        EVENTS.append(dict(classification='Counterexample candidate',actual_step=step,accepted=passed,corrections=len(rows)-1,
            seconds=time.monotonic()-began,rows=rows))
        write(OUT/f'long-double-{first.TAG[0]}.json',dict(classification='Counterexample candidate',events=EVENTS,
            method='Flexible GMRES in long double with the double-precision preconditioner, only after both double solvers failed.'))
        if not passed:raise AssertionError(('Four-moment linear residual',relative,[float(v) for v in moments]))
        linear.append(dict(initial_info=1,corrections=len(rows)-1,extended_residual=relative,long_double_fgmres=True))
        m.max_residual=max(m.max_residual,relative);return sol
    def solve(m,op,P,rhs,guess):
        try:return original(m,op,P,rhs,guess)
        except AssertionError as exc:
            if not failed(exc):raise
            tb=exc.__traceback__;last=None
            while tb:
                if tb.tb_frame.f_code is full.__code__ and 'sol' in tb.tb_frame.f_locals:last=tb.tb_frame.f_locals['sol']
                tb=tb.tb_next
            return refine(m,op,P,rhs,guess if last is None else last)
    Model.run=bind(Model.run,stages=bind(stage,solve=solve));return Model


first.initialize=initialize


def prepare():
    stop=read(HGATE/'coarse-receipt.json');assert "('Four-moment linear residual'" in stop['error']
    assert read(HGATE/'capture-64.json')['actual_steps']==118 and read(HGATE/'pipeline-status.json')['state']=='failed'
    probe=read(PROBE);assert probe['reproduced_history'] and probe['reproduced_true'] and probe['linear_gates_passed'] and probe['stage_passed']
    ld1=read(LD1/'fine-receipt.json');assert 'is not JSON serializable' in ld1['error'] and read(LD1/'pipeline-status.json')['state']=='failed'
    assert all(read(LD1/f'restart-regression-{n}.json')['passed'] for n in [64,128])
    test=read(Path('.phase265-ld-refine-test.json'));assert test['passed'] and test['json_written']
    assert not OUT.exists();OUT.mkdir()
    files=list((HGATE/'sweep-0').rglob('*.npz'))+[p for part in ['gr','metric'] for p in (HGATE/part).iterdir() if p.is_file()]
    files += [HGATE/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','symbolic.json','endpoint-check.json']]
    files += [HGATE/f'sweep-1/photons/interval-{j}-{n}{ext}' for n in [64,128] for j in range(1,BEGIN+1) for ext in ['.npz','.json']]
    for p in files:
        dst=OUT/p.relative_to(HGATE);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    (OUT/'sweep-1/material').mkdir(parents=True,exist_ok=True)
    counts={}
    for n in [64,128]:
        # The live run rewrites recovered-{n}.npz, so the prefix capture is a separate copy, never a link.
        src=HGATE/f'recovered-{n}.npz';dst=OUT/f'prefix-recovered-{n}.npz';shutil.copyfile(src,dst);assert sha(src)==sha(dst)
        with np.load(dst) as r,np.load(OUT/f'sweep-1/photons/interval-{BEGIN}-{n}.npz') as z:
            assert np.array_equal(r['times'],z['joint_stage_times']);counts[str(n)]=len(z['actual_step_edges'])-1
        files += [src,dst]
    for name in ['charge-reader.py','exterior-reader.py']:
        text=(HGATE/name).read_text()
        for a,b in [('native-photon-boundary-hgate-charge265-work',str(CHARGE)),('native-photon-boundary-hgate-exterior265-work',str(EXTERIOR)),
                    ('native-photon-boundary-hgate265-work',str(OUT))]:text=text.replace(a,b)
        assert 'hgate' not in text,name
        (OUT/name).write_text(text);files += [HGATE/name,OUT/name]
    (OUT/'expanded-hgate-initialize.py').write_text(hgate.PATCHED[0])
    evidence=[HGATE/n for n in ['coarse-receipt.json','right-calls-coarse.json','right-calls-fine.json','last-short-exhaustion-coarse.json',
        'last-short-exhaustion-fine.json','pipeline-status.json','plan.json','capture-64.json','capture-128.json','last-pair-64.npz','boundary-check.json']]
    evidence=[p for p in evidence if p.exists()]+[PROBE,Path('.phase265-ld-probe.py'),Path(ld.__file__),Path(__file__),OUT/'expanded-hgate-initialize.py']
    evidence += [LD1/n for n in ['fine-receipt.json','pipeline-status.json','plan.json','capture-64.json','capture-128.json','right-calls-coarse.json','right-calls-fine.json',
        'restart-regression-64.json','restart-regression-128.json']]+[Path('.phase265-ld-refine-test.json'),Path('.phase265-ld-refine-test.py')]
    files += evidence
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    rows=probe['ld_refinement'];plan=read(HGATE/'plan.json')
    plan.update(classification='Conjectural',
        claim='Complete the one-return photon-boundary pair (user-approved H material gate2e-13) and read its own charge, with a long-double flexible GMRES as the last linear fallback.',
        hgate_attempt=dict(directory=str(HGATE),error=stop['error'][:400],coarse_accepted_steps=118,fine_accepted_steps=read(HGATE/'capture-128.json')['actual_steps'],
            double_refinement_ratio_range=[0.0117,450.67]),
        first_long_double_attempt=dict(directory=str(LD1),error=ld1['error'],coarse_accepted_steps=read(LD1/'capture-64.json')['actual_steps'],
            fine_accepted_steps=read(LD1/'capture-128.json')['actual_steps'],
            cause='Implementation error in the fallback log (numpy.bool_ passed to json) at the first fallback, fine step231. The bool came from a correction row whose vector residual was below 1e-14 while a physical moment was not, so at least one long-double correction had run; the write raised before the refine returned or raised, so the outcome of that fallback was not recorded. Its restart checks passed. Fixed with bool(); the production refine was then run on the rebuilt coarse step119 system in a scratch directory, including its log write (the first such test stopped on a KeyError in its own last bookkeeping line after the refine had returned; preserved, harness line fixed and rerun). A scratch replay of the full fine step231 stage runs as a separate diagnostic.'),
        long_double_fallback=dict(method=f'Flexible GMRES({LD_RESTART}) in long double with the existing double-precision preconditioner and iterative refinement on the long-double true residual; at most {LD_CORRECTIONS} corrections of {LD_CYCLES} cycles, inner relative target {float(LD_INNER):.0e}.',
            order='short right-preconditioned proposals, then the original double solver, then this fallback; systems accepted earlier are unchanged.',
            probe=dict(system='coarse step119 of the H-gate run, rebuilt from its saved last pair; first GMRES call reproduced bitwise',
                corrections=len(rows)-1,final=rows[-1],stage_passed=probe['stage_passed'],stage_seconds=probe['stage_seconds'])),
        reuse=f'The fallback never acts on a system the double solvers accept, so every accepted step of the H-gate run is also accepted here. Reuse its canonical intervals1..{BEGIN} of both clocks ({counts} actual steps) after an exact zero-step restart and a replay of all saved native rates; evolve only interval16. All sources, fields, the shifted metric and the primary history are reused unchanged.',
        restart_interval=BEGIN,reused_actual_steps=counts,
        forecast='H-gate run: coarse interval16 took33min for7steps before the failing step, fine late steps cost up to7min each (259: last16fine solves about29min). Expect about35..60min for both clocks in parallel, a few minutes per long-double fallback (probe about2.4min), checks about10min, then readers about35min. Caps: coarse6h, fine8h.',
        decision='Continue the actual pair through its own compact, frozen-exterior, mass and complete photon-boundary readout. Any original gate failure, an H residual at or above2e-13, another material component at or above1e-13 or a failed long-double fallback stops the run and is preserved; no further gate, clock, grid, path or tolerance change without a new decision.',
        budgets=CAPS,CPU_affinities=dict(coarse=4,fine=6),scientific_gates_changed=False,internal_acceptance_changed_for_H=True,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    plan['bindings']={str(p):sha(p) for p in dict.fromkeys(files)}
    write(OUT/'plan.json',plan)


def boundary():
    """The H-gate module's shifted-metric and stage-driver check, in its own process (one initialization)."""
    first.check()


def check(n):
    Model=initialize();m=Model(n);label=f'interval-{BEGIN}-{n}'
    row=m.run(n,f'identity-{n}',BEGIN*n//16,label);assert row['passed']
    with np.load(OUT/f'sweep-1/photons/{label}.npz') as a,np.load(OUT/f'sweep-1/photons/identity-{n}.npz') as b:
        assert set(a.files)==set(b.files)
        for k in a.files:assert np.array_equal(a[k],b[k]),k
        arrays=len(a.files)
    del m;gc.collect();native=[];branches=[];errors=[];last=0;intervals=[]
    for j in range(1,BEGIN+1):
        m=Model(n);p=OUT/f'sweep-1/photons/interval-{j}-{n}.npz'
        with np.load(p) as z:
            for t,q,rate in zip(z['joint_stage_times'][last:],z['joint_stage_conserved_scaled'][last:],z['joint_native_rates_scaled'][last:]):
                g=base.actual.bridge.recovery.prior.restored_gas(m,q);actual=m.native(t,g)*m.units
                error=float(np.max(np.sum(abs(actual-rate),axis=0)/np.maximum(np.sum(abs(rate),axis=0),LD('1e-290'))))
                errors.append(error);assert error<1e-12,(j,t,error)
            last=len(z['joint_stage_times'])
        native+=m.anchor_checks;branches+=m.branch_checks;intervals.append(read(p.with_suffix('.json')));del m;gc.collect()
    write(OUT/f'prefix-seed-{n}.json',dict(anchor_checks=native,branch_checks=branches,intervals=intervals))
    write(OUT/f'restart-regression-{n}.json',dict(classification='Counterexample candidate',passed=True,every_saved_array_exact=True,
        arrays=arrays,physical_steps=0,reused_stage_count=last,replayed_stage_native_rates=len(errors),max_saved_native_relative=max(errors),
        scope='Exact zero-step restart of the H-gate canonical interval and replay of every saved native rate; no new physical step.'))


def evolve(n):
    assert read(OUT/f'restart-regression-{n}.json')['passed']
    source=inspect.getsource(base.evolve)
    changes=[('interval-14-{n}',f'interval-{BEGIN}-{{n}}'),("recovered=np.load(OLD/f'recovered-{n}.npz')","recovered=np.load(OUT/f'prefix-recovered-{n}.npz')"),
             ("oldrow=read(OLD/f'run-{n}.json')","oldrow=read(OUT/f'prefix-seed-{n}.json')"),
             ("oldrow['intervals'][:14]",f"oldrow['intervals'][:{BEGIN}]"),('for j in [15,16]:',f'for j in range({BEGIN+1},17):')]
    for a,b in changes:
        assert a in source,a;source=source.replace(a,b)
    ns=dict(base.evolve.__globals__,OUT=OUT,OLD=OUT,initialize=initialize)
    exec(compile(source,__file__,'exec'),ns);(OUT/f'expanded-restart-evolve-{n}.py').write_text(source);ns['evolve'](n)
    row=read(OUT/f'run-{n}.json');row.update(full_declared_period=True,source_knot_side_corrected=True,restart_interval=BEGIN,prefix_directory=str(HGATE))
    write(OUT/f'run-{n}.json',row)


def audit():
    hgate.audit();r=read(OUT/'result.json');events={}
    for tag in ['coarse','fine']:
        p=OUT/f'long-double-{tag}.json';events[tag]=read(p)['events'] if p.exists() else []
    assert all(e['accepted'] for v in events.values() for e in v)
    r.update(long_double_fallback_events={k:[dict(actual_step=e['actual_step'],accepted=e['accepted'],corrections=e['corrections'],
        seconds=e['seconds'],final=e['rows'][-1]) for e in v] for k,v in events.items()},restart_interval=BEGIN,
        hgate_prefix_reused=True,hgate_attempt_preserved=True,hgate_attempt_directory=str(HGATE),
        first_long_double_attempt_preserved=True,first_long_double_attempt_directory=str(LD1))
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None;first.TAG[0]=action
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    first.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            assert not (OUT/f'{action}-receipt.json').exists()
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action=='prepare':prepare()
        elif action=='boundary':boundary()
        elif action in ['check64','check128']:check(int(action[5:]))
        elif action=='audit':audit()
        elif action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:raise SystemExit(action)
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),
            hgate_source_sha256=sha(hgate.__file__),long_double_source_sha256=sha(ld.__file__),long_double_events=len(EVENTS),
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
