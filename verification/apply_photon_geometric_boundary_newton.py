"""Retry the one-return photon-boundary pair with a practical Newton budget.

Counterexample candidate. The first attempt (apply_photon_geometric_boundary,
native-photon-boundary265-work) stopped at coarse step 3. Its vector residual
passed (3.8e-16) and its physical moments passed, but the separate phase-257
per-component material gate held the stage: H relative 1.0039, 1.0030 and
1.0023e-13 after the three Newton proposals, against 1e-13. The same step
passed at 7.46e-14 in phase 259. Following the 2026-09-24 instruction on
generous iteration budgets (precedent resume_full_native_newton, 8 proposals),
this retry allows 12 Newton proposals per stage. Every accuracy gate, including
the 1e-13 material gate, is unchanged, and a stage that converges within three
proposals is computed exactly as before. The first attempt, its rejected stage
and its receipts are preserved; nothing is reused from it.
"""
from pathlib import Path
import inspect,resource,sys,time
import apply_photon_geometric_boundary as first
import finish_returned_material_accuracy as material

NEWTON=12
FIRST=first.OUT
# The first retry directory (native-photon-boundary-newton265-work) stopped before any physical step:
# its stage source named NEWTON, which the stage's own globals do not contain (NameError). Preserved.
SECOND=Path('native-photon-boundary-newton265-work')
OUT=first.OUT=Path('native-photon-boundary-newton265b-work')
CHARGE=first.CHARGE=Path('native-photon-boundary-newton-charge265b-work')
EXTERIOR=first.EXTERIOR=Path('native-photon-boundary-newton-exterior265b-work')
read,write,sha=first.read,first.write,first.sha
CAPS=first.CAPS


def patch():
    """Replace only the stage Newton count inside the phase-257 initializer's stage source.

    The count is written as a literal: the stage source is executed in the stage's own
    globals, not in this initializer's.
    """
    if getattr(material.initialize,'newton_budget',None)==NEWTON:return
    source=inspect.getsource(material.initialize)
    old="source=(OUT/'expanded-full-stages.py').read_text()"
    new=old+(f";source=base.base.replace(base.base.replace(source,'    for newton in range(3):','    for newton in range({NEWTON}):'),"
             f"'if newton==2:','if newton=={NEWTON-1}:')")
    assert source.count(old)==1
    ns=dict(material.initialize.__globals__);exec(compile(source.replace(old,new),material.__file__,'exec'),ns)
    fn=ns['initialize'];fn.newton_budget=NEWTON;material.initialize=fn


patch()


def prepare():
    failed=read(FIRST/'coarse-receipt.json');assert 'True native joint Radau equation' in failed['error']
    rejected=read(FIRST/'rejected-joint-stage.json');status=read(FIRST/'pipeline-status.json')
    assert status['state']=='failed' and not (FIRST/'result.json').exists()
    material_last=[max(v['material_relative']) for v in rejected['equations']]
    assert len(material_last)==3 and min(material_last)>1e-13
    broken=[read(SECOND/f'{a}-receipt.json') for a in ['coarse','fine']]
    assert all('NameError' in v['error'] for v in broken) and not (SECOND/'sweep-1/photons/interval-1-64.npz').exists()
    first.prepare();plan=read(OUT/'plan.json')
    evidence=[FIRST/n for n in ['coarse-receipt.json','rejected-joint-stage.json','rejected-joint-stage.npz','pipeline-status.json','plan.json','boundary-check.json']]
    evidence += [SECOND/n for n in ['coarse-receipt.json','fine-receipt.json','pipeline-status.json','plan.json']]
    plan['bindings'].update({str(p):sha(p) for p in evidence+[Path(__file__)]})
    plan.update(claim=plan['claim']+' Retry with a12-proposal Newton budget after the first attempt stopped at coarse step3.',
        first_attempt=dict(directory=str(FIRST),error=failed['error'][:400],material_relative_per_proposal=material_last,
            time=rejected['time'],step=rejected['step'],phase259_same_step_material_relative=7.461873222732249e-14),
        second_attempt=dict(directory=str(SECOND),error=broken[0]['error'],physical_steps=0,
            cause='The budget was written as a name that the stage source, executed in its own globals, could not resolve. Implementation error before any physical step; now a literal count.'),
        resource_change='Allow12Newton proposals per physical stage instead of3 (2026-09-24user instruction; precedent resume_full_native_newton used8). All accuracy gates unchanged, including the separate1e-13material gate; stages converging within3proposals are identical.',
        newton_proposals=NEWTON,scientific_gates_changed=False)
    write(OUT/'plan.json',plan)


def audit():
    first.audit();r=read(OUT/'result.json')
    r.update(newton_proposals=NEWTON,first_attempt_failure_preserved=True,first_attempt_directory=str(FIRST),
        second_attempt_failure_preserved=True,second_attempt_directory=str(SECOND))
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
        elif action=='audit':audit()
        elif action in ['coarse','fine']:first.evolve(64 if action=='coarse' else 128)
        else:getattr(first,action)()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),
            first_source_sha256=sha(first.__file__),newton_proposals=NEWTON,peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
