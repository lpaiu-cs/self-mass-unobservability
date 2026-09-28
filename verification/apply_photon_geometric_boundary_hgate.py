"""Evolve the one-return photon-boundary pair with the user-approved H material gate.

Counterexample candidate. Three earlier attempts are preserved:
native-photon-boundary265-work stopped at coarse step 3 on the phase-257
per-component material gate (H 1.0039/1.0030/1.0023e-13 against 1e-13 after
three Newton proposals); native-photon-boundary-newton265-work stopped before
any physical step on an implementation error (NameError); and
native-photon-boundary-newton265b-work kept H at 1.00274..1.00282e-13 over
twelve proposals while the vector residual fell to 5.4e-17.

On 2026-09-26 the user approved setting the internal H gate to 2e-13 (option
1 of the report). Only that component changes, in both places where the
phase-257 gate acts: the linear proposal and the true nonlinear stage
acceptance. Etilde, B and S stay at 1e-13. The original vector,
physical-moment, native constitutive, shared-face balance, port, 2 percent
paired-time and the reader's dense 1e-12 source gates are unchanged. The
twelve-proposal Newton budget of the previous retry is kept; a stage that
converges within three proposals is computed as before.
"""
from pathlib import Path
from collections import Counter
import inspect,resource,sys,time
import apply_photon_geometric_boundary as first
import finish_returned_material_accuracy as material

NEWTON=12
GAS_GATES=(1e-13,2e-13,1e-13,1e-13)  # Etilde, H (neutral), B, S; the same for both Radau stages
H=1
ATTEMPTS=dict(first=first.OUT,name_error=Path('native-photon-boundary-newton265-work'),newton12=Path('native-photon-boundary-newton265b-work'))
OUT=first.OUT=Path('native-photon-boundary-hgate265-work')
CHARGE=first.CHARGE=Path('native-photon-boundary-hgate-charge265-work')
EXTERIOR=first.EXTERIOR=Path('native-photon-boundary-hgate-exterior265-work')
read,write,sha=first.read,first.write,first.sha
CAPS=first.CAPS
PATCHED=[None]


def gas_gate(values):
    """Per-component material gate on gas_relative(...): [stage0 Etilde,H,B,S, stage1 Etilde,H,B,S]."""
    assert len(values)==2*len(GAS_GATES),values
    return all(v<GAS_GATES[k%len(GAS_GATES)] for k,v in enumerate(values))


def patch():
    """Rewrite only the gate expressions and the Newton count inside the phase-257 initializer."""
    if getattr(material.initialize,'hgate',False):return
    source=inspect.getsource(material.initialize)
    read_stage="source=(OUT/'expanded-full-stages.py').read_text()"
    changes=[
        (read_stage,read_stage+f";source=base.base.replace(base.base.replace(source,'    for newton in range(3):','    for newton in range({NEWTON}):'),'if newton==2:','if newton=={NEWTON-1}:')",1),
        ('max(gas_relative(m,residual,sol))<1e-13','gas_gate(gas_relative(m,residual,sol))',1),
        ('max(gas_relative(m,defect,sol))<1e-13','gas_gate(gas_relative(m,defect,sol))',1),
        ('gas_relative=gas_relative)','gas_relative=gas_relative,gas_gate=gas_gate)',2)]
    for a,b,count in changes:
        assert source.count(a)==count,(a,source.count(a));source=source.replace(a,b)
    ns=dict(material.initialize.__globals__,gas_gate=gas_gate)
    exec(compile(source,f'{material.__file__}#hgate-initialize','exec'),ns)
    fn=ns['initialize'];fn.hgate=True;material.initialize=fn;PATCHED[0]=source


patch()


def worst(equations):
    rows=[e['material_relative'] for e in equations]
    return [max(v[H],v[H+4]) for v in rows],[max(x for k,x in enumerate(v) if k%4!=H) for v in rows]


def prepare():
    stop=read(ATTEMPTS['first']/'coarse-receipt.json');assert 'True native joint Radau equation' in stop['error']
    broken=[read(ATTEMPTS['name_error']/f'{a}-receipt.json') for a in ['coarse','fine']]
    assert all('NameError' in v['error'] for v in broken)
    retry=read(ATTEMPTS['newton12']/'coarse-receipt.json');assert 'True native joint Radau equation' in retry['error']
    rejected={k:read(ATTEMPTS[k]/'rejected-joint-stage.json') for k in ['first','newton12']}
    hvalues={k:worst(v['equations'])[0] for k,v in rejected.items()};others=[x for v in rejected.values() for x in worst(v['equations'])[1]]
    assert len(hvalues['first'])==3 and len(hvalues['newton12'])==NEWTON
    assert min(min(v) for v in hvalues.values())>1e-13 and max(max(v) for v in hvalues.values())<GAS_GATES[H] and max(others)<1e-13
    first.prepare();plan=read(OUT/'plan.json')
    (OUT/'expanded-hgate-initialize.py').write_text(PATCHED[0])
    candidates=[ATTEMPTS['first']/n for n in ['coarse-receipt.json','rejected-joint-stage.json','rejected-joint-stage.npz','pipeline-status.json','plan.json','boundary-check.json']]
    candidates += [ATTEMPTS['name_error']/n for n in ['coarse-receipt.json','fine-receipt.json','pipeline-status.json','plan.json']]
    candidates += [ATTEMPTS['newton12']/n for n in ['coarse-receipt.json','rejected-joint-stage.json','rejected-joint-stage.npz','pipeline-status.json','plan.json']]
    evidence=[p for p in candidates if p.exists()]
    plan['bindings'].update({str(p):sha(p) for p in evidence+[Path(__file__),OUT/'expanded-hgate-initialize.py']})
    plan.update(claim=plan['claim']+' Internal H material gate2e-13 by user approval; all other gates unchanged.',
        attempts={k:str(v) for k,v in ATTEMPTS.items()},
        attempt_H_material_relative={k:v for k,v in hvalues.items()},
        attempt_errors=dict(first=stop['error'][:300],name_error=broken[0]['error'],newton12=retry['error'][:300]),
        gate_change=dict(component='H (neutral) per-component material residual of the phase-257 internal gate',old=1e-13,new=GAS_GATES[H],
            applies='linear proposal and true nonlinear stage acceptance',
            unchanged='Etilde/B/S material1e-13; vector1e-14 linear and1e-12 nonlinear; physical moments1e-13; native constitutive0.2percent; shared-face balance1e-8; port1e-12; paired time2percent; reader dense source1e-12',
            reason='With the new boundary the accepted-quality coarse step3 H residual stays at1.0027..1.0039e-13 (a long-double floor of the stiff H rate; twelve proposals do not reduce it) while the vector residual reaches5.4e-17. Phase259 accepted H up to8.79e-14. The reader dense gate keeps a5x margin.',
            user_approval='2026-09-26: the user chose option1 of the report (H internal gate2e-13).'),
        internal_material_gate=dict(zip(['Etilde','H','B','S'],GAS_GATES)),newton_proposals=NEWTON,
        solver='Unchanged259solver (short right-preconditioned GMRES80/1 proposals, original257solver as fallback) with at most12Newton proposals. Gates: vector1e-14 linear, true nonlinear1e-12, physical moments1e-13, material Etilde/B/S1e-13 and H2e-13 (linear and nonlinear), native constitutive0.2percent, balance1e-8, paired time2percent.',
        forecast='265b measured about20s per fine step and about70s for the first coarse steps including initialization; 259measured coarse111steps3615s and fine231steps6056s. Expect coarse about60..80min and fine about80..110min in parallel, then readers about35min (259:1880s). Caps unchanged: coarse6h, fine8h, readers as registered.',
        decision='Continue the actual pair through its own compact, frozen-exterior, mass and complete photon-boundary readout. Any original gate failure, an H residual at or above2e-13 or another material component at or above1e-13 stops the run and is preserved; no further gate, clock, grid, path or tolerance change without a new decision.',
        scientific_gates_changed=False,internal_acceptance_changed_for_H=True)
    write(OUT/'plan.json',plan)


def audit():
    first.audit();r=read(OUT/'result.json');rows={}
    for n in [64,128]:
        steps=read(OUT/f'checks-{n}.json')['newton'];hv,ov=worst([s[-1] for s in steps])
        assert max(hv)<GAS_GATES[H] and max(ov)<1e-13,(n,max(hv),max(ov))
        rows[str(n)]=dict(steps=len(steps),proposals_used={str(k):v for k,v in sorted(Counter(len(s) for s in steps).items())},
            maximum_H=max(hv),steps_with_H_above_1e_13=int(sum(v>=1e-13 for v in hv)),maximum_other_material=max(ov))
    r.update(internal_material_gate=dict(zip(['Etilde','H','B','S'],GAS_GATES)),internal_acceptance_changed_for_H=True,
        H_gate_rows=rows,newton_proposals=NEWTON,attempts={k:str(v) for k,v in ATTEMPTS.items()},attempt_failures_preserved=True)
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
            first_source_sha256=sha(first.__file__),newton_proposals=NEWTON,internal_material_gate=list(GAS_GATES),
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
