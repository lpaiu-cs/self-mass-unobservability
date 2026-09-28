"""Continue the unchanged caloric operator after the failed 16/32 time gate."""
import json,shutil,sys
import gr_caloric_increment as base
import gr_caloric_resume as resume
import verify_gr_caloric_increment as auditor

g=base.g;PRIOR=g.OUT/'gr-caloric-refinement';OUT=g.OUT/'gr-caloric-64'


def activate():
    base.OUT=OUT;resume.OUT=OUT;auditor.OUT=OUT/'exact-audit'


def prepare():
    auditor.OUT=PRIOR/'exact-audit';auditor.verify()
    record=PRIOR/'exact-audit/path-32.json';r=json.loads(record.read_text())
    assert r['passed'] and not r['finite_time_refinement']['passed']
    assert not OUT.exists();OUT.mkdir()
    old=PRIOR/'path-32.npz';target=OUT/'path-32.npz';shutil.copy2(old,target);assert g.c.sha(old)==g.c.sha(target)
    source,changed,replacements=resume.transformed_source()
    (OUT/'original-run.py').write_text(source);(OUT/'resumed-run.py').write_text(changed)
    p=json.loads((PRIOR/'plan.json').read_text())
    assert p['inherited_algorithm_sha256']==g.c.sha(g.ROOT/'verification/gr_caloric_increment.py')
    p.update(checkpoint='40892f8',step_counts=[64],processes=8,
        continuation='The original exact 16/32 maximum stored lnT difference exceeded the unchanged 1e-4 time gate. Compute 64 steps from the same initial state with the identical duration, EOS, opacity, face operator, Newton method, arithmetic and tolerances. Retain all previous failures.',
        resume='Reuse the already checked reversible endpoint replay overlay. On first execution it restores zero steps. A later resume must exactly replay every saved accepted endpoint; no accepted state is reset or duplicated.',
        source_overlay=replacements,physical_EOS_certified=False,full_GR_evolution=False)
    paths=[g.ROOT/'verification/gr_caloric_64.py',g.ROOT/'verification/gr_caloric_resume.py',
        PRIOR/'plan.json',record,old,target,OUT/'original-run.py',OUT/'resumed-run.py']
    p['bindings'].update({path.relative_to(g.ROOT).as_posix():g.c.sha(path) for path in paths})
    activate();base.save('plan.json',p);base.symbolic();auditor.prepare()
    base.save('imported-path-32.json',dict(classification='Imported from prior work',
        copied_bitwise=True,not_recomputed=True,
        bindings={path.relative_to(g.ROOT).as_posix():g.c.sha(path) for path in [record,old,target]},
        previous_time_gate_passed=False))
    print('PREPARED unchanged 64-step caloric path; failed 16/32 time gate retained',flush=True)


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in {**p['bindings'],**p['diagnostic_binding']}.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return p


def run():
    bindings();activate();source,changed,_=resume.transformed_source()
    assert source==(OUT/'original-run.py').read_text() and changed==(OUT/'resumed-run.py').read_text()
    namespace=dict(base.run.__globals__,OUT=OUT,initialize_path=resume.initialize_path,save_path=resume.save_path)
    exec(compile(changed,str(OUT/'resumed-run.py'),'exec'),namespace)
    namespace['run']();auditor.audit(64);verify()


def verify():
    bindings();activate();auditor.verify()
    r=json.loads((auditor.OUT/'path-64.json').read_text());assert r['passed'] and 'finite_time_refinement' in r
    print('PASS saved 64-step caloric audit; finite time verdict:',r['finite_time_refinement'],flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
