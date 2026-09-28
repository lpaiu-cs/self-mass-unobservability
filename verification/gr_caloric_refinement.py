"""Continue the unchanged caloric algorithm after failed finite time refinement."""
import json,shutil,sys
import gr_caloric_increment as base
import verify_gr_caloric_increment as auditor

g=base.g;INITIAL=base.OUT;OUT=g.OUT/'gr-caloric-refinement'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    auditor.verify()
    failure=json.loads((auditor.OUT/'path-2.json').read_text())
    assert failure['passed'] and not failure['finite_time_refinement']['passed']
    plan=json.loads((INITIAL/'plan.json').read_text())
    for rel,digest in {**plan['bindings'],**plan['diagnostic_binding']}.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    plan.update(checkpoint='ed0f495',step_counts=[8,16,32],processes=8,
        continuation='Use the unchanged direct-EOS caloric driver, duration, physical operator, Newton method and all original tolerances. The original 1/2 time refinement failed. Compute three finer paths independently from the same initial state; preserve their passes and failures separately.',
        inherited_algorithm_sha256=g.c.sha(g.ROOT/'verification/gr_caloric_increment.py'),
        full_GR_evolution=False,physical_EOS_certified=False)
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_caloric_refinement.py',INITIAL/'plan.json',auditor.OUT/'path-2.json']})
    (OUT/'plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')
    base.OUT=OUT;auditor.OUT=OUT/'exact-audit';auditor.prepare()
    base.symbolic()
    print('PREPARED unchanged caloric paths 8/16/32; original time gate 1e-4 retained',flush=True)


def run():
    base.OUT=OUT;base.run()
    for steps in [8,16,32]:audit(steps)


def audit(steps):
    base.OUT=OUT;auditor.OUT=OUT/'exact-audit'
    if steps==8:
        previous=INITIAL/'path-4.npz';record=INITIAL/'path-4.json'
        assert json.loads(record.read_text())['completed']
        target=OUT/'path-4.npz'
        if not target.exists():shutil.copy2(previous,target)
        assert g.c.sha(target)==g.c.sha(previous)
        link=OUT/'imported-path-4.json'
        binding=dict(classification='Imported from prior work',not_recomputed=True,
            bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [previous,record,target]})
        if link.exists():assert json.loads(link.read_text())==binding
        else:link.write_text(json.dumps(binding,indent=2)+'\n')
    auditor.audit(steps);auditor.verify()


if __name__=='__main__':
    if sys.argv[1]=='audit':audit(int(sys.argv[2]))
    else:globals()[sys.argv[1]]()
