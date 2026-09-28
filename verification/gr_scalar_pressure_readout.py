"""Keep the pressure-root GR state; compare its audited density replay separately."""
import ast,inspect,json,sys
import gr_scalar_native_audit as audit

original=audit.original;g=original.g;OUT=g.OUT/'gr-scalar-pressure-readout'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def source_overlay():
    source=inspect.getsource(original.run);tree=ast.parse(source)
    calls=[node for node in ast.walk(tree) if isinstance(node,ast.Call) and isinstance(node.func,ast.Name)
        and node.func.id=='save' and isinstance(node.args[0],ast.Constant) and node.args[0].value=='plan.json']
    assert len(calls)==1;plan_write=ast.get_source_segment(source,calls[0])
    changes={'assert not OUT.exists();OUT.mkdir()':'assert OUT.exists()',
        plan_write:"plan=json.loads((OUT/'plan.json').read_text())",
        "assert np.array_equal(native[:,2],state['u_W'])":
        "assert json.loads((audit.OUT/'result.json').read_text())['existing_energy_budget_passed']"}
    changed=source
    for old,new in changes.items():assert changed.count(old)==1;changed=changed.replace(old,new)
    reverse=changed
    for old,new in reversed(list(changes.items())):reverse=reverse.replace(new,old,1)
    assert reverse==source
    return source,changed,changes


def prepare():
    assert not OUT.exists();OUT.mkdir()
    for rel,digest in json.loads((audit.OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest
    assert json.loads((audit.OUT/'result.json').read_text())['existing_energy_budget_passed']
    source,changed,changes=source_overlay()
    (OUT/'original-run.py').write_text(source);(OUT/'pressure-run.py').write_text(changed)
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='5090e6f',change='Preserve the original pressure-root GR state and every scalar IVP coefficient/formula/gate. Replace only the failed exact equality with the independently completed, existing 2 erg/g or 32-ULP energy comparison. The two rounded EOS calls remain distinct; no state entry is overwritten.',
        substitutions=changes,original_bitwise_entry_failed=True,original_GR_state_preserved=True)
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_scalar_pressure_readout.py',original.OUT/'manifest.json',audit.OUT/'manifest.json',OUT/'pressure-run.py']})
    save('plan.json',plan)


def run():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    source,changed,_=source_overlay();assert changed==(OUT/'pressure-run.py').read_text()
    namespace=dict(original.run.__globals__,OUT=OUT,save=save,verify=verify,audit=audit)
    exec(compile(changed,str(OUT/'pressure-run.py'),'exec'),namespace);namespace['run']()


def verify():
    from types import FunctionType
    FunctionType(original.verify.__code__,dict(original.verify.__globals__,OUT=OUT))()


if __name__=='__main__':globals()[sys.argv[1]]()
