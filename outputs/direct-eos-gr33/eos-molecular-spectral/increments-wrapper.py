def increments():
    assert not THERMO.exists();THERMO.mkdir();plan=json.loads((thermo.original.OUT/'plan.json').read_text())
    _,changed,replacements=thermo.overlay();(THERMO/'candidate-run.py').write_text(changed)
    plan.update(checkpoint='432ffb0',substitutions=replacements,
        intervention='The existing material increment experiment with the new H2/H2+ spectral EOS and declared chemical-energy anchors, retained molecules. All 55 paths, analytic controls, endpoint and quadrature gates are unchanged.',
        physical_scope=json.loads((OUT/'plan.json').read_text())['boundary'])
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/eos_molecular_spectral.py',OUT/'plan.json',OUT/'controls.json',THERMO/'candidate-run.py']})
    (THERMO/'plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')
    def save_thermo(name,value):(THERMO/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    def verify_thermo():FunctionType(thermo.original.verify.__code__,dict(thermo.original.verify.__globals__,OUT=THERMO))()
    namespace=dict(thermo.original.run.__globals__,OUT=THERMO,save=save_thermo,verify=verify_thermo,Continuation=EOS)
    exec(compile(changed,str(THERMO/'candidate-run.py'),'exec'),namespace);namespace['run']()
