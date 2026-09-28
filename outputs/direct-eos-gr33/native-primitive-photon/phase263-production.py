"""Reuse the unchanged packet/stress owner with the verified primitive solver."""
from pathlib import Path
source=Path('.phase262-production.py').read_text()
changes={
    'import propagate_returned_moments as adapter':'import propagate_primitive_characteristics as adapter',
    'deadline(14400)':'deadline(28800)',
    "assert b.read(out/'pilot.json')['eligible']":"assert b.read(out/'direct-reference/result.json')['passed']\n    assert b.read(out/'method-comparison.json')['passed']",
    'a,t,g=p.SETTINGS[action]; b.Moments=adapter.Moments; p.Metric=b.Metric':'a,t,g=p.SETTINGS[action]; adapter.install()',
    'p.initialize(); m=p.Photons(a,g)':'p.initialize(); m=adapter.Photons(a,g)',
    "row.update(reference_port=port,component=":"row.update(work_identity_independent=False,reference_port=port,component=",
    "passed=all(np.max(v)<.002 for v in errors.values())":"boundary={k:np.asarray([[r['photon_J_source_cm'],r['photon_lapse_particular']] for r in v],p.LD) for k,v in rows.items()}\n        bnorm=np.maximum(np.max(abs(boundary['fine']),axis=0),p.LD('1e-290'))\n        boundary_controls={k:np.asarray(np.max(abs(v-boundary['fine']),axis=0)/bnorm,float).tolist() for k,v in boundary.items() if k!='fine'}\n        passed=all(np.max(v)<.002 for v in list(errors.values())+list(boundary_controls.values()))",
    'passed=passed,controls=errors,':'passed=passed,controls=errors,boundary_controls=boundary_controls,',
}
for old,new in changes.items():
    assert source.count(old)==1,(old,source.count(old));source=source.replace(old,new)
exec(compile(source,__file__,'exec'))
