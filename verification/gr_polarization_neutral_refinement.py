"""Propagate the narrower certified inverse derivatives through the same bounds."""
from types import SimpleNamespace
import inspect,json,runpy,sys
import gr_polarization_neutral_runner as previous
import gr_electron_density_refined_runner as refined

g=previous.g;OUT=g.OUT/'gr-polarization-neutral-refined'
ns=runpy.run_path(previous.base.__file__)['run'].__globals__
ns.update(OUT=OUT,export=SimpleNamespace(**refined.export_namespace()))
source=inspect.getsource(ns['symbolic']);assert source.count(previous.old)==1
exec(compile(source.replace(previous.old,previous.new),__file__,'exec'),ns)


def prepare():
    refined.verify();previous.verify();assert not json.loads((previous.OUT/'result.json').read_text())['budget_passed']
    ns['prepare']();p=json.loads((OUT/'plan.json').read_text())
    paths=[g.ROOT/'verification/gr_polarization_neutral_refinement.py',
        g.ROOT/'verification/gr_electron_density_refined_runner.py',previous.OUT/'manifest.json',
        previous.OUT/'result.json',refined.OUT/'manifest.json']
    p['bindings'].update({x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in paths})
    p['refinement']='Use the certified narrower eta-box inverse-derivative intervals. Reuse the original numerical center, neutral coefficient candidates, all-Q bounds and 2e-7 field-error budget without alteration.'
    ns['save']('plan.json',p)


run=ns['run'];verify=ns['verify']
if __name__=='__main__':globals()[sys.argv[1]]()
