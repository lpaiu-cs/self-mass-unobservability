"""Tighten the certified eta box without changing inputs or the Gauss rule."""
from fractions import Fraction as F
from pathlib import Path
from types import SimpleNamespace
import inspect,json,runpy,shutil,sys
import gr_electron_density_certificate as original
import verify_electron_density_export as exported

g=original.g;OUT=g.OUT/'gr-electron-density-refined';CACHE=g.CACHE/'electron-density-refined'
ns=runpy.run_path(original.__file__)['run'].__globals__
ns.update(OUT=OUT,CACHE=CACHE)
OVERLAYS={
    'analyze':('radius=I(1)/2**40','radius=I(1)/2**46'),
    'build':("str(g.ROOT/'verification/gr_electron_density_certificate.cpp')","str(OUT/'density-refined.cpp')")}
for name,(old,new) in OVERLAYS.items():
    source=inspect.getsource(ns[name]);assert source.count(old)==1
    exec(compile(source.replace(old,new),__file__,'exec'),ns)


def prepare():
    exported.verify();assert not OUT.exists();OUT.mkdir();CACHE.mkdir(parents=True,exist_ok=True)
    roots=json.loads((exported.OUT/'result.json').read_text())['rows']
    maximum=max(F(row['root_error_from_center_upper']) for row in roots)
    assert maximum<F(1,2**46)
    failure=g.OUT/'gr-polarization-neutral-defined/result.json'
    assert not json.loads(failure.read_text())['budget_passed']
    cpp=(g.ROOT/'verification/gr_electron_density_certificate.cpp').read_text()
    old='MI radius=MI(1)/MI(1099511627776.0)';new='MI radius=MI(1)/MI(70368744177664.0)'
    assert cpp.count(old)==1;(OUT/'density-refined.cpp').write_text(cpp.replace(old,new))
    names=['candidates.npz','targets.json','states.tsv','candidate-result.json','candidate-manifest.json']
    for name in names:shutil.copy2(original.OUT/name,OUT/name)
    p=json.loads((original.OUT/'plan.json').read_text())
    p.update(checkpoint='35988e4',root_half_width_power=-46,processes=4,
        refinement='The previous exact root intervals all lie within each original center +/- 2^-46. Reevaluate the same density partials on this 64-times narrower box with unchanged 128-bit arithmetic, 8192-panel Gauss rule and uniform remainder. No candidate, composition, constant, quadrature tolerance or scientific budget changes.',
        prior_center_error_upper_exact=str(maximum),source_overlays=OVERLAYS,
        cpp_overlay=dict(original=old,replacement=new))
    p['root_certificate']=p['root_certificate'].replace('this bracket','this narrower bracket')
    p['exact_rule']=p['exact_rule'].replace('2^-40','2^-46')
    paths=[g.ROOT/'verification/gr_electron_density_refinement.py',OUT/'density-refined.cpp',
        original.OUT/'manifest.json',exported.OUT/'manifest.json',exported.OUT/'result.json',failure,
        *[OUT/name for name in names]]
    p['bindings'].update({str(path):g.c.sha(path) for path in paths})
    ns['save']('plan.json',p)
    print('PREPARED unchanged inputs with certified 64-times narrower eta boxes',flush=True)


def export_namespace():
    target=runpy.run_path(exported.__file__)['run'].__globals__
    target.update(original=SimpleNamespace(**ns),OUT=g.OUT/'gr-electron-density-refined-export')
    return target


build=ns['build'];controls=ns['controls'];run=ns['run'];verify=ns['verify']
def export_run():export_namespace()['run']()
def export_verify():export_namespace()['verify']()
if __name__=='__main__':globals()[sys.argv[1]]()
