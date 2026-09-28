"""Retain the failed 64/128-point self integral and refine the outer rule."""
from types import FunctionType
import json, shutil, sys
import gr_polarization_self_energy as original

g=original.g;OUT=g.OUT/'gr-polarization-self-refinement'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=original.finite.previous.original.previous.scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.bindings()
    for rel,digest in json.loads((original.OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert not json.loads((original.OUT/'result.json').read_text())['passed']
    shutil.copy2(g.ROOT/'outputs/gr-polarization-self-energy33-run.log',OUT/'original-error.log')
    paths=[g.ROOT/'verification/gr_polarization_self_refinement.py',original.OUT/'plan.json',
        original.OUT/'manifest.json',original.OUT/'result.json',OUT/'original-error.log']
    save('plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},gauss_orders=[256,512],
        original_failure='64/128 difference 8.841853974472791e-5 exceeded the fixed 2e-7 gate. Independent adaptive controls at cells4933 and5734 also failed. This is an interior quadrature-resolution failure, not a failure of the separate analytic UV bound.',
        change='Change only the two Gauss-Legendre orders to 256/512. Recompute all states and controls. Equations, inner tolerance, cut, comparison gates and analytic tail/global bounds are unchanged. Require bitwise-identical symbolic and bound records.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((original.OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert not json.loads((original.OUT/'result.json').read_text())['passed']
    return dict(original.bindings(),gauss_orders=plan['gauss_orders'])


def run():
    ns=dict(original.__dict__,OUT=OUT,save=save,bindings=bindings,verify=verify)
    for name in ['symbolic','tail_bounds','run']:ns[name]=FunctionType(original.__dict__[name].__code__,ns)
    ns['run']()


def verify():
    ns=dict(original.__dict__,OUT=OUT,bindings=bindings)
    FunctionType(original.verify.__code__,ns)()
    for name in ['symbolic.json','tail-bounds.json']:assert (original.OUT/name).read_bytes()==(OUT/name).read_bytes(),name
    print('PASS refined outer self integral; failed64/128 and identical analytic bounds retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
