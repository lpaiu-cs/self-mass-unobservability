"""Preserve an unevaluated temperature slot in the symbolic chain-rule check."""
from types import FunctionType
import inspect, json, shutil, sys
import gr_plasma_photon_model as original

g=original.g;OUT=g.OUT/'gr-plasma-photon-defined'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();plan=json.loads((original.OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert not (original.OUT/'result.json').exists()
    shutil.copy2(g.ROOT/'outputs/gr-plasma-photon33-run.log',original.OUT/'failure.log')
    (original.OUT/'failure.json').write_text(json.dumps(dict(classification='Counterexample candidate',failed=True,
        reason='SymPy Subs simplified a partial temperature derivative after inserting temperature-dependent b. An independent dummy first argument retains the partial-derivative slot; direct chain-rule expression simplifies to zero.'),indent=2)+'\n')
    (original.OUT/'failure-manifest.json').write_text(json.dumps(dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in original.OUT.iterdir() if p.is_file()}),indent=2)+'\n')
    source=inspect.getsource(original.symbolic)
    before='assert sp.simplify(dT-sp.Subs(sp.Derivative(F(T,b),T),b,Bt)-sp.Subs(sp.Derivative(F(T,b),b),b,Bt)*sp.diff(Bt,T))==0'
    after="dummy=sp.Dummy('temperature')\n    assert sp.simplify(dT-sp.Subs(sp.diff(F(dummy,Bt),dummy),dummy,T)-sp.Subs(sp.Derivative(F(T,b),b),b,Bt)*sp.diff(Bt,T))==0"
    assert source.count(before)==1;source=source.replace(before,after);(OUT/'corrected-symbolic.py').write_text(source)
    paths=[g.ROOT/'verification/gr_plasma_photon_runner.py',original.OUT/'failure-manifest.json',OUT/'corrected-symbolic.py']
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths})
    plan['symbolic_correction']='Represent the first argument derivative with an independent dummy temperature; retain its value substitution until after differentiation. Same chain rule, thermodynamic/ray model, numerical controls, tolerances and physical boundaries. Original symbolic implementation failure preserved.'
    save('plan.json',plan)


def namespace():
    env=dict(vars(original),OUT=OUT)
    for name in ['save','run','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    exec(compile((OUT/'corrected-symbolic.py').read_text(),str(OUT/'corrected-symbolic.py'),'exec'),env)
    return env


def run():namespace()['run']()


def verify():
    for rel,digest in json.loads((original.OUT/'failure-manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    namespace()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
