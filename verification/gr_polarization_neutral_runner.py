"""Use simplified-difference equality in the frozen predictor-bound check."""
import inspect,json,sys
import gr_polarization_neutral_error as base

g=base.g;ORIGINAL=base.OUT;OUT=g.OUT/'gr-polarization-neutral-defined'
source=inspect.getsource(base.symbolic)
old='assert sp.factor(1/(4*s)-x/(x+s)**2)==(s-x)**2/(4*s*(s+x)**2)'
new='assert sp.simplify(1/(4*s)-x/(x+s)**2-(s-x)**2/(4*s*(s+x)**2))==0'
assert source.count(old)==1
namespace=dict(base.symbolic.__globals__)
exec(compile(source.replace(old,new),__file__,'exec'),namespace)
base.symbolic=namespace['symbolic'];base.OUT=OUT


def prepare():
    base.prepare();p=json.loads((OUT/'plan.json').read_text())
    paths=[g.ROOT/'verification/gr_polarization_neutral_runner.py',ORIGINAL/'plan.json',ORIGINAL/'failure.json']
    p['bindings'].update({x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in paths})
    p['symbolic_repair']=dict(original=old,replacement=new,science_or_budget_changed=False)
    base.save('plan.json',p)


run=base.run;verify=base.verify
if __name__=='__main__':globals()[sys.argv[1]]()
