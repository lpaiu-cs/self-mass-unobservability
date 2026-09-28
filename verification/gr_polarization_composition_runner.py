"""Check the frozen composition identities with an arbitrary C2 Taylor jet."""
import inspect,json,sys
import gr_polarization_composition as base

g=base.g;ORIGINAL=base.OUT;OUT=g.OUT/'gr-polarization-composition-defined'
source=inspect.getsource(base.symbolic);OLD="fun=sp.Function('g')"
NEW="q0,T0=sp.symbols('q0 T0');jet_coeff=sp.symbols('g0 gn gT gnn gnT gTT');fun=lambda x,y:jet_coeff[0]+jet_coeff[1]*(x-q0)+jet_coeff[2]*(y-T0)+jet_coeff[3]*(x-q0)**2/2+jet_coeff[4]*(x-q0)*(y-T0)+jet_coeff[5]*(y-T0)**2/2"
assert source.count(OLD)==1
namespace=dict(base.symbolic.__globals__)
exec(compile(source.replace(OLD,NEW),__file__,'exec'),namespace)
base.symbolic=namespace['symbolic'];base.OUT=OUT


def check():
    saved=namespace['save'];namespace['save']=lambda *args:None
    try:base.symbolic()
    finally:namespace['save']=saved
    print('PASS arbitrary C2 Taylor jet composition identities',flush=True)


def prepare():
    base.prepare();p=json.loads((OUT/'plan.json').read_text())
    paths=[g.ROOT/'verification/gr_polarization_composition_runner.py',ORIGINAL/'plan.json',ORIGINAL/'symbolic-failure.log']
    p['bindings'].update({x.relative_to(g.ROOT).as_posix():g.c.sha(x) for x in paths})
    p['symbolic_repair']=dict(original=OLD,replacement=NEW,
        reason='Uninterpreted mixed derivatives remained in different Subs/Derivative forms after simplification. Check the same identities using six arbitrary coefficients of a second-order Taylor polynomial about arbitrary n_e,T. Every C2 function has such a jet at every point, so the first/second derivative identities retain their generality.',
        science_or_numerical_inputs_changed=False)
    base.save('plan.json',p)


run=base.run;verify=base.verify
if __name__=='__main__':globals()[sys.argv[1]]()
