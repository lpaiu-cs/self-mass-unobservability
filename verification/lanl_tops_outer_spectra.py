"""Preserve outer physical inputs while isolating the spectral output failure."""
from types import FunctionType
import json, shutil, sys
import lanl_tops_boundary_control as control

g=control.g;original=control.original;OUT=g.OUT/'lanl-tops-outer-spectra'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();control.verify()
    shutil.copy2(original.OUT/'lanl-tops-form.html',OUT/'lanl-tops-form.html')
    base=next(x for x in json.loads((control.audit.retrieval.OUT/'requests.json').read_text()) if x['cell']==0)
    changes=[('outer1500on',dict(temps='0.0015')),('outer2000on',dict(temps='0.002')),
             ('outerpairoff',dict(plasnu='off')),('outerpairboth',dict(datype='contgrup'))]
    rows=[dict(name=name,cell=0,fields=dict(base['fields'],mixname=name,**change)) for name,change in changes]
    save('requests.json',rows)
    paths=[g.ROOT/'verification/lanl_tops_outer_spectra.py',control.OUT/'manifest.json',
           control.audit.retrieval.OUT/'manifest.json',OUT/'requests.json',OUT/'lanl-tops-form.html']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='f305247',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        changes='The exact same outer composition and density are kept. Try each of the original two temperatures individually, the original pair with cutoff OFF, and the original pair requesting both spectral and default group output. Each query recomputes through /submit and /results; original failed pair remains frozen.',
        acceptance='Retrieval is not acceptance. Check all echoes/warnings and compare spectra with the saved gray means. A successful alternate output mode does not prove non-LTE atmospheric physics or full GR evolution.'))


def namespace():
    env=dict(vars(original),OUT=OUT)
    for name in ['save','bindings','fetch','run','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    return env


def run():namespace()['run']()


def verify():namespace()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
