"""Separate the outer-cell output failure and the central plasma-cutoff rule."""
from types import FunctionType
import json, shutil, sys
import lanl_tops_stellar_audit as audit

original=audit.retrieval.original;g=audit.g;OUT=g.OUT/'lanl-tops-boundary-control'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();audit.verify()
    shutil.copy2(original.OUT/'lanl-tops-form.html',OUT/'lanl-tops-form.html')
    prior=json.loads((audit.retrieval.OUT/'requests.json').read_text());rows=[]
    for cutoff in ['on','off']:
        base=next(x for x in prior if x['cell']==0)
        name='outergray'+cutoff;data=dict(base['fields'],mixname=name,datype='gray',plasnu=cutoff)
        rows.append(dict(name=name,cell=0,plasma_cutoff=cutoff,fields=data))
        base=next(x for x in prior if x['cell']==5734 and x['plasma_cutoff']==cutoff)
        name='coregroups'+cutoff
        data=dict(base['fields'],mixname=name,datype='groups',temps='1.75',egrid='specific',
                  energies='.01 .1 1 3 5 7 9 12 20 40')
        rows.append(dict(name=name,cell=5734,plasma_cutoff=cutoff,fields=data))
    save('requests.json',rows)
    paths=[g.ROOT/'verification/lanl_tops_boundary_control.py',audit.OUT/'manifest.json',
           audit.retrieval.OUT/'manifest.json',OUT/'requests.json',OUT/'lanl-tops-form.html']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='361c8e0',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        outer='Keep the exact original composition, density and both temperatures. Change only requested output to gray, plus a separately marked cutoff-OFF diagnostic. A gray success does not erase the original spectral failure, nor prove LTE in the outer atmosphere.',
        centre='At the same actual central mixture/density, request nine energy groups at the lower existing table temperature 1.75 keV, for both cutoff options. Check documented below-cutoff sentinel 1e10 against the actual group output; preserve unchanged monochromatic arrays and changed gray means from the prior audit.',
        acceptance='Capture every response. Follow-up must validate composition, density, temperature, warnings and group boundaries before physical use. No inference of full GR evolution, full opacity/EOS certification or an actual scattering redistribution kernel.'))


def namespace():
    env=dict(vars(original),OUT=OUT)
    for name in ['save','bindings','fetch','run','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    return env


def run():namespace()['run']()


def verify():namespace()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
