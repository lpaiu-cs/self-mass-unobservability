"""Correct the EOS energy normalization; preserve the rejected predecessor."""
from pathlib import Path
from types import FunctionType
import argparse
import inspect
import json
import numpy as np
import def_free_surface_response as old

h=old.h
OUT=old.OUT/'normalized'


def eos_cx(X):
    # EOS neutral-element normalization and isotope rest-mass CX are distinct.
    eos=h.molecular.model.EOS()
    return float(((X/h.molecular.g.c.A)@eos.mapping)@eos.weights)


text=inspect.getsource(old.build)
before="u0=-cx0*a[-1,11]"
after="u0=-eos_cx(X)*a[-1,11]"
assert text.count(before)==1
text=text.replace(before,after).replace('files=[Path(__file__),','files=[Path(__file__),Path(old.__file__),h.molecular.g.d.OUT/\'model-data.json\',h.ROOT/\'verification/direct_ion_eos.py\',')
namespace=dict(vars(old),old=old,OUT=OUT,__file__=__file__,eos_cx=eos_cx)
exec(compile(text,__file__,'exec'),namespace)
for name in ['sample','operators','solve','pilot','run']:
    namespace[name]=FunctionType(getattr(old,name).__code__,namespace)


def build():
    namespace['build']()
    data,_=h.inputs();cx=eos_cx(data['X'][0]);rest=float(data['CX'][0])
    h.write(OUT/'normalization-correction.json',dict(classification='Counterexample candidate',
        rest_CX=rest,EOS_CX=cx,old_source_sha256=h.digest(Path(old.__file__)),corrected_source_sha256=h.digest(Path(__file__)),
        old_build_rejected=True,old_surface_radius_input_withdrawn=True,
        cause='Phase48 multiplied the native H2 binding output (per EOS neutral-element gram) by isotope-rest CX instead of the EOS normalization computed by direct_ion_eos.EOS.__call__. The warm enthalpy differences remain invariant; the inferred zero-temperature endpoint shifts.',
        substitution=dict(old=before,new=after),surface_radius_m=float(np.load(OUT/'background.npz')['r'][-1]*np.load(OUT/'background.npz')['R']),
        physical_EOS_certified=False,full_dynamic_charge_solved=False))


def correct_surface():
    original=old.surface
    source=inspect.getsource(original.main)
    a='h0=cx-cx*mp.mpf(float(a[11]))/(c*c)'
    b="h0=cx-mp.mpf(eos_cx(data['X'][0]))*mp.mpf(float(a[11]))/(c*c)"
    assert source.count(a)==1;source=source.replace(a,b)
    source=source.replace('paths=[Path(__file__),','paths=[Path(__file__),Path(original.__file__),h.molecular.g.d.OUT/\'model-data.json\',h.ROOT/\'verification/direct_ion_eos.py\',')
    source=source.replace('multiplied by CX for baryon-gram h0','multiplied by the EOS neutral-element normalization, distinct from isotope-rest CX, for baryon-gram h0')
    ns=dict(vars(original),original=original,OUT=OUT/'surface-enclosure-corrected',__file__=__file__,eos_cx=eos_cx)
    exec(compile(source,__file__,'exec'),ns);ns['main']()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['build','correct_surface','pilot','run']);action=p.parse_args().action
    globals()[action]() if action in ['build','correct_surface'] else namespace[action]()
