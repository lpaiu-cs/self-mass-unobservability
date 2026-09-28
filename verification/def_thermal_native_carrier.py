"""Reuse positive archived luminosities for the native evaluation carrier.

Luminosities initialize the disposable MESA carrier only; they are neither the
new GR heat flux nor a prescribed physical boundary. All source inputs are held.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import shutil
import numpy as np
import def_free_surface_thermal as old

OUT=old.OUT


def inputs():
    bg,raw,state=old.inputs();reference,_=old.h.inputs()
    state['L']=reference['L'].copy()
    assert np.all(state['L']>0)
    return bg,raw,state


def prepare():
    assert not (OUT/'carrier-plan.json').exists()
    log=old.CACHE/'pilot/execution.log'
    shutil.copy2(log,OUT/'pilot-zero-luminosity-failure.log')
    old.h.write(OUT/'carrier-plan.json',dict(classification='Counterexample candidate',
        bindings={str(p.relative_to(old.h.ROOT)):old.h.digest(p) for p in [Path(__file__),Path(old.__file__),old.h.OLD/'molecular-state-17-8.npz']},
        failure='Zero-luminosity 32-cell carrier reached the native net evaluations but finish_load_model failed in set_vars before a complete independent profile was written; no source output accepted.',
        intervention='Only the carrier luminosity is replaced by the stored positive reference luminosity. rho,T,X and the requested EOS auxiliaries stay fixed. It is not the physical new-state luminosity.',
        new_labels='Append -positive-L; preserve the entire failed native directory and original plan.'))


ns=dict(vars(old),inputs=inputs)
value_source=FunctionType(old.source.__code__,ns)


def source(label,indices,arguments):
    return value_source(label+'-positive-L',indices,arguments)


ns['source']=source
pilot=FunctionType(old.pilot.__code__,ns)
sources=FunctionType(old.sources.__code__,ns)
# The source profile name must track the carrier label; scientific equations
# and rate-value controls are exactly the frozen parent implementation.
import inspect
text=inspect.getsource(old.sources).replace("OUT/'actual-profile.data.gz'","OUT/'actual-positive-L-profile.data.gz'")
exec(compile(text,__file__,'exec'),ns)
sources=ns['sources']

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot','sources'])
    globals()[p.parse_args().action]()
