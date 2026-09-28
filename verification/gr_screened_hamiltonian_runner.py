"""Retain the first numpy-bool serialization failure; change only JSON output."""
from types import FunctionType
import json, shutil, sys
import numpy as np
import gr_screened_hamiltonian_matching as original

g=original.g;OUT=g.OUT/'gr-screened-hamiltonian-defined'


def scalar(value):
    if isinstance(value,np.generic):return value.item()
    raise TypeError(f'Unsupported JSON value: {type(value).__name__}')


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2,default=scalar)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.bindings()
    error=g.ROOT/'outputs/gr-screened-hamiltonian-matching33-run.log'
    assert 'not JSON serializable' in error.read_text() and not (original.OUT/'result.json').exists()
    shutil.copy2(error,OUT/'original-error.log')
    paths=[g.ROOT/'verification/gr_screened_hamiltonian_runner.py',original.OUT/'plan.json',
        original.OUT/'symbolic.json',original.OUT/'states.npz',OUT/'original-error.log']
    save('plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        change='Use np.generic.item() only at the JSON serialization boundary. Original equations, controls, input order and gates are unchanged. Preserve original symbolic/state outputs and require bitwise-equal state replay.',
        original_failure='First run computed states.npz but failed to serialize a numpy.bool_ in the kernel control record. No result.json or manifest.json was produced.'))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    return original.bindings()


def run():
    assert json.loads(json.dumps({'v':np.bool_(True)},default=scalar))=={'v':True}
    ns=dict(original.run.__globals__,OUT=OUT,save=save,bindings=bindings,verify=verify)
    ns['symbolic']=FunctionType(original.symbolic.__code__,ns)
    FunctionType(original.run.__code__,ns)()


def verify():
    ns=dict(original.verify.__globals__,OUT=OUT,bindings=bindings)
    FunctionType(original.verify.__code__,ns)()
    before=dict(np.load(original.OUT/'states.npz'));after=dict(np.load(OUT/'states.npz'))
    assert before.keys()==after.keys()
    for key in before:assert np.array_equal(before[key],after[key]),key
    assert (original.OUT/'symbolic.json').read_bytes()==(OUT/'symbolic.json').read_bytes()
    print('PASS serialization-only replay; original state arrays and symbolic identities unchanged',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
