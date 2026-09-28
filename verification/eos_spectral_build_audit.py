"""Preserve and repair a build-metadata spill; frozen libraries are never rebuilt."""
import json, subprocess, sys
from types import FunctionType
import eos_h2plus_spectral as previous

g=previous.g;OUT=g.OUT/'eos-spectral-build-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def build(candidate):
    """Redirect the writer as well as its caller's output path."""
    path=previous.switch.OUT/'bridge.json';before=path.read_bytes()
    fn=previous.switch.e.build
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=candidate.OUT,CACHE=candidate.CACHE,
        NAME=candidate.NAME,LIB=candidate.LIB,save=candidate.save),closure=fn.__closure__)()
    assert path.read_bytes()==before,'parent metadata was modified'


def repair():
    assert not OUT.exists();OUT.mkdir()
    path=previous.switch.OUT/'bridge.json';rel=path.relative_to(g.ROOT).as_posix()
    original=subprocess.check_output(['git','show','15cb34b:'+rel],cwd=g.ROOT)
    expected=previous.switch.read('manifest.json')['sha256'][rel]
    import hashlib
    assert hashlib.sha256(original).hexdigest()==expected
    wrong=path.read_bytes();assert wrong!=original
    record=json.loads(wrong);assert record['library_sha256']==g.c.sha(previous.LIB)
    assert record['bridge_sha256']==g.c.sha(previous.CACHE/'excitation.so')
    (OUT/'misdirected-h2plus-bridge.json').write_bytes(wrong)
    (OUT/'original-parent-bridge.json').write_bytes(original)
    save('repair.json',dict(classification='Counterexample candidate',checkpoint='15cb34b',
        cause='The previous FunctionType build overlay redirected OUT but retained the parent save callable. Only bridge.json was misdirected; both libraries and their bridge binaries retain their frozen SHA values.',
        restored_file=rel,restored_sha256=expected,
        correction='Preserve the misdirected record here, restore parent bytes from its frozen commit and manifest, and redirect save explicitly in future candidate builds. Do not rerun the frozen faulty build wrapper.',
        original_runtime=previous.switch.read('manifest.json')['runtime'],
        candidate_runtime=json.loads((previous.OUT/'manifest.json').read_text())['runtime']))
    path.write_bytes(original)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in
        [*OUT.iterdir(),g.ROOT/'verification/eos_spectral_build_audit.py']}))
    verify()


def verify():
    record=json.loads((OUT/'repair.json').read_text())
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(g.ROOT/record['restored_file'])==record['restored_sha256']
    for key in ['original_runtime','candidate_runtime']:
        for path,digest in record[key].items():assert g.c.sha(path)==digest,path
    print('PASS spectral build metadata restoration; all four frozen binaries unchanged',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
