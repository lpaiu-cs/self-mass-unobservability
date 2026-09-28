"""Load immutable NPZ members and constants once for the existing full audits."""
from types import FunctionType,SimpleNamespace
from functools import lru_cache
import json,sys
import gr_dense_plasma_defined as d

g=d.g;OUT=d.OUT/'eager'


def prepare():
    assert not OUT.exists();OUT.mkdir();d.bindings();changes={}
    for name,var in [('return-run.py','a'),('screen-run.py','saved')]:
        before=(d.OUT/name).read_text();old=var+"=np.load(d.OUT/'stellar-comparison.npz')";new=var+"=dict(np.load(d.OUT/'stellar-comparison.npz'))"
        assert before.count(old)==1;after=before.replace(old,new);assert after.replace(new,old)==before
        (OUT/name).write_text(after);changes[name]=dict(old=old,new=new)
    paths=[g.ROOT/'verification/gr_dense_plasma_eager.py',d.OUT/'plan.json']+list(OUT.iterdir())
    (OUT/'plan.json').write_text(json.dumps(dict(classification='Counterexample candidate',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},substitutions=changes,
        interrupted='Read-only return audit PID8298 terminated while repeatedly reopening stellar-comparison.npz; no return-audit result existed. The shell && chain stopped before starting full screening replay. Long-running caloric/GR jobs were untouched.',
        reason='NpzFile indexing was decompressing entire members for every cell. Materialize identical stored members once. Cache only the source-derived constant dictionary in mixture(); no arithmetic, inputs, source fit, state order or acceptance threshold changes.'),indent=2)+'\n')


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel


def namespace(module,branch):
    ns=d.namespace(module,branch)
    original=d.d.mixture
    constants=lru_cache(maxsize=1)(d.d.plasma.constants)
    fn=FunctionType(original.__code__,dict(original.__globals__,plasma=SimpleNamespace(constants=constants)))
    ns['d']=SimpleNamespace(**{**vars(d.d),'mixture':fn})
    return ns


def run():
    bindings()
    for module,branch,name in [(d.audit,'return','return-run.py'),(d.screening,'screening','screen-run.py')]:
        ns=namespace(module,branch);exec(compile((OUT/name).read_text(),str(OUT/name),'exec'),ns);ns['run']()
    d.save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in d.OUT.rglob('*') if p.is_file()}))
    verify()


def verify():bindings();d.verify();print('PASS same inputs/criteria with eager immutable-array replay',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
