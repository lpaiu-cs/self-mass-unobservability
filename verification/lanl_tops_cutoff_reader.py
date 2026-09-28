"""Read TOPS's returned lower/upper energy fields without changing audit gates."""
from types import FunctionType
import inspect, json, shutil, sys
import lanl_tops_cutoff_audit as original

g=original.g;OUT=g.OUT/'lanl-tops-cutoff-reader'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    plan=json.loads((original.OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert not (original.OUT/'result.json').exists()
    shutil.copy2(g.ROOT/'outputs/lanl-tops-cutoff33-run.log',original.OUT/'failure.log')
    (original.OUT/'failure.json').write_text(json.dumps(dict(classification='Counterexample candidate',
        failed=True,reason="Reader expected returned 'energies', but actual form uses egplow/egphigh; no scientific gate evaluated."),indent=2)+'\n')
    (original.OUT/'failure-manifest.json').write_text(json.dumps(dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in original.OUT.iterdir() if p.is_file()}),indent=2)+'\n')
    source=inspect.getsource(original.header)
    before="assert second['energies'].split()==submitted"
    after="assert F(Decimal(second['egplow']))==F(Decimal(submitted[0])) and F(Decimal(second['egphigh']))==F(Decimal(submitted[-1]))"
    assert source.count(before)==1;source=source.replace(before,after)
    (OUT/'corrected-header.py').write_text(source)
    for name in ['lanl-opacity-tables2025.pdf','lanl-opacity-tables2025.txt']:
        shutil.copy2(original.OUT/name,OUT/name)
    paths=[g.ROOT/'verification/lanl_tops_cutoff_reader.py',original.OUT/'failure-manifest.json',OUT/'corrected-header.py']
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths})
    plan['reader_correction']='The full submitted grid is checked by the returned lower-edge rows plus the actual returned egplow/egphigh. Compare exact decimal numbers; only a nonexistent field lookup changes. Original failure is frozen. Numerical and physical acceptance gates are unchanged.'
    save('plan.json',plan)


def namespace():
    env=dict(vars(original),OUT=OUT)
    exec(compile((OUT/'corrected-header.py').read_text(),str(OUT/'corrected-header.py'),'exec'),env)
    for name in ['save','run','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    return env


def run():namespace()['run']()


def verify():
    for rel,digest in json.loads((original.OUT/'failure-manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    namespace()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
