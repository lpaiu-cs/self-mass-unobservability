"""Separate the documented exact cutoff sentinel from large physical opacity."""
from types import FunctionType
import inspect, json, shutil, sys
import lanl_tops_outer_groups as original

g=original.g;OUT=g.OUT/'lanl-tops-outer-group-reader'


def prepare():
    assert not OUT.exists();OUT.mkdir();original.verify()
    plan=json.loads((original.OUT/'plan.json').read_text())
    source=inspect.getsource(original.analyze);before='np.all(group[:,1:]<1e10)';after='np.all(group[:,1:]!=1e10)'
    assert source.count(before)==1;source=source.replace(before,after);(OUT/'corrected-analyze.py').write_text(source)
    for path in original.OUT.iterdir():
        if path.name=='requests.json' or path.name.endswith(('-table.txt','-results-request.json')):shutil.copy2(path,OUT/path.name)
    paths=[g.ROOT/'verification/lanl_tops_outer_group_reader.py',original.OUT/'analysis-manifest.json']+list(OUT.iterdir())
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths})
    plan['reader_correction']='The earlier implementation used opacity<1e10, an undocumented ceiling, and rejected finite values 2.8e10 to 2.5e11 even with cutoff OFF. Provider documentation specifies the exact sentinel 1e10, not a physical upper bound. A separately classified reader rejects equality to the sentinel while retaining all input, positivity, finite 1e-3 integration and refinement gates. Preserve the original failed analysis and distinguish this source-informed gate correction from passing its stricter old ceiling.'
    (OUT/'plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')


def namespace():
    env=dict(vars(original),OUT=OUT)
    for name in ['save','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    exec(compile((OUT/'corrected-analyze.py').read_text(),str(OUT/'corrected-analyze.py'),'exec'),env)
    return env


def run():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    namespace()['analyze']()


def verify():original.verify();namespace()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
