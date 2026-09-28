"""Preserve the failed H-only NCRSE assumption; validate each element's index."""
from types import FunctionType
import inspect, json, shutil, sys
import op_planck_connection as original

g=original.g;OUT=g.OUT/'op-planck-reader'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def source():
    before=inspect.getsource(original.records)
    old='assert ntot==10000 and ncoarse==1 and jstep==2'
    new="assert ntot==10000 and ncoarse==int((SOURCE/'mono'/f'm{z:02}.index').read_text().splitlines()[3].split()[0]) and jstep==2"
    assert before.count(old)==1;after=before.replace(old,new);assert after.replace(new,old)==before
    return after


def prepare():
    assert not OUT.exists();OUT.mkdir();original.bindings()
    shutil.copy2(g.ROOT/'outputs/op-planck-connection33-run.log',OUT/'original-layout-failure.log')
    for name in ['source-files.json','native-build.json']+[f'native-m{z:02}.smry' for z in original.ELEMENTS]:
        shutil.copy2(original.OUT/name,OUT/name)
    (OUT/'candidate-records.py').write_text(source())
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan['record_layout_repair']='The frozen reader incorrectly required every element to have the hydrogen coarse-grid count NCRSE=1. The actual 91 H files have 1, and the other 1456 files have 100. Validate NCRSE against each element\'s bound index file; NTOT=10000 and density-index spacing 2 stay unchanged. Reuse the completed unchanged Fortran reader outputs. No opacity or target-domain criterion is relaxed.'
    paths=[g.ROOT/'verification/op_planck_reader_runner.py',OUT/'candidate-records.py',OUT/'original-layout-failure.log',
           original.OUT/'plan.json']+[p for p in OUT.iterdir() if p.is_file()]
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths})
    save('plan.json',plan)


def namespace():
    env=dict(vars(original),OUT=OUT)
    for name in ['save','bindings','headers','connect','run','verify']:
        fn=getattr(original,name);env[name]=FunctionType(fn.__code__,env,argdefs=fn.__defaults__)
    exec(compile(source(),str(OUT/'candidate-records.py'),'exec'),env)
    assert (OUT/'candidate-records.py').read_text()==source()
    def native_control():
        build=json.loads((OUT/'native-build.json').read_text())
        assert g.c.sha(build['executable'])==build['sha256']
        for z in original.ELEMENTS:assert (OUT/f'native-m{z:02}.smry').read_bytes()==(original.OUT/f'native-m{z:02}.smry').read_bytes()
    env['native_control']=native_control
    return env


def run():namespace()['run']()


def verify():namespace()['verify']()


if __name__=='__main__':globals()[sys.argv[1]]()
