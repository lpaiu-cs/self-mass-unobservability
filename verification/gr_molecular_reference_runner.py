"""Use the already-defined 21-output Python EOS boundary, preserving the failed guard."""
from types import FunctionType
import inspect, json, shutil, sys
import gr_molecular_reference as original

g=original.g;OUT=g.OUT/'gr-molecular-reference-defined'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    before=inspect.getsource(original.block)
    old="a=snap['eos'];assert a.shape==(22,)\n        a=np.delete(a,20);b=old(2,r,t,X);assert a.shape==b.shape==(21,)"
    new="a=snap['eos'];b=old(2,r,t,X);assert a.shape==b.shape==(21,)"
    assert before.count(old)==1;after=before.replace(old,new);assert after.replace(new,old)==before
    (OUT/'candidate-block.py').write_text(after)
    shutil.copy2(g.ROOT/'outputs/gr-molecular-reference33-run.log',OUT/'original-shape-failure.log')
    plan=json.loads((original.OUT/'plan.json').read_text())
    plan.update(checkpoint='226f237',correction='The live snapshot and EOS __call__ both return 21 defined outputs, inherited from direct_eos_gr.EOS. The frozen failed guard confused this Python boundary with its internal 22-output C buffer. Remove the extra deletion and compare the two measured 21-output Python APIs. No physical model, data or acceptance-gate change.')
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/gr_molecular_reference_runner.py',OUT/'candidate-block.py',OUT/'original-shape-failure.log',original.OUT/'plan.json']})
    save('plan.json',plan)


def bindings():return FunctionType(original.bindings.__code__,dict(original.bindings.__globals__,OUT=OUT))()


def block(start):
    namespace=dict(original.block.__globals__,OUT=OUT,bindings=bindings,save=save)
    exec(compile((OUT/'candidate-block.py').read_text(),str(OUT/'candidate-block.py'),'exec'),namespace)
    return namespace['block'](start)


def run():
    fn=original.run
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,bindings=bindings,save=save,block=block,verify=verify))()


def verify():
    fn=original.verify
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,bindings=bindings))()


if __name__=='__main__':globals()[sys.argv[1]]()
