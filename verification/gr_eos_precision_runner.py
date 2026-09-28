"""Supply the frozen build helper's omitted default, retaining the failed wrapper."""
import inspect, json, sys
import gr_eos_extended_precision as original

g=original.g;OUT=original.OUT


def build():
    assert not (OUT/'runner.json').exists() and not (OUT/'build.json').exists()
    before=inspect.getsource(original.build)
    old='FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,save=save))'
    new='FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,save=save),argdefs=fn.__defaults__)'
    assert before.count(old)==1;after=before.replace(old,new);assert after.replace(new,old)==before
    (OUT/'candidate-build.py').write_text(after)
    original.save('runner.json',dict(classification='Counterexample candidate',
        preserved_prebuild_error="TypeError: build_at() missing 1 required positional argument: 'prefix'",
        correction='Carry over the existing build_at defaults when cloning it. The first call failed before any compiler/build command. No precision, physical option, experimental input or acceptance gate changes.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/gr_eos_precision_runner.py',OUT/'plan.json',OUT/'candidate-build.py']}))
    namespace=dict(original.build.__globals__)
    exec(compile(after,str(OUT/'candidate-build.py'),'exec'),namespace);namespace['build']()


def run():original.run();verify()


def verify():
    original.verify()
    for rel,digest in json.loads((OUT/'runner.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS precision runner provenance; frozen prebuild error retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
