"""Extended arithmetic for the same declared two-spectrum molecular EOS."""
from types import FunctionType
import ctypes, json, shutil, sys
import mpmath as mp
import numpy as np
import eos_molecular_spectral as model
import gr_eos_precision_runner as previous
import gr_molecular_reference_runner as reference

q=previous.original;g=q.g;OUT=g.OUT/'gr-molecular-precision';CACHE=g.CACHE/'molecular-precision'
NAME='free_eos_direct24_molecular_extended';BRIDGE=CACHE/'split-input.so'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    model.verify();previous.verify();reference.verify()
    source=CACHE/'source';parent=model.CACHE/'source';shutil.copytree(parent,source)
    path=source/'src/CMakeLists.txt';before=path.read_text();after=before
    pairs=[('set(supported_fp_precision 15 CACHE','set(supported_fp_precision 33 CACHE'),
        ('set(supported_fp_exponent 300 CACHE','set(supported_fp_exponent 4000 CACHE'),
        ('OUTPUT_NAME '+model.NAME,'OUTPUT_NAME '+NAME)]
    for a,b in pairs:assert after.count(a)==1;after=after.replace(a,b)
    reverse=after
    for a,b in reversed(pairs):assert reverse.count(b)==1;reverse=reverse.replace(b,a)
    assert reverse==before;path.write_text(after)
    shutil.copy2(parent/'src/CMakeLists.txt',OUT/'before-CMakeLists.txt');shutil.copy2(path,OUT/path.name)
    shutil.copy2(q.OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    shutil.copy2(q.OUT/'candidate-build.py',OUT/'candidate-build.py')
    save('source-tree.json',dict(before={p.relative_to(parent).as_posix():g.c.sha(p) for p in parent.rglob('*') if p.is_file()},
        after={p.relative_to(source).as_posix():g.c.sha(p) for p in source.rglob('*') if p.is_file()},substitutions=pairs))
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='cf062eb',
        cells=[0,1175,1176,2972,5734],log_temperature_offsets=[-.001,0.,.001],
        root_iterations=15,known_root_log_temperature_tolerance=1e-10,
        intervention='Use exactly the new two-spectrum molecular source and chemical anchors with upstream 33-digit arithmetic, the already validated split-input bridge, and retain_molecules=True once per initialized library. Retain physical coefficients, Fermi tolerances and binary64 LAPACK. Do not use the original-model extended EOS in a new-model root.',
        gate='At the same target entropy and pressure, require the original max(2 erg/g,32 ulp(abs(H))) energy-unit gate and known-root lnT error below 1e-10. Record arithmetic differences at identical binary64 inputs without turning them into physical error bounds.',
        boundary='Extended arithmetic and finite known-root tests do not certify continuous native, mixed-precision switch, Fermi, plasma, occupation or physical EOS error.',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/gr_molecular_precision.py',g.ROOT/'verification/gr_eos_extended_precision.py',
            model.OUT/'manifest.json',q.OUT/'manifest.json',reference.OUT/'manifest.json',
            OUT/'source-tree.json',OUT/'direct_ion_bridge.f90',OUT/'candidate-build.py']}))


def build():
    namespace=dict(q.build.__globals__,OUT=OUT,CACHE=CACHE,NAME=NAME,BRIDGE=BRIDGE,save=save)
    exec(compile((OUT/'candidate-build.py').read_text(),str(OUT/'candidate-build.py'),'exec'),namespace)
    namespace['build']()


class EOS(q.EOS):
    def __init__(self):
        fn=q.EOS.__init__;FunctionType(fn.__code__,dict(fn.__globals__,BRIDGE=BRIDGE))(self)
        ctypes.c_int.in_dll(self.lib,'__mod_free_eos_MOD_retain_molecules').value=1


def run():
    mp.mp.dps=80;plan=json.loads((OUT/'plan.json').read_text())
    data=dict(np.load(reference.OUT/'reference-state.npz'));eos=EOS();binary=model.EOS();rows=[]
    for i in plan['cells']:
        for delta in plan['log_temperature_offsets']:
            lp=mp.mpf(float(data['logP'][i]));t=mp.mpf(float(data['lnT'][i]+delta));X=data['X'][i]
            a=eos(lp,t,X);b=binary(1,float(lp),float(t),X)
            best,trace=q.solve(eos,lp,a[3],X,t+mp.mpf('.01'),plan['root_iterations'])
            error=abs(best[1]-t);passed=best[0]<=1 and error<plan['known_root_log_temperature_tolerance']
            rows.append(dict(cell=i,delta=delta,score=float(best[0]),lnT_error=mp.nstr(error,50),trace=trace,
                new_minus_binary_at_same_input=[mp.nstr(v-mp.mpf(float(w)),40) for v,w in zip(a,b)],passed=bool(passed)))
            print('MOLECULAR EXTENDED ROOT',i,delta,float(best[0]),float(error),flush=True)
    save('result.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        all_passed=all(r['passed'] for r in rows),physical_EOS_certified=False,continuous_errors_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in json.loads((OUT/'build.json').read_text())['sha256'].items():assert g.c.sha(path)==digest,path
    tree=json.loads((OUT/'source-tree.json').read_text())
    for key,folder in [('before',model.CACHE/'source'),('after',CACHE/'source')]:
        for rel,digest in tree[key].items():assert g.c.sha(folder/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS molecular extended evaluator provenance and finite controls; no physical certificate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
