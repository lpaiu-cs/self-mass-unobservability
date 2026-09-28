"""Accumulate positive residual bounds with exact outward dyadic rounding."""
from fractions import Fraction as Q
from types import FunctionType
import inspect, json, shutil, sys
import gr_scalar_residual_certificate as original

g=original.g;OUT=g.OUT/'gr-scalar-residual-dyadic'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def ceil_bound(value):
    assert value>=0
    scale=1<<128
    result=Q((value.numerator*scale+value.denominator-1)//value.denominator,scale)
    assert value<=result<value+Q(1,scale)
    return result


def prepare():
    assert not OUT.exists()
    fn=original.prepare;FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,save=save))()
    before=inspect.getsource(original.run);old='error+=err';new='error+=ceil_bound(err)'
    assert before.count(old)==1;after=before.replace(old,new);assert after.replace(new,old)==before
    (OUT/'candidate-run.py').write_text(after)
    shutil.copy2(g.ROOT/'outputs/gr-scalar-residual-certificate33-run.log',OUT/'original-serialization-failure.log')
    plan=json.loads((OUT/'plan.json').read_text())
    plan.update(accumulation='The initial exact sum of rational positive bounds acquired a denominator too large for the 4300-decimal-digit output guard. Preserve that run and its piecewise bounds. In this separate certificate, round each nonnegative piece bound UP to the next multiple of 2^-128 using integer ceiling, then sum exactly. The added bound is strictly below 5736*2^-128. No guard, physical coefficient, residual polynomial or acceptance criterion is weakened.',
        rounded_bound_bits=128,checkpoint='427ebf8')
    paths=[g.ROOT/'verification/gr_scalar_residual_runner.py',OUT/'candidate-run.py',
        OUT/'original-serialization-failure.log',original.OUT/'plan.json',original.OUT/'controls.json',original.OUT/'piece-residuals.json']
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths});save('plan.json',plan)
    assert ceil_bound(Q(0))==0 and ceil_bound(Q(1,3))>=Q(1,3) and ceil_bound(Q(1,2))==Q(1,2)


def run():
    source=(OUT/'candidate-run.py').read_text()
    assert source.replace('error+=ceil_bound(err)','error+=err')==inspect.getsource(original.run)
    namespace=dict(original.run.__globals__,OUT=OUT,save=save,ceil_bound=ceil_bound,verify=verify)
    exec(compile(source,str(OUT/'candidate-run.py'),'exec'),namespace);namespace['run']()


def verify():
    fn=original.verify;FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT))()


if __name__=='__main__':globals()[sys.argv[1]]()
