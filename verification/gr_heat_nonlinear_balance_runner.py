"""Retain the failed coordinate normalization check and correct only its factor."""
from types import ModuleType
import json,sys
import gr_heat_nonlinear_balance as original

OUT=original.OUT.parent/'gr-heat-nonlinear-balance-defined'


def source():
    before=(original.g.ROOT/'verification/gr_heat_nonlinear_balance.py').read_text()
    changes={"OUT=g.OUT/'gr-heat-nonlinear-balance'":"OUT=g.OUT/'gr-heat-nonlinear-balance-defined'",
        'divergence_covariant(tensor,0)+N*BE/a':'divergence_covariant(tensor,0)+BE/a'}
    after=before
    for old,new in changes.items():assert after.count(old)==1;after=after.replace(old,new)
    return after,changes


def module(text):
    obj=ModuleType('gr_heat_nonlinear_balance_defined');exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def run():
    assert not OUT.exists();text,changes=source();obj=module(text);old_save=obj.save
    def save(name,value):
        if name=='plan.json':
            (OUT/'candidate.py').write_text(text)
            error=original.g.ROOT/'outputs/gr-heat-nonlinear-balance33-run.log';assert 'AssertionError' in error.read_text()
            paths=[original.g.ROOT/'verification/gr_heat_nonlinear_balance_runner.py',OUT/'candidate.py',original.OUT/'plan.json',error]
            value.update(checkpoint='1609b7e7',substitutions=changes,
                correction='The coordinate covector divergence is nabla_mu T^mu_t=-BE/a. The first check inserted an extra lapse N after the metric determinant had already canceled it. Correct the normalization of that equality only; the proposed energy balance BE=0, stress tensor, metric, other equations and all symbolic gates are unchanged. The original failure is preserved.')
            value['bindings'].update({p.relative_to(original.g.ROOT).as_posix():original.g.c.sha(p) for p in paths})
        old_save(name,value)
    obj.save=save;obj.run()


def verify():
    text,_=source();assert text==(OUT/'candidate.py').read_text();module(text).verify()


if __name__=='__main__':globals()[sys.argv[1]]()
