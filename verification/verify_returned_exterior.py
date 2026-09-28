"""Independently check each mass/charge component and its original controls."""
from pathlib import Path
from decimal import Decimal, localcontext
import json, sys
import numpy as np
import read_returned_exterior as run


def audit(scope):
    out=run.ROOT/scope;result=run.read(out/'result.json');plan=run.read(out/'plan.json')
    assert result['passed'] and run.read(run.ROOT/f'{scope}-receipt.json')['error'] is None
    assert run.read(run.ROOT/f'{scope}-receipt.json')['source_sha256']==run.sha(Path(run.__file__))
    for p,h in plan['bindings'].items():assert run.sha(p)==h,p
    rows=result['rows'];controls={}
    for component in ['high','low']:
        values={(r['clock'],r['angular'],r['radial']):r for r in rows if r['component']==component}
        reference=values[128,8,8];controls[component]={}
        for name,setting in [('time',(64,8,8)),('angular',(128,4,8)),('radial',(128,8,4))]:
            changes={}
            for key in ['normalized_standalone','homogeneous_mass_cm','compact']:
                a=Decimal(str(values[setting][key]));b=Decimal(str(reference[key]))
                changes[key]=float(abs(a-b)/max(abs(b),Decimal('1e-290')))
            controls[component][name]=changes
    errors=[];low_errors=[]
    high=run.common.OUT if scope=='common' else run.full.OUT
    with np.load(high/'gr/source-128.npz') as d,localcontext() as ctx:
        ctx.prec=160
        def exact(v):
            a,b=np.longdouble(v).as_integer_ratio();return Decimal(a)/Decimal(b)
        alpha=-exact(d['K_cm'])/exact(d['M_cm'])
        for n in [64,128]:
            h,l=[next(v for v in rows if v['component']==part and (v['clock'],v['angular'],v['radial'])==(n,8,8)) for part in ['high','low']]
            hs,ls=Decimal(h['scalar_numerator']),Decimal(l['scalar_numerator'])
            hm=Decimal(h['kappa'])-Decimal(h['epsilon']);lm=Decimal(l['kappa'])-Decimal(l['epsilon'])
            qh=(alpha+hs)/(1+hm)-alpha
            total=(alpha+hs+ls)/(1+hm+lm)-alpha
            r=next(v for v in result['components'] if v['clock']==n)
            errors.append(float(abs(total-Decimal(r['total']))/abs(total)))
            low_errors.append(float(abs((total-qh)-Decimal(r['same_solution_low_increment']))/abs(total-qh)))
        # A positive mass port must decrease alpha; emitted energy has the
        # opposite sign. These finite controls expose a swapped convention.
        a=Decimal('.004');k=Decimal('.001')
        assert a/(1+k)-a<0 and a/(1-k)-a>0
    assert max(errors+low_errors)<1e-12,(errors,low_errors)
    passed=all(max(v.values())<plan['gates'][name] for group in controls.values() for name,v in group.items())
    files=[Path(__file__),Path(run.__file__),out/'plan.json',out/'result.json',run.ROOT/f'{scope}-receipt.json']
    verdict=dict(classification='Counterexample candidate',passed=passed,controls=controls,
        exact_rational_total_relative=errors,exact_rational_low_increment_relative=low_errors,
        pure_mass_and_arrival_sign_controls=True,
        conditional_compact_sign_survives_frozen_exterior=all(v['nominal_compact_sign_unchanged'] for v in result['components']) if passed else None,
        no_new_physical_steps=True,physical_energy_conversion_and_dynamic_exterior_complete=False,
        background_mass_remap_complete=False,physical_final_charge_solved=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):run.sha(p) for p in files})
    run.write(out/'audit.json',verdict);print(json.dumps(verdict),flush=True);assert passed,verdict


if __name__=='__main__':audit(sys.argv[1])
