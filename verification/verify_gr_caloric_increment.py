"""Exact saved-number audit; separate from native/continuous EOS certification."""
from fractions import Fraction as F
from decimal import Decimal, localcontext, ROUND_CEILING
import json,sys
import numpy as np
import gr_caloric_increment as run

g=run.g;OUT=g.OUT/'gr-caloric-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def rational(value):return F(*value.as_integer_ratio())
def upper(value):
    with localcontext() as ctx:
        ctx.prec=40;ctx.rounding=ROUND_CEILING
        return str(Decimal(value.numerator)/Decimal(value.denominator))


def prepare():
    assert not OUT.exists();OUT.mkdir()
    inputs=[g.ROOT/'verification/verify_gr_caloric_increment.py',run.OUT/'plan.json',g.OUT/'initial-state-17-4.npz']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='a4fb34c',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in inputs},
        already_available_paths_at_prepare=[n for n in [1,2,4] if (run.OUT/f'path-{n}.json').exists()],
        scope='Convert stored binary64 and long-double values by exact as_integer_ratio, with no cast to binary64 first. Audit the caloric coordinate and fixed reconstructed binary64 lapse weights; no native primitive, continuum exp, quadrature-remainder or physical error certificate.',
        corruption_control='Defined after the first path was seen: adding one percent of the exact exchanged energy to the total increment must fail the original global 1e-8 gate. Does not mutate the original output or count as a preregistered native fault injection.'))


def audit(steps):
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    record=json.loads((run.OUT/f'path-{steps}.json').read_text());assert record['completed']
    source_plan=json.loads((run.OUT/'plan.json').read_text());assert steps in source_plan['step_counts']
    for rel,digest in {**source_plan['bindings'],**source_plan['diagnostic_binding']}.items():assert g.c.sha(g.ROOT/rel)==digest,rel
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));dm=state['dm'];weight=dm*np.exp(state['nu'])
    masses=list(map(rational,weight));baryons=list(map(rational,dm))
    duration=json.loads((g.OUT/'gr-nonlinear-thermal/duration.json').read_text())['coordinate_seconds'];dt=F(duration/steps)
    sums=[];denominators=[];rows=[];energy_total=np.zeros(len(dm),dtype=np.longdouble);shift_total=np.zeros_like(energy_total)
    paths=[run.OUT/f'path-{steps}.json',run.OUT/f'path-{steps}.npz',g.OUT/'gr-nonlinear-thermal/duration.json']
    for i in range(steps):
        path=run.OUT/f'endpoint-{steps}-{i}.npz';data=dict(np.load(path));paths.append(path)
        U=list(map(rational,data['caloric_increment']));low=list(map(rational,data['two_point_increment']))
        L=[F(0),*map(rational,data['interior_Linf']),F(0)];target=[dt*(b-a) for a,b in zip(L,L[1:])]
        assert sum(target)==0
        energy=[m*u for m,u in zip(masses,U)];residual=[e-t for e,t in zip(energy,target)]
        den=sum(map(abs,target));assert den>0
        total=sum(energy);ratio=abs(total)/den;local=max(abs(r)/(m*rational(c)) for r,m,c in zip(residual,masses,data['eos'][:,10]))
        quadrature=sum(m*abs(u-v) for m,u,v in zip(masses,U,low))/den
        entropy=sum(b*rational(s) for b,s in zip(baryons,data['entropy_increment']))
        assert ratio<F(str(source_plan['global_energy_relative_to_exchange_tolerance']))
        assert local<F(str(source_plan['local_energy_residual_scaled_tolerance']))
        assert quadrature<F(str(source_plan['finite_quadrature_difference_tolerance'])) and entropy>=0
        corrupted=abs(total+den/100)/den;assert corrupted>F('1e-8')
        energy_total+=weight.astype(np.longdouble)*data['caloric_increment'];shift_total+=data['step_shift'].astype(np.longdouble)
        sums.append(total);denominators.append(den)
        rows.append(dict(step=i,exact_frozen_global_energy_ratio_upper=upper(ratio),
            exact_frozen_local_scaled_residual_upper=upper(local),
            exact_frozen_two_four_point_difference_ratio_upper=upper(quadrature),
            exact_frozen_entropy_nonnegative=entropy>=0,post_result_corruption_control_rejected=True))
    final=dict(np.load(run.OUT/f'path-{steps}.npz'))
    assert np.array_equal(energy_total,final['total_energy_increment'])
    assert np.array_equal(shift_total,final['lnT_shift'])
    value=dict(classification='Proven',steps=steps,passed=True,rows=rows,
        exact_frozen_path_global_energy_ratio_upper=upper(abs(sum(sums))/sum(denominators)),
        final_stored_energy_and_temperature_increment_replayed=True,
        no_native_continuous_or_physical_certificate=True,
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths})
    if steps>1 and (run.OUT/f'path-{steps//2}.npz').exists():
        previous=run.OUT/f'path-{steps//2}.npz';old=dict(np.load(previous))
        error=max(abs(rational(a)-rational(b)) for a,b in zip(final['lnT_shift'],old['lnT_shift']))
        value['finite_time_refinement']=dict(classification='Counterexample candidate',
            exact_saved_logT_difference_upper=upper(error),passed=error<F(str(source_plan['finite_time_refinement_logT_tolerance'])),
            rigorous_time_error_enclosure=False)
        value['bindings'][previous.relative_to(g.ROOT).as_posix()]=g.c.sha(previous)
    save(f'path-{steps}.json',value)
    print('EXACT CALORIC AUDIT',steps,value['exact_frozen_path_global_energy_ratio_upper'],value.get('finite_time_refinement'),flush=True)


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path in sorted(OUT.glob('path-*.json')):
        record=json.loads(path.read_text());assert record['passed']
        for rel,digest in record['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    print('PASS saved caloric path bindings',flush=True)


if __name__=='__main__':
    if sys.argv[1]=='audit':audit(int(sys.argv[2]))
    else:globals()[sys.argv[1]]()
