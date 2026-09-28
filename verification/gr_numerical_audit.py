"""Read-only independent saved-number audits of GR mass additions and roots."""
import json,sys
from fractions import Fraction as F
import numpy as np
import gr_heat_primitive_newton as newton
import gr_mass_increment_replay as replay

g=newton.g;OUT=g.OUT/'gr-numerical-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir();newton.verify()
    paths=[g.ROOT/'verification/gr_numerical_audit.py',newton.OUT/'manifest.json',
        newton.original.OUT/'manifest.json',newton.scaling.OUT/'manifest.json',replay.OUT/'plan.json',
        *[replay.OUT/f'cell-{i}-baseline.npz' for i in range(3)]]
    save('plan.json',dict(classification='Proven',checkpoint='a2d50b8',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        scope='Exact saved mass additions and independent finite root/trace bookkeeping. No fresh physical EOS, continuum or native root error certificate.'))
    reference=np.load(g.OUT/'reference-state.npz');B=float(sum(reference['dm']))*.001*g.c.gr.G/g.c.gr.C**2
    mass_rows=[]
    for i in range(3):
        data=dict(np.load(replay.OUT/f'cell-{i}-baseline.npz'));trace=data['mass_additions']
        assert trace[0,0]==data['start'][1] and trace[-1,2]==data['end'][1]
        assert np.array_equal(trace[:-1,2],trace[1:,0])
        increments=sum((F(float(v)) for v in trace[:,1]),F(0))
        stored=F(float(data['end'][1]))-F(float(data['start'][1]))
        errors=sum((F(float(after))-F(float(before))-F(float(dy)) for before,dy,after in trace),F(0))
        assert stored-increments==errors and increments<0
        missing=errors/(-increments)
        mass_rows.append(dict(classification='Proven',cell=i,additions=len(trace),
            unchanged_additions=int(np.sum(trace[:,0]==trace[:,2])),
            proposed_decrease_geom_cm=float(-increments*F(B)*100),
            actual_decrease_geom_cm=float(-stored*F(B)*100),
            exact_missing_fraction=str(missing),missing_fraction_display=float(missing),
            exact_addition_roundoff_geom_cm=str(errors*F(B)*100),
            scope='RK addition stage only, excluding seed conversion and output multiplication. Positive missing fraction means a smaller accumulated decrease; negative means excess decrease.'))
    original=json.loads((newton.original.OUT/'result.json').read_text());scaled=json.loads((newton.scaling.OUT/'result.json').read_text())
    fixed=json.loads((newton.OUT/'result.json').read_text());plan=json.loads((newton.OUT/'plan.json').read_text())
    discrepancies=[]
    for row in fixed['rows']:
        data=dict(np.load(newton.OUT/f"cell-{row['cell']}-case-{row['case']}.npz"))
        errors=abs(data['true_primitive']-data['recovered_primitive'])
        assert errors[0]<=plan['root_log_density_tolerance'] and errors[1]<=plan['root_logT_tolerance'] and errors[2]<=plan['root_velocity_tolerance']
        assert max(abs(data['residual']))<=plan['scaled_residual_tolerance']
        discrepancies.append(errors)
    fd=json.loads((newton.scaling.OUT/'root-diagnostics.json').read_text())['rows']
    traces=json.loads((newton.OUT/'Newton-traces.json').read_text())['rows'];rejections=[]
    for trace in traces:
        history=trace['history'];assert trace['success']
        for current,following in zip(history,history[1:]):
            assert current['merit']==max(map(abs,current['residual']))
            picked=[x for x in current['trials'] if x['backtrack']==current['accepted_backtrack']][0]
            assert picked['x']==following['x'] and picked['merit']==following['merit']<current['merit']
            trial=current['trials'][0]
            if 'merit' in trial and trial['merit']>=current['merit']:
                rejections.append(dict(case_index=trace['case_index'],iteration=current['iteration'],
                    full_step_merit_ratio=trial['merit']/current['merit'],accepted_backtrack=current['accepted_backtrack']))
        assert history[-1]['merit']<=plan['Newton_scaled_residual']
    summary=dict(classification='Counterexample candidate',
        original_passed=sum(r['passed'] for r in original['rows']),unit_scaling_passed=sum(r['passed'] for r in scaled['rows']),
        decreasing_Newton_passed=len(discrepancies),total=len(fixed['rows']),
        maximum_recovery_errors_logrho_logT_v=map(float,np.max(discrepancies,axis=0)),
        maximum_scaled_residual=max(r['scaled_residual'] for r in fixed['rows']),
        maximum_initial_Jacobian_finite_row_relative_error=max(d['row_relative_error'] for a in fd for d in a['finite_differences']),
        maximum_initial_Jacobian_condition=max(a['initial_Jacobian_condition'] for a in fd),
        maximum_accepted_iterations=max(len(t['history'])-1 for t in traces),
        full_Newton_step_rejections=rejections,
        interpretation='The unchanged equations/targets recover under decreasing coupled Newton. Unit scaling alone did not repair the failures. Observed full-step residual increases establish nonlinear step difficulty; the internal MINPACK failure mechanism is not completely identified.')
    summary['maximum_recovery_errors_logrho_logT_v']=list(summary['maximum_recovery_errors_logrho_logT_v'])
    save('mass-additions.json',dict(classification='Proven',rows=mass_rows,exact_saved_identities_passed=True))
    save('primitive-summary.json',summary)
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    print('PASS mass addition identities and 27 unchanged-target nonlinear primitive recoveries',flush=True);verify()


def verify():
    for rel,digest in json.loads((OUT/'plan.json').read_text())['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'mass-additions.json').read_text())['exact_saved_identities_passed']
    report=json.loads((OUT/'primitive-summary.json').read_text());assert report['decreasing_Newton_passed']==report['total']==27
    print('PASS independent GR numerical audit bindings',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
