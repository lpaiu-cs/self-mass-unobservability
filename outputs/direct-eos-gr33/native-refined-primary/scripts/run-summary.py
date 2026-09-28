"""Phase267: summary of the refined 64-clock run (segments, gates, fallbacks, approved exceptions)."""
import glob, json, os
W = '/home/lpaiu/work/native-refined267-runtime/primary267-refined64-work/'
rows = []; seconds = 0.; linear_exc = []; nonlinear_exc = []; fallbacks = 0
for f in sorted(glob.glob(W + 'seg-*-driver.json'), key=lambda p: int(p.split('seg-')[1].split('-')[0])):
    d = json.load(open(f)); r = d.get('row', {}); seconds += d['seconds']; fallbacks += d.get('fallbacks') or 0
    linear_exc += d.get('vector_gate_exceptions') or []; nonlinear_exc += d.get('nonlinear_exceptions') or []
    rows.append(dict(segment=os.path.basename(f)[:6], passed=d['passed'], seconds=round(d['seconds']), actual_steps=r.get('actual_completed_steps'),
                     linear=r.get('linear_relative'), energy=r.get('energy_balance_relative'), species=r.get('species_balance_relative'),
                     krylov=r.get('max_Krylov_iterations'), jet=r.get('velocity_jet_relative'), fallbacks=d.get('fallbacks'),
                     linear_exceptions=len(d.get('vector_gate_exceptions') or []), nonlinear_exceptions=len(d.get('nonlinear_exceptions') or [])))
summary = dict(segments=rows, total_seconds=seconds, total_fallbacks=fallbacks,
               linear_exceptions=dict(count=len(linear_exc), max_vector=max((e['vector_relative'] for e in linear_exc), default=None),
                                      max_physical=max((e['physical_moment_max'] for e in linear_exc), default=None), max_material=max((e['material_component_max'] for e in linear_exc), default=None)),
               nonlinear_exceptions=nonlinear_exc)
open(W + 'run-summary.json', 'w').write(json.dumps(summary, indent=1) + '\n')
for r in rows: print(r)
print('total seconds', round(seconds), 'fallbacks', fallbacks, 'linear exceptions', summary['linear_exceptions'])
for e in nonlinear_exc: print('nonlinear', e)
