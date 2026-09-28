"""Replay saved coupled comparisons without any EOS call or evolution."""
from pathlib import Path
import hashlib
import json
import numpy as np

out=Path(__file__).resolve().parent
read=lambda name:json.loads((out/name).read_text())
for record in [read('plan.json'),read('input.json'),read('rates/result.json')]:
    for name,sha in record['bindings'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==sha,name
names=['fine-16','fine-32','fine-64','coarse-64','quadrature4-64','zero-atomic-64']
saved={name:dict(np.load(out/(name+'.npz'))) for name in names}
result=read('response.json');fine=saved['fine-64'];Cm=float(fine['Cm'])
def distance(x,y):
    assert np.array_equal(x['u'],y['u']) and np.array_equal(x['Ci'],y['Ci'])
    assert x['Cm']==y['Cm'] and x['duration']==y['duration']
    return float(np.sqrt(abs(x['T']-y['T'])**2+np.linalg.norm(x['E']-y['E'])**2)/np.sqrt(Cm))
errors=[distance(saved['fine-16'],saved['fine-32']),distance(saved['fine-32'],fine)]
assert errors==result['time_differences_initial']
assert float(np.log2(errors[0]/errors[1]))==result['time_order']
assert result['time_order']>=1.8 and errors[-1]<.001
for name in ['coarse','quadrature4','zero-atomic']:
    assert distance(fine,saved[name+'-64'])==result['differences_initial'][name]
tail=read('input.json')['line_integration'][-1]['omitted_wing_response_bound']
assert result['source_grid_plus_two_tails']==result['differences_initial']['coarse']+2*tail
assert result['source_grid_plus_two_tails']<.001
assert result['differences_initial']['quadrature4']<.0001
assert result['differences_initial']['zero-atomic']>100*result['source_grid_plus_two_tails']
for name in names:
    row=read(name+'.json')
    assert row['balance']<1e-9 and row['energy_residual']<1e-9
    assert row['solver_residual']<1e-11 and row['entropy_growth']<1e-10
    assert max(row['algebra']['energy_number_null_relative'])<1e-10
assert result['actual_coupled_evolution'] and result['numerical_response_passed']
assert not any(result[k] for k in ['original_complete_input_failure_resolved',
    'physical_opacity_certified','full_GR_photon_feedback_evolved','full_dynamic_charge_solved'])
print('PASS six saved coupled paths and unchanged numerical gates; original full-input failure remains open')
