"""Independent endpoint and thermodynamic replay; no EOS or time evolution."""
from pathlib import Path
import hashlib
import json
import numpy as np

out=Path(__file__).resolve().parent
read=lambda name:json.loads((out/name).read_text())
source=Path('verification/def_photon_population_coupled.py')
assert hashlib.sha256((out/'executed-source.py').read_bytes()).hexdigest()==read('execution.json')['source_sha256']
for record in [read('build.json'),read('input.json')]:
    for path,sha in record['bindings'].items():assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==sha,path
names=['fine-16','fine-32','fine-64','coarse-64','quadrature4-64']
saved={name:dict(np.load(out/(name+'.npz'))) for name in names}
t=dict(np.load(out/'thermo.npz'));r=read('response.json');fine=saved['fine-64'];Ceq=float(t['Ceq'])
assert np.allclose(t['U'].T@t['U'],float(t['temperature'])**2*1.380649e-16*np.linalg.inv(t['S']),rtol=1e-10,atol=1e-30)
assert abs(float(t['Cf'])+t['g']@t['g']-Ceq)<1e-14
def distance(x,y):return float(np.sqrt(np.linalg.norm(x['x']-y['x'])**2+np.linalg.norm(x['E']-y['E'])**2)/np.sqrt(Ceq))
errors=[distance(saved['fine-16'],saved['fine-32']),distance(saved['fine-32'],fine)]
assert errors==r['time_differences_initial']
assert float(np.log2(errors[0]/errors[1]))==r['time_order']
for name in ['coarse','quadrature4']:assert distance(fine,saved[name+'-64'])==r['input_differences_initial'][name]
bank=dict(np.load(out.parent/'def-photon-hhe-coupled/bank.npz'))
Cr=bank['Ci'][bank['u']<=60].sum()
ratio=(1+Cr/float(t['Cf']))/(1+Cr/Ceq)
prior=json.loads((out.parent/'def-photon-hhe-coupled/input.json').read_text())
bound=prior['line_integration'][-1]['omitted_wing_response_bound']*ratio
assert r['line_tail_capacity_rescale']==ratio
assert r['source_plus_two_tails']==r['input_differences_initial']['coarse']+2*bound
assert r['time_order']>=1.8 and errors[-1]<.001 and r['source_plus_two_tails']<.001
assert r['input_differences_initial']['quadrature4']<.0001
for name in names:
    row=read(name+'.json')
    assert row['balance']<1e-9 and row['solver_residual']<1e-11 and row['entropy_growth']<1e-10
    assert np.allclose(saved[name]['weights'],np.r_[np.sqrt(t['Cf']),t['g']],rtol=1e-14)
old=dict(np.load(out.parent/'def-photon-hhe-coupled/fine-64.npz'))
temperature=old['T']/np.sqrt(Ceq)
lift=dict(x=fine['weights']*temperature,E=old['E'])
assert distance(fine,lift)==r['finite_vs_LTE_difference']
assert r['actual_coupled_evolution'] and r['numerical_response_passed']
assert not any(r[k] for k in ['physical_opacity_certified','original_complete_input_failure_resolved','full_GR_photon_feedback_evolved','full_dynamic_charge_solved'])
result=dict(classification='Counterexample candidate',passed=True,endpoint_comparisons=5,
    thermo_without_double_counting=True,new_EOS_calls=0,new_evolution_paths=0,
    finite_vs_LTE_difference=r['finite_vs_LTE_difference'],
    bindings={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,*[out/(n+'.npz') for n in names],out/'thermo.npz',out/'moments.npz']})
(out/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print('PASS saved finite-population evolution, unchanged gates and native thermodynamic split')
