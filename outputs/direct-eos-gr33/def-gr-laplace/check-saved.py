"""Independent saved-state/readout replay; no new GR resolvents or EOS calls."""
from pathlib import Path
import json
import numpy as np
import def_gr_laplace as b

o=b.OUT
for name in ['plan.json','pilot-plan.json']:
    plan=json.loads((o/name).read_text())
    for path,sha in plan['bindings'].items():
        path=Path(path)
        if name=='pilot-plan.json' and path.name=='def_gr_laplace.py':path=o/'pilot-source.py'
        assert b.prior.digest(path)==sha,path
p=json.loads((o/'dehoog-plan.json').read_text())
assert p['source']==b.prior.digest(o/'dehoog.py')
assert p['input']==b.prior.digest(o/'fine-transform.npy')
z=np.load(o/'fine-transform.npy').reshape(513,-1,4)
assert np.array_equal(z[:16].reshape(16,-1),np.load(o/'pilot.npy'))
problem=b.Problem(b.prior.OUT/'fine-bank.npz');r=problem.bg.grid
for cutoff in [128,256,512]:
    row=json.loads((o/f'fine-{cutoff}.json').read_text())
    saved=dict(np.load(o/f'fine-{cutoff}.npz'))
    s=6+2j*np.pi*np.arange(cutoff+1)/4
    # Independent direct Fourier sum at the final time, including k=0 half.
    factor=np.exp(2j*np.pi*np.arange(cutoff+1)/4);factor[0]*=.5
    u=np.exp(6)/2*np.real(np.einsum('k,kij->ij',factor,z[:cutoff+1]))
    v=np.exp(6)/2*np.real(np.einsum('k,kij->ij',factor*s,z[:cutoff+1]))
    u[:,2]+=r*problem.bg.nodes['v']*u[:,0]
    v[:,2]+=r*problem.bg.nodes['v']*v[:,0]
    u-=problem.heat.lift(1.,problem.bg.nodes).reshape(-1,4)
    v-=problem.heat.lift(1.,problem.bg.nodes,True).reshape(-1,4)
    for actual,field in [(u,'response'),(v,'velocity')]:
        assert np.max(abs(actual-saved[field]))<2e-12*np.max(abs(saved[field])),field
    cv=np.interp(problem.native,r,problem.speed*saved['velocity'][:,0])
    cf=np.interp(problem.native,r,saved['Eulerian_scalar'])
    direct=[np.sqrt(problem.weights@(cv*cv)),np.sqrt(problem.weights@(cf*cf))]
    direct += [np.sqrt(problem.weights[m]@cv[m]**2/problem.weights[m].sum()) for m in problem.masks]
    assert np.allclose(direct,[row['history'][-1][f] for f in b.FIELDS],rtol=2e-12,atol=0)
for accelerated in [False,True]:
    result=json.loads((o/('dehoog-result.json' if accelerated else 'result.json')).read_text())
    for j,field in enumerate(b.FIELDS):
        series=[]
        for n in [128,256,512]:
            if accelerated:
                d=dict(np.load(o/f'dehoog-{n}.npz'))
                speed=d['velocity_m_s'];scalar=d['Eulerian_scalar']
                expected=np.column_stack([np.sqrt((speed*speed)@problem.weights),np.sqrt((scalar*scalar)@problem.weights)]+
                    [np.sqrt((speed[:,m]**2)@problem.weights[m]/problem.weights[m].sum()) for m in problem.masks])
                assert np.allclose(expected,d['readouts'],rtol=2e-14,atol=0)
                series.append(d['readouts'][:,j])
            else:series.append(np.array([x[field] for x in json.loads((o/f'fine-{n}.json').read_text())['history']]))
        a,c,d=series;norm=max(abs(d).max(),1e-100)
        before=max(abs(a-c))/norm;last=max(abs(c-d))/norm;order=np.log2(before/last)
        stored=result['comparisons'][field]
        assert np.allclose([before,last,order],[stored['previous' if accelerated else 'cutoff_previous'],stored['last' if accelerated else 'cutoff_last'],stored['order']],rtol=2e-13,atol=1e-14)
    assert not result['cutoff_gates_passed' if accelerated else 'passed']
control=json.loads((o/'dehoog-control.json').read_text())
assert control['errors']['300.0'][-1]>.9 and control['errors']['600.0'][-1]>.9
b.write(o/'verification.json',dict(classification='Counterexample candidate',saved_replay_passed=True,
    source_bindings_checked=True,original_initial_pilot_reused=True,failed_inversions_replayed=2,
    new_GR_resolvents=0,new_EOS_calls=0,original_failure_resolved=False))
print('PASS: exact saved resolvents, endpoint states, four readouts and both failed verdicts replayed')
