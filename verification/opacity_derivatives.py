"""Finite native opacity checks along explicitly linear electron-input paths.

Counterexample candidate. These are derivatives of the saved opacity model at
fixed nuclear composition, not certified physical opacity or full transport.
"""
import json, sys
import numpy as np
import native_opacity as o

g=o.g;OUT=o.OUT/'derivatives'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert json.loads((o.OUT/'identity-control.json').read_text())['passed']
    assert not OUT.exists();OUT.mkdir()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='9d30f86',
        inputs_sha256={str(p.relative_to(g.ROOT)):g.c.sha(p) for p in [
            o.OUT/'baseline-type1-captured.npz',o.OUT/'identity-control.json',
            g.OUT/'reference-state.npz',g.ROOT/'verification/native_opacity.py',
            g.ROOT/'verification/opacity_derivatives.py']},
        log_steps=[1e-3,5e-4],relative_tolerance=1e-3,absolute_floor=1,
        path='At each original cell hold composition fixed, vary ln rho or ln T, and set lnfree_e to its baseline value plus the corresponding supplied baseline slope times the actual log displacement. The derivative of this path at zero is exactly the supplied electron derivative; no assumption about finite-step EOS curvature is needed.',
        comparison='Centered differences of ln opacity versus the returned baseline logarithmic derivatives. Compare both steps and their mutual difference, recording every failed cell without relaxing the gate.',
        scope='Finite perturbations of the unchanged native Type1 opacity model, including its input rounding and table interpolation. These are not continuous derivative enclosures or new-EOS finite curves.',
        physical_opacity_certified=False,full_GR_transport=False))


def case(axis,h,sign):
    label=f'derivative-{axis}-{h:g}-'+('plus' if sign>0 else 'minus')
    trace=o.OUT/(label+'-trace.json')
    if trace.exists():
        assert json.loads(trace.read_text())['profile_passed'];return label
    data=dict(np.load(g.OUT/'reference-state.npz'))
    baseline=dict(np.load(o.OUT/'baseline-type1-captured.npz'))
    data[axis]=data[axis]+sign*h
    g.c.native.OUT=o.OUT;g.c.native.CACHE=o.CACHE
    g.c.native.setup(label,data,species=g.c.NAMES,
        network=(g.c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
    j=0 if axis=='lnd' else 1
    displacement=data[axis]-baseline['parameters'][:,3+j]*np.log(10)
    replacement=baseline['used'].copy();replacement[:,0]+=replacement[:,1+j]*displacement
    o.trace(label,replacement)
    actual=dict(np.load(o.OUT/(label+'-captured.npz')))
    assert np.array_equal(actual['used'],replacement)
    assert np.array_equal(actual['X'],baseline['X'])
    assert np.max(abs(actual['parameters'][:,3+j]*np.log(10)-data[axis]))<1e-12
    return label


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['inputs_sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    baseline=dict(np.load(o.OUT/'baseline-type1-captured.npz'));rows=[];slopes=[];bindings={}
    for j,axis in enumerate(['lnd','lnT']):
        direction=[]
        for h in plan['log_steps']:
            labels=[case(axis,h,sign) for sign in [-1,1]]
            pair=[dict(np.load(o.OUT/(label+'-captured.npz'))) for label in labels]
            for label in labels:
                for suffix in ['-captured.npz','-trace.json','-input.npz','-inputs.json']:
                    p=o.OUT/(label+suffix);bindings[str(p.relative_to(g.ROOT))]=g.c.sha(p)
            step=(pair[1]['parameters'][:,3+j]-pair[0]['parameters'][:,3+j])*np.log(10)
            assert np.all(abs(step/(2*h)-1)<1e-10)
            numeric=np.log(pair[1]['outputs'][:,0]/pair[0]['outputs'][:,0])/step
            analytic=baseline['outputs'][:,1+j]
            score=abs(numeric-analytic)/np.maximum(plan['absolute_floor'],abs(analytic))
            direction.append(numeric);worst=int(np.argmax(score))
            rows.append(dict(axis=axis,step=h,maximum_score=float(score[worst]),worst_cell=worst,
                returned_derivative=float(analytic[worst]),finite_derivative=float(numeric[worst]),
                failed_cells=np.flatnonzero(score>plan['relative_tolerance']).tolist(),
                passed=bool(np.all(score<=plan['relative_tolerance']))))
        slopes.append(direction)
    slopes=np.array(slopes);mutual=abs(slopes[:,0]-slopes[:,1])/np.maximum(1,abs(slopes[:,1]))
    np.savez_compressed(OUT/'finite-slopes.npz',slopes=slopes,mutual_scores=mutual)
    passed=all(row['passed'] for row in rows) and np.all(mutual<=plan['relative_tolerance'])
    save('result.json',dict(classification='Counterexample candidate',passed=bool(passed),rows=rows,
        maximum_two_step_difference=float(mutual.max()),sha256=bindings,
        physical_or_continuous_derivative_certificate=False))
    print('NATIVE OPACITY DERIVATIVES',bool(passed),rows,'two-step',float(mutual.max()),flush=True)
    assert passed,'Preserved failed native opacity derivative gate'


if __name__=='__main__': globals()[sys.argv[1]]()
