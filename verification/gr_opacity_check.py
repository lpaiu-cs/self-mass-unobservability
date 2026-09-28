"""Frozen opacity checks on the new EOS GR state, before inspecting its opacity.

Counterexample candidate. Native-implementation values and the repaired spline
tangent have separate controls; none is physical opacity calibration.
"""
import json, sys
import numpy as np
import opacity_tables as t

o=t.o;g=t.g;OUT=g.OUT/'gr-opacity'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    plan=json.loads((o.OUT/'plan.json').read_text())
    plan['new_GR_opacity_protocol']=dict(classification='Counterexample candidate',
        model_sha256={str(p.relative_to(g.ROOT)):g.c.sha(p) for p in [
            g.ROOT/'verification/native_opacity.py',g.ROOT/'verification/opacity_tables.py',
            g.ROOT/'verification/opacity_cubic.py',g.ROOT/'verification/gr_opacity_check.py',
            o.OUT/'tables-enriched-tables.npz']},
        state='outputs/direct-eos-gr33/initial-state-17-4.npz',
        state_file_existed_at_protocol_creation=(g.OUT/'initial-state-17-4.npz').exists(),
        source_state_not_read_by_prepare=True,
        native_value_control='Match the native Type1 values using the otherwise identical double evaluator with the documented zero conduction mixed field. This identifies implementation agreement; the buggy native derivative verdict remains failed.',
        repaired_tangent_control='For the restored mixed-field polynomial, check both returned log derivatives at natural-log steps 5e-5 and 2.5e-5, and compare the two steps. Its intentional native-value shift is reported separately.',
        value_relative_tolerance=1e-4,derivative_relative_tolerance=1e-3,derivative_floor=1,
        finite_log_steps=[5e-5,2.5e-5],
        scope='A separately frozen check at the new GR state, not a new observational population or rescue of any original failed gate.',
        physical_EOS_or_opacity_certified=False,full_GR_evolution=False)
    save('plan.json',plan)


def run():
    protocol=json.loads((OUT/'plan.json').read_text())['new_GR_opacity_protocol']
    for rel,digest in protocol['model_sha256'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((g.OUT/'initial-GR.json').read_text())['completed']
    state=dict(np.load(g.ROOT/protocol['state']))
    repaired=t.Opacity();native_model=t.Opacity(include_mixed=False)
    o.OUT=OUT;o.CACHE=g.CACHE/'gr-opacity';o.CACHE.mkdir(exist_ok=True)
    g.c.native.OUT=OUT;g.c.native.CACHE=o.CACHE
    g.c.native.setup('new-GR',state,species=g.c.NAMES,network=(g.c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
    o.trace('new-GR',inspect_internal=True)
    native=dict(np.load(OUT/'new-GR-captured.npz'));parameters=native['parameters']
    values=[];compatibility=[];scores=[];mutual=[]
    for i,p in enumerate(parameters):
        a=repaired(p);values.append(a);compatibility.append(native_model(p)[0]);directions=[]
        for j in [0,1]:
            estimates=[]
            for h in protocol['finite_log_steps']:
                pair=[]
                for sign in [-1,1]:
                    shifted=p.copy();shifted[3+j]+=sign*h/np.log(10);pair.append(repaired(shifted)[0])
                estimates.append(np.log(pair[1]/pair[0])/(2*h))
            directions.append(estimates)
        directions=np.array(directions)
        scores.append(abs(directions-a[1:,None])/np.maximum(1,abs(a[1:,None])))
        mutual.append(abs(directions[:,0]-directions[:,1])/np.maximum(1,abs(directions[:,1])))
        if i%1000==0: print('NEW GR OPACITY',i,flush=True)
    values=np.array(values);compatibility=np.array(compatibility);scores=np.array(scores);mutual=np.array(mutual)
    native_error=abs(compatibility/native['outputs'][:,0]-1)
    passed=bool(native_error.max()<protocol['value_relative_tolerance'] and scores.max()<protocol['derivative_relative_tolerance'] and mutual.max()<protocol['derivative_relative_tolerance'])
    np.savez_compressed(OUT/'evaluation.npz',values=values,native_compatible_values=compatibility,derivative_scores=scores,mutual_scores=mutual)
    save('result.json',dict(classification='Counterexample candidate',passed=passed,cells=len(values),
        new_GR_state_sha256=g.c.sha(g.ROOT/protocol['state']),
        maximum_native_implementation_value_error=float(native_error.max()),
        maximum_repaired_tangent_score=float(scores.max()),maximum_two_step_difference=float(mutual.max()),
        maximum_intentional_repaired_vs_native_value_change=float(abs(values[:,0]/native['outputs'][:,0]-1).max()),
        prior_failed_gates_preserved=True,physical_EOS_or_opacity_certified=False,full_GR_evolution=False))
    print('NEW GR OPACITY COMPLETE',passed,float(native_error.max()),float(scores.max()),float(mutual.max()),flush=True)
    assert passed,'Preserve the failed new-GR opacity protocol'


if __name__=='__main__': globals()[sys.argv[1]]()
