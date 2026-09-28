"""Counterexample candidate: common EOS inputs on the new material GR state.

Reuse the frozen native tracers and weak correction. Compute electron and
reaction auxiliaries together, in isolated Fortran processes, without changing
the historical experiments or equating a local input check with evolution.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, sys
import numpy as np
import source_retry as retry
import native_opacity as opacity

g=retry.g;OUT=g.OUT/'gr-microphysics';CACHE=g.CACHE/'gr-microphysics'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();CACHE.mkdir(exist_ok=True)
    paths=[g.ROOT/'verification'/(n+'.py') for n in [
        'gr_microphysics','direct_eos_gr','source_retry','native_opacity',
        'closure_precision','conservative_star']]
    paths += [g.OUT/'source-retry/plan.json',
        g.OUT/'source-retry/materialized_native_eos_bridge_v2.py',
        g.d.OUT/'full-integral-build.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='84eeb6c',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        state='outputs/direct-eos-gr33/initial-state-17-4.npz',
        state_existed_at_prepare=(g.OUT/'initial-state-17-4.npz').exists(),
        state_not_read_at_prepare=True,processes=8,block_cells=128,
        finite_log_steps=[5e-5,2.5e-5],finite_derivative_score_tolerance=1e-3,
        native_input_log_and_composition_tolerance=1e-12,
        serial_control_cells=[0,2972,5734],
        source_columns=['free_e_per_baryon','eta','deta_dlnT','deta_dlnrho'],
        opacity_columns=['lnfree_e_per_baryon','dlnfree_e_dlnrho','dlnfree_e_dlnT'],
        method='Nine evaluations of the same declared EOS per GR cell supply both interfaces. Preserve binary64 state inputs; quantify native import rounding. Identity replay precedes each new source intervention. Finite derivative checks are not continuous enclosures.',
        physical_EOS_certified=False,full_GR_evolution=False))


def cell(eos,r,t,x,steps):
    seen=[]
    def sample(lr,lt):
        a=eos(2,lr,lt,x);assert a[13]>0
        seen.append(a[12]);return a,np.array([a[12],np.log(a[13])-lr])
    a,value=sample(r,t);slopes=[]
    for h in steps:
        slopes.append(np.array([(sample(r+h,t)[1]-sample(r-h,t)[1])/(2*h),
            (sample(r,t+h)[1]-sample(r,t-h)[1])/(2*h)]).T)
    slopes=np.array(slopes)
    scores=np.max(abs(slopes[1]-slopes[0])/np.maximum(1,abs(slopes[1])),axis=1)
    source=np.array([a[13]/np.exp(r),a[12],slopes[1,0,1],slopes[1,0,0]])
    electron=np.r_[value[1],slopes[1,1]]
    return a,source,electron,scores,np.array([min(seen),max(seen)])


def block(start):
    plan=json.loads((OUT/'plan.json').read_text());path=g.ROOT/plan['state']
    state=dict(np.load(path));stop=min(start+plan['block_cells'],len(state['X']))
    target=OUT/f'block-{start}.npz';record=OUT/f'block-{start}.json'
    binding=g.c.sha(path);protocol=g.c.sha(OUT/'plan.json')
    if target.exists() and record.exists():
        old=json.loads(record.read_text())
        assert old['state_sha256']==binding and old['plan_sha256']==protocol
        assert old['output_sha256']==g.c.sha(target);return old
    eos=g.EOS();rows=[cell(eos,state['lnd'][i],state['lnT'][i],state['X'][i],
        plan['finite_log_steps']) for i in range(start,stop)]
    arrays={key:np.array([r[j] for r in rows]) for j,key in enumerate([
        'eos','source','electron','derivative_scores','eta_ranges'])}
    np.savez_compressed(target,**arrays);score=float(arrays['derivative_scores'].max())
    result=dict(classification='Counterexample candidate',start=start,stop=stop,
        maximum_score=score,passed=score<plan['finite_derivative_score_tolerance'],
        state_sha256=binding,plan_sha256=protocol,output_sha256=g.c.sha(target))
    save(record.name,result);assert result['passed'],result;return result


def compute():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    gr=json.loads((g.OUT/'initial-GR.json').read_text());assert gr['completed']
    path=g.ROOT/plan['state'];assert g.c.sha(path)==gr['state_sha256']
    state=dict(np.load(path));n=len(state['X']);starts=list(range(0,n,plan['block_cells']));records=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for done in as_completed([pool.submit(block,k) for k in starts]):
            records.append(done.result());save('progress.json',dict(
                completed_cells=sum(r['stop']-r['start'] for r in records),total_cells=n))
            print('GR COMMON AUXILIARIES',len(records),'/',len(starts),flush=True)
    parts=[dict(np.load(OUT/f'block-{k}.npz')) for k in starts]
    arrays={key:np.concatenate([p[key] for p in parts]) for key in parts[0]}
    eos=g.EOS()
    for i in plan['serial_control_cells']:
        row=cell(eos,state['lnd'][i],state['lnT'][i],state['X'][i],plan['finite_log_steps'])
        assert all(np.array_equal(a,arrays[k][i]) for k,a in zip(arrays,row)),i
    np.savez_compressed(OUT/'auxiliaries.npz',**arrays)
    save('auxiliaries.json',dict(classification='Counterexample candidate',passed=True,cells=n,
        state_sha256=g.c.sha(path),output_sha256=g.c.sha(OUT/'auxiliaries.npz'),
        maximum_derivative_scores=arrays['derivative_scores'].max(0).tolist(),
        independent_serial_controls_bitwise=True,physical_or_continuous_certificate=False))
    print('GR COMMON AUXILIARIES COMPLETE',n,arrays['derivative_scores'].max(0),flush=True)


def inputs():
    plan=json.loads((OUT/'plan.json').read_text());result=json.loads((OUT/'auxiliaries.json').read_text())
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert result['passed'] and g.c.sha(g.ROOT/plan['state'])==result['state_sha256']
    assert g.c.sha(OUT/'auxiliaries.npz')==result['output_sha256']
    return dict(np.load(g.ROOT/plan['state'])),dict(np.load(OUT/'auxiliaries.npz'))


def source():
    state,aux=inputs();retry.configure();g.s.OUT=OUT;g.s.CACHE=CACHE
    g.c.native.OUT=OUT;g.c.native.CACHE=CACHE;g.c.native.context()
    shim=g.s.v.shim;shim.OUT=OUT;shim.CACHE=CACHE;shim.save=save
    g.c.native.setup('baseline',state,species=g.c.NAMES,
        network=(g.c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
    shim.trace('baseline');baseline=dict(np.load(OUT/'baseline-native.npz'))
    g.s.auxiliary=lambda data:baseline['aux_used']
    identity,_=g.s.evaluate('identity',state)
    equal={key:bool(np.array_equal(value,identity[key])) for key,value in baseline.items()}
    save('identity.json',dict(classification='Counterexample candidate',passed=all(equal.values()),bitwise_equal=equal))
    assert all(equal.values())
    g.s.auxiliary=lambda data:aux['source']
    native,corrected=g.s.evaluate('common-EOS',state)
    assert np.array_equal(native['aux_used'],aux['source'])
    imported=np.column_stack([np.log(native['rho']),np.log(native['T'])])
    error=float(abs(imported-np.column_stack([state['lnd'],state['lnT']])).max())
    xerror=float(abs(native['X']-state['X']).max());assert max(error,xerror)<1e-12
    previous=dict(np.load(g.OUT/'reference-cached-corrected.npz'))
    changes={key:float(np.max(abs(corrected[key]-previous[key])/
        np.maximum(1e-30,np.maximum(abs(corrected[key]),abs(previous[key])))))
        for key in ['dxdt','heat','neutrino']}
    save('source.json',dict(classification='Counterexample candidate',passed=True,
        cells=len(state['X']),new_EOS_inputs_used_bitwise=True,native_identity_bitwise=True,
        maximum_native_import_log_error=error,maximum_native_composition_error=xerror,
        model_replacement_change_from_fixed_reference=changes,
        physical_EOS_certified=False,time_evolution=False))
    print('NEW GR SOURCE COMPLETE',len(state['X']),changes,flush=True)


def apply_opacity():
    state,aux=inputs();target=g.OUT/'gr-opacity'
    assert json.loads((target/'result.json').read_text())['passed']
    opacity.OUT=target;opacity.CACHE=g.CACHE/'gr-opacity'
    g.c.native.OUT=target;g.c.native.CACHE=opacity.CACHE
    g.c.native.setup('common-EOS',state,species=g.c.NAMES,
        network=(g.c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())
    opacity.trace('common-EOS',aux['electron'])
    actual=dict(np.load(target/'common-EOS-captured.npz'));base=dict(np.load(target/'new-GR-captured.npz'))
    assert np.array_equal(actual['used'],aux['electron'])
    assert np.array_equal(actual['parameters'][:,:5],base['parameters'][:,:5])
    assert np.array_equal(actual['X'],base['X'])
    error=float(abs(actual['parameters'][:,3:5]*np.log(10)-np.column_stack([state['lnd'],state['lnT']])).max())
    assert error<1e-12
    save('opacity.json',dict(classification='Counterexample candidate',passed=True,
        cells=len(state['X']),same_EOS_electron_inputs_used_bitwise=True,
        maximum_native_import_log_error=error,
        maximum_opacity_relative_change=float(abs(actual['outputs'][:,0]/base['outputs'][:,0]-1).max()),
        maximum_log_derivative_change=abs(actual['outputs'][:,1:]-base['outputs'][:,1:]).max(0).tolist(),
        physical_opacity_certified=False,full_GR_evolution=False))
    print('NEW GR OPACITY EOS INPUTS COMPLETE',flush=True)


def selfcheck():
    def manufactured(mode,r,t,x):
        a=np.ones(21);a[12]=r*r+3*t;a[13]=np.exp(r+.2*r-.4*t);return a
    a,s,e,score,_=cell(manufactured,.3,.2,None,[5e-5,2.5e-5])
    assert np.allclose(s,[np.exp(-.02),.69,3,.6],rtol=0,atol=1e-10)
    assert np.allclose(e,[-.02,.2,-.4],rtol=0,atol=1e-10) and score.max()<1e-9
    print('PASS shared EOS auxiliary column and derivative control',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
