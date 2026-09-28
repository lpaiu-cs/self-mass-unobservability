"""Reuse saved EOS auxiliaries and materialize the native tracer's input once.

Counterexample candidate. Keep the original 180-second limit and all numerical
checks; preserve the frozen historical tracer and the timed-out attempt.
"""
import importlib.util, json, shutil, sys, time
import numpy as np
import direct_eos_gr as g

OUT=g.OUT/'source-retry'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    source=g.ROOT/'verification/native_eos_bridge.py';text=source.read_text()
    old="inp=np.load(OUT/(label+'-input.npz'));"
    new="inp=dict(np.load(OUT/(label+'-input.npz')));"
    assert text.count(old)==1 and text.count('time.monotonic()-begin<180')==2
    shutil.copy2(source,OUT/'original-native_eos_bridge.py')
    (OUT/'materialized_native_eos_bridge.py').write_text(text.replace(old,new))
    shutil.copy2(g.ROOT/'outputs/direct-eos-gr33-source.log',OUT/'original-timeout.log')
    assert json.loads((g.OUT/'reference-auxiliary.json').read_text())['passed']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='2438b36',
        original_source_sha256=g.c.sha(source),replacement_source_sha256=g.c.sha(OUT/'materialized_native_eos_bridge.py'),
        reference_auxiliary_sha256=g.c.sha(g.OUT/'reference-auxiliary.npy'),
        reference_state_sha256=g.c.sha(g.OUT/'reference-state.npz'),
        original_failure='The reference native trace exceeded 180 seconds after the new EOS auxiliary calculation passed. The original trace repeatedly decompresses lnd, lnT and X from an NPZ object inside each cell.',
        change='Materialize those immutable input arrays once. Preserve all ABI, scientific inputs, output checks and the 180-second limit. Replay the previous initial source bitwise before using the saved new EOS inputs.',
        physical_EOS_certified=False,full_GR_evolution=False))
    path=g.OUT/'reference-input.npz';rounds=32
    with np.load(path) as archive:
        start=time.perf_counter()
        for i in range(rounds): a=[archive[k] for k in ['lnd','lnT','X']]
        repeated=time.perf_counter()-start
        data=dict(archive);start=time.perf_counter()
        for i in range(rounds): b=[data[k] for k in ['lnd','lnT','X']]
        cached=time.perf_counter()-start
    assert all(np.array_equal(x,y) for x,y in zip(a,b))
    save('input-access-control.json',dict(classification='Counterexample candidate',passed=True,
        rounds=rounds,repeated_decompression_s=repeated,materialized_access_s=cached,
        scope='Input-access measurement, not proof that every contributor to the timeout is identified.'))


def configure():
    plan=json.loads((OUT/'plan.json').read_text());path=OUT/plan.get('replacement_source_filename','materialized_native_eos_bridge.py')
    assert g.c.sha(path)==plan['replacement_source_sha256']
    assert g.c.sha(g.ROOT/'verification/native_eos_bridge.py')==plan['original_source_sha256']
    spec=importlib.util.spec_from_file_location('request33_materialized_eos_trace',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    g.s.v.shim=module;g.s.OUT=g.OUT;g.s.CACHE=g.CACHE
    # Keep g.s.v.OUT at the frozen weak-table stage, as in the original source path.
    return plan


def identity():
    plan=configure();data=dict(np.load(g.OLD/'initial-state.npz'));aux=np.load(g.OLD/'initial-replacement.npy')
    g.s.auxiliary=lambda state:aux
    raw,_=g.s.evaluate(plan.get('identity_label','cached-identity'),data);previous=dict(np.load(g.OLD/'initial-native.npz'))
    equal={key:bool(np.array_equal(value,raw[key])) for key,value in previous.items()}
    save('identity-control.json',dict(classification='Counterexample candidate',passed=all(equal.values()),bitwise_equal=equal))
    assert all(equal.values()),equal;print('MATERIALIZED TRACE IDENTITY',len(raw['X']),'bitwise equal',flush=True)


def reference():
    plan=configure();assert json.loads((OUT/'identity-control.json').read_text())['passed']
    assert g.c.sha(g.OUT/'reference-state.npz')==plan['reference_state_sha256']
    assert g.c.sha(g.OUT/'reference-auxiliary.npy')==plan['reference_auxiliary_sha256']
    data=dict(np.load(g.OUT/'reference-state.npz'));aux=np.load(g.OUT/'reference-auxiliary.npy')
    g.s.auxiliary=lambda state:aux
    raw,corrected=g.s.evaluate('reference-cached',data);assert np.array_equal(raw['aux_used'],aux)
    previous=dict(np.load(g.OLD/'initial-corrected.npz'))
    changes={name:float(np.max(abs(corrected[name]-previous[name])/np.maximum(1e-30,np.maximum(abs(corrected[name]),abs(previous[name]))))) for name in ['dxdt','heat','neutrino']}
    g.save('reference-source.json',dict(classification='Counterexample candidate',completed=True,
        label='reference-cached',new_EOS_inputs_confirmed_bitwise=True,fixed_rho_T_X=True,
        relative_model_changes=changes,native_source_sha256=g.c.sha(g.OUT/'reference-cached-native.npz'),
        corrected_source_sha256=g.c.sha(g.OUT/'reference-cached-corrected.npz'),
        original_timeout_preserved=True,timeout_seconds_unchanged=180,physical_EOS_certified=False,GR_time_evolution=False))
    print('DIRECT EOS SOURCE CACHED',len(raw['X']),changes,flush=True)


def repair_profile_read():
    plan=json.loads((OUT/'plan.json').read_text());path=OUT/'materialized_native_eos_bridge.py'
    assert g.c.sha(path)==plan['replacement_source_sha256']
    text=path.read_text();old='except (ValueError,IndexError,KeyError,OSError): pass'
    assert text.count(old)==1
    new='except (ValueError,IndexError,KeyError,OSError,StopIteration): pass'
    shutil.copy2(OUT/'plan.json',OUT/'before-profile-plan.json')
    shutil.copy2(g.ROOT/'outputs/direct-eos-gr33-source-retry-identity.log',OUT/'partial-header-failure.log')
    target=OUT/'materialized_native_eos_bridge_v2.py';target.write_text(text.replace(old,new))
    plan['replacement_source_filename']=target.name;plan['replacement_source_sha256']=g.c.sha(target)
    plan['identity_label']='cached-identity-v2'
    plan['profile_read_repair']='The faster trace reached the newly created profile before its six header lines were complete. Treat StopIteration as an incomplete poll, alongside existing incomplete-file exceptions. Keep the full-profile success condition and the same 180-second deadline.'
    save('plan.json',plan)


if __name__=='__main__': globals()[sys.argv[1]]()
