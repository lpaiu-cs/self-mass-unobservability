"""Apply the existing same-return source and charge readers to the full period."""
from pathlib import Path
import fcntl,inspect,json,os,resource,sys,time
import numpy as np
import read_compensated_return_charge as charge

endpoint=charge.prior.prior;dense=charge.prior;base=charge.base
read,write,sha,bind=charge.read,charge.write,charge.sha,charge.bind
CHECK=len(sys.argv)>1 and sys.argv[1].startswith('check_')
ROOT=Path('native-full-charge251-work');OUT=ROOT/('check-retry' if CHECK else 'full')
ACTUAL=Path('native-stage-metric227-work' if CHECK else 'native-full-return249-work/full')
PRIMARY=Path('native-complete-radau224-check-work' if CHECK else 'native-full-captured244-work')
HIGH=Path('native-stage-metric227-work' if CHECK else 'native-retarded-extension248-work')
POLYNOMIAL=Path('.phase247-polynomial-audit.py')
CAPS=dict(endpoint_prepare=600,endpoint_source=5400,dense_prepare=600,dense_geometry=3600,
    dense_source=7200,charge_prepare=600,charge_polynomial=600,charge_field1288=10800,
    charge_field648=10800,charge_field1284=10800,charge_collect=300,charge_audit=900,
    charge_compare=300,regression=600)

# Change paths in this process only. The running producers remain immutable.
endpoint.INPUT=ACTUAL;endpoint.OUT=OUT/'endpoint'
dense.INPUT=endpoint.OUT;dense.OUT=OUT/'dense';dense.PRIMARY=PRIMARY
charge.INPUT=dense.OUT;charge.OUT=OUT/'charge';charge.HIGH=HIGH
charge.Response.setup=bind(endpoint.returned.Response.setup,INPUT=dense.OUT,bridge=charge.bridge)

# The original reader fixes its old horizon. Keep its complete provenance and
# all scientific gates, changing only the declared full-period step counts.
source=inspect.getsource(endpoint.prepare)
source=base.replace(source,'[111,215]','[119,231]')
ns=dict(endpoint.prepare.__globals__);exec(compile(source,__file__,'exec'),ns)
endpoint_prepare=ns['prepare']

# High fields live in248; the actual coupled stage receipts belong to249.
source=inspect.getsource(charge.prepare)
assert source.count("HIGH/f'run-{n}.json'")==2
source=source.replace("HIGH/f'run-{n}.json'","ACTUAL/f'run-{n}.json'")
ns=dict(charge.prepare.__globals__,ACTUAL=ACTUAL);exec(compile(source,__file__,'exec'),ns)
charge_prepare=ns['prepare']


def polynomial():
    s=POLYNOMIAL.read_text()
    s=base.replace(s,"folder=Path('native-dense-returned246-work')/('check' if check else 'full')",f'folder=Path({str(dense.OUT)!r})')
    s=base.replace(s,"out=Path('native-compensated-charge247-work')/('check' if check else 'full')",f'out=Path({str(charge.OUT)!r})')
    exec(compile(s,str(POLYNOMIAL),'exec'),dict(__file__=str(POLYNOMIAL)))


def regression():
    assert CHECK;rows=[];files=[Path(__file__),POLYNOMIAL]
    pairs=[(endpoint.OUT,Path('native-returned-source245-work/check'),[f'gr/endpoint-{n}.npz' for n in [64,128]]),
        (dense.OUT,Path('native-dense-returned246-work/check'),['geometry.npz']+[f'gr/source-{n}.npz' for n in [64,128]]),
        (charge.OUT,Path('native-compensated-charge247-work/check'),[f'gr/fields-{n}-g{q}.npz' for n,q in charge.FIELDS.values()])]
    for actual,old,names in pairs:
        for name in names:
            with np.load(actual/name) as a,np.load(old/name) as b:
                assert set(a.files)==set(b.files),(name,'keys')
                for k in a.files:assert np.array_equal(a[k],b[k]),(name,k)
                rows.append(dict(file=str(actual/name),reference=str(old/name),exact_arrays=len(a.files)))
            files += [actual/name,old/name]
    r=read(charge.OUT/'result.json');old=read(Path('native-compensated-charge247-work/check/result.json'))
    assert r['charge_comparison_admitted'] and r['components']==old['components']
    assert read(charge.OUT/'polynomial-audit.json')['passed']
    write(ROOT/'regression.json',dict(classification='Counterexample candidate',passed=True,
        rows=rows,componentwise_charge_exact=True,new_physical_steps=0,
        bindings={str(p):sha(p) for p in files},
        scope='Full reader route reproduces the completed short same-return endpoints, geometry, source, three GR fields and charge. No physical or full-period error certificate.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False))


def field(action):
    n,q=charge.FIELDS[action];worker=charge.OUT/action
    # Legacy initialization writes retained-motion-return151-work generated
    # code. Serialize that shared side effect, then run the fields in parallel.
    with (ROOT/'initialization.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        bind(base.endpoint.initialize,OUT=worker)()
    fn=bind(base.base.gr.base.Response.run,OUT=worker/'gr');fn(charge.Response(),n,q)


def execute(action):
    if action=='regression':return regression()
    part,verb=action.split('_',1);module=dict(endpoint=endpoint,dense=dense,charge=charge)[part]
    if verb!='prepare':
        for p,h in read(module.OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
    elif not CHECK:
        assert read(ROOT/'regression.json')['passed']
        assert read(ROOT/'regression-receipt.json')['source_sha256']==sha(__file__)
    if action=='endpoint_prepare':endpoint_prepare()
    elif action=='charge_prepare':charge_prepare()
    elif action=='charge_polynomial':polynomial()
    elif part=='charge' and verb in charge.FIELDS:field(verb)
    else:getattr(module,verb)()
    if verb=='prepare':
        r=read(module.OUT/'plan.json');r.update(
            claim='Read the identical actual full-period249returned gas, photons, boundary and applied metric through the previously validated245/246/247operator.' if not CHECK else 'Reproduce the completed short actual227source and charge through the new path adapter.',
            reuse='No physical integration. Full: actual249return,244primary mass-constraint source and248high GR. Check: unchanged227/224short inputs.',
            scope='One compensated return; sampled Hermite representation. Full EOS, uniform derivative/time, spatial/boundary, nonlinear/selfGR, static/observational and infinity-normalization obligations remain.',
            forecast='Endpoint about8minutes; geometry about3minutes; dense source about12minutes from measured245/246long receipts. Fields use the measured247output-count times source-cut scaling with a2x allowance plus5minutes. Caps90/60/120minutes for endpoint/geometry/source and3hours each field;16GiB each. Later marginal cost is not measured.',
            budgets=CAPS,actual_solution_directory=str(ACTUAL),primary_source_directory=str(PRIMARY),high_field_directory=str(HIGH),
            full_declared_period=not CHECK,source_adapter_sha256=sha(__file__),
            final_charge_conclusion='unadjudicated',full_goal_complete=False)
        r['bindings'][str(POLYNOMIAL)]=sha(POLYNOMIAL);write(module.OUT/'plan.json',r)
    if action=='charge_compare':
        r=read(charge.OUT/'result.json');r.update(full_declared_period=not CHECK,
            actual_solution_directory=str(ACTUAL),primary_source_directory=str(PRIMARY),high_field_directory=str(HIGH))
        write(charge.OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1].removeprefix('check_');assert action in CAPS
    OUT.mkdir(parents=True,exist_ok=True)
    part,verb=action.split('_',1) if '_' in action else ('','')
    folder=dict(endpoint=endpoint.OUT,dense=dense.OUT,charge=charge.OUT).get(part,ROOT)
    receipt=folder/(verb+'-receipt.json' if verb else 'regression-receipt.json');assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,16*1024**3))
    base.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(CAPS[action]);start=time.monotonic();error=None
    try:execute(action)
    except BaseException as exc:error=repr(exc);raise
    finally:
        if folder.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
