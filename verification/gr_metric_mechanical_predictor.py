"""Apply the tested target-only predictor without changing conservative targets."""
from types import ModuleType
import json,subprocess,sys
import gr_metric_predictor_diagnostic as diagnostic

prior=diagnostic.prior;ROOT=prior.ROOT;OUT=prior.OUT.parent/'gr-metric-mechanical-predictor';sha=prior.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def source():
    before,_=prior.source();old,new=diagnostic.transformed();assert before.count(old)==1
    after=before.replace(old,new)
    a='row=dict(cell=index,nodes=number,case=case,solver_success=success,errors=error.tolist(),'
    b='row=dict(cell=index,nodes=number,case=case,solver_success=success,errors=error.tolist(),predictor_history=cell.predictor_history,'
    assert after.count(a)==1;after=after.replace(a,b)
    assert after.replace(b,a).replace(new,old)==before
    return after,[dict(old=old,new=new),dict(old=a,new=b)]


def prepare():
    assert not OUT.exists();diagnostic.verify()
    assert json.loads((diagnostic.OUT/'result.json').read_text())['all_passed']
    OUT.mkdir();candidate,changes=source();(OUT/'candidate.py').write_text(candidate)
    plan=json.loads((prior.OUT/'plan.json').read_text())
    plan.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        pilot_cells=[0,1,2,76,77,78,82,83,86,2688,2972,5734],substitutions=changes,
        predictor='Use the frozen diagnostic predictor from the same [0,.003,0] initial guess. Only the initial search path changes: up to eight baryon/momentum iterations at fixed initial temperature followed by the original full Newton iterations. No true root is supplied to the predictor; no target, metric/EOS equation or final gate changes. Mechanical and full residual descent are separate phase criteria.',
        scope='Post-failure numerical candidate. Its 72-case preflight includes every original representative case and every case in the six currently frozen failed cells. The 24 failed original cases remain failed in their original run. A later full run is a distinct population evaluation, not a replacement of the old formal verdict. No physical EOS certification or GR trajectory.')
    files=[ROOT/'verification/gr_metric_mechanical_predictor.py',ROOT/'verification/gr_metric_predictor_diagnostic.py',
        OUT/'candidate.py',diagnostic.OUT/'manifest.json',prior.OUT/'preflight.json',
        prior.OUT/'block-0064/result.json',prior.OUT/'block-0080/result.json']
    plan['bindings'].update({p.relative_to(ROOT).as_posix():sha(p) for p in files})
    save('plan.json',plan)
    for name in ['symbolic.json','quadrature-controls.json']:(OUT/name).write_bytes((prior.OUT/name).read_bytes())


def engine():
    source_text,_=source();assert source_text==(OUT/'candidate.py').read_text()
    name='gr_metric_mechanical_predictor_candidate';obj=ModuleType(name);obj.__file__=str(OUT/'candidate.py')
    sys.modules[name]=obj;obj.mechanical_predictor=diagnostic.mechanical_predictor
    exec(compile(source_text,obj.__file__,'exec'),obj.__dict__);obj.OUT=OUT;return obj


def archived(cell,n,case):
    if cell in [76,77,78,82,83,86]:return diagnostic.archived(cell,n,case)
    row=next(c for c in json.loads((prior.OUT/'preflight.json').read_text())['rows'] if c['cell']==cell)
    return next(r for r in row['rows'] if r['nodes']==n and r['case']==case)


def compare():
    after=json.loads((OUT/'preflight.json').read_text());assert after['passed']
    records=[]
    for cell in after['rows']:
        for row in cell['rows']:
            old=archived(row['cell'],row['nodes'],row['case'])
            assert row['target_changes_exact']==old['target_changes_exact'] and row['passed']
            assert row['predictor_history'][0]['x']==[0.,.003,0.]
            for key,merit in [('predictor_history','mechanical_merit'),('history','merit')]:
                assert all(b[merit]<a[merit] for a,b in zip(row[key],row[key][1:]))
            records.append(dict(cell=row['cell'],nodes=row['nodes'],case=row['case'],original_passed=old['passed'],passed=row['passed']))
    assert len(records)==72 and sum(not r['original_passed'] for r in records)==24
    save('comparison.json',dict(classification='Counterexample candidate',passed=True,rows=records,
        original_failed_cases=24,all_target_dyadics_identical=True,physical_EOS_certified=False,full_GR_evolution=False))
    files=[OUT/'preflight-manifest.json',OUT/'comparison.json',diagnostic.OUT/'manifest.json']
    save('comparison-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}))


def preflight():engine().preflight();compare();verify_preflight()


def verify_preflight():
    engine().verify_preflight()
    for rel,digest in json.loads((OUT/'comparison-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'comparison.json').read_text());assert r['passed'] and len(r['rows'])==72 and r['original_failed_cases']==24
    print('PASS same72 metric-coupled targets with target-only predictor; all24 original failures retained',flush=True)


def run():verify_preflight();engine().run()


def verify():verify_preflight();engine().verify()


if __name__=='__main__':globals()[sys.argv[1]]()
