"""Remove the unintended cumulative log box; retain the failed frozen run."""
from types import ModuleType
import json,sys
import gr_metric_coupled_subcell as original

ROOT=original.ROOT;OUT=original.OUT.parent/'gr-metric-coupled-subcell-newton';sha=original.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def source():
    before=(ROOT/'verification/gr_metric_coupled_subcell.py').read_text()
    old='if max(abs(proposed[:2]))>.1 or abs(proposed[2])>=.5:continue'
    new='if abs(proposed[2])>=.5:continue'
    assert before.count(old)==1;after=before.replace(old,new)
    assert after.replace(new,old)==before;return after,dict(old=old,new=new)


def prepare():
    assert not OUT.exists();plan=original.bindings()
    failed=json.loads((original.OUT/'preflight.json').read_text());assert not failed['passed'] and not failed['failures']
    assert all(c['all_Jacobian_controls_passed'] and c['finite_quadrature_passed'] for c in failed['rows'])
    for rel,digest in json.loads((original.OUT/'preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    bad=[r for c in failed['rows'] for r in c['rows'] if not r['passed']]
    assert bad and all(not r['solver_success'] and max(abs(x) for x in r['history'][-1]['x'][:2])>.099 for r in bad)
    OUT.mkdir();text,change=source();(OUT/'candidate.py').write_text(text)
    plan.update(classification='Counterexample candidate',checkpoint='0beea1d2',substitution=change,
        original_failed_preflight=True,original_failed_cases=len(bad),
        correction='The original per-step log caps were also applied as an unintended total |log shift|<=0.1 box. The manufactured root is inside the box, but the recorded residual-decreasing Newton path stalls at its boundary. Remove only that cumulative box. Keep all per-step caps, the velocity interval, positive-EOS/metric requirements, targets, initial guess, Jacobian, quadratures and scientific acceptance gates unchanged.',
        scope='Separate post-failure numerical repair for the same metric-coupled cell moments. This does not relax a physical EOS domain, certify global root uniqueness, evolve Q, integrate a GR time path, solve an atmosphere or complete observations.')
    files=[ROOT/'verification/gr_metric_coupled_subcell_newton.py',OUT/'candidate.py',
        original.OUT/'preflight.json',original.OUT/'preflight-manifest.json']
    plan['bindings'].update({p.relative_to(ROOT).as_posix():sha(p) for p in files});save('plan.json',plan)
    for name in ['symbolic.json','quadrature-controls.json']:(OUT/name).write_bytes((original.OUT/name).read_bytes())


def engine():
    text,_=source();assert text==(OUT/'candidate.py').read_text()
    name='gr_metric_coupled_subcell_newton_candidate';obj=ModuleType(name);obj.__file__=str(OUT/'candidate.py')
    sys.modules[name]=obj;exec(compile(text,obj.__file__,'exec'),obj.__dict__);obj.OUT=OUT;return obj


def comparison():
    before=json.loads((original.OUT/'preflight.json').read_text());after=json.loads((OUT/'preflight.json').read_text())
    assert not before['passed'] and after['passed'];a=[r for c in before['rows'] for r in c['rows']];b=[r for c in after['rows'] for r in c['rows']]
    assert len(a)==len(b)==36
    for x,y in zip(a,b):
        assert (x['cell'],x['nodes'],x['case'],x['target_changes_exact'])==(y['cell'],y['nodes'],y['case'],y['target_changes_exact'])
        assert y['passed'];assert all(h['merit']>j['merit'] for h,j in zip(y['history'],y['history'][1:]))
    maximum=max(abs(h['x'][1]) for r in b for h in r['history'])
    assert maximum>.1
    save('repair-comparison.json',dict(classification='Counterexample candidate',passed=True,cases=36,
        original_failed_cases=sum(not r['passed'] for r in a),repaired_passed_cases=sum(r['passed'] for r in b),
        all_target_dyadics_identical=True,maximum_temporary_logT_shift=maximum,
        all_scientific_gates_unchanged=True,full_GR_evolution=False))
    files=[OUT/'preflight-manifest.json',OUT/'repair-comparison.json',original.OUT/'preflight-manifest.json']
    save('repair-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}));verify_preflight()


def preflight():engine().preflight();comparison()


def verify_preflight():
    engine().verify_preflight()
    for rel,digest in json.loads((OUT/'repair-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'repair-comparison.json').read_text());assert r['passed'] and r['cases']==36 and r['all_target_dyadics_identical']
    print('PASS same36 metric-coupled targets without the unintended cumulative box; original failure retained',flush=True)


def run():verify_preflight();engine().run()


def verify():verify_preflight();engine().verify()


if __name__=='__main__':globals()[sys.argv[1]]()
