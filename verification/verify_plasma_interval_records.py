"""Bind every exported residual interval to an exact-rational root-error budget."""
from fractions import Fraction as F
import json, sys
from pathlib import Path
import gr_plasma_interval as original

g=original.g;OUT=g.OUT/'gr-plasma-interval-record-audit'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();original.bindings()
    paths=[g.ROOT/'verification/verify_plasma_interval_records.py',original.OUT/'controls-manifest.json',original.OUT/'plan.json']
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='fbcfefe',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
        expected_cells=5735,expected_tests_per_cell=11,threshold='1e-10',
        rule='Wait for the original complete frozen manifest, verify runtime libraries and exact row coverage, then use max(original rational error, absolute outward exported residual endpoints) as the reported bound. Preserve any original endpoint-budget discrepancies and report additional padding. No tolerance or physical parameter changes.',
        physical_EOS_certified=False))


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    original.verify();runtime=json.loads((original.OUT/'runtime.json').read_text())
    for path,digest in runtime['libraries'].items():assert g.c.sha(Path(path))==digest,path
    assert g.c.sha(Path(runtime['binary']))==runtime['binary_sha256']
    result=json.loads((original.OUT/'result.json').read_text());rows=[];failed=0;max_padding=F(0);maximum=F(0)
    assert sorted(r['cell'] for r in result['rows'])==list(range(plan['expected_cells']))
    for cell in result['rows']:
        assert sorted((t['kind'],t['index']) for t in cell['tests'])==sorted([('T',i) for i in range(6)]+[('L',i) for i in range(5)])
        for test in cell['tests']:
            a,b=map(F,test['residual'][1:-1].split(','));assert a<=b
            before=F(test['h_error_upper_rational']);bound=max(before,abs(a),abs(b));assert bound<F(plan['threshold'])
            padding=bound-before;failed+=padding>0;max_padding=max(max_padding,padding);maximum=max(maximum,bound)
            rows.append(dict(cell=cell['cell'],kind=test['kind'],index=test['index'],h_error_upper_rational=str(bound)))
    save('result.json',dict(classification='Proven',passed=True,cells=plan['expected_cells'],root_tests=len(rows),rows=rows,
        original_exported_interval_budget_discrepancies=int(failed),maximum_export_padding_rational=str(max_padding),
        maximum_certified_h_error_rational=str(maximum),display_only=dict(maximum_h_error=float(maximum),maximum_export_padding=float(max_padding)),
        original_result_sha256=g.c.sha(original.OUT/'result.json'),original_manifest_sha256=g.c.sha(original.OUT/'manifest.json'),
        threshold_unchanged=True,physical_parameter_uncertainty_included=False,physical_EOS_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    original.verify()
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['passed'] and result['root_tests']==5735*11
    assert g.c.sha(original.OUT/'manifest.json')==result['original_manifest_sha256']
    print('PASS 63085 exported plasma root certificates with exact rational budgets',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
