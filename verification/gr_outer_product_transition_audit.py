"""Freeze coverage and preservation evidence for the stopped-to-adaptive transition."""
from fractions import Fraction as F
from pathlib import Path
import json,sys
import gr_outer_product_adaptive_table as table

ROOT=table.ROOT;sha=table.sha;OUT=table.OUT.parent/'gr-outer-product-transition-audit'


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def run():
    assert not OUT.exists();table.verify_preflight();table.original.pilot.verify();inventory=table.check_inventory()
    receipt=json.loads((table.OUT/'producer-stop.json').read_text());assert receipt['stopped'] and receipt['partial_files_preserved']
    assert sha(ROOT/'verification/gr_outer_product_transition.py')==receipt['transition_source_sha256']
    for rel,digest in receipt['bindings'].items():assert sha(ROOT/rel)==digest,rel
    rows=[];files=[Path(__file__),table.OUT/'inventory.json',table.OUT/'producer-stop.json',
        table.OUT/'preflight-manifest.json',table.original.pilot.OUT/'manifest.json']
    for record in inventory['records']:
        i=record['position'];path=ROOT/record['path'];files.append(path)
        result=json.loads(path.read_text()) if record['kind']=='original_pilot' else table.original.verify_cell(i)
        rows.append(result)
    i=inventory['preflight_position'];rows.append(table.worker().verify_cell(i));rows.sort(key=lambda r:r['position'])
    assert len(rows)==len({r['cell'] for r in rows})==len({r['position'] for r in rows})==12
    maximum=F(0)
    for r in rows:
        assert r['passed'] and r['complete_finite_interval'] and r['point_budget_passed']
        assert F(r['response_midpoint_error_upper'])<F('1e-15') and len(r['field_error_upper'])==7
        for interval,point,error in zip(r['complete_field_enclosures'],r['field_approximation'],r['field_error_upper'],strict=True):
            lo,hi=table.original.pilot.cusp.endpoints(interval)
            assert max(abs(lo-F.from_float(point)),abs(hi-F.from_float(point)))==F(error)<F('2e-7')
            maximum=max(maximum,F(error))
    OUT.mkdir();save('plan.json',dict(classification='Proven',bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        scope='Retained byte identity of every stopped original producer file and the finite complete-state inventory. No verdict is assigned to an interrupted partial state. The four new adaptive workers are outside these frozen directories.'))
    save('result.json',dict(classification='Proven',passed=True,unique_states=len(rows),physical_outputs=7*len(rows),
        actual_middle_Q_nodes=sum(r['nodes'] for r in rows),maximum_field_error_upper=str(maximum),
        completed_positions=[r['position'] for r in rows],remaining_states=len(inventory['remaining_positions']),
        preserved_original_files=len(receipt['bindings']),preserved_partial_folders=inventory['original_partial_folders'],
        all3206_complete=False,physical_EOS_certified=False,nonlinear_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert sha(ROOT/rel)==digest,rel
    receipt=json.loads((table.OUT/'producer-stop.json').read_text())
    for rel,digest in receipt['bindings'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text())
    assert r['passed'] and r['unique_states']==12 and r['physical_outputs']==84 and r['actual_middle_Q_nodes']==174592
    assert r['remaining_states']==3194 and F(r['maximum_field_error_upper'])<F('2e-7') and not r['all3206_complete']
    print('PASS transition:12 complete states,84 fields,174592 Q nodes and all89 original files retained',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
