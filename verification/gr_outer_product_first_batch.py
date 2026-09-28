"""Bind the first four scheduled full-table completions without changing the run."""
from fractions import Fraction as F
import json,sys
import gr_outer_product_table as table

ROOT=table.ROOT;OUT=table.OUT.parent/'gr-outer-product-first-batch';sha=table.sha


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def run():
    positions=[2,3,4,5]
    for i in positions:assert (table.OUT/f'cell-{i:04d}'/'manifest.json').exists()
    assert not OUT.exists();OUT.mkdir();table.bindings()
    files=[ROOT/'verification/gr_outer_product_first_batch.py',table.OUT/'plan.json',table.OUT/'candidate.py',table.OUT/'preflight-manifest.json',table.pilot.OUT/'manifest.json']
    files += [table.OUT/f'cell-{i:04d}'/'manifest.json' for i in positions]
    save('plan.json',dict(classification='Proven',positions=positions,bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        target='Audit the first four positions scheduled by the frozen whole-table producer, in original order. No node,threshold,physical model or live source is changed.',
        boundary='Completed finite set of declared ideal-electron RPA self terms only. Whole3206-state H controls/finalization, full physical EOS, GR evolution and observations remain open.'))
    rows=[table.verify_cell(i) for i in positions];table.verify_preflight();table.pilot.verify()
    previous=[json.loads((table.pilot.OUT/f'cell-{i:04d}-result.json').read_text()) for i in [0,1603,3205]]+[table.verify_cell(1)]
    allrows=sorted(previous+rows,key=lambda x:x['position']);assert len({r['cell'] for r in allrows})==8
    save('result.json',dict(classification='Proven',passed=True,positions=[r['position'] for r in allrows],states=8,outputs=56,
        new_positions=positions,new_nodes=sum(r['nodes'] for r in rows),total_nodes=sum(r['nodes'] for r in allrows),
        new_maximum_field_error=str(max(F(e) for r in rows for e in r['field_error_upper'])),
        total_maximum_field_error=str(max(F(e) for r in allrows for e in r['field_error_upper'])),rows=rows,whole_table_complete=False,physical_EOS_certified=False))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()


def verify():
    for name,key in [('plan.json','bindings'),('manifest.json','sha256')]:
        for rel,digest in json.loads((OUT/name).read_text())[key].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'result.json').read_text());assert r['passed'] and r['states']==8 and r['outputs']==56 and r['new_positions']==[2,3,4,5]
    assert F(r['total_maximum_field_error'])<F('2e-7') and not r['whole_table_complete'] and not r['physical_EOS_certified']
    print('PASS first four scheduled full-table cells; eight unique full-Q states and56 physical outputs bound at the original gate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
