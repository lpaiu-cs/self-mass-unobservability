"""Continue the finite 3206-state inventory using the certified degree bound."""
from concurrent.futures import ProcessPoolExecutor,as_completed
from fractions import Fraction as F
from types import FunctionType,SimpleNamespace
import gzip,json,subprocess,sys
import gr_outer_product_table as original
import gr_outer_product_adaptive_pilot as pilot
import gr_outer_path_comparison as paths

ROOT=original.ROOT;sha=original.sha
OUT=original.OUT.parent/'gr-outer-product-adaptive-table'
CACHE=original.CACHE.parent/'outer-product-adaptive-table'


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();pilot.verify();original.verify_preflight()
    OUT.mkdir();CACHE.mkdir();text,changes=original.source();(OUT/'candidate.py').write_text(text)
    files=[ROOT/'verification/gr_outer_product_adaptive_table.py',ROOT/'verification/gr_outer_product_table.py',
        OUT/'candidate.py',original.OUT/'plan.json',original.OUT/'preflight-manifest.json',
        pilot.OUT/'manifest.json',ROOT/'verification/gr_outer_product_adaptive_pilot.py',
        ROOT/'verification/gr_response_product_adaptive.py',pilot.adaptive.OUT/'manifest.json',
        pilot.adaptive.OUT/'runtime.json']
    plan=original.bindings();bound=dict(plan['bindings']);bound.update(pilot.bindings()['bindings']);bound.update({p.relative_to(ROOT).as_posix():sha(p) for p in files})
    plan.update(checkpoint=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        classification='Counterexample candidate',bindings=bound,
        runtime_change='Reuse the exact original evaluate and verify_cell code objects with only OUT,CACHE,bindings and the adaptive native runtime changed. The original scientific candidate source and exact schedule are unchanged.',
        preflight_position=1,processes=4,
        preflight_gate='One additional full-Q state outside the three pilot states, with byte-exact Q replay and seven physical-field intersections with the original full-degree result. All original 1e-15 response and 2e-7 field gates retained.',
        reuse='After successful preflight and stopping the identified original producer, freeze all completed original per-cell manifests. Reuse only independently verified completed states. Retain every original partial file; incomplete cells are computed in this separate directory. Any recorded scientific failure prevents transition.',
        final_gate='All 3206 unique positions and cells, original per-state gates, and the completed H-table certificate and independent H controls. No missing or failed state may be skipped.',
        boundary='Specified finite ideal-electron RPA self contribution. No complete correlated or partially ionized physical EOS, self-consistent GR trajectory, atmosphere or observation inference.')
    save('plan.json',plan)


def bindings():
    plan=json.loads((OUT/'plan.json').read_text());text,_=original.source();assert text==(OUT/'candidate.py').read_text()
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return plan


def worker():
    facade=SimpleNamespace(**dict(original.pilot.__dict__,engine=pilot.adaptive.engine))
    ns=dict(original.__dict__,OUT=OUT,CACHE=CACHE,bindings=bindings,pilot=facade)
    for name in ['evaluate','verify_cell']:
        fn=getattr(original,name);ns[name]=FunctionType(fn.__code__,ns,name,fn.__defaults__)
    return SimpleNamespace(evaluate=ns['evaluate'],verify_cell=ns['verify_cell'])


def evaluate(position):return worker().evaluate(position)


def preflight():
    plan=bindings();i=plan['preflight_position'];new=evaluate(i);old=original.verify_cell(i)
    with gzip.open(original.OUT/f'cell-{i:04d}'/f'cell-{i:04d}-q.tsv.gz','rb') as a,gzip.open(OUT/f'cell-{i:04d}'/f'cell-{i:04d}-q.tsv.gz','rb') as b:
        assert a.read()==b.read()
    comparison=paths.compare(new,old);assert comparison['passed']
    save('preflight-comparison.json',dict(classification='Counterexample candidate',passed=True,
        exact_Q_byte_replay=True,comparison=comparison,physical_EOS_certified=False))
    files=[OUT/'plan.json',OUT/'candidate.py',OUT/'preflight-comparison.json',OUT/f'cell-{i:04d}'/'manifest.json']
    save('preflight-manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}));verify_preflight()


def verify_preflight():
    plan=bindings()
    for rel,digest in json.loads((OUT/'preflight-manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    result=worker().verify_cell(plan['preflight_position']);control=json.loads((OUT/'preflight-comparison.json').read_text())
    assert control['passed'] and control['exact_Q_byte_replay'] and control['comparison']['passed']
    print('PASS additional adaptive full-Q preflight',result['position'],result['nodes'],flush=True)


def inventory():
    """Read a separately frozen stop receipt, then bind the completed inventory."""
    verify_preflight();assert not (OUT/'inventory.json').exists()
    receipt=json.loads((OUT/'producer-stop.json').read_text());assert receipt['stopped'] and receipt['partial_files_preserved']
    progress=json.loads((original.OUT/'progress.json').read_text());assert not progress['failures']
    plan=bindings();records=[]
    for i in plan['pilot_positions']:
        path=original.pilot.OUT/f'cell-{i:04d}-result.json';r=json.loads(path.read_text());assert r['passed']
        records.append(dict(position=i,kind='original_pilot',path=path.relative_to(ROOT).as_posix(),sha256=sha(path)))
    for folder in sorted(original.OUT.glob('cell-*')):
        if not (folder/'manifest.json').exists():continue
        i=int(folder.name.split('-')[1]);r=original.verify_cell(i)
        if i==plan['preflight_position']:continue
        path=folder/'manifest.json';records.append(dict(position=i,kind='original_table',path=path.relative_to(ROOT).as_posix(),sha256=sha(path)))
    assert len({r['position'] for r in records})==len(records)
    covered={r['position'] for r in records}|{plan['preflight_position']};assert set(progress['completed_positions'])<=covered
    partial=[f.name for f in sorted(original.OUT.glob('cell-*')) if not (f/'manifest.json').exists()]
    files=[OUT/'producer-stop.json',original.OUT/'progress.json']
    save('inventory.json',dict(classification='Proven',records=records,original_partial_folders=partial,
        bindings={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        remaining_positions=[i for i in plan['positions'] if i not in covered],preflight_position=plan['preflight_position']))


def check_inventory():
    inv=json.loads((OUT/'inventory.json').read_text())
    for rel,digest in inv['bindings'].items():assert sha(ROOT/rel)==digest,rel
    for r in inv['records']:assert sha(ROOT/r['path'])==r['sha256'],r['path']
    positions=[r['position'] for r in inv['records']]+[inv['preflight_position']]+inv['remaining_positions']
    assert sorted(positions)==bindings()['positions'];return inv


def run():
    verify_preflight();inv=check_inventory();completed=sorted([r['position'] for r in inv['records']]+[inv['preflight_position']]);failures=[]
    with ProcessPoolExecutor(max_workers=bindings()['processes']) as pool:
        jobs={pool.submit(evaluate,i):i for i in inv['remaining_positions']}
        for future in as_completed(jobs):
            i=jobs[future]
            try:result=future.result();completed.append(i);print('ADAPTIVE FULL PRODUCT',i,len(completed),result['passed'],flush=True)
            except Exception as error:failures.append(dict(position=i,type=type(error).__name__,message=str(error)));print('ADAPTIVE FULL PRODUCT FAILURE',i,repr(error),flush=True)
            save('progress.json',dict(classification='Counterexample candidate',completed_positions=sorted(completed),failures=failures))
    assert not failures


def all_results():
    inv=check_inventory();records={r['position']:r for r in inv['records']};rows=[]
    for i in bindings()['positions']:
        if i in records:
            record=records[i]
            result=json.loads((ROOT/record['path']).read_text()) if record['kind']=='original_pilot' else original.verify_cell(i)
        else:result=worker().verify_cell(i)
        assert result['position']==i and result['passed'];rows.append(result)
    assert len({r['cell'] for r in rows})==3206;return rows


def finalize():
    verify_preflight();original.H_table.verify();rows=all_results()
    save('result.json',dict(classification='Proven',passed=True,cells=3206,rows=rows,
        maximum_field_error_upper=str(max(F(e) for r in rows for e in r['field_error_upper'])),physical_EOS_certified=False))
    files=[OUT/'plan.json',OUT/'candidate.py',OUT/'inventory.json',OUT/'producer-stop.json',OUT/'result.json',
        OUT/'preflight-manifest.json',original.H_table.OUT/'manifest.json']+list(OUT.glob('cell-*/manifest.json'))
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in files}));verify()


def verify():
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    result=json.loads((OUT/'result.json').read_text());assert result['passed'] and result['cells']==3206
    assert F(result['maximum_field_error_upper'])<F('2e-7') and result['rows']==all_results()
    print('PASS all3206 full-Q states with bounded-degree continuation; full physical EOS, GR and observations remain separate',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
