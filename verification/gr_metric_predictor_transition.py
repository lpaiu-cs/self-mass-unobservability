"""Preserve the failed original inventory before changing its numerical solver."""
from pathlib import Path
import json,os,signal,sys,time
import gr_metric_mechanical_predictor as predictor
from gr_outer_product_transition import process,snapshot,descendants,same

ROOT=predictor.ROOT;OLD=predictor.prior.OUT;OUT=OLD.parent/'gr-metric-predictor-transition';sha=predictor.sha
SCRIPT='verification/gr_metric_coupled_subcell_newton.py'


def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')


def prepare():
    assert not OUT.exists();predictor.verify_preflight();OUT.mkdir()
    paths=[Path(__file__),ROOT/'verification/gr_outer_product_transition.py',
        predictor.OUT/'preflight-manifest.json',predictor.OUT/'comparison-manifest.json',predictor.OUT/'comparison.json']
    save('plan.json',dict(classification='Counterexample candidate',bindings={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        reason='After the 72-case same-target predictor preflight, stop the exact original producer tree. Preserve every original byte, completed failure and partial block. The original full run is failed and incomplete. Prioritize executing coupled fluid/heat/metric evolution instead of launching another manufactured full5735-cell inventory. Do not reuse old successful rows as new-solver evidence.',
        original_script=SCRIPT,population=5735,physical_EOS_certified=False,full_GR_evolution=False))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert sha(ROOT/rel)==digest,rel
    return plan


def inspect():
    rows=snapshot();candidates={pid:r for pid,r in rows.items() if r['cwd']==str(ROOT.resolve())
        and len(r['command'])==3 and Path(r['command'][0]).name=='python3' and r['command'][1:]==[SCRIPT,'run']}
    roots=[pid for pid,r in candidates.items() if r['ppid'] not in candidates];assert len(roots)==1,roots
    print(json.dumps(dict(root=roots[0],processes=[rows[p] for p in sorted(descendants(roots[0],rows))]),indent=2),flush=True)


def run():
    plan=bindings();assert not (OUT/'producer-stop.json').exists()
    rows=snapshot();root_path=str(ROOT.resolve())
    candidates={pid for pid,r in rows.items() if r['cwd']==root_path and len(r['command'])==3
        and Path(r['command'][0]).name=='python3' and r['command'][1:]==[SCRIPT,'run']}
    roots=[pid for pid in candidates if rows[pid]['ppid'] not in candidates];assert len(roots)==1,roots
    root=roots[0];owned=descendants(root,rows);assert owned==candidates and 1<=len(owned)<=3
    stopped={};killed=False
    try:
        for pid in sorted(owned):
            record=rows[pid];assert same(record);os.kill(pid,signal.SIGSTOP);stopped[pid]=record
        for _ in range(100):
            if all((r:=process(pid)) is None or r['state'] in ['T','Z'] for pid in stopped):break
            time.sleep(.05)
        else:raise RuntimeError('Original inverse producer did not quiesce')
        files=[p for p in OLD.rglob('*') if p.is_file()]+[ROOT/'outputs/gr-metric-coupled-subcell-newton33-run.log']
        assert all(p.exists() for p in files)
        hashes={p.relative_to(ROOT).as_posix():sha(p) for p in files}
        completed=[];invalid=[];bad=[];exceptions=[];quadrature_failures=[]
        for folder in sorted(OLD.glob('block-*')):
            result=folder/'result.json'
            if not result.exists():continue
            try:value=json.loads(result.read_text())
            except (ValueError,UnicodeError) as error:invalid.append(dict(path=result.relative_to(ROOT).as_posix(),error=str(error)));continue
            assert len(value['cells'])==len(set(value['cells']))
            completed+=value['cells']
            exceptions+=value['failures']
            quadrature_failures += [c['cell'] for c in value['rows'] if not c['finite_quadrature_passed']]
            bad += [dict(cell=r['cell'],nodes=r['nodes'],case=r['case']) for c in value['rows'] for r in c['rows'] if not r['passed']]
        assert len(completed)==len(set(completed)) and len(bad)>=24
        assert all(0<=i<plan['population'] for i in completed)
        receipt=dict(classification='Counterexample candidate',stopped=False,partial_files_preserved=True,
            original_verdict='failed-and-incomplete',process_records=list(stopped.values()),bindings=hashes,
            completed_cells=sorted(completed),unevaluated_or_partial_cells=sorted(set(range(plan['population']))-set(completed)),
            completed_cells_definition='Slots in completed block inventories, including explicitly retained evaluation exceptions; not a list of certified states.',
            failed_inverse_cases=bad,evaluation_exceptions=exceptions,failed_quadrature_cells=quadrature_failures,
            partial_result_files=invalid,full_GR_evolution=False,physical_EOS_certified=False)
        save('producer-stop.json',receipt)
        for pid,record in sorted(stopped.items(),key=lambda item:item[0]==root):
            if same(record):os.kill(pid,signal.SIGKILL)
        killed=True
        for _ in range(100):
            if all((r:=process(pid)) is None or r['start']!=record['start'] or r['state']=='Z' for pid,record in stopped.items()):break
            time.sleep(.05)
        else:raise RuntimeError('Identified original producer is still active')
        for rel,digest in hashes.items():assert sha(ROOT/rel)==digest,rel
        receipt['stopped']=True;save('producer-stop.json',receipt)
        save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():sha(p) for p in OUT.iterdir() if p.is_file()}));verify()
    finally:
        if not killed:
            for pid,record in stopped.items():
                if same(record):os.kill(pid,signal.SIGCONT)


def verify():
    plan=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert sha(ROOT/rel)==digest,rel
    r=json.loads((OUT/'producer-stop.json').read_text());assert r['stopped'] and r['partial_files_preserved']
    for rel,digest in r['bindings'].items():assert sha(ROOT/rel)==digest,rel
    assert sorted(r['completed_cells']+r['unevaluated_or_partial_cells'])==list(range(plan['population']))
    assert r['original_verdict']=='failed-and-incomplete' and len(r['failed_inverse_cases'])>=24
    print('PASS original inverse bytes and failed/incomplete population preserved:',len(r['completed_cells']),len(r['failed_inverse_cases']),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
