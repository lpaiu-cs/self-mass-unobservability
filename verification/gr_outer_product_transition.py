"""Stop only this repository's old producer after the new full preflight passes."""
from pathlib import Path
import json,os,signal,sys,time
import gr_outer_product_adaptive_table as table


def process(pid):
    try:
        directory=Path('/proc')/str(pid);fields=(directory/'stat').read_text().rsplit(')',1)[1].split()
        return dict(pid=pid,ppid=int(fields[1]),state=fields[0],start=fields[19],
            command=(directory/'cmdline').read_bytes().decode().strip('\0').split('\0'),cwd=str((directory/'cwd').resolve(strict=True)))
    except (FileNotFoundError,ProcessLookupError,PermissionError):return None


def snapshot():return {int(p.name):r for p in Path('/proc').iterdir() if p.name.isdigit() and (r:=process(int(p.name))) is not None}


def descendants(root,rows):
    owned={root}
    while True:
        extended=owned|{pid for pid,r in rows.items() if r['ppid'] in owned}
        if extended==owned:return owned
        owned=extended


def same(record):
    current=process(record['pid']);return current is not None and current['start']==record['start']


def inspect():
    rows=snapshot();candidates={pid:r for pid,r in rows.items() if r['cwd']==str(table.ROOT.resolve())
        and len(r['command'])==3 and r['command'][1:]==['verification/gr_outer_product_table.py','run']}
    roots=[pid for pid,r in candidates.items() if r['ppid'] not in candidates];assert len(roots)==1,roots
    print(json.dumps(dict(root=roots[0],processes=[rows[pid] for pid in sorted(descendants(roots[0],rows))],
        allowed_native=str(table.original.pilot.CACHE/'product')),indent=2),flush=True)


def run():
    table.verify_preflight();path=table.OUT/'producer-stop.json';assert not path.exists()
    assert not json.loads((table.original.OUT/'progress.json').read_text())['failures']
    rows=snapshot();root_path=str(table.ROOT.resolve())
    candidates={pid for pid,r in rows.items() if r['cwd']==root_path and len(r['command'])==3
        and Path(r['command'][0]).name=='python3' and r['command'][1:]==['verification/gr_outer_product_table.py','run']}
    roots=[pid for pid in candidates if rows[pid]['ppid'] not in candidates];assert len(roots)==1,roots
    root=roots[0];stopped={};killed=False
    try:
        for _ in range(100):
            rows=snapshot();owned=descendants(root,rows)
            for pid in sorted(owned):
                if pid in stopped:continue
                record=rows[pid];assert record['cwd']==root_path
                assert pid in candidates or record['command'][0]==str(table.original.pilot.CACHE/'product'),record
                assert same(record);os.kill(pid,signal.SIGSTOP);stopped[pid]=record
            current=[process(pid) for pid in stopped]
            if all(r is None or r['state'] in ['T','Z'] for r in current):
                again=snapshot()
                if descendants(root,again)<=set(stopped):break
            time.sleep(.05)
        else:raise RuntimeError('Owned process tree did not quiesce')
        progress=json.loads((table.original.OUT/'progress.json').read_text());assert not progress['failures']
        for result_path in table.original.OUT.glob('cell-*/cell-*-result.json'):
            assert json.loads(result_path.read_text())['passed'],str(result_path)
        files=[p for p in table.original.OUT.rglob('*') if p.is_file()]
        bindings={p.relative_to(table.ROOT).as_posix():table.sha(p) for p in files}
        receipt=dict(classification='Counterexample candidate',stopped=False,partial_files_preserved=True,
            reason='Switch the remaining inventory to the separately certified adaptive degree implementation after its additional full-Q preflight; do not erase or finalize the interrupted original producer.',
            process_records=list(stopped.values()),original_completed_positions=progress['completed_positions'],
            bindings=bindings,transition_source_sha256=table.sha(Path(__file__)),
            new_preflight_manifest_sha256=table.sha(table.OUT/'preflight-manifest.json'))
        table.save('producer-stop.json',receipt)
        # Every target was observed as a descendant, identity-checked and stopped.
        # Kill without resuming so no partial scientific output can change.
        for pid,record in sorted(stopped.items(),key=lambda item:item[0]==root):
            if same(record):os.kill(pid,signal.SIGKILL)
        killed=True
        for _ in range(100):
            if all((r:=process(pid)) is None or r['start']!=record['start'] or r['state']=='Z' for pid,record in stopped.items()):break
            time.sleep(.05)
        else:raise RuntimeError('An identified producer process is still active')
        for rel,digest in bindings.items():assert table.sha(table.ROOT/rel)==digest,rel
        receipt['stopped']=True;table.save('producer-stop.json',receipt)
        table.inventory()
        print('PASS old producer stopped, all original bytes retained, completed inventory frozen',flush=True)
    finally:
        if not killed:
            for pid,record in stopped.items():
                if same(record):os.kill(pid,signal.SIGCONT)


if __name__=='__main__':globals()[sys.argv[1]]()
