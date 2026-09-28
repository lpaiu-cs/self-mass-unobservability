"""Reversible CPU handoff; never infer ownership from a stale PID alone."""
from pathlib import Path
import json, os, signal, sys
from gr_outer_product_transition import snapshot, descendants, same, process

receipt=Path(__file__).with_suffix('.json')
rootpath=str(Path(__file__).resolve().parents[3])
if sys.argv[1]=='pause':
    assert not receipt.exists()
    rows=snapshot()
    owned={p for p,r in rows.items() if r['cwd']==rootpath and len(r['command'])==3
           and Path(r['command'][0]).name=='python3'
           and r['command'][1:]==['verification/gr_caloric_64.py','run']}
    roots=[p for p in owned if rows[p]['ppid'] not in owned]
    assert len(roots)==1 and descendants(roots[0],rows)==owned and len(owned)==9
    records=[rows[p] for p in sorted(owned)]
    receipt.write_text(json.dumps(dict(reason='Temporarily prioritize actual coupled evolution',processes=records),indent=2)+'\n')
    for r in records:
        assert same(r)
        os.kill(r['pid'],signal.SIGSTOP)
    print('PAUSED CALORIC',len(records),flush=True)
elif sys.argv[1]=='resume':
    saved=json.loads(receipt.read_text())
    for r in saved['processes']:
        assert same(r),r['pid']
        os.kill(r['pid'],signal.SIGCONT)
    print('RESUMED CALORIC',len(saved['processes']),flush=True)
else:
    raise ValueError(sys.argv[1])
