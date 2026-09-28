from pathlib import Path
import json,os
work=Path.cwd().resolve();out=(work/'native-full-return249-work').resolve()
archive=(work/'native-full-charge251-work/clock-failure/existing-seed-attempt').resolve()
assert work in out.parents and work in archive.parents
read=lambda p:json.loads(p.read_text())
entry=read(out/'controller-start.json');state=read(out/'controller-status.json')
assert state['state']=='failed' and state['completed']==[dict(action='prepare',returncode=1)]
proc=Path('/proc')/str(entry['pid']);assert not proc.exists() or proc.joinpath('stat').read_text().split()[21]!=entry['process_start_ticks']
assert not (out/'full/plan.json').exists() and not (out/'metric-receipt.json').exists()
archive.mkdir()
for name in ['controller-start.json','controller-status.json','prepare-receipt.json','prepare.stdout.log','prepare.stderr.log','full']:
 src=(out/name).resolve();dst=(archive/name).resolve();assert out in src.parents and archive in dst.parents;os.rename(src,dst)
