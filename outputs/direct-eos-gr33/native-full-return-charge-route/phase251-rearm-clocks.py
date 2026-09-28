from pathlib import Path
import hashlib,json,os,shutil
root=Path('native-full-charge251-work/clock-failure');root.mkdir()
read=lambda p:json.loads(Path(p).read_text())
for folder,names in [('native-full-return249-work',['controller-start.json','controller-status.json','prepare-receipt.json','prepare.stdout.log','prepare.stderr.log']),
 ('native-full-charge251-work/full',['controller-start.json','controller-status.json'])]:
 p=Path(folder);entry=read(p/'controller-start.json');state=read(p/'controller-status.json');assert state['state']=='failed'
 proc=Path('/proc')/str(entry['pid']);assert not proc.exists() or proc.joinpath('stat').read_text().split()[21]!=entry['process_start_ticks']
 assert state['completed']==([dict(action='prepare',returncode=1)] if '249' in folder else [])
 dst=root/folder.replace('/','-');dst.mkdir()
 for name in names:os.rename(p/name,dst/name)
for name in ['verification/complete_returned_period.py','.phase249-followthrough.py']:
 shutil.copyfile(name,root/Path(name).name)
