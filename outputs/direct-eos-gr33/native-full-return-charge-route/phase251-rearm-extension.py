from pathlib import Path
import hashlib,json,os,shutil
root=Path('native-full-charge251-work/extension-failure');root.mkdir()
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
bindings={}
for folder,names in [
 ('native-retarded-extension248-work',['controller-start.json','controller-status.json','prepare-receipt.json','prepare.stdout.log','prepare.stderr.log']),
 ('native-full-return249-work',['controller-start.json','controller-status.json'])]:
 p=Path(folder);entry=read(p/'controller-start.json');assert read(p/'controller-status.json')['state']=='failed'
 proc=Path('/proc')/str(entry['pid'])
 assert not proc.exists() or proc.joinpath('stat').read_text().split()[21]!=entry['process_start_ticks'],'Original controller still live'
 assert read(p/'controller-status.json')['completed']==([] if '249' in folder else [dict(action='prepare',returncode=1)])
 dst=root/folder;dst.mkdir()
 for name in names:
  src=p/name;bindings[str(src)]=sha(src);os.rename(src,dst/name);assert sha(dst/name)==bindings[str(src)]
for name in ['verification/extend_retarded_history.py','.phase248-followthrough.py','.phase249-followthrough.py']:
 p=Path(name);shutil.copyfile(p,root/p.name);bindings[name]=sha(p)
(root/'preserved.json').write_text(json.dumps(dict(classification='Counterexample candidate',
 original_failure_preserved=True,controllers_verified_terminal=True,new_physical_steps=0,bindings=bindings),indent=2)+'\n')
