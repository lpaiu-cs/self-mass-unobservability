"""Resume only uncommitted reader actions after a verified WSL boot change."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time,shutil
job=sys.argv[1];root=Path('native-full-charge251-work')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
configs={
 '247':('native-compensated-charge247-work','full','full-','read_compensated_return_charge.py','', [('field1288',0),('field648',8),('field1284',12)],['collect','audit','compare']),
 '248':('native-retarded-extension248-work','','','extend_unchanged_source_prefix.py','', [('field1288',1),('field648',5),('field1284',9)],['collect','audit']),
 '251':('native-full-charge251-work','check-retry','','read_full_return_charge.py','check_', [('charge_field1288',4),('charge_field648',10),('charge_field1284',14)],['charge_collect','charge_audit','charge_compare','regression'])}
folder,sub,status_prefix,name,prefix,fields,tail=configs[job];folder=Path(folder);out=folder/sub
control=out if job=='251' else folder;producer=Path('verification')/name;launcher=Path('.phase251-reader-launch.py')
status_name=status_prefix+'controller-status.json';start_name=status_prefix+'controller-start.json'
old_start=read(control/start_name);old_status=read(control/status_name);boot=Path('/proc/sys/kernel/random/boot_id').read_text().strip()
assert old_start['boot_id']!=boot,'Not a verified boot-change recovery'
assert old_status['state']=='running'
archive=root/'interruption-0605'/job;archive.mkdir(parents=True)
for p in [control/start_name,control/status_name,producer,launcher]:shutil.copyfile(p,archive/p.name)
bound={str(producer):sha(producer),str(launcher):sha(launcher)}
expected=old_start.get('producer_sha256',old_start.get('bindings',{}).get(str(producer)))
assert expected==bound[str(producer)],'Producer changed since interrupted run'
def write(name,value):
 p=control/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(json.dumps(value,indent=2)+'\n');os.replace(q,p)
def receipt(action):
 if job=='251':
  if action=='regression':return root/'regression-receipt.json'
  part,verb=action.split('_',1);return out/part/(verb+'-receipt.json')
 return out/(action+'-receipt.json')
def passed(action):
 p=receipt(action)
 if not p.exists():return False
 r=read(p);assert r['error'] is None,(action,r);assert r['source_sha256']==bound[str(producer)],action
 return True
start=time.monotonic();completed=[];workers={}
write(start_name,dict(pid=os.getpid(),process_start_ticks=Path('/proc/self/stat').read_text().split()[21],boot_id=boot,
 started_unix=time.time(),controller_sha256=sha(__file__),bindings=bound,producer_sha256=bound[str(producer)],
 resumption_after_boot_change=True,interrupted_start=old_start,archive=str(archive),
 policy='Keep completed receipts and scientific gates. Resume only missing actions, preserve unfinished files/logs. Serialize shared generated-code initialization only.',
 final_charge_conclusion='unadjudicated',full_goal_complete=False))
def launch(action,cpu=6):
 for p,h in bound.items():assert sha(p)==h,p
 for suffix in ['stdout.log','stderr.log']:
  p=(out if job=='251' else folder)/(('full-' if job=='247' else '')+action+'.'+suffix)
  if p.exists():shutil.copyfile(p,archive/p.name)
 # Raw partial source files are retained before the same pending source action.
 if job=='251' and action=='dense_source':
  for p in (out/'dense/gr').glob('source-*.npz'):shutil.copyfile(p,archive/p.name)
 with (archive/(action+'.resumed.stdout.log')).open('w') as stdout,(archive/(action+'.resumed.stderr.log')).open('w') as stderr:
  return subprocess.Popen(['taskset','-c',str(cpu),sys.executable,str(launcher),str(producer),prefix+action],stdout=stdout,stderr=stderr)
def sequential(action):
 if passed(action):completed.append(dict(action=action,reused=True));return
 child=launch(action);workers[action]=child
 write(status_name,dict(state='running',action=action,child_pid=child.pid,completed=completed,elapsed_seconds=time.monotonic()-start))
 code=child.wait();completed.append(dict(action=action,returncode=code));assert code==0,(action,code)
try:
 if job=='251':
  for action in ['endpoint_prepare','endpoint_source','dense_prepare','dense_geometry','dense_source','charge_prepare','charge_polynomial']:sequential(action)
 else:
  assert passed('prepare');assert read(out/'plan.json')['final_charge_conclusion']=='unadjudicated'
  forecast=read(folder/('full-dispatch-forecast.json' if job=='247' else 'dispatch-forecast.json'));assert forecast['passed']
 workers={}
 for action,cpu in fields:
  if passed(action):completed.append(dict(action=action,reused=True))
  else:workers[action]=launch(action,cpu)
 write(status_name,dict(state='running',action='parallel_fields',workers={k:v.pid for k,v in workers.items()},completed=completed,elapsed_seconds=time.monotonic()-start))
 while any(v.poll() is None for v in workers.values()):
  failed={k:v.returncode for k,v in workers.items() if v.poll() not in [None,0]};assert not failed,failed;time.sleep(5)
 completed.extend(dict(action=k,returncode=v.returncode) for k,v in workers.items());assert not any(v.returncode for v in workers.values()),completed
 for action in tail:sequential(action)
 write(status_name,dict(state='completed',completed=completed,elapsed_seconds=time.monotonic()-start,resumed_after_boot_change=True,
  final_charge_conclusion='unadjudicated',full_goal_complete=False))
except BaseException as exc:
 for child in workers.values():
  if child.poll() is None:child.terminate()
 for child in workers.values():
  if child.poll() is None:child.wait()
 write(status_name,dict(state='failed',completed=completed,error=repr(exc),elapsed_seconds=time.monotonic()-start,final_charge_conclusion='unadjudicated',full_goal_complete=False));raise
