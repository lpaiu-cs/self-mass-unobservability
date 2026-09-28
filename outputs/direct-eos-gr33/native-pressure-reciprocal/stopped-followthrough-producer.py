"""One-shot continuation of the live175 process; no retry or extra sweep."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

out=Path('native-pressure-reciprocal175-work');state=out/'followthrough-state.json'
guard=out/'followthrough.lock';fd=os.open(guard,os.O_CREAT|os.O_EXCL|os.O_WRONLY);os.close(fd)
read=lambda p:json.loads(Path(p).read_text())
sources=['verification/return_native_pressure_reciprocity.py','verification/complete_native_pressure_reciprocity.py','verification/read_native_pressure_charge.py']
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
frozen={p:sha(p) for p in sources};started=time.time()
identity=read(out/'status.json');assert identity['pid']==324 and identity['process_start_ticks']=='373'
assert identity['source_sha256']==frozen[sources[0]]

def save(**values):
    tmp=state.with_suffix('.tmp');tmp.write_text(json.dumps(dict(pid=os.getpid(),started_unix=started,source_bindings=frozen,**values),indent=2));tmp.replace(state)

try:
    save(state='waiting',dependency_pid=324,dependency_start_ticks='373')
    while Path('/proc/324/stat').exists():
        assert Path('/proc/324/stat').read_text().split()[21]=='373','Dependency PID reused'
        assert time.time()-started<4800,'Dependency observation cap reached; do not restart it'
        time.sleep(30)
    assert read(out/'run-receipt.json')['error'] is None
    assert read(out/'photon-result.json')['passed'] and read(out/'status.json')['state']=='completed'
    actions=[(sources[1],a) for a in ['prepare','material_pilot','material_production','residual','block']]
    actions.extend((sources[2],a) for a in ['prepare','compact','packets','mass','audit'])
    for source,action in actions:
        assert all(sha(p)==h for p,h in frozen.items()),'Continuation source changed'
        save(state='running',source=source,action=action)
        print(json.dumps(dict(event='start',source=source,action=action,unix=time.time())),flush=True)
        subprocess.run([sys.executable,source,action],check=True)
    save(state='completed',final_charge_conclusion=read('native-pressure-charge176-work/result.json')['final_charge_conclusion'])
except BaseException as exc:
    save(state='failed',error=repr(exc));raise
