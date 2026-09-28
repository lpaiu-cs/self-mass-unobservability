"""Record the complete bounded process output and terminal status."""
from pathlib import Path
import json,subprocess,sys,time
out=Path('outputs/direct-eos-gr33/native-retained-completion')
action=sys.argv[1];cap={'evolve':1500,'charge':150}[action];start=time.monotonic()
assert not (out/(action+'-receipt.json')).exists()
with (out/(action+'-external.log')).open('w') as stream:
    p=subprocess.run(['timeout','--signal=TERM','--kill-after=5s',str(cap)+'s',sys.executable,
        'verification/continue_native_retained_tail.py',action],stdout=stream,stderr=subprocess.STDOUT)
result=dict(action=action,returncode=p.returncode,seconds=time.monotonic()-start,external_cap_seconds=cap)
(out/(action+'-receipt.json')).write_text(json.dumps(result,indent=2)+'\n')
print((out/(action+'-external.log')).read_text(),end='');print(json.dumps(result),flush=True)
sys.exit(p.returncode)
