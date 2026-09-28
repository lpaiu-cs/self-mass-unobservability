"""Persist the external process exit, stdout and stderr without filtering."""
from pathlib import Path
import json,subprocess,sys,time
out=Path('outputs/direct-eos-gr33/native-retained-tail/supported-temperature')
action=sys.argv[1];cap=sys.argv[2];start=time.monotonic()
with (out/(action+'-external.log')).open('w') as stream:
    p=subprocess.run(['timeout','--signal=TERM','--kill-after=5s',cap+'s',sys.executable,
        sys.argv[3] if len(sys.argv)>3 else 'verification/evolve_native_retained_tail.py',action],stdout=stream,stderr=subprocess.STDOUT)
result=dict(action=action,returncode=p.returncode,seconds=time.monotonic()-start,external_cap_seconds=float(cap))
(out/(action+'-receipt.json')).write_text(json.dumps(result,indent=2)+'\n')
print((out/(action+'-external.log')).read_text(),end='');print(json.dumps(result),flush=True)
sys.exit(p.returncode)
