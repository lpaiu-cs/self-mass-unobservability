"""Resolve executed versions without rewriting their original plans."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path('E:/lab/self-mass-unobservability') if sys.platform=='win32' else Path('/mnt/e/lab/self-mass-unobservability')
OUT=ROOT/'outputs/direct-eos-gr33/native-retained-tail'
cache={}
def sha(p):
    key=str(p.resolve())
    if key not in cache:
        h=hashlib.sha256()
        with p.open('rb') as stream:
            for b in iter(lambda:stream.read(1024*1024),b''):h.update(b)
        cache[key]=h.hexdigest()
    return cache[key]
snapshots={sha(p):p for p in OUT.rglob('*.py')}
rows=[];missing=[];external=[]
for p in sorted(OUT.rglob('*.json')):
    if 'plan' not in p.name:continue
    plan=json.loads(p.read_text())
    for name,h in plan.get('bindings',{}).items():
        rel=name.replace('/home/lpaiu/work/native-retained-tail-runtime/','').replace('/mnt/e/lab/self-mass-unobservability/','')
        q=ROOT/rel
        if name.startswith('/home/') and rel==name:
            if sys.platform=='win32':external.append(dict(path=name,sha256=h));continue
            q=Path(name)
        if not q.is_file() or sha(q)!=h:q=snapshots.get(h)
        row=dict(plan=p.relative_to(ROOT).as_posix(),declared_path=name,sha256=h)
        if q is None:missing.append(row)
        else:
            assert sha(q)==h
            row['retained_path']=str(q);rows.append(row)
result=dict(verified=rows,unavailable=missing,external_checks_required=external,
    all_historical_input_bindings_verified=not missing and not external)
(OUT/'bindings-audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(verified=len(rows),unavailable=missing,external_checks_required=len(external))))
