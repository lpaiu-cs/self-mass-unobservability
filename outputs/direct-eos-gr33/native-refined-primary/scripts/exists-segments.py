"""Every existence check made by the refined segments must agree with the same relative path in the original runtime."""
import glob, json, os
new = '/home/lpaiu/work/native-refined267-runtime/'; orig = '/home/lpaiu/work/native-retained-tail-runtime/'
paths = {}
for f in sorted(glob.glob(new + '.phase267-refined64-logs/exists-seg-*.json')) + glob.glob(new + '.phase267-readout-refined64-*-record.json'):
    d = json.load(open(f)); d = d.get('checked', d)
    for p, r in d.items(): paths.setdefault(p, set()).add(r)
bad = []
for p, rs in sorted(paths.items()):
    if not p.startswith(new): continue
    rel = p[len(new):]
    if rel.startswith(('primary267-', 'readout267-', 'verification/', '.')): continue
    o = os.path.exists(orig + rel)
    if rs != {o}: bad.append((sorted(rs), o, rel))
print(len(paths), 'paths checked;', len(bad), 'disagree with the original runtime')
for row in bad: print(*row)
