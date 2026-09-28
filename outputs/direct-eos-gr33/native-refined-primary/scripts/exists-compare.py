import json, os
root = '/home/lpaiu/work/native-retained-tail-runtime/'; new = '/home/lpaiu/work/native-refined267-runtime/'
d = json.load(open(root + '.phase267-exists-original.json'))
rows = []
for p, r in sorted(d.items()):
    if not p.startswith(root): continue
    rel = p[len(root):]
    if rel.startswith('primary267-') or rel.startswith('verification/') or rel.startswith('.'): continue
    nr = os.path.exists(new + rel)
    if r != nr: rows.append((r, nr, os.path.getsize(p) if r and os.path.isfile(p) else -1, rel))
print(len(d), 'checked;', len(rows), 'differ (original, refined, size, path):')
for row in rows: print(*row)
