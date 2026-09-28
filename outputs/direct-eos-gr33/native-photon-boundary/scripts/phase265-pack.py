"""Pack the per-item queue values of the 575-time boundary and verify them against the merged file (WSL)."""
from pathlib import Path
import hashlib,json,os
import numpy as np
geom=Path('native-geometric-clock265-work')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
assert read(geom/'result.json')['passed'] and not (geom/'queue-items.npz').exists()
items=sorted(int(p.stem.split('-')[1]) for p in geom.glob('queue-[0-9]*.json'))
rows=[read(geom/f'queue-{i}.json') for i in items];zs=[np.load(geom/f'queue-{i}.npz') for i in items]
values=np.array([z['value'] for z in zs]);launch=np.array([float(z['launch']) for z in zs])
merged=np.load(geom/'boundary-575.npz')
for i,v,l in zip(items,values,launch):
    assert np.array_equal(merged['photon_geometric_mass_cm'][i],v[0]) and np.array_equal(merged['photon_geometric_lapse'][i],v[1])
    assert np.array_equal(merged['lapse_energy_radius_angle_parts'][i],v[2:]) and merged['physical_launch_energy_increment_erg'][i]==l
q=read(geom/'queue.json');assert sorted(v['index'] for v in q['items'])==items and len(q['claimed'])==len(items)
np.savez_compressed(geom/'queue-items.npz',index=np.array(items),value=values,launch=launch)
write(geom/'queue-items.json',dict(items=rows))
write(geom/'queue-items-check.json',dict(classification='Counterexample candidate',passed=True,items=len(items),
    merged_values_bitwise=True,queue_items_npz_sha256=sha(geom/'queue-items.npz'),queue_items_json_sha256=sha(geom/'queue-items.json'),
    boundary_sha256=sha(geom/'boundary-575.npz'),source_sha256=sha(__file__)))
print(json.dumps(read(geom/'queue-items-check.json')))
