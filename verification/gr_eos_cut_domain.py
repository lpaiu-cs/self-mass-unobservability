"""Locate the actual native temperature switch and classify saved affine paths."""
import json,sys
import numpy as np
import gr_material_thermo_continuation as continuation

g=continuation.g;switch=continuation.switch;OUT=g.OUT/'gr-eos-cut-domain'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def run():
    assert not OUT.exists();OUT.mkdir();continuation.verify()
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));eos=switch.EOS();i=1175
    center=float(np.log(1e6));temperatures=[center]
    left=right=center
    for _ in range(8):
        left=float(np.nextafter(left,-np.inf));right=float(np.nextafter(right,np.inf));temperatures.extend([left,right])
    probes=[]
    for t in sorted(temperatures):
        snap,_=eos.sample(state['lnd'][i],t,state['X'][i],False)
        probes.append(dict(lnT=t,lnT_hex=t.hex(),active=bool(snap['mol_flags'][3])))
    cold=max(r['lnT'] for r in probes if r['active']);hot=min(r['lnT'] for r in probes if not r['active'])
    assert np.nextafter(cold,np.inf)==hot
    assert all(r['active']==(r['lnT']<hot) for r in probes)
    rows=[];paths=[];base=state['lnT'].astype(np.longdouble)
    for folder,steps in [('gr-caloric-increment',1),('gr-caloric-increment',2),('gr-caloric-increment',4),('gr-caloric-refinement',8)]:
        count=0;minimum=np.inf;crossed=[]
        for step in range(steps):
            path=g.OUT/folder/f'endpoint-{steps}-{step}.npz';data=dict(np.load(path));paths.append(path)
            low=base+data['old_shift'];high=low+data['step_shift'].astype(np.longdouble)
            low64=np.asarray(low,float);high64=np.asarray(high,float)
            crossing=(low64<hot)!=(high64<hot);indices=np.flatnonzero(crossing)
            count+=len(indices);crossed.extend(dict(step=step,cell=int(j)) for j in indices)
            minimum=min(minimum,float(np.min(abs(low-np.longdouble(hot)))),float(np.min(abs(high-np.longdouble(hot)))))
        rows.append(dict(folder=folder,steps=steps,checked_affine_segments=steps*len(base),native_switch_crossings=count,
            crossed=crossed,minimum_endpoint_logT_distance=minimum))
    table_path=g.OUT/'initial-adiabats-17.npz';table=np.load(table_path)['values'][:,:,1]
    selected=np.flatnonzero((table.min(axis=1)<hot)&(table.max(axis=1)>=hot))
    save('result.json',dict(classification='Proven',completed=True,native_flag_probes=probes,
        last_active_binary64_lnT_hex=cold.hex(),first_inactive_binary64_lnT_hex=hot.hex(),rows=rows,
        reference_pressure_tables_straddling_temperature_cut=selected.tolist(),
        scope='Exact comparisons of saved binary temperatures and monotone-rounded affine caloric line paths with the observed flag threshold. Pressure-table ranges spanning the threshold are not a smooth native EOS certificate. Other native state-dependent switches, true time paths, continuum errors and physical validity are not covered.'))
    sources=[g.ROOT/'verification/gr_eos_cut_domain.py',g.OUT/'initial-state-17-4.npz',table_path,
        continuation.OUT/'manifest.json',OUT/'result.json',*paths]
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in sources}))
    verify();print('EOS CUT DOMAIN',rows,'pressure tables',selected.tolist(),flush=True)


def verify():
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS saved affine paths and observed native cut threshold bindings',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
