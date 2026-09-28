"""Reevaluate the current star with both molecular spectra before any GR refit."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import json, sys
import numpy as np
import eos_molecular_spectral as model

g=model.g;OUT=g.OUT/'gr-molecular-reference'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();model.verify()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='226f237',
        cells=5735,block_size=128,workers=3,inventory_tolerance=1e-10,charge_tolerance=1e-10,
        intervention='Hold every current baryon cell mass, nuclear fraction, baryon density and temperature fixed; reevaluate both the original GR EOS and the new two-spectrum/chemical-anchor EOS. The new EOS defines new reference entropies for a subsequent GR solve. Old geometry is an initial guess only; no gravitational mass or baryon mass fitting occurs here.',
        comparison='Pair old and new evaluations at the same saved rho,T,X. Do not demand bit equality to the saved old pressure-inverse energy, which followed a different rounded inverse path. No fitted normalization or after-result smallness gate for model differences.',
        scope='Full-grid finite EOS/source connection and conserved inventory check, not a GR solution, time evolution, physical EOS error certificate or observational inference.',
        runtime=json.loads((model.OUT/'manifest.json').read_text())['runtime'],
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [g.ROOT/'verification/gr_molecular_reference.py',
            model.OUT/'manifest.json',g.OUT/'initial-state-17-4.npz']}))


def bindings():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in plan['runtime'].items():assert g.c.sha(path)==digest,path
    return plan


def block(start):
    plan=bindings();state=dict(np.load(g.OUT/'initial-state-17-4.npz'));eos=model.EOS();old=g.EOS()
    stop=min(start+plan['block_size'],plan['cells']);rows=[];values=[];previous=[];molecules=[]
    assert not (OUT/f'block-{start:04}.json').exists()
    for i in range(start,stop):
        r,t,X=state['lnd'][i],state['lnT'][i],state['X'][i]
        snap=eos.snapshot(r,t,X);report=model.previous.switch.molecular.p.c.s.check(snap,X,r,eos)
        a=snap['eos'];assert a.shape==(22,)
        a=np.delete(a,20);b=old(2,r,t,X);assert a.shape==b.shape==(21,)
        passed=not report['missing_nonzero_elements'] and report['inventory_error']<plan['inventory_tolerance'] and report['charge_error']<plan['charge_tolerance'] and a[0]>0 and a[1]>0 and a[10]>0
        rows.append(dict(cell=i,**report,passed=bool(passed)));values.append(a);previous.append(b)
        molecules.append(snap['molecular_H_fractions'])
    target=OUT/f'block-{start:04}.npz'
    np.savez_compressed(target,cells=np.arange(start,stop),new=np.array(values),old=np.array(previous),molecules=np.array(molecules))
    record=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,all_passed=all(r['passed'] for r in rows),
        plan_sha256=g.c.sha(OUT/'plan.json'),output_sha256=g.c.sha(target))
    save(f'block-{start:04}.json',record);print('NEW MOLECULAR REFERENCE',start,stop,record['all_passed'],flush=True)
    return record


def run():
    plan=bindings();starts=list(range(0,plan['cells'],plan['block_size']))
    assert not any(OUT.glob('block-*.json')),'preserve all started outputs'
    with ProcessPoolExecutor(max_workers=plan['workers']) as pool:
        for done in as_completed([pool.submit(block,i) for i in starts]):done.result()
    arrays=[dict(np.load(OUT/f'block-{i:04}.npz')) for i in starts]
    a=np.concatenate([p['new'] for p in arrays]);b=np.concatenate([p['old'] for p in arrays])
    old=dict(np.load(g.OUT/'initial-state-17-4.npz'))
    new={**old,'logP':np.log(a[:,1]),'u_W':a[:,2],'s_B':a[:,3]}
    np.savez_compressed(OUT/'reference-state.npz',**new)
    assert all(np.array_equal(new[k],old[k]) for k in ['dm','X','lnd','lnT','CX'])
    records=[json.loads((OUT/f'block-{i:04}.json').read_text()) for i in starts]
    save('result.json',dict(classification='Counterexample candidate',completed=True,cells=len(a),
        all_finite_inventory_gates_passed=all(r['all_passed'] for r in records),
        maximum_absolute_EOS_change=abs(a-b).max(0).tolist(),
        maximum_scaled_EOS_change=(abs(a-b)/np.maximum(abs(b),1)).max(0).tolist(),
        same_baryon_nuclear_density_temperature_inputs=True,already_GR_solution=False,
        physical_EOS_certified=False,full_GR_evolution=False))
    save('manifest.json',dict(sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.iterdir() if p.is_file()}))
    verify()


def verify():
    plan=bindings()
    for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    cells=[]
    for path in sorted(OUT.glob('block-*.json')):
        r=json.loads(path.read_text());assert r['plan_sha256']==g.c.sha(OUT/'plan.json')
        assert r['output_sha256']==g.c.sha(path.with_suffix('.npz'));cells.extend(x['cell'] for x in r['rows'])
    assert cells==list(range(plan['cells']))
    assert json.loads((OUT/'result.json').read_text())['completed']
    print('PASS new molecular EOS whole-star reference bindings; geometry has not been refitted',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
