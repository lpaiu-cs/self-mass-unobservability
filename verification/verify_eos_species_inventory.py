"""Independent full-state inventory replay against the prior EOS table."""
from pathlib import Path
import json, sys
import numpy as np
import eos_species_inventory as s


def read(path): return json.loads(path.read_text())


def run():
    plan=read(s.OUT/'plan.json');result=read(s.OUT/'result.json');assert result['passed']
    for rel,digest in plan['bindings'].items(): assert s.g.c.sha(s.g.ROOT/rel)==digest,rel
    for name,digest in plan['mapping_reuse']['preserved'].items(): assert s.g.c.sha(s.OUT/name)==digest,name
    assert s.g.c.sha(s.LIB)==plan['library_sha256']
    assert s.g.c.sha(s.CACHE/'inventory.so')==read(s.OUT/'bridge.json')['bridge_sha256']
    state=dict(np.load(s.g.ROOT/plan['state']));table=dict(np.load(s.g.OUT/'initial-adiabats-17.npz'))
    assert s.g.c.sha(s.g.OUT/'initial-adiabats-17.npz')==read(s.g.OUT/'table-control.json')['output_sha256']
    assert all(np.array_equal(state[k],table[k]) for k in ['X','lnT','lnd'])
    records=[];parts=[];stop=0
    for start in range(0,plan['cells'],plan['block_cells']):
        row=read(s.OUT/f'block-{start}.json');path=s.OUT/f'block-{start}.npz'
        assert row['passed'] and row['start']==stop and row['state_sha256']==s.g.c.sha(s.g.ROOT/plan['state'])
        assert row['plan_sha256']==s.g.c.sha(s.OUT/'plan.json') and row['output_sha256']==s.g.c.sha(path)
        stop=row['stop'];parts.append(dict(np.load(path)));records.append(row)
    assert stop==len(state['X'])==plan['cells']
    data={k:np.concatenate([p[k] for p in parts]) for k in parts[0]}
    assert np.array_equal(data['eos'],table['reference']),'All 21 outputs must reproduce the independently saved EOS table'
    for i in plan['control_cells']:
        control=dict(np.load(s.OUT/f'control-{i}.npz'));assert all(np.array_equal(data[k][i],control[k]) for k in control)
    X=state['X'].astype(np.longdouble)
    # Group nuclear counts directly by charge, without the driver's map.
    ym=np.column_stack([np.sum(X[:,s.g.c.Z==z]/s.g.c.A[s.g.c.Z==z],axis=1) for z in s.g.d.CHARGES])
    weights=np.array(read(s.g.d.OUT/'model-data.json')['atomic_weights'],dtype=np.longdouble)
    cx=ym@weights;counts=data['number_fractions'].astype(np.longdouble)*cx[:,None,None]
    active=data['active'];assert active.dtype==bool and np.all(np.isfinite(counts)) and np.all(counts>=0)
    assert np.all(counts[~active]==0)
    for i,z in enumerate(s.g.d.CHARGES):
        assert not active[:,i,z+1:].any()
        assert np.all(active[:,i,:z+1]==active[:,i,0:1])
    present=ym>0;assert np.all(active.any(2)[present])
    molecules=data['molecular_H_fractions'].astype(np.longdouble)
    assert np.all(molecules>=0) and np.all(molecules<=1)
    actual=counts.sum(2);actual[:,0]+=ym[:,0]*molecules.sum(1)
    inventory=np.zeros_like(actual);inventory[present]=abs(actual[present]/ym[present]-1)
    charge=(counts*np.arange(29)).sum((1,2))+ym[:,0]*molecules[:,1]/2
    # Use the stored native rho_B to avoid importing an independently rounded exp.
    expected=data['eos'][:,13].astype(np.longdouble)/data['eos'][:,0].astype(np.longdouble)
    charge_error=abs(charge/expected-1)
    assert inventory.max()<plan['inventory_relative_tolerance'] and charge_error.max()<plan['charge_relative_tolerance']
    fractions=np.zeros_like(counts);np.divide(counts,ym[:,:,None],out=fractions,where=ym[:,:,None]>0)
    checked=np.stack([fractions[:,0,1],fractions[:,1,1],fractions[:,1,2]],axis=1)
    exported=data['eos'][:,[14,15,16]].astype(np.longdouble)
    valid=np.column_stack([present[:,0],present[:,1],present[:,1]])
    ferr=abs(checked-exported);assert ferr[valid].max()<plan['reported_fraction_absolute_tolerance']
    defined=active&present[:,:,None];zeros=defined&(counts==0);positive=fractions[defined&(counts>0)]
    report=dict(classification='Counterexample candidate',passed=True,cells=len(X),
        all_21_outputs_bitwise_equal_prior_independent_table=True,
        original_serial_controls_bitwise_equal_parallel_census=True,
        independent_longdouble_inventory_relative_max=float(inventory.max()),
        independent_longdouble_charge_relative_max=float(charge_error.max()),
        independent_reported_fraction_absolute_max=float(ferr[valid].max()),
        defined_stage_entries=int(defined.sum()),exact_zero_defined_stage_entries=int(zeros.sum()),
        cells_with_zero_defined_stage=int(zeros.any((1,2)).sum()),minimum_positive_stage_fraction=float(positive.min()),
        maximum_H2_and_H2plus_nuclear_fractions=np.max(molecules,axis=0).astype(float).tolist(),
        boundary='A defined native stage can be numerically zero. Positive-interior entropy Hessian formulas cannot use these zeros as certified physical zero abundances; logarithmic populations or a separately enclosed boundary/underflow treatment are required.',
        native_species_fractions_available=True,nonideal_Hessian_certified=False,physical_EOS_certified=False,
        table_sha256=s.g.c.sha(s.g.OUT/'initial-adiabats-17.npz'))
    assert not (s.OUT/'audit.json').exists();s.save('audit.json',report)
    paths=[p for p in s.OUT.rglob('*') if p.is_file()]+[s.g.ROOT/'verification/eos_species_inventory.py',
        s.g.ROOT/'verification/verify_eos_species_inventory.py']
    s.save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(s.g.ROOT).as_posix():s.g.c.sha(p) for p in paths},
        state_sha256=s.g.c.sha(s.g.ROOT/plan['state']),library_sha256=plan['library_sha256'],physical_EOS_certified=False))
    print('PASS FULL SPECIES AUDIT',len(X),'cells; 21 outputs bitwise; zero defined stages',int(zeros.sum()),
        '; minimum positive fraction',float(positive.min()),flush=True)


def verify():
    manifest=read(s.OUT/'manifest.json')
    for rel,digest in manifest['sha256'].items(): assert s.g.c.sha(s.g.ROOT/rel)==digest,rel
    assert s.g.c.sha(s.LIB)==manifest['library_sha256']
    print('PASS SPECIES INVENTORY',len(manifest['sha256']),'SHA bindings',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
