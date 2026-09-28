"""Read actual native ionization inventories without modifying the EOS library."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import ctypes, json, shutil, subprocess, sys
import numpy as np
import direct_eos_gr as g

OUT=g.OUT/'species-inventory';CACHE=g.CACHE/'species-inventory'
LIB=g.d.CACHE/'full-integral-build/src/libfree_eos_direct24_integral_full.so.1.0.0'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists();OUT.mkdir();CACHE.mkdir(exist_ok=True)
    paths=[g.ROOT/'verification/eos_species_inventory.py',g.ROOT/'verification/direct_eos_gr.py',
        g.ROOT/'verification/direct_ion_eos.py',g.OUT/'reference-state.npz',g.d.OUT/'full-integral-build.json',
        g.d.OUT/'direct_ion_bridge.f90']
    for name in ['mod_nuvar.f90','ionize.f90','free_eos_detailed.f90']:
        target=OUT/name;shutil.copy2(g.d.CACHE/'full-integral-source/src'/name,target);paths.append(target)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='04bddb6',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},library_sha256=g.c.sha(LIB),
        state='outputs/direct-eos-gr33/reference-state.npz',cells=5735,processes=4,block_cells=128,
        control_cells=[0,2972,5734],inventory_relative_tolerance=1e-10,
        charge_relative_tolerance=1e-10,reported_fraction_absolute_tolerance=1e-12,
        history_bitwise_required=True,old_EOS_outputs_bitwise_required=True,
        method='Read exported mod_nuvar symbols immediately after a deterministic EOS call. A thin wrapper additionally returns the already-computed H2/H2+ fractions. Compare all 21 previously defined outputs and active species after intervening calls. Only source-defined active slots are read; derivatives in mod_nuvar are not assumed current.',
        scope='Actual same-model ion/neutral/molecular inventory on the fixed reference state. It is not the new GR solution, a nonideal Hessian certificate or a physical EOS calibration.',
        physical_EOS_certified=False,full_GR_evolution=False))


def build():
    plan=json.loads((OUT/'plan.json').read_text());assert g.c.sha(LIB)==plan['library_sha256']
    original=(g.d.OUT/'direct_ion_bridge.f90').read_text()
    text=original.replace('direct_ion_eos','ionization_inventory').replace('res(22)','res(24)')
    assert text!=original and 'res(13:) = [' in text
    text=text.replace('res(13:) = [','res(13:22) = [')
    text=text.replace('end subroutine ionization_inventory','  res(23:24) = [h2rat,h2plusrat]\nend subroutine ionization_inventory')
    assert 'res(23:24)' in text
    source=OUT/'inventory_bridge.f90';source.write_text(text)
    root=LIB.parent;module=next(root.parent.rglob('mod_nuvar.mod')).parent
    command=['gfortran','-O2','-fPIC','-shared','-I'+str(module),str(source),'-L'+str(root),
        '-Wl,-rpath,'+str(root),'-lfree_eos_direct24_integral_full','-o',str(CACHE/'inventory.so')]
    result=subprocess.run(command,capture_output=True,text=True);(OUT/'build.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr;assert g.c.sha(LIB)==plan['library_sha256']
    save('bridge.json',dict(classification='Counterexample candidate',command=command,
        bridge_sha256=g.c.sha(CACHE/'inventory.so'),source_sha256=g.c.sha(source),original_library_unchanged=True))


class InventoryEOS(g.EOS):
    def __init__(self):
        super().__init__();self.inventory_lib=ctypes.CDLL(str(CACHE/'inventory.so'))
        call=self.inventory_lib.ionization_inventory
        call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,
            np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),ctypes.POINTER(ctypes.c_int)]
        call.restype=None
        def capture(mode,value,t,eps,out,info):
            raw=np.full(24,np.nan);call(mode,value,t,eps,raw,info)
            out[:]=raw[:22];out[20]=0.;self.molecules=raw[22:].copy()
        self.call=capture

    def snapshot(self,r,t,x):
        a=self(2,r,t,x);lib=self.inventory_lib
        n=ctypes.c_int.in_dll(lib,'__mod_nuvar_MOD_nuvar_nelements').value;assert 0<=n<=24
        def integers(name): return np.ctypeslib.as_array((ctypes.c_int*24).in_dll(lib,'__mod_nuvar_MOD_'+name))[:n].copy()
        indices=integers('nuvar_index_element')-1;charges=integers('nuvar_atomic_number')
        assert len(set(indices))==n and np.all((indices>=0)&(indices<24))
        assert np.array_equal(charges,g.d.CHARGES[indices])
        raw=np.ctypeslib.as_array((ctypes.c_double*(29*24)).in_dll(lib,'__mod_nuvar_MOD_nuvar')).reshape((29,24),order='F')
        values=np.zeros((24,29));active=np.zeros((24,29),dtype=bool)
        for j,(index,z) in enumerate(zip(indices,charges)):
            values[index,:z+1]=raw[:z+1,j];active[index,:z+1]=True
        assert np.all(np.isfinite(values)) and np.all(values>=0) and np.all(np.isfinite(self.molecules))
        return dict(eos=a,number_fractions=values,active=active,molecular_H_fractions=self.molecules.copy())


def check(snapshot,x,r):
    values=snapshot['number_fractions'];ym=(x/g.c.A)@g.d.EOS().mapping
    weights=np.array(json.loads((g.d.OUT/'model-data.json').read_text())['atomic_weights'])
    cx=float(ym@weights);eps=ym/cx
    present=eps>0;missing=present&~snapshot['active'].any(1)
    counts=values.sum(1);h2,h2plus=snapshot['molecular_H_fractions']
    counts[0]+=eps[0]*(h2+h2plus)
    inventory_error=float(np.max(abs(counts[present]/eps[present]-1)))
    charge=float((values*np.arange(29)).sum()+eps[0]*h2plus/2)
    expected=snapshot['eos'][13]/(np.exp(r)*cx)
    charge_error=float(abs(charge/expected-1))
    fractions=[]
    if eps[0]>0: fractions.append(abs(values[0,1]/eps[0]-snapshot['eos'][14]))
    if eps[1]>0:
        fractions.extend([abs(values[1,1]/eps[1]-snapshot['eos'][15]),abs(values[1,2]/eps[1]-snapshot['eos'][16])])
    return dict(missing_nonzero_elements=np.flatnonzero(missing).tolist(),inventory_error=inventory_error,
        charge_error=charge_error,reported_fraction_error=float(max(fractions,default=0)),
        active_states=int(snapshot['active'].sum()),molecular_H_fractions=[float(h2),float(h2plus)])


def validate(row,plan):
    return not row['missing_nonzero_elements'] and row['inventory_error']<plan['inventory_relative_tolerance'] and row['charge_error']<plan['charge_relative_tolerance'] and row['reported_fraction_error']<plan['reported_fraction_absolute_tolerance']


def control():
    plan=json.loads((OUT/'plan.json').read_text());state=dict(np.load(g.ROOT/plan['state']));new=InventoryEOS();old=g.EOS();rows=[]
    for k,i in enumerate(plan['control_cells']):
        x=state['X'][i];r=state['lnd'][i];t=state['lnT'][i]
        a=new.snapshot(r,t,x);baseline=old(2,r,t,x)
        j=plan['control_cells'][(k+1)%3];new.snapshot(state['lnd'][j],state['lnT'][j],state['X'][j])
        repeated=new.snapshot(r,t,x)
        row=dict(cell=i,**check(a,x,r),old_outputs_bitwise=bool(np.array_equal(a['eos'],baseline)),
            history_bitwise=all(np.array_equal(a[key],repeated[key]) for key in a))
        row['passed']=validate(row,plan) and row['old_outputs_bitwise'] and row['history_bitwise'];rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',**a);print('SPECIES CONTROL',row,flush=True)
    save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows,
        physical_EOS_certified=False));assert all(r['passed'] for r in rows)


def block(start):
    plan=json.loads((OUT/'plan.json').read_text());state=dict(np.load(g.ROOT/plan['state']));stop=min(start+plan['block_cells'],len(state['X']))
    eos=InventoryEOS();records=[];arrays=[]
    for i in range(start,stop):
        s=eos.snapshot(state['lnd'][i],state['lnT'][i],state['X'][i]);row=check(s,state['X'][i],state['lnd'][i]);row['cell']=i
        row['passed']=validate(row,plan);records.append(row);arrays.append(s)
    target=OUT/f'block-{start}.npz';np.savez_compressed(target,**{k:np.array([s[k] for s in arrays]) for k in arrays[0]})
    result=dict(classification='Counterexample candidate',start=start,stop=stop,rows=records,
        passed=all(r['passed'] for r in records),output_sha256=g.c.sha(target),state_sha256=g.c.sha(g.ROOT/plan['state']),plan_sha256=g.c.sha(OUT/'plan.json'))
    save(f'block-{start}.json',result);assert result['passed'],start;return result


def run():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(LIB)==plan['library_sha256']
    assert g.c.sha(CACHE/'inventory.so')==json.loads((OUT/'bridge.json').read_text())['bridge_sha256']
    assert json.loads((OUT/'control.json').read_text())['passed']
    starts=list(range(0,plan['cells'],plan['block_cells']));records=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for done in as_completed([pool.submit(block,i) for i in starts]):
            records.append(done.result());save('progress.json',dict(completed_cells=sum(r['stop']-r['start'] for r in records),total_cells=plan['cells']))
            print('SPECIES INVENTORY',len(records),'/',len(starts),flush=True)
    records.sort(key=lambda r:r['start']);rows=[row for r in records for row in r['rows']]
    save('result.json',dict(classification='Counterexample candidate',passed=True,cells=len(rows),
        maximum_inventory_error=max(r['inventory_error'] for r in rows),maximum_charge_error=max(r['charge_error'] for r in rows),
        maximum_reported_fraction_error=max(r['reported_fraction_error'] for r in rows),
        active_states_range=[min(r['active_states'] for r in rows),max(r['active_states'] for r in rows)],
        state_sha256=g.c.sha(g.ROOT/plan['state']),physical_EOS_certified=False,nonideal_Hessian_certified=False))
    print('SPECIES INVENTORY COMPLETE',len(rows),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
