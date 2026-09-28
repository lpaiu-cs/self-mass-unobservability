"""Counterexample candidate: isolate the actual 1e6 K molecular model switch.

Retaining the existing Taylor partition extension is a comparison model, not a
physical repair. The original EOS, all thresholds, and failures are preserved.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
import ctypes, json, shutil, sys
import numpy as np
import eos_molecular_stationarity as molecular

e=molecular.e;g=molecular.g;Base=e.EOS
PREVIOUS_OUT=e.OUT;PREVIOUS_CACHE=e.CACHE;PREVIOUS_NAME=e.NAME
OUT=g.OUT/'molecular-switch';CACHE=g.CACHE/'molecular-switch'
NAME='free_eos_direct24_molecular_switch'
e.OUT=OUT;e.CACHE=CACHE;e.NAME=NAME;e.LIB=CACHE/'build/src'/('lib'+NAME+'.so.1.0.0')


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def read(name): return json.loads((OUT/name).read_text())


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    source=CACHE/'source';shutil.copytree(PREVIOUS_CACHE/'source',source)
    path=source/'src/mod_free_eos.f90';raw=path.read_text()
    old='    if(tl.ge.tllim) then\n       ! turn off molecules'
    assert raw.count(old)==1 and raw.count('  public free_eos, version')==1
    assert 'tllim = log(1.e6_fp_kind)' in raw
    text=raw.replace('  public free_eos, version','  public free_eos, version\n  logical, save, public :: retain_molecules=.false.')
    text=text.replace(old,'    if(tl.ge.tllim.and..not.retain_molecules) then\n       ! turn off molecules')
    path.write_text(text)
    cmake=source/'src/CMakeLists.txt';rawcmake=cmake.read_text()
    assert rawcmake.count('OUTPUT_NAME '+PREVIOUS_NAME)==1
    cmake.write_text(rawcmake.replace('OUTPUT_NAME '+PREVIOUS_NAME,'OUTPUT_NAME '+NAME))
    for name in ['mod_free_eos.f90','CMakeLists.txt']:
        shutil.copy2(PREVIOUS_CACHE/'source/src'/name,OUT/('before-'+name))
        shutil.copy2(source/'src'/name,OUT/name)
    shutil.copy2(PREVIOUS_OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='8ab926a',
        original_cutoff_K=1e6,controls=[0,1175,1176,2972,5734],processes=4,block_cells=128,
        boundary_cells=[1175,1176],boundary_log_steps=[1e-3,1e-4,1e-5,1e-6,1e-7,1e-8],
        inventory_tolerance=1e-10,charge_tolerance=1e-10,original_outputs_bitwise_required=True,
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/eos_molecular_switch.py',PREVIOUS_OUT/'manifest.json',
            g.OUT/'reference-state.npz',g.OUT/'initial-state-17-4.npz']},
        original_library_sha256=g.c.sha(molecular.p.c.s.LIB),
        intervention='Only the final common ifh2/ifh2plus=0 temperature switch can be bypassed. All other options and partition-function Taylor continuation are unchanged. Default false must reproduce the original EOS.',
        policy='Report every finite change and failure, with no after-result smallness gate. A forced molecular continuation is not a validated physical partition function. Samples and one-sided finite sequences are not uniform bounds.',
        physical_EOS_certified=False,continuous_EOS_certified=False,full_GR_evolution=False))


class EOS(Base):
    def sample(self,r,t,x,retain):
        ctypes.c_int.in_dll(self.inventory_lib,'__mod_free_eos_MOD_retain_molecules').value=int(retain)
        snap=self.snapshot(r,t,x)
        flags=np.ctypeslib.as_array((ctypes.c_int*4).in_dll(self.inventory_lib,'__mod_nuvar_MOD_mol_flags')).copy()
        snap.update(mol_flags=flags,mol_logs=self.array('mod_nuvar','mol_logs',3))
        report=molecular.p.c.s.check(snap,x,r,self)
        assert not report['missing_nonzero_elements']
        assert report['inventory_error']<1e-10 and report['charge_error']<1e-10,report
        return snap,report


def controls():
    state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();original=g.EOS();rows=[]
    for i in read('plan.json')['controls']:
        r,t,x=state['lnd'][i],state['lnT'][i],state['X'][i]
        a,ra=eos.sample(r,t,x,False);baseline=original(2,r,t,x)
        b,rb=eos.sample(r,t,x,True);again,_=eos.sample(r,t,x,False)
        same=np.array_equal(a['eos'],baseline) and all(np.array_equal(v,again[k]) for k,v in a.items())
        assert same,i
        if t<np.log(1e6):assert all(np.array_equal(v,b[k]) for k,v in a.items()),i
        row=dict(cell=i,original_bitwise=same,original_molecules=a['molecular_H_fractions'].tolist(),
            retained_molecules=b['molecular_H_fractions'].tolist(),retained_flags=b['mol_flags'].tolist(),
            inventory_error=rb['inventory_error'],charge_error=rb['charge_error'])
        rows.append(row);print('MOLECULAR SWITCH CONTROL',row,flush=True)
    save('control.json',dict(classification='Counterexample candidate',passed=True,rows=rows))


def boundary():
    plan=read('plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[];samples=[]
    t0=np.log(1e6)
    for i in plan['boundary_cells']:
        r,x=state['lnd'][i],state['X'][i]
        for h in plan['boundary_log_steps']:
            values=[];molecules=[];flags=[]
            for retain in [False,True]:
                for sign in [-1,1]:
                    snap,_=eos.sample(r,t0+sign*h,x,retain)
                    values.append(snap['eos']);molecules.append(snap['molecular_H_fractions']);flags.append(snap['mol_flags'])
            values=np.array(values);mol=np.array(molecules)
            assert np.array_equal(values[0],values[2])
            assert flags[0][3]==1 and flags[1][3]==0 and flags[2][3]==flags[3][3]==1
            samples.append(dict(cell=i,step=h,values=values,molecules=mol))
            rows.append(dict(cell=i,log_step=h,
                original_pressure_crossing=float((values[1,1]-values[0,1])/values[0,1]),
                retained_pressure_crossing=float((values[3,1]-values[2,1])/values[2,1]),
                pressure_switch_difference=float((values[1,1]-values[3,1])/values[3,1]),
                energy_switch_difference_erg_g=float(values[1,2]-values[3,2]),
                entropy_switch_difference_erg_g_K=float(values[1,3]-values[3,3]),
                high_side_retained_molecules=mol[3].tolist()))
    np.savez_compressed(OUT/'boundary.npz',cells=np.array([s['cell'] for s in samples]),
        steps=np.array([s['step'] for s in samples]),values=np.array([s['values'] for s in samples]),
        molecules=np.array([s['molecules'] for s in samples]))
    save('boundary.json',dict(classification='Counterexample candidate',completed=True,rows=rows,
        continuous_limit_certified=False,physical_partition_extension_certified=False))
    print('MOLECULAR BOUNDARY COMPLETE',rows[-1],flush=True)


def block(start):
    plan=read('plan.json');state=dict(np.load(g.OUT/'reference-state.npz'))
    stop=min(start+plan['block_cells'],len(state['X']));eos=EOS();original=g.EOS();rows=[];arrays=[]
    for i in range(start,stop):
        r,t,x=state['lnd'][i],state['lnT'][i],state['X'][i]
        base,rb=eos.sample(r,t,x,False);assert np.array_equal(base['eos'],original(2,r,t,x)),i
        new,rn=eos.sample(r,t,x,True)
        if t<np.log(1e6):assert all(np.array_equal(v,new[k]) for k,v in base.items()),i
        rows.append(dict(cell=i,inventory_error=rn['inventory_error'],charge_error=rn['charge_error']))
        arrays.append(dict(original=base['eos'],retained=new['eos'],original_molecules=base['molecular_H_fractions'],
            retained_molecules=new['molecular_H_fractions'],retained_flags=new['mol_flags'],retained_logs=new['mol_logs']))
    target=OUT/f'block-{start}.npz';np.savez_compressed(target,**{k:np.array([a[k] for a in arrays]) for k in arrays[0]})
    result=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,
        plan_sha256=g.c.sha(OUT/'plan.json'),output_sha256=g.c.sha(target),passed=True)
    save(f'block-{start}.json',result);return result


def run():
    plan=read('plan.json');assert read('control.json')['passed']
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    starts=list(range(0,5735,plan['block_cells']));rows=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for done in as_completed([pool.submit(block,i) for i in starts]):
            rows.append(done.result());print('MOLECULAR SWITCH',len(rows),'/',len(starts),flush=True)
    parts=[dict(np.load(OUT/f'block-{i}.npz')) for i in starts]
    a={k:np.concatenate([p[k] for p in parts]) for k in parts[0]};before=a['original'];after=a['retained']
    changed=np.any(before!=after,axis=1)
    result=dict(classification='Counterexample candidate',completed=True,cells=5735,
        all_original_EOS_outputs_bitwise=True,changed_cells=int(changed.sum()),
        maximum_absolute_EOS_changes=abs(after-before).max(0).tolist(),
        maximum_scaled_EOS_changes=(abs(after-before)/np.maximum(1,abs(before))).max(0).tolist(),
        maximum_retained_molecular_fractions=a['retained_molecules'].max(0).tolist(),
        physical_partition_extension_certified=False,continuous_EOS_certified=False,physical_EOS_certified=False)
    save('result.json',result);print('MOLECULAR SWITCH COMPLETE',result,flush=True)


def verify():
    plan=read('plan.json')
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(molecular.p.c.s.LIB)==plan['original_library_sha256']
    assert read('control.json')['passed'] and read('result.json')['completed']
    seen=[]
    for path in sorted(OUT.glob('block-*.json')):
        row=json.loads(path.read_text());assert row['passed']
        assert row['output_sha256']==g.c.sha(path.with_suffix('.npz'))
        assert row['plan_sha256']==g.c.sha(OUT/'plan.json');seen.extend(range(row['start'],row['stop']))
    assert sorted(seen)==list(range(5735))
    if (OUT/'manifest.json').exists():
        for rel,digest in read('manifest.json')['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    else:
        paths=[p for p in OUT.iterdir() if p.is_file()]+[g.ROOT/'verification/eos_molecular_switch.py']
        save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},
            runtime={str(p):g.c.sha(p) for p in [e.LIB,CACHE/'excitation.so']}))
    for path,digest in read('manifest.json')['runtime'].items():assert g.c.sha(g.d.Path(path))==digest,path
    print('PASS MOLECULAR SWITCH',len(read('manifest.json')['sha256']),'artifact SHA',flush=True)


build=e.build
if __name__=='__main__':globals()[sys.argv[1]]()
