"""Read current ionization log weights from an isolated, output-controlled build."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import ctypes, json, shutil, sys
import numpy as np
import eos_species_inventory as s

g=s.g;OUT=g.OUT/'current-mask';CACHE=g.CACHE/'current-mask'
NAME='free_eos_direct24_current_mask'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def read(path): return json.loads(path.read_text())


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    original=g.d.CACHE/'full-integral-source';source=CACHE/'source';shutil.copytree(original,source)
    changed={}
    def replace(rel,pairs):
        path=source/rel;before=path.read_text();text=before
        for old,new in pairs:
            assert text.count(old)==1,(rel,old,text.count(old));text=text.replace(old,new)
        (OUT/('before-'+path.name)).write_text(before);path.write_text(text)
        shutil.copy2(path,OUT/path.name);changed[rel]=dict(before=g.c.sha(original/rel),after=g.c.sha(path))
    replace('src/mod_nuvar.f90', [('  implicit none', '''  implicit none
  public audit_logweight, audit_zero, audit_reused, audit_ifnr, audit_calls
  real(fp_kind), save :: audit_logweight(29,24)=0._fp_kind
  integer, save :: audit_zero(29,24)=-1, audit_reused=0, audit_ifnr=-1, audit_calls=0''')])
    replace('src/ionize.f90', [
        ('  use mod_nuvar, only:', '  use mod_nuvar, only: audit_logweight, audit_zero, audit_reused, audit_ifnr, audit_calls\n  use mod_nuvar, only:'),
        ('  n_partial_elements = size(partial_elements) - 2', '''  audit_logweight=0._fp_kind
  audit_zero=-1
  audit_reused=merge(1,0,ifsame_under)
  audit_ifnr=ifnr
  audit_calls=audit_calls+1
  n_partial_elements = size(partial_elements) - 2'''),
        ('        arg = -ln_fractmax - dvzero(ielement)', '        arg = -ln_fractmax - dvzero(ielement)\n        audit_logweight(1,index_element)=arg'),
        ('        logsum0 = arg', '        audit_zero(1,index_element)=merge(1,0,ifneutral_zero(ielement))\n        logsum0 = arg'),
        ('           arg = fract(jndex_f) - ln_fractmax', '           arg = fract(jndex_f) - ln_fractmax\n           audit_logweight(jndex_f+1,index_element)=arg'),
        ('           if(fract(jndex_f).gt.fractmax) then', '           audit_zero(jndex_f+1,index_element)=merge(1,0,ifion_zero(ion0+jndex_f))\n           if(fract(jndex_f).gt.fractmax) then')])
    replace('src/CMakeLists.txt',[('OUTPUT_NAME free_eos_direct24_integral_full','OUTPUT_NAME '+NAME)])
    shutil.copy2(s.OUT/'inventory_bridge.f90',OUT/'direct_ion_bridge.f90')
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='d9c3c1b',
        cells=5735,processes=4,block_cells=128,control_cells=[0,2972,5734],
        log_population_absolute_tolerance=1e-10,output_bitwise_required=True,
        source_changes=changed,bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in
            [g.ROOT/'verification/eos_current_mask.py',s.OUT/'manifest.json',g.OUT/'reference-state.npz',OUT/'direct_ion_bridge.f90']},
        original_library_sha256=g.c.sha(s.LIB),
        measured='Final ionize call log weights before exponentiation and masks after applying fresh or reused decisions. Compare reconstructed retained atomic softmax with native populations, normalized within each atomic element; hydrogen molecular rescaling cancels in this conditional atomic distribution.',
        controls='The isolated instrumentation must reproduce all original 21 EOS outputs, atomic populations, active stages and molecule fractions bitwise, including an intervening state call and parallel blocks.',
        scope='Current frozen log weights and masks only. Their exponentials and current omitted-to-retained ratio can be enclosed separately. This does not enclose native log-weight calculation, molecular equilibrium or self-consistent nonideal root derivatives.',
        physical_EOS_certified=False,full_GR_evolution=False))


def build():
    old=g.d.OUT
    try:
        g.d.OUT=OUT
        g.d.build_at(CACHE/'source',CACHE/'build',NAME,CACHE/'inventory.so')
    finally: g.d.OUT=old
    assert g.c.sha(s.LIB)==read(OUT/'plan.json')['original_library_sha256']


class EOS(s.InventoryEOS):
    def __init__(self):
        super().__init__();self.inventory_lib=ctypes.CDLL(str(CACHE/'inventory.so'))
        call=self.inventory_lib.ionization_inventory
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
        call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_int)];call.restype=None
        def capture(mode,value,t,eps,out,info):
            raw=np.full(24,np.nan);call(mode,value,t,eps,raw,info)
            out[:]=raw[:22];out[20]=0.;self.molecules=raw[22:].copy()
        self.call=capture

    def snapshot(self,r,t,x):
        result=super().snapshot(r,t,x);lib=self.inventory_lib;prefix='__mod_nuvar_MOD_'
        n=ctypes.c_int.in_dll(lib,prefix+'nuvar_nelements').value
        indices=np.ctypeslib.as_array((ctypes.c_int*24).in_dll(lib,prefix+'nuvar_index_element'))[:n]-1
        def raw(name,dtype): return np.ctypeslib.as_array((dtype*(29*24)).in_dll(lib,prefix+name)).reshape((24,29)).copy()
        logs=np.zeros((24,29));mask=np.zeros((24,29),dtype=bool)
        a=raw('audit_logweight',ctypes.c_double);b=raw('audit_zero',ctypes.c_int)
        for j,index in enumerate(indices):
            active=result['active'][index];assert np.all((b[j,active]==0)|(b[j,active]==1))
            logs[index,active]=a[j,active];mask[index,active]=b[j,active]==1
        result.update(logweights=logs,mask=mask,
            reused=np.int64(ctypes.c_int.in_dll(lib,prefix+'audit_reused').value),
            ifnr=np.int64(ctypes.c_int.in_dll(lib,prefix+'audit_ifnr').value))
        return result


def logsum(v):
    v=np.asarray(v,dtype=np.longdouble);m=v.max();return m+np.log(np.exp(v-m).sum())


def check(a):
    error=0.;ratios=[];mismatches=0
    for j,active in enumerate(a['active']):
        if not active.any(): continue
        w=a['logweights'][j,active];mask=a['mask'][j,active];n=a['number_fractions'][j,active]
        assert np.all(np.isfinite(w)) and (~mask).any()
        mismatches+=int(np.count_nonzero(mask!=(n==0)))
        positive=n>0
        actual=np.log(n[positive].astype(np.longdouble)/np.sum(n,dtype=np.longdouble))
        expected=w[positive]-logsum(w[~mask]);error=max(error,float(abs(actual-expected).max()))
        if mask.any(): ratios.append(float(logsum(w[mask])-logsum(w[~mask])))
    return dict(log_population_error=error,mask_zero_mismatches=mismatches,
        maximum_log_omitted_retained_ratio=max(ratios) if ratios else None,
        reused=int(a['reused']),ifnr=int(a['ifnr']))


def gates(row,plan): return row['mask_zero_mismatches']==0 and row['log_population_error']<plan['log_population_absolute_tolerance']


def control():
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[]
    for k,i in enumerate(plan['control_cells']):
        args=[state['lnd'][i],state['lnT'][i],state['X'][i]];a=eos.snapshot(*args)
        old=dict(np.load(s.OUT/f'control-{i}.npz'))
        same=all(np.array_equal(a[key],old[key]) for key in old)
        j=plan['control_cells'][(k+1)%3];eos.snapshot(state['lnd'][j],state['lnT'][j],state['X'][j]);again=eos.snapshot(*args)
        history=all(np.array_equal(a[key],again[key]) for key in a)
        row=dict(cell=i,**check(a),original_outputs_bitwise=same,history_bitwise=history)
        row['passed']=same and history and gates(row,plan);rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',**a);print('CURRENT MASK CONTROL',row,flush=True)
    save('control.json',dict(classification='Counterexample candidate',rows=rows,passed=all(r['passed'] for r in rows)))
    assert all(r['passed'] for r in rows)


def block(start):
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();samples=[];rows=[]
    stop=min(start+plan['block_cells'],plan['cells']);old=dict(np.load(s.OUT/f'block-{start}.npz'))
    for i in range(start,stop):
        a=eos.snapshot(state['lnd'][i],state['lnT'][i],state['X'][i]);row=dict(cell=i,**check(a))
        row['original_outputs_bitwise']=all(np.array_equal(a[k],old[k][i-start]) for k in old)
        row['passed']=row['original_outputs_bitwise'] and gates(row,plan);rows.append(row);samples.append(a)
    path=OUT/f'block-{start}.npz';np.savez_compressed(path,**{k:np.array([a[k] for a in samples]) for k in samples[0]})
    result=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,passed=all(r['passed'] for r in rows),
        output_sha256=g.c.sha(path),plan_sha256=g.c.sha(OUT/'plan.json'))
    save(f'block-{start}.json',result);assert result['passed'],start;return result


def run():
    plan=read(OUT/'plan.json');assert read(OUT/'control.json')['passed']
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in read(OUT/'build.json')['sha256'].items(): assert g.c.sha(s.Path(path))==digest,path
    rows=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for future in as_completed([pool.submit(block,i) for i in range(0,plan['cells'],plan['block_cells'])]):
            rows+=future.result()['rows'];print('CURRENT MASK',len(rows),'/',plan['cells'],flush=True)
    ratios=[r['maximum_log_omitted_retained_ratio'] for r in rows if r['maximum_log_omitted_retained_ratio'] is not None]
    save('result.json',dict(classification='Counterexample candidate',passed=True,cells=len(rows),
        original_EOS_and_species_outputs_bitwise=True,maximum_log_population_error=max(r['log_population_error'] for r in rows),
        mask_zero_mismatches=sum(r['mask_zero_mismatches'] for r in rows),reused_mask_cells=sum(r['reused'] for r in rows),
        ifnr_counts={str(k):sum(r['ifnr']==k for r in rows) for k in sorted(set(r['ifnr'] for r in rows))},
        maximum_log_omitted_retained_ratio=max(ratios),current_ratio_interval_certified=False,physical_EOS_certified=False))
    print('CURRENT MASK COMPLETE',len(rows),'max log ratio',max(ratios),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
