"""Connect the sourced 423-level sum to an isolated FreeEOS candidate."""
from types import FunctionType
import ctypes, json, shutil, sys
import numpy as np
import h2plus_mold_data as data
import gr_material_thermo_continuation as thermo

g=data.g;switch=data.base.switch
OUT=g.OUT/'eos-h2plus-spectral';CACHE=g.CACHE/'h2plus-spectral'
NAME='free_eos_direct24_h2plus_spectral';LIB=CACHE/'build/src'/('lib'+NAME+'.so.1.0.0')
THERMO=OUT/'material-increments'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir();data.verify();thermo.verify()
    original=switch.CACHE/'source';source=CACHE/'source';shutil.copytree(original,source)
    rows=json.loads((data.OUT/'levels.json').read_text())['rows'];assert len(rows)==423
    energy=np.array([float(r['excitation_cm_inverse']) for r in rows])
    kelvin=energy*(6.62607015e-34*299792458*100/1.380649e-23)
    weights=np.array([(2*r['N']+1)*(1 if r['N']%2==0 else 3)/2 for r in rows])
    def array(name,values):
        return '    real(fp_kind), parameter :: '+name+'(423)=[ &\n'+', &\n'.join('      '+format(v,'.17e')+'_fp_kind' for v in values)+' ]\n'
    declaration=array('spectral_kelvin',kelvin)+array('spectral_weight',weights)+'''    real(fp_kind) :: sx(423), sw(423), sumw, meanx
'''
    tail='''    ! Fixed 423-level ground-electronic-state comparison; no plasma occupation error bound.
    if(ifh2plus.eq.2) then
       sx=spectral_kelvin/exp(tl)
       sw=spectral_weight*exp(-sx)
       sumw=sum(sw)
       meanx=sum(sw*sx)/sumw
       qh2plus=log(sumw)
       qh2plust=meanx
       qh2plustt=sum(sw*(sx-meanx)**2)/sumw-meanx
    endif
'''
    changes={};replacements={
        'src/mod_molecular_hydrogen.f90':[
            ('    real(fp_kind), intent(out) :: qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt',
             '    real(fp_kind), intent(out) :: qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt\n'+declaration),
            ('  end subroutine molecular_hydrogen',tail+'  end subroutine molecular_hydrogen')],
        'src/CMakeLists.txt':[('OUTPUT_NAME '+switch.NAME,'OUTPUT_NAME '+NAME)]}
    for rel,pairs in replacements.items():
        path=source/rel;before=path.read_text();after=before
        for old,new in pairs:assert after.count(old)==1;after=after.replace(old,new)
        reverse=after
        for old,new in reversed(pairs):assert reverse.count(new)==1;reverse=reverse.replace(new,old)
        assert reverse==before;path.write_text(after)
        shutil.copy2(original/rel,OUT/('before-'+path.name));shutil.copy2(path,OUT/path.name)
        changes[rel]=dict(before=g.c.sha(original/rel),after=g.c.sha(path),substitutions=pairs)
    shutil.copy2(switch.OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    np.savez_compressed(OUT/'compiled-level-inputs.npz',kelvin=kelvin,weights=weights)
    save('source-tree.json',dict(before={p.relative_to(original).as_posix():g.c.sha(p) for p in original.rglob('*') if p.is_file()},
        after={p.relative_to(source).as_posix():g.c.sha(p) for p in source.rglob('*') if p.is_file()}))
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='ddfe050',changes=changes,
        partition_temperatures_K=[1000,2000,3150,5040,8400,8999,9000,9001,12600,16800,25200,100000,999999,1000000,1000001,32000000],
        native_partition_absolute_tolerance=1e-10,controls=[0,1175,1176,2972,5734],
        inventory_tolerance=1e-10,charge_tolerance=1e-10,
        intervention='One shared molecular_hydrogen routine supplies lnQ, DlnQ, D2lnQ for H2+ option 2 using the 423 sourced levels. All three callers (eos_calc, excitation_sum, excitation_pi) retain the same common provider. H2 partition, atomic/excitation/nonideal physics and zero-point binding constants are unchanged. The existing retain_molecules switch is set True once in the new candidate; no within-library partition switching or cache-invalidating flag is introduced.',
        boundary='A sourced ground-electronic-state partition substitution, not a complete physical EOS repair. H2 still uses the previously rejected high-T Taylor extension. Physical level uncertainties, plasma occupation and binding-energy consistency with the original chemical zero are not certified. Existing full GR/reference/time paths and native libraries are not replaced.',
        controls_scope='Direct partition replay, original H2 bit equality, finite density/temperature coupled mixture increments and inventory conservation. These are finite implementation checks, not continuous/native/libm/root/physical error bounds.',
        runtime_original=switch.read('manifest.json')['runtime'],
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
            g.ROOT/'verification/eos_h2plus_spectral.py',data.OUT/'manifest.json',switch.OUT/'manifest.json',
            OUT/'source-tree.json',OUT/'compiled-level-inputs.npz',OUT/'direct_ion_bridge.f90',
            g.OUT/'initial-state-17-4.npz']}))


def build():
    fn=switch.e.build
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,CACHE=CACHE,NAME=NAME,LIB=LIB),closure=fn.__closure__)()


class EOS(switch.EOS):
    def __init__(self):
        fn=switch.EOS.__init__
        FunctionType(fn.__code__,dict(fn.__globals__,CACHE=CACHE),closure=fn.__closure__)(self)
        ctypes.c_int.in_dll(self.inventory_lib,'__mod_free_eos_MOD_retain_molecules').value=1


def partition_controls():
    plan=json.loads((OUT/'plan.json').read_text());eos=EOS();fn=eos.inventory_lib.__mod_molecular_hydrogen_MOD_molecular_hydrogen
    fn.restype=None;fn.argtypes=[ctypes.POINTER(ctypes.c_int)]*3+[ctypes.POINTER(ctypes.c_double)]*7
    old=data.base.native_function();a=dict(np.load(OUT/'compiled-level-inputs.npz'));rows=[]
    for T in plan['partition_temperatures_K']:
        ints=[ctypes.c_int(i) for i in [0,3,2]];t=ctypes.c_double(np.log(float(T)));values=[ctypes.c_double() for _ in range(6)]
        fn(*[ctypes.byref(i) for i in ints],ctypes.byref(t),*[ctypes.byref(v) for v in values]);v=np.array([x.value for x in values]).reshape(2,3)
        x=a['kelvin']/np.exp(t.value);w=a['weights']*np.exp(-x);p=w/w.sum();mean=p@x
        expected=np.array([np.log(w.sum()),mean,p@((x-mean)**2)-mean]);error=float(abs(v[1]-expected).max())
        same=np.array_equal(v[0],old(T)[0]);passed=same and error<plan['native_partition_absolute_tolerance']
        rows.append(dict(T_K=T,H2_bitwise_original=same,H2plus_jets=v[1].tolist(),maximum_absolute_replay_error=error,passed=passed))
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));inventories=[]
    for i in plan['controls']:
        r,t,X=state['lnd'][i],state['lnT'][i],state['X'][i];snap=eos.snapshot(r,t,X)
        report=switch.molecular.p.c.s.check(snap,X,r,eos)
        passed=not report['missing_nonzero_elements'] and report['inventory_error']<plan['inventory_tolerance'] and report['charge_error']<plan['charge_tolerance']
        inventories.append(dict(cell=i,**report,passed=passed))
    save('controls.json',dict(classification='Counterexample candidate',partition=rows,inventory=inventories,
        all_passed=all(r['passed'] for r in rows+inventories)))
    assert all(r['passed'] for r in rows+inventories)
    print('SPECTRAL EOS CONTROLS',max(r['maximum_absolute_replay_error'] for r in rows),len(inventories),flush=True)


def increments():
    assert not THERMO.exists();THERMO.mkdir();plan=json.loads((thermo.original.OUT/'plan.json').read_text())
    _,changed,replacements=thermo.overlay();(THERMO/'candidate-run.py').write_text(changed)
    plan.update(checkpoint='ddfe050',substitutions=replacements,
        intervention='The existing material increment experiment with the new sourced H2+ spectral EOS, retained molecules. All 55 paths, analytic controls, endpoint and quadrature gates are unchanged.',
        physical_scope=json.loads((OUT/'plan.json').read_text())['boundary'])
    plan['bindings'].update({p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [
        g.ROOT/'verification/eos_h2plus_spectral.py',OUT/'plan.json',OUT/'controls.json',THERMO/'candidate-run.py']})
    (THERMO/'plan.json').write_text(json.dumps(plan,ensure_ascii=False,indent=2)+'\n')
    def save_thermo(name,value):(THERMO/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
    def verify_thermo():FunctionType(thermo.original.verify.__code__,dict(thermo.original.verify.__globals__,OUT=THERMO))()
    namespace=dict(thermo.original.run.__globals__,OUT=THERMO,save=save_thermo,verify=verify_thermo,Continuation=EOS)
    exec(compile(changed,str(THERMO/'candidate-run.py'),'exec'),namespace);namespace['run']()


def verify():
    plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    for path,digest in plan['runtime_original'].items():assert g.c.sha(path)==digest,path
    trees=json.loads((OUT/'source-tree.json').read_text())
    for name,folder in [('before',switch.CACHE/'source'),('after',CACHE/'source')]:
        for rel,digest in trees[name].items():assert g.c.sha(folder/rel)==digest,rel
    FunctionType(thermo.original.verify.__code__,dict(thermo.original.verify.__globals__,OUT=THERMO))()
    assert json.loads((OUT/'controls.json').read_text())['all_passed']
    if (OUT/'manifest.json').exists():
        for rel,digest in json.loads((OUT/'manifest.json').read_text())['sha256'].items():assert g.c.sha(g.ROOT/rel)==digest,rel
    else:
        save('manifest.json',dict(classification='Counterexample candidate',sha256={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in OUT.rglob('*') if p.is_file()},
            runtime={str(p):g.c.sha(p) for p in [LIB,CACHE/'excitation.so']}))
    for path,digest in json.loads((OUT/'manifest.json').read_text())['runtime'].items():assert g.c.sha(path)==digest,path
    print('PASS sourced H2+ native candidate bindings; read material-increments/result.json for gates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
