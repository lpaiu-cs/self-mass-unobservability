"""Use the two sourced molecular spectra and an explicit closed chemical-energy cycle."""
from types import FunctionType
import ctypes, inspect, json, shutil, sys
import numpy as np
import sympy as sp
import h2_spectre_refinement as h2
import eos_h2plus_spectral as previous
import eos_spectral_build_audit as build_audit

g=previous.g;OUT=g.OUT/'eos-molecular-spectral';CACHE=g.CACHE/'molecular-spectral'
NAME='free_eos_direct24_molecular_spectral';LIB=CACHE/'build/src'/('lib'+NAME+'.so.1.0.0')
THERMO=OUT/'material-increments'


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not OUT.exists() and not CACHE.exists();OUT.mkdir();CACHE.mkdir()
    h2.verify();previous.verify();build_audit.verify()
    source=CACHE/'source';original=previous.CACHE/'source';shutil.copytree(original,source)
    rows=json.loads((h2.OUT/'levels.json').read_text())['rows'];assert len(rows)==302
    energy=np.array([float(r['excitation_cm_inverse']) for r in rows])
    kelvin=energy*(6.62607015e-34*299792458*100/1.380649e-23)
    weights=np.array([(2*r['J']+1)*(1 if r['J']%2==0 else 3)/4 for r in rows])
    def array(name,values):
        return '    real(fp_kind), parameter :: '+name+'(302)=[ &\n'+', &\n'.join('      '+format(v,'.17e')+'_fp_kind' for v in values)+' ]\n'
    declaration=array('h2_spectral_kelvin',kelvin)+array('h2_spectral_weight',weights)+'''    real(fp_kind) :: h2sx(302), h2sw(302), h2sumw, h2meanx
'''
    tail='''    ! Declared H2 X-state sum. Physical plasma/continuum errors are not enclosed.
    if(ifh2.eq.3) then
       h2sx=h2_spectral_kelvin/exp(tl)
       h2sw=h2_spectral_weight*exp(-h2sx)
       h2sumw=sum(h2sw)
       h2meanx=sum(h2sw*h2sx)/h2sumw
       qh2=log(h2sumw)
       qh2t=h2meanx
       qh2tt=sum(h2sw*(h2sx-h2meanx)**2)/h2sumw-h2meanx
    endif
'''
    h2_D=next(r['dissociation_cm_inverse'] for r in rows if (r['v'],r['J'])==(0,0))
    babb=previous.data.previous
    h2p_D=next(r['binding_cm_inverse'] for r in json.loads((babb.OUT/'levels.json').read_text())['rows'] if (r['v'],r['N'])==(0,0))
    replacements={
        'src/mod_molecular_hydrogen.f90':[
            ('    real(fp_kind), intent(out) :: qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt',
             '    real(fp_kind), intent(out) :: qh2, qh2t, qh2tt, qh2plus, qh2plust, qh2plustt\n'+declaration),
            ('  end subroutine molecular_hydrogen',tail+'  end subroutine molecular_hydrogen')],
        'src/mod_ionization_data.f90':[
            ('  ! H2 I.P. from Huber and Herzberg, footnote b\n  real(fp_kind), parameter :: h2_ip = 124417.2_fp_kind',
             '  ! Declared mixed-source cycle: H2SPECTRE D(H2), Babb D(H2+), original I(H).\n  real(fp_kind), parameter :: h2_ip = '+h2_D+'_fp_kind + monatomic_ip(1) - '+h2p_D+'_fp_kind'),
            ('  ! H2 dissociation energy cm^{-1} taken from Huber and Herzberg\n  ! footnote a.\n  real(fp_kind), parameter :: h2diss = 36118.3_fp_kind',
             '  ! H2SPECTRE 7.4 computed D00; finite grid controls, no hard physical error bound.\n  real(fp_kind), parameter :: h2diss = '+h2_D+'_fp_kind')],
        'src/CMakeLists.txt':[('OUTPUT_NAME '+previous.NAME,'OUTPUT_NAME '+NAME)]}
    changes={}
    for rel,pairs in replacements.items():
        path=source/rel;before=path.read_text();after=before
        for old,new in pairs:assert after.count(old)==1;after=after.replace(old,new)
        reverse=after
        for old,new in reversed(pairs):assert reverse.count(new)==1;reverse=reverse.replace(new,old)
        assert reverse==before;path.write_text(after)
        shutil.copy2(original/rel,OUT/('before-'+path.name));shutil.copy2(path,OUT/path.name)
        changes[rel]=dict(before=g.c.sha(original/rel),after=g.c.sha(path),substitutions=pairs)
    shutil.copy2(previous.OUT/'direct_ion_bridge.f90',OUT/'direct_ion_bridge.f90')
    np.savez_compressed(OUT/'compiled-level-inputs.npz',kelvin=kelvin,weights=weights)
    save('source-tree.json',dict(before={p.relative_to(original).as_posix():g.c.sha(p) for p in original.rglob('*') if p.is_file()},
        after={p.relative_to(source).as_posix():g.c.sha(p) for p in source.rglob('*') if p.is_file()}))
    D,Dp,I=sp.symbols('D Dp I');Ip=D+I-Dp;Ipp=D+2*I-Ip
    assert sp.simplify(D+I-Ip-Dp)==0 and sp.simplify(Ipp-Dp-I)==0
    save('chemical-cycle.json',dict(classification='Proven',H2_D_cm_inverse=h2_D,H2plus_D_cm_inverse=h2p_D,
        identities=['I(H2)=D(H2)+I(H)-D(H2+)','I(H2+)=D(H2+)+I(H)'],
        scope='Algebraic energy-cycle consistency for the declared anchors. The H2+ excitation spectrum is MOL-D and its D00 anchor is Babb, a deliberately mixed-source model rather than one Hamiltonian or a certified physical error enclosure. Original atomic ionization energies are retained.'))
    wrapper=inspect.getsource(previous.increments).replace("checkpoint='ddfe050'","checkpoint='432ffb0'").replace(
        'new sourced H2+ spectral EOS','new H2/H2+ spectral EOS and declared chemical-energy anchors').replace(
        "g.ROOT/'verification/eos_h2plus_spectral.py'","g.ROOT/'verification/eos_molecular_spectral.py'")
    (OUT/'increments-wrapper.py').write_text(wrapper)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='432ffb0',changes=changes,
        partition_temperatures_K=json.loads((previous.OUT/'plan.json').read_text())['partition_temperatures_K'],
        native_partition_absolute_tolerance=1e-10,controls=[0,1175,1176,2972,5734],
        inventory_tolerance=1e-10,charge_tolerance=1e-10,
        intervention='Replace H2 option 3 in the same shared provider with the 302-level sum. Keep the H2+ 423-level sum bitwise. Align H2 dissociation to computed H2SPECTRE D00 and the H2+ ground anchor to the explicit Babb D00; derive both molecular ionization energies by the original chemical cycle. Retain molecules once on initialization. Reuse the identical original 55 material increment paths and gates.',
        boundary='The H2+ relative spectrum and its dissociation anchor are from different calculations. This declared mixed-source model has a consistent energy cycle, not a common-Hamiltonian certification. Excited electronic, quasibound, continuum and plasma occupation/physical errors remain. No running GR state/path or native library is replaced.',
        runtime_original=json.loads((previous.OUT/'manifest.json').read_text())['runtime'],
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in [g.ROOT/'verification/eos_molecular_spectral.py',
            h2.OUT/'manifest.json',previous.OUT/'manifest.json',build_audit.OUT/'manifest.json',babb.OUT/'manifest.json',
            OUT/'source-tree.json',OUT/'compiled-level-inputs.npz',OUT/'chemical-cycle.json',OUT/'increments-wrapper.py',
            OUT/'direct_ion_bridge.f90',g.OUT/'initial-state-17-4.npz']}))


def build():build_audit.build(sys.modules[__name__])


class EOS(previous.EOS):
    def __init__(self):
        fn=previous.EOS.__init__
        FunctionType(fn.__code__,dict(fn.__globals__,CACHE=CACHE))(self)


def controls():
    plan=json.loads((OUT/'plan.json').read_text());eos=EOS();parent=previous.EOS();rows=[]
    a=dict(np.load(OUT/'compiled-level-inputs.npz'))
    def evaluate(eos,T):
        fn=eos.inventory_lib.__mod_molecular_hydrogen_MOD_molecular_hydrogen
        fn.restype=None;fn.argtypes=[ctypes.POINTER(ctypes.c_int)]*3+[ctypes.POINTER(ctypes.c_double)]*7
        ints=[ctypes.c_int(i) for i in [0,3,2]];t=ctypes.c_double(np.log(float(T)));values=[ctypes.c_double() for _ in range(6)]
        fn(*[ctypes.byref(i) for i in ints],ctypes.byref(t),*[ctypes.byref(v) for v in values])
        return np.array([v.value for v in values]).reshape(2,3),t.value
    for T in plan['partition_temperatures_K']:
        value,t=evaluate(eos,T);old,_=evaluate(parent,T)
        x=a['kelvin']/np.exp(t);w=a['weights']*np.exp(-x);p=w/w.sum();mean=p@x
        expected=np.array([np.log(w.sum()),mean,p@((x-mean)**2)-mean])
        error=float(abs(value[0]-expected).max());same=np.array_equal(value[1],old[1])
        rows.append(dict(T_K=T,H2plus_bitwise_parent=same,H2_jets=value[0].tolist(),maximum_absolute_replay_error=error,
            passed=same and error<plan['native_partition_absolute_tolerance']))
    state=dict(np.load(g.OUT/'initial-state-17-4.npz'));inventories=[]
    for i in plan['controls']:
        r,t,X=state['lnd'][i],state['lnT'][i],state['X'][i];snap=eos.snapshot(r,t,X)
        report=previous.switch.molecular.p.c.s.check(snap,X,r,eos)
        passed=not report['missing_nonzero_elements'] and report['inventory_error']<plan['inventory_tolerance'] and report['charge_error']<plan['charge_tolerance']
        inventories.append(dict(cell=i,**report,passed=passed))
    save('controls.json',dict(classification='Counterexample candidate',partition=rows,inventory=inventories,
        all_passed=all(r['passed'] for r in rows+inventories)))
    assert all(r['passed'] for r in rows+inventories)
    print('MOLECULAR SPECTRAL CONTROLS',max(r['maximum_absolute_replay_error'] for r in rows),flush=True)


def increments():
    namespace=dict(previous.increments.__globals__,OUT=OUT,THERMO=THERMO,EOS=EOS)
    exec(compile((OUT/'increments-wrapper.py').read_text(),str(OUT/'increments-wrapper.py'),'exec'),namespace)
    namespace['increments']()


def verify():
    fn=previous.verify
    FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT,CACHE=CACHE,LIB=LIB,THERMO=THERMO,
        save=save,switch=previous))()
    build_audit.verify()
    assert json.loads((THERMO/'result.json').read_text())['all_passed']
    print('PASS both native spectra, declared chemical cycle and 55 material increment gates',flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
