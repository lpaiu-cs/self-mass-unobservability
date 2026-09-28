"""Direct Li/Be/B/F species in a separate FreeEOS model, not a physics certificate.

Counterexample candidate. The original library/runtime and ongoing GR paths are
unchanged. New elements inherit the declared FreeEOS metal approximations.
"""
from pathlib import Path
import ctypes, json, re, shutil, subprocess, sys
import numpy as np
import structured_enthalpy as e

OUT=e.OUT/'direct-ions';CACHE=Path('/home/lpaiu/work/direct-ions32')
CHARGES=np.r_[e.c.EZ,[3,4,5,9]];SYMBOLS=['Li','Be','B','F']


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def extend_array(text,name,addition,before_molecules=False):
    pattern=r'(\b'+name+r'\([^)]*\)\s*=\s*\[&)(.*?)(\])'
    matches=list(re.finditer(pattern,text,re.S));assert len(matches)==1,(name,len(matches))
    match=matches[0];body=match[2]
    if before_molecules:
        boundary=re.search(r'\n\s*!\s*H2(?:, H2\+)?\n',body);assert boundary,name
        body=body[:boundary.start()]+'\n'+addition+',&'+body[boundary.start():]
    else: body=body.rstrip()+',&\n'+addition
    return text[:match.start(2)]+body+text[match.end(2):]


def prepare():
    assert not (OUT/'plan.json').exists() and not (OUT/'source-patch.json').exists()
    OUT.mkdir(exist_ok=True);CACHE.mkdir(exist_ok=True)
    # Upstream website symlinks include unavailable old figures. They are not
    # scientific/library inputs; retain the original failure log separately.
    source=CACHE/'source';shutil.copytree(e.c.gr.SOURCE,source,dirs_exist_ok=True,ignore=shutil.ignore_patterns('www'))
    data=json.loads((e.s.v.OLD/'nist-trace-data.json').read_text());rows=data['rows']
    assert [r['Z'] for r in rows]==[3,4,5,9]
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='3d41c64',
        objective='Remove the positive charge-remapping obstruction by representing Li, Be, B and F directly in the same declared free-energy model.',
        elements=CHARGES.tolist(),ionization_stages=316,
        model='Original EOS1 options (3,1,-2), unchanged supported-element inputs. Append four elements and all their ionization stages; use archived NIST ground-state energies, the existing isoelectronic constant-term weight approximation and the existing generic metal pressure-ionization prescription.',
        limits='This is a declared extension of an approximate physical model, not an experimental EOS calibration or an error bound for plasma physics. Isotope shifts, missing nuclear/atomic partition information and model error remain separate.',
        checks='Zero-new-element agreement on all actual initial states; direct pure-element low-density electron-count limits; actual full-mixture thermodynamic derivative and stability checks.',
        supported_comparison_relative=1e-9,thermodynamic_derivative_relative=1e-5,
        archived_NIST_data_sha256=e.c.sha(e.s.v.OLD/'nist-trace-data.json'),
        original_library_sha256=e.c.sha(e.c.gr.CACHE/'build/src/libfree_eos.so.1.0.0'),
        original_runtime_modified=False,physical_EOS_certified=False))
    old={p.relative_to(source).as_posix():p.read_bytes() for p in (source/'src').glob('*.f90')}
    changed={}
    for rel,raw in old.items():
        text=raw.decode()
        text=re.sub(r'(\b(?:nelements\w*|neps\w*)\s*=\s*)20\b',r'\g<1>24',text)
        text=re.sub(r'(\bnions\w*\s*=\s*)295\b',r'\g<1>316',text)
        text=re.sub(r'(\bnions\w*\s*=\s*)297\b',r'\g<1>318',text)
        if text.encode()!=raw: changed[rel]=text
    def get(name): return changed.get('src/'+name,old['src/'+name].decode())
    def put(name,text): changed['src/'+name]=text
    name='free_eos_detailed.f90';text=get(name)
    text=extend_array(text,'iatomic_number','       3,4,5,9')
    text=extend_array(text,'iftracemetal','       0,0,0,0');put(name,text)
    name='mod_isotopic_mass_data.f90';text=get(name)
    text=extend_array(text,'isotopic_mass','       '+',&\n       '.join(f'isotopic_mass_2d(1,{z})' for z in [3,4,5,9]));put(name,text)
    name='mod_ionization_data.f90';text=get(name)
    addition='\n'.join('       ! '+sym+'\n       '+','.join(map(str,range(1,row['Z']+1)))+(',&' if k<3 else '')
        for k,(sym,row) in enumerate(zip(SYMBOLS,rows)))
    text=extend_array(text,'nion',addition,True)
    conversion=e.c.EV/(6.62607015e-27*(e.c.gr.C*100))
    addition='\n'.join('       ! '+sym+'\n       '+',&\n       '.join(f'{val*conversion:.17e}_fp_kind' for val in row['energies_eV'])+(',&' if k<3 else '')
        for k,(sym,row) in enumerate(zip(SYMBOLS,rows)))
    text=extend_array(text,'monatomic_ip',addition);put(name,text)
    name='mod_statistical_weight_data.f90';text=get(name)
    host=re.search(r'! Ne\n(.*?)! Na\n',text.split('iqion(nions_stat)',1)[1],re.S)[1]
    electron_weights=list(reversed([int(n) for n in re.findall(r'\d+',host)]))+[1]
    assert electron_weights==[1,2,1,2,1,6,9,4,9,6,1]
    text=extend_array(text,'iqneutral','       '+','.join(str(electron_weights[z]) for z in [3,4,5,9]))
    addition='\n'.join('       ! '+sym+'\n       '+','.join(str(electron_weights[z-q]) for q in range(1,z+1))+(',&' if k<3 else '')
        for k,(sym,z) in enumerate(zip(SYMBOLS,[3,4,5,9])))
    text=extend_array(text,'iqion',addition);put(name,text)
    name='mod_pi_fit.f90';text=get(name)
    text=text.replace('nelements_pi_fitp2 = 22','nelements_pi_fitp2 = 26')
    text=text.replace('spread(0._fp_kind,1,18)','spread(0._fp_kind,1,22)')
    addition='\n'.join('       ! '+sym+'\n       0._fp_kind,0._fp_kind,spread(6._fp_kind,1,'+str(z-2)+')'+(',&' if k<3 else '')
        for k,(sym,z) in enumerate(zip(SYMBOLS,[3,4,5,9])))
    text=extend_array(text,'pi_fit_ion_ln',addition,True);put(name,text)
    # Different library name/SONAME prevents the loader from substituting the
    # installed 20-element library. No installation or toolchain change.
    rel='src/CMakeLists.txt';original=(source/rel).read_bytes();old[rel]=original
    changed[rel]=original.decode()+'\nset_target_properties(${WRITEABLE_TARGET}free_eos PROPERTIES OUTPUT_NAME free_eos_direct24)\n'
    for rel,text in changed.items():
        (source/rel).write_text(text);dest=OUT/'sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text)
    weights=json.loads((e.s.v.OLD/'freeeos-weights.json').read_text())['atomic_weights']
    original=old['src/mod_isotopic_mass_data.f90'].decode()
    for z in [3,4,5,9]:
        block=original.split(f'isotopic_mass_{z:03d}(min_nmz:max_nmz) = [&',1)[1].split(']',1)[0]
        values=[float(v) for v in re.findall(r'([0-9]+\.[0-9]*)_fp_kind',block)];assert len(values)==73
        weights.append(values[9])
    save('model-data.json',dict(classification='Counterexample candidate',atomic_weights=weights,
        electron_count_term_weights=electron_weights,NIST_data=data,
        model_parameters='The added neutral metal radius factors are 1. Ion factors are 1 for the first two stages and exp(2) thereafter, exactly the same generic prescription as supported metals. This is an explicit model assumption, not an independent calibration.',
        isotope_entropy='This first direct-element control retains the existing element-group entropy convention. Isotope-resolved entropy terms need a separately declared consistent correction.'))
    save('source-patch.json',dict(classification='Proven',files={rel:dict(before=e.hashlib.sha256(old[rel]).hexdigest(),
        after=e.c.sha(source/rel)) for rel in changed},original_source_unchanged=True))
    bridge=(e.ROOT/'verification/common_eos_bridge.f90').read_text().replace('common_eos_aux','direct_ion_eos').replace('eps(20)','eps(24)')
    (OUT/'direct_ion_bridge.f90').write_text(bridge)
    # A control mixture restricted to the old supported elements, evaluated
    # before the new library is loaded. It is not a physical replacement star.
    base=dict(np.load(e.OUT/'initial-state.npz'));X=base['X'].copy()
    X[:,~np.isin(e.c.Z,e.c.EZ)]=0;X/=X.sum(1)[:,None];eos=e.s.v.ColdEOS()
    values=np.array([eos(2,r,t,x) for r,t,x in zip(base['lnd'],base['lnT'],X)])
    np.savez_compressed(OUT/'supported-reference.npz',lnd=base['lnd'],lnT=base['lnT'],X=X,values=values)
    print('PREPARED DIRECT IONS',len(changed),'source files',len(X),'reference states',flush=True)


def repair_dimensions():
    """Preserve the failed candidate and fix the missed neps_local dimension."""
    assert not (OUT/'dimension-fix.json').exists()
    for name in ['build.json','source-patch.json']:
        shutil.copy2(OUT/name,OUT/('before-dimension-fix-'+name))
    shutil.copy2(e.OUT/'direct-ions-supported.log',OUT/'dimension-failure.log')
    patch=json.loads((OUT/'source-patch.json').read_text());changes={}
    for path in (CACHE/'source/src').glob('*.f90'):
        before=OUT/('before-dimension-fix-'+path.name)
        raw=before.read_text() if before.exists() else path.read_text()
        text=re.sub(r'(\bneps\w*\s*=\s*)20\b',r'\g<1>24',raw)
        if text==raw: continue
        rel='src/'+path.name
        if not before.exists(): shutil.copy2(path,before)
        changes[rel]=dict(before=e.c.sha(before))
        path.write_text(text);(OUT/'sources'/rel).write_text(text)
        changes[rel]['after']=e.c.sha(path);patch['files'][rel]['after']=e.c.sha(path)
    assert changes.keys()=={'src/free_eos_detailed.f90','src/mod_free_eos.f90'}
    save('source-patch.json',patch)
    save('dimension-fix.json',dict(classification='Proven',changes=changes,
        cause='The first parameter replacement matched neps but missed neps_local, so the input-size guard rejected all 24-element queries. The revised preparation also covers suffixed neps parameters.'))


def build():
    build_at(CACHE/'source',CACHE/'build','free_eos_direct24',CACHE/'direct_ion_bridge.so')


def build_at(source,folder,name,bridge,prefix=''):
    commands=[['cmake','-S',str(source),'-B',str(folder),'-DCMAKE_BUILD_TYPE=Release',
        '-DBUILD_SHARED_LIBS=ON','-DBUILD_TEST=ON','-DBUILD_DOC=OFF','-DBUILD_DOX_DOC=OFF'],
        ['cmake','--build',str(folder),'--target','free_eos','-j','2']]
    records=[]
    for i,cmd in enumerate(commands):
        result=subprocess.run(cmd,capture_output=True,text=True);(OUT/f'{prefix}build-{i}.log').write_text(result.stdout+result.stderr)
        records.append(dict(command=cmd,returncode=result.returncode));save(prefix+'build.json',dict(classification='Counterexample candidate',commands=records,completed=False))
        assert result.returncode==0,(i,result.stderr[-2000:])
    module=next(folder.rglob('mod_free_eos.mod')).parent
    cmd=['gfortran','-O2','-fPIC','-shared','-I'+str(module),str(OUT/'direct_ion_bridge.f90'),
        '-L'+str(folder/'src'),'-Wl,-rpath,'+str(folder/'src'),'-l'+name,'-o',str(bridge)]
    result=subprocess.run(cmd,capture_output=True,text=True);(OUT/(prefix+'bridge-build.log')).write_text(result.stdout+result.stderr)
    records.append(dict(command=cmd,returncode=result.returncode));assert result.returncode==0,result.stderr
    save(prefix+'build.json',dict(classification='Counterexample candidate',commands=records,completed=True,
        sha256={str(p):e.c.sha(p) for p in [bridge,folder/'src'/('lib'+name+'.so.1.0.0')]},
        original_runtime_modified=False,physical_EOS_certified=False))
    print('BUILT SEPARATE LIBRARY',name,flush=True)


def build_integral():
    source=CACHE/'integral-source';assert not source.exists()
    save('integral-plan.json',dict(classification='Counterexample candidate',checkpoint='6dd7317',
        change='Replace only the default electron Fermi evaluation morder=13 (CT plus fitted relativistic correction) by existing morder=21 direct numerical integrals. Preserve EOS1 radiation, exchange, Coulomb, ionization and all compositions.',
        controls='Original finite mixture derivative tolerance 1e-5 and steps 1e-4,5e-5. First check the observed boundary and representatives, then all initial cells.',
        scope='A changed electron numerical evaluator, not a new physical plasma validation or rigorous quadrature error certificate.'))
    shutil.copytree(CACHE/'source',source);changes={}
    for rel in ['src/mod_free_eos.f90','src/CMakeLists.txt']:
        path=source/rel;changes[rel]=dict(before=e.c.sha(path));text=path.read_text()
        if rel.endswith('.f90'):
            text,count=re.subn(r'^    morder = 13$',r'    morder = 21',text,flags=re.M);assert count==1
        else: text=text.replace('OUTPUT_NAME free_eos_direct24)','OUTPUT_NAME free_eos_direct24_integral)')
        path.write_text(text);dest=OUT/'integral-sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text)
        changes[rel]['after']=e.c.sha(path)
    save('integral-source-patch.json',dict(classification='Proven',files=changes))
    build_at(source,CACHE/'integral-build','free_eos_direct24_integral',CACHE/'direct_ion_integral_bridge.so','integral-')


def build_tight_integral():
    source=CACHE/'tight-integral-source';assert not source.exists()
    save('tight-integral-plan.json',dict(classification='Counterexample candidate',
        change='Tighten only direct Fermi quadrature fderr from 1e-9 to 1e-12. Keep the original finite derivative gates and both earlier candidates.',
        reason='The initial direct-integral representative derivative score is 8.60e-6 against a 1e-5 gate. The upstream local quadrature error estimator is not a rigorous global bound, and its noise can be amplified by differencing.',
        tolerance=1e-12,physical_or_rigorous_quadrature_certificate=False))
    shutil.copytree(CACHE/'integral-source',source);changes={}
    for rel in ['src/fermi_dirac_direct.f90','src/CMakeLists.txt']:
        path=source/rel;changes[rel]=dict(before=e.c.sha(path));text=path.read_text()
        if rel.endswith('.f90'):
            assert text.count('fderr = 1.e-09_fp_kind')==1;text=text.replace('fderr = 1.e-09_fp_kind','fderr = 1.e-12_fp_kind')
        else: text=text.replace('OUTPUT_NAME free_eos_direct24_integral)','OUTPUT_NAME free_eos_direct24_integral_tight)')
        path.write_text(text);dest=OUT/'tight-integral-sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text)
        changes[rel]['after']=e.c.sha(path)
    save('tight-integral-source-patch.json',dict(classification='Proven',files=changes))
    build_at(source,CACHE/'tight-integral-build','free_eos_direct24_integral_tight',CACHE/'direct_ion_integral_tight_bridge.so','tight-integral-')


def build_full_integral():
    source=CACHE/'full-integral-source';assert not source.exists()
    save('full-integral-plan.json',dict(classification='Counterexample candidate',checkpoint='792a8e2',
        change='Keep morder=21 and fderr=1e-12, and replace the common fermi_dirac_ct evaluator for the exchange K/I call as well. Use F_3/2 prime = 3/2 F_1/2 and F_3/2 double prime = 3/4 F_-1/2, with the existing direct integrator for higher eta derivatives.',
        root_cause='The EOS1 exchange_gcpf call independently uses CT derivatives in its K integral. Changing the main Fermi morder alone leaves this branch in the actual flow.',
        controls='Independent 60-digit polylog values for orders 0 through 5, the observed branch scan, original representative and whole-mixture derivative controls.',
        supported_scope='The registered EOS1 option path. Alternate FreeEOS options and the unused CT entropy-difference approximation are not claimed validated.',
        physical_or_continuous_certificate=False))
    shutil.copytree(CACHE/'tight-integral-source',source);changes={}
    replacement='''function fermi_dirac_ct(eta, nderiv, ifsimple)
  use mod_free_eos_constants, only: pi
  real(fp_kind) fermi_dirac_ct
  real(fp_kind), intent(in) :: eta
  integer, intent(in) :: nderiv, ifsimple
  real(fp_kind), save :: last_eta = -huge(1._fp_kind), cached(0:5) = 0._fp_kind
  logical, save :: have_zero = .false., have_middle = .false., have_five = .false.
  real(fp_kind) answer(7)
  if(nderiv.lt.0.or.nderiv.gt.5) error stop 'direct CT replacement: bad derivative order'
  if(ifsimple.eq.1) then
     fermi_dirac_ct = 0.75_fp_kind*sqrt(pi)*exp(eta)
     return
  endif
  ! The existing EOS library is serial/stateful; this cache has the same scope.
  if(eta.ne.last_eta) then
     last_eta = eta
     have_zero = .false.
     have_middle = .false.
     have_five = .false.
  endif
  if(nderiv.eq.0.and..not.have_zero) then
     call fermi_dirac_direct(1.5_fp_kind,eta,0._fp_kind,answer(:1))
     cached(0) = answer(1)
     have_zero = .true.
  elseif(1.le.nderiv.and.nderiv.le.4.and..not.have_middle) then
     call fermi_dirac_direct(0.5_fp_kind,eta,0._fp_kind,answer)
     cached(1:4) = 1.5_fp_kind*[answer(1),answer(2),answer(4),answer(7)]
     have_middle = .true.
  elseif(nderiv.eq.5.and..not.have_five) then
     call fermi_dirac_direct(-0.5_fp_kind,eta,0._fp_kind,answer)
     cached(5) = 0.75_fp_kind*answer(7)
     have_five = .true.
  endif
  fermi_dirac_ct = cached(nderiv)
end function fermi_dirac_ct'''
    for rel in ['src/fermi_dirac_ct.f90','src/CMakeLists.txt']:
        path=source/rel;changes[rel]=dict(before=e.c.sha(path));text=path.read_text()
        if rel.endswith('.f90'):
            text,count=re.subn(r'^function fermi_dirac_ct\(eta, nderiv, ifsimple\).*?^end function fermi_dirac_ct$',replacement,text,flags=re.M|re.S);assert count==1
        else: text=text.replace('OUTPUT_NAME free_eos_direct24_integral_tight)','OUTPUT_NAME free_eos_direct24_integral_full)')
        path.write_text(text);dest=OUT/'full-integral-sources'/rel;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(text)
        changes[rel]['after']=e.c.sha(path)
    save('full-integral-source-patch.json',dict(classification='Proven',files=changes))
    build_at(source,CACHE/'full-integral-build','free_eos_direct24_integral_full',CACHE/'direct_ion_integral_full_bridge.so','full-integral-')


class EOS:
    def __init__(self,evaluator='ct'):
        names={'ct':'direct_ion_bridge.so','integral':'direct_ion_integral_bridge.so','tight':'direct_ion_integral_tight_bridge.so','full':'direct_ion_integral_full_bridge.so'}
        self.lib=ctypes.CDLL(str(CACHE/names[evaluator]));self.call=self.lib.direct_ion_eos
        self.call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,
            np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS'),ctypes.POINTER(ctypes.c_int)]
        self.call.restype=None;self.weights=np.array(json.loads((OUT/'model-data.json').read_text())['atomic_weights'])
        self.mapping=np.array([[float(z==q) for q in CHARGES] for z in e.c.Z])
        assert np.all(self.mapping.sum(1)==1) and np.array_equal(self.mapping@CHARGES,e.c.Z)
        assert np.array_equal(self.mapping@(CHARGES**2),e.c.Z**2)

    def __call__(self,mode,value,t,x):
        assert x.shape==(26,) and np.all(x>=0) and abs(x.sum()-1)<1e-11
        ym=(x/e.c.A)@self.mapping;cx=float(ym@self.weights);eps=np.ascontiguousarray(ym/cx)
        seed=np.zeros(24);seed[2 if x[e.c.NAMES.index('c12')]<.5 else 0]=1
        seed/=seed@self.weights;out=np.full(22,np.nan);info=ctypes.c_int(-999)
        self.call(0,-20.,np.log(1e6),seed,out,ctypes.byref(info));assert info.value==0
        self.call(mode,float(value+np.log(cx) if mode==2 else value),float(t),eps,out,ctypes.byref(info))
        assert info.value==0 and np.all(np.isfinite(out)),(info.value,mode,value,t,x)
        out[0]/=cx;out[2:4]*=cx;out[9:11]*=cx
        return out


def supported_control():
    ref=dict(np.load(OUT/'supported-reference.npz'));eos=EOS()
    actual=np.array([eos(2,r,t,x)[:12] for r,t,x in zip(ref['lnd'],ref['lnT'],ref['X'])])
    scale=np.maximum(abs(ref['values']),1.);difference=abs(actual-ref['values'])/scale
    maximum=float(difference.max());passed=maximum<1e-9
    np.savez_compressed(OUT/'supported-control.npz',actual=actual,relative_difference=difference)
    save('supported-control.json',dict(classification='Counterexample candidate',passed=passed,
        cells=len(actual),maximum_relative_difference=maximum,worst_index=list(map(int,np.unravel_index(difference.argmax(),difference.shape))),
        physical_or_continuous_certificate=False))
    assert passed,maximum;print('SUPPORTED CONTROL',maximum,flush=True)


def control_plan():
    assert not (OUT/'control-plan.json').exists()
    save('control-plan.json',dict(classification='Counterexample candidate',
        density_g_cm3=1e-5,temperature_K=3e6,pure_free_electron_relative=1e-3,
        pure_species=[n for n,z in zip(e.c.NAMES,e.c.Z) if z in [3,4,5,9]],
        finite_log_steps=[1e-4,5e-5],derivative_relative=1e-5,
        scope='Every actual initial cell at fixed physical baryon density, temperature and 26-species composition. Finite derivative checks and a dilute near-stripped limit are numerical controls, not continuous or physical uncertainty bounds.'))


def old_full_reference():
    # A separate process never loads the new library for this reference.
    assert not (OUT/'old-full-reference.npz').exists()
    base=dict(np.load(e.OUT/'initial-state.npz'));eos=e.s.v.ColdEOS()
    values=np.array([eos(2,r,t,x) for r,t,x in zip(base['lnd'],base['lnT'],base['X'])])
    np.savez_compressed(OUT/'old-full-reference.npz',values=values,
        state_sha256=e.c.sha(e.OUT/'initial-state.npz'),library_sha256=e.c.sha(e.c.gr.CACHE/'build/src/libfree_eos.so.1.0.0'))
    print('OLD FULL REFERENCE',len(values),flush=True)


def pure_control():
    plan=json.loads((OUT/'control-plan.json').read_text());eos=EOS();rows=[];raw=[]
    for name in plan['pure_species']:
        k=e.c.NAMES.index(name);x=np.eye(26)[k]
        a=eos(2,np.log(plan['density_g_cm3']),np.log(plan['temperature_K']),x)
        # rmue is electron number divided by N_A per cm^3; baryon rho recovers Ye.
        actual=a[13]/a[0];target=e.c.Z[k]/e.c.A[k]
        error=float(abs(actual/target-1));raw.append(a)
        rows.append(dict(species=name,free_Ye=float(actual),fully_stripped_Ye=float(target),
            relative_difference=error,passed=error<plan['pure_free_electron_relative']))
    np.savez_compressed(OUT/'pure-control.npz',values=np.array(raw))
    passed=all(r['passed'] for r in rows)
    save('pure-control.json',dict(classification='Counterexample candidate',rows=rows,passed=passed,
        model_physics_certified=False,plan_sha256=e.c.sha(OUT/'control-plan.json')))
    print('PURE CONTROL',rows,flush=True);assert passed


def mixture_control(evaluator='ct'):
    plan=json.loads((OUT/'control-plan.json').read_text());base=dict(np.load(e.OUT/'initial-state.npz'))
    eos=EOS(evaluator);values=[];errors=[];prefix={'ct':'','integral':'integral-','tight':'tight-integral-','full':'full-integral-'}[evaluator]
    columns=['P_T','P_rho','u_T','entropy_T','Maxwell','entropy_rho']
    for i,(r,t,x) in enumerate(zip(base['lnd'],base['lnT'],base['X'])):
        a=eos(2,r,t,x);values.append(a);cell=[]
        for step in plan['finite_log_steps']:
            dT=(eos(2,r,t+step,x)-eos(2,r,t-step,x))/(2*step)
            dr=(eos(2,r+step,t,x)-eos(2,r-step,t,x))/(2*step)
            cell.append([dT[1]/(a[1]*a[6])-1,dr[1]/(a[1]*a[5])-1,
                dT[2]/a[10]-1,dT[3]*np.exp(t)/a[10]-1,
                (dr[2]-(a[1]-dT[1])/a[0])/max(abs(a[10]),1.),
                (dr[3]+dT[1]/(a[0]*np.exp(t)))/max(abs(a[3]),1.)])
        errors.append(cell)
        if i%1000==0: print('DIRECT MIXTURE',i,flush=True)
    values=np.array(values);errors=np.array(errors)
    old=dict(np.load(OUT/'old-full-reference.npz'))
    assert str(old['state_sha256'])==e.c.sha(e.OUT/'initial-state.npz')
    difference=values[:,:12]-old['values'];sound=values[:,21]/(e.c.gr.C*100)**2
    worst=np.max(abs(errors),axis=(0,1));finite=bool(np.all(np.isfinite(errors)))
    stable=bool(np.all(values[:,10]>0) and np.all(values[:,5]>0) and np.all(sound>0) and np.all(sound<1))
    passed=bool(finite and stable and max(worst)<plan['derivative_relative'])
    np.savez_compressed(OUT/(prefix+'mixture-control.npz'),values=values,errors=errors,old_model_difference=difference,
        columns=columns,state_sha256=e.c.sha(e.OUT/'initial-state.npz'))
    save(prefix+'mixture-control.json',dict(classification='Counterexample candidate',cells=len(values),passed=passed,
        finite=finite,positive_capacity_compressibility_subluminal=stable,
        worst_relative=dict(zip(columns,map(float,worst))),
        worst_index=list(map(int,np.unravel_index(abs(errors).argmax(),errors.shape))),
        sound2_over_c2_range=[float(sound.min()),float(sound.max())],
        old_model_pressure_relative_max=float(max(abs(difference[:,1]/values[:,1]))),
        old_model_energy_over_cvT_max=float(max(abs(difference[:,2]/values[:,10]))),
        plan_sha256=e.c.sha(OUT/'control-plan.json'),physical_or_continuous_certificate=False))
    print('MIXTURE CONTROL',passed,dict(zip(columns,worst)),flush=True);assert passed


def integral_mixture_control(): mixture_control('integral')


def tight_integral_mixture_control(): mixture_control('tight')


def full_integral_mixture_control(): mixture_control('full')


def failure_probe(old=False):
    """Retain the original failure and compare step sizes in each library."""
    saved=dict(np.load(OUT/'mixture-control.npz'));base=dict(np.load(e.OUT/'initial-state.npz'))
    selected=np.flatnonzero(np.max(abs(saved['errors']),axis=(1,2))>=1e-5)
    eos=e.s.v.ColdEOS() if old else EOS();rows=[];raw=[]
    steps=[1e-4,5e-5,1e-5,5e-6,1e-6,5e-7,1e-7]
    for i in selected:
        r,t,x=base['lnd'][i],base['lnT'][i],base['X'][i];a=eos(2,r,t,x);cell=[]
        for h in steps:
            ap=eos(2,r,t+h,x);am=eos(2,r,t-h,x);dT=(ap-am)/(2*h)
            cell.append([dT[1]/(a[1]*a[6])-1,dT[2]/a[10]-1,dT[3]*np.exp(t)/a[10]-1])
            raw.append([a,ap,am])
        rows.append(dict(cell=int(i),rho_g_cm3=float(np.exp(r)),T_K=float(np.exp(t)),errors=cell))
    label='old' if old else 'new'
    np.savez_compressed(OUT/f'mixture-failure-{label}.npz',values=np.array(raw),cells=selected,steps=steps)
    save(f'mixture-failure-{label}.json',dict(classification='Counterexample candidate',rows=rows,
        steps=steps,columns=['P_T','u_T','entropy_T'],original_failure_retained=True))
    print('FAILURE PROBE',label,'cells',len(rows),flush=True)


def old_failure_probe(): failure_probe(True)


def branch_probe(evaluator='ct'):
    base=dict(np.load(e.OUT/'initial-state.npz'));eos=EOS(evaluator);i=2972;prefix='full-integral-' if evaluator=='full' else ''
    offsets=np.linspace(-1.5e-4,1.5e-4,121);r,t,x=base['lnd'][i],base['lnT'][i],base['X'][i]
    values=np.array([eos(2,r,t+dt,x) for dt in offsets])
    derivatives=np.column_stack([values[:,1]*values[:,6],values[:,10],values[:,10]/np.exp(t+offsets)])
    residual=np.diff(values[:,[1,2,3]],axis=0)-np.diff(offsets)[:,None]*(derivatives[1:]+derivatives[:-1])/2
    score=abs(residual)/np.maximum(abs(values[:-1,[1,2,3]]),1.)
    k=int(np.unravel_index(score.argmax(),score.shape)[0])
    np.savez_compressed(OUT/(prefix+'branch-probe.npz'),offsets=offsets,values=values,residual=residual)
    save(prefix+'branch-probe.json',dict(classification='Counterexample candidate',cell=i,
        eta_range=[float(values[:,12].min()),float(values[:,12].max())],
        largest_local_residual_lnT_interval=offsets[k:k+2].tolist(),
        eta_at_interval=values[k:k+2,12].tolist(),
        maximum_integrated_derivative_residual_relative=score.max(axis=0).tolist(),
        original_mixture_failure_retained=True,
        scope='Finite scan and trapezoidal consistency diagnostic, not a rigorous discontinuity or derivative bound.'))
    print('BRANCH PROBE',json.loads((OUT/(prefix+'branch-probe.json')).read_text()),flush=True)


def ct_boundary_probe():
    """Read native exchange state and evaluate its unchanged CT piece boundary."""
    base=dict(np.load(e.OUT/'initial-state.npz'));scan=dict(np.load(OUT/'branch-probe.npz'));eos=EOS();i=2972
    primed=ctypes.c_double.in_dll(eos.lib,'__mod_master_exchange_data_MOD_psiprime')
    rows=[]
    for offset in scan['offsets']:
        a=eos(2,base['lnd'][i],base['lnT'][i]+offset,base['X'][i]);rows.append([offset,a[12],primed.value])
    rows=np.array(rows);cross=np.flatnonzero((rows[:-1,2]-1)*(rows[1:,2]-1)<0)
    call=eos.lib.__mod_fermi_dirac_MOD_fermi_dirac_ct
    call.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int)]
    call.restype=ctypes.c_double;zero=ctypes.c_int(0);limits=[]
    for boundary in [1.,4.]:
        vals=[]
        for point in [np.nextafter(boundary,-np.inf),np.nextafter(boundary,np.inf)]:
            eta=ctypes.c_double(point)
            vals.append([call(ctypes.byref(eta),ctypes.byref(ctypes.c_int(k)),ctypes.byref(zero)) for k in range(6)])
        limits.append(dict(eta=boundary,left=vals[0],right=vals[1],
            relative_difference=((np.array(vals[1])-vals[0])/np.array(vals[0])).tolist()))
    np.savez_compressed(OUT/'ct-boundary-probe.npz',exchange=rows)
    save('ct-boundary-probe.json',dict(classification='Counterexample candidate',cell=i,
        primed_eta_one_crossing_intervals=[rows[k:k+2].tolist() for k in cross],native_CT_limits=limits,
        scope='Native values adjacent to the piece boundary and the actual final exchange argument. This localizes the observed finite failure; arbitrary-precision coefficient bounds are separate.'))
    print('CT BOUNDARY',json.loads((OUT/'ct-boundary-probe.json').read_text()),flush=True)


def ct_interval():
    """A nonzero jump for the stated real CT formulas, not a stellar bound."""
    from mpmath import iv
    iv.dps=60
    path=CACHE/'source/src/mod_fermi_dirac.f90';text=path.read_text();tables=[]
    for name in ['p0','q0']:
        block=text.split(name+'(0:nx,ntab) = reshape([&',1)[1].split('], shape',1)[0]
        block=re.sub(r'!.*','',block)
        coefficients=re.findall(r'([+-]?(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?)_fp_kind',block)
        assert len(coefficients)==15
        # Enclose each actual binary64 coefficient and its decimal literal.
        # The width also covers conversion, rather than assuming a decimal
        # coefficient is identical to a compiled floating point constant.
        intervals=[]
        for value in coefficients:
            val=float(value)
            lo=np.nextafter(val,-np.inf);hi=np.nextafter(val,np.inf)
            intervals.append(iv.mpf([lo,hi]))
        tables.append([intervals[k:k+5] for k in [0,5,10]])
    def poly(coeff,x):
        result=iv.mpf(0)
        for value in reversed(coeff): result=result*x+value
        return result
    def ratio(j,x):
        p,q=tables[0][j],tables[1][j];P,Q=poly(p,x),poly(q,x)
        dp=poly([k*p[k] for k in range(1,5)],x);dq=poly([k*q[k] for k in range(1,5)],x)
        return P/Q,(dp*Q-P*dq)/Q**2
    y=iv.exp(1);A=.75*np.sqrt(np.pi);constant=iv.mpf([A-8*np.spacing(A),A+8*np.spacing(A)])
    mathematical_constant=iv.mpf('.75')*iv.sqrt(iv.pi)
    assert mathematical_constant.a>constant.a and mathematical_constant.b<constant.b
    R,D=ratio(0,y);left=[y*(constant+y*R),y*constant+2*y*y*R+y**3*D]
    R,D=ratio(1,iv.mpf(1));right=[R,D]
    jump=[b-a for a,b in zip(left,right)]
    assert jump[0].a>0 and jump[1].b<0
    save('ct-interval.json',dict(classification='Proven',passed=True,precision_decimal_digits=60,
        eta=1,function_jump=str(jump[0]),first_derivative_jump=str(jump[1]),
        source_sha256=e.c.sha(path),
        statement='The two declared real rational CT branches, including one-ulp coefficient enclosures and eight-ulp root-pi constant enclosure, have unequal limits at eta=1. Therefore this piecewise approximation is not C0 and its returned first derivative also has a nonzero jump.',
        limits='This is a coefficient-level mathematical statement. It does not enclose every floating point EOS operation, prove the exact global stellar jump, or bound physical EOS error.'))
    print('CT INTERVAL',jump,flush=True)


def integral_control(evaluator='integral'):
    base=dict(np.load(e.OUT/'initial-state.npz'));old=np.load(OUT/'mixture-control.npz')['values'];eos=EOS(evaluator)
    prefix={'integral':'integral-','tight':'tight-integral-','full':'full-integral-'}[evaluator]
    selected=np.unique(np.r_[np.linspace(0,5734,17).astype(int),np.argmax(base['X'],axis=0),2972])
    rows=[];raw=[]
    for i in selected:
        r,t,x=base['lnd'][i],base['lnT'][i],base['X'][i];a=eos(2,r,t,x);raw.append(a);errors=[]
        for h in [1e-4,5e-5]:
            dT=(eos(2,r,t+h,x)-eos(2,r,t-h,x))/(2*h)
            dr=(eos(2,r+h,t,x)-eos(2,r-h,t,x))/(2*h)
            errors.append([dT[1]/(a[1]*a[6])-1,dr[1]/(a[1]*a[5])-1,
                dT[2]/a[10]-1,dT[3]*np.exp(t)/a[10]-1,
                (dr[2]-(a[1]-dT[1])/a[0])/max(abs(a[10]),1.),
                (dr[3]+dT[1]/(a[0]*np.exp(t)))/max(abs(a[3]),1.)])
        rows.append(dict(cell=int(i),errors=errors,pressure_difference_from_CT=float(a[1]/old[i,1]-1),
            energy_difference_over_cvT=float((a[2]-old[i,2])/a[10])))
        print('INTEGRAL CONTROL',int(i),float(np.max(abs(np.array(errors)))),flush=True)
    worst=float(max(np.max(abs(np.array(row['errors']))) for row in rows));passed=worst<1e-5
    np.savez_compressed(OUT/(prefix+'control.npz'),values=np.array(raw),cells=selected)
    save(prefix+'control.json',dict(classification='Counterexample candidate',rows=rows,passed=passed,
        maximum_derivative_error=worst,physical_or_continuous_certificate=False))
    assert passed,worst


def tight_integral_control(): integral_control('tight')


def full_integral_control(): integral_control('full')


def full_branch_probe(): branch_probe('full')


def fermi_polylog_control():
    import mpmath as mp
    mp.mp.dps=60;eos=EOS('full');call=eos.lib.__mod_fermi_dirac_MOD_fermi_dirac_ct
    call.argtypes=[ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_int)]
    call.restype=ctypes.c_double;zero=ctypes.c_int(0);rows=[]
    for eta in [-20.,-2.,0.,1.-1e-9,1.,1.+1e-9,4.-1e-9,4.,4.+1e-9,10.,40.,100.]:
        value=ctypes.c_double(eta);actual=[];reference=[];errors=[]
        for k in range(6):
            a=call(ctypes.byref(value),ctypes.byref(ctypes.c_int(k)),ctypes.byref(zero))
            ref=-mp.gamma(mp.mpf('2.5'))*mp.polylog(mp.mpf('2.5')-k,-mp.exp(mp.mpf(eta)))
            assert abs(mp.im(ref))<mp.mpf('1e-45');ref=mp.re(ref)
            actual.append(a);reference.append(str(ref));errors.append(float(abs(mp.mpf(a)-ref)/max(abs(ref),mp.mpf('1e-8'))))
        rows.append(dict(eta=eta,actual=actual,reference=reference,relative_errors=errors))
        print('POLYLOG CONTROL',eta,max(errors),flush=True)
    worst=max(max(r['relative_errors']) for r in rows);passed=worst<1e-8
    save('fermi-polylog-control.json',dict(classification='Counterexample candidate',rows=rows,passed=passed,
        maximum_scaled_error=worst,scope='Independent high-precision numerical reference, not interval-certified quadrature or full plasma physics.'))
    assert passed,worst


def fermi_tail_control():
    """Uniform upper-tail bound for nonrelativistic Fermi eta derivatives."""
    import sympy as sp
    from mpmath import iv
    iv.dps=60;q=sp.symbols('q');polynomial=q;coefficients=[]
    for k in range(6):
        p=sp.Poly(polynomial,q)
        assert p.eval(0)==0 and p.degree()==k+1
        coefficients.append([int(p.nth(j)) for j in range(1,k+2)])
        polynomial=sp.expand(q*(1-q)*sp.diff(polynomial,q))
    assert coefficients[1]==[1,-1] and coefficients[2]==[1,-3,2]
    rows=[];width=iv.mpf('0.001')
    for eta in [-20,-2,0,1,4,10,40,100]:
        B=iv.mpf(eta)+width;L=40+abs(iv.mpf(eta))+width;delta=iv.exp(B-L)
        for nu in ['-0.5','0.5','1.5']:
            order=iv.mpf(nu);denominator=1-max(float(nu),0)/L;bounds=[]
            assert denominator.a>0 and delta.b<1
            for cs in coefficients:
                factor=sum(abs(c)*delta**j for j,c in enumerate(cs))
                bound=factor*delta*L**order/denominator
                bounds.append(str(bound))
            rows.append(dict(eta_center=eta,eta_half_width='0.001',nu=nu,L=str(L),tail_absolute_bounds=bounds))
    save('fermi-tail-control.json',dict(classification='Proven',passed=True,
        definition='F_nu(eta)=integral_0^infinity x^nu/(1+exp(x-eta)) dx, nu>-1.',
        smoothness='On each compact eta interval, every eta derivative is dominated by an integrable polynomial times exp(B-x) at large x and by x^nu near zero. Thus F_nu is C infinity.',
        derivative_polynomials=coefficients,
        bound='For eta<=B and L>max(nu,0), delta=exp(B-L)<=1, write d_eta^k q=sum_j a_j q^j. The absolute tail is at most [sum_j |a_j| delta^(j-1)] delta L^nu/[1-max(nu,0)/L].',
        proof='Use q<=exp(B-x)<=delta on the tail, bound the polynomial by q times its absolute coefficients at delta, and use log(x/L)<=x/L-1 for nu>=0 or x^nu<=L^nu for nu<0; integrate the resulting exponential majorant.',
        rows=rows,precision_decimal_digits=60,
        limits='Only the omitted upper tail of the declared nonrelativistic integrals and eta derivatives. The finite-interval quadrature, relativistic beta derivatives, other EOS components and stellar propagation remain outside this bound.'))
    print('FERMI TAIL',len(rows),'parameter boxes, orders 0 through 5',flush=True)


def isotope_no_go():
    import sympy as sp
    M,N,me,q=sp.symbols('M N me q',positive=True)
    difference=(M-(q+1)*me)*(N-q*me)-(N-(q+1)*me)*(M-q*me)
    assert sp.simplify(difference-me*(M-N))==0
    constants=CACHE/'source/src/mod_free_eos_constants.f90';source=constants.read_text()
    assert 'ifreducedmass = 1' in (CACHE/'source/src/mod_free_eos.f90').read_text()
    assert 'isotopic_mass(ielement) - real((nion(ion)),fp_kind)*electron_mass' in (CACHE/'source/src/ionize.f90').read_text()
    # Reproduce the declared source constants for a finite magnitude census.
    light=2.99792458e10;h=6.62607015e-27;charge=1.602176634e-19*.1*light
    electron=e.c.NA*109737.31568160*(h/(charge*charge))*(h/charge)**2*light/(2*np.pi*np.pi)
    weights=np.array(json.loads((OUT/'model-data.json').read_text())['atomic_weights']);rows=[]
    for i,name in enumerate(e.c.NAMES):
        j=int(np.flatnonzero(CHARGES==e.c.Z[i])[0]);mass=e.c.W[i];reference=weights[j]
        if mass==reference: continue
        stages=np.arange(int(e.c.Z[i]));assert min(mass,reference)>e.c.Z[i]*electron
        ratio=((mass-(stages+1)*electron)/(mass-stages*electron))/((reference-(stages+1)*electron)/(reference-stages*electron))
        errors=ratio**1.5-1
        rows.append(dict(species=name,isotope_atomic_mass=float(mass),reference_atomic_mass=float(reference),
            maximum_adjacent_stage_mass_factor_relative=float(max(abs(errors)))))
    save('isotope-no-go.json',dict(classification='Proven',passed=True,
        identity='(M-(q+1)me)(N-q me)-(N-(q+1)me)(M-q me)=me(M-N)',
        assumptions='Classical translational partition factors for two unequal isotope masses M,N > Z me, me > 0; otherwise equal electronic level data; freely varying charge populations at fixed nuclear isotope inventories.',
        consequence='The adjacent-ionization translational mass ratios differ. A free-energy addition depending only on temperature and conserved isotope inventories is constant under charge redistribution and cannot change their equilibrium Saha ratio. Therefore an element-group EOS plus a composition-only entropy offset is not an exact isotope-resolved partial-ionization EOS.',
        boundary='A separately assumed fixed charge state, or deliberately identical isotope masses, removes this specific obstruction. Neither removes other physical EOS uncertainties.',
        symbolic_check=True,source_sha256={str(constants):e.c.sha(constants),str(CACHE/'source/src/ionize.f90'):e.c.sha(CACHE/'source/src/ionize.f90')},
        finite_source_constant_electron_mass_amu=electron,finite_factor_census=rows,
        census_is_physical_pressure_bound=False,physical_EOS_certified=False))
    print('ISOTOPE NO-GO',len(rows),'mass differences; largest finite factor',max(r['maximum_adjacent_stage_mass_factor_relative'] for r in rows),flush=True)


def provenance_control():
    original=json.loads((e.s.OLD/'provenance.json').read_text())
    for path,digest in original['executable_and_libraries'].items(): assert e.c.sha(Path(path))==digest,path
    p29=json.loads((e.s.v.OLD/'provenance.json').read_text())
    for rel,digest in p29['FreeEOS_source_sha256'].items(): assert e.c.sha(e.c.gr.SOURCE/rel)==digest,rel
    libraries={};sources={}
    variants=[('', 'source'),('integral-','integral-source'),('tight-integral-','tight-integral-source'),('full-integral-','full-integral-source')]
    for prefix,folder in variants:
        build=json.loads((OUT/(prefix+'build.json')).read_text());assert build['completed']
        for path,digest in build['sha256'].items(): assert e.c.sha(Path(path))==digest,path
        libraries.update(build['sha256'])
        root=CACHE/folder
        sources[str(root)]={p.relative_to(root).as_posix():e.c.sha(p) for p in root.rglob('*') if p.is_file()}
    # Each successive candidate changes only the recorded files.
    previous=e.c.gr.SOURCE
    for prefix,folder in variants:
        patch=json.loads((OUT/(prefix+'source-patch.json')).read_text())['files'];root=CACHE/folder
        for rel,digest in sources[str(root)].items():
            if rel in patch:
                assert e.c.sha(previous/rel)==patch[rel]['before'] and digest==patch[rel]['after'],rel
            else: assert e.c.sha(previous/rel)==digest,rel
        previous=root
    save('provenance.json',dict(classification='Proven',passed=True,libraries=libraries,sources=sources,
        unchanged_executable_and_libraries=original['executable_and_libraries'],
        unchanged_upstream_source_count=len(p29['FreeEOS_source_sha256']),
        compiler=subprocess.check_output(['gfortran','--version'],text=True).splitlines()[0],
        no_original_library_runtime_or_data_overwrite=True))
    print('PROVENANCE',len(libraries),'new library hashes;',sum(map(len,sources.values())),'candidate source hashes;',
        len(p29['FreeEOS_source_sha256']),'unchanged upstream sources',flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
