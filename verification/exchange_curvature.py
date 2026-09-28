"""Actual nonlinear exchange transform and saved-support component curvature."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import ctypes, json, shutil, subprocess, sys
import numpy as np
import sympy as sp
import coulomb_curvature as c

g=c.g;OUT=g.OUT/'exchange-curvature';CACHE=g.CACHE/'exchange-curvature'
NAME='free_eos_direct24_exchange_trace'
LIB=CACHE/'build/src'/('lib'+NAME+'.so.1.0.0')
FIELDS=['free','mu_over_kT','dmu_over_kT_dlnf','dmu_over_kT_dlnT',
        'ne','dlnne_dlnf','kT','free_f','exchange_ne_H_over_kT',
        'ideal_plus_exchange_ne_H_over_kT','psi','psi_prime','dpsi_prime_dlnf',
        'Legendre_mu_residual','Legendre_mu_derivative_residual','iforder']


def save(name,value): (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def read(path): return json.loads(path.read_text())


def prepare():
    assert not OUT.exists();OUT.mkdir();CACHE.mkdir(exist_ok=True)
    paths=[g.ROOT/'verification/exchange_curvature.py',c.OUT/'manifest.json',g.OUT/'reference-state.npz']
    for name in ['mod_exchange.f90','master_exchange.f90','exchange_gcpf.f90']:
        path=OUT/name;shutil.copy2(g.d.CACHE/'full-integral-source/src'/name,path);paths.append(path)
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='fe58f46',
        bindings={p.relative_to(g.ROOT).as_posix():g.c.sha(p) for p in paths},library_sha256=g.c.sha(c.s.LIB),
        cells=5735,processes=4,block_cells=128,control_cells=[0,2972,5734],
        settings=dict(morder=21,ifexchange_in=14,nonlinear_transform_iforder=2),
        old_output_bitwise_required=True,original_exchange_relative_tolerance=1e-10,
        normalized_Legendre_tolerance=1e-10,finite_log_steps=[1e-4,5e-5],finite_derivative_tolerance=1e-6,
        method='Read the cached original EOS exchange state before an independent call to the same nonlinear master_exchange and exchange_free routines. Verify the free value, chemical potential and derivative, ideal electron density, nonlinear Legendre identity and finite differences. Add the resulting scalar electron curvature to the previously frozen Coulomb and ideal-ion support.',
        omitted=['Pressure ionization','Density-dependent partition/excitation','Unresolved populations','Native continuous and physical EOS error'],
        full_EOS_Hessian_certified=False,physical_EOS_certified=False,full_GR_evolution=False))


def trace_build():
    source=CACHE/'source';assert not source.exists()
    original=g.d.CACHE/'full-integral-source';shutil.copytree(original,source);changes={}
    patches={
        'mod_exchange.f90': [('  public&\n       iforder', '  real(fp_kind), save, public :: last_fl_input=0._fp_kind, last_tl_input=0._fp_kind\n  public&\n       iforder')],
        'master_exchange.f90': [('  use mod_master_exchange_data, only:&', '  use mod_master_exchange_data, only: last_fl_input,last_tl_input\n  use mod_master_exchange_data, only:&'),
            ('  fl = max(min(fl_in, ln_overflow_limit), ln_underflow_limit)', '  last_fl_input=fl_in;last_tl_input=tl\n  fl = max(min(fl_in, ln_overflow_limit), ln_underflow_limit)')],
        'CMakeLists.txt': [('OUTPUT_NAME free_eos_direct24_integral_full','OUTPUT_NAME '+NAME)]}
    for name,pairs in patches.items():
        path=source/'src'/name;text=path.read_text()
        for old,new in pairs:
            assert old in text;(oldcount:=text.count(old));text=text.replace(old,new,1)
        path.write_text(text);shutil.copy2(path,OUT/('traced-'+name))
        changes[name]=dict(before=g.c.sha(original/'src'/name),after=g.c.sha(path))
    # Existing build helper only; output is confined to this new candidate.
    shutil.copy2(c.s.OUT/'inventory_bridge.f90',OUT/'direct_ion_bridge.f90');old=g.d.OUT
    try:
        g.d.OUT=OUT;g.d.build_at(source,CACHE/'build',NAME,CACHE/'unused-inventory.so',prefix='traced-')
    finally: g.d.OUT=old
    plan=read(OUT/'plan.json')
    plan['bindings']['verification/exchange_curvature.py']=g.c.sha(g.ROOT/'verification/exchange_curvature.py')
    plan['trace_revision']=dict(source_changes=changes,library_sha256=g.c.sha(LIB),
        preserved={p.name:g.c.sha(p) for p in OUT.glob('*cached-grouping*')},
        reason='Same returned fl is not guaranteed to equal the last input used by the cached EOS exchange state. Record actual internal inputs, compare the independent exchange call at those inputs with the unchanged original thresholds, and separately retain the original returned-input mismatch. The grouping-only change did not fix the original control failure.',
        original_failed_control_preserved=True,original_tolerances_unchanged=True)
    save('plan.json',plan)


def build():
    source=OUT/'exchange_bridge.f90'
    source.write_text('''subroutine exchange_cached(res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,cpe
  use mod_master_exchange_data, only: iforder,fexprime,n_e,t,p_e,pstarprime, &
       psiprime,psi,muex,muex2,muexf,muex2f,last_fl_input,last_tl_input
  implicit none
  real(c_double), intent(out) :: res(7)
  real(c_double) :: fex2
  fex2=boltzmann*t*(psiprime-psi)*n_e-(cpe*pstarprime(1)-p_e)
  res=[fexprime+fex2, &
       muex+muex2,muexf+muex2f,n_e,real(iforder,c_double),last_fl_input,last_tl_input]
end subroutine exchange_cached

subroutine exchange_probe(fl,tl,res) bind(C)
  use iso_c_binding
  use mod_free_eos_constants, only: boltzmann,c_e
  use mod_exchange, only: master_exchange,exchange_free
  use mod_master_exchange_data, only: iforder,psi,psiprime,dpsiprimedf,flprimef
  implicit none
  real(c_double), value :: fl,tl
  real(c_double), intent(out) :: res(16)
  real(c_double) :: r(9),p(9),ss(3),uu(3),dve,dvef,dvet,free,freef,ne,kt,psif
  call master_exchange(0,fl,tl,r,p,ss,uu,21,14,dve,dvef,dvet)
  call exchange_free(r,p,free,freef)
  ne=c_e*r(1);kt=boltzmann*exp(tl);psif=sqrt(1.d0+exp(fl))
  res=[free,-dve,-dvef,-dvet,ne,r(2),kt,freef, &
       -dvef/r(2),(psif-dvef)/r(2),psi,psiprime,dpsiprimedf*flprimef, &
       -dve-(psiprime-psi),psif-dvef-dpsiprimedf*flprimef,real(iforder,c_double)]
end subroutine exchange_probe
''')
    command=['gfortran','-O2','-fPIC','-shared','-I'+str(LIB.parent),str(source),str(c.OUT/'coulomb_bridge.f90'),
        '-L'+str(LIB.parent),'-Wl,-rpath,'+str(LIB.parent),'-l'+NAME,'-o',str(CACHE/'exchange.so')]
    result=subprocess.run(command,capture_output=True,text=True);(OUT/'build.log').write_text(result.stdout+result.stderr)
    assert result.returncode==0,result.stderr
    save('bridge.json',dict(classification='Counterexample candidate',command=command,source_sha256=g.c.sha(source),bridge_sha256=g.c.sha(CACHE/'exchange.so')))


class EOS(c.EOS):
    def __init__(self):
        super().__init__();self.exchange_lib=ctypes.CDLL(str(CACHE/'exchange.so'))
        array=np.ctypeslib.ndpointer(np.float64,flags='C_CONTIGUOUS')
        self.inventory_lib=self.exchange_lib;call=self.inventory_lib.coulomb_inventory
        call.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_double,array,array,ctypes.POINTER(ctypes.c_int)];call.restype=None
        def capture(mode,value,t,eps,out,info):
            raw=np.full(25,np.nan);call(mode,value,t,eps,raw,info)
            out[:]=raw[:22];out[20]=0.;self.molecules=raw[22:24].copy();self.fl=float(raw[24]);self.rho_native=float(raw[0])
        self.call=capture;self.probe_call=self.inventory_lib.coulomb_probe
        self.probe_call.argtypes=[ctypes.c_double]*4+[array];self.probe_call.restype=None
        self.exchange_lib.exchange_cached.argtypes=[array];self.exchange_lib.exchange_cached.restype=None
        self.exchange_lib.exchange_probe.argtypes=[ctypes.c_double,ctypes.c_double,array];self.exchange_lib.exchange_probe.restype=None

    def probe_exchange(self,fl,tl):
        result=np.empty(16);self.exchange_lib.exchange_probe(float(fl),float(tl),result)
        assert np.all(np.isfinite(result));return result

    def full(self,r,t,x):
        snap,a,H,arg,pops,mols,row=self.actual(r,t,x)
        old=np.empty(7);self.exchange_lib.exchange_cached(old)
        matched=self.probe_exchange(old[5],old[6])
        b=self.probe_exchange(arg[0],t);assert old[4]==b[15]==2
        replay=float(np.max(abs(matched[[0,1,2,4]]-old[:4])/np.maximum(abs(old[:4]),1e-100)))
        legendre=max(abs(b[13]),abs(b[14]))/max(1.,abs(b[1]),abs(b[12]))
        free_derivative=abs(b[7]/(b[6]*b[4]*b[5])-b[1])/max(1.,abs(b[1]))
        other=a.copy();other[13]=a[13]+b[2]
        margins,G=c.curvature(pops,mols,other)
        row.update(exchange_relative_error=replay,normalized_Legendre_error=float(legendre),
            normalized_free_derivative_error=float(free_derivative),electron_ne_H_over_kT=float(b[9]),
            actual_cached_input_fl_difference=float(old[5]-arg[0]),
            actual_cached_input_T_difference=float(old[6]-t),
            returned_input_exchange_relative_error=float(np.max(abs(b[[0,1,2,4]]-old[:4])/np.maximum(abs(old[:4]),1e-100))),
            same_ideal_electron_density=bool(b[4]==a[10]))
        return snap,a,b,arg,old,margins,G,row


def gates(row,plan):
    return row['exchange_relative_error']<plan['original_exchange_relative_tolerance'] and row['same_ideal_electron_density'] and max(row['normalized_Legendre_error'],row['normalized_free_derivative_error'])<plan['normalized_Legendre_tolerance']


def symbolic():
    n,kT=sp.symbols('n kT',positive=True);psi=sp.Function('psi')(n);P=sp.Function('P')
    F=n*kT*psi-P(psi)
    first=sp.diff(F,n).subs(sp.diff(P(psi),psi),n*kT)
    assert sp.simplify(first-kT*psi)==0
    second=sp.diff(first,n);assert sp.simplify(second-kT*sp.diff(psi,n))==0
    save('symbolic.json',dict(classification='Proven',passed=True,
        assumption='At fixed positive T, the corrected grand-canonical pressure Pstar(psi)=Pideal(psi)-Fexchange_grand(psi) is C2 with Pstar_second>0, and the canonical density equation Pstar_prime=kT*n has an interior root.',
        conclusion='The Legendre Helmholtz density F=n*kT*psi-Pstar(psi) has F_prime=kT*psi and F_second=(kT)^2/Pstar_second>0. Subtracting the ideal electron Legendre pair gives the canonical exchange correction; its own Hessian need not be positive.',
        scope='Conditional exact transform. The native fit, numerical transform residual, continuous positivity and physical exchange model error require separate certification.'))


def control():
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));eos=EOS();rows=[];symbolic()
    for i in plan['control_cells']:
        snap,a,b,arg,old,margins,G,row=eos.full(state['lnd'][i],state['lnT'][i],state['X'][i])
        prior=dict(np.load(c.OUT/f'control-{i}.npz'));assert np.array_equal(a,prior['values']) and np.array_equal(G,prior['Gram'])
        differences=[]
        for h in plan['finite_log_steps']:
            lo=eos.probe_exchange(arg[0]-h,arg[1]);hi=eos.probe_exchange(arg[0]+h,arg[1])
            differences.append([(hi[0]-lo[0])/(2*h)/(b[6]*b[4]*b[5]),(hi[1]-lo[1])/(2*h)])
        errors=abs(np.array(differences)-b[[1,2]])/np.maximum(1.,abs(b[[1,2]]))
        row.update(cell=i,normalized_finite_derivative_error=float(errors.max()),combined_margin=float(margins[1]))
        row['passed']=gates(row,plan) and row['normalized_finite_derivative_error']<plan['finite_derivative_tolerance'];rows.append(row)
        np.savez_compressed(OUT/f'control-{i}.npz',values=b,original=old,arguments=arg,finite=differences,margins=margins)
        print('EXCHANGE CONTROL',row,flush=True)
    save('control.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows));assert all(r['passed'] for r in rows)


def block(start):
    plan=read(OUT/'plan.json');state=dict(np.load(g.OUT/'reference-state.npz'));old=dict(np.load(c.OUT/f'block-{start}.npz'));eos=EOS()
    rows=[];values=[];margins=[];arguments=[];original=[];stop=min(start+128,plan['cells'])
    for i in range(start,stop):
        snap,a,b,arg,native,margin,G,row=eos.full(state['lnd'][i],state['lnT'][i],state['X'][i])
        assert np.array_equal(a,old['values'][i-start]) and np.array_equal(snap['eos'],old['eos'][i-start]) and np.array_equal(G,old['Gram'][i-start])
        row.update(cell=i,passed=gates(row,plan));rows.append(row);values.append(b);margins.append(margin);arguments.append(arg[:2]);original.append(native)
    path=OUT/f'block-{start}.npz';np.savez_compressed(path,values=values,margins=margins,arguments=arguments,original=original)
    record=dict(classification='Counterexample candidate',start=start,stop=stop,rows=rows,passed=all(r['passed'] for r in rows),
        output_sha256=g.c.sha(path),plan_sha256=g.c.sha(OUT/'plan.json'))
    save(f'block-{start}.json',record);assert record['passed'],start;return record


def run():
    plan=read(OUT/'plan.json');assert read(OUT/'control.json')['passed']
    for rel,digest in plan['bindings'].items(): assert g.c.sha(g.ROOT/rel)==digest,rel
    assert g.c.sha(c.s.LIB)==plan['library_sha256'] and g.c.sha(CACHE/'exchange.so')==read(OUT/'bridge.json')['bridge_sha256']
    assert g.c.sha(LIB)==plan['trace_revision']['library_sha256']
    for name,digest in plan['trace_revision']['preserved'].items(): assert g.c.sha(OUT/name)==digest,name
    rows=[];records=[]
    with ProcessPoolExecutor(max_workers=plan['processes']) as pool:
        for future in as_completed([pool.submit(block,i) for i in range(0,plan['cells'],plan['block_cells'])]):
            rec=future.result();rows+=rec['rows'];records.append(rec);print('EXCHANGE',len(rows),'/',plan['cells'],flush=True)
    records.sort(key=lambda r:r['start']);parts=[dict(np.load(OUT/f"block-{r['start']}.npz")) for r in records]
    a=np.concatenate([p['values'] for p in parts]);margins=np.concatenate([p['margins'] for p in parts])
    save('result.json',dict(classification='Counterexample candidate',passed=True,cells=len(rows),
        original_EOS_Coulomb_and_Gram_bitwise=True,maximum_exchange_replay_relative=max(r['exchange_relative_error'] for r in rows),
        maximum_normalized_Legendre_error=max(r['normalized_Legendre_error'] for r in rows),
        maximum_normalized_free_derivative_error=max(r['normalized_free_derivative_error'] for r in rows),
        electron_ne_H_over_kT_range=[float(a[:,9].min()),float(a[:,9].max())],
        minimum_ideal_ions_Coulomb_electron_exchange_margin=float(margins[:,1].min()),
        full_EOS_Hessian_certified=False,physical_EOS_certified=False,full_GR_evolution=False))
    print('EXCHANGE COMPLETE',len(rows),'minimum margin',margins[:,1].min(),flush=True)


if __name__=='__main__': globals()[sys.argv[1]]()
