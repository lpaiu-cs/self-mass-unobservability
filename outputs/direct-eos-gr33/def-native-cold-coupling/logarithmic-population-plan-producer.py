"""Repair the saved cold native state before continuing its coupled trajectory."""
from pathlib import Path
import ctypes
import json
import signal
import sys
import time
import shutil
import subprocess
from types import FunctionType, MethodType
import inspect
import textwrap
import numpy as np
import def_native_two_way_atmosphere as prior

OUT=prior.OUT.parent/'def-native-cold-coupling'
write=prior.write;sha=prior.sha
COLD=prior.optical.ex.old.cold
CACHE=COLD.CACHE.parent/'native-cold-coupling'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='1297eba71',
        claim='Repair the exact native EOS state blocking the real two-way photon/moving-atmosphere fine trajectory, then apply the validated input to its remaining original interval.',
        reuse='Original448/64 and448/128 completed paths, fine104 complete snapshot and fine108 history, original918 native states and the saved failed conservative cell. No full-path replay or added spatial/time grid.',
        first_decision='Use exposed native molecular diagnostics and constrained states to determine why finite requested density underflows internally at160K. Preserve failed attempts. No ideal-EOS substitution or temperature clamp.',
        first_budget=dict(native_calls=60,seconds=20,CPU_threads=1,memory_GB=1,fluid_steps=0),
        forecast='Phase112 seven native endpoint calls took1.75s; this bounded single-cell diagnostic is expected1-10s. If native/build/support changes are needed, assess and register those before execution.',
        gates=dict(constitutive=.002,population=1e-12,relative_H=1e-8,primitive=2e-11,energy=1e-8,space_trace=.02,space_mass=.02),
        stop='Do not enlarge the horizon, add full flow paths, relax gates or silently extrapolate an EOS or spectrum. A recovered state must actually enter the original coupled evolution before calling its blocking condition solved.',
        bindings={str(p):sha(p) for p in [Path(__file__),prior.OUT/'root-failure-state.npz',prior.OUT/'cells-896-steps-128.npz',prior.OUT/'result.json',prior.OUT/'spectrum-bank.npz',prior.optical.ex.OUT/'repaired-bank.npz']}))


def diagnose():
    assert not (OUT/'diagnostic.json').exists();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(20);start=time.monotonic()
    native=prior.optical.ex.Native(cap=60);f=prior.Flow(896);m=f.base;U=np.load(prior.OUT/'root-failure-state.npz')['U'][:,416]
    D=U[0];tau=(U[2]-(m.a[416]-m.a0)*f.eos.cx*D)/m.a[416];v=U[1]/(f.eos.cx*D+tau);rho=D*np.sqrt(1-v*v);x=float(np.log(rho));y=float(U[3]/D)
    rows=[];fields=[];lib=native.ion.gas.gas_lib
    def array(name,n):return np.ctypeslib.as_array((ctypes.c_double*n).in_dll(lib,'__mod_nuvar_MOD_'+name)).copy().tolist()
    try:
        for T in [240.,200.,170.,160.]:
            failure=None
            try:state=native.state(x,float(np.log(T)),y)
            except Exception as exc:failure=repr(exc)
            row=dict(T=T,failure=failure,native_calls=native.ion.calls,molecular_log_populations=array('mol_logs',3),molecular_equilibria=array('mol_eq',2),molecular_fields=array('mol_dv',2),fields_H_H2_H2plus=native.ion.fields[[0,316,317]].tolist())
            if failure is None:row.update(pressure=float(state['raw'][1]),energy=float(state['raw'][2]),population_error=state['population_error'])
            rows.append(row);fields.append(native.ion.fields.copy());print(json.dumps(row),flush=True)
    finally:
        np.savez_compressed(OUT/'diagnostic-states.npz',fields=fields,temperatures=[r['T'] for r in rows],x=x,y=y,rho=rho,U=U)
        write(OUT/'diagnostic.json',dict(classification='Counterexample candidate',checks=rows,native_calls=native.ion.calls,seconds=time.monotonic()-start,source_sha256=sha(__file__)))
        signal.alarm(0)


def diagnostic_build():
    assert not CACHE.exists();CACHE.mkdir()
    write(OUT/'partition-build-plan.json',dict(classification='Counterexample candidate',
        evidence='At200K the molecular equilibrium and internal H2plus field become infinite while its external constrained field remains finite. The failure is upstream of density recovery.',
        decision='Read the exponential and partition operands at their owner; correct only their demonstrated arithmetic cause and test old supported states before using new support.',
        budget=dict(build_seconds=60,diagnostic_native_calls=8,diagnostic_seconds=20),
        forecast='The same excitation module previously compiled in3.59s. Reuse all unchanged native objects and build a separate diagnostic library; no old library is overwritten.',
        source_sha256=sha(__file__)))
    for name in ['mod_excitation.f90','excitation_pi.f90','excitation_sum.f90']:shutil.copyfile(COLD.CACHE/name,CACHE/name)
    p=CACHE/'excitation_pi.f90';src=p.read_text();anchor='                 mu(ion) = exparg*qstar(nmin_s,izqstar)';assert src.count(anchor)==1
    src=src.replace(anchor,"                 if(tl.lt.log(240._fp_kind)) write(*,*) 'COLD_EXP',iz,nmin_s,izqstar,exparg,qstar(nmin_s,izqstar),c2t*bion(ion),plop(ion),qh2plus\n"+anchor);p.write_text(src)
    build_library('diagnostic')


def build_library(label):
    old=json.loads((COLD.OUT/'stable-build.json').read_text())['receipts'];start=time.monotonic();receipts=[]
    library='free_eos_native_cold_coupling_'+label
    compile_cmd=old[0]['command'].copy();compile_cmd=[v.replace(str(COLD.CACHE),str(CACHE)).replace('mod_excitation-stable.o','mod_excitation-'+label+'.o') for v in compile_cmd]
    link_cmd=old[1]['command'].copy();link_cmd=[v.replace(str(COLD.CACHE/'libfree_eos_native_cold_stable.so'),str(CACHE/('lib'+library+'.so'))).replace('libfree_eos_native_cold_stable.so','lib'+library+'.so').replace('mod_excitation-stable.o','mod_excitation-'+label+'.o') for v in link_cmd]
    (CACHE/label).mkdir(exist_ok=True)
    bridge=COLD.old.native.OUT.parent/'gr-radiation-eos-split/gas-bridge.f90'
    bridge_cmd=['gfortran','-O2','-fPIC','-shared','-I'+str(COLD.old.BUILD),str(bridge),'-L'+str(CACHE),'-Wl,-rpath,'+str(CACHE),'-l'+library,'-o',str(CACHE/label/'gas.so')]
    for i,cmd in enumerate([compile_cmd,link_cmd,bridge_cmd]):
        row=subprocess.run(cmd,cwd=CACHE,capture_output=True,text=True,timeout=max(1,60-(time.monotonic()-start)))
        (OUT/f'{label}-build-{i}.log').write_text(row.stdout+row.stderr);receipts.append(dict(command=cmd,returncode=row.returncode));assert row.returncode==0,row.stderr
    write(OUT/(label+'-build.json'),dict(classification='Counterexample candidate',seconds=time.monotonic()-start,receipts=receipts,
        bindings={str(p):sha(p) for p in [CACHE/'excitation_pi.f90',CACHE/'excitation_sum.f90',CACHE/(label+'/gas.so'),CACHE/('lib'+library+'.so')]}))


def native_variant(label,cap):
    native=prior.optical.ex.Native(cap=cap)
    ion=object.__new__(prior.optical.ex.old.InventoryIons)
    init=COLD.old.Ions.__init__
    FunctionType(init.__code__,dict(init.__globals__,CACHE=CACHE/label))(ion,cap-1)
    # The physical baseline is unchanged; retain the initial constructor call
    # in resource accounting and verify the replacement surface evaluation.
    baseline=ion.snapshot(native.lr,np.log(native.fan.T),np.zeros(318))
    assert np.max(abs(baseline['eos']-native.base['eos'])/np.maximum(abs(native.base['eos']),1.))<1e-10
    native.ion=ion;native.variant_initial_calls=1
    return native


def diagnostic_native():
    signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(20);start=time.monotonic();native=native_variant('diagnostic',8)
    d=np.load(OUT/'diagnostic-states.npz');failure=None
    try:native.state(float(d['x']),np.log(200.),float(d['y']))
    except Exception as exc:failure=repr(exc)
    write(OUT/'partition-diagnostic.json',dict(classification='Counterexample candidate',failure=failure,native_calls=native.ion.calls+native.variant_initial_calls,seconds=time.monotonic()-start));signal.alarm(0)


def repair_build():
    import difflib
    assert not (OUT/'scaled-build.json').exists()
    write(OUT/'scaled-partition-plan.json',dict(classification='Counterexample candidate',
        failure='The200K diagnostic measured the H2plus Rydberg partition as infinity while its multiplying Boltzmann factor was2.60517e-150. Evaluating exp(a/n^2) before multiplying the compensating factor overflows even though the physical product is finite.',
        repair='Represent each cold partition tail and every first/second partial derivative by the same positive exponential scale. Combine its log scale with the caller Boltzmann factor before exponentiating. Leave ordinary supported evaluations bitwise unchanged where no scale is required. Retain the same native Rydberg levels, occupation law and constraints.',
        derivative='The returned jets are exp(-s) times the derivatives of the original partition at the evaluation point, not derivatives of the scaled representation; derivative ratios and physical products therefore recover the original jets exactly in real arithmetic.',
        summation='Use the existing direct quantum-number sum for tails whose leading exponent exceeds680. Each requested lower level has its own scale; never zero a sibling tail merely because another tail is larger. Use the original termination test in log units and the same occupation cutoff.',
        budget=dict(build_seconds=60,controls_native_calls=100,controls_seconds=30),
        gates=dict(warm_constitutive=1e-10,cold_population=1e-12,cold_relative_H=1e-8),
        stop='No physical EOS parameter, quantum level, temperature floor or acceptance tolerance changes. Stop before support extension if the exact cold states or old controls fail.',source_sha256=sha(__file__)))
    for name in ['mod_excitation.f90','excitation_pi.f90','excitation_sum.f90']:shutil.copyfile(COLD.CACHE/name,CACHE/name)
    for name in ['qstar_calc.f90','qryd_calc.f90']:shutil.copyfile(COLD.old.SOURCE/name,CACHE/name)
    changes={}
    p=CACHE/'mod_excitation.f90';old=p.read_text();src=old.replace('contains','  real(fp_kind), save :: qstar_logscale(10,29)=0._fp_kind, qryd_logscale=0._fp_kind\ncontains',1)
    helpers='''
  function partition_product(logfactor, partition, logscale) result(value)
    real(fp_kind), intent(in) :: logfactor, partition, logscale
    real(fp_kind) :: value
    if(partition.le.0._fp_kind) then
       value=0._fp_kind
    elseif(logscale.eq.0._fp_kind.and.abs(logfactor).lt.700._fp_kind) then
       value=exp(logfactor)*partition
    else
       value=exp(logfactor+logscale+log(partition))
    endif
  end function partition_product

  subroutine scaled_summand(rn2,arg,lnocc,value,first,second)
    real(fp_kind), intent(in) :: rn2,arg,lnocc
    real(fp_kind), intent(out) :: value,first,second
    real(fp_kind) :: factor,em
    if(arg.gt.40._fp_kind) then
       factor=exp(arg+lnocc-qryd_logscale);em=exp(-arg)
       value=2._fp_kind*rn2*factor*(1._fp_kind-(1._fp_kind+arg)*em)
       first=2._fp_kind*factor*(1._fp_kind-em)
       second=2._fp_kind*factor/rn2
    else
       call plsummand_normalized(rn2,arg,value,first,second)
       factor=exp(lnocc-qryd_logscale)
       value=value*factor;first=first*factor;second=second*factor
    endif
  end subroutine scaled_summand
'''
    src=src.replace('end module mod_excitation',helpers+'end module mod_excitation');p.write_text(src);changes[p.name]=(old,src)
    for name in ['excitation_pi.f90','excitation_sum.f90']:
        p=CACHE/name;old=p.read_text();src=old.replace('exparg = exp(exparg)','! Keep the Boltzmann factor logarithmic until the partition product.')
        for index in ['iz','izqstar']:
            src=src.replace(f'exparg*qstar(nmin_s,{index})',f'partition_product(exparg,qstar(nmin_s,{index}),qstar_logscale(nmin_s,{index}))')
        src=src.replace('exparg*qmhd_he1','partition_product(exparg,qmhd_he1,0._fp_kind)')
        for index in ['iz','izhi+1']:
            anchor=f'qstarx2(:,:,:nmin_max_max,{index}))';assert src.count(anchor)==1
            src=src.replace(anchor,f'qstarx2(:,:,:nmin_max_max,{index}),qstar_logscale(:nmin_max_max,{index}))')
        assert 'exparg*qstar' not in src and src!=old;p.write_text(src);changes[name]=(old,src)
    p=CACHE/'qstar_calc.f90';old=p.read_text();src=old
    anchor='     qstar, qstart, qstarx, qstart2, qstartx, qstarx2)';assert src.count(anchor)==1
    src=src.replace(anchor,anchor[:-1]+', logscale)')
    src=src.replace('  ! Internal variables','  real(fp_kind), intent(out), optional :: logscale(:)\n\n  ! Internal variables',1)
    src=src.replace('  integer nmax_reached, n, ib, nx, nmin_max','  integer nmax_reached, n, ib, nx, nmin_max, ntarget')
    src=src.replace('  real(fp_kind), allocatable :: b(:), bconst(:)','  real(fp_kind), allocatable :: b(:), bconst(:), sq(:), st(:), sx(:,:), stt(:), stx(:,:), sxx(:,:,:)')
    anchor='  if(ifapprox) then';assert src.count(anchor)==1
    src=src.replace(anchor,'''  if(present(logscale)) logscale=0._fp_kind
  qryd_logscale=0._fp_kind
  if(present(logscale).and.a/real(nmin*nmin,fp_kind).gt.680._fp_kind) then
     allocate(sq(nmin_max),st(nmin_max),sx(nb,nmin_max),stt(nmin_max),stx(nb,nmin_max),sxx(nb,nb,nmin_max))
     qstar=0._fp_kind;qstart=0._fp_kind;qstarx=0._fp_kind
     qstart2=0._fp_kind;qstartx=0._fp_kind;qstarx2=0._fp_kind
     do ntarget=nmin,nmin_max
        logscale(ntarget)=max(0._fp_kind,a/real(ntarget*ntarget,fp_kind)-400._fp_kind)
        qryd_logscale=logscale(ntarget)
        call qryd_calc(ifpl,ifmhd,ifneutral,eps_factor,ntarget,nmax,a,b,sq,st,sx,stt,stx,sxx,nmax_reached)
        qstar(ntarget)=sq(ntarget);qstart(ntarget)=st(ntarget);qstarx(:,ntarget)=sx(:,ntarget)
        qstart2(ntarget)=stt(ntarget);qstartx(:,ntarget)=stx(:,ntarget);qstarx2(:,:,ntarget)=sxx(:,:,ntarget)
     enddo
     qryd_logscale=0._fp_kind
  elseif(ifapprox) then''')
    p.write_text(src);changes[p.name]=(old,src)
    p=CACHE/'qryd_calc.f90';old=p.read_text();src=old
    src=src.replace('summand.gt.eps/eps_factor',"(summand.gt.0._fp_kind.and.log(max(summand,tiny(summand)))+qryd_logscale.gt.log(eps)-log(eps_factor))")
    src=src.replace('     if(ifpl) then\n        call plsummand_normalized', '''     if(qryd_logscale.gt.0._fp_kind) then
        if(ifpl) then
           call scaled_summand(rn2,arg,0._fp_kind,summand,summanda,summanda2)
        else
           summand=2._fp_kind*rn2*exp(arg-qryd_logscale)
           summanda=summand/rn2;summanda2=summanda/rn2
        endif
     elseif(ifpl) then
        call plsummand_normalized''')
    anchor='''        occupation = exp(lnoccupation)
        summand = summand*occupation
        summanda = summanda*occupation
        summanda2 = summanda2*occupation''';assert src.count(anchor)==1
    src=src.replace(anchor,'''        if(qryd_logscale.gt.0._fp_kind) then
           if(ifpl) then
              call scaled_summand(rn2,arg,lnoccupation,summand,summanda,summanda2)
           else
              summand=2._fp_kind*rn2*exp(arg+lnoccupation-qryd_logscale)
              summanda=summand/rn2;summanda2=summanda/rn2
           endif
        else
'''+anchor+'''
        endif''')
    p.write_text(src);changes[p.name]=(old,src)
    for name,(old,src) in changes.items():
        (OUT/(name+'.patch')).write_text(''.join(difflib.unified_diff(old.splitlines(True),src.splitlines(True),fromfile='prior/'+name,tofile='scaled/'+name)))
        shutil.copyfile(CACHE/name,OUT/('scaled-'+name))
    build_library('scaled')


def repaired_controls():
    assert not (OUT/'scaled-controls.json').exists();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(30);start=time.monotonic()
    native=native_variant('scaled',100);d=np.load(OUT/'diagnostic-states.npz');rows=[];failure=None
    try:
        for T in [240.,200.,170.,160.,120.,80.]:
            state=native.state(float(d['x']),float(np.log(T)),float(d['y']));rows.append(dict(T=T,pressure=float(state['raw'][1]),energy=float(state['raw'][2]),population_error=state['population_error']))
            print(json.dumps(rows[-1]),flush=True)
    except Exception as exc:failure=repr(exc)
    result=dict(classification='Counterexample candidate',passed=failure is None,checks=rows,failure=failure,native_calls=native.ion.calls+native.variant_initial_calls,seconds=time.monotonic()-start)
    write(OUT/'scaled-controls.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def logarithmic_build():
    assert not (OUT/'logarithmic-build.json').exists()
    write(OUT/'logarithmic-population-plan.json',dict(classification='Counterexample candidate',
        evidence='Scaled partition evaluation recovered240,200,170,160,120K. At80K the EOS constraint succeeds but the Python optical reader asserts because native excitation diagnostics omit a smaller-than-representable H excitation term.',
        repair='Export the H excited/ground log partition ratio before native floating-point suppression. Reconstruct logarithmic finite-level fractions and combine absorption/emission exponents before exponentiation. Keep the actual constrained EOS and its physical free energy unchanged.',
        budget=dict(build_seconds=60,remaining_control_native_calls=68,control_seconds=30),
        acceptance='Old retained states agree within1e-10; actual failed native energy root within2e-11; full constrained inventories1e-12 and relative H1e-8. A log optical representation is required before any cold table extension.',
        source_sha256=sha(__file__)))
    p=CACHE/'mod_excitation.f90';src=p.read_text();src=src.replace('contains\n','  real(fp_kind), save, public :: native_h_log_ratio=0._fp_kind\ncontains\n',1);p.write_text(src)
    for name,quantity in [('excitation_pi.f90','mu(ion)'),('excitation_sum.f90','qratio')]:
        p=CACHE/name;src=p.read_text();anchor=f'if({quantity}.gt.exp_lim) then';at=src.index(anchor)
        statement='''if(ion.eq.1.and.qstar(nmin_s,iz).gt.0._fp_kind) native_h_log_ratio= &
                exparg+qstar_logscale(nmin_s,iz)+log(qstar(nmin_s,iz))+ &
                log(real(iqion(ion),fp_kind)/real(iqneutral(ielement),fp_kind))-exparg_shift
           '''
        src=src[:at]+statement+src[at:];p.write_text(src)
    for name in ['mod_excitation.f90','excitation_pi.f90','excitation_sum.f90','qstar_calc.f90','qryd_calc.f90']:shutil.copyfile(CACHE/name,OUT/('logarithmic-'+name))
    build_library('logarithmic')


def log_metadata(native,a,err,x,lt,y):
    lib=native.ion.gas.gas_lib
    lr=ctypes.c_double.in_dll(lib,'__mod_excitation_MOD_native_h_log_ratio').value
    xaux=np.ctypeslib.as_array((ctypes.c_double*5).in_dll(lib,'__mod_excitation_block_MOD_x')).copy()
    terms=np.zeros(30);native.levels(1,float(lt),xaux,terms);tail=terms.reshape(3,10)[0];terms=tail-np.r_[tail[1:],0.]
    assert np.all(terms>=0) and terms[1:].sum()>0 and np.isfinite(lr)
    L=float(np.logaddexp(0.,lr));logfraction=np.full(10,-np.inf);logfraction[0]=-L
    logfraction[1:]=lr-L+np.log(terms[1:])-np.log(terms[1:].sum())
    fraction=np.exp(logfraction);assert abs(fraction.sum()-1)<1e-13
    return dict(raw=a['eos'],fraction=fraction,log_fraction=logfraction,affinity=float(native.ion.fields[0]),population_error=err,x=x,lt=lt,y=y,L=L)


def logarithmic_native(cap):
    native=native_variant('logarithmic',cap)
    src=textwrap.dedent(inspect.getsource(prior.optical.ex.Native.state));src=src[:src.index('    lib=self.ion.gas.gas_lib')]+"    return log_metadata(self,a,err,x,lt,y)\n"
    ns=dict(vars(prior.optical.ex),log_metadata=log_metadata);exec(compile(src,__file__,'exec'),ns)
    native.state=MethodType(ns['state'],native);return native


def logarithmic_controls():
    from scipy.optimize import brentq
    assert not (OUT/'logarithmic-controls.json').exists();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(30);start=time.monotonic()
    n=logarithmic_native(68);d=np.load(OUT/'diagnostic-states.npz');f=prior.Flow(896);m=f.base;U=d['U'];D=U[0];tau=(U[2]-(m.a[416]-m.a0)*f.eos.cx*D)/m.a[416]
    rows=[];failure=None;root=None
    try:
        for T in [240.,160.,120.,80.]:
            state=n.state(float(d['x']),np.log(T),float(d['y']));rows.append(dict(T=T,population_error=state['population_error'],log_fraction=state['log_fraction'].tolist(),pressure=float(state['raw'][1]),energy=float(state['raw'][2])))
        reference=json.loads((OUT/'diagnostic.json').read_text())['checks'][0]
        warm=max(abs(rows[0]['pressure']/reference['pressure']-1),abs(rows[0]['energy']/reference['energy']-1));assert warm<1e-10
        def energy(T):
            pp=0.
            for _ in range(2):
                v=U[1]/(f.eos.cx*D+tau+pp);r=np.sqrt(1-v*v);rho=D*r;wm=v*v/(r*(1+r));state=n.state(float(np.log(rho)),np.log(T),float(d['y']));pp=state['raw'][1]/(f.eos.rho0*prior.C**2)
            u=state['raw'][2]/prior.C**2
            return float((f.eos.cx*D*wm+(rho*u+pp)/(1-v*v)-pp-tau)/tau)
        root=float(brentq(energy,160,240,xtol=1e-5));root_error=abs(energy(root));assert root_error<2e-11
    except Exception as exc:failure=repr(exc)
    result=dict(classification='Counterexample candidate',passed=failure is None,failure=failure,native_root_temperature_K=root,checks=rows,
        root_relative_energy_residual=root_error if root else None,native_calls=n.ion.calls+n.variant_initial_calls,seconds=time.monotonic()-start)
    write(OUT/'logarithmic-controls.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
