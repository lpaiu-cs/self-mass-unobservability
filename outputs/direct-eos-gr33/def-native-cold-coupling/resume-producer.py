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


def support():
    assert json.loads((OUT/'logarithmic-controls.json').read_text())['passed']
    assert not (OUT/'repaired-bank.npz').exists()
    write(OUT/'support-plan.json',dict(classification='Counterexample candidate',
        decision='Apply the demonstrated native arithmetic repair to the original failed fine evolution. Extend only cold constitutive and log optical support; keep every old node and use exactly the old interpolators at and above240K.',
        new_temperature_K=[80,120,170,200],new_states=136,reused_states=918,
        budget=dict(native_calls=1000,seconds=60,CPU_threads=1,memory_GB=1),
        forecast='58 constrained native calls took2.10s.136 states at4-6 calls/state imply20-30s;60s cap includes serialization. Dense cold states remain unmeasured.',
        gates=dict(constitutive=.002,spectrum=.002,population=1e-12,relative_H=1e-8),
        stop='Stop on native or independent control failure. No floor clamp, further automatic support extension, or change to original evolution paths.',source_sha256=sha(__file__)))
    signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(60);start=time.monotonic();n=logarithmic_native(1000)
    d=dict(np.load(prior.optical.ex.OUT/'repaired-bank.npz'));s=dict(np.load(prior.OUT/'spectrum-bank.npz'));lt=np.log([80.,120.,170.,200.])
    raw=np.zeros((2,len(d['x']),4,21));lf=np.zeros((*raw.shape[:3],10));af=np.zeros(raw.shape[:3]);done=np.zeros(af.shape,bool);failure=None
    try:
        for iy,y in enumerate(d['ys']):
            for ix,x in enumerate(d['x']):
                for it,t in enumerate(lt):
                    a=n.state(float(x),float(t),float(y));raw[iy,ix,it]=a['raw'];lf[iy,ix,it]=a['log_fraction'];af[iy,ix,it]=a['affinity'];done[iy,ix,it]=True
        d['lt']=np.r_[lt,d['lt']];d['raw']=np.concatenate([raw,d['raw']],axis=2)
        # Only inherited constructor storage: prescribed-bath rates are never
        # used by the actual-photon evolution and are not cold rate evidence.
        d['rates']=np.concatenate([np.ones((*raw.shape[:3],3,2)),d['rates']],axis=2)
        np.savez_compressed(OUT/'repaired-bank.npz',**d)
        s['lt']=d['lt'];s['log_fraction']=np.concatenate([lf,np.log(s.pop('fraction'))],axis=2);s['affinity']=np.concatenate([af,s['affinity']],axis=2)
        np.savez_compressed(OUT/'spectrum-bank.npz',**s)
    except Exception as exc:failure=repr(exc)
    finally:
        np.savez_compressed(OUT/'cold-native-states.npz',**{k:np.array([r[k] for r in n.ion.states]) for k in n.ion.states[0]})
        np.savez_compressed(OUT/'cold-nodes.npz',raw=raw,log_fraction=lf,affinity=af,done=done,lt=lt)
        write(OUT/'support.json',dict(classification='Counterexample candidate',passed=failure is None,failure=failure,states=int(done.sum()),native_calls=n.ion.calls+1,seconds=time.monotonic()-start));signal.alarm(0)
    print((OUT/'support.json').read_text(),flush=True);assert failure is None


class ColdEOS(prior.optical.ex.EOS):
    def __init__(self):
        self.warm=prior.optical.ex.EOS();self.__dict__.update(self.warm.__dict__)
        self.cold=object.__new__(prior.optical.ex.EOS)
        fn=prior.optical.ex.EOS.__init__;FunctionType(fn.__code__,dict(fn.__globals__,OUT=OUT))(self.cold)
        self.lt=self.cold.lt

    def evaluate(self,rho,lt):
        cold=lt<self.warm.lt[0];result=np.zeros((7,len(rho)));y=np.broadcast_to(self.y,rho.shape)
        for mask,model in [(cold,self.cold),(~cold,self.warm)]:
            if mask.any():model.y=y[mask];result[:,mask]=np.array(model.evaluate(rho[mask],lt[mask]))
        return tuple(result)

    def reactions(self,*_):raise AssertionError('Use the actual coupled photon spectrum')


class ColdSpectrum(prior.optical.Spectrum):
    def __init__(self):
        super().__init__();self.cold=object.__new__(prior.optical.Spectrum);self.cold.__dict__.update(self.__dict__)
        c=self.cold;c.d=d=np.load(OUT/'spectrum-bank.npz');c.f=[];c.rev=[]
        for k,y in enumerate(self.ys):
            T=np.exp(d['lt'])[None,:,None];lf=d['log_fraction'][k]
            pref=lf+(self.binding[0]-self.binding)[None,None,:]/(prior.optical.ex.K*T)
            reverse=lf+np.log(y/(1-y))+d['affinity'][k,:,:,None]-self.binding[None,None,:]/(prior.optical.ex.K*T)
            c.f.append([prior.optical.RectBivariateSpline(d['x'],d['lt'],pref[:,:,j]) for j in range(10)])
            c.rev.append([prior.optical.RectBivariateSpline(d['x'],d['lt'],reverse[:,:,j]) for j in range(10)])

    def levels(self,rho,lt,y):
        cold=lt<self.d['lt'][0];f=np.zeros((len(rho),10));r=np.zeros_like(f)
        for mask,model in [(cold,self.cold),(~cold,self)]:
            if mask.any():f[mask],r[mask]=prior.optical.Spectrum.levels(model,rho[mask],lt[mask],y[mask])
        return f,r


def support_controls():
    assert not (OUT/'support-controls.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(30)
    write(OUT/'support-controls-plan.json',dict(classification='Counterexample candidate',native_calls=180,seconds=30,
        states=[[-17.7,90,2e-14],[-16.5,140,5e-9],[-9,180,.0003],[-.1,220,.0007]],
        quantity='Native P,u and separately positive absorption/spontaneous/stimulated emission integrated against the retained frequency measures and three test spectra. Also check the exact saved failure primitive and old-support bitwise reuse.',
        gate=.002,source_sha256=sha(__file__)))
    n=logarithmic_native(180);e=ColdEOS();s=ColdSpectrum();d=np.load(prior.optical.photons.OUT/'bank-16-8.npz');a=d['face_a'][-1];E=d['Einf']/a;num=d['num']/a**3;rows=[];failure=None
    try:
        for x,T,y in [(-17.7,90,2e-14),(-16.5,140,5e-9),(-9,180,.0003),(-.1,220,.0007)]:
            t=np.log(T);z=n.state(x,t,y);e.y=np.array([y]);p,u,*_=e(np.array([np.exp(x)]),np.array([t]));constitutive=max(abs(p[0]*e.rho0*prior.C**2/z['raw'][1]-1),abs(u[0]*prior.C**2/z['raw'][2]-1))
            ab=np.zeros_like(E);em=np.zeros_like(E)
            for j in range(10):
                sigma=s.cross(E,j+1);ab+=np.exp(np.log(y)+z['log_fraction'][j])*sigma
                good=sigma>0
                em[good]+=np.exp(np.log(y)+z['log_fraction'][j]+z['affinity']-E[good]/(prior.optical.ex.K*T))*sigma[good]
            aa,ee=s.coefficients(np.array([e.rho0*np.exp(x)]),np.array([t]),np.array([y]),E[None,:]);errors=[]
            for Trad in [13400,20000,80000]:
                occ=1/np.expm1(E/(prior.optical.ex.K*Trad))
                for true,est,field in [(ab,aa[0],occ),(em,ee[0],occ),(em,ee[0],np.ones_like(E))]:
                    for factor in [num,num*E]:
                        weight=factor*field;errors.append(float(np.sum(abs(est-true)*weight)/max(float(true@weight),1e-280)))
            rows.append(dict(x=x,T=T,y=y,constitutive=float(constitutive),spectrum=max(errors)))
        assert max(max(z['constitutive'],z['spectrum']) for z in rows)<.002
        # The saved complete fine snapshot must still yield identical warm
        # primitives, while the exact previously failed state must now recover.
        f=prior.Flow(896);saved=np.load(prior.OUT/'cells-896-steps-128.npz');old=f.primitive(saved['snapshot_U'][-1]);f.eos=ColdEOS();new=f.primitive(saved['snapshot_U'][-1]);warm=max(float(np.max(abs(a-b))) for a,b in zip(old,new));assert warm==0
        failed=np.load(prior.OUT/'root-failure-state.npz');V=f.primitive(failed['U']);Troot=float(np.exp(V[2][416]));assert 160<Troot<200
    except Exception as exc:failure=repr(exc)
    result=dict(classification='Counterexample candidate',passed=failure is None,failure=failure,checks=rows,native_calls=n.ion.calls+1,seconds=time.monotonic()-start)
    if failure is None:result.update(warm_primitive_bitwise=True,failed_cell_temperature_K=Troot)
    write(OUT/'support-controls.json',result);signal.alarm(0);print(json.dumps(result),flush=True);assert failure is None


def resumed_runner():
    src=textwrap.dedent(inspect.getsource(prior.Coupled.run))
    start=src.index('    assert not ');end=src.index('    p0=')
    src=src[:start]+'''    assert not (OUT/(label+'.npz')).exists();start=time.monotonic();b=self.bulk;f=self.flow;m=self.m;h=prior.END/steps;count=steps
    saved=np.load(prior.OUT/'cells-896-steps-128.npz');assert saved['snapshot_t'][-1]==104*h
    U=saved['snapshot_U'][-1].copy();I=saved['snapshot_I'][-1].copy();xb=saved['snapshot_bulk_I'][-1]/b.scale
    theta=saved['snapshot_theta'][-1].copy();eta=saved['snapshot_eta'][-1].copy();u=b.eos.gas(theta,eta)[1]
    refU=U.copy();refI=I.copy();refb=xb.copy();refu=u.copy()
    ledger=np.zeros(6);discard=np.zeros(4);boundary=np.zeros(4);balance=0.;gas_balance=0.;local_balance=0.;iterations=0;residual=0.;substeps=0;snapshots=[];failure=None
    keys=['t','bulk_trace','atmosphere_trace','outside_mass','maximum_speed','surface_luminosity_per_frequency']
    history=[{k:saved[k][j] for k in keys} for j in range(104)]
    np.savez_compressed(OUT/'restart104.npz',U=refU,I=refI,bulk_I=refb*b.scale,u=refu,theta=theta,eta=eta,t=104*h)
''' + src[end:]
    src=src.replace('record(0.)','record(104*h)')
    src=src.replace('for j in range(count):','for j in range(104,count):')
    src=src.replace('(xb*b.scale-b.initial)','((xb-refb)*b.scale)').replace('(I.sum(0)-self.initial_I)','(I.sum(0)-refI.sum(0))')
    src=src.replace('b.gas_weight@(u-b.u0)','b.gas_weight@(u-refu)')
    for k in [0,2,3]:src=src.replace(f'(U[{k}]-f.initial[{k}])',f'(U[{k}]-refU[{k}])')
    src=src.replace('np.sum(I[1]*self.energy_weight)','np.sum((I[1]-refI[1])*self.energy_weight)')
    src=src.replace("                if (j+1)%max", "                if (j+1)%max") if False else src
    anchor='            if (j+1)%max'
    src=src.replace(anchor,'''            if j==108:
                elapsed=time.monotonic()-start;forecast=elapsed+(127-j)*(elapsed/5)*1.5+15
                write(OUT/'remaining-measured-budget.json',dict(first_five_seconds=elapsed,forecast_total_seconds=forecast,cap_seconds=190,eligible=forecast<190))
                assert forecast<190,'Remaining measured trajectory budget'
            if (j+1)%4==0 or j==108:print(json.dumps(dict(step=j+1,Tmin_K=float(np.exp(f.primitive(U)[2][U[0]>=f.eos.floor]).min()),seconds=time.monotonic()-start)),flush=True)
''' + anchor)
    src=src.replace('coupled_photons_and_moving_atmosphere=failure is None,','restart_step=104,ledger_scope="Segment104-128 only; earlier conservation remains separately certified",coupled_photons_and_moving_atmosphere=failure is None,')
    (OUT/'expanded-resume.py').write_text(src)
    ns=dict(vars(prior),OUT=OUT,prior=prior);exec(compile(src,__file__,'exec'),ns);return ns['run']


def resume():
    assert json.loads((OUT/'support-controls.json').read_text())['passed']
    assert not (OUT/'resume-plan.json').exists()
    write(OUT/'resume-plan.json',dict(classification='Counterexample candidate',
        decision='Finish the originally registered fine1024-fluid-cell trajectory using the repaired native cold EOS and logarithmic optical support, then adjudicate original2percent spatial gates against retained coarse results.',
        restart='Only complete snapshot104 is used. Raw failed terminal state is not a checkpoint. Histories0-104 are reused;104-128 is one new segment with its own conservative baseline and ledgers. Earlier accepted prefix certificate is retained separately.',
        budget=dict(seconds=190,CPU_threads=1,memory_GB=3,macroscopic_steps=24,native_evolution_calls=4),
        forecast='The saved104-108 replay took21.78s, about5.45s per step.24 steps predict131s, with15s serialization and1.3x execution margin below190s. First5 steps are retained and measured; stop if their forecast exceeds cap.',
        overlap='Compare recomputed105-108 bulk/atmosphere/mass with old saved histories before interpreting the new suffix. Gate1e-8 of each old peak response; old interpolators remain exact for warm states.',
        gates=json.loads((prior.OUT/'connection-plan.json').read_text())['gates'],
        stop='No rerun of successful coarse paths, new grid, extended horizon or post-result gate relaxation. Stop on positivity, primitive, support, conservation, overlap or budget failure.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'repaired-bank.npz',OUT/'spectrum-bank.npz',prior.OUT/'cells-896-steps-128.npz']}))
    signal.signal(signal.SIGALRM,prior.optical.timeout);signal.alarm(190)
    model=prior.Coupled(896,8);model.flow.eos=ColdEOS();model.spectrum=ColdSpectrum()
    row=resumed_runner()(model,128,'resumed-896-128');signal.alarm(0)
    old=np.load(prior.OUT/'cells-896-steps-128.npz');fine=np.load(OUT/'resumed-896-128.npz');coarse=np.load(prior.OUT/'cells-448-steps-128.npz');overlap={};errors={}
    for k in ['bulk_trace','atmosphere_trace','outside_mass']:
        last=min(109,len(fine[k]));overlap[k]=float(np.max(abs(fine[k][104:last]-old[k][104:last]))/max(float(np.max(abs(old[k]))),1.))
    if row['passed']:
        for k in ['atmosphere_trace','outside_mass']:errors[k]=float(np.max(abs(fine[k]-coarse[k]))/max(float(np.max(abs(fine[k]))),1.))
    result=dict(classification='Counterexample candidate',passed=bool(row['passed'] and max(overlap.values())<1e-8 and max(errors.values(),default=1)<.02),
        continuation=row,overlap=overlap,space_comparison=errors,reused_time_comparison=json.loads((prior.OUT/'time-comparison.json').read_text()),
        cold_failure_fixed_in_actual_coupled_evolution=bool(row['completed_steps']>109 and row['failure'] is None),full_GR_scalar_feedback=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
