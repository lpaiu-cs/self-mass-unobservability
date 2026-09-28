"""Request27: direct native reaction outputs, with unchanged scientific binary.

Counterexample candidate: debugger extraction is checked against untraced exports;
it is not a new EOS, a GR evolution solver, or a global error certificate.
"""
from pathlib import Path
import ctypes as ct
import json, os, shutil, signal, struct, subprocess, sys, time
import re
import numpy as np
import fresh_microphysics as fresh
from thermal_restart import sha
from thermal_wd import mesa

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/native-closure27'
CACHE=Path('/home/lpaiu/work/native-closure27')
OLD=ROOT/'outputs/remaining-closure26'
ISOS=list(fresh.ISOS)

def save(name, obj):
    (OUT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

def context():
    fresh.OUT=OUT;fresh.CACHE=CACHE;fresh.ISOS=list(ISOS)

def prepare():
    assert not OUT.exists() and not CACHE.exists()
    OUT.mkdir();CACHE.mkdir();context()
    save('plan.json',dict(classification='Counterexample candidate',checkpoint='3b70803',
        previous_manifest_sha256=sha(OLD/'manifest.json'),
        stages=['Direct native dxdt/Jacobian extraction at the original GR input',
                'Independent comparison with channel reconstruction and untraced profile',
                'Generated-fluorine common native EOS and composition continuation controls'],
        tolerances=dict(input_absolute=1e-12,profile_relative=1e-10,
                        vector_relative_to_zone_max=1e-8,baryon_relative=1e-12),
        constraints='No executable file changes, no runtime rebuild, no new empirical observations. '
        'A disposable single-thread child is stopped at the net_get entry/return; only breakpoint '
        'instruction bytes and instruction pointer are restored. Numerical argument/result memory is read only. '
        'Retain failed gates. Initial state evaluation is not stellar evolution.'))
    setup('direct',dict(np.load(ROOT/'outputs/fresh-microphysics25/gr-input.npz')))

def setup(label,data,extra='',species=None,network=None):
    context();h,_=mesa(fresh.SOURCE)
    if species is not None: fresh.ISOS=list(species)
    fresh.setup_run(label,data,float(h['star_age']))
    folder=CACHE/label
    if species is not None:
        save(label+'-species.json',species)
        p=folder/'input.mod';p.write_text(re.sub(r'(?m)^(\s*species\s+)22$',lambda m:m[1]+str(len(species)),p.read_text()))
        fresh.gzcopy(p,OUT/'inputs'/label/'input.mod.gz')
    if network is not None:
        (folder/'cno_extras.net').write_text(network)
        shutil.copy2(folder/'cno_extras.net',OUT/'inputs'/label/'cno_extras.net')
    if extra:
        p=folder/'inlist1';p.write_text(p.read_text().replace('&controls','&controls\n'+extra+'\n',1))
    p=folder/'profile_columns.list';names={line.split('!')[0].strip() for line in p.read_text().splitlines()}
    with p.open('a') as f:
        f.write('\n'+'\n'.join(k for k in ['energy','pressure','entropy','cv','cp','gamma1','chiRho','chiT','grada',*fresh.ISOS] if k not in names)+'\n')
    for name in ['inlist1','profile_columns.list']:
        shutil.copy2(folder/name,OUT/'inputs'/label/name)
    fresh.save(label+'-inputs.json',dict(classification='Proven',sha256={p.name:sha(p) for p in folder.iterdir() if p.is_file()}))

# ponytail: exact x86-64 SysV/gfortran ABI for the SHA-bound executable only.
# A different executable/compiler requires its own ABI validation, not a guessed layout.
class Registers(ct.Structure):
    _fields_=[(k,ct.c_ulonglong) for k in
        'r15 r14 r13 r12 rbp rbx r11 r10 r9 r8 rax rcx rdx rsi rdi orig_rax rip cs eflags rsp ss fs_base gs_base ds es fs gs'.split()]

def trace(label):
    context();folder=CACHE/label
    if (OUT/(label+'-species.json')).exists(): fresh.ISOS=json.loads((OUT/(label+'-species.json')).read_text())
    expected=json.loads((OUT/(label+'-inputs.json')).read_text())['sha256']
    for name,digest in expected.items(): assert sha(folder/name)==digest,name
    assert sha(folder/'binary')==json.loads((ROOT/'outputs/fresh-microphysics25/runtime.json').read_text())['sha256']
    symbols=subprocess.check_output(['nm',str(folder/'binary')],text=True)
    line=[s for s in symbols.splitlines() if s.endswith(' __net_lib_MOD_net_get')]
    assert len(line)==1;entry=int(line[0].split()[0],16)
    libc=ct.CDLL(None,use_errno=True);libc.ptrace.restype=ct.c_long
    libc.ptrace.argtypes=[ct.c_ulong,ct.c_ulong,ct.c_void_p,ct.c_void_p]
    def ptrace(request,pid,addr=0,data=0):
        ct.set_errno(0);result=libc.ptrace(request,pid,addr,data)
        if result==-1 and ct.get_errno(): raise OSError(ct.get_errno(),os.strerror(ct.get_errno()))
        return result
    def child(): ptrace(0,0)
    env=os.environ.copy();runtime=fresh.MESA.parent
    env.update(MESA_DIR=str(fresh.MESA),LD_LIBRARY_PATH=str(runtime/'mesasdk/lib')+':'+str(runtime/'mesasdk/lib64'),OMP_NUM_THREADS='1')
    assert not (folder/'execution.log').exists()
    begin=time.monotonic();records=[];mem=None;reg=Registers()
    with (folder/'execution.log').open('w') as log:
        proc=subprocess.Popen(['./binary'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=child)
        pid=proc.pid
        def wait():
            while time.monotonic()-begin<180:
                done,status=os.waitpid(pid,os.WNOHANG)
                if done:
                    assert os.WIFSTOPPED(status),(label,status,len(records))
                    if os.WSTOPSIG(status)==signal.SIGCHLD:
                        ptrace(7,pid,0,signal.SIGCHLD);continue
                    assert os.WSTOPSIG(status)==signal.SIGTRAP,(label,status)
                    return
                time.sleep(.0001)
            raise TimeoutError(label)
        def read(addr,n):
            result=os.pread(mem,n,addr);assert len(result)==n;return result
        def u64(addr): return struct.unpack('<Q',read(addr,8))[0]
        def scalar(addr): return struct.unpack('<d',read(addr,8))[0]
        def integer(addr): return struct.unpack('<i',read(addr,4))[0]
        def array(desc,shape):
            header=struct.unpack('<'+'q'*(3+3*len(shape)),read(desc,8*(3+3*len(shape))))
            assert header[3]==1,header
            for i,n in enumerate(shape):
                assert header[5+3*i]-header[4+3*i]+1==n,(header,shape)
                if i: assert header[3+3*i]==int(np.prod(shape[:i])),header
            return np.frombuffer(read(header[0],8*int(np.prod(shape))),dtype='<f8').copy().reshape(shape,order='F')
        def getregs(): ptrace(12,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def rewind(addr):
            reg.rip=addr;ptrace(13,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def breakpoint(addr):
            word=u64(addr);ptrace(4,pid,addr,(word & ~255)|0xcc);return word
        try:
            wait();mem=os.open('/proc/'+str(pid)+'/mem',os.O_RDONLY)
            original=breakpoint(entry);ptrace(7,pid)
            inp=np.load(OUT/(label+'-input.npz'));nzone=len(inp['dm']);ns=inp['X'].shape[1]
            while len(records)<nzone:
                wait();getregs();assert reg.rip==entry+1,hex(reg.rip)
                args=[reg.rdi,reg.rsi,reg.rdx,reg.rcx,reg.r8,reg.r9]+[u64(reg.rsp+8*(i-5)) for i in range(6,38)]
                assert integer(args[3])==ns and integer(args[1])==0
                row=dict(T=scalar(args[6]),rho=scalar(args[8]),X=array(args[5],(ns,)))
                ret=u64(reg.rsp)
                ptrace(4,pid,entry,original);rewind(entry);ptrace(9,pid);wait()
                retword=breakpoint(ret);ptrace(7,pid);wait();getregs();assert reg.rip==ret+1
                assert integer(args[37])==0
                row.update(heat=scalar(args[23]),heat_rho=scalar(args[24]),heat_T=scalar(args[25]),
                    heat_X=array(args[26],(ns,)),dxdt=array(args[27],(ns,)),
                    dxdt_rho=array(args[28],(ns,)),dxdt_T=array(args[29],(ns,)),
                    jacobian=array(args[30],(ns,ns)),neutrino=scalar(args[34]))
                assert all(np.all(np.isfinite(a)) for a in row.values())
                records.append(row)
                ptrace(4,pid,ret,retword);rewind(ret)
                if len(records)<nzone: breakpoint(entry)
                else: ptrace(17,pid);break
                ptrace(7,pid)
            complete=False
            while time.monotonic()-begin<180:
                path=folder/'LOGS1/profile1.data'
                if path.exists():
                    try:
                        h,d=mesa(path)
                        complete=int(h['model_number'])==1 and len(d['zone'])==nzone and d['zone'][-1]==nzone
                    except (ValueError,IndexError,KeyError,OSError): pass
                    if complete: break
                if proc.poll() is not None: break
                time.sleep(.05)
            assert complete,(label,'profile incomplete')
        finally:
            if mem is not None: os.close(mem)
            # This is our disposable child, never another session's process.
            if proc.poll() is None:
                os.kill(pid,signal.SIGKILL)
                try: os.waitpid(pid,0)
                except ChildProcessError: pass
    assert sha(folder/'binary')==expected['binary']
    np.savez_compressed(OUT/(label+'-native.npz'),**{k:np.array([r[k] for r in records]) for k in records[0]})
    shutil.copy2(folder/'execution.log',OUT/(label+'-execution.log'));fresh.collect(label)
    save(label+'-trace.json',dict(classification='Proven',entry_hex=hex(entry),calls=len(records),
        elapsed_s=time.monotonic()-begin,threads=1,executable_unchanged=True,
        numerical_memory_written=False,independent_profile_control_pending=True))
    print('TRACE',label,len(records),'native calls',flush=True)

def direct(): trace('direct-v2')

def retry_pilot():
    save('instrumentation-pilot.json',dict(classification='Proven',
        failure='The first tracer stopped on a normal SIGCHLD from a native initialization helper. '
        'Forward that signal unchanged; no scientific result was produced.'))
    shutil.copy2(CACHE/'direct/execution.log',OUT/'instrumentation-pilot.log')
    setup('direct-v2',dict(np.load(ROOT/'outputs/fresh-microphysics25/gr-input.npz')))

def analyze():
    d=np.load(OUT/'direct-v2-native.npz');state=np.load(OUT/'direct-v2-input.npz')
    _,profile=mesa(OUT/'direct-v2-profile.data.gz');_,old=mesa(OLD/'full-profile.data.gz')
    # Native loop order is checked, not inferred from the number of records.
    control=dict(logT=float(max(abs(np.log(d['T'])-state['lnT']))),
        logrho=float(max(abs(np.log(d['rho'])-state['lnd']))),
        X=float(np.max(abs(d['X']-state['X']))))
    assert max(control.values())<1e-12,control
    profile_error={k:float(max(abs(profile[k]-old[k])/np.maximum(1,abs(old[k])))) for k in
        ['rho','logT','eps_nuc','eps_nuc_neu_total','energy','pressure','opacity']}
    direct_error={k:float(max(abs(d[key]-profile[k])/np.maximum(1,abs(profile[k])))) for key,k in
        [('heat','eps_nuc'),('neutrino','eps_nuc_neu_total')]}
    assert max([*profile_error.values(),*direct_error.values()])<1e-10,(profile_error,direct_error)
    rebuilt=np.load(OLD/'reconstructed-composition-source.npz')
    scale=np.maximum(1e-30,abs(d['dxdt']).max(axis=1));err=abs(d['dxdt']-rebuilt['dX_dproper_s'])
    relative=err/scale[:,None];worst=np.unravel_index(relative.argmax(),relative.shape)
    baryon=float(max(abs(d['dxdt'].sum(axis=1))/np.maximum(1e-30,abs(d['dxdt']).sum(axis=1))))
    jac_baryon=float(np.max(abs(d['jacobian'].sum(axis=1))/np.maximum(1e-30,abs(d['jacobian']).sum(axis=1))))
    fluor=[fresh.ISOS.index(k) for k in ['f17','f18','f19']]
    fluor_generation=d['dxdt'][:,fluor].sum(axis=1)
    result=dict(classification='Counterexample candidate',zones=len(scale),species=22,
        trace_input_errors=control,untraced_profile_errors=profile_error,direct_output_errors=direct_error,
        vector_max_relative_to_zone_max=float(relative.max()),worst_zone=int(worst[0]+1),worst_species=fresh.ISOS[worst[1]],
        vector_max_absolute=float(err.max()),baryon_relative=baryon,jacobian_baryon_relative=jac_baryon,
        max_F_generation_per_s=float(fluor_generation.max()),positive_F_generation_cells=int(sum(fluor_generation>0)),
        independently_validated=bool(relative.max()<1e-8 and baryon<1e-12),
        meaning='Direct dxdt returned by the unchanged native net_get, independently compared with the previous channel-energy reconstruction. '
        'The native Jacobian is exported but not yet certified; time evolution, EOS consistency and GR structure remain separate gates.')
    save('direct-vector-validation.json',result);print(json.dumps(result,indent=2),flush=True)

def reaction_plan():
    assert json.loads((OUT/'direct-vector-validation.json').read_text())['independently_validated']
    save('reaction-continuation-plan.json',dict(classification='Counterexample candidate',
        duration_proper_s=1.6294*86400,step_counts=[1,2,4,8],
        method='Backward Euler with native dxdt and native approximate Jacobian, reevaluating '
        'the complete native EOS/network/opacity at every Newton iterate. Keep each zone rho,T fixed. '
        'Initial trace result provides the first guess; use fresh result for acceptance.',
        residual_tolerance=1e-7,abundance_floor_for_residual=1e-24,refinement_tolerance=1e-3,
        gates='Nonnegative 22-species fractions; baryon sum error <1e-12; residual per species '
        'normalized to max(current change, dt*initial source, 1e-24) <1e-7. '
        'Compare successive end states relative to max(total change,1e-20).',
        interpretation='Externally thermostated, fixed-density reaction control. All 22 species '
        'including generated F enter the same native EOS and reaction evaluator. This is not '
        'GR thermal evolution or proof of common EOS thermodynamic consistency.'))

def reactions(arithmetic=False):
    plan_name='reaction-refinement-plan.json' if arithmetic==3 else 'reaction-absolute-plan.json' if arithmetic==2 else 'reaction-arithmetic-plan.json' if arithmetic else 'reaction-continuation-plan.json'
    plan=json.loads((OUT/plan_name).read_text());base=dict(np.load(OUT/'direct-v2-input.npz'))
    first=dict(np.load(OUT/'direct-v2-native.npz'));runs=[];previous=None
    if arithmetic==3: previous=np.load(OUT/'burn-8-end.npz')['X']
    for steps in plan['step_counts']:
        X=base['X'].copy();current=first;dt=plan['duration_proper_s']/steps;records=[]
        for step in range(steps):
            initial=X.copy();initial_f=current['dxdt'].copy()
            for iteration in range(12):
                residual=X-initial-dt*current['dxdt']
                scale=np.maximum.reduce([abs(X-initial),dt*abs(initial_f),np.full_like(X,1e-24)])
                error=float(np.max(abs(residual)/scale))
                budget=plan['residual_tolerance']*scale
                if arithmetic: budget+=64*np.finfo(float).eps*np.maximum(abs(X),abs(initial))
                if arithmetic>=2: budget+=1e-24
                score=float(np.max(abs(residual)/budget))
                if iteration>0 and score<1: break
                delta=np.linalg.solve(np.eye(22)[None,:,:]-dt*current['jacobian'],-residual[...,None])[...,0]
                # Mass-conserving roundoff correction, with its size recorded.
                largest=np.argmax(X,axis=1);zones=np.arange(len(X))
                mass_correction=-delta.sum(axis=1);delta[zones,largest]+=mass_correction
                trial=X+delta
                # Only arithmetic-scale negatives can be projected; record every correction.
                negative=float(np.min(trial));assert negative>-1e-22,('negative abundance',steps,step,iteration,negative)
                projection=np.maximum(-trial,0);trial+=projection;trial[zones,largest]-=projection.sum(axis=1)
                assert trial.min()>=0 and max(abs(trial.sum(axis=1)-1))<1e-12
                label=f'heatbath-{steps}-{step}-{iteration}' if arithmetic else f'burn-{steps}-{step}-{iteration}'
                if arithmetic>=2: label=f'controlled-{steps}-{step}-{iteration}'
                data={**base,'X':trial}
                if (OUT/(label+'-native.npz')).exists():
                    cached=np.load(OUT/(label+'-input.npz'))
                    assert all(np.array_equal(cached[k],v) for k,v in data.items()),label
                    assert json.loads((OUT/(label+'-import-control.json')).read_text())['passed']
                else:
                    if (OUT/(label+'-inputs.json')).exists():
                        partial=label;label+='-resume'
                        log=CACHE/partial/'execution.log'
                        if log.exists(): shutil.copy2(log,OUT/(partial+'-interrupted.log'))
                        save(partial+'-interrupted.json',dict(classification='Proven',
                            reason='User interruption before a complete trace artifact; restart from unchanged input in a new folder.',
                            resumed_label=label))
                    setup(label,data);trace(label)
                X=trial;current=dict(np.load(OUT/(label+'-native.npz')))
                # Validate every traced state's exact matching independently of solver residual.
                assert np.max(abs(current['X']-X))<1e-12
                assert max(abs(np.log(current['T'])-base['lnT']))<1e-12
                assert max(abs(np.log(current['rho'])-base['lnd']))<1e-12
                records.append(dict(step=step,iteration=iteration,label=label,pre_update_residual=error,pre_update_budget_score=score,
                    mass_correction_max=float(max(abs(mass_correction))),negative_projection_max=float(projection.max())))
            else: raise AssertionError(('Newton not converged',steps,step,error))
            records[-1]['accepted_residual']=error;records[-1]['accepted_budget_score']=score
            print('BURN ACCEPT',steps,step,iteration,error,flush=True)
        _,final=mesa(OUT/(records[-1]['label']+'-profile.data.gz'))
        fluor=[fresh.ISOS.index(k) for k in ['f17','f18','f19']]
        denominator=np.maximum.reduce([abs(X-base['X']),128*np.finfo(float).eps*abs(X),np.full_like(X,1e-20)])
        final_difference=None if previous is None else float(np.max(abs(X-previous)/denominator))
        change=X-base['X'];np.savez_compressed(OUT/f'burn-{steps}-end.npz',X=X,dxdt=current['dxdt'])
        row=dict(steps=steps,records=records,end_label=records[-1]['label'],
            max_absolute_change=float(abs(change).max()),max_F_fraction=float(X[:,fluor].sum(axis=1).max()),
            baryon_max_absolute=float(max(abs(X.sum(axis=1)-1))),
            refinement_max_relative_change=final_difference,
            max_relative_heat_change=float(max(abs(final['eps_nuc']-first['heat'])/np.maximum(1,abs(first['heat'])))))
        output_name='reaction-controlled-16.json' if arithmetic==3 else 'reaction-controlled.json' if arithmetic==2 else 'reaction-continuation.json'
        runs.append(row);save(output_name,dict(classification='Counterexample candidate',runs=runs,
            completed_all_refinements=len(runs)==len(plan['step_counts']),
            full_GR_thermal_evolution=False,common_EOS_thermodynamics_certified=False))
        previous=X.copy();print('BURN END',steps,{k:v for k,v in row.items() if k!='records'},flush=True)

def arithmetic_plan():
    original=json.loads((OUT/'reaction-continuation-plan.json').read_text())
    original.update(original_plan_sha256=sha(OUT/'reaction-continuation-plan.json'),
        residual_budget='1e-7*max(abs(X-X0),dt*abs(f0),1e-24) + 64*eps*max(abs(X),abs(X0))',
        refinement_denominator='max(abs(total change),128*eps*abs(X),1e-20)',
        scope='Separate arithmetic-aware solve; preserve the strict relative-residual failure. '
        'The rounding term covers storage/subtraction scale only, not EOS or physical-model error.')
    save('reaction-arithmetic-plan.json',original)
    base=np.load(OUT/'direct-v2-input.npz');first=np.load(OUT/'direct-v2-native.npz');rows=[]
    for p in sorted(OUT.glob('burn-1-0-*-native.npz')):
        d=np.load(p);dt=original['duration_proper_s'];r=d['X']-base['X']-dt*d['dxdt']
        scale=np.maximum.reduce([abs(d['X']-base['X']),dt*abs(first['dxdt']),np.full_like(r,1e-24)])
        i=np.unravel_index((abs(r)/scale).argmax(),r.shape)
        rows.append(dict(file=p.name,relative_residual=float((abs(r)/scale).max()),worst_zone=int(i[0]+1),
            worst_species=fresh.ISOS[i[1]],absolute_residual_at_worst=float(r[i]),abundance_at_worst=float(d['X'][i])))
    save('strict-residual-failure.json',dict(classification='Counterexample candidate',passed=False,runs=rows,
        cause='A relative change residual can demand absolute accuracy far below the spacing of a stored nonzero abundance. '
        'The later arithmetic-aware acceptance is a separately defined gate; this original gate remains failed.'))

def heatbath(): reactions(arithmetic=True)

def absolute_plan():
    previous=json.loads((OUT/'reaction-arithmetic-plan.json').read_text())
    previous.update(step_counts=[8],previous_plan_sha256=sha(OUT/'reaction-arithmetic-plan.json'),
        absolute_residual_tolerance=1e-24,
        correction='Add an explicit absolute reaction-balance tolerance of 1e-24 in mass fraction. '
        'The prior arithmetic-only residual failed near F17 balance with an absolute residual of about 5e-30. '
        'Preserve that gate and all its iterates; do not describe this tolerance as a physical error certificate.')
    save('reaction-absolute-plan.json',previous)
    save('arithmetic-residual-failure.json',dict(classification='Counterexample candidate',passed=False,
        step_count=8,failed_step_zero_based=4,iterations=12,
        witness='heatbath-8-4-11-native.npz',worst_zone=2255,worst_species='f17',
        residual_about=-4.95057050648201e-30,budget_about=1.7864769868129497e-30,
        conclusion='A separate absolute reaction-balance tolerance is necessary near cancellation.'))

def controlled(): reactions(arithmetic=2)

def refinement_plan():
    p=json.loads((OUT/'reaction-absolute-plan.json').read_text());p.update(step_counts=[16],
        previous_result_sha256=sha(OUT/'reaction-controlled.json'),
        additional_mixed_refinement_budget='1e-16 + 1e-3*abs(total change) per isotope and zone',
        original_relative_change_gate_preserved=True)
    save('reaction-refinement-plan.json',p)

def controlled16(): reactions(arithmetic=3)

def eos_plan():
    save('EOS-continuation-plan.json',dict(classification='Counterexample candidate',
        steps=[.02,.01,.005],local_derivative_tolerance=1e-3,
        purpose='Larger symmetric stencils test whether native single-precision table outputs '
        'permit a usable local derivative check. The previous small-step failures remain frozen. '
        'No finite-difference agreement is a global error certificate.',
        tests=['pressure derivatives versus chiT and chiRho','dE/dlnT versus T*cv',
               'T*dS/dlnT versus T*cv','dE/dlnrho - T*dS/dlnrho - P/rho'],
        native_HELM_density_note='Source helm_electron_positron.dek deliberately uses an independently '
        'interpolated dP/drho instead of the derivative of its interpolated free-energy pressure. '
        'Inspect the source convention separately from the former numerical discrepancy.'))

def eos_runs():
    context();plan=json.loads((OUT/'EOS-continuation-plan.json').read_text())
    base=dict(np.load(OUT/'direct-v2-input.npz'))
    for h in plan['steps']:
        for var in ['lnT','lnd']:
            for sign in [-1,1]:
                data={**base,var:base[var]+sign*h};label=f'EOS-{h}-{var}-{sign}'
                setup(label,data);fresh.run(label);fresh.collect(label)
    eos_analysis()

def eos_analysis():
    plan=json.loads((OUT/'EOS-continuation-plan.json').read_text())
    _,d=mesa(OUT/'direct-v2-profile.data.gz');T=10**d['logT'];rho=d['rho'];P=d['pressure']
    unit=1.3806504e-16*6.02214179e23;rows=[];previous=None
    for h in plan['steps']:
        fd={}
        for var in ['lnT','lnd']:
            _,p=mesa(OUT/f'EOS-{h}-{var}-1-profile.data.gz');_,m=mesa(OUT/f'EOS-{h}-{var}--1-profile.data.gz')
            fd[var]={k:(p[k]-m[k])/(2*h) for k in ['pressure','energy','entropy']}
        pairs={'chiT':(fd['lnT']['pressure']/P,d['chiT']),
            'chiRho':(fd['lnd']['pressure']/P,d['chiRho']),
            'cv':(fd['lnT']['energy']/T,d['cv']),
            'entropy_T_cv':(unit*fd['lnT']['entropy'],d['cv'])}
        errors={k:float(max(abs(a-b)/np.maximum(1,abs(b)))) for k,(a,b) in pairs.items()}
        errors['first_law_rho']=float(max(abs(fd['lnd']['energy']-T*unit*fd['lnd']['entropy']-P/rho)/(P/rho)))
        drift=None if previous is None else {v+'_'+k:float(max(abs(fd[v][k]-previous[v][k])/
            np.maximum(abs(previous[v][k]), {'pressure':P,'energy':T*d['cv'],'entropy':d['cv']/unit}[k])))
            for v in fd for k in fd[v]}
        rows.append(dict(step=h,errors=errors,halving_drift=drift,passed=max(errors.values())<plan['local_derivative_tolerance']))
        previous=fd
    save('EOS-continuation.json',dict(classification='Counterexample candidate',rows=rows,
        common_EOS_certified=False,original_failures_promoted=False,
        interpretation='Bigger stencils reduce roundoff amplification but mix interpolation truncation, '
        'branch transitions and differences between tabulated values and reported derivatives. '
        'This diagnostic alone cannot turn the common native evaluator into a certified thermodynamic EOS.'))
    print('EOS',json.dumps(rows),flush=True)

def boundary_plan():
    save('network-boundary-plan.json',dict(classification='Counterexample candidate',
        intervention='At original peak-heating zone 2592 only, transfer its entire H1 fraction to He4. '
        'Keep rho,T and all other abundances unchanged. Evaluate the original 22-species network.',
        positivity_gate='At X_H1=0, the returned dX_H1/dt must be nonnegative for an ODE invariant on the abundance simplex.',
        reason='The effective PP-II/III stoichiometry consumes H1 although its initial He3+He4 rate omits H1. '
        'This boundary probe checks the validity domain of the eliminated-intermediate network; it is not the actual GR initial state.'))
    base=dict(np.load(OUT/'direct-v2-input.npz'));base['X']=base['X'].copy();i=2591
    base['X'][i,2]+=base['X'][i,0];base['X'][i,0]=0
    setup('H-zero-boundary',base)

def boundary():
    trace('H-zero-boundary');d=np.load(OUT/'H-zero-boundary-native.npz');i=2591
    save('network-boundary.json',dict(classification='Counterexample candidate',zone=i+1,
        rho=float(d['rho'][i]),T=float(d['T'][i]),H_fraction=float(d['X'][i,0]),H_derivative=float(d['dxdt'][i,0]),
        simplex_inward=bool(d['dxdt'][i,0]>=0),
        meaning='A negative derivative at zero abundance is a counterexample to global simplex invariance '
        'of this reduced network, not evidence that the actual initial-star path has reached this boundary.'))
    print('H BOUNDARY',d['X'][i,0],d['dxdt'][i,0],flush=True)

def quantization_plan():
    base=np.load(OUT/'direct-v2-input.npz');i=2336
    log10T=np.float32(base['lnT'][i]/np.log(10));neighbor=np.nextafter(log10T,np.float32(np.inf))
    boundary=(float(log10T)+float(neighbor))/2*np.log(10)
    save('EOS-quantization-plan.json',dict(classification='Counterexample candidate',zone=i+1,
        logarithmic_boundary=boundary,half_widths=[1e-7,1e-9,1e-11],
        purpose='Cross a known REAL(log10T) rounding boundary in the default EOS source. '
        'Check whether the finite output jump shrinks with the bracketing interval. '
        'This concerns the evaluator and is not proof that the continuum physical EOS is discontinuous.'))

def quantization():
    plan=json.loads((OUT/'EOS-quantization-plan.json').read_text());i=plan['zone']-1
    base=dict(np.load(OUT/'direct-v2-input.npz'));rows=[]
    for h in plan['half_widths']:
        pair=[]
        for sign in [-1,1]:
            data={**base,'lnT':base['lnT'].copy()};data['lnT'][i]=plan['logarithmic_boundary']+sign*h
            label=f'rounding-{h}-{sign}';setup(label,data);fresh.run(label);fresh.collect(label)
            _,d=mesa(OUT/(label+'-profile.data.gz'));pair.append(d)
        m,p=pair
        rows.append(dict(half_width=h,actual_lnT_separation=float((p['logT'][i]-m['logT'][i])*np.log(10)),
            relative_jumps={k:float((p[k][i]-m[k][i])/max(abs(m[k][i]),1e-30)) for k in ['pressure','energy','entropy']}))
    save('EOS-quantization.json',dict(classification='Counterexample candidate',zone=i+1,rows=rows,
        continuum_physical_EOS_discontinuity_claimed=False));print('ROUNDING',rows,flush=True)

def explicit_plan():
    extras=['h2','li7','be7','b8'];species=['h1','h2','he3','he4','li7','be7','b8']+ISOS[3:]
    network=(ROOT/'outputs/thermal-closure22/inputs/cno_extras.net').read_text()+"\ninclude 'add_pp_extras'\n"
    save('explicit-PP-plan.json',dict(classification='Counterexample candidate',species=species,extras=extras,
        method='Use the existing native add_pp_extras network definition: resolve H2, Li7, Be7 and B8, '
        'remove their effective PP replacement reactions. Initialize the added inventories to zero explicitly. '
        'Compare the H=0 boundary and extract the conditional intermediate-state Jacobian.',
        gates='All input rho,T,X retained exactly. Test H boundary inwardness, baryon source/Jacobian sums. '
        'Any relaxation spectrum is local at frozen bulk composition, density, temperature and EOS auxiliaries; '
        'it is not an observable scalar pole or an independently initialized stellar model.'))
    for label,origin in [('explicit-PP','direct-v2'),('explicit-PP-Hzero','H-zero-boundary')]:
        base=dict(np.load(OUT/(origin+'-input.npz')));old=base['X'];base['X']=np.zeros((len(old),len(species)))
        for j,name in enumerate(species):
            if name in ISOS: base['X'][:,j]=old[:,ISOS.index(name)]
        setup(label,base,species=species,network=network)

def explicit_pp():
    for label in ['explicit-PP','explicit-PP-Hzero']: trace(label)
    explicit_analysis()

def explicit_analysis():
    plan=json.loads((OUT/'explicit-PP-plan.json').read_text());species=plan['species']
    d=np.load(OUT/'explicit-PP-native.npz');zero=np.load(OUT/'explicit-PP-Hzero-native.npz')
    extra=[species.index(k) for k in plan['extras']];A=d['jacobian'][:,extra,:][:,:,extra]
    b=d['dxdt'][:,extra];eigen=np.linalg.eigvals(A)
    times=np.divide(-1,eigen.real,out=np.full_like(eigen.real,np.inf),where=eigen.real<0)
    # A is the chemostatted intermediate subsystem; do not label it the full stellar Jacobian.
    stable=np.all(eigen.real<0,axis=1);i=2591
    equilibrium=np.linalg.solve(A[stable],-b[stable,:,None])[...,0]
    rates=np.load(OUT/'direct-v2-native.npz')
    np.savez_compressed(OUT/'explicit-PP-local-subsystem.npz',A=A,b=b,eigenvalues=eigen,
        times_s=times,stable=stable,conditional_equilibrium=equilibrium)
    result=dict(classification='Counterexample candidate',species=len(species),zones=len(A),
        peak_zone=i+1,peak_subsystem_times_s=sorted(times[i].tolist()),
        peak_zero_H_derivative=float(zero['dxdt'][i,species.index('h1')]),
        zero_H_boundary_inward=bool(zero['dxdt'][i,species.index('h1')]>=0),
        stable_intermediate_subsystem_zones=int(stable.sum()),
        intermediate_equilibrium_max=float(np.max(equilibrium)),
        max_baryon_source_relative=float(np.max(abs(d['dxdt'].sum(axis=1))/np.maximum(1e-30,abs(d['dxdt']).sum(axis=1)))),
        zero_intermediate_heat_max_normalized_difference=float(np.max(abs(d['heat']-rates['heat'])/np.maximum(1,abs(rates['heat'])))),
        state_domain='Zero added intermediate abundances; fixed bulk chemostat. The equilibrium solve is a local affine diagnostic, not a conserved-composition stellar initialization.',
        actual_scalar_observable_pole=False,full_stellar_evolution=False)
    save('explicit-PP.json',result);print('EXPLICIT PP',json.dumps(result),flush=True)

def pp_response_plan():
    save('PP-response-plan.json',dict(classification='Counterexample candidate',
        domain='Original 896 zones with abs(eps_nuc)>1 erg/g/s and stable added-intermediate block.',
        initialization='Solve four intermediate source equations with their added baryon mass removed from He4. '
        'This is a declared changed-composition state; no GR rematching is asserted.',
        temperature_drive='External delta ln T at fixed baryon density and held bulk abundances. '
        'Intermediate changes exchange mass with the He4 reservoir. Heat is the output.',
        periods_days=[1.6294,327.26],derivative_steps=[1e-4,5e-5],
        derivative_tolerance=1e-5,steady_relative_tolerance=1e-6,
        limits='A conditional nuclear heat-transfer pole does not establish a scalar charge, '
        'a free-fall force, or SEP violation. At the GR zero-scalar branch linear thermal-to-scalar coupling remains zero.'))

def pp_setup(label,data):
    plan=json.loads((OUT/'explicit-PP-plan.json').read_text())
    net=(OUT/'inputs/explicit-PP/cno_extras.net').read_text()
    setup(label,data,species=plan['species'],network=net)

def pp_steady():
    plan=json.loads((OUT/'explicit-PP-plan.json').read_text());species=plan['species']
    extra=[species.index(k) for k in plan['extras']];he4=species.index('he4')
    data=dict(np.load(OUT/'explicit-PP-input.npz'));initial=data['X'].copy()
    d=dict(np.load(OUT/'explicit-PP-native.npz'));sub=np.load(OUT/'explicit-PP-local-subsystem.npz')
    old=np.load(OUT/'direct-v2-native.npz');domain=(abs(old['heat'])>1)&sub['stable'];indices=np.flatnonzero(domain)
    assert len(indices)==896
    rows=[]
    for iteration in range(3):
        A=d['jacobian'][:,extra,:][:,:,extra]-d['jacobian'][:,extra,he4][:,:,None]
        delta=np.linalg.solve(A[domain],-d['dxdt'][domain][:,extra,None])[...,0]
        X=data['X'].copy();X[np.ix_(indices,extra)]+=delta;X[domain,he4]-=delta.sum(axis=1)
        assert X.min()>=0 and max(abs(X.sum(axis=1)-1))<1e-12
        data={**data,'X':X};label=f'PP-steady-{iteration}';pp_setup(label,data);trace(label)
        d=dict(np.load(OUT/(label+'-native.npz')))
        scale=np.maximum(1e-30,np.sum(abs(A[domain]*X[np.ix_(indices,extra)][:,None,:]),axis=2))
        residual=float(np.max(abs(d['dxdt'][domain][:,extra])/scale))
        rows.append(dict(label=label,relative_residual=residual,max_update=float(abs(delta).max())))
        print('PP STEADY',iteration,residual,flush=True)
    np.savez_compressed(OUT/'PP-steady-state.npz',**data,domain=domain)
    save('PP-steady.json',dict(classification='Counterexample candidate',rows=rows,domain_cells=int(domain.sum()),
        baryon_max_absolute=float(max(abs(data['X'].sum(axis=1)-1))),
        max_composition_change=float(abs(data['X']-initial).max()),passed=residual<1e-6,
        GR_structure_rematched=False))

def pp_derivatives():
    data=dict(np.load(OUT/'PP-steady-state.npz'));data.pop('domain')
    plan=json.loads((OUT/'PP-response-plan.json').read_text())
    for h in plan['derivative_steps']:
        for sign in [-1,1]:
            label=f'PP-T-{h}-{sign}';pp_setup(label,{**data,'lnT':data['lnT']+sign*h});trace(label)
    pp_response()

def pp_response():
    plan=json.loads((OUT/'PP-response-plan.json').read_text());species=json.loads((OUT/'explicit-PP-plan.json').read_text())['species']
    extra=[species.index(k) for k in ['h2','li7','be7','b8']];he4=species.index('he4')
    state=np.load(OUT/'PP-steady-state.npz');domain=state['domain'];d=np.load(OUT/'PP-steady-2-native.npz')
    A=d['jacobian'][:,extra,:][:,:,extra]-d['jacobian'][:,extra,he4][:,:,None]
    C=d['heat_X'][:,extra]-d['heat_X'][:,he4,None]
    b=d['dxdt_T'][:,extra]*d['T'][:,None];D=d['heat_T']*d['T'];checks=[];previous=None
    for h in plan['derivative_steps']:
        p=np.load(OUT/f'PP-T-{h}-1-native.npz');m=np.load(OUT/f'PP-T-{h}--1-native.npz')
        fd=(p['dxdt'][:,extra]-m['dxdt'][:,extra])/(2*h);heat=(p['heat']-m['heat'])/(2*h)
        scale=np.maximum(1e-30,abs(b).max(axis=1));err=np.max(abs(fd-b)/scale[:,None],axis=1)
        checks.append(dict(step=h,source_derivative_error=float(err[domain].max()),
            heat_derivative_error=float(np.max((abs(heat-D)/np.maximum(1,abs(D)))[domain])),
            source_step_halving=None if previous is None else float(np.max((abs(fd-previous)/scale[:,None])[domain]))))
        previous=fd
    eigen=np.linalg.eigvals(A[domain]);i=2591;results=[]
    metric=np.load(ROOT/'outputs/fresh-microphysics25/gr-input.npz');N=np.exp(metric['nu'])
    weight=metric['dm']*N*N/fresh.LSUN
    for period in plan['periods_days']:
        omega=2*np.pi/(period*86400)
        response=D[domain]+np.einsum('ij,ij->i',C[domain],np.linalg.solve(1j*(omega/N[domain])[:,None,None]*np.eye(4)[None,:,:]-A[domain],b[domain,:,None])[...,0])
        zero=D[domain]+np.einsum('ij,ij->i',C[domain],np.linalg.solve(-A[domain],b[domain,:,None])[...,0])
        peak=int(np.flatnonzero(np.flatnonzero(domain)==i)[0]);total=weight[domain]@response;static=weight[domain]@zero
        results.append(dict(period_days=period,peak_heat_response_real=float(response[peak].real),
            peak_heat_response_imag=float(response[peak].imag),peak_static_response=float(zero[peak]),
            domain_weighted_response_real_Lsun=float(total.real),domain_weighted_response_imag_Lsun=float(total.imag),
            domain_weighted_static_response_Lsun=float(static)))
    np.savez_compressed(OUT/'PP-response-matrices.npz',A=A,b=b,C=C,D=D,domain=domain,eigenvalues=eigen)
    save('PP-response.json',dict(classification='Counterexample candidate',checks=checks,results=results,
        peak_intermediate_times_s=sorted((-1/np.linalg.eigvals(A[i]).real).tolist()),
        native_derivative_check_passed=all(max(c['source_derivative_error'],c['heat_derivative_error'])<plan['derivative_tolerance'] for c in checks),
        input='Externally imposed delta ln T harmonic in asymptotic coordinate time; omega_proper=omega_infinity/N. Fixed density, bulk reservoir and metric weights.',
        output='Local nuclear heat. The weighted sum is a source diagnostic, not emergent surface luminosity.',
        actual_scalar_charge_or_force_transfer=False,physical_stellar_initialization_certified=False))
    print('PP RESPONSE',checks,results,flush=True)

def pp_reduction():
    a=np.load(OUT/'PP-response-matrices.npz');i=2591
    eigen,V=np.linalg.eig(a['A'][i]);assert np.max(abs(eigen.imag))<1e-12 and np.all(eigen.real<0)
    residue=(a['C'][i]@V)*np.linalg.solve(V,a['b'][i]);assert max(abs(residue.imag))<1e-12
    tau=-1/eigen.real;amplitude=(residue*tau).real;slow=int(tau.argmax());fast=np.arange(4)!=slow
    N=float(np.exp(np.load(ROOT/'outputs/fresh-microphysics25/gr-input.npz')['nu'][i]))
    Dfast=float(a['D'][i]+amplitude[fast].sum());rows=[]
    for period in [1.6294,327.26]:
        w=2*np.pi/(period*86400*N);exact=a['D'][i]+sum(amplitude/(1+1j*w*tau))
        reduced=Dfast+amplitude[slow]/(1+1j*w*tau[slow])
        bound=float(sum(abs(amplitude[fast])*w*tau[fast]))
        assert abs(exact-reduced)<=bound*(1+1e-6)+1e-9
        rows.append(dict(period_days=period,absolute_reduction_error=float(abs(exact-reduced)),
            error_bound=bound,relative_error=float(abs(exact-reduced)/abs(exact))))
    save('PP-one-state-reduction.json',dict(classification='Counterexample candidate',zone=i+1,
        tau_proper_days=float(tau[slow]/86400),tau_coordinate_days=float(tau[slow]/(86400*N)),lapse=N,
        tau_fast_proper_s=sorted(tau[fast].tolist()),
        direct_plus_fast_heat_per_dlnT=Dfast,slow_heat_amplitude_per_dlnT=float(amplitude[slow]),
        all_static_modal_amplitudes=amplitude.tolist(),rows=rows,
        equation='tau_coordinate*dz/dt_infinity + z = delta ln T; delta q_nuc = D_fast*delta ln T + a_slow*z',
        coefficient_origin='Local native reaction Jacobian at the separately initialized chemostatted intermediate equilibrium.',
        mass_sensitivity_or_scalar_coupling_derived=False,
        conclusion='A concrete one-state nuclear heat response emerges in this conditional reduction. '
        'It supplies an internal-time candidate, not the missing map from orbital/scalar drive to mass/force readout.'))
    print('ONE STATE',tau[slow]/86400,amplitude[slow],rows,flush=True)

def time_convention():
    for name in ['PP-response.json','PP-one-state-reduction.json']:
        dest=OUT/name.replace('.json','-proper-period-pilot.json')
        assert not dest.exists();shutil.copy2(OUT/name,dest)
    save('PP-time-convention.json',dict(classification='Proven',
        clarification='The original pilot used the same proper-period in every zone. '
        'For a common asymptotic time harmonic, use omega_proper=omega_infinity/N in every local block. '
        'Keep both pilot outputs; the final weighted diagnostic uses the latter convention.'))
    pp_response();pp_reduction()

def pp_composition_plan():
    save('PP-composition-derivative-plan.json',dict(classification='Counterexample candidate',zone=2592,
        fractional_step=1e-4,relative_tolerance=1e-5,
        method='At the steady peak zone, perturb each positive intermediate abundance by +/-1e-4 of itself, '
        'compensate He4, and reevaluate native EOS/rates. Check both the constrained intermediate Jacobian '
        'and the heat readout C. No untested native composition derivative is promoted to a full certificate.'))

def pp_composition():
    data=dict(np.load(OUT/'PP-steady-state.npz'));data.pop('domain');plan=json.loads((OUT/'PP-composition-derivative-plan.json').read_text())
    species=json.loads((OUT/'explicit-PP-plan.json').read_text())['species'];i=plan['zone']-1;he4=species.index('he4');rows=[]
    extra=[species.index(k) for k in ['h2','li7','be7','b8']];mats=np.load(OUT/'PP-response-matrices.npz')
    for j,name in enumerate(['h2','li7','be7','b8']):
        ix=species.index(name);h=data['X'][i,ix]*plan['fractional_step'];assert h>0;pair=[]
        for sign in [-1,1]:
            X=data['X'].copy();X[i,ix]+=sign*h;X[i,he4]-=sign*h;label=f'PP-X-{name}-{sign}'
            pp_setup(label,{**data,'X':X});trace(label);pair.append(np.load(OUT/(label+'-native.npz')))
        m,p=pair;fd=(p['dxdt'][i,extra]-m['dxdt'][i,extra])/(2*h)
        Cfd=(p['heat'][i]-m['heat'][i])/(2*h)
        Aerror=float(np.max(abs(fd-mats['A'][i,:,j]))/max(1e-30,np.max(abs(mats['A'][i,:,j]))))
        Cerror=float(abs(Cfd-mats['C'][i,j])/max(1e-30,abs(mats['C'][i,j])))
        rows.append(dict(species=name,step=h,source_column_error=Aerror,heat_readout_error=Cerror,
            passed=max(Aerror,Cerror)<plan['relative_tolerance']))
    save('PP-composition-derivatives.json',dict(classification='Counterexample candidate',zone=i+1,rows=rows,
        passed=all(r['passed'] for r in rows),global_composition_derivatives_certified=False))
    print('PP COMPOSITION',rows,flush=True)

def composition_summary():
    base=np.load(OUT/'direct-v2-input.npz')['X'];rows=[];previous=None
    for steps in [1,2,4,8,16]:
        X=np.load(OUT/f'burn-{steps}-end.npz')['X'];change=X-base
        row=dict(steps=steps,max_absolute_change=float(abs(change).max()))
        if previous is not None:
            diff=abs(X-previous);raw=diff/np.maximum.reduce([abs(change),128*np.finfo(float).eps*abs(X),np.full_like(X,1e-20)])
            mixed=diff/(1e-16+1e-3*abs(change))
            worst=np.unravel_index(mixed.argmax(),mixed.shape)
            row.update(max_absolute_refinement=float(diff.max()),relative_to_global_max_change=float(diff.max()/abs(change).max()),
                original_component_relative_max=float(raw.max()),original_relative_gate_passed=bool(raw.max()<1e-3),
                additional_mixed_budget_score=float(mixed.max()),additional_mixed_gate_passed=bool(mixed.max()<1),
                mixed_worst_zone=int(worst[0]+1),mixed_worst_species=ISOS[worst[1]])
        rows.append(row);previous=X
    old=np.load(OUT/'direct-v2-native.npz');label=json.loads((OUT/'reaction-controlled-16.json').read_text())['runs'][0]['end_label']
    _,final=mesa(OUT/(label+'-profile.data.gz'));_,initial=mesa(OUT/'direct-v2-profile.data.gz')
    state=np.load(OUT/'direct-v2-input.npz');weights=state['dm']*np.exp(2*state['nu'])/fresh.LSUN
    changes={k:float(max(abs(final[k]-initial[k])/np.maximum(1e-30,abs(initial[k])))) for k in ['energy','pressure','opacity','cv']}
    save('composition-continuation-summary.json',dict(classification='Counterexample candidate',rows=rows,
        original_strict_Newton_gate=False,original_arithmetic_only_Newton_gate=False,
        final_mixed_time_refinement_passed=rows[-1]['additional_mixed_gate_passed'],
        final_native_EOS_fractional_changes=changes,
        initial_redshift_weighted_heat_Lsun=float(weights@old['heat']),
        final_redshift_weighted_heat_Lsun=float(weights@final['eps_nuc']),
        meaning='Fixed-rho,T 22-species reaction evolution, with all generated species passed through native EOS/rates/opacity. '
        'Source powers use fixed prior GR weights and are not emergent luminosities or a GR evolution solution.'))
    print('COMPOSITION SUMMARY',rows,changes,flush=True)

def exponential_plan():
    save('reaction-exponential-plan.json',dict(classification='Counterexample candidate',steps=[1,2,4],
        prior_BE_failure_sha256=sha(OUT/'composition-continuation-summary.json'),duration_proper_s=1.6294*86400,
        method='Exponential Rosenbrock-Euler: delta=dt*phi1(dt*J)*f, calculated by the installed scipy.linalg.expm '
        'on an augmented matrix. Reevaluate the full native EOS/rates/J at every step. '
        'This treats the stiff linearized reaction transient exactly and leaves a nonlinear discretization error.',
        mixed_refinement_budget='1e-16 + 1e-3*abs(total change) per isotope and zone',
        constraints='Same fixed rho,T external heatbath and 22-species model. No EOS/GR physical closure claim. '
        'Mass correction and any tiny negative projection are recorded. Preserve both failed BE refinement gates.'))

def exponential(compensated=False):
    from scipy.linalg import expm
    plan=json.loads((OUT/('reaction-compensated-plan.json' if compensated else 'reaction-exponential-plan.json')).read_text())
    # Independent two-species irreversible conversion has a known analytic flow.
    rate=.7;dt=2.;J=np.array([[-rate,0],[rate,0]]);x=np.array([.8,.2]);aug=np.zeros((3,3))
    aug[:2,:2]=J;aug[:2,2]=J@x;got=x+expm(dt*aug)[:2,2]
    assert np.max(abs(got-np.array([.8*np.exp(-rate*dt),1-.8*np.exp(-rate*dt)])))<1e-14
    base=dict(np.load(OUT/'direct-v2-input.npz'));first=dict(np.load(OUT/'direct-v2-native.npz'));rows=[];previous=None
    for steps in plan['steps']:
        X=base['X'].copy();d=first;dt=plan['duration_proper_s']/steps;records=[]
        accumulated=np.zeros_like(X,dtype=np.longdouble)
        for step in range(steps):
            aug=np.zeros((len(X),23,23));aug[:,:22,:22]=d['jacobian'];aug[:,:22,22]=d['dxdt']
            delta=expm(dt*aug)[:,:22,22]
            mass=-delta.sum(axis=1);largest=np.argmax(X,axis=1);zones=np.arange(len(X));delta[zones,largest]+=mass
            if compensated:
                accumulated+=delta.astype(np.longdouble)
                trial=np.asarray(base['X'].astype(np.longdouble)+accumulated,dtype=float)
            else: trial=X+delta
            assert trial.min()>-1e-22,('exponential negative',steps,step,trial.min())
            projection=np.maximum(-trial,0);trial+=projection;trial[zones,largest]-=projection.sum(axis=1)
            if compensated:
                accumulated+=projection.astype(np.longdouble)
                accumulated[zones,largest]-=projection.sum(axis=1).astype(np.longdouble)
            assert trial.min()>=0 and max(abs(trial.sum(axis=1)-1))<1e-12
            prefix='expc' if compensated else 'exp';label=f'{prefix}-{steps}-{step}';setup(label,{**base,'X':trial});trace(label)
            d=dict(np.load(OUT/(label+'-native.npz')));X=trial
            records.append(dict(label=label,mass_correction_max=float(abs(mass).max()),projection_max=float(projection.max())))
        difference=None if previous is None else float(np.max(abs(X-previous)/(1e-16+1e-3*abs(X-base['X']))))
        np.savez_compressed(OUT/f'{prefix}-{steps}-end.npz',X=X,dxdt=d['dxdt'])
        rows.append(dict(steps=steps,records=records,max_absolute_change=float(abs(X-base['X']).max()),
            mixed_refinement_score=difference,passed=None if difference is None else difference<1,
            max_F_fraction=float(X[:,[ISOS.index(k) for k in ['f17','f18','f19']]].sum(axis=1).max()),
            baryon_max_absolute=float(max(abs(X.sum(axis=1)-1)))))
        previous=X;save('reaction-compensated.json' if compensated else 'reaction-exponential.json',dict(classification='Counterexample candidate',rows=rows,
            analytic_linear_conversion_control=True,completed=len(rows)==len(plan['steps']),
            prior_BE_failed_gates_promoted=False,full_GR_evolution=False))
        print('EXPONENTIAL END',steps,difference,flush=True)

def compensated_plan():
    plan=json.loads((OUT/'reaction-exponential-plan.json').read_text())
    plan.update(previous_exponential_result_sha256=sha(OUT/'reaction-exponential.json'),
        change='Accumulate only the abundance increments in numpy.longdouble, then add to the original state '
        'before converting to native float64 inputs. This avoids repeated rounding of a tiny increment '
        'into a large background abundance. All native function evaluations remain float64.',
        same_mixed_refinement_gate=True)
    save('reaction-compensated-plan.json',plan)

def compensated(): exponential(compensated=True)

def symbolic():
    import sympy as sp
    t,tau,w=sp.symbols('t tau w',positive=True,real=True)
    z=(sp.cos(w*t)+w*tau*sp.sin(w*t))/(1+(w*tau)**2)
    assert sp.simplify(tau*sp.diff(z,t)+z-sp.cos(w*t))==0
    x=sp.symbols('x',positive=True,real=True)
    err=1/(1+sp.I*x)-1
    assert sp.simplify(sp.expand_complex(err*sp.conjugate(err))-x*x/(1+x*x))==0
    r33,r34=sp.symbols('r33 r34',nonnegative=True)
    hdot=2*r33-r34;assert hdot.subs({r33:0,r34:1})==-1
    beta,phi0,dphi,T0,dT=sp.symbols('beta phi0 dphi T0 dT');eps=sp.symbols('eps')
    linear=sp.diff(beta*(phi0+eps*dphi)*(T0+eps*dT),eps).subs(eps,0)
    assert sp.expand(linear.subs(phi0,0))==beta*T0*dphi
    save('symbolic.json',dict(classification='Proven',single_state_monochromatic_identity=True,
        finite_frequency_fast_mode_error='For real tau_j>0, error <= sum_fast abs(a_j)*omega*tau_j.',
        adiabatic_collapse='If omega*tau_slow <<1 the slow state reduces to an instantaneous coefficient at leading order; if its drive or readout residue is zero it is unobservable.',
        reduced_PP_simplex_boundary=str(hdot),
        jump_lemma='If f has one-sided jump J and g is continuous, sup|f-g| >= |J|/2. '
        'Triangle inequality at the two one-sided limits proves this; the measured native jump is a separate numerical result.',
        scalar_zero_branch_linear_thermal_coupling=False,
        scope='Algebraic model statements under their explicit assumptions; not certificates for native EOS, complete GR dynamics, or astronomical observations.'))
    print('PASS symbolic single-state reduction, fast-mode bound, positivity boundary and scalar decoupling',flush=True)

def maintain():
    docs=['model-definition','observable-targets','adiabatic-limit','nonadiabatic-regime','failure-ledger-dynamic-chi']
    old=json.loads((OLD/'manifest.json').read_text())['sha256'];dest=OUT/'request26-notes';dest.mkdir()
    bindings={}
    for rel in ['docs/'+k+'.md' for k in docs]+['paper/revision-manifest.json']:
        assert sha(ROOT/rel)==old[rel],rel
        snap=dest/Path(rel).name;shutil.copy2(ROOT/rel,snap)
        bindings[rel]=dict(snapshot=snap.relative_to(ROOT).as_posix(),sha256=old[rel],historical_manifest=rel.startswith('paper/'))
    save('historical-note-bindings.json',bindings)
    bodies=[
        '분류: Counterexample candidate. Request27은 변경하지 않은 native 반응 함수에서 22종 변화율과 Jacobian을 직접 추출해, 기존 채널 복원 벡터를 독립 검증했다. 생성된 불소를 포함하여 같은 native EOS·반응망에 조성을 되먹이는 고정 밀도·온도 대조를 수행했다. 기존 축약 PP 망은 H=0 경계에서 바깥 방향 변화율을 보였고, 중간 핵종 4개를 명시한 26종 망은 같은 경계 대조를 통과했다.',
        '분류: Counterexample candidate. 별도로 초기화한 국소 반응 부분계에서 외부 온도 구동→핵 가열 응답의 느린 단일 상태를 분리했다. 최대 가열 구역의 고유시간은 약 5.248988일이며, 고정 GR 배경의 좌표시간으로 약 5.249074일이다. 입출력 미분을 native 재평가로 대조했다. 이는 핵 가열 전달함수이며 scalar 전하·자유낙하 힘·표면 광도 전달함수는 도출하지 않았다.',
        '분류: Proven. 안정한 실수 모드들의 응답이 D+Σaⱼ/(1+iωτⱼ)이면 빠른 모드를 정적 계수로 치환한 오차는 Σfast|aⱼ|ωτⱼ 이하이다. 느린 모드도 ωτ≪1에서는 선도 정적 계수로 붕괴하며, 그 구동 또는 출력 잔차가 0이면 관측되지 않는다. 영 scalar 배경의 선형 열 분리 경계는 새 핵 반응 상태의 존재로 해제되지 않는다.',
        '분류: Counterexample candidate. 26종 망의 명시적 중간 핵종에서 일 단위의 핵 가열 응답을 계산했다. 지정 두 주파수에서 4상태→1상태 축약의 상대 오차는 약 1.49e-6과 7.32e-9이며 빠른 모드 오차식 안에 든다. 전 주파수의 동등성이나 유한 carrier에서 자유로운 정적 고차 계수로 흡수되지 않는다는 주장은 하지 않는다.',
        '분류: Counterexample candidate. 직접 반응 벡터 검증의 미해결 항목은 지정 초기 상태에서 해소했다. 엄격 상대 잔차, 반올림만 포함한 잔차, 후진 Euler의 조성별 시간 정밀화 실패를 각각 보존했다. EOS 차분은 큰 간격에서도 실패했고 REAL 입력 경계에 걸친 출력 점프를 실제 측정했다. 이는 물리적 연속 EOS의 불연속성이나 제1법칙 위반의 증명이 아니다.\n\n분류: Conjectural. 생성 불소와 PP 중간 핵종을 포함한 공통 EOS의 정칙성·오차 제어, 보존된 항성 초기화·열/유체/metric 진화 및 비영 scalar 배경의 구동/전하 연결이 남는다. 국소 핵 가열 상태를 완전한 관측 추론으로 승격하지 않는다.'
    ]
    rev=json.loads((ROOT/'paper/revision-manifest.json').read_text())
    for stem,body in zip(docs,bodies):
        rel='docs/'+stem+'.md'
        with (ROOT/rel).open('ab') as f: f.write(('\n\n## Request 27 직접 반응 벡터와 명시적 PP 상태\n\n'+body+'\n\n세부 근거: [한글 보고서](../notes/REQUEST27_NATIVE_CLOSURE_KO.md).\n').encode())
        rev['sha256'][rel]=sha(ROOT/rel)
    rev['request27_supporting_note_update']=dict(evidence_manifest='outputs/native-closure27/manifest.json',
        historical_notes='outputs/native-closure27/historical-note-bindings.json',
        status='직접 반응 벡터 독립 검증 및 명시적 핵 가열 상태 후보; 공통 EOS와 실제 GR/관측 폐쇄 미완료',
        artifact_status='원고 PDF/ZIP은 Request12 역사 산출물이며 이번 연구 노트 갱신과 구분한다.')
    (ROOT/'paper/revision-manifest.json').write_text(json.dumps(rev,ensure_ascii=False,indent=2)+'\n')

def source_bindings():
    for rel in ['net/public/net_lib.f90','star/private/net.f90','star/private/profile_getval.f90',
        'eos/private/eosdt_eval.f90','eos/private/helm_electron_positron.dek',
        'data/net_data/nets/add_pp_extras','data/rates_data/reactions.list','const/public/const_def.f90']:
        p=OUT/'sources'/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(fresh.MESA/rel,p)
    save('source-bindings.json',dict(classification='Imported from prior work',
        source_sha256={p.relative_to(OUT/'sources').as_posix():sha(p) for p in (OUT/'sources').rglob('*') if p.is_file()},
        literature=[dict(title='Paxton et al. MESA (2011)',url='https://arxiv.org/abs/1009.1622'),
            dict(title='MESA reaction-network documentation',url='https://docs.mesastar.org/en/latest/net/nets.html')],
        note='Current documentation supplies general context only; exact old network definitions, source and unchanged executable establish the numerical boundary.'))

def audit():
    counts=dict(native_profiles=0,untraced_profiles=0,native_zone_calls=0);worst=0.
    for p in sorted(OUT.glob('*-import-control.json')):
        label=p.name[:-len('-import-control.json')];assert json.loads(p.read_text())['passed']
        inp=np.load(OUT/(label+'-input.npz'));species=ISOS
        if (OUT/(label+'-species.json')).exists(): species=json.loads((OUT/(label+'-species.json')).read_text())
        _,prof=mesa(OUT/(label+'-profile.data.gz'));assert len(prof['zone'])==len(inp['dm'])
        path=OUT/(label+'-native.npz')
        if not path.exists(): counts['untraced_profiles']+=1;continue
        d=np.load(path);ns=len(species);assert d['jacobian'].shape==(len(inp['dm']),ns,ns)
        assert np.max(abs(np.log(d['T'])-inp['lnT']))<1e-12
        assert np.max(abs(np.log(d['rho'])-inp['lnd']))<1e-12
        assert np.max(abs(d['X']-inp['X']))<1e-12
        for key,col in [('heat','eps_nuc'),('neutrino','eps_nuc_neu_total')]:
            err=float(np.max(abs(d[key]-prof[col])/np.maximum(1,abs(prof[col]))));worst=max(worst,err);assert err<1e-12
        counts['native_profiles']+=1;counts['native_zone_calls']+=len(inp['dm'])
    save('trace-audit.json',dict(classification='Proven',**counts,max_direct_vs_profile_relative=worst,
        all_complete_traces_checked=True,abi_scope='SHA-bound x86-64 SysV/gfortran descriptor layout only.'))
    print('AUDIT',counts,'max direct/profile',worst,flush=True)

def seal():
    audit();exp=json.loads((OUT/'reaction-compensated.json').read_text())
    assert exp['completed']
    save('gates.json',dict(classification='Proven',
        native_22_vector_independently_validated=json.loads((OUT/'direct-vector-validation.json').read_text())['independently_validated'],
        common_22_species_evaluation_path=True,
        exponential_fixed_rho_T_refinement_passed=exp['rows'][-1]['passed'],
        original_uncompensated_exponential_gate=False,
        explicit_PP_boundary_control_passed=json.loads((OUT/'explicit-PP.json').read_text())['zero_H_boundary_inward'],
        PP_steady_control_passed=json.loads((OUT/'PP-steady.json').read_text())['passed'],
        PP_temperature_derivatives_checked=json.loads((OUT/'PP-response.json').read_text())['native_derivative_check_passed'],
        PP_peak_composition_derivatives_checked=json.loads((OUT/'PP-composition-derivatives.json').read_text())['passed'],
        common_EOS_thermodynamics_certified=False,global_derivative_certificate=False,
        full_GR_thermal_fluid_metric_evolution=False,scalar_drive_charge_map=False,
        complete_nonlinear_observation=False,final_submission_package_updated=False,
        classification_detail='Theorem progress: fast-mode reduction error, positivity boundary and zero-field decoupling. '
        'Loophole progress: native vector validation, generated-species continuation and explicit nuclear heat-memory candidate.'))
    save('provenance.json',dict(classification='Proven',checkpoint='3b70803',previous_manifest_sha256=sha(OLD/'manifest.json'),
        executable_sha256=sha(fresh.BINARY),runtime_data_binding='outputs/fresh-microphysics25/runtime-data-bindings.json',
        historical_policies='Preserve all previous failures and hypotheses. No new build, runtime installation, or empirical pulsar/LLR work.',
        source_version_note='The reused data/version_number says 7623; the source folder is named r7624. No independent build certificate is asserted.'))
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json']
    files += [ROOT/'verification/native_closure.py',ROOT/'notes/REQUEST27_NATIVE_CLOSURE_KO.md',ROOT/'.gitattributes']
    files += [ROOT/k for k in json.loads((OUT/'historical-note-bindings.json').read_text())]
    save('manifest.json',dict(classification='Proven',sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sorted(files)}))

def verify():
    stages=['validated-variational','remaining-levers15','nbody-readout16','nonzero-drive17','thermal-wd18',
        'thermal-restart19','thermal-robustness20','gr-mass21','thermal-closure22','baryon-entropy23',
        'reactive-energy24','fresh-microphysics25','remaining-closure26','native-closure27']
    histories={f'outputs/{a}/manifest.json':f'outputs/{b}/historical-note-bindings.json' for a,b in zip(stages,stages[1:])};count=0
    for label in ['outputs/research-remediation/manifest.json',*histories,'paper/revision-manifest.json','outputs/native-closure27/manifest.json']:
        old=json.loads((ROOT/histories[label]).read_text()) if label in histories else {}
        for name,expected in json.loads((ROOT/label).read_text())['sha256'].items():
            path=ROOT/name
            if name in old:
                bind=old[name];path=ROOT/bind['snapshot'];assert expected==bind['sha256']
                if bind.get('historical_manifest'):
                    before=json.loads(path.read_text());after=json.loads((ROOT/name).read_text())
                    for k,v in before.items():
                        if k!='sha256': assert after[k]==v,k
                    for k,v in before['sha256'].items():
                        if k not in old: assert after['sha256'][k]==v,k
                else: assert (ROOT/name).read_bytes().startswith(path.read_bytes())
            assert sha(path)==expected,(label,name);count+=1
    assert sha(fresh.BINARY)==json.loads((OUT/'provenance.json').read_text())['executable_sha256']
    data=json.loads((ROOT/'outputs/fresh-microphysics25/runtime-data-bindings.json').read_text())
    for path,digest in data['sha256'].items(): assert sha(Path(data['root'])/path)==digest,path
    for name in ['strict-residual-failure.json','arithmetic-residual-failure.json']:
        assert not json.loads((OUT/name).read_text())['passed']
    assert not json.loads((OUT/'network-boundary.json').read_text())['simplex_inward']
    assert not any(r['passed'] for r in json.loads((OUT/'EOS-continuation.json').read_text())['rows'])
    gates=json.loads((OUT/'gates.json').read_text())
    for key in ['common_EOS_thermodynamics_certified','global_derivative_certificate','full_GR_thermal_fluid_metric_evolution',
        'scalar_drive_charge_map','complete_nonlinear_observation','final_submission_package_updated']: assert not gates[key]
    assert sha(OUT/'composition-continuation-summary.json')==json.loads((OUT/'reaction-exponential-plan.json').read_text())['prior_BE_failure_sha256']
    print('PASS',count,'artifact/history SHA;',len(data['sha256']),'MESA data SHA; original failed gates retained',flush=True)

def recheck():
    import tempfile, contextlib, io
    global OUT,CACHE
    original,cache=OUT,CACHE
    labels=['direct-v2','H-zero-boundary','explicit-PP-Hzero','PP-steady-2','PP-X-be7-1',
        'rounding-1e-11--1','rounding-1e-11-1','expc-4-3']
    with tempfile.TemporaryDirectory(prefix='native-closure27-') as temp:
        OUT=Path(temp)/'out';CACHE=Path(temp)/'runs';CACHE.mkdir();shutil.copytree(original,OUT);context()
        try:
            worst=0.
            for label in labels:
                folder=CACHE/label;folder.mkdir();shutil.copy2(fresh.BINARY,folder/'binary')
                for p in (original/'inputs'/label).iterdir():
                    if p.name=='input.mod.gz':
                        with fresh.gzip.open(p,'rb') as src,(folder/'input.mod').open('wb') as dst: shutil.copyfileobj(src,dst)
                    else: shutil.copy2(p,folder/p.name)
                traced=(original/(label+'-native.npz')).exists()
                with contextlib.redirect_stdout(io.StringIO()):
                    if traced: trace(label)
                    else: context();fresh.run(label);fresh.collect(label)
                if traced:
                    a=np.load(OUT/(label+'-native.npz'));b=np.load(original/(label+'-native.npz'))
                    for key in a.files:
                        err=float(np.max(abs(a[key]-b[key])/np.maximum(1e-100,abs(b[key]))));worst=max(worst,err)
                        assert err<1e-10,(label,key,err)
                _,a=mesa(OUT/(label+'-profile.data.gz'));_,b=mesa(original/(label+'-profile.data.gz'))
                for key in ['rho','logT','eps_nuc','eps_nuc_neu_total','energy','pressure']:
                    assert np.max(abs(a[key]-b[key])/np.maximum(1e-100,abs(b[key])))<1e-10,(label,key)
                print('RECHECK native state',label,flush=True)
            with contextlib.redirect_stdout(io.StringIO()):
                analyze();eos_analysis();explicit_analysis();pp_response();pp_reduction();symbolic()
            print('PASS 8 actual native state reruns and symbolic/response recomputation; max difference',worst,flush=True)
        finally: OUT,CACHE=original,cache;context()

if __name__=='__main__': globals()[sys.argv[1]]()
