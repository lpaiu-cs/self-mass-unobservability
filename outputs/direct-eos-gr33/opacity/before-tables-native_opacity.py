"""Read/change explicit opacity inputs in a disposable SHA-bound stellar child.

Counterexample candidate. Reuses the existing native model setup, collection and
x86-64 register layout. No executable file or returned result is overwritten.
"""
import ctypes as ct
import json, os, shutil, signal, struct, subprocess, sys, time
import numpy as np
import direct_eos_gr as g
from native_closure import Registers

OUT=g.OUT/'opacity';CACHE=g.CACHE/'opacity'


def save(name,value):
    (OUT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')


def prepare():
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True);CACHE.mkdir(exist_ok=True)
    save('plan.json',dict(classification='Counterexample candidate',
        executable_sha256=json.loads((g.ROOT/'outputs/fresh-microphysics25/runtime.json').read_text())['sha256'],
        call='__opacities_MOD_eval_kap_type1, 17 by-reference gfortran arguments; use the explicit cell index and composition descriptor to bind each call.',
        controls='First read-only baseline: actual rho,T,26-species X, species sums, and returned opacity against the collected profile. A later identity-input replay must reproduce outputs before changing EOS electron inputs.',
        argument_slots=dict(cell=1,zbar=2,XH=3,Z=4,q=5,Pgas_div_P=6,
            log10_rho=7,log10_T=8,composition=9,lnfree_e=10,free_e_rho=11,free_e_T=12,
            opacity=13,opacity_rho=14,opacity_T=15,ierr=16),
        input_log_tolerance=1e-12,composition_absolute_tolerance=1e-12,profile_opacity_relative_tolerance=1e-12,
        source_contract='Total free electron plus positron number per nucleon, not merely nuclear net Ye. Original total opacity already includes conduction.',
        physical_opacity_certified=False,nonlinear_transport_completed=False))


def setup(label):
    g.c.native.OUT=OUT;g.c.native.CACHE=CACHE
    g.c.native.setup(label,dict(np.load(g.OUT/'reference-state.npz')),species=g.c.NAMES,
        network=(g.c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())


def trace(label,replacement=None,inspect_internal=False,capture_cubic=False):
    folder=CACHE/label;plan=json.loads((OUT/'plan.json').read_text())
    expected=json.loads((OUT/(label+'-inputs.json')).read_text())['sha256']
    for name,digest in expected.items(): assert g.c.sha(folder/name)==digest,name
    assert expected['binary']==plan['executable_sha256']
    symbols=subprocess.check_output(['nm',str(folder/'binary')],text=True)
    lines=[line for line in symbols.splitlines() if line.endswith(' __opacities_MOD_eval_kap_type1')]
    assert len(lines)==1;entry=int(lines[0].split()[0],16)
    symbol_addresses={line.split()[-1]:int(line.split()[0],16) for line in symbols.splitlines() if len(line.split())==3}
    inner_entry=symbol_addresses['__kap_eval_MOD_combine_rad_with_compton_and_conduction']
    vector_entry=symbol_addresses['__interp_1d_lib_MOD_interpolate_vector']
    libc=ct.CDLL(None,use_errno=True);libc.ptrace.restype=ct.c_long
    libc.ptrace.argtypes=[ct.c_ulong,ct.c_ulong,ct.c_void_p,ct.c_void_p]
    def ptrace(request,pid,address=0,value=0):
        ct.set_errno(0);result=libc.ptrace(request,pid,address,value)
        if result==-1 and ct.get_errno(): raise OSError(ct.get_errno(),os.strerror(ct.get_errno()))
        return result
    def child(): ptrace(0,0)
    runtime=g.c.fresh.MESA.parent;env=os.environ.copy()
    env.update(MESA_DIR=str(g.c.fresh.MESA),LD_LIBRARY_PATH=str(runtime/'mesasdk/lib')+':'+str(runtime/'mesasdk/lib64'),OMP_NUM_THREADS='1')
    assert not (folder/'execution.log').exists()
    data=dict(np.load(OUT/(label+'-input.npz')));n=len(data['dm']);ns=data['X'].shape[1]
    if replacement is not None: assert replacement.shape==(n,3) and np.all(np.isfinite(replacement))
    begin=time.monotonic();records=[];mem=None;reg=Registers();global_controls={};cubic=[]
    with (folder/'execution.log').open('w') as log:
        proc=subprocess.Popen(['./binary'],cwd=folder,env=env,stdout=log,stderr=subprocess.STDOUT,preexec_fn=child);pid=proc.pid
        def wait():
            while time.monotonic()-begin<180:
                done,status=os.waitpid(pid,os.WNOHANG)
                if done:
                    assert os.WIFSTOPPED(status),(label,status,len(records))
                    if os.WSTOPSIG(status)==signal.SIGCHLD: ptrace(7,pid,0,signal.SIGCHLD);continue
                    assert os.WSTOPSIG(status)==signal.SIGTRAP,(label,status);return
                time.sleep(.0001)
            raise TimeoutError(label)
        def read(address,size):
            result=os.pread(mem,size,address);assert len(result)==size;return result
        def u64(address): return struct.unpack('<Q',read(address,8))[0]
        def real(address): return struct.unpack('<d',read(address,8))[0]
        def single(address): return struct.unpack('<f',read(address,4))[0]
        def i64(address): return struct.unpack('<q',read(address,8))[0]
        def integer(address): return struct.unpack('<i',read(address,4))[0]
        def vector(address,n):
            descriptor=struct.unpack('<6q',read(address,48))
            assert descriptor[3]==1 and descriptor[5]-descriptor[4]+1==n,descriptor
            return np.frombuffer(read(descriptor[0],8*n),dtype='<f8').copy()
        def getregs(): ptrace(12,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def rewind(address): reg.rip=address;ptrace(13,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def breakpoint(address):
            word=u64(address);ptrace(4,pid,address,(word & ~255)|0xcc);return word
        try:
            wait();mem=os.open('/proc/'+str(pid)+'/mem',os.O_RDWR if replacement is not None else os.O_RDONLY)
            original=breakpoint(entry);ptrace(7,pid)
            while len(records)<n:
                wait();getregs();assert reg.rip==entry+1,hex(reg.rip)
                args=[reg.rdi,reg.rsi,reg.rdx,reg.rcx,reg.r8,reg.r9]+[u64(reg.rsp+8*(k-5)) for k in range(6,17)]
                i=integer(args[1])-1;assert i==len(records),(i,len(records))
                descriptor=struct.unpack('<6q',read(args[9],48))
                assert descriptor[3]==1 and descriptor[5]-descriptor[4]+1==ns,descriptor
                x=np.frombuffer(read(descriptor[0],8*ns),dtype='<f8').copy()
                parameters=np.array([real(args[k]) for k in [2,3,4,7,8,10,11,12]])
                assert np.all(np.isfinite(parameters))
                assert max(abs(parameters[3:5]*np.log(10)-[data['lnd'][i],data['lnT'][i]]))<plan['input_log_tolerance']
                assert max(abs(x-data['X'][i]))<plan['composition_absolute_tolerance']
                saved=[]
                if replacement is not None:
                    for slot,value in zip([10,11,12],replacement[i]):
                        saved.append((args[slot],read(args[slot],8)))
                        assert os.pwrite(mem,struct.pack('<d',float(value)),args[slot])==8
                used=np.array([real(args[k]) for k in [10,11,12]])
                if replacement is not None: assert np.array_equal(used,replacement[i])
                ret=u64(reg.rsp);ptrace(4,pid,entry,original);rewind(entry);ptrace(9,pid);wait()
                retword=breakpoint(ret)
                if inspect_internal: innerword=breakpoint(inner_entry)
                if capture_cubic: vectorword=breakpoint(vector_entry)
                ptrace(7,pid);wait();getregs()
                while capture_cubic and reg.rip==vector_entry+1:
                    vector_args=[reg.rdi,reg.rsi,reg.rdx,reg.rcx,reg.r8,reg.r9]+[u64(reg.rsp+8*(k-5)) for k in range(6,11)]
                    assert integer(vector_args[0])==4 and integer(vector_args[2])==1
                    name=read(vector_args[9],19).decode('ascii')
                    assert name in ['Get_Kap_for_X_cubic','Get_Kap_for_Z_cubic'],name
                    item=dict(cell=i,axis=name[12],x=vector(vector_args[1],4),
                        at=float(vector(vector_args[3],1)[0]),values=vector(vector_args[4],4))
                    vector_ret=u64(reg.rsp);ptrace(4,pid,vector_entry,vectorword);rewind(vector_entry)
                    ptrace(9,pid);wait();vector_retword=breakpoint(vector_ret);ptrace(7,pid);wait();getregs()
                    assert reg.rip==vector_ret+1 and integer(vector_args[10])==0
                    item['result']=float(vector(vector_args[5],1)[0]);cubic.append(item)
                    ptrace(4,pid,vector_ret,vector_retword);rewind(vector_ret);breakpoint(vector_entry)
                    ptrace(7,pid);wait();getregs()
                if inspect_internal:
                    assert reg.rip==inner_entry+1,hex(reg.rip)
                    inner_args=[reg.rdi,reg.rsi,reg.rdx,reg.rcx,reg.r8,reg.r9]+[u64(reg.rsp+8*(k-5)) for k in range(6,16)]
                    inner=np.array([single(inner_args[k]) for k in range(1,6)]+[real(inner_args[k]) for k in range(6,12)])
                    assert np.all(np.isfinite(inner)) and inner[8]>0
                    if i==0:
                        # Fixed SHA-bound gfortran layout, verified from kap_def
                        # and get_logT_Compton_blend_hi disassembly: 64/192 bytes.
                        for name in ['kap_z_tables','kap_lowt_z_tables']:
                            address=symbol_addresses['__kap_def_MOD_'+name]
                            z=u64(address)+64*(i64(address+8)+i64(address+24))
                            xfirst=u64(z+16)+192*(i64(z+24)+i64(z+40))
                            global_controls[name]=dict(logT_min=single(xfirst+80),logT_max=single(xfirst+84))
                            count=i64(address+40)-i64(address+32)+1;assert 1<=count<=100
                            global_controls[name]['Z']=[single(z+64*j+4) for j in range(count)]
                        for name in ['kap_blend_logt_lower_bdy','kap_blend_logt_upper_bdy']:
                            global_controls[name]=real(symbol_addresses['__kap_def_MOD_'+name])
                        global_controls['compton_blend_hi']=global_controls['kap_z_tables']['logT_max']-.01
                        request=u64(inner_args[0])
                        global_controls['choices']={name:integer(request+4*j) for j,name in enumerate([
                            'cubic_interpolation_in_X','cubic_interpolation_in_Z','include_electron_conduction'])}
                        assert all(v in [0,1] for v in global_controls['choices'].values())
                    ptrace(4,pid,inner_entry,innerword);rewind(inner_entry);ptrace(7,pid);wait();getregs()
                assert reg.rip==ret+1
                if capture_cubic: ptrace(4,pid,vector_entry,vectorword)
                assert integer(args[16])==0
                result=np.array([real(args[k]) for k in [13,14,15]])
                assert np.all(np.isfinite(result)) and result[0]>0
                for address,value in saved:
                    assert os.pwrite(mem,value,address)==8 and read(address,8)==value
                records.append(dict(cell=i,parameters=parameters,X=x,used=used,outputs=result))
                if inspect_internal: records[-1]['inner']=inner
                ptrace(4,pid,ret,retword);rewind(ret)
                if len(records)==n: ptrace(17,pid);break
                breakpoint(entry);ptrace(7,pid)
            complete=False
            while time.monotonic()-begin<180:
                profile=folder/'LOGS1/profile1.data'
                if profile.exists():
                    try:
                        h,p=g.c.mesa(profile);complete=int(h['model_number'])==1 and len(p['zone'])==n and p['zone'][-1]==n
                    except (ValueError,IndexError,KeyError,OSError): pass
                    if complete: break
                if proc.poll() is not None: break
                time.sleep(.05)
            assert complete,(label,'profile incomplete')
        finally:
            if mem is not None: os.close(mem)
            if proc.poll() is None:
                os.kill(pid,signal.SIGKILL)
                try: os.waitpid(pid,0)
                except ChildProcessError: pass
            if records: np.savez_compressed(OUT/(label+'-captured.npz'),**{k:np.array([r[k] for r in records]) for k in records[0]})
            if cubic: np.savez_compressed(OUT/(label+'-cubic.npz'),**{k:np.array([r[k] for r in cubic]) for k in cubic[0]})
            shutil.copy2(folder/'execution.log',OUT/(label+'-execution.log'))
    assert g.c.sha(folder/'binary')==expected['binary']
    g.c.fresh.OUT=OUT;g.c.fresh.CACHE=CACHE;g.c.fresh.ISOS=list(g.c.NAMES);g.c.fresh.collect(label)
    _,profile=g.c.mesa(OUT/(label+'-profile.data.gz'))
    values=np.array([r['outputs'][0] for r in records]);score=float(abs(values/profile['opacity']-1).max())
    save(label+'-trace.json',dict(classification='Counterexample candidate',calls=len(records),entry_hex=hex(entry),
        executable_unchanged=True,inputs_restored_after_each_call=True,returned_results_overwritten=False,
        replacement_used=replacement is not None,profile_opacity_relative_difference=score,
        internal_inputs_captured=inspect_internal,global_controls=global_controls,
        cubic_calls=len(cubic),
        profile_passed=score<plan['profile_opacity_relative_tolerance'],physical_or_nonlinear_certificate=False))
    assert score<plan['profile_opacity_relative_tolerance'];print('OPACITY TRACE',label,len(records),score,flush=True)


def baseline(): setup('baseline-type1');trace('baseline-type1')


def identity():
    baseline=dict(np.load(OUT/'baseline-type1-captured.npz'));assert json.loads((OUT/'baseline-type1-trace.json').read_text())['profile_passed']
    setup('identity-type1');trace('identity-type1',baseline['used'])
    replay=dict(np.load(OUT/'identity-type1-captured.npz'))
    equal={key:bool(np.array_equal(baseline[key],replay[key])) for key in baseline}
    save('identity-control.json',dict(classification='Counterexample candidate',passed=all(equal.values()),bitwise_equal=equal))
    assert all(equal.values()),equal


def composition_control():
    values=dict(np.load(OUT/'baseline-type1-captured.npz'));x=values['X'];p=values['parameters']
    nuclei=x/g.c.A;zbar=(nuclei@g.c.Z)/nuclei.sum(1)
    X=x[:,g.c.Z==1].sum(1);He=x[:,g.c.Z==2].sum(1);Z=np.clip(1-X-He,0,.1)
    errors=dict(zbar=float(abs(p[:,0]-zbar).max()),hydrogen=float(abs(p[:,1]-X).max()),
        metallicity=float(abs(p[:,2]-Z).max()))
    passed=max(errors.values())<1e-12
    save('composition-control.json',dict(classification='Counterexample candidate',passed=passed,
        absolute_errors=errors,scope='Actual native Type1 composition inputs versus captured 26-species inventories. The nuclear mean charge does not certify a partially ionized conduction model.'))
    assert passed,errors;print('OPACITY COMPOSITION',errors,flush=True)


def repair_call_selection():
    plan=json.loads((OUT/'plan.json').read_text())
    source=g.c.fresh.MESA/'star/defaults/controls.defaults'
    assert 'use_Type2_opacities = .false.' in source.read_text()
    assert 'use_Type2_opacities' not in (OUT/'inputs/baseline/inlist1').read_text()
    plan['call']='__opacities_MOD_eval_kap_type1; actual Type1 path, 17 by-reference arguments.'
    plan['argument_slots']=dict(cell=1,zbar=2,XH=3,Z=4,q=5,Pgas_div_P=6,
        log10_rho=7,log10_T=8,composition=9,lnfree_e=10,free_e_rho=11,free_e_T=12,
        opacity=13,opacity_rho=14,opacity_T=15,ierr=16)
    plan['retained_failure']='Type2 breakpoint had zero hits before normal native exit; it was not the active opacity path.'
    plan['before_plan_sha256']=g.c.sha(OUT/'before-type1-plan.json')
    shutil.copy2(source,OUT/'controls.defaults');plan['controls_source_sha256']=g.c.sha(OUT/'controls.defaults')
    save('plan.json',plan)


if __name__=='__main__': globals()[sys.argv[1]]()
