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
        executable_sha256=g.c.sha(g.c.fresh.CACHE/'binary') if (g.c.fresh.CACHE/'binary').exists() else
            json.loads((g.ROOT/'outputs/fresh-microphysics25/runtime.json').read_text())['sha256'],
        call='__opacities_MOD_eval_kap_type2, 21 by-reference gfortran arguments; use the explicit cell index and composition descriptor to bind each call.',
        controls='First read-only baseline: actual rho,T,26-species X, species sums, and returned opacity against the collected profile. A later identity-input replay must reproduce outputs before changing EOS electron inputs.',
        argument_slots=dict(cell=1,zbar=2,XH=3,Z=4,Zbase=5,XC=6,XN=7,XO=8,XNe=9,
            log10_rho=10,log10_T=11,composition=12,lnfree_e=13,free_e_rho=14,free_e_T=15,
            fraction_type2=16,opacity=17,opacity_rho=18,opacity_T=19,ierr=20),
        input_log_tolerance=1e-12,composition_absolute_tolerance=1e-12,profile_opacity_relative_tolerance=1e-12,
        source_contract='Total free electron plus positron number per nucleon, not merely nuclear net Ye. Original total opacity already includes conduction.',
        physical_opacity_certified=False,nonlinear_transport_completed=False))


def setup(label):
    g.c.native.OUT=OUT;g.c.native.CACHE=CACHE
    g.c.native.setup(label,dict(np.load(g.OUT/'reference-state.npz')),species=g.c.NAMES,
        network=(g.c.cell.OLD/'inputs/explicit-PP/cno_extras.net').read_text())


def trace(label,replacement=None):
    folder=CACHE/label;plan=json.loads((OUT/'plan.json').read_text())
    expected=json.loads((OUT/(label+'-inputs.json')).read_text())['sha256']
    for name,digest in expected.items(): assert g.c.sha(folder/name)==digest,name
    assert expected['binary']==plan['executable_sha256']
    symbols=subprocess.check_output(['nm',str(folder/'binary')],text=True)
    lines=[line for line in symbols.splitlines() if line.endswith(' __opacities_MOD_eval_kap_type2')]
    assert len(lines)==1;entry=int(lines[0].split()[0],16)
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
    begin=time.monotonic();records=[];mem=None;reg=Registers()
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
        def integer(address): return struct.unpack('<i',read(address,4))[0]
        def getregs(): ptrace(12,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def rewind(address): reg.rip=address;ptrace(13,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def breakpoint(address):
            word=u64(address);ptrace(4,pid,address,(word & ~255)|0xcc);return word
        try:
            wait();mem=os.open('/proc/'+str(pid)+'/mem',os.O_RDWR if replacement is not None else os.O_RDONLY)
            original=breakpoint(entry);ptrace(7,pid)
            while len(records)<n:
                wait();getregs();assert reg.rip==entry+1,hex(reg.rip)
                args=[reg.rdi,reg.rsi,reg.rdx,reg.rcx,reg.r8,reg.r9]+[u64(reg.rsp+8*(k-5)) for k in range(6,21)]
                i=integer(args[1])-1;assert i==len(records),(i,len(records))
                descriptor=struct.unpack('<6q',read(args[12],48))
                assert descriptor[3]==1 and descriptor[5]-descriptor[4]+1==ns,descriptor
                x=np.frombuffer(read(descriptor[0],8*ns),dtype='<f8').copy()
                parameters=np.array([real(args[k]) for k in list(range(2,12))+list(range(13,16))])
                assert np.all(np.isfinite(parameters))
                assert max(abs(parameters[8:10]*np.log(10)-[data['lnd'][i],data['lnT'][i]]))<plan['input_log_tolerance']
                assert max(abs(x-data['X'][i]))<plan['composition_absolute_tolerance']
                saved=[]
                if replacement is not None:
                    for slot,value in zip([13,14,15],replacement[i]):
                        saved.append((args[slot],read(args[slot],8)))
                        assert os.pwrite(mem,struct.pack('<d',float(value)),args[slot])==8
                used=np.array([real(args[k]) for k in [13,14,15]])
                if replacement is not None: assert np.array_equal(used,replacement[i])
                ret=u64(reg.rsp);ptrace(4,pid,entry,original);rewind(entry);ptrace(9,pid);wait()
                retword=breakpoint(ret);ptrace(7,pid);wait();getregs();assert reg.rip==ret+1
                assert integer(args[20])==0
                result=np.array([real(args[k]) for k in [16,17,18,19]])
                assert np.all(np.isfinite(result)) and result[1]>0
                for address,value in saved:
                    assert os.pwrite(mem,value,address)==8 and read(address,8)==value
                records.append(dict(cell=i,parameters=parameters,X=x,used=used,outputs=result))
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
            shutil.copy2(folder/'execution.log',OUT/(label+'-execution.log'))
    assert g.c.sha(folder/'binary')==expected['binary']
    g.c.fresh.OUT=OUT;g.c.fresh.CACHE=CACHE;g.c.fresh.ISOS=list(g.c.NAMES);g.c.fresh.collect(label)
    _,profile=g.c.mesa(OUT/(label+'-profile.data.gz'))
    values=np.array([r['outputs'][1] for r in records]);score=float(abs(values/profile['opacity']-1).max())
    save(label+'-trace.json',dict(classification='Counterexample candidate',calls=len(records),entry_hex=hex(entry),
        executable_unchanged=True,inputs_restored_after_each_call=True,returned_results_overwritten=False,
        replacement_used=replacement is not None,profile_opacity_relative_difference=score,
        profile_passed=score<plan['profile_opacity_relative_tolerance'],physical_or_nonlinear_certificate=False))
    assert score<plan['profile_opacity_relative_tolerance'];print('OPACITY TRACE',label,len(records),score,flush=True)


def baseline(): setup('baseline');trace('baseline')


def identity():
    baseline=dict(np.load(OUT/'baseline-captured.npz'));assert json.loads((OUT/'baseline-trace.json').read_text())['profile_passed']
    setup('identity');trace('identity',baseline['used'])
    replay=dict(np.load(OUT/'identity-captured.npz'))
    equal={key:bool(np.array_equal(baseline[key],replay[key])) for key in baseline}
    save('identity-control.json',dict(classification='Counterexample candidate',passed=all(equal.values()),bitwise_equal=equal))
    assert all(equal.values()),equal


if __name__=='__main__': globals()[sys.argv[1]]()
