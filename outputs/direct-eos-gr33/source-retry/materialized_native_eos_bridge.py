"""Request29: explicit EOS-input shim for the SHA-bound native net_get.

Counterexample candidate. Changes four EOS INPUTS (optionally also three
composition moments) in our disposable child, records and restores them.
No returned numerical results are overwritten.
"""
from pathlib import Path
import ctypes as ct
import json, os, shutil, signal, struct, subprocess, time
import numpy as np
import fresh_microphysics as fresh
from native_closure import Registers
from thermal_restart import sha
from thermal_wd import mesa
from common_eos import ROOT, OUT, CACHE, save

def context():
    fresh.OUT=OUT;fresh.CACHE=CACHE

def trace(label, replacement=None, weak_blend=None, capture_weak_table=False):
    context();folder=CACHE/label
    if (OUT/(label+'-species.json')).exists(): fresh.ISOS=json.loads((OUT/(label+'-species.json')).read_text())
    expected=json.loads((OUT/(label+'-inputs.json')).read_text())['sha256']
    for name,digest in expected.items(): assert sha(folder/name)==digest,name
    assert sha(folder/'binary')==json.loads((ROOT/'outputs/fresh-microphysics25/runtime.json').read_text())['sha256']
    symbols=subprocess.check_output(['nm',str(folder/'binary')],text=True)
    def address(name):
        found=[s for s in symbols.splitlines() if s.endswith(' '+name)]
        assert len(found)==1,name
        return int(found[0].split()[0],16)
    line=[s for s in symbols.splitlines() if s.endswith(' __net_lib_MOD_net_get')]
    assert len(line)==1;entry=int(line[0].split()[0],16)
    blend_addresses=[]
    if weak_blend is not None:
        assert 0<weak_blend[0]<weak_blend[1]
        for suffix in ['off','on']:
            name='__rates_def_mic_MOD_t9_weaklib_full_'+suffix
            found=[s for s in symbols.splitlines() if s.endswith(' '+name)]
            assert len(found)==1,name
            blend_addresses.append(int(found[0].split()[0],16))
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
        def table_array(name,rank,dtype):
            desc=address('__rates_def_mic_MOD_'+name)
            header=struct.unpack('<'+'q'*(3+3*rank),read(desc,8*(3+3*rank)))
            shape=tuple(header[5+3*i]-header[4+3*i]+1 for i in range(rank))
            assert all(0<n<=100000 for n in shape) and np.prod(shape)<10000000,(name,header)
            assert header[3]==1
            for i in range(1,rank): assert header[3+3*i]==int(np.prod(shape[:i])),(name,header)
            return np.frombuffer(read(header[0],np.dtype(dtype).itemsize*int(np.prod(shape))),dtype=dtype).copy().reshape(shape,order='F')
        def getregs(): ptrace(12,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def rewind(addr):
            reg.rip=addr;ptrace(13,pid,0,ct.cast(ct.byref(reg),ct.c_void_p))
        def breakpoint(addr):
            word=u64(addr);ptrace(4,pid,addr,(word & ~255)|0xcc);return word
        try:
            wait();mem=os.open('/proc/'+str(pid)+'/mem',os.O_RDWR)
            original=breakpoint(entry);ptrace(7,pid)
            inp=dict(np.load(OUT/(label+'-input.npz')));nzone=len(inp['dm']);ns=inp['X'].shape[1]
            while len(records)<nzone:
                wait();getregs();assert reg.rip==entry+1,hex(reg.rip)
                args=[reg.rdi,reg.rsi,reg.rdx,reg.rcx,reg.r8,reg.r9]+[u64(reg.rsp+8*(i-5)) for i in range(6,38)]
                assert integer(args[3])==ns and integer(args[1])==0
                row=dict(T=scalar(args[6]),rho=scalar(args[8]),X=array(args[5],(ns,)))
                assert integer(args[31])==2, 'Only extended screening is audited'
                saved=[]
                if weak_blend is not None:
                    row['weak_blend_before']=np.array([scalar(addr) for addr in blend_addresses])
                    for addr,value in zip(blend_addresses,weak_blend):
                        saved.append((addr,read(addr,8)))
                        assert os.pwrite(mem,struct.pack('<d',float(value)),addr)==8
                    row['weak_blend_used']=np.array([scalar(addr) for addr in blend_addresses])
                    assert np.array_equal(row['weak_blend_used'],weak_blend)
                row['aux_before']=np.array([scalar(args[k]) for k in [13,14,15,16]])
                columns=[13,14,15,16]
                if replacement is not None and replacement.shape==(nzone,7):
                    columns=[10,11,12,13,14,15,16]
                    row['moments_before']=np.array([scalar(args[k]) for k in [10,11,12]])
                if replacement is not None:
                    assert replacement.shape==(nzone,len(columns))
                    wanted=np.r_[inp['lnd'][len(records)],inp['lnT'][len(records)]]
                    assert max(abs(np.log([row['rho'],row['T']])-wanted))<1e-12
                    assert max(abs(row['X']-inp['X'][len(records)]))<1e-12
                    for k,value in zip(columns,replacement[len(records)]):
                        saved.append((args[k],read(args[k],8)))
                        assert os.pwrite(mem,struct.pack('<d',float(value)),args[k])==8
                row['aux_used']=np.array([scalar(args[k]) for k in [13,14,15,16]])
                if len(columns)==7:
                    row['moments_used']=np.array([scalar(args[k]) for k in [10,11,12]])
                ret=u64(reg.rsp)
                ptrace(4,pid,entry,original);rewind(entry);ptrace(9,pid);wait()
                retword=breakpoint(ret);ptrace(7,pid);wait();getregs();assert reg.rip==ret+1
                assert integer(args[37])==0
                if capture_weak_table and not records:
                    table={name:table_array(name,rank,dtype) for name,rank,dtype in [
                        ('weak_reactions_data',5,'<f4'),('weak_lhs_nuclide_id',1,'<i4'),
                        ('weak_rhs_nuclide_id',1,'<i4'),('weak_reaclib_id',1,'<i4'),
                        ('weak_lowt_rate',1,'<f8')]}
                    for name,size in [('weak_reaction_t9s',12),('weak_reaction_lyerhos',11)]:
                        table[name]=np.frombuffer(read(address('__rates_def_mic_MOD_'+name),size*4),dtype='<f4').copy()
                    table['weak_bicubic']=integer(address('__rates_def_mic_MOD_weak_bicubic'))
                    assert table['weak_reactions_data'].shape[:4]==(4,12,11,3)
                    np.savez_compressed(OUT/(label+'-weak-table.npz'),**table)
                row.update(heat=scalar(args[23]),heat_rho=scalar(args[24]),heat_T=scalar(args[25]),
                    heat_X=array(args[26],(ns,)),dxdt=array(args[27],(ns,)),
                    dxdt_rho=array(args[28],(ns,)),dxdt_T=array(args[29],(ns,)),
                    jacobian=array(args[30],(ns,ns)),neutrino=scalar(args[34]))
                assert all(np.all(np.isfinite(a)) for a in row.values())
                for addr,value in saved:
                    assert os.pwrite(mem,value,addr)==8
                    assert read(addr,8)==value
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
        numerical_memory_written=replacement is not None or weak_blend is not None,
        weak_blend_inputs=weak_blend,
        weak_table_read_only_capture=capture_weak_table,
        argument_indices_zero_based=columns,
        argument_meanings=(['abar','zbar','z2bar'] if len(columns)==7 else [])+['free_electron_abundance','eta','deta_dlnT','deta_dlnrho'],
        inputs_restored_after_each_call=True,returned_results_overwritten=False,
        replacement_sha256=sha(OUT/(label+'-replacement.npy')) if replacement is not None else None,
        independent_profile_control_pending=True))
    print('TRACE',label,len(records),'native calls',flush=True)
