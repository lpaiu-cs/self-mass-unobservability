"""Read the actual gas-EOS excitation state before connecting atomic opacity."""
from pathlib import Path
import argparse
import ctypes
import json
import shutil
import signal
import time
import numpy as np
import def_photon_actual_opacity as a

OUT=a.OUT.parent/'def-photon-eos-populations'
SOURCE=Path('/home/lpaiu/work/direct-eos-gr33/molecular-spectral/source/src')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=['mod_excitation_block.f90','excitation_sum.f90','qstar_calc.f90','qryd_calc.f90',
        'mod_free_eos_constants.f90','mod_ionization_data.f90','qmhd_calc.f90']
    for name in files:shutil.copyfile(SOURCE/name,OUT/name)
    a.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Use the actual same-EOS excitation/occupation state and atomic-source partitions to resolve the microscopic population boundary before changing opacity. Do not fit line strengths or relabel an injected population as a consistent LTE EOS.',
        decision='Identify which physical levels and partition corrections exist in the frozen EOS; reconstruct them and connect compatible microscopic absorption. If incompatible state sets require an EOS extension, state it rather than silently adding opacity-only excited states.',
        budget=dict(native_EOS_calls=4,seconds=60,native_atomic_runs=0,coupled_runs=0,new_stellar_steps=0),
        gates=dict(history_bitwise=True,inventory_relative=1e-10,saved_state_relative=1e-12),
        bindings={str(p):a.digest(p) for p in [Path(__file__),a.previous.OUT/'inventory.npz',a.OUT/'profile-integral/populations.json']+[OUT/n for n in files]}))


def capture():
    assert not (OUT/'eos-state.npz').exists();signal.alarm(60);begin=time.monotonic()
    prior=a.previous;old=prior.old;d,p=prior.base.inputs();gas=old.previous.matter.old.GasEOS()
    gas.inventory_lib=gas.gas_lib;native=gas.gas_lib.ionization_inventory
    def call(mode,value,t,eps,out,info):
        raw=np.full(24,np.nan);native(mode,value,t,eps,raw,info)
        out[:]=raw[:22];out[20]=0.;gas.molecules=raw[22:].copy()
    gas.call=call;lib=gas.gas_lib
    def arr(name,n,dtype=ctypes.c_double):
        return np.ctypeslib.as_array((dtype*n).in_dll(lib,'__mod_excitation_block_MOD_'+name)).copy()
    T=float(np.load(old.OUT/'bank.npz')['T']);snaps=[];checks=[]
    for temperature in [T,T*np.exp(-1e-4),T*np.exp(1e-4),T]:
        snap=prior.inventory_reader.InventoryEOS.snapshot(gas,float(d['lnd'][0]),np.log(temperature),d['X'][0])
        checks.append(prior.inventory_reader.check(snap,d['X'][0],float(d['lnd'][0]),gas))
        count=int(arr('extrace_count',1,ctypes.c_int)[0])
        snap.update(T=temperature,ids=arr('extrace_ids',636,ctypes.c_int).reshape(318,2)[:count],
            value=arr('extrace_value',318*6).reshape(318,6)[:count],scale=arr('extrace_scale',318)[:count],
            x=arr('x',5),qstar=arr('qstar',290).reshape(29,10)[:2,:3],
            flags=np.array([int(arr(n,1,ctypes.c_int)[0]) for n in ['extrace_mode','extrace_nr','ifpi_fit_old_excitation','ifpl_logical_old_excitation','ifmhd_logical_old_excitation']]))
        assert all(np.all(np.isfinite(v)) for v in snap.values());snaps.append(snap)
    assert all(np.array_equal(snaps[0][key],snaps[-1][key]) for key in snaps[0])
    assert all(x['inventory_error']<1e-10 and x['charge_error']<1e-10 for x in checks)
    saved=np.load(prior.OUT/'inventory.npz');assert np.allclose(snaps[0]['eos'],saved['full_EOS'],rtol=1e-12,atol=0)
    np.savez_compressed(OUT/'eos-state.npz',**{k:np.asarray([s[k] for s in snaps]) for k in snaps[0]})
    info=dict(classification='Counterexample candidate',passed=True,seconds=time.monotonic()-begin,
        checks=checks,history_bitwise=True,saved_state_matched=True,ids=snaps[0]['ids'].tolist(),
        flags=snaps[0]['flags'].tolist(),excitation_log_partition=(snaps[0]['value'][:,3]/snaps[0]['scale']).tolist())
    a.write(OUT/'capture.json',info);print(info,flush=True);signal.alarm(0)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','capture']);globals()[p.parse_args().action]()
