def station():
    assert not (OUT/'station.npz').exists();signal.alarm(10);begin=time.monotonic()
    prior=a.p.a.previous;gas_module=prior.old.previous.matter.old
    gas=object.__new__(gas_module.GasEOS);init=gas_module.GasEOS.__init__
    FunctionType(init.__code__,dict(init.__globals__,BRIDGE=a.CACHE/'gas.so'),closure=init.__closure__)(gas)
    gas.inventory_lib=gas.gas_lib;native=gas.gas_lib.ionization_inventory
    def call(mode,value,t,eps,out,info):
        raw=np.full(24,np.nan);native(mode,value,t,eps,raw,info)
        out[:]=raw[:22];out[20]=0.;gas.molecules=raw[22:].copy()
    gas.call=call
    d,_=prior.base.inputs();r=float(d['lnd'][0]);X=d['X'][0]
    T=float(np.load(a.p.OUT/'eos-state.npz')['T'][0])
    snap=prior.inventory_reader.InventoryEOS.snapshot(gas,r,float(np.log(T)),X)
    saved=np.load(a.OUT/'eos-state.npz')
    assert all(np.array_equal(saved[k][0],snap[k]) for k in snap)
    values={}
    for key,n in [('dv',316),('zero',24),('ce',316),('binding',316),('plop',318),('tc2',1),('flags',3),('active',24)]:
        dtype=ctypes.c_int if key in ['flags','active'] else ctypes.c_double
        values[key]=np.ctypeslib.as_array((dtype*n).in_dll(gas.gas_lib,'__mod_nuvar_MOD_station_'+key)).copy()
    values['T']=np.array(T);values['rho']=np.array(np.exp(r))
    values['ground_logw']=np.ctypeslib.as_array((ctypes.c_double*318).in_dll(gas.gas_lib,'__mod_excitation_MOD_shared_ground_logw')).copy()
    np.savez_compressed(OUT/'station.npz',**values)
    write(OUT/'station.json',dict(classification='Counterexample candidate',same_state_bitwise=True,native_EOS_calls=1,seconds=time.monotonic()-begin))
    signal.alarm(0);print('STATION',time.monotonic()-begin,flush=True)
