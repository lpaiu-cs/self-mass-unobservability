emission={}
def inputs(n):
    if n in emission:return emission[n]
    with np.load(EV/'accepted-ports-128.npz') as z:lum=z['angular_luminosity'].astype(LD);bh=float(z['h']);h=bh*128/n
    with np.load(PHOTON/f'steps-{n}-reference-128.npz') as z:
        t=z['accepted_angular_times'];delta=z['accepted_angular_luminosity'].astype(LD);ports=z['radial_ports'];canonical=z['t']
    gamma=1-1/np.sqrt(2);assert lum.shape==(128,4) and delta.shape==(2*n,4)
    assert np.max(abs(t-(h*(np.arange(n)[:,None]+[gamma,1.])).ravel()))<1e-18
    assert abs(h*n-m.T)<1e-18 and np.min(lum)>=0
    w=np.tile([1-gamma,gamma],n).astype(LD)*LD(h);packets=w[:,None]*delta
    q=np.vstack([np.zeros(4,LD),np.cumsum(LD(bh)*lum,axis=0)])
    r=np.vstack([np.zeros(4,LD),np.cumsum(LD(bh)*q[:-1]+LD(bh)**2/2*lum,axis=0)])
    dq=np.vstack([np.zeros(4,LD),np.cumsum(packets,axis=0)])
    dtq=np.vstack([np.zeros(4,LD),np.cumsum(packets*t[:,None],axis=0)])
    aw=np.arange(1,8,2,dtype=LD)/32
    with np.load(EV/'source-128.npz') as z:baseport=z['outer_cumulative_energy_erg'];bt=z['t']
    err0=float(np.max(abs(q@aw-baseport))/max(abs(baseport[-1]),1.))
    ids=np.array([int(np.argmin(abs(bt-v))) for v in canonical]);assert np.max(abs(bt[ids]-canonical))<1e-18
    err1=float(np.max(abs((dq[2*np.rint(canonical/h).astype(int)]@aw)-ports[:,1,1]))/max(np.sum(abs(packets)@aw),1.))
    assert max(err0,err1)<1e-12,(err0,err1)
    def primitives(at):
        x=np.clip(np.asarray(at),0,m.T);j=np.clip((x/bh).astype(int),0,127);s=(x-j*bh).astype(LD)[...,None]
        hb=q[j]+s*lum[j];hhb=r[j]+s*q[j]+s*s/2*lum[j]
        k=np.searchsorted(t,x,side='right');hd=dq[k];hhd=x[...,None]*dq[k]-dtq[k]
        return hb,hhb,hd,hhd
    full=primitives(m.T);assert np.max(abs(full[0]-q[-1]))/np.max(q[-1])<1e-14
    emitted=float((q[-1]+dq[-1])@aw);absolute=float(q[-1]@aw+np.sum(abs(packets)@aw))
    row=dict(steps=n,background_port_relative=err0,response_port_relative=err1,emitted_energy_erg=emitted,
             absolute_emission_energy_erg=absolute,signed_increment_energy_erg=float(dq[-1]@aw))
    emission[n]=SimpleNamespace(primitives=primitives,row=row)
    return emission[n]

def evaluate(n,a,r):
    begin=time.monotonic();source=inputs(n);k=m.kernel(a,r);mass=[];stress=[];arrival=[]
    for t in m.t:
        hb,hhb,hd,hhd=source.primitives(np.maximum(t-k['delay'],0.));ids=np.arange(len(hb));bins=k['node_bins']
        mass.append([-G/C**3*float(np.sum(k['weights']*k['mass']*v[ids,bins],dtype=LD)) for v in [hhb,hhd]])
        stress.append([-G/(2*C**4)*float(np.sum(k['weights']*k['stress']*v[ids,bins],dtype=LD)) for v in [hb,hd]])
        hb,_,hd,_=source.primitives(np.maximum(t-k['infinity'],0.));ids=np.arange(len(hb))
        arrival.append([float(np.sum(k['mw']*k['mu']*v[ids,k['bins']],dtype=LD)) for v in [hb,hd]])
    mass=np.asarray(mass);stress=np.asarray(stress);arrival=np.asarray(arrival)
    d=dict(t=m.t,mass_U=mass,stress_U=stress,normalized_exterior=-(mass+stress).sum(1)/m.M,
           arrived_energy_erg=arrival.sum(1),arrived_parts_erg=arrival,normalized_exterior_parts=-(mass+stress)/m.M)
    row=dict(source.row,angular=a,radial=r,seconds=time.monotonic()-begin,kernel_seconds=k['seconds'],
        delay_inverse_seconds=k['inversion_error'],endpoint_exterior=float(d['normalized_exterior'][-1]),
        endpoint_arrived_energy_erg=float(d['arrived_energy_erg'][-1]),endpoint_increment_exterior=float(d['normalized_exterior_parts'][-1,1]))
    label=f'exterior-{n}-a{a}-r{r}';np.savez_compressed(OUT/f'{label}.npz',**d);write(OUT/f'{label}.json',row)
    return d,row

