def solve(steps,bank,neutrino=False,outer=2,label=None):
    start=time.monotonic();ray=dict(np.load(task.coupled.OUT/'fine-rays.npz'));radiation=task.Radiation(ray,bank,neutrino)
    bg=task.coupled.Background(radiation,outer);fn,_=task.coupled.reactive.symbolic();K,D,base=task.coupled.assemble(bg,fn)
    dt=1/steps;step=Step(K,D,dt);y=np.zeros(K.shape[0]);v=y.copy();history=[]
    native=bg.native['radius_cm'][::-1]/(100*radiation.geometry.R);weights=bg.native['dm'][::-1];weights/=weights.sum()
    r=bg.grid;N,a=radiation.geometry.metric(r)
    native_faces=bg.native['faces_cm']/(100*radiation.geometry.R)
    width=538000/radiation.geometry.R
    masks=[abs(native-native_faces[i])<=width for i in [4123,1723]]
    assert all(np.any(mask) for mask in masks)

    def output(t):
        q=(y-radiation.heat.lift(t,bg.nodes)).reshape(-1,4);p=(v-radiation.heat.lift(t,bg.nodes,True)).reshape(-1,4)
        speed=np.interp(native,r,a/N*r*p[:,0]*task.h.gr.C);scalar=np.interp(native,r,q[:,2]-r*q[:,0]*bg.nodes['v'])
        return dict(tau=t,velocity_mass_RMS_m_s=float(np.sqrt(weights@(speed*speed))),scalar_mass_RMS=float(np.sqrt(weights@(scalar*scalar))),**{name:float(np.sqrt((weights[mask]@speed[mask]**2)/weights[mask].sum())) for name,mask in zip(['old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s'],masks)}),q,p
    forcing=lambda t:base(t)+K@radiation.heat.lift(t,bg.nodes)
    history.append(output(0.)[0]);began=time.monotonic()
    for j in range(steps):
        y,v=step.advance(y,v,j*dt,forcing);history.append(output((j+1)*dt)[0])
    step_seconds=time.monotonic()-began;_,q,p=output(1.);flux,energy=radiation.heat.faces(1.)
    balance=float(abs(np.sum(-np.diff(energy),dtype=np.longdouble))/max(abs(energy).max(),1e-100))
    row=dict(classification='Counterexample candidate',steps=steps,neutrino=neutrino,outer=outer,linear_residual=step.error,
        seconds=time.monotonic()-start,step_seconds=step_seconds,heat_telescoping=balance,history=history)
    if label:
        np.savez_compressed(OUT/(label+'.npz'),grid=bg.grid,response=q,velocity=p,heat_luminosity_erg_s=flux,heat_cumulative_energy_erg=energy)
        task.core.ex.write(OUT/(label+'.json'),row)
    print('RADAU',label,steps,'SECONDS',row['seconds'],'END',history[-1],flush=True)
    return row
