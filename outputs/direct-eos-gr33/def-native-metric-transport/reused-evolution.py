def evolve(self,seconds,steps,label):
    assert not (OUT/f'{label}.npz').exists()
    start=time.monotonic();gamma=1-1/np.sqrt(np.longdouble(2))
    times=np.longdouble(seconds/self.bg.tc)*(np.arange(steps+1,dtype=np.longdouble)/steps)**2
    state=tuple(np.zeros(n,np.longdouble) for n in [self.size,self.size,self.n,self.n])
    self.history_t=[0.];self.history_e=[0.];self.history_d=[0.]
    self.history_z=[0.];self.history_zd=[0.];self.max_closure_error=0.;self.max_closure_iterations=0
    self.max_error=self.max_heat_error=0.
    rows=[];temp=[];vel=[];scalar=[];surface=[];weights=self.dm/self.dm.sum()
    p=self.bg.sample(self.native);b=1-2*p['m']/self.native
    speed=old.C/100*self.native/(p['N']*np.sqrt(b))
    for j,t in enumerate(times):
        if j:
            dt=t-times[j-1];stage=self.stage(gamma*dt)
            first=stage(state,times[j-1]+gamma*dt)
            base=tuple(a+(1-gamma)/gamma*(c-a) for a,c in zip(state,first));state=stage(base,t)
            self.history_t.append(float(t));self.history_e.append(float(state[2][-1]));self.history_d.append(float(state[3][-1]))
            self.history_z.append(float((self.surfaceV[0]@state[0])[0]));self.history_zd.append(float((self.surfaceV[0]@state[1])[0]))
        q,v,e,d=state;T=self.Tq@q+t*self.TE0+self.TE@e
        velocity=speed*(self.nativeV[0]@v);s=self.nativeV[1]@q
        surf=[float((a@z)[0]) for z in [q,v] for a in self.surfaceV]
        row=dict(t=float(t*self.bg.tc),maximum_delta_lnT=float(max(abs(T))),base_cell_delta_lnT=float(T[self.core_count-1]),
            velocity_RMS_m_s=float(np.sqrt(weights@velocity**2)),maximum_velocity_m_s=float(max(abs(velocity))),
            scalar_RMS=float(np.sqrt(weights@s**2)),surface_displacement=surf[0],surface_scalar=surf[1],
            surface_velocity_coordinate=surf[2],outgoing_luminosity_relative=float(d[-1]/self.f0[-1]))
        self.closure(state,t)
        row.update(surface_lapse=float(self.last_lapse_surface),maximum_lapse=float(max(abs(self.last_lapse_native))))
        rows.append(row);temp.append(T);vel.append(velocity);scalar.append(s);surface.append(surf)
        assert row['maximum_delta_lnT']<.05,('Tangent window exceeded',row)
    dmgeom=self.native**2*b*p['v']*s-4*np.pi*self.native**3*np.exp(-8*p['phi']**2)*(p['e']+p['p'])*(self.nativeV[0]@q)+t*self.J0+self.Jmap@e
    np.savez_compressed(OUT/f'{label}.npz',temperature=temp,velocity=vel,scalar=scalar,surface=surface,
        q=q,v=v,e=e,d=d,f0=self.f0,radius=self.native,edges=self.edges,grid=self.grid,indices=self.indices,
        Eulerian_mass_geom_increment_cm=dmgeom*self.bg.R,emission_times=self.history_t,
        emission_energy_deviation=self.history_e,emission_flux_deviation=self.history_d)
    result=dict(classification='Counterexample candidate',steps=steps,degree=4,graded=True,history=rows,
        seconds=time.monotonic()-start,setup_seconds=self.setup_seconds,max_linear_residual=self.max_error,
        max_local_deviation_heat_identity=self.max_heat_error,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        max_closure_relative=self.max_closure_error,max_closure_iterations=self.max_closure_iterations,
        moving_surface_solved=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/f'{label}.json',result);print('BALANCED',label,result['seconds'],rows[-1],flush=True)
    return result
