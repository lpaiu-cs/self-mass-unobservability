def run(self,steps,label,limit=None,restart=None):
    assert not (OUT/f'{label}.npz').exists();start=time.monotonic();h=self.t[-1]/steps;gamma=old.prior.GAMMA
    lu=splu(sparse.eye(self.n*self.q,format='csc')-gamma*h*self.A)
    x=np.zeros_like(self.I[0]);g=np.zeros((self.n,2));ledger=np.zeros(2);escape=np.zeros(3);impulse=np.zeros(self.n);ports=np.zeros((2,2));port_history=[];photon_history=[];gas_history=[];transfer_history=[]
    error=0.;species_error=0.;records=[];times=[];begin=0;moment=0.
    def record(t):
        k=int(np.argmin(abs(self.t-t))); pmap=self.point(k)['pressure_map']
        assert abs(self.t[k]-t)<1e-18 or t==count*h,'Only canonical pressure snapshots'
        # Noncanonical pilot endpoint gets an interpolated pressure map.
        if abs(self.t[k]-t)>1e-18:
            j=max(0,min(np.searchsorted(self.t,t)-1,15));f=(t-self.t[j])/(self.t[j+1]-self.t[j]);pmap=(1-f)*self.point(j)['pressure_map']+f*self.point(j+1)['pressure_map']
        pressure=np.einsum('nj,nj->n',pmap,g)*self.volume
        port_history.append(ports.copy()*AMPLITUDE);photon_history.append(x*self.scale*AMPLITUDE);gas_history.append(g.copy()*AMPLITUDE);transfer_history.append(g*np.stack([self.eu,self.nu],axis=-1)*AMPLITUDE)
        times.append(t);records.append(np.stack([np.sum(x*self.Eweight,axis=(1,2)),g[:,0]*self.eu,g[:,1]*self.nu,impulse,np.sum(abs(x)*self.Eweight,axis=(1,2)),np.sum(x*self.Eweight*self.model.bulk.mu2[None,:,None],axis=(1,2)),pressure])*AMPLITUDE)
    count=steps if limit is None else limit
    if restart is None:record(0.)
    else:
        z=np.load(OUT/f'{restart}.npz');previous=json.loads((OUT/f'{restart}.json').read_text());begin=previous['completed_steps']
        x=z['delta_packet_scaled_occupation']/(self.scale*AMPLITUDE);g=z['delta_material']/AMPLITUDE
        ledger=z['ledger']/AMPLITUDE;escape=z['escape']/AMPLITUDE;impulse=z['moments'][-1,3]/AMPLITUDE
        ports=z['radial_ports'][-1]/AMPLITUDE;port_history=list(z['radial_ports']);photon_history=list(z['photon_history_scaled_occupation']);gas_history=list(z['material_history']);transfer_history=list(z['collision_transfer']);self.max_residual=previous['linear_relative'];self.max_iterations=previous['max_Krylov_iterations'];self.angular_times=list(z['accepted_angular_times']);self.angular=list(z['accepted_angular_luminosity']);times=z['t'].tolist();records=list(z['moments']);error=previous['energy_balance_relative'];species_error=previous['species_balance_relative']
    def stream(xx):return (self.A@xx.reshape(self.n*self.q,self.nf)).reshape(xx.shape)
    for k in range(begin,count):
        t=k*h;c=self.local(t+gamma*h);inverse=self.inverse(c,gamma*h)
        def L(v):
            xx,gg=self.unpack(v);p,q,*_=self.collision(c,xx,gg);return self.pack(stream(xx)+p,q)
        def mat(v):return v-gamma*h*L(v)
        def pre(v):
            xx,gg=self.unpack(v);xx=lu.solve(xx.reshape(self.n*self.q,self.nf)).reshape(xx.shape)
            return inverse(self.pack(xx,gg))
        op=LinearOperator((self.size+2*self.n,)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
        def solve(rhs,guess):
            iterations=[];sol,info=gmres(op,rhs,x0=guess,M=P,rtol=1e-12,atol=0.,restart=20,maxiter=5,callback=iterations.append,callback_type='pr_norm')
            err=float(np.linalg.norm(mat(sol)-rhs)/max(np.linalg.norm(rhs),1e-300));self.max_residual=max(self.max_residual,err);self.max_iterations=max(self.max_iterations,len(iterations))
            assert info==0 and err<1e-12,('Monolithic stage',info,err,len(iterations))
            return self.unpack(sol)
        qg=self.gas(c['q'],c['qb'],c['qe'])
        s1,l1,e1=self.source(t+gamma*h);s1=s1/(self.scale*AMPLITUDE);l1=l1/AMPLITUDE
        y,gy=solve(self.pack(x+gamma*h*(s1+c['q']),g+gamma*h*qg),self.pack(x,g))
        p1,g1,es1,_=self.collision(c,y,gy,True);f1=stream(y)+p1+s1
        c=self.local(t+h);inverse=self.inverse(c,gamma*h);qg=self.gas(c["q"],c["qb"],c["qe"])
        s2,l2,e2=self.source(t+h);s2=s2/(self.scale*AMPLITUDE);l2=l2/AMPLITUDE
        z,gz=solve(self.pack(x+(1-gamma)*h*f1+gamma*h*(s2+c['q']),g+(1-gamma)*h*g1+gamma*h*qg),self.pack(y,gy))
        p2,g2,es2,_=self.collision(c,z,gz,True)
        ports+=h*((1-gamma)*self.boundary_ports(t+gamma*h,y)+gamma*self.boundary_ports(t+h,z))
        esc=h*((1-gamma)*es1+gamma*es2);escape+=esc.sum(1)+h*((1-gamma)*l1[:3]+gamma*l2[:3])
        ledger+=h*((1-gamma)*(self.port(y*self.scale)+l1[3:])+gamma*(self.port(z*self.scale)+l2[3:]))-esc[:2].sum(1)
        impulse-=h*np.sum(((1-gamma)*p1+gamma*p2)*self.Eweight*self.mu[None,:,None],axis=(1,2))+esc[2]
        x,g=z,gz;moment=max(moment,e1,e2)
        photon=self.moments(x*self.scale);total=np.array([photon[0]-np.sum(g[:,1]*self.nu),photon[1]+np.sum(g[:,0]*self.eu)])
        norm=np.array([max(np.sum(abs(x)*self.Nweight),np.sum(abs(g[:,1])*self.nu),1.),max(np.sum(abs(x)*self.Eweight),np.sum(abs(g[:,0])*self.eu),1.)])
        defect=abs(total-ledger)/np.maximum(norm,abs(ledger));species_error=max(species_error,float(defect[0]));error=max(error,float(defect[1]))
        if (k+1)%max(1,steps//16)==0 or k+1==count:
            record((k+1)*h)
            path=OUT/f'{label}-checkpoint.npz';tmp=path.with_suffix('.tmp')
            with tmp.open('wb') as handle:np.savez_compressed(handle,x=x,g=g,ledger=ledger,escape=escape,impulse=impulse,completed_steps=k+1,t=times,moments=records,ports=ports,photon_history=photon_history,gas_history=gas_history,port_history=port_history,transfer_history=transfer_history,accepted_angular_times=self.angular_times,accepted_angular_luminosity=self.angular)
            tmp.replace(path);write(OUT/f'{label}-progress.json',dict(completed_steps=k+1,planned_steps=steps,seconds=time.monotonic()-start))
    solve_seconds=time.monotonic()-start
    np.savez_compressed(OUT/f'{label}.npz',t=times,moments=records,delta_packet_scaled_occupation=x*self.scale*AMPLITUDE,delta_material=g*AMPLITUDE,material_energy_units=self.eu,material_neutral_units=self.nu,ledger=ledger*AMPLITUDE,escape=escape*AMPLITUDE,radius_E=self.r,radial_ports=port_history,photon_history_scaled_occupation=photon_history,material_history=gas_history,collision_transfer=transfer_history)
    row=dict(classification='Counterexample candidate',reference=self.reference,steps=steps,completed_steps=count,new_steps=count-begin,seconds=time.monotonic()-start,operator_point_seconds=self.point_seconds,operator_points=self.point_count,stepping_seconds=solve_seconds-self.point_seconds,energy_balance_relative=error,species_balance_relative=species_error,linear_relative=self.max_residual,max_Krylov_iterations=self.max_iterations,frequency_moment_relative=moment,endpoint_photon_reference_energy_erg=float(np.sum(x*self.Eweight)*AMPLITUDE),endpoint_material_reference_energy_erg=float(np.sum(g[:,0]*self.eu)*AMPLITUDE),endpoint_material_energy_L1_erg=float(np.sum(abs(g[:,0])*self.eu)*AMPLITUDE),additional_material_motion_evolved=False,full_GR_feedback=False,final_charge_solved=False)
    row['passed']=max(error,species_error)<1e-8 and self.max_residual<1e-12
    write(OUT/f'{label}.json',row);print(json.dumps(row),flush=True);return row
