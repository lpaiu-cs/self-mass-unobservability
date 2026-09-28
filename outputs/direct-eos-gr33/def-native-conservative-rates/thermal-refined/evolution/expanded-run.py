def run_capture(self,steps,count=None,resume=False):
    started=time.monotonic();b=self.bulk;f=self.flow;m=self.m;h=flow.old.END/steps
    checkpoint=OUT/f'checkpoint-{steps}.npz';sidecar=OUT/f'history-{steps}.npz'
    U=f.initial.copy();I=np.stack([self.initial_I,np.zeros_like(self.initial_I)]);xb=b.initial/b.scale
    u=b.u0.copy();eta=np.zeros(b.n);theta=b.eos.recover(u,eta,np.zeros(b.n));ledger=np.zeros(6);discard=np.zeros(4);ports=np.zeros(2);begin=0;owner_port=0.;steps_local=0
    if resume:
        z=np.load(checkpoint);U=z['U'];I=z['I'];xb=z['xb'];u=z['u'];theta=z['theta'];eta=z['eta'];begin=int(z['completed'])
        self.Pi=z['Pi'];self.h=z['h'];self.j=z['j'];self.mass=self.mass0-self.h[1:]+self.h[:-1]
        for key in SCALARS:setattr(self,key,float(z['scalar_'+key]))
        ledger=z['ledger'];discard=z['discard'];ports=z['ports'];owner_port=float(z['owner_port']);self.set_material(begin*h)
        saved=np.load(sidecar);self.rows=list(saved['moments']);self.times=list(saved['t']);self.port_rows=list(saved['ports']);self.discard_rows=list(saved['discard']);self.escape_rows=list(saved['escape'])
    else:
        self.rows=[self.moment(U,I,xb,u,theta,eta)];self.times=[0.];self.port_rows=[ports.copy()];self.discard_rows=[discard.copy()];self.escape_rows=[0.]
    def save(completed):
        # Only compact moments persist after success; checkpoint permits
        # restart without another full identical-path replay.
        for path,values in [(sidecar,dict(t=self.times,moments=self.rows,ports=self.port_rows,discard=self.discard_rows,escape=self.escape_rows)),
            (checkpoint,dict(U=U,I=I,xb=xb,u=u,theta=theta,eta=eta,Pi=self.Pi,h=self.h,j=self.j,ledger=ledger,discard=discard,ports=ports,owner_port=owner_port,completed=completed,
                **{'scalar_'+key:getattr(self,key,0.) for key in SCALARS}))]:
            tmp=path.with_suffix('.tmp')
            with tmp.open('wb') as stream:np.savez_compressed(stream,**values)
            tmp.replace(path)
    last=steps if count is None else count
    for k in range(begin,last):
        U,I,u,theta,eta,ll,dd,ss=self.local_joint(U,I,u,theta,eta,k*h,h/2);ledger+=ll;discard+=dd;steps_local+=ss
        xb,I,u,theta,eta,port,it,err=self.radiate(xb,I,u,theta,eta,h)
        actual=h*self.port_energy(xb,I);assert abs((actual[0]-actual[1])-port[0])<1e-12*max(abs(port[0]),1.)
        ports+=actual;owner_port+=port[0]
        U,I,u,theta,eta,ll,dd,ss=self.local_joint(U,I,u,theta,eta,k*h+h/2,h/2);ledger+=ll;discard+=dd;steps_local+=ss
        self.rows.append(self.moment(U,I,xb,u,theta,eta));self.times.append((k+1)*h);self.port_rows.append(ports.copy());self.discard_rows.append(discard.copy());self.escape_rows.append(float(ledger[5]+getattr(self,'deep_escape',0.)))
        if (k+1)%16==0:save(k+1);print(json.dumps(dict(steps=steps,completed=k+1,seconds=time.monotonic()-started)),flush=True)
    save(last)
    source=self.template.copy();arr=np.asarray(self.rows)
    for i,key in enumerate(KEYS):source[key]=arr[:,i]
    source.update(t=np.array(self.times),inner_cumulative_energy_erg=np.array(self.port_rows)[:,0],outer_cumulative_energy_erg=np.array(self.port_rows)[:,1],
        metric_stress_erg=arr[:,0]*LD(self.cx)*LD(C)**2+arr[:,3]+arr[:,5]-arr[:,6])
    # Old diagnostic fields have a different time grid; never leave them
    # looking like the newly captured exact accepted-stage histories.
    for key in ['Killing_cell_energy_erg','inner_luminosity','outer_luminosity','port_mismatch_erg','port_error_bound_erg','spectral_escape_erg','luminosity_per_mu','velocity']:source.pop(key,None)
    energy=(arr[:,0]*LD(self.cx)*LD(C)**2+arr[:,1]+arr[:,5])*source['a']
    disc=np.array(self.discard_rows);loss=(disc[:,2]+self.m.a0*self.cx*disc[:,0])*self.gas_scale+np.array(self.escape_rows)
    defect=np.sum(energy,axis=1,dtype=LD)+loss-np.diff(-np.array(self.port_rows),axis=1)[:,0]
    scale=max(np.max(abs(np.array(self.port_rows))),1.);balance=float(np.max(abs(defect))/scale)
    reference=np.load(flow.OUT/f'coupled-{steps}.npz');idx=np.argmin(abs(reference['snapshot_t']-self.times[-1]));assert abs(reference['snapshot_t'][idx]-self.times[-1])<1e-18
    state_errors={}
    for key,value,base in [('U',U,f.initial),('I',I,np.stack([self.initial_I,np.zeros_like(self.initial_I)])),('bulk_I',xb*b.scale,b.initial),('u',u,b.u0),('h',self.h,np.zeros_like(self.h)),('Pi',self.Pi,np.zeros_like(self.Pi))]:
        old=reference['snapshot_'+key][idx];state_errors[key]=float(np.sum(abs(value-old),dtype=LD)/max(np.sum(abs(old-base),dtype=LD),LD('1e-100')))
    matched=[];legacy=np.load(gr.base.OUT/f'source-{steps}.npz')
    for oldt in legacy['t']:
        if oldt>self.times[-1]+1e-18:continue
        ii=int(np.argmin(abs(source['t']-oldt)));jj=int(np.argmin(abs(legacy['t']-oldt)))
        for key in KEYS:
            denom=max(np.max(np.sum(abs(legacy[key]),axis=1)),LD('1e-100'))
            matched.append(float(np.sum(abs(source[key][ii]-legacy[key][jj]),dtype=LD)/denom))
    row=dict(classification='Counterexample candidate',steps=steps,completed=last,seconds=time.monotonic()-started,local_steps=steps_local,
        conserved_Killing_balance=balance,prior_state_change=state_errors,prior_source_change=max(matched,default=0.),
        exact_owner_port_relative=abs(ports[0]-ports[1]-owner_port)/scale,
        complete=last==steps,passed=bool(balance<1e-8 and abs(ports[0]-ports[1]-owner_port)/scale<1e-10),
        full_source_error_enclosed=False,final_charge_solved=False)
    np.savez_compressed(OUT/f'source-{steps}.npz',**source)
    write(OUT/f'{"result" if last==steps else "pilot"}-{steps}.json',row);print(json.dumps(row),flush=True);return row
