"""Recover missing stage flux/source moments by bounded identical-path replay.

No EOS/model/clock change. The existing sparse snapshots cannot recover the
accepted implicit-stage flux integral; record it where the owner returns it.
"""
from pathlib import Path
import json,signal,sys,time
import numpy as np
import def_native_energy_source_reconciliation as prior

flow=prior.flow;gr=prior.gr;OUT=prior.OUT/'stage-history';C=flow.C;LD=np.longdouble;write=flow.write;sha=flow.sha
KEYS=['baryon_g','gas_nonrest_energy_erg','nonrest_trace_erg','nonrest_stress_erg','pressure_volume_erg','photon_energy_erg','photon_radial_pressure_erg']
SCALARS=['mechanical_energy','radiation_work','boundary_species','max_frame','max_density','maximum_inner_iterations','mechanical_steps','deep_escape','join_mass','join_momentum','join_energy','join_neutral']


class Capture(flow.Coupled):
    def __init__(self):
        super().__init__();f=self.flow;m=self.m;b=self.bulk
        self.template=dict(np.load(gr.base.OUT/'source-128.npz'))
        self.volumes=np.r_[b.volume,4*np.pi*m.RJ*m.RJ*m.vol]
        old_y=f.eos.y;f.eos.y=f.eos.y0
        self.ip=f.eos(f.initial[0],f.initial_temperature)[0];f.eos.y=old_y
        self.tau0=(f.initial[2]-(m.a-m.a0)*f.eos.cx*f.initial[0])/m.a
        self.rows=[];self.port_rows=[];self.discard_rows=[];self.escape_rows=[];self.times=[]

    def moment(self,U,I,xb,u,theta,eta):
        f=self.flow;m=self.m;b=self.bulk;old_seed=f.seed.copy();old_y=np.copy(f.eos.y)
        try:
            rho,v,lt,y=f.primitive(U);p=f.eos(rho,lt)[0]
            tau=(U[2]-(m.a-m.a0)*f.eos.cx*U[0])/m.a
            pp=b.eos.gas(theta,eta)[0];dm=-np.diff(self.h);beta=self.velocity()
            de=dm*b.u0+self.mass*(u-b.u0);ktrace=-self.mass*(self.cx*C*C+u)*beta**2/(1+np.sqrt(1-beta**2))
            dpb=(pp-self.f0['p0'])*b.volume;unit=f.eos.rho0*C*C*self.volumes[b.n:];dp=p-self.ip
            baryon=np.r_[dm,(U[0].astype(LD)-f.initial[0])*f.eos.rho0*self.volumes[b.n:]]
            gas=np.r_[de+self.kinetic(),(tau-self.tau0)*unit]
            trace=np.r_[de+ktrace-3*dpb,((tau-self.tau0)-U[1]*v-3*dp)*unit]
            stress=np.r_[de+ktrace-dpb,((tau-self.tau0)-U[1]*v-dp)*unit]
            pressure=np.r_[dpb,dp*unit]
            # Difference first, as in the actual conservation owner.
            ib=xb*b.scale-b.initial;ia=I.sum(0)-self.initial_I;en=b.d['num']*b.d['Einf']
            ph=np.r_[np.einsum('iqf,q,f->i',ib,b.w,en)/b.d['a']**4,np.einsum('iqf,q,f->i',ia,b.w,en)/m.a**4]*self.volumes
            pr=np.r_[np.einsum('iqf,q,f->i',ib,b.w*b.mu2,en)/b.d['a']**4,np.einsum('iqf,q,f->i',ia,b.w*b.mu2,en)/m.a**4]*self.volumes
            return np.array([baryon,gas,trace,stress,pressure,ph,pr],dtype=LD)
        finally:f.seed=old_seed;f.eos.y=old_y

    def port_energy(self,xb,I):
        b=self.bulk;inner=np.where((b.mu>0)[:,None],b.incoming,xb[0]*b.scale)
        outer=np.where((b.mu>0)[:,None],I.sum(0)[-1],0.)
        weights=b.d['num']*b.d['Einf']
        return np.array([4*np.pi*C*b.area[0]*((b.w*b.mu@inner)*weights).sum(),
            4*np.pi*C*self.area[-1]*((b.w*b.mu@outer)*weights).sum()])

    def run_capture(self,steps,count=None,resume=False):
        started=time.monotonic();b=self.bulk;f=self.flow;m=self.m;h=flow.old.END/steps
        checkpoint=OUT/f'checkpoint-{steps}.npz';sidecar=OUT/f'history-{steps}.npz'
        U=f.initial.copy();I=np.stack([self.initial_I,np.zeros_like(self.initial_I)]);xb=b.initial/b.scale
        u=b.u0.copy();theta=np.zeros(b.n);eta=np.zeros(b.n);ledger=np.zeros(6);discard=np.zeros(4);ports=np.zeros(2);begin=0;owner_port=0.;steps_local=0
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
            conserved_Killing_balance=balance,replay_state_relative=state_errors,legacy_source_relative=max(matched,default=0.),
            exact_owner_port_relative=abs(ports[0]-ports[1]-owner_port)/scale,
            complete=last==steps,passed=balance<1e-8 and max(state_errors.values())<1e-8 and max(matched,default=0.)<1e-8,
            full_source_error_enclosed=False,final_charge_solved=False)
        np.savez_compressed(OUT/f'source-{steps}.npz',**source)
        write(OUT/f'{"result" if last==steps else "pilot"}-{steps}.json',row);print(json.dumps(row),flush=True);return row


def prepare():
    assert not OUT.exists();OUT.mkdir()
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='414e242ad',
        claim='Recover exact accepted-stage inner/outer energy integrals and every original-step GR source moment, then apply them to the final charge. Resolve the sparse17-knot port defect without fitting a correction.',
        reason_for_replay='Saved17-knot snapshots and final total port cannot identify separate cumulative inner/outer ports at intermediate steps. The mismatch is mostly a posteriori trapezoid integration of implicit Euler stage rates, not a failure of conserved energy conversion. A bounded identical-path replay is necessary to recover the missing histories.',
        reuse='Same accepted64/128 paths, fixed source equations,EOS banks,531 cells,8angles,152frequencies,3.434ms. No new physical resolution. Store compact moments at each already existing step; preserve old artifacts.',
        changes='Capture exact radiation owner stage ports and stable photon increments; primitive readout restores its seed/species state so it cannot steer the next step. Record actual floor/escape ledgers separately; do not erase lost matter by an energy renormalization.',
        gates=dict(energy=1e-8,replay=1e-8,source=1e-8,charge_time=.02,charge_quadrature=.002),
        budget=dict(pilot_seconds=35,production_seconds=600,readout_seconds=90,CPU_threads=1,memory_GB=3),
        forecast='Original two-path run cost about350s. Pilot two full steps on each clock; compare per-step estimate to both saved full-run timings and require1.5x the larger estimate plus20s to fit600s. No automatic clock/grid/horizon expansion.',
        stop='Halt on conservation, replay, source, forecast or wall-time gate. Each16-step checkpoint supports resume. Do not rerun a successful path or loosen failed limits.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(flow.__file__),prior.OUT/'diagnosis.json',flow.OUT/'coupled-64.npz',flow.OUT/'coupled-128.npz',gr.base.OUT/'source-64.npz',gr.base.OUT/'source-128.npz']}))


def pilot():
    assert not (OUT/'pilot.json').exists();signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(35);rows=[]
    for steps in [64,128]:
        rows.append(Capture().run_capture(steps,2))
        if not rows[-1]['passed']:break
    estimate=sum(r['seconds']/2*(r['steps']-2) for r in rows)
    prior_cost=sum(json.loads((flow.OUT/f'coupled-{s}.json').read_text())['seconds'] for s in [64,128])
    upper=1.5*max(estimate,prior_cost)+20
    result=dict(classification='Counterexample candidate',paths=rows,measured_step_estimate=estimate,prior_full_path_seconds=prior_cost,upper_seconds=upper,
        eligible=len(rows)==2 and all(r['passed'] for r in rows) and upper<600)
    write(OUT/'pilot.json',result);signal.alarm(0);print(json.dumps(result),flush=True)
    if result['eligible']:write(OUT/'execution-plan.json',dict(budget_seconds=600,upper_seconds=upper,bindings={str(p):sha(p) for p in [Path(__file__),OUT/'plan.json',OUT/'pilot.json']}))


def production():
    assert not (OUT/'result.json').exists();p=json.loads((OUT/'execution-plan.json').read_text())
    for file,h in p['bindings'].items():assert sha(file)==h,file
    signal.signal(signal.SIGALRM,flow.old.optical.timeout);signal.alarm(p['budget_seconds']);start=time.monotonic();rows=[]
    try:
        for steps in [64,128]:
            row=Capture().run_capture(steps,resume=True);rows.append(row)
            if not row['passed']:break
        result=dict(classification='Counterexample candidate',passed=len(rows)==2 and all(r['passed'] for r in rows),paths=rows,seconds=time.monotonic()-start,
            exact_stage_port_history_captured=True,applied_to_GR_charge=False,final_charge_solved=False,full_goal_complete=False)
        write(OUT/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:write(OUT/'failure.json',dict(error=repr(exc),paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
