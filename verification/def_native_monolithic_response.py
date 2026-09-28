"""Counterexample candidate: remove transport/collision time partition.

Same linear response, physical mesh, native banks and acceptance gates as
the rejected split run. Krylov solves the full coupled stage equation.
"""
from pathlib import Path
import json
import signal
import sys
import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu,LinearOperator,gmres
import def_native_collision_response as old

OUT=old.OUT/'monolithic';write=old.write;sha=old.sha;AMPLITUDE=old.AMPLITUDE


class Response(old.CoupledResponse):
    def run(self,steps,label,limit=None,restart=None):
        assert not (OUT/f'{label}.npz').exists();start=time.monotonic();h=self.t[-1]/steps;gamma=old.prior.GAMMA
        lu=splu(sparse.eye(self.n*self.q,format='csc')-gamma*h*self.A)
        x=np.zeros_like(self.I[0]);g=np.zeros((self.n,2));ledger=np.zeros(2);escape=np.zeros(3);impulse=np.zeros(self.n)
        error=0.;species_error=0.;records=[];times=[];begin=0;moment=0.
        def record(t):
            k=int(np.argmin(abs(self.t-t))); pmap=self.point(k)['pressure_map']
            assert abs(self.t[k]-t)<1e-18 or t==count*h,'Only canonical pressure snapshots'
            # Noncanonical pilot endpoint gets an interpolated pressure map.
            if abs(self.t[k]-t)>1e-18:
                j=max(0,min(np.searchsorted(self.t,t)-1,15));f=(t-self.t[j])/(self.t[j+1]-self.t[j]);pmap=(1-f)*self.point(j)['pressure_map']+f*self.point(j+1)['pressure_map']
            pressure=np.einsum('nj,nj->n',pmap,g)*self.volume
            times.append(t);records.append(np.stack([np.sum(x*self.Eweight,axis=(1,2)),g[:,0]*self.eu,g[:,1]*self.nu,impulse,np.sum(abs(x)*self.Eweight,axis=(1,2)),np.sum(x*self.Eweight*self.model.bulk.mu2[None,:,None],axis=(1,2)),pressure])*AMPLITUDE)
        count=steps if limit is None else limit
        if restart is None:record(0.)
        else:
            z=np.load(OUT/f'{restart}.npz');previous=json.loads((OUT/f'{restart}.json').read_text());begin=previous['completed_steps']
            x=z['delta_packet_scaled_occupation']/(self.scale*AMPLITUDE);g=z['delta_material']/AMPLITUDE
            ledger=z['ledger']/AMPLITUDE;escape=z['escape']/AMPLITUDE;impulse=z['moments'][-1,3]/AMPLITUDE
            times=z['t'].tolist();records=list(z['moments']);error=previous['energy_balance_relative'];species_error=previous['species_balance_relative']
        def stream(xx):return (self.A@xx.reshape(self.n*self.q,self.nf)).reshape(xx.shape)
        for k in range(begin,count):
            t=k*h;c=self.local(t+h/2);inverse=self.inverse(c,gamma*h)
            def L(v):
                xx,gg=self.unpack(v);p,q,*_=self.collision(c,xx,gg);return self.pack(stream(xx)+p,q)
            def mat(v):return v-gamma*h*L(v)
            def pre(v):
                xx,gg=self.unpack(v);xx=lu.solve(xx.reshape(self.n*self.q,self.nf)).reshape(xx.shape)
                return inverse(self.pack(xx,gg))
            op=LinearOperator((self.size+2*self.n,)*2,mat,dtype=float);P=LinearOperator(op.shape,pre,dtype=float)
            def solve(rhs,guess):
                iterations=[];sol,info=gmres(op,rhs,x0=guess,M=P,rtol=1e-10,atol=0.,restart=20,maxiter=5,callback=iterations.append,callback_type='pr_norm')
                err=float(np.linalg.norm(mat(sol)-rhs)/max(np.linalg.norm(rhs),1e-300));self.max_residual=max(self.max_residual,err);self.max_iterations=max(self.max_iterations,len(iterations))
                assert info==0 and err<1e-10,('Monolithic stage',info,err,len(iterations))
                return self.unpack(sol)
            qg=self.gas(c['q'],c['qb'],c['qe'])
            s1,l1,e1=self.source(t+gamma*h);s1=s1/(self.scale*AMPLITUDE);l1=l1/AMPLITUDE
            y,gy=solve(self.pack(x+gamma*h*(s1+c['q']),g+gamma*h*qg),self.pack(x,g))
            p1,g1,es1,_=self.collision(c,y,gy,True);f1=stream(y)+p1+s1
            s2,l2,e2=self.source(t+h);s2=s2/(self.scale*AMPLITUDE);l2=l2/AMPLITUDE
            z,gz=solve(self.pack(x+(1-gamma)*h*f1+gamma*h*(s2+c['q']),g+(1-gamma)*h*g1+gamma*h*qg),self.pack(y,gy))
            p2,g2,es2,_=self.collision(c,z,gz,True)
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
                with tmp.open('wb') as handle:np.savez_compressed(handle,x=x,g=g,ledger=ledger,escape=escape,impulse=impulse,completed_steps=k+1,t=times,moments=records)
                tmp.replace(path);write(OUT/f'{label}-progress.json',dict(completed_steps=k+1,planned_steps=steps,seconds=time.monotonic()-start))
        solve_seconds=time.monotonic()-start
        np.savez_compressed(OUT/f'{label}.npz',t=times,moments=records,delta_packet_scaled_occupation=x*self.scale*AMPLITUDE,delta_material=g*AMPLITUDE,material_energy_units=self.eu,material_neutral_units=self.nu,ledger=ledger*AMPLITUDE,escape=escape*AMPLITUDE,radius_E=self.r)
        row=dict(classification='Counterexample candidate',reference=self.reference,steps=steps,completed_steps=count,new_steps=count-begin,seconds=time.monotonic()-start,operator_point_seconds=self.point_seconds,operator_points=self.point_count,stepping_seconds=solve_seconds-self.point_seconds,energy_balance_relative=error,species_balance_relative=species_error,linear_relative=self.max_residual,max_Krylov_iterations=self.max_iterations,frequency_moment_relative=moment,endpoint_photon_reference_energy_erg=float(np.sum(x*self.Eweight)*AMPLITUDE),endpoint_material_reference_energy_erg=float(np.sum(g[:,0]*self.eu)*AMPLITUDE),endpoint_material_energy_L1_erg=float(np.sum(abs(g[:,0])*self.eu)*AMPLITUDE),additional_material_motion_evolved=False,full_GR_feedback=False,final_charge_solved=False)
        row['passed']=max(error,species_error)<1e-8 and self.max_residual<1e-10
        write(OUT/f'{label}.json',row);print(json.dumps(row),flush=True);return row


def prepare():
    assert not OUT.exists();OUT.mkdir();prior=json.loads((old.OUT/'result.json').read_text());assert not prior['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='e71ed7eba',
        failure='The completed split evolution conserves energy and H, but64/128 photon energy and momentum impulse differences are3.7903percent and4.6305percent, exceeding the unchanged2percent gate.',
        repair='Solve streaming, actual moving collisions, thermal energy and H in a single SDIRK stage. The split inverse is only a Krylov preconditioner; require the residual of the full equation below1e-10.',
        reuse='All corrected backgrounds, EOS derivative banks,531 cells,8 angles,152 frequencies,64/128 clocks and3.434ms horizon. Old split trajectories remain rejected. Only the time equation changes.',
        coefficients='Freeze the linearized local operator at each full-step midpoint, with both geometric forcing stage times retained. This is a second-order time approximation tested on the unchanged clocks.',
        saved_sources='Save full time histories of photon radial pressure and material pressure, as well as energy, neutral inventory and delivered momentum. The previous split history did not include photon radial pressure and is not a complete GR source export.',
        limits='Additional density, velocity, advected inventory and GR feedback remain open. No final physical charge.',
        budget=dict(pilot_seconds=55,production_seconds=480,CPU_threads=1,memory_GB=3,new_native_bank_calls=0),
        forecast='Measure both actual four-step prefixes; separate17-point operator preparation from time stepping, include three model/output overheads and2x total margin. Do not dispatch if it exceeds480s. Late Krylov cost is unmeasured.',
        gates=dict(time=.02,background_time=.02,energy=1e-8,species=1e-8,linear=1e-10),
        stop='Preserve the old failure and all new prefixes. No finer mesh or clocks, longer horizon, or weakened gates. Stop on any failed path or480s; no automatic enlargement.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.__file__),old.OUT/'result.json',old.OUT/'bank-result.json',old.OUT/'execution-plan.json']}))


def pilot():
    assert not (OUT/'pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.flow.old.optical.timeout);signal.alarm(55)
    rows=[]
    try:
        for n in [64,128]:rows.append(Response(128).run(n,f'pilot-{n}',4))
        point=max(r['operator_point_seconds']/r['operator_points'] for r in rows)
        estimate=point*51+rows[0]['stepping_seconds']/4*60+rows[1]['stepping_seconds']/4*252+15
        result=dict(classification='Counterexample candidate',rows=rows,forecast_seconds=estimate,upper_seconds=2*estimate,eligible=all(r['passed'] for r in rows) and 2*estimate<480,seconds=time.monotonic()-start)
        write(OUT/'pilot.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)
    except Exception as exc:
        write(OUT/'pilot-failure.json',dict(classification='Counterexample candidate',error=repr(exc),completed_paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


def warm_pilot():
    assert not (OUT/'warm-pilot.json').exists();start=time.monotonic();signal.signal(signal.SIGALRM,old.flow.old.optical.timeout);signal.alarm(35)
    rows=[Response(128).run(n,f'warm-pilot-{n}',8,f'pilot-{n}') for n in [64,128]]
    point=max(r['operator_point_seconds']/r['operator_points'] for r in rows)
    estimate=point*51+rows[0]['stepping_seconds']/4*56+rows[1]['stepping_seconds']/4*248+15
    result=dict(classification='Counterexample candidate',rows=rows,forecast_seconds=estimate,upper_seconds=2*estimate,eligible=all(r['passed'] for r in rows) and 2*estimate<480,seconds=time.monotonic()-start,
        optimization='Start each Krylov solve at the existing state or prior stage; retain full1e-10 residual. Continue actual saved prefixes; same equation and gates.',
        previous_dispatch_eligible=False)
    write(OUT/'warm-pilot.json',result);print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True);signal.alarm(0)


def budget():
    p=json.loads((OUT/'warm-pilot.json').read_text()); assert all(r['passed'] for r in p['rows']) and p['upper_seconds']<500
    assert not (OUT/'execution-plan.json').exists()
    files=[Path(__file__),Path(old.__file__),OUT/'first-pilot-producer.py',OUT/'warm-pilot-producer.py',old.OUT/'execution-plan.json',old.OUT/'bank-result.json']
    files+=list(old.OUT.glob('bank-*/*.npz'))+[OUT/f'warm-pilot-{n}.npz' for n in [64,128]]
    write(OUT/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=510,original_cap_seconds=480,original_dispatch_eligible=False,
        forecast_seconds=p['forecast_seconds'],upper_seconds=p['upper_seconds'],
        reassessment='Actual prefixes require6-8 full-equation Krylov iterations and predict249s, or498s with the retained2x margin. Reusing prior stages as guesses reduces the predicted bound from509s to498s. A single510s cap covers the same necessary three paths; no production was dispatched under the insufficient480s plan.',
        cheaper_alternatives='The160.66s split calculation completed but failed time accuracy; reuse its corrected backgrounds and coefficient banks, not its inaccurate trajectory. Continue the accepted monolithic8-step prefixes. An extra finer clock would cost more and would leave the partition problem unresolved.',
        limits='Numerical finite-dimensional first variation at prescribed background material motion. Additional mechanics, deep/exterior sources and GR closure remain open.',
        paths=[[64,128,'warm-pilot-64'],[128,128,'warm-pilot-128'],[128,64,None]],
        gates=dict(time=.02,background_time=.02,energy=1e-8,species=1e-8,linear=1e-10),
        stop='Stop on any failed path or510s. Atomic checkpoints retain work. No finer clocks or automatic further budget change.',bindings={str(f):sha(f) for f in files}))


def production():
    assert not (OUT/'result.json').exists();plan=json.loads((OUT/'execution-plan.json').read_text())
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    start=time.monotonic();signal.signal(signal.SIGALRM,old.flow.old.optical.timeout);signal.alarm(plan['hard_cap_seconds']);rows=[]
    try:
        for steps,reference,restart in plan['paths']:
            row=Response(reference).run(steps,f'steps-{steps}-reference-{reference}',restart=restart);rows.append(row)
            if not row['passed']:break
        comparisons={};indices=[0,1,2,3,5,6]
        if len(rows)==3 and all(r['passed'] for r in rows):
            def history(steps,reference):
                z=np.load(OUT/f'steps-{steps}-reference-{reference}.npz');t=np.linspace(0,old.flow.old.END,17);ids=[int(np.argmin(abs(z['t']-tt))) for tt in t]
                assert np.max(abs(z['t'][ids]-t))<1e-18
                return z['moments'][ids][:,indices]
            fine=history(128,128);norm=np.maximum(np.max(np.sum(abs(fine),axis=2),axis=0),1.)
            for key,value in [('time',history(64,128)),('background_time',history(128,64))]:comparisons[key]=np.max(np.sum(abs(value-fine),axis=2)/norm,axis=0).tolist()
        passed=len(rows)==3 and all(r['passed'] for r in rows) and max([max(v) for v in comparisons.values()],default=1)<.02
        result=dict(classification='Counterexample candidate',passed=passed,paths=rows,comparisons=comparisons,comparison_order=['photon_energy','material_energy','neutral_count','momentum_impulse','photon_radial_pressure','material_pressure'],seconds=time.monotonic()-start,
            actual_monolithic_radiation_thermal_H_response_evolved=len(rows)==3,additional_material_motion_evolved=False,full_GR_feedback=False,final_charge_solved=False)
        write(OUT/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'production-failure.json',dict(classification='Counterexample candidate',error=repr(exc),completed_paths=rows,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


if __name__=='__main__':globals()[sys.argv[1]]()
