"""Counterexample candidate: carry actual Radau collisions into free matter.

One bounded replay records missing stage inputs and must reproduce177 exactly.
Then integrate the same native material flux with those inputs, before any
full-horizon calculation. Canonical-linear interpolation failure stays frozen.
"""
from pathlib import Path
from types import FunctionType
import gc,json,resource,shutil,sys,time
import numpy as np
import sympy as sp
import couple_native_neutral_transport as prior

OUT=Path('native-stage-collisions178-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha;LD=prior.LD;AMP=prior.AMP
CAPS=dict(prepare=10,capture=220,material=60)


def prepare():
    assert not OUT.exists() and not read(OLD/'material-prefix-result.json')['passed'];OUT.mkdir()
    p=read(OLD/'pilot-result.json');assert p['passed']
    for src in list(OLD.glob('face-*.npz'))+[OLD/'normalization.json',OLD/'photon-conservation-plan.json']:
        shutil.copyfile(src,OUT/src.name)
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    for src in (OLD/'sweep-0').rglob('*.npz'):
        dst=OUT/src.relative_to(OLD);shutil.copyfile(src,dst);assert sha(src)==sha(dst)
    files=[Path(__file__),Path(prior.__file__),OLD/'pilot-result.json',OLD/'material-prefix-result.json']
    files.extend(OLD/f'sweep-1/photons/pilot-{n}.npz' for n in [64,128])
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='4a63b5900',
        claim='Replace canonical17-knot collision interpolation by the actual accepted Radau stage history, then test original0.2percent neutral-transport consistency in the same short coupled response.',
        hypothesis='The177joint H equation passes but the free material forced by a canonical-linear collision history differs66percent in H-C. The stage history of its strong collision input is missing from that consumer. No unique-cause claim before actual return.',
        replay_reason='177saved neutral transport stages but not gas/collision/impulse stages. These missing same-solution inputs cannot be recovered from its two canonical prefix snapshots. Replay ONLY the acceptedT/16prefix, with unchanged equations and stage times; require bit-identical original physical arrays. No full path, new waveform, clock or resolution.',
        material_method='Use the actual two-stage collocation polynomial C(t) on each accepted Radau interval. Evolve U=z-C(t) with native F(z), z=U+C(t), using explicit SSPRK3 on the original64/128macro clocks and native hydro CFL, cut at recorded source-stage edges. Record floor projection and all conservation ledgers. This source removal is in the forced free-material equation only; the stiff photon/thermal/H Radau equation is unchanged. No H overwrite.',
        method_choice='SSPRK3 is fixed before the trial so a quadratic collocation input is sampled at its midpoint as well as endpoints. No subsequent order/clock/mesh ladder if it fails.',
        gates=dict(replay_identity=0.,time=.02,transport=.002,directional=.002,conservation=1e-8,owner=1e-8,branch=.01),
        forecast=dict(previous_capture_equivalent_seconds=read(OLD/'pilot_retry-receipt.json')['seconds'],
            capture_range_seconds=[150,220],previous_material_pair_seconds=read(OLD/'material-prefix-retry-receipt.json')['seconds'],
            material_range_seconds=[15,60],assumption='Capture adds small O(stage*cell) arrays; material adds one SSP stage and exact source-edge cuts. Costs are not guarantees.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=3,full_horizon_authorized=False,
        stop='Any physical gate, replay difference, future input or action cap rejects the trial. Preserve177failure and all gates. No automatic full path or new refinement.',
        bindings={str(p):sha(p) for p in files}))
    x=sp.symbols('x');w=sp.Matrix([sp.Rational(3,2)*x-sp.Rational(3,4)*x*x,sp.Rational(3,4)*x*x-x/2])
    assert w.subs(x,1)==sp.Matrix([sp.Rational(3,4),sp.Rational(1,4)])
    assert w.subs(x,sp.Rational(1,3))==sp.Matrix([sp.Rational(5,12),-sp.Rational(1,12)])
    C,U,F=sp.symbols('C U F',cls=sp.Function);t=sp.symbols('t');z=U(t)+C(t)
    assert sp.simplify(sp.diff(z,t).subs(sp.diff(U(t),t),F(z))-(F(z)+sp.diff(C(t),t)))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Collocation integrated basis reproduces both Radau stage rows. z=U+C converts z_dot=F(z)+C_dot to U_dot=F(U+C). No stiff photon transformation or full numerical error theorem.'))


def initialize():
    prior.OUT=OUT;prior.initialize_coupled()


def capture():
    initialize();model=prior.Response;original=model.run.__globals__['stages'];rows=[]
    def stages(m,t,h,x,g,lus):
        result,mech=original(m,t,h,x,g,lus)
        for j,(_,gg,p,q,e,_,_) in enumerate(result):
            transport=m.neutral_rates[-2+j]/LD(AMP)
            rates=np.zeros((4,m.n),LD)
            rates[1]=(-np.sum(p*m.Eweight*m.mu[None,:,None],axis=(1,2))-e[2])/m.a
            rates[2]=(q[:,0]-mech[:,0])*m.eu
            rates[3]=q[:,1]*m.nu-transport
            m.stage_collisions.append(rates)
            m.stage_gas.append(gg*np.stack([m.eu,m.nu],axis=-1))
        return result,mech
    model.run=FunctionType(model.run.__code__,dict(model.run.__globals__,stages=stages),argdefs=model.run.__defaults__)
    for n in [64,128]:
        m=model(n);m.stage_collisions=[];m.stage_gas=[];row=m.run(n,f'pilot-{n}',n//16)
        file=OUT/f'sweep-1/photons/pilot-{n}.npz';p=dict(np.load(file));old=dict(np.load(OLD/f'sweep-1/photons/pilot-{n}.npz'))
        keys=['t','moments','material_history','collision_transfer','radial_ports','accepted_angular_times',
            'accepted_angular_luminosity','accepted_angular_quadrature_weights','actual_step_edges','delta_material','delta_packet_scaled_occupation']
        for k in keys:assert np.array_equal(p[k],old[k]),('Physical replay changed',n,k)
        p.update(native_neutral_stage_times=np.array(m.neutral_times),native_neutral_stage_weights=np.array(m.neutral_weights),
            native_neutral_stage_rates=np.array(m.neutral_rates),collision_stage_rates_scaled=np.array(m.stage_collisions),
            gas_stage_scaled=np.array(m.stage_gas))
        weights=p['accepted_angular_quadrature_weights'].astype(LD)
        total=np.sum(weights[:,None,None]*p['collision_stage_rates_scaled'],axis=0,dtype=LD)
        expected=np.stack([np.zeros(m.n),p['moments'][-1,3]/m.a,p['collision_transfer'][-1,:,0],p['collision_transfer'][-1,:,1]])/LD(AMP)
        err=(np.sum(abs(total-expected),axis=1,dtype=LD)/np.maximum(np.sum(abs(expected),axis=1,dtype=LD),1.)).astype(float).tolist()
        assert max(err)<1e-12,(n,err)
        row.update(bit_identical_physical_replay=True,stage_collision_integral_relative=err)
        np.savez_compressed(file,**p);write(file.with_suffix('.json'),row);rows.append(row)
        del m,p,old;gc.collect()
    write(OUT/'capture-result.json',dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows),rows=rows,
        actual_original_physics_unchanged=True,full_horizon_completed=False,final_charge_conclusion='unadjudicated'))


def dense(p):
    edges=p['actual_step_edges'].astype(LD);h=np.diff(edges);rates=p['collision_stage_rates_scaled'].reshape(-1,2,4,len(p['radius_E']))
    increments=h[:,None,None]*(LD('.75')*rates[:,0]+LD('.25')*rates[:,1])
    cumulative=np.concatenate([np.zeros((1,4,rates.shape[-1]),LD),np.cumsum(increments,axis=0,dtype=LD)])
    def C(t):
        assert 0<=t<=edges[-1]+1e-18,('Uncomputed stage input',t)
        k=int(np.clip(np.searchsorted(edges,t,side='left')-1,0,len(h)-1));s=(LD(t)-edges[k])/h[k]
        return cumulative[k]+h[k]*((LD('1.5')*s-LD('.75')*s*s)*rates[k,0]+(LD('.75')*s*s-LD('.5')*s)*rates[k,1])
    return C,edges


def material():
    assert read(OUT/'capture-result.json')['passed'];initialize();owner=prior.previous.run.c;old_transfer=owner.transfers
    def transfers(m,path,n):
        if Path(path)!=prior.paths(1)[0]:return old_transfer(m,path,n)
        # Collisions enter the exact known lift C(t); the flux owner adds none.
        m.transfer=np.zeros((len(m.t),4,m.n),LD)
    owner.transfers=transfers;rows=[];histories=[]
    for n in [64,128]:
        m=owner.Material(128,n);p=dict(np.load(OUT/f'sweep-1/photons/pilot-{n}.npz'));C,edges=dense(p)
        end=edges[-1];u=np.zeros((4,m.n),LD);ledger=np.zeros(4,LD);discard=ledger.copy();norm=ledger.copy();t=LD(0);count=0
        times=[t];values=[u.copy()];maximum=0.;started=time.monotonic()
        def rhs(t,u,closing=False,probe=1.):
            assert 0<=t<=end+1e-18
            m.closing_stage=closing
            try:return m.rhs(t,u+C(t),probe)
            finally:m.closing_stage=False
        while t<end-1e-18:
            r1,l1,cfl=rhs(t,u);next_edge=edges[np.searchsorted(edges,t+LD('1e-18'),side='right')]
            next_macro=m.t[-1]*(np.floor(float(t/(m.t[-1]/n))+1e-10)+1)/n
            h=min(LD(cfl)*64/n,next_edge-t,LD(next_macro)-t,end-t);assert h>0
            for attempt in range(10):
                v=u+h*r1;r2,l2,cfl2=rhs(t+h,v,True)
                w=LD('.75')*u+LD('.25')*(v+h*r2);r3,l3,cfl3=rhs(t+h/2,w)
                if h<=min(cfl2,cfl3)*64/n*(1+1e-10):break
                h=min(h/2,LD(min(cfl2,cfl3))*64/n)
            else:raise AssertionError('Native SSP3 CFL')
            u=u/3+LD(2)/3*(w+h*r3)
            ledger+=h*(l1.astype(LD)/6+l2.astype(LD)/6+LD(2)/3*l3.astype(LD))
            norm+=h*(np.sum(abs(r1),axis=1,dtype=LD)/6+np.sum(abs(r2),axis=1,dtype=LD)/6+LD(2)/3*np.sum(abs(r3),axis=1,dtype=LD))
            t+=h;full=u+C(t);active=m.active(t);discard+=np.sum(full[:,~active],axis=1,dtype=LD);full[:,~active]=0;u=full-C(t)
            count+=1;assert np.isfinite(u).all() and count<1000
            balance=float(np.max(abs(np.sum(u,axis=1,dtype=LD)+discard-ledger)/np.maximum(norm,1.)));maximum=max(maximum,balance)
            assert balance<1e-8,('Lifted material conservation',balance)
            times.append(t);values.append(u.copy())
        times=np.array(times);values=np.array(values);z=u+C(end)
        rates=[rhs(end,u,True,v)[0] for v in [.5,1.,2.]]
        probe=[(np.sum(abs(v-rates[1]),axis=1)/np.maximum(np.sum(abs(rates[1]),axis=1),1.)).astype(float).tolist() for v in [rates[0],rates[2]]]
        joint=np.sum(p['native_neutral_stage_weights'][:,None].astype(LD)*p['native_neutral_stage_rates'],axis=0,dtype=LD)/LD(AMP)
        mismatch=float(np.sum(abs(u[3]-joint),dtype=LD)/max(np.sum(abs(u[3]),dtype=LD),1.))
        differences=[];norms=[];gas_errors=[]
        for st,rate,gas in zip(p['native_neutral_stage_times'],p['native_neutral_stage_rates'],p['gas_stage_scaled']):
            j=int(np.clip(np.searchsorted(times,st,side='left')-1,0,len(times)-2));w=(st-times[j])/(times[j+1]-times[j])
            state=(1-w)*values[j]+w*values[j+1]+C(st)
            gas_errors.append(np.sum(abs(state[3]-gas[:,1]),dtype=LD))
            k,v,_,_=m.fields(st);F=sum(a*prior.flux(m,i,state) for i,a in [(k,1-v),(k+1,v)] if a)
            actual=-np.diff(F)*AMP;differences.append(np.sum(abs(actual-rate),dtype=LD));norms.append(np.sum(abs(actual),dtype=LD))
        rate_error=float(max(differences)/max(max(norms),LD('1e-290')))
        gas_error=float(max(gas_errors)/max(np.max(np.sum(abs(p['gas_stage_scaled'][:,:,1]),axis=1,dtype=LD)),1.))
        paired=owner.relative((p['moments'][-1,[1,2]]/AMP)[None],z[[2,3]][None])
        row=dict(classification='Counterexample candidate',steps=n,substeps=count,seconds=time.monotonic()-started,
            balance_relative=maximum,mechanical_H_relative=mismatch,native_H_rate_relative=rate_error,
            actual_stage_H_relative=gas_error,paired_E_H=paired,endpoint_probe=probe,
            owner=max(v['owner_error'] for v in m.cache.values()),physical_branch_ratio=m.physical_branch_ratio)
        row['passed']=bool(max([mismatch,rate_error]+paired)<.002 and np.max(probe)<.002 and row['owner']<1e-8 and row['physical_branch_ratio']<.01)
        rows.append(row);histories.append(z[None]);write(OUT/f'material-{n}.json',row)
        np.savez_compressed(OUT/f'material-{n}.npz',t=times,lifted_history_scaled=values,delta_scaled=z,ledger_scaled=ledger,discard_scaled=discard,norm_scaled=norm)
        print(json.dumps(row),flush=True);del m,p;gc.collect()
    errors=owner.relative(*histories);result=dict(classification='Counterexample candidate',passed=all(r['passed'] for r in rows) and max(errors)<.02,
        rows=rows,time_comparison=errors,full_horizon_completed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'material-result.json',result);assert result['passed'],result


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));prior.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
