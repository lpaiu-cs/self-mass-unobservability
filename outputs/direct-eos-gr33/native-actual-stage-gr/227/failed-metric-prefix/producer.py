"""Evaluate the same dense GR history at actual coupled stages, without secants."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import gc,inspect,json,os,resource,sys,time
import numpy as np
from scipy.interpolate import PPoly
import continue_dense_gr_return as previous

prior=previous.prior;bridge=prior.bridge;geometry=prior.geometry
OUT=Path('native-stage-metric227-work');INPUT=prior.INPUT
read,write,sha,LD,C,AMP=bridge.read,bridge.write,bridge.sha,bridge.LD,bridge.C,bridge.AMP
CAPS=dict(prepare=180,field1288=1800,field648=1800,field1284=1800,metric=1200,coarse=2400,fine=4800,audit=180)


def representation(n):
    raw=dict(np.load(INPUT/f'gr/source-{n}.npz'));knots,co=bridge.coefficients(raw)
    return raw,{k:PPoly(np.asarray(v[::-1],float),knots) for k,v in co.items()}


def prepare():
    failed=read(previous.OUT/'result.json');assert not failed['passed']
    assert all(read(previous.OUT/f'run-{n}.json')['passed'] for n in [64,128])
    assert not (OUT/'plan.json').exists();OUT.mkdir(exist_ok=True);files=[]
    for part in ['sweep-0','sweep-1/photons','sweep-1/material','gr','metric']:(OUT/part).mkdir(parents=True,exist_ok=True)
    for src in list((previous.OUT/'sweep-0').rglob('*.npz'))+[previous.OUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(previous.OUT);dst.parent.mkdir(parents=True,exist_ok=True)
        if dst.exists():assert sha(dst)==sha(src)
        else:os.link(src,dst)
        files += [src,dst]
    candidates=[0.]
    for n in [64,128]:
        path=previous.OUT/f'sweep-1/photons/return-{n}.npz';p=np.load(path);candidates+=list(p['joint_stage_times'])+list(p['actual_step_edges']);files.append(path)
    clock=[]
    for t in sorted(candidates):
        if not clock or t-clock[-1]>1e-18:clock.append(t)
    clock=np.array(clock)
    reconstruction=[]
    for n in [64,128]:
        raw,poly=representation(n);d=dict(raw,t=clock,original_clock=np.array(n))
        d.update({k:f(clock) for k,f in poly.items()})
        indices=np.array([np.argmin(abs(clock-t)) for t in raw['t']]);assert np.max(abs(clock[indices]-raw['t']))<1e-18
        jumps={k:float(np.max(abs(d[k][indices]-raw[k]))/max(np.max(abs(raw[k])),LD('1e-290'))) for k in poly}
        # The dense final polynomial closes BEFORE the saved floor operation.
        # Preserve the original post-floor endpoint, not a smoothed polynomial.
        for k in poly:d[k][indices]=raw[k]
        assert all(np.array_equal(d[k][indices],np.asarray(raw[k],float)) for k in poly)
        np.savez_compressed(OUT/f'gr/source-{n}.npz',**d)
        reconstruction.append(dict(clock=n,original_endpoint_values_exact=True,pre_floor_to_saved_endpoint_relative=jumps))
        files += [INPUT/f'gr/source-{n}.npz',prior.saved(n)]
    files += [previous.OUT/n for n in ['result.json','audit-receipt.json','coarse-receipt.json','fine-receipt.json']]+[OUT/'rate-alias.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Remove measured aliasing of the returned metric derivative at actual coupled stages, then evolve and adjudicate the same original two physical clocks.',
        evidence='226both physical paths pass but paired time max5.03484percent fails2percent. Actual Radau sampling of225piecewise-linear metric misintegrates u by5.05290percent and lambda by3.54655percent on64, while128is exact to2.36e-16. This diagnoses a concrete input quadrature error, not a proof that it explains all coupled error.',
        method='Reuse224dense source polynomials and characteristic GR solver. Evaluate both source histories on the union of EXISTING64/128stage and endpoint times. At original endpoints preserve their actual post-floor values; the final dense polynomial is pre-floor and must not overwrite that operation. Use directly propagated U_t and differentiate the exact-center mass constraint inside its same source segments. Floor jumps remain explicit events, not a uniform derivative theorem. Stage driver performs exact-time lookup; no metric secant or interpolated derivative. Same EOS, physical grid, amplitude, clocks, fronts and native stage equations.',
        packets='Use the same original two Radau angular luminosities as their linear collocation polynomial on each emission interval. Integrate this same dense emission along the null rays at requested intermediate times. Its whole-step integral is exactly the original3/4,1/4packet measure. Check energy and4/8ray quadrature; this replaces the old endpoint-only discrete packet approximation explicitly.',
        decision='Original2percent field/metric including pointwise rates,0.2percent quadrature,1e-12energy ports and all unchanged coupled stage/physical/native/constitutive/balance/branch/time gates. Failed226and188verdicts stay. No final charge or nonlinear/self-GR closure claim.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=12,
        forecast='224fine8GR took120s for30times; union count below. Three fields roughly8..20minutes if output count scaling holds; allow30minutes EACH. Metric20minutes.226actual15/29steps334/640s; allow40/80minutes. No new physical resolution or period; old physics must be rerun because its applied GR representation changes.',
        output_times=len(clock),actual_physical_steps=[15,29],capture='Copy actual accepted stage photon moments and boundary returns during evolution, avoiding a later missing-photon recovery solve.',
        stop='Any original acceptance gate, native anchor, provenance or generous wall cap. No automatic further refinement.218and223/224sources and plans remain unchanged.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as sp
    x,a,b=sp.symbols('x a b');lum=sp.Rational(3,2)*(1-x)*a+(sp.Rational(3,2)*x-sp.Rational(1,2))*b
    assert sp.integrate(lum,(x,0,1))==sp.Rational(3,4)*a+sp.Rational(1,4)*b
    f,e,k,q=sp.symbols('f e k q');assert sp.diff(k*f+q*e,f)==k and sp.diff(k*f+q*e,e)==q
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Same Radau emission integral and static-background linear mass-constraint derivative identities only. No uniform physical error theorem.'))
    write(OUT/'source-endpoint-check.json',dict(classification='Counterexample candidate',passed=True,rows=reconstruction))


class Response(bridge.Response):
    def setup(self,d,order):
        raw=dict(np.load(INPUT/f"gr/source-{int(d['original_clock'])}.npz"));knots,co=bridge.coefficients(raw)
        shape=co['baryon_g'].shape[:2];samples=dict(raw,t=np.zeros(shape[0]*shape[1]))
        samples.update({k:v.reshape((-1,)+raw[k].shape[1:]) for k,v in co.items()})
        bridge.base.gr.Response.setup(self,samples,order)
        source=self.source.reshape(*shape,-1).copy();direct=self.direct.reshape(*shape,-1).copy()
        bridge.base.gr.Response.setup(self,d,order)
        self.source=PPoly(source[::-1],knots);self.direct=PPoly(direct[::-1],knots)


def field(n,order):
    path=OUT/f'gr/fields-{n}-g{order}.npz'
    if path.exists():
        for p in [path,path.with_suffix('.json')]:assert sha(p)==read(OUT/'field-reuse.json')[str(p)]
        row=read(path.with_suffix('.json'))
    else:
        FunctionType(bridge.endpoint.initialize.__code__,dict(bridge.endpoint.initialize.__globals__,OUT=OUT))()
        m=Response();fn=FunctionType(bridge.base.gr.base.Response.run.__code__,dict(bridge.base.gr.base.Response.run.__globals__,OUT=OUT/'gr'))
        row=fn(m,n,order)
    old=np.load(INPUT/f'gr/fields-{n}-g{order}.npz');new=np.load(path)
    controls={k:bridge.endpoint.aligned(old,new,k) for k in ['U','U_t','U_x']}
    result=dict(classification='Counterexample candidate',passed=max(controls.values())<1e-9,original_output_reproduction=controls)
    write(OUT/f'field-{n}-g{order}-check.json',result);assert result['passed'],result
    return row


class Lapse(geometry.metric.Lapse):
    def boundary(self,d,n,order):
        p=np.load(prior.saved(n));edges=p['actual_step_edges'];lum=p['accepted_angular_luminosity'].reshape(-1,2,4)
        end=d['t'][-1];mu,mw,ids,solution,invariant=self.rays(order,end)
        gx,gw=np.polynomial.legendre.leggauss(order);photon=[];energy=[];norm=[]
        for now in d['t']:
            left=edges[:-1];right=np.minimum(edges[1:],now);mask=right>left
            lo=left[mask];hi=right[mask];tt=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*gx/2;ww=(hi-lo)[:,None]*gw/2
            theta=(tt-lo[:,None])/np.diff(edges)[mask,None]
            ll=1.5*(1-theta[:,:,None])*lum[mask,0,None]+(1.5*theta[:,:,None]-.5)*lum[mask,1,None]
            weight=ww.ravel()[:,None]*ll.reshape(-1,4)[:,ids]*mu*mw
            if not len(weight):photon.append(0.);energy.append(0.);norm.append(0.);continue
            ages=(now-tt.ravel())/end;assert np.min(ages)>=0
            state=solution.sol(ages).reshape(3,len(mu),len(ages)).transpose(0,2,1)
            rp=self.r0*state[0];mp=state[1];cc=self.bg.metric(rp.ravel()/self.model.m.R)[3].reshape(rp.shape)
            kernel=state[2]/self.r0-mp*mp/(rp*cc)
            photon.append(float(geometry.metric.G/C**4*np.sum(weight*kernel,dtype=LD)))
            energy.append(float(np.sum(weight,dtype=LD)));norm.append(float(np.sum(abs(weight),dtype=LD)))
        energy=np.array(energy);error=float(np.max(abs(energy-d['outer_cumulative_energy_erg']))/max(max(norm),1e-290));assert error<1e-12,error
        z=(gx+1)/2;_,N,b,_,_=self.bg.metric(self.r0/(z*self.model.m.R));kernel=float(np.sum(gw/(2*self.r0*N*b**1.5)))
        return np.array(photon),energy,dict(ray_invariant=invariant,emitted_energy_relative=error,unit_ADM_lapse_kernel_per_cm=kernel,same_Radau_dense_emission=True)


def metric():
    FunctionType(bridge.endpoint.initialize.__code__,dict(bridge.endpoint.initialize.__globals__,OUT=OUT))()
    m=Lapse();scope=dict(geometry.metric.Lapse.run.__globals__,OUT=OUT/'metric',prior=SimpleNamespace(OUT=OUT/'gr',new=SimpleNamespace(OUT=OUT/'gr'),centers=geometry.constraints.centers))
    fn=FunctionType(geometry.metric.Lapse.run.__code__,scope);rows=[]
    for n,q in [(128,8),(64,8),(128,4)]:
        rows.append(fn(m,n,q));path=OUT/f'metric/metric-{n}-g{q}.npz';d=dict(np.load(OUT/f'gr/source-{n}.npz'));z=dict(np.load(path));f=np.load(OUT/f'gr/fields-{n}-g8.npz')
        _,poly=representation(n);dot=dict(d);dot.update({k:p.derivative()(d['t']) for k,p in poly.items()})
        derivative=dict(delta_phi=f['U_t']/f['radius_E'],delta_Phi=np.zeros_like(f['U_t']))
        exact=geometry.constraints.centers(m.response,dot,derivative,q)
        z['actual_delta_lambda_rate']=exact['delta_lambda'][:,:-1];np.savez_compressed(path,**z)
    fine=dict(np.load(OUT/'metric/metric-128-g8.npz'));controls={}
    keys=list(geometry.KEYS)+['delta_u_t','actual_delta_lambda_rate'];keys.remove('delta_lambda_rate')
    for label,n,q in [('time',64,8),('quadrature',128,4)]:
        other=dict(np.load(OUT/f'metric/metric-{n}-g{q}.npz'))
        controls[label]={k:bridge.endpoint.aligned(other,fine,k) for k in keys}
    a,b=[np.load(OUT/f'gr/fields-128-g{q}.npz') for q in [4,8]]
    controls['field_quadrature']={k:bridge.endpoint.aligned(a,b,k) for k in ['U','U_t','U_x']}
    result=dict(classification='Counterexample candidate',passed=max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002 and max(controls['field_quadrature'].values())<.002 and max(r['ray_invariant'] for r in rows)<1e-10,controls=controls,rows=rows,actual_stage_fields=True,final_charge_conclusion='unadjudicated')
    write(OUT/'metric-result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


class StageDriver(geometry.ReturnOnly):
    def returned(self,t):
        j=int(np.argmin(abs(self.clock-t)));assert abs(self.clock[j]-t)<1e-18,('Uncomputed actual GR time',t)
        row={k:v[j].astype(LD) for k,v in self.g.items() if k.startswith('delta_') and v.shape==(len(self.clock),self.n)}
        row['delta_lambda_rate']=self.g['actual_delta_lambda_rate'][j].astype(LD)
        row['delta_log_volume']=3*row['delta_u']+row['delta_lambda'];row['delta_log_areal_radius']=row['delta_u']
        return row


def evolve(n):
    assert read(OUT/'metric-result.json')['passed'];geometry.ReturnOnly=StageDriver
    Model=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))()
    run=Model.run;stage=run.__globals__['stages'];boundary=Model.boundary_ports;moments=[];ports=[];native=[];branches=[];intervals=[]
    def observed(m,t,h,x,g,lus):
        pair,mechanical=stage(m,t,h,x,g,lus)
        for v in pair:
            xx=v[0];moments.append(np.array([np.sum(xx*m.Eweight,axis=(1,2)),np.sum(xx*m.Eweight*m.model.bulk.mu2[None,:,None],axis=(1,2)),np.sum(xx*m.Nweight,axis=(1,2))])*AMP)
        np.savez_compressed(OUT/f'last-pair-{n}.npz',time=t,step=h,photons=[v[0] for v in pair],gas=[v[1] for v in pair],x_initial=x,g_initial=g)
        write(OUT/f'capture-{n}.json',dict(actual_steps=len(moments)//2,time=float(m.stage_t[-1])))
        return pair,mechanical
    def port(m,t,x):
        value=boundary(m,t,x);ports.append(value.copy()*AMP);return value
    Model.boundary_ports=port;Model.run=FunctionType(run.__code__,dict(run.__globals__,stages=observed),argdefs=run.__defaults__)
    restart=None
    for j in [1,2]:
        m=Model(n);label=f'interval-{j}-{n}';row=m.run(n,label,j*n//16,restart);assert row['passed'];intervals.append(dict(row))
        native+=m.anchor_checks;branches+=m.branch_checks;path=OUT/f'sweep-1/photons/{label}.npz';z=np.load(path)
        np.savez_compressed(OUT/f'recovered-{n}.npz',times=z['joint_stage_times'],weights=z['joint_stage_weights'],photon_moments=moments,radial_ports=ports,collision_rates=z['joint_collision_rates_scaled'],angular=z['accepted_angular_luminosity'],endpoint_occupation=z['restart_x']*m.scale*AMP)
        checks=dict(newton=m.newton_iterations,stages=m.stage_log);restart=label;del m;gc.collect()
    dst=OUT/f'sweep-1/photons/return-{n}.npz';os.link(path,dst);os.link(path.with_suffix('.json'),dst.with_suffix('.json'))
    audit,_,_,_=geometry.prior.run.verify(dst,dst)
    weights=z['joint_stage_weights'];expected=np.sum(weights[:,None,None]*ports,axis=0,dtype=LD);old=z['radial_ports'][-1]
    error=float(np.max(abs(expected-old)/np.maximum(abs(old),LD('1e-290'))));assert error<1e-12,error
    anchor=np.load(prior.saved(n));assert np.array_equal(z['actual_step_edges'],anchor['actual_step_edges']) and np.array_equal(z['joint_stage_times'],anchor['joint_stage_times'])
    row.update(audit=audit,intervals=intervals,anchor_checks=native,branch_checks=branches,actual_stage_photon_moments_captured=True,captured_radial_port_relative=error,
        maximum_true_stage=max(v[-1]['relative'] for v in checks['newton']),maximum_true_physical_stage=max(max(v[-1]['moments']) for v in checks['newton']),same_saved_stage_equation=True,actual_return_time_evolved=True,self_GR_return_closed=False,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/f'run-{n}.json',row);write(OUT/f'checks-{n}.json',dict(classification='Counterexample candidate',**checks))


def audit():FunctionType(previous.audit.__code__,dict(previous.audit.__globals__,OUT=OUT))()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(12*1024**3,12*1024**3));prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action.startswith('field'):field(64 if action=='field648' else 128,4 if action=='field1284' else 8)
        elif action in ['coarse','fine']:evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
