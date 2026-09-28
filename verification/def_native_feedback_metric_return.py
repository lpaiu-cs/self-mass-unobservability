"""Apply the saved additional GR field as a compensated transport increment."""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,resource,signal,sys,time
import numpy as np
import def_native_updated_joint_feedback as previous
import verify_native_updated_joint_feedback as returned

prior=previous.prior;metric=prior.metric;C=prior.C;G=prior.G
OUT=previous.OUT/'metric-return';LAPSE=OUT/'lapse';RESPONSE=OUT/'response'
write=prior.write;sha=prior.sha;PATHS=previous.PATHS


def positive_increment(old,delta):
    """max(old+delta,0)-max(old,0), without losing a sub-ulp delta."""
    return np.where(old>=0,np.where(delta>=-old,delta,-old),np.where(delta>-old,old+delta,0.))


def prepare():
    assert not OUT.exists();LAPSE.mkdir(parents=True);RESPONSE.mkdir()
    assert json.loads((returned.OUT/'audit.json').read_text())['passed']
    for x in [-2.,0.,2.]:
        for dx in [-4.,-.5,0.,.5,4.]:
            assert positive_increment(np.array(x),np.array(dx))==max(x+dx,0)-max(x,0)
    assert positive_increment(np.array(2.),np.array(1e-30))==1e-30
    write(OUT/'plan.json',dict(classification='Counterexample candidate',before_checkpoint='60877e975',
        claim='Apply the actual returned material/photon GR field and SDIRK-stage outgoing photons to asymptotic lapse, then to actual simultaneous photons/energy/H. Return those incremental transfers to matter/GR and measure the next finite correction.',
        reuse='Preserve completed background,primary GR,Phase134 motion and all large response histories. Evolve only the new compensated GR forcing on the same531 cells,8 angles,152 frequencies,3.434ms and original three paths. No new native states or nonlinear physical trajectory.',
        compensation='Never add the sub-ulp metric increment to the primary metric. Streaming and collision drives are affine. For signed frequency upwinding use the exact positive-part increment about the existing primary frequency drift, including sign crossings; do not choose donors from the correction alone.',
        boundary='Propagate actual signed outgoing angular luminosity with its accepted SDIRK stage times and weights, preserving the exact cumulative numerical energy port. This is a discrete-stage representation, not a continuum emission certificate.',
        budgets=dict(lapse_seconds=120,response_pilot_seconds=65,response_production_seconds=650,material_pilot_seconds=45,material_production_seconds=350,source_seconds=75,GR_seconds=120,CPU_processes=3,threads_each=1,response_total_virtual_GiB=9,material_total_virtual_GiB=6),
        forecast='Lapse owner previously completed three paths within90s; four controls allowed120s. Response uses concurrent4/8/4 prefixes and measured17-point/remaining-step cost with2x margin under650s. Material uses previous full late-CFL raw-call counts and a measured concurrent prefix,2x margin under350s.',
        gates=dict(lapse_quadrature=.002,time=.02,background=.02,energy=1e-8,H=1e-8,linear=1e-12,angular_port=1e-12,ray=1e-10,velocity_jet=1e-4,primitive=1e-10,directional=.002,pressure=.002,GR_independent=1e-9),
        stop='Stop at a failed gate or measured budget. No automatic extra clock,mesh,horizon,EOS support or repeated waveform sweep. Preserve failures and accepted prefixes.',
        limits='Finite compensated response about the declared updated background. Phase134 material/photon mismatch is a separate remaining error. A small new GR correction is not a contraction theorem, nonlinear GR or final-charge enclosure.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(previous.__file__),Path(returned.__file__),Path(metric.__file__),returned.OUT/'audit.json',returned.OUT/'closure-audit.json',prior.LAPSE/'result.json']}))
    write(OUT/'signed-frequency-check.json',dict(classification='Proven',passed=True,
        scope='For real x,d, the branch formula computes max(x+d,0)-max(x,0). Represented scalar cases cover both signs,zero,crossings; a1e-30 increment beside2 is retained. This identity does not certify a continuum frequency discretization.'))


lapse_prior=SimpleNamespace(**dict(vars(metric.prior),OUT=returned.GR,new=SimpleNamespace(**dict(vars(prior.wave),OUT=returned.GR))))


class Lapse(metric.Lapse):
    run=FunctionType(metric.Lapse.run.__code__,dict(vars(metric),OUT=LAPSE,prior=lapse_prior))

    def boundary(self,d,label,order):
        n,ref=map(int,str(label).split('-reference-'));p=np.load(returned.photon_path(n,ref))
        times=p['accepted_angular_times'];lum=p['accepted_angular_luminosity'];h=d['t'][-1]/n
        gamma=1-1/np.sqrt(2);weights=np.tile([1-gamma,gamma],n)*h
        assert len(times)==len(lum)==2*n
        mu,mw,ids,solution,invariant=self.rays(order,d['t'][-1]);photon=[];energy=[]
        for now in d['t']:
            count=int(round(now/h));assert abs(count*h-now)<1e-18
            tt=times[:2*count];packet=weights[:2*count,None]*lum[:2*count,ids]*mu*mw
            if count==0:photon.append(0.);energy.append(0.);continue
            assert np.max(tt)<=now+1e-18
            z=solution.sol((now-tt)/d['t'][-1]).reshape(3,len(mu),len(tt)).transpose(0,2,1)
            rp=self.r0*z[0];cc=self.bg.metric(rp.ravel()/self.model.m.R)[3].reshape(rp.shape)
            kernel=z[2]/self.r0-z[1]**2/(rp*cc)
            photon.append(float(G/C**4*np.sum(packet*kernel,dtype=np.longdouble)))
            energy.append(float(np.sum(packet,dtype=np.longdouble)))
        energy=np.array(energy);err=float(np.max(abs(energy-d['outer_cumulative_energy_erg']))/max(np.max(abs(energy)),1e-300))
        assert err<1e-12
        gx,gw=np.polynomial.legendre.leggauss(order);z=(gx+1)/2
        _,N,b,_,_=self.bg.metric(self.r0/(z*self.model.m.R))
        return np.array(photon),energy,dict(ray_invariant=invariant,emitted_energy_relative=err,
            unit_ADM_lapse_kernel_per_cm=float(np.sum(gw/(2*self.r0*N*b**1.5))),actual_accepted_SDIRK_angular_emission=True)


def lapse():
    assert not (LAPSE/'result.json').exists();prior.configure();start=time.monotonic()
    signal.signal(signal.SIGALRM,prior.run.flow.old.optical.timeout);signal.alarm(120)
    m=Lapse();paths=[('128-reference-128',8),('128-reference-128',4),('64-reference-128',8),('128-reference-64',8)]
    rows=[m.run(label,q) for label,q in paths];fine=np.load(LAPSE/'metric-128-reference-128-g8.npz')
    keys=['delta_log_lapse','delta_lambda','delta_u','delta_log_speed'];errors={}
    for name,label,q in [('quadrature','128-reference-128',4),('time','64-reference-128',8),('background','128-reference-64',8)]:
        d=np.load(LAPSE/f'metric-{label}-g{q}.npz')
        errors[name]={k:float(np.max(abs(d[k]-fine[k]))/max(np.max(abs(fine[k])),1e-300)) for k in keys}
    old=np.load(prior.LAPSE/'corrected/metric-128-g8.npz')
    ratios={k:float(np.max(abs(fine[k]))/max(np.max(abs(old[k])),1e-300)) for k in keys}
    result=dict(classification='Counterexample candidate',passed=max(errors['quadrature'].values())<.002 and max(v for name,e in errors.items() if name!='quadrature' for v in e.values())<.02,
        errors=errors,increment_over_primary_metric=ratios,paths=rows,seconds=time.monotonic()-start,
        actual_returned_GR_lapse_computed=True,actual_transport_applied=False,final_charge_solved=False)
    result['passed']=result['passed'] and max(r['ray_invariant'] for r in rows)<1e-10
    write(LAPSE/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


class Response(previous.Response):
    def __init__(self,reference,steps=128):
        super().__init__(reference);self.primary_metric=self.g
        self.g=dict(np.load(LAPSE/f'metric-{steps}-reference-{reference}-g8.npz'))
        assert np.array_equal(self.g['t'],self.t) and np.array_equal(self.g['radius_E'],self.r)
        self.motion.fill(0.);self.mechanical.fill(0.);self.xi.fill(0.);self.energy_offset.fill(0.)
        self.branch_crossings=0;self.zero_primary_drifts=0;self.max_drift_ratio=0.

    def source(self,t):
        j=max(0,min(np.searchsorted(self.t,t,side='left')-1,len(self.t)-2));h=self.t[j+1]-self.t[j];a=(t-self.t[j])/h
        g=self.primary_metric;blend=lambda k:(1-a)*g[k][j]+a*g[k][j+1]
        self.primary_omega=-C*self.cc[:,None]*self.mu[None,:]*(blend('delta_nu_prime')+blend('delta_u_prime'))[:,None]-(g['delta_u'][j+1]-g['delta_u'][j])[:,None]/h-self.mu[None,:]**2*g['delta_lambda_interval_rate'][j,:,None]
        return super().source(t)

    def frequency(self,I,omega):
        base=self.primary_omega;up=positive_increment(base,omega);down=up-omega
        self.branch_crossings+=int(np.count_nonzero(((base>0)&(omega<-base))|((base<0)&(omega>-base))))
        self.zero_primary_drifts+=int(np.count_nonzero((base==0)&(omega!=0)))
        self.max_drift_ratio=max(self.max_drift_ratio,float(np.max(np.divide(abs(omega),abs(base),out=np.zeros_like(omega),where=base!=0))))
        count=I*self.num;u=count*up[:,:,None]*self.E/(self.extended[2:]-self.E)
        d=count*down[:,:,None]*self.E/(self.E-self.extended[:-2]);change=-u-d
        change[:,:,1:]+=u[:,:,:-1];change[:,:,:-1]+=d[:,:,1:]
        en=d[:,:,0]+u[:,:,-1];ee=d[:,:,0]*self.extended[0]+u[:,:,-1]*self.extended[-1]
        work=(count*self.E*omega[:,:,None]).sum(2)
        scale=np.max((abs(u)+abs(d)).sum(2));escale=np.max((abs(u)+abs(d))@self.E)
        error=max(float(np.max(abs(change.sum(2)+en))/max(scale,1e-300)),float(np.max(abs(change@self.E+ee-work))/max(escale,1e-300)))
        return change/self.num,en,ee,work,error


runner_ns=dict(previous.namespace,OUT=RESPONSE)
exec(compile(previous.runner,__file__,'exec'),runner_ns);Response.run=runner_ns['run']
dispatch=FunctionType(previous.parallel.dispatch.__code__,dict(vars(previous.parallel),OUT=RESPONSE,__file__=__file__))


def worker(n,ref,label,limit,restart):
    cap=3*1024**3;resource.setrlimit(resource.RLIMIT_AS,(cap,cap));prior.configure();start=time.monotonic();cpu=time.process_time()
    model=Response(ref,n);row=model.run(n,label,limit,restart)
    earlier=json.loads((RESPONSE/f'{restart}.json').read_text()) if restart else {}
    row.update(worker_wall_seconds=time.monotonic()-start,worker_CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        frequency_branch_crossings=model.branch_crossings+earlier.get('frequency_branch_crossings',0),
        zero_primary_drifts=model.zero_primary_drifts+earlier.get('zero_primary_drifts',0),
        maximum_increment_over_primary_drift=max(model.max_drift_ratio,earlier.get('maximum_increment_over_primary_drift',0.)))
    write(RESPONSE/f'{label}.json',row);assert row['passed']


def pilot():
    assert json.loads((LAPSE/'result.json').read_text())['passed'];assert not (RESPONSE/'pilot.json').exists()
    rows,elapsed=dispatch([(64,128,'pilot-64-128',4,None),(128,128,'pilot-128-128',8,None),(128,64,'pilot-128-64',4,None)],65)
    estimates=[r['operator_point_seconds']/r['operator_points']*17+r['stepping_seconds']/r['new_steps']*(r['steps']-r['completed_steps'])+20 for r in rows]
    p=dict(classification='Counterexample candidate',rows=rows,forecast_each_seconds=estimates,upper_seconds=2*max(estimates),eligible=all(r['passed'] for r in rows) and 2*max(estimates)<650,seconds=elapsed)
    write(RESPONSE/'pilot.json',p);print(json.dumps(p),flush=True)
    if p['eligible']:write(RESPONSE/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,hard_cap_seconds=650,
        paths=[[n,r,f'pilot-{n}-{r}'] for n,r in PATHS],bindings={str(f):sha(f) for f in [Path(__file__),OUT/'plan.json',LAPSE/'result.json',RESPONSE/'pilot.json']}))


production_source=prior.replace(inspect.getsource(previous.parallel.production),
    'norm=np.maximum(np.max(np.sum(abs(fine),axis=2),axis=0),1.)',
    'norm=np.maximum(np.max(np.sum(abs(fine),axis=2),axis=0),1e-300)')
production_scope=dict(vars(previous.parallel),OUT=RESPONSE,dispatch=dispatch)
exec(compile(production_source,__file__,'exec'),production_scope);production=production_scope['production']


if __name__=='__main__':
    if sys.argv[1]=='worker':worker(int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],None if sys.argv[5]=='None' else int(sys.argv[5]),None if sys.argv[6]=='None' else sys.argv[6])
    else:globals()[sys.argv[1]]()
