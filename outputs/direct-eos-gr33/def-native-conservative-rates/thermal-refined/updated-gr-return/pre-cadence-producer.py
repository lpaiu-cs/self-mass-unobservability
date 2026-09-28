"""Updated physical histories -> GR/lapse -> simultaneous photons/heat/H."""
from pathlib import Path
from types import FunctionType
import inspect,textwrap,json,signal,sys,time
import numpy as np
import sympy as sp
import def_native_updated_gr_return as run
import def_native_characteristic_gr as wave
import verify_native_anisotropic_gr as constraints
import def_native_dynamic_lapse as metric
import def_native_monolithic_response as mono

OUT=run.OUT;GR=run.GR;LAPSE=run.LAPSE;C=run.C;G=run.flow.G
write=run.write;sha=run.sha;old=mono.old;physical=run.physical


class Model(physical.Coupled):
    def __init__(self):
        super().__init__();run.cache_inputs(self)


def configure():
    # Process-local routing only. Historical producers and outputs stay intact.
    run.flow.OUT=run.BG;run.flow.Coupled=Model
    wave.base.OUT=GR;wave.OUT=GR;constraints.OUT=GR
    metric.OUT=LAPSE;old.OUT=run.COLL


class Fields(wave.Response):
    run=FunctionType(wave.Response.run.__code__,dict(wave.Response.run.__globals__,OUT=GR))


def fields():
    assert json.loads((OUT/'capture-result.json').read_text())['passed'];assert not (GR/'result.json').exists()
    write(GR/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply corrected-EOS actual conserved sources to the same initial first-variation GR operator; retain full spatial fields for actual transport return.',
        gates=dict(time=.02,quadrature=.002,cadence=.002,full_history_outer=.002),budget_seconds=180,
        controls='Original64/128 backgrounds;4/8 radial rules;17/33 source knots; compare the outer free component against the previously accepted all129-stage readout. No extra physical time integration.',
        limits='Initial GR coefficients retain the declared frozen representation; compact fields plus bounded-potential iteration, not nonlinear GR.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(wave.__file__),Path(wave.base.__file__),OUT/'capture-result.json',GR/'source-128.npz',GR/'source-64.npz']}))
    start=time.monotonic();signal.signal(signal.SIGALRM,run.flow.old.optical.timeout);signal.alarm(180);configure()
    d=dict(np.load(physical.EV/'source-128.npz'));selected=np.arange(0,129,4)
    for key in d:
        if key=='t' or key in physical.capture.KEYS+['metric_stress_erg','inner_cumulative_energy_erg','outer_cumulative_energy_erg']:d[key]=d[key][selected]
    np.savez_compressed(GR/'source-128-cadence33.npz',**d)
    model=Fields();rows=[model.run(n,q) for n,q in [(128,8),(64,8),(128,4),('128-cadence33',8)]]
    fine=np.load(GR/'fields-128-g8.npz');norm=float(np.max(abs(fine['U'])))
    errors={k:float(np.max(abs(fine['U']-np.load(GR/f'fields-{label}.npz')['U'][::stride]))/norm)
        for k,label,stride in [('time','64-g8',1),('quadrature','128-g4',1),('cadence','128-cadence33-g8',2)]}
    previous=np.load(physical.GR/'wave-128-g8.npz');mass=float(d['M_cm'])
    free=-fine['direct_and_mass_stress_U'][:,-1]/mass
    errors['full_history_outer']=float(np.max(abs(free-previous['free_scalar'][::8]))/np.max(abs(previous['free_scalar'])))
    result=dict(classification='Counterexample candidate',passed=all(v<(.02 if k=='time' else .002) for k,v in errors.items()),
        errors=errors,paths=rows,seconds=time.monotonic()-start,actual_updated_source_fields=True,full_GR_feedback=False,final_charge_solved=False)
    write(GR/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


class Lapse(metric.Lapse):
    run=FunctionType(metric.Lapse.run.__code__,dict(vars(metric),OUT=LAPSE/'corrected'))

    def boundary(self,d,steps,order):
        t=d['t'];end=t[-1];mu,mw,ids,solution,invariant=self.rays(order,end)
        ports=np.load(run.BG/f'accepted-ports-{steps}.npz');lum=ports['angular_luminosity'];h=float(ports['h'])
        gx,gw=np.polynomial.legendre.leggauss(order);photon=[];energy=[]
        for now in t:
            count=int(round(now/h));assert abs(count*h-now)<1e-18
            if count==0:photon.append(0.);energy.append(0.);continue
            tt=(np.arange(count)[:,None]+(gx+1)/2)*h;tt=tt.ravel()
            packet=(h*gw/2)[None,:,None]*lum[:count,None,ids]*mu*mw;packet=packet.reshape(len(tt),len(mu))
            state=solution.sol((now-tt)/end).reshape(3,len(mu),len(tt)).transpose(0,2,1)
            rp=self.r0*state[0];mp=state[1];cc=self.bg.metric(rp.ravel()/self.model.m.R)[3].reshape(rp.shape)
            kernel=state[2]/self.r0-mp*mp/(rp*cc)
            photon.append(float(G/C**4*np.sum(packet*kernel,dtype=np.longdouble)));energy.append(float(np.sum(packet,dtype=np.longdouble)))
        energy=np.array(energy);photon=np.array(photon)
        error=float(np.max(abs(energy-d['outer_cumulative_energy_erg']))/max(energy[-1],1.));assert error<1e-12
        z=(gx+1)/2;_,N,b,_,_=self.bg.metric(self.r0/(z*self.model.m.R))
        kernel=float(np.sum(gw/(2*self.r0*N*b**1.5)))
        return photon,energy,dict(ray_invariant=invariant,emitted_energy_relative=error,unit_ADM_lapse_kernel_per_cm=kernel,
            accepted_stage_angular_emission=True)


def lapse():
    assert json.loads((GR/'result.json').read_text())['passed'];assert not (LAPSE/'result.json').exists();configure()
    write(LAPSE/'plan.json',dict(classification='Counterexample candidate',budget_seconds=90,
        claim='Use actual accepted-stage angular photon emission to normalize updated lapse at infinity, then export compensated radial/angular/frequency GR forcing.',
        correction='Piecewise-constant stage luminosity integrates to the saved exact accepted port. Do not restore the previously rejected endpoint-trapezoid energy mismatch.',
        gates=dict(accepted_port=1e-12,ray=1e-10,quadrature=.002,time=.02),
        limits='Same finite angular-bin and initial linear GR model; discarded-material transport and exterior-generated/scattered scalar are not newly closed.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(metric.__file__),GR/'result.json',run.BG/'accepted-ports-128.npz',run.BG/'accepted-ports-64.npz']}))
    start=time.monotonic();signal.signal(signal.SIGALRM,run.flow.old.optical.timeout);signal.alarm(90)
    mu,a,b=sp.symbols('mu a b');assert sp.simplify(sp.integrate(mu,(mu,a,b))-(b-a)*(a+b)/2)==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='For angular-bin constant occupation and accepted-stage constant luminosity, the exact angular flux and time integral equal the shared numerical photon port. This is not a continuum angular or time certificate.'))
    model=Lapse();rows=[model.run(n,q) for n,q in [(128,8),(128,4),(64,8)]]
    fine=np.load(LAPSE/'corrected/metric-128-g8.npz');keys=['delta_log_lapse','delta_lambda','delta_u','delta_log_speed']
    errors={name:{k:float(np.max(abs(fine[k]-np.load(LAPSE/f'corrected/metric-{label}.npz')[k]))/max(np.max(abs(fine[k])),1e-300)) for k in keys}
        for name,label in [('quadrature','128-g4'),('time','64-g8')]}
    result=dict(classification='Counterexample candidate',passed=max(errors['quadrature'].values())<.002 and max(errors['time'].values())<.02 and max(r['ray_invariant'] for r in rows)<1e-10 and max(r['emitted_energy_relative'] for r in rows)<1e-12,
        errors=errors,paths=rows,seconds=time.monotonic()-start,actual_accepted_photon_lapse=True,final_charge_solved=False)
    write(LAPSE/'result.json',result);signal.alarm(0);print(json.dumps(result),flush=True)


def bank():
    assert json.loads((LAPSE/'result.json').read_text())['passed'];configure()
    write(run.COLL/'plan.json',dict(classification='Counterexample candidate',budget_seconds=90,new_native_states=0,
        claim='Rebuild actual moving thermal/H/photon coefficient jets on corrected-EOS states using the existing owner and unchanged representative derivative gate1e-4.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(old.__file__),LAPSE/'result.json',physical.OUT/'bank.npz',run.BG/'coupled-128.npz',run.BG/'coupled-64.npz']}))
    old.bank()


def replace(source,before,after):
    assert source.count(before)==1,before
    return source.replace(before,after)


# Same simultaneous equation; evaluate its nonautonomous coefficients at the
# actual SDIRK stages and retain complete states for the next material return.
runner=textwrap.dedent(inspect.getsource(mono.Response.run))
runner=replace(runner,'t=k*h;c=self.local(t+h/2);','t=k*h;c=self.local(t+gamma*h);')
runner=replace(runner,'s2,l2,e2=self.source(t+h);',
    'c=self.local(t+h);inverse=self.inverse(c,gamma*h);qg=self.gas(c["q"],c["qb"],c["qe"])\n        s2,l2,e2=self.source(t+h);')
runner=runner.replace('rtol=1e-10','rtol=1e-12').replace('err<1e-10','err<1e-12').replace('self.max_residual<1e-10','self.max_residual<1e-12')
runner=replace(runner,'impulse=np.zeros(self.n)',
    'impulse=np.zeros(self.n);ports=np.zeros((2,2));port_history=[];photon_history=[];gas_history=[];transfer_history=[]')
runner=replace(runner,'times.append(t);records.append(',
    'port_history.append(ports.copy()*AMPLITUDE);photon_history.append(x*self.scale*AMPLITUDE);gas_history.append(g.copy()*AMPLITUDE);transfer_history.append(g*np.stack([self.eu,self.nu],axis=-1)*AMPLITUDE)\n        times.append(t);records.append(')
runner=replace(runner,"times=z['t'].tolist();records=list(z['moments']);",
    "ports=z['radial_ports'][-1]/AMPLITUDE;port_history=list(z['radial_ports']);photon_history=list(z['photon_history_scaled_occupation']);gas_history=list(z['material_history']);transfer_history=list(z['collision_transfer']);self.max_residual=previous['linear_relative'];self.max_iterations=previous['max_Krylov_iterations'];times=z['t'].tolist();records=list(z['moments']);")
runner=replace(runner,'p2,g2,es2,_=self.collision(c,z,gz,True)',
    'p2,g2,es2,_=self.collision(c,z,gz,True)\n        ports+=h*((1-gamma)*self.boundary_ports(t+gamma*h,y)+gamma*self.boundary_ports(t+h,z))')
runner=replace(runner,'completed_steps=k+1,t=times,moments=records)',
    'completed_steps=k+1,t=times,moments=records,ports=ports,photon_history=photon_history,gas_history=gas_history,port_history=port_history,transfer_history=transfer_history)')
runner=replace(runner,'radius_E=self.r)',
    'radius_E=self.r,radial_ports=port_history,photon_history_scaled_occupation=photon_history,material_history=gas_history,collision_transfer=transfer_history)')
namespace=dict(vars(mono),OUT=run.RESPONSE);exec(compile(runner,__file__,'exec'),namespace)


class Response(mono.Response):
    run=namespace['run']

    def boundary_ports(self,t,x):
        k=max(0,min(np.searchsorted(self.t,t,side='left')-1,15));f=(t-self.t[k])/(self.t[k+1]-self.t[k])
        I=(1-f)*self.I[k]+f*self.I[k+1]
        speed=((1-f)*self.g['delta_log_speed'][k]+f*self.g['delta_log_speed'][k+1])/old.AMPLITUDE
        actual=x*self.scale+I*speed[:,None,None];ports=[]
        for j,area,mask in [(0,self.area[0],self.mu<0),(-1,self.area[-1],self.mu>0)]:
            values=4*np.pi*C*area*np.sum(actual[j,mask]*(self.w*self.mu)[mask,None]*self.num,axis=0)
            ports.append([float(values.sum()),float(values@self.E)])
        return np.array(ports)


def pilot():
    assert json.loads((run.COLL/'bank-result.json').read_text())['passed'];configure();folder=run.RESPONSE
    assert not (folder/'pilot.json').exists()
    write(folder/'plan.json',dict(classification='Counterexample candidate',
        claim='Apply actual updated-source GR/lapse fields to simultaneous moving photon/thermal/H equations across all531 physical cells.',
        reuse='Same64/128 clocks,17 background knots,8 angles,152 frequencies and3.434ms. Actual SDIRK stage-time coefficients; GMRES/full residual1e-12. Store full angular response,thermal/H and radial ports for the next actual material return.',
        limits='Prescribed background baryon/momentum/inventory; this sweep does not solve their extra mechanical response or a GR fixed point.',
        budgets=dict(pilot_seconds=65,production_seconds=650,CPU_threads=1,memory_GB=3),
        gates=dict(time=.02,background_time=.02,energy=1e-8,species=1e-8,linear=1e-12),
        forecast='Equal-horizon4/8 steps;51 operator points plus remaining60/248 steps and25s setup/export,2x margin; no dispatch above650s. Later Krylov cost extrapolated.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(mono.__file__),Path(old.__file__),run.COLL/'bank-result.json',LAPSE/'result.json']}))
    (folder/'expanded-run.py').write_text(runner)
    start=time.monotonic();signal.signal(signal.SIGALRM,run.flow.old.optical.timeout);signal.alarm(65);rows=[]
    try:
        rows=[Response(128).run(n,f'pilot-{n}',k) for n,k in [(64,4),(128,8)]]
        a=np.load(folder/'pilot-64.npz')['moments'][-1:,[0,1,2,3,5,6]]
        b=np.load(folder/'pilot-128.npz')['moments'][-1:,[0,1,2,3,5,6]]
        comparison=(np.sum(abs(a-b),axis=2)/np.maximum(np.sum(abs(b),axis=2),1.)).ravel().tolist()
        point=max(r['operator_point_seconds']/r['operator_points'] for r in rows)
        estimate=point*51+rows[0]['stepping_seconds']/4*60+rows[1]['stepping_seconds']/8*248+25
        result=dict(classification='Counterexample candidate',rows=rows,equal_horizon_comparison=comparison,
            forecast_seconds=estimate,upper_seconds=2*estimate,eligible=all(r['passed'] for r in rows) and max(comparison)<.02 and 2*estimate<650,seconds=time.monotonic()-start)
        write(folder/'pilot.json',result);print(json.dumps(result),flush=True)
        if result['eligible']:write(folder/'execution-plan.json',dict(classification='Counterexample candidate',eligible=True,
            hard_cap_seconds=650,paths=[[64,128,'pilot-64'],[128,128,'pilot-128'],[128,64,None]],
            bindings={str(p):sha(p) for p in [Path(__file__),Path(run.__file__),folder/'plan.json',folder/'pilot.json',run.COLL/'bank-result.json',LAPSE/'result.json']}))
    except Exception as exc:
        write(folder/'pilot-failure.json',dict(error=repr(exc),completed=rows,seconds=time.monotonic()-start));raise
    finally:signal.alarm(0)


production_owner=FunctionType(mono.production.__code__,dict(vars(mono),OUT=run.RESPONSE,Response=Response))
def production():
    configure();production_owner()


def audit():
    folders=[OUT/'capture-result.json',GR/'result.json',LAPSE/'result.json',run.COLL/'bank-result.json',run.RESPONSE/'result.json']
    results=[json.loads(p.read_text()) for p in folders];assert all(r['passed'] for r in results)
    replay=[json.loads((run.BG/f'replay-{n}.json').read_text()) for n in [64,128]]
    assert all(r['passed'] and max(r['state_response_relative'].values())<1e-12 and max(r['source_relative'].values())<1e-12 for r in replay)
    capture_plan=json.loads((OUT/'capture-execution-plan.json').read_text())
    capture_seconds=capture_plan.get('prior_spent_seconds',0.)+results[0]['seconds'];assert capture_seconds<650
    checked=0
    for p in [GR/'plan.json',LAPSE/'plan.json',run.COLL/'plan.json',run.RESPONSE/'plan.json',run.RESPONSE/'execution-plan.json']:
        for path,h in json.loads(p.read_text())['bindings'].items():assert sha(path)==h,path;checked+=1
    plan=json.loads((OUT/'plan.json').read_text())
    for path,h in plan['bindings'].items():
        target=OUT/'uncached-producer.py' if Path(path).name==Path(run.__file__).name else Path(path)
        assert sha(target)==h,str(target);checked+=1
    assert sha(run.__file__)==json.loads((OUT/'capture-execution-plan.json').read_text())['source_sha256']
    histories=[]
    for n,ref in [(64,128),(128,128),(128,64)]:
        p=run.RESPONSE/f'steps-{n}-reference-{ref}.npz';d=np.load(p)
        ids=[int(np.argmin(abs(d['t']-t))) for t in np.linspace(0,run.flow.old.END,17)]
        assert np.max(abs(d['t'][ids]-np.linspace(0,run.flow.old.END,17)))<1e-18
        assert len(d['t'])==len(d['photon_history_scaled_occupation'])==len(d['material_history'])
        units=np.stack([d['material_energy_units'],d['material_neutral_units']],axis=-1)
        # With no additional mechanical source in this sweep, all integrated
        # energy/H changes are the paired collision transfer by definition.
        error=float(np.max(abs(d['collision_transfer']-d['material_history']*units))/max(np.max(abs(d['collision_transfer'])),1e-300))
        assert error<1e-14
        row=json.loads(p.with_suffix('.json').read_text());assert row['linear_relative']<1e-12
        histories.append(dict(steps=n,reference=ref,canonical_states=17,paired_transfer_relative=error,
            material_energy_erg=row['endpoint_material_reference_energy_erg'],photon_energy_erg=row['endpoint_photon_reference_energy_erg']))
    result=dict(classification='Counterexample candidate',passed=True,bindings_checked=checked,paths=histories,capture_total_production_seconds=capture_seconds,
        maximum_replay_state_response_relative=max(max(r['state_response_relative'].values()) for r in replay),
        maximum_replay_source_relative=max(max(r['source_relative'].values()) for r in replay),
        physical_replay_bitwise=False,physical_replay_registered_relative_gate=True,updated_GR_actually_applied_to_photon_thermal_H=True,
        complete_angular_and_material_response_histories_saved=True,
        additional_material_motion_evolved=False,coupled_fixed_point_verified=False,
        full_EOS_history_error_enclosed=False,discarded_material_transport_closed=False,
        nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
    write(OUT/'audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':globals()[sys.argv[1]]()
