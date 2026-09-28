"""Counterexample candidate: return the reciprocal incident source's own GR.

Keep a separately normalized increment on the existing evolving background.
The finite residual is not a continuum bound or nonlinear Einstein closure.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json,math,resource,sys,time
import numpy as np
import sympy as sp
from scipy.sparse.linalg import gmres
import solve_native_incident_reciprocal as coupled
import read_native_incident_response as readout
import def_native_incident_deep_tangent as deep
import def_native_characteristic_gr as wave
import verify_native_anisotropic_gr as constraints
import def_native_dynamic_lapse as metric

OUT=Path('native-incident-self-gr157-work')
BEFORE=Path('native-incident-reciprocal156-work/sweep-2')
FIELDS=OUT/'fields';METRIC=OUT/'metric'
read,write,sha=coupled.read,coupled.write,coupled.sha
AMP=coupled.AMP;LD=np.longdouble;C=wave.C;G=wave.base.G
CAPS=dict(prepare=30,fields=180,lapse=90,check=90,photon_pilot=120,
          photon_production=1100,material_pilot=60,material_production=400,
          residual=60,compact=180,metric_residual=180)
CAPS.update(linear_plan=30,photon_repair=120)
CAPS.update(state_precision_plan=30,photon_extended_pilot=120)
CAPS['export_plan']=30
original_initialize=coupled.initialize
LINEAR=[]


def json_text(value):
    def scalar(v):
        if isinstance(v,np.bool_):return bool(v)
        if isinstance(v,np.integer):return int(v)
        if isinstance(v,np.floating):return float(v)
        raise TypeError(type(v).__name__)
    return json.dumps(value,default=scalar)


def export_plan():
    pilot=read(paths(1)[0]/'pilot.json');assert pilot['passed'] and pilot['eligible']
    write(OUT/'export-plan.json',dict(classification='Counterexample candidate',source_sha256=sha(__file__),
        error=read(OUT/'photon_extended_pilot-1-receipt.json')['error'],
        repair='Both physical prefixes and their forecast were saved successfully. The subsequent console json.dumps rejected a NumPy bool/longdouble. Convert NumPy scalars only in console serialization, matching the existing saved-JSON writer. Reuse both accepted prefixes byte-for-byte; no physical run is repeated.',
        prefix_forecast_seconds=pilot['upper_remaining_seconds'],production_cap=CAPS['photon_production'],
        bindings={str(p):sha(p) for p in [OUT/'extended-state-producer.py',OUT/'state-precision-plan.json',OUT/'photon_extended_pilot-1-receipt.json',paths(1)[0]/'pilot.json']+
                  [paths(1)[0]/f'pilot-{n}.npz' for n in [64,128]]}))


def refined_gmres(op,rhs,**options):
    args=dict(options,x0=None if options['x0'] is None else np.asarray(options['x0'],float))
    sol,info=gmres(op,np.asarray(rhs,float),**args);original_info=info
    sol=sol.astype(LD)
    for correction in range(4):
        # Same operator, wider residual arithmetic. The correction solve itself
        # remains binary64 with the same sparse/thermochemical preconditioner.
        residual=rhs.astype(LD)-op.matvec(sol.astype(LD))
        relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290))
        if relative<1e-14:
            LINEAR.append(dict(initial_info=int(original_info),corrections=correction,extended_residual=relative))
            return sol,0
        if correction==3:break
        args=dict(options,x0=None,rtol=1e-12)
        delta,_=gmres(op,np.asarray(residual,float),**args)
        sol+=delta
    raise AssertionError(('Extended stage residual',relative,original_info))


def state_precision_plan():
    failure=read(OUT/'photon_repair-1-receipt.json');assert failure['error']
    write(OUT/'state-precision-plan.json',dict(classification='Counterexample candidate',
        failure=failure['error'],source_sha256=sha(__file__),
        repair='The residual-only refinement still rounds each correction into a binary64 stage state. Retain the accumulated solution and SDIRK stage states in longdouble; cast only the Krylov residual equation and initial guess to binary64. Same operator, up to3 corrections and unchanged1e-14/1e-12/1e-8 gates.',
        budget='Same120s prefixes,1100s production admission and4600s total including both failures. No physical/grid/EOS change.',
        bindings={str(p):sha(p) for p in [OUT/'linear-plan.json',OUT/'residual-only-producer.py',OUT/'rejected-residual-only-run.py',OUT/'photon_repair-1-receipt.json']}))


def linear_plan():
    failure=read(OUT/'photon_pilot-1-receipt.json');assert failure['error']
    write(OUT/'linear-plan.json',dict(classification='Counterexample candidate',
        failure=failure['error'],source_sha256=sha(__file__),
        repair='Keep the same finite stage, preconditioner and1e-14 internal target. Form its true residual in longdouble, solve a binary64 correction equation with zero initial guess, and repeat at most3 times. Require extended relative residual<1e-14 and the original binary64 full-equation1e-12 plus physical1e-8 conservation gates.',
        evidence='The first short photon stage stagnated at binary64 residual3.42913e-14 with GMRES info5. No macro prefix or full trajectory completed. Preserve its receipt and expanded owner.',
        budget='Repeat only equal-horizon4/8 prefixes under120s, then require the original1100s forecast admission. Count rejected and repaired actions in4600s.',
        bindings={str(p):sha(p) for p in [OUT/'registered-producer.py',OUT/'rejected-expanded-run.py',OUT/'photon_pilot-1-receipt.json']}))


def paths(sweep):
    p=OUT/f'sweep-{sweep}'
    return p/'photons-precise',p/'material-analytic'


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(BEFORE/'result.json')['passed']
    for p in [FIELDS,METRIC]+[p for i in range(4) for p in paths(i)]:p.mkdir(parents=True)
    files=[Path(__file__),Path(coupled.__file__),Path(coupled.prior.__file__),
           Path(coupled.prior.lift.__file__),Path(coupled.base.__file__),Path(deep.__file__),
           Path(coupled.feedback.__file__),Path(readout.__file__),Path(wave.__file__),
           Path(wave.base.__file__),Path(constraints.__file__),Path(metric.__file__),
           BEFORE/'result.json']
    for n in [64,128]:
        source=BEFORE/f'gr/source-{n}-reference-128.npz'
        (FIELDS/f'source-{n}.npz').write_bytes(source.read_bytes())
        pp=BEFORE/f'photons-precise/steps-{n}-reference-128.npz'
        mm=BEFORE/f'material-analytic/steps-{n}-reference-128.npz'
        files += [source,pp,mm,BEFORE/f'gr/wave-{n}-g8.npz']
        # Exact zero initial waveform for the separate metric-induced increment.
        with np.load(mm) as d:
            np.savez_compressed(paths(0)[1]/f'steps-{n}-reference-128.npz',
                t=d['t'],history_scaled=np.zeros_like(d['history_scaled']))
        with np.load(pp) as p:
            np.savez_compressed(paths(0)[0]/f'steps-{n}-reference-128.npz',
                t=p['t'],moments=np.zeros_like(p['moments']),collision_transfer=np.zeros_like(p['collision_transfer']))
    z,R,H,b,s=sp.symbols('z R H b s',nonzero=True)
    assert sp.expand((z+H)-(R*b+R*H)-(z-R*b+H-R*H))==0
    assert sp.simplify((R*(s*b))/s-R*b)==0
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='2fb64643e',
        claim='Compute scalar/mass/lapse from the accepted reciprocal incident source and its actual signed angular packets; apply that metric to photons and free material, iterate their finite17-knot waveform residual, and read the compact charge correction.',
        decision='Whether the original compact response survives actual self-GR return, and whether the returned source changes the applied self metric by less than the declared finite residual. A metric-only diagnostic is not completion of this stage.',
        separation='Solve a separate zero-initial increment. Use no original incident pulse or affine pulse lift in this increment. Linearity on the fixed saved background permits one final addition to Phase156 only.',
        normalization='After physical lapse construction, multiply all driving metric increments by a single positive power of two so their largest dimensionless metric component is at most1e-30 and greater than5e-31. Solve and audit this normalized finite response; divide every reported physical correction by that exact factor. Do not add sub-ulp increments to a background.',
        equations='Same corrected evolving EOS background,531 cells,8 angles,152 frequencies,3.434ms,64/128 clocks,17 coefficient knots. Etilde coordinate, paired collision/mechanical partition, actual-stage SDIRK and left-closing SSP remain. Deep tangent is analytic; shared HLL/atmospheric directional owner remains.',
        boundary='Only accepted signed Phase156 angular packets source the incremental exterior lapse. Keep scalar-vacuum rays, leading compact outgoing scalar lapse and separately measured ADM residual. No background luminosity double count or outer lapse clamp.',
        gates=dict(time=.02,quadrature=.002,waveform=.002,metric_residual=.002,linear=1e-12,
                   conservation=1e-8,source_identity=1e-12,independent_GR=1e-9,port=1e-12,
                   ray=1e-10,directional=.002,mapping=1e-10,velocity=1e-4,branch=.01),
        budget=dict(per_action_seconds=CAPS,total_action_seconds=4600,CPU_threads=1,virtual_GiB=3,max_material_photon_sweeps=3,max_self_GR_returns=1),
        forecast='Previous two reciprocal sweeps plus readout used1128.78 action seconds. Measure first17-time field and equal-horizon4/8 photon prefixes plus two-step material prefixes. Reuse accepted prefixes;2x forecasts and previous late-step floors must fit1100s photons and400s material. New forcing and late Krylov costs remain extrapolated.',
        stop='Any original gate, nondecreasing reciprocal residual after sweep2, three sweeps, metric residual>=0.002 or4600 action seconds. No automatic new GR return, finer clock/grid/domain, longer horizon, background replay or threshold relaxation.',
        scope='One self-GR return and its actual finite equation residual on retained linear maps. Initial-slice GR operator, full native/continuum derivative and error enclosure, full exterior scalar scattering, orbital/static-observable matching and nonlinear Einstein star remain open.',
        bindings={str(p):sha(p) for p in files}))
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Linearity permits a separated, homogeneously normalized first return. This does not prove a contraction.',
        constraints=wave.base.symbolic(),null_transport=metric.symbolic(),characteristic=wave.check()))


def geometry_initialize(folder):
    coupled.base.old.fields.initialize()
    wave.OUT=folder;wave.base.OUT=folder;constraints.OUT=folder


def gr_fields(folder,source_readout=None):
    geometry_initialize(folder);m=wave.Response();start=time.monotonic()
    fn=FunctionType(wave.Response.run.__code__,dict(wave.Response.run.__globals__,OUT=folder))
    rows=[fn(m,128,8)];first=time.monotonic()-start
    forecast=4*first+10
    write(folder/'pilot.json',dict(classification='Counterexample candidate',first_seconds=first,
        remaining_upper_seconds=forecast,eligible=forecast<CAPS['fields']-first))
    assert forecast<CAPS['fields']-first
    for n,q in [(64,8),(128,4)]:rows.append(fn(m,n,q))
    fine=np.load(folder/'fields-128-g8.npz');controls={}
    for name,n,q in [('time',64,8),('quadrature',128,4)]:
        other=np.load(folder/f'fields-{n}-g{q}.npz')
        controls[name]=float(np.max(abs(fine['U']-other['U']))/max(np.max(abs(fine['U'])),1e-290))
    if source_readout:
        direct=np.load(source_readout)['free_scalar']
        actual=-fine['direct_and_mass_stress_U'][:,-1]/float(m.data['M_cm'])
        controls['compact_reproduction']=float(np.max(abs(actual-direct))/max(np.max(abs(direct)),1e-290))
        assert controls['compact_reproduction']<1e-9
    result=dict(classification='Counterexample candidate',passed=controls['time']<.02 and controls['quadrature']<.002,controls=controls,rows=rows)
    write(folder/'result.json',result);assert result['passed'],result
    return result


def fields():return gr_fields(FIELDS,BEFORE/'gr/wave-128-g8.npz')


class Lapse(metric.Lapse):
    def __init__(self,photon):super().__init__();self.photon=photon
    def boundary(self,d,steps,order):
        t=d['t'];end=t[-1];mu,mw,ids,solution,invariant=self.rays(order,end)
        p=np.load(self.photon/f'steps-{steps}-reference-128.npz');gamma=1-1/np.sqrt(2)
        packets=end/steps*np.tile([1-gamma,gamma],steps)[:,None]*p['accepted_angular_luminosity']
        stage=p['accepted_angular_times'];photon=[];energy=[];norm=[]
        for now in t:
            n=int(round(now/end*steps))*2
            if not n:photon.append(0.);energy.append(0.);norm.append(0.);continue
            weight=packets[:n,ids]*mu*mw
            ages=np.maximum((now-stage[:n])/end,0.)
            state=solution.sol(ages).reshape(3,len(mu),n).transpose(0,2,1)
            rp=self.r0*state[0];mp=state[1]
            cc=self.bg.metric(rp.ravel()/self.model.m.R)[3].reshape(rp.shape)
            kernel=state[2]/self.r0-mp*mp/(rp*cc)
            photon.append(float(G/C**4*np.sum(weight*kernel,dtype=LD)))
            energy.append(float(np.sum(weight,dtype=LD)));norm.append(float(np.sum(abs(weight),dtype=LD)))
        energy=np.array(energy);error=float(np.max(abs(energy-d['outer_cumulative_energy_erg']))/max(max(norm),1e-290))
        assert error<1e-12,error
        gx,gw=np.polynomial.legendre.leggauss(order);z=(gx+1)/2
        _,N,b,_,_=self.bg.metric(self.r0/(z*self.model.m.R));kernel=float(np.sum(gw/(2*self.r0*N*b**1.5)))
        return np.array(photon),energy,dict(ray_invariant=invariant,emitted_energy_relative=error,
            unit_ADM_lapse_kernel_per_cm=kernel,only_actual_signed_incremental_packets=True)


METRIC_KEYS=['delta_u','delta_lambda','delta_log_lapse','delta_log_speed','delta_u_prime','delta_nu_prime']


def gr_lapse(folder,destination,photon):
    assert read(folder/'result.json')['passed'];geometry_initialize(folder)
    destination.mkdir(exist_ok=True);m=Lapse(photon)
    fn=FunctionType(metric.Lapse.run.__code__,dict(metric.Lapse.run.__globals__,OUT=destination))
    rows=[fn(m,n,q) for n,q in [(128,8),(64,8),(128,4)]]
    fine=np.load(destination/'metric-128-g8.npz');controls={}
    for name,n,q in [('time',64,8),('quadrature',128,4)]:
        other=np.load(destination/f'metric-{n}-g{q}.npz')
        controls[name]={k:float(np.max(abs(fine[k]-other[k]))/max(np.max(abs(fine[k])),1e-290)) for k in METRIC_KEYS}
    result=dict(classification='Counterexample candidate',passed=max(controls['time'].values())<.02 and max(controls['quadrature'].values())<.002 and max(r['ray_invariant'] for r in rows)<1e-10,
        controls=controls,rows=rows,actual_metric_applied_to_evolution=False)
    write(destination/'result.json',result);assert result['passed'],result
    return result


def lapse():
    result=gr_lapse(FIELDS,METRIC,BEFORE/'photons-precise')
    d=np.load(METRIC/'metric-128-g8.npz')
    magnitude=max(float(np.max(abs(d[k]))) for k in METRIC_KEYS[:4]);assert magnitude>0
    exponent=math.floor(math.log2(1e-30/magnitude));factor=math.ldexp(1.,exponent)
    write(OUT/'normalization.json',dict(classification='Counterexample candidate',factor=factor,power_two_exponent=exponent,
        physical_metric_maximum=magnitude,normalized_metric_maximum=magnitude*factor,
        reason='Homogeneous arithmetic normalization only; physical metric and reported charge remain unscaled.'))
    return result


class SavedMetric:
    def __init__(self,n):
        self.factor=read(OUT/'normalization.json')['factor']
        self.g={k:(v*self.factor if k.startswith('delta_') else v) for k,v in dict(np.load(METRIC/f'metric-{n}-g8.npz')).items()}
    def view(self,t,clock,side='left'):
        assert np.array_equal(clock,self.g['t'])
        return self.g


def initialize(sweep):
    coupled.OUT=OUT;coupled.paths=paths;original_initialize(sweep)
    Photon,Matter=coupled.Response,coupled.Material
    class Response(Photon):
        def __init__(self,n):
            super().__init__(n);self.driver=SavedMetric(n);self.g=self.driver.g
            self.material.driver=self.driver
        def lift(self,t):
            return np.zeros_like(self.I[0]),np.zeros_like(self.I[0]),np.zeros(3),np.zeros(3),0.
    runner=(OUT/f'sweep-{sweep}/expanded-run.py').read_text()
    runner=coupled.base.old.fields.previous.prior.updated.replace(runner,"    row['passed']=",
        "    row.update(maximum_extended_stage_residual=max(v['extended_residual'] for v in LINEAR),maximum_refinement_steps=max(v['corrections'] for v in LINEAR),initial_GMRES_stagnations=sum(v['initial_info']!=0 for v in LINEAR))\n    row['passed']=")
    scope=dict(Photon.run.__globals__,gmres=refined_gmres,LINEAR=LINEAR)
    exec(compile(runner,__file__,'exec'),scope);Response.run=scope['run']
    (OUT/f'sweep-{sweep}/expanded-refined-run.py').write_text(runner)
    class Material(Matter):
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps);self.driver=SavedMetric(steps)
        deep_tangent=deep.deep_tangent
        directional_rhs=deep.rhs(coupled.base.old.aligned)[0]
    coupled.Response=Response;coupled.Material=Material


def check(sweep):
    coupled.check(sweep)
    # Reuse the already initialized object graph for the zero/odd/scaling test.
    m=coupled.Response(128);rows=[]
    for t in [m.t[-1]*.31,m.t[-1]*.63]:
        samples={}
        for scale in [0.,1.,-1.,2.]:
            m.drive_scale=scale;s,l,e=m.source(t)
            assert e<1e-12;samples[scale]=[s,l]
        assert all(np.count_nonzero(v)==0 for v in samples[0.])
        error=max(float(np.max(abs(samples[a][j]-a*samples[1.][j]))/max(np.max(abs(samples[1.][j])),1e-290)) for a in [-1.,2.] for j in range(2))
        assert error<1e-12;rows.append(dict(time=t,zero_exact=True,odd_scaling_relative=error))
    p=OUT/f'sweep-{sweep}/check.json';r=read(p);r['metric_transport_controls']=rows;write(p,r)


def compact(sweep):
    out=OUT/f'sweep-{sweep}';photon,material=paths(sweep);gr=out/'gr'
    r=read(out/'residual.json');assert r['passed'] and r['finite_reciprocal_residual_passed'];gr.mkdir()
    files=[Path(__file__),out/'residual.json',OUT/'normalization.json']
    files += [p/f'steps-{n}-reference-128.npz' for p in [photon,material] for n in [64,128]]
    write(out/'readout-plan.json',dict(classification='Counterexample candidate',bindings={str(p):sha(p) for p in files},
        scope='Normalized self-metric response only, physical correction is divided by normalization.factor.'))
    def setup():initialize(sweep);coupled.base.Material=coupled.Material
    fn=FunctionType(readout.compact.__code__,dict(readout.compact.__globals__,OUT=out,PHOTON=photon,MATERIAL=material,GR=gr,fixed=SimpleNamespace(initialize=setup,write=write)))
    fn();result=read(gr/'result.json');factor=read(OUT/'normalization.json')['factor']
    baseline=read(BEFORE/'result.json')['compact_return_endpoint'];correction=result['compact_return_endpoint']/factor
    result.update(normalized_compact_endpoint=result['compact_return_endpoint'],normalization_factor=factor,
        physical_self_GR_correction_endpoint=correction,baseline_compact_endpoint=baseline,
        signed_correction_over_baseline=correction/abs(baseline),total_compact_endpoint=baseline+correction,
        finite_material_photon_residual=r['maximum_residual'],actual_self_GR_returned=True,
        full_goal_complete=False,full_null_infinity_charge=False)
    write(out/'result.json',result);return result


def metric_residual(sweep):
    out=OUT/f'sweep-{sweep}';gr=out/'gr';returned=out/'returned-fields';returned.mkdir()
    for n in [64,128]:(returned/f'source-{n}.npz').write_bytes((gr/f'source-{n}-reference-128.npz').read_bytes())
    gr_fields(returned,gr/'wave-128-g8.npz');dest=out/'returned-metric'
    gr_lapse(returned,dest,paths(sweep)[0]);factor=read(OUT/'normalization.json')['factor'];rows=[]
    for n in [64,128]:
        before=np.load(METRIC/f'metric-{n}-g8.npz');after=np.load(dest/f'metric-{n}-g8.npz')
        relative={k:float(np.max(abs(after[k]))/factor/max(np.max(abs(before[k])),1e-290)) for k in METRIC_KEYS}
        rows.append(dict(steps=n,returned_over_applied_metric=relative))
    maximum=max(max(r['returned_over_applied_metric'].values()) for r in rows)
    result=read(out/'result.json');result.update(metric_residual_rows=rows,maximum_self_metric_residual=maximum,
        finite_self_metric_residual_passed=maximum<.002,passed=result['passed'] and maximum<.002,
        full_nonlinear_fixed_point=False,complete_reciprocal_fixed_point=False,
        scope='Actual self-GR application with finite material/photon and metric residuals; not a contraction, continuum error bound or final null-infinity charge.')
    write(OUT/'result.json',result);assert result['passed'],result
    return result


if __name__=='__main__':
    action=sys.argv[1];sweep=int(sys.argv[2]) if len(sys.argv)>2 else 0
    assert action in CAPS
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));coupled.base.drive.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None;receipt=OUT/f'{action}-{sweep}-receipt.json'
    assert not receipt.exists()
    try:
        if action not in ['prepare','linear_plan','state_precision_plan','export_plan']:
            plan=read(OUT/'plan.json')
            for p,h in plan['bindings'].items():
                target=OUT/'registered-producer.py' if Path(p).name==Path(__file__).name else p
                assert sha(target)==h,p
            repair=read(OUT/'linear-plan.json');assert sha(OUT/'residual-only-producer.py')==repair['source_sha256']
            revision=read(OUT/'state-precision-plan.json');assert sha(OUT/'extended-state-producer.py')==revision['source_sha256']
            export=read(OUT/'export-plan.json');assert sha(__file__)==export['source_sha256']
            for p,h in dict(repair['bindings'],**revision['bindings'],**export['bindings']).items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.rglob('*-receipt.json'))+CAPS[action]<=plan['budget']['total_action_seconds']
            if sweep>1:assert read(OUT/f'sweep-{sweep-1}/residual.json')['next_sweep_allowed']
        coupled.OUT=OUT;coupled.paths=paths;coupled.initialize=initialize;coupled.CAPS=CAPS
        coupled.json=SimpleNamespace(dumps=json_text)
        if action.startswith(('photon_','material_')):
            getattr(coupled,action.split('_')[0])(sweep,action.endswith(('pilot','repair')))
        elif action=='residual':coupled.residual(sweep)
        else:
            result=globals()[action](sweep) if sweep else globals()[action]()
            if result is not None:print(json.dumps(result),flush=True)
    except Exception as exc:error=repr(exc);raise
    finally:
        write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
