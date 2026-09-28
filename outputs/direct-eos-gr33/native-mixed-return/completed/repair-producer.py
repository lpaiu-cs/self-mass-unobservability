"""Return the computed mixed compact field through actual coupled transport.

Counterexample candidate. The compact-source linear exterior component is
matched; the additional exterior bilinear source is a separate open component.
"""
from pathlib import Path
from types import FunctionType, SimpleNamespace
import inspect, json, math, resource, shutil, sys, time
import numpy as np
from numpy.polynomial import legendre as leg
import couple_native_mixed_gr as mixed
import solve_native_incident_self_gr as reuse
import audit_native_incident_self_gr as audit_owner

OUT=Path('native-mixed-return162-work');METRIC=OUT/'metric'
OLD=reuse.OUT;wave=mixed.wave;inf=mixed.inf;LD=np.longdouble;C=mixed.C;G=mixed.G
read,write,sha=mixed.read,mixed.write,mixed.sha
CAPS=dict(prepare=40,metric=90,check=100,photon_pilot=120,photon_production=1250,
          material_pilot=70,material_production=450,block=300,compact=180,infinity=30,repair=480)
TOTAL=4200


def paths(sweep):
    return OUT/f'sweep-{sweep}/photons',OUT/f'sweep-{sweep}/material'


def prepare():
    assert not OUT.exists();OUT.mkdir();METRIC.mkdir()
    files=[Path(__file__),Path(mixed.__file__),Path(reuse.__file__),Path(reuse.coupled.__file__),
           Path(reuse.readout.__file__),Path(reuse.deep.__file__),Path(audit_owner.__file__),Path(inf.__file__),
           mixed.OUT/'result.json',mixed.OUT/'audit.json',mixed.OUT/'applied-charge.npz',
           inf.OUT/'packet-kernels.npz',inf.BACKGROUND]
    for n,q in mixed.SETTINGS:
        files += [mixed.OUT/f'field-{n}-g{q}.npz',mixed.previous.FIELDS/f'source-{n}.npz']
    for s in range(4):
        for p in paths(s):p.mkdir(parents=True)
    for n in [64,128]:
        mm=OLD/f'sweep-2/material-common-gr/steps-{n}-reference-128.npz'
        pp=OLD/f'sweep-2/photons-common-gr/steps-{n}-reference-128.npz';files += [mm,pp,pp.with_suffix('.json'),mm.with_suffix('.json')]
        with np.load(mm) as d:np.savez_compressed(paths(0)[1]/mm.name,t=d['t'],history_scaled=np.zeros_like(d['history_scaled']))
        with np.load(pp) as p:np.savez_compressed(paths(0)[0]/pp.name,t=p['t'],moments=np.zeros_like(p['moments']),collision_transfer=np.zeros_like(p['collision_transfer']))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b27561096',
        previous_turn='Progress: Phase161 applied the missing compact reciprocal GR source and stored its actual field, preserving the initial-source derivative failure and repair.',
        claim='Match that computed compact-source field to its initial-vacuum exterior component, apply its full scalar/mass/lapse drive to actual photons and free material, close the declared finite block input, then read the resulting charge using actual signed emission.',
        decision='Measure the previously omitted transport feedback of the mixed field. A small input norm alone is not accepted as its final-charge bound.',
        component='Linearity of the initial unknown-field operator separates compact and exterior bilinear forcing. This action returns the compact component. Additional exterior mixed forcing and physical mixed ADM normalization remain open; the auxiliary mass source is not reclassified as emitted energy.',
        boundary='No zero-lapse clamp. Match the generated outgoing scalar leading continuation and the vacuum homogeneous J component with lapse fixed at infinity. Keep the auxiliary asymptotic mass component explicit and do not claim a complete dynamical Einstein exterior.',
        drive='Fine33-knot field is the same input for64/128 transport clocks. Query its piecewise-linear values and left/right interval rates at actual stages through the existing17-knot background owner. No17-knot downsampling of the drive; no original incident affine lift.',
        normalization='One positive power of two scales the separate metric drive to5e-31..1e-30 for arithmetic. Divide all physical return quantities by that exact factor.',
        equations='Existing531cells,8angles,152frequencies,3.434ms,SDIRK photon/totalE-H and SSP2 free matter, analytic deep tangent, extended stage refinement and source-knot alignment. No background replay or new EOS roots.',
        gates=dict(time=.02,quadrature=.002,block=.002,source=1e-12,conservation=1e-8,port=1e-12,
                   stage=1e-12,directional=.002,branch=.01,observer=1e-10),
        convergence='Use the actual lagged B/S/xi/M inputs and stage source defect, with M=(Etilde,H)-collision transfer. Record whole-state changes separately. Maximum3 sweeps; stop if the actual block defect is nondecreasing after sweep2.',
        budget=dict(actions=CAPS,total_action_seconds=TOTAL,CPU_threads=1,virtual_GiB=3,max_sweeps=3),
        forecast='Prior self-GR work took2659s including repairs. Measure equal-horizon4/8 photon prefixes and two-step material prefixes. Reuse them;2x remaining forecast with the completed Phase157 late cost floor must fit1250s photons and450s material. Any third sweep still requires the total4200s admission.',
        stop='Stop on failed gates or admission; preserve failures. No automatic finer clock, mesh, longer interval, extra GR iteration or threshold relaxation.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    (OUT/'registered-producer.py').write_bytes(Path(__file__).read_bytes())


def integrate_centers(q,values,radius):
    """Existing Jordan nodal polynomial, exactly integrated to each center."""
    nc=len(q.h);order=len(q.w);v=values.reshape(len(values),nc,order)
    whole=q.h*(v@q.w);prefix=np.cumsum(whole,axis=1,dtype=LD)-whole
    xx=(radius-q.r[:,0]+q.h*(1+q.x[0]))/q.h-1
    Q=np.column_stack([leg.legval(xx,leg.legint(np.eye(order)[j]))-leg.legval(-1,leg.legint(np.eye(order)[j])) for j in range(order)])
    co=v@np.linalg.inv(leg.legvander(q.x,order-1)).T
    return np.column_stack([prefix+q.h*np.sum(co*Q,axis=2),np.sum(whole,axis=1,dtype=LD)])


def metric_route(n,order):
    folder=OUT/'corrected-fields' if (OUT/'corrected-fields').exists() else mixed.OUT
    d=dict(np.load(mixed.previous.FIELDS/f'source-{n}.npz'));f=np.load(folder/f'field-{n}-g{order}.npz')
    m=wave.Response();m.setup(d,order);t=d['t'];q=wave.base.flow.initial.Quadrature(d['edges'],order)
    r=m.r;z=m.z;a=z['lapse']*np.sqrt(z['b']);jac=1/(np.exp(-2*z['phi']**2)*(1+z['alpha']*r*z['Phi']))
    test=integrate_centers(q,np.ones((1,len(r))),d['radius'])[0]
    exact=np.r_[d['radius'],d['edges'][-1]]-d['edges'][0]
    assert np.max(abs(test-exact))/max(exact)<1e-12
    integral=integrate_centers(q,a*r*f['q_lambda']*jac,d['radius'])
    rc=m.tr;zc=m.tz;J=f['J_centers'] if 'J_centers' in f.files else np.asarray(integral*np.sqrt(zc['b'])/zc['lapse'],float)
    phi=f['delta_phi'];lam=rc*zc['Phi']*phi+J/(rc*zc['b'])
    fn=np.array([np.interp(r,rc,v) for v in phi]);ln=r*z['Phi']*fn+f['J_mixed']/(r*z['b'])
    zr=(2/(r*z['b'])+4*np.pi*r*z['A4']/z['b']*(3*z['Pg']-z['Eg']-z['Kg']-z['Er']+z['R4']))*ln
    zr+=4*np.pi*r*z['A4']/z['b']*z['alpha']*(7*z['Pg']-z['Eg']-3*z['Kg'])*fn+f['q_nu']-f['q_lambda']
    zinteg=integrate_centers(q,zr*jac,d['radius'])
    r0=rc[-1];a0=zc['lapse'][-1]*np.sqrt(zc['b'][-1]);driver=mixed.inc.Driver(order)
    gx,gw=leg.leggauss(order);rr=r0/((gx+1)/2);vac=m.coeff(rr)
    kernel=float(np.sum(gw/(2*r0*vac['lapse']*vac['b']**1.5)))
    # Same leading outgoing scalar continuation as the existing lapse owner.
    outer_scalar=[]
    for i,now in enumerate(t):
        cuts=driver.inverse(C*(now-t[:i+1][::-1]))
        if not i:outer_scalar.append(0.);continue
        radius=(cuts[:-1,None]+np.diff(cuts)[:,None]*(gx+1)/2).ravel();weight=(np.diff(cuts)[:,None]*gw/2).ravel()
        bg=m.coeff(radius);delay=driver.optical(radius)/C
        uu=np.interp(np.clip(now-delay,0,now),t,f['U'][:,-1])
        outer_scalar.append(float(-np.sum(weight*2*bg['Phi']*uu/(radius*bg['b']),dtype=LD)))
    outer_scalar=np.array(outer_scalar);mass=J[:,-1]*zc['lapse'][-1]/np.sqrt(zc['b'][-1])
    outer_J=-mass*(kernel+1/(r0*a0));outer=outer_scalar+outer_J
    zet=np.asarray(outer[:,None]-(zinteg[:,-1,None]-zinteg),float);nu=lam+zet
    scalar_nu=zc['Phi'][-1]*f['U'][:,-1]+outer_scalar;nu_out=scalar_nu-mass*kernel
    edge=float(np.max(abs(nu[:,-1]-nu_out))/max(np.max(abs(nu_out)),1e-290));assert edge<1e-12,edge
    qn_center=f['q_nu_centers'] if 'q_nu_centers' in f.files else np.array([np.interp(rc[:-1],r,v) for v in f['q_nu']])
    volume=3*zc['alpha'][:-1]*phi[:,:-1]+lam[:,:-1]
    dp=-zc['Kg'][:-1]*volume-4*zc['alpha'][:-1]*zc['Pr'][:-1]*phi[:,:-1]-(3*zc['Pr'][:-1]-zc['R4'][:-1])*lam[:,:-1]
    P=zc['Pg'][:-1]+zc['Pr'][:-1];r1=rc[:-1];b=zc['b'][:-1];A4=zc['A4'][:-1]
    nup=(1+8*np.pi*r1*r1*A4*P)*lam[:,:-1]/(r1*b)+4*np.pi*r1*A4/b*(dp+4*zc['alpha'][:-1]*P*phi[:,:-1])
    nup+=r1*zc['Phi'][:-1]*f['delta_Phi'][:,:-1]+qn_center
    alpha=zc['alpha'][:-1];u=alpha*phi[:,:-1];up=-4*zc['Phi'][:-1]*phi[:,:-1]+alpha*f['delta_Phi'][:,:-1]
    data=dict(t=t,radius_E=rc[:-1],delta_phi=phi[:,:-1],delta_u=u,delta_u_prime=up,
        delta_nu=nu[:,:-1],delta_nu_prime=nup,delta_lambda=lam[:,:-1],delta_log_speed=zet[:,:-1],
        delta_log_lapse=nu[:,:-1]+u,delta_log_radial_length=lam[:,:-1]+u,delta_log_areal_radius=u,
        outer_scalar_zeta=outer_scalar,outer_J_zeta=outer_J,auxiliary_asymptotic_mass_cm=mass,
        J_centers=J,interface_relative=edge)
    np.savez_compressed(METRIC/f'metric-{n}-g{order}.npz',**data)
    return data


def metric():
    mixed.previous.initialize();fields=[metric_route(n,q) for n,q in mixed.SETTINGS];fine=fields[0];controls={}
    for name,other in zip(['quadrature','source_cadence'],fields[1:]):
        controls[name]={k:float(np.max(abs((fine[k][::2] if name=='source_cadence' else fine[k])-other[k]))/max(np.max(abs(fine[k])),1e-290)) for k in reuse.METRIC_KEYS}
    mag=max(float(np.max(abs(fine[k]))) for k in reuse.METRIC_KEYS[:4]);power=math.floor(math.log2(1e-30/mag));factor=math.ldexp(1.,power)
    write(OUT/'normalization.json',dict(factor=factor,power_two_exponent=power,physical_metric_maximum=mag,normalized_maximum=factor*mag))
    result=dict(classification='Counterexample candidate',passed=max(controls['quadrature'].values())<.002 and max(controls['source_cadence'].values())<.02,
        controls=controls,maximum_interface_relative=max(float(f['interface_relative']) for f in fields),
        compact_source_exterior_component_matched=True,additional_exterior_mixed_source_closed=False,
        auxiliary_mass_not_reclassified_as_photon_energy=True,full_goal_complete=False)
    write(METRIC/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


def repair():
    failed=read(METRIC/'result.json');assert not failed['passed'];assert not (OUT/'repair-plan.json').exists()
    shutil.copytree(METRIC,OUT/'rejected-metric')
    paths_bound=[Path(__file__),Path(__file__).with_name('resolve_native_mixed_return.py'),OUT/'metric-0-receipt.json']+list((OUT/'rejected-metric').iterdir())
    write(OUT/'repair-plan.json',dict(classification='Counterexample candidate',failure=failed,
        correction='Resolve the original representation: integrate the primary pulse analytically in time; split spatial constraint integrals at existing field centers and pulse fronts; evaluate the constraint at actual material centers. Keep all existing physical cells,33/17 source clocks,4/8 coefficient representations and acceptance gates.',
        field_consistency='Repropagate the affected local primary trace and J-source changes through the same wave/potential operator before constructing the corrected metric. Do not repair only the readout or force table.',
        exclusions='The same additional exterior mixed source and full physical errors remain open. No permission to run transport until the unchanged metric gates pass.',
        budget_seconds=480,total_action_seconds=TOTAL,first_route_admission='Four times the first corrected route plus15s must fit the remaining480s. No automatic increase.',
        original_producer_sha256=sha(OUT/'registered-producer.py'),bindings={str(p):sha(p) for p in paths_bound}))
    import resolve_native_mixed_return
    resolve_native_mixed_return.main();metric()


class StageMetric:
    def __init__(self,n):
        factor=read(OUT/'normalization.json')['factor'];self.factor=factor
        z=dict(np.load(METRIC/'metric-128-g8.npz'));self.clock=z['t']
        self.values={k:v*factor for k,v in z.items() if k.startswith('delta_')}
        self.g={k:v[::2].copy() for k,v in self.values.items()}
    def view(self,t,clock,side='left'):
        t=reuse.canonical_time(t,self.clock)
        k=np.clip(np.searchsorted(self.clock,t,side=side)-1,0,len(self.clock)-2)
        dt=self.clock[k+1]-self.clock[k];w=(t-self.clock[k])/dt
        rows={key:(1-w)*v[k]+w*v[k+1] for key,v in self.values.items()}
        result={key:np.broadcast_to(v,(len(clock),len(v))).copy() for key,v in rows.items()}
        j=np.clip(np.searchsorted(clock,t,side=side)-1,0,len(clock)-2);h=clock[j+1]-clock[j];a=(t-clock[j])/h
        for key in ['delta_u','delta_lambda','delta_nu']:
            rate=(self.values[key][k+1]-self.values[key][k])/dt
            result[key][j]=rows[key]-a*h*rate;result[key][j+1]=rows[key]+(1-a)*h*rate
            result[key+'_interval_rate']=np.broadcast_to(rate,(len(clock)-1,len(rate)))
        return result


def initialize(sweep):
    reuse.OUT=OUT;reuse.METRIC=METRIC;reuse.paths=paths;reuse.SavedMetric=StageMetric
    reuse.initialize(sweep)


def install():
    c=reuse.coupled;c.OUT=OUT;c.paths=paths;c.initialize=initialize;c.CAPS=CAPS
    c.json=SimpleNamespace(dumps=reuse.json_text)
    def cost_read(p):
        p=Path(p)
        if p.name.startswith('steps-') and p.suffix=='.json':
            if p.parent==c.prior.lift.PHOTON:p=OLD/'sweep-2/photons-common-gr'/p.name
            elif p.parent==c.prior.MATERIAL:p=OLD/'sweep-2/material-common-gr'/p.name
        return read(p)
    c.read=cost_read


def check(sweep):
    assert read(METRIC/'result.json')['passed'];install();reuse.check(sweep)
    driver=StageMetric(128);clock=driver.clock[::2];errors=[]
    for t in driver.clock[1:-1]:
        for side in ['left','right']:
            g=driver.view(t,clock,side);j=np.clip(np.searchsorted(clock,t,side=side)-1,0,len(clock)-2);w=(t-clock[j])/(clock[j+1]-clock[j])
            for key in ['delta_u','delta_lambda','delta_nu']:
                got=(1-w)*g[key][j]+w*g[key][j+1];k=int(np.argmin(abs(driver.clock-t)))
                errors.append(float(np.max(abs(got-driver.values[key][k]))/max(np.max(abs(driver.values[key])),1e-290)))
                kk=k-1 if side=='left' else k;rate=np.diff(driver.values[key][kk:kk+2],axis=0)[0]/np.diff(driver.clock[kk:kk+2])[0]
                actual=(g[key][j+1]-g[key][j])/(clock[j+1]-clock[j])
                errors.append(float(np.max(abs(actual-rate))/max(np.max(abs(rate)),1e-290)))
    assert max(errors)<1e-12,errors
    p=OUT/f'sweep-{sweep}/check.json';r=read(p);r['actual33_knot_values_max_relative']=max(errors);write(p,r)


def block(sweep):
    # Reuse the established actual-input/stage audit, without its historical
    # post-failure prerequisites. This criterion is registered before this run.
    s=inspect.getsource(audit_owner.block)
    s=s.replace("OUT/'block-plan.json'", "OUT/f'sweep-{sweep}/block-plan.json'")
    s=s.replace("OUT/'block-result.json'", "OUT/f'sweep-{sweep}/block-result.json'")
    s=s.replace("assert max(list(inputs.values())+source+paired)<.002,rows[-1]", "# A nonconverged first sweep is recorded, not accepted.")
    s=s.replace('passed=True,rows=rows,','passed=maximum<.002,rows=rows,')
    namespace=dict(audit_owner.block.__globals__,OUT=OUT,run=SimpleNamespace(initialize=initialize,coupled=reuse.coupled,paths=paths,AMP=reuse.AMP,C=C))
    path=OUT/f'sweep-{sweep}';assert not (path/'block-plan.json').exists()
    write(path/'block-plan.json',dict(sweep=sweep,bindings={},scope='Actual finite lagged input and stage-source consistency, preregistered in this phase; not a contraction or continuum error bound.'))
    # Store the usual physical ledgers/time comparison before the block audit.
    install();reuse.coupled.residual(sweep)
    exec(compile(s,__file__,'exec'),namespace);namespace['block'](sweep)
    result=read(path/'block-result.json');result.pop('original_stop_preserved');result.pop('original_waveform_change_test_passed')
    before=read(OUT/f'sweep-{sweep-1}/block-result.json')['maximum_block_defect'] if sweep>1 else None
    result.update(criteria_preregistered=True,next_sweep_allowed=not result['passed'] and sweep<3 and (before is None or result['maximum_block_defect']<before))
    write(path/'block-result.json',result)


def compact(sweep):
    install();out=OUT/f'sweep-{sweep}';photon,material=paths(sweep);gr=out/'gr';gr.mkdir()
    assert read(out/'block-result.json')['passed']
    files=[out/'block-result.json',OUT/'normalization.json']+[p/f'steps-{n}-reference-128.npz' for p in [photon,material] for n in [64,128]]
    write(out/'readout-plan.json',dict(bindings={str(p):sha(p) for p in files}))
    def setup():initialize(sweep);reuse.coupled.base.Material=reuse.coupled.Material
    fn=FunctionType(reuse.readout.compact.__code__,dict(reuse.readout.compact.__globals__,OUT=out,PHOTON=photon,MATERIAL=material,GR=gr,fixed=SimpleNamespace(initialize=setup,write=write)))
    fn()


def infinity(sweep):
    folder=OUT/f'sweep-{sweep}';photon,_=paths(sweep);bg=np.load(inf.BACKGROUND);clock=bg['t'];T=clock[-1]
    z=np.load(inf.OUT/'packet-kernels.npz');kernels={tuple(k):v for k,v in zip(z['keys'],z['values'])}
    source=np.load(folder/'gr/source-128-reference-128.npz');M=LD(source['M_cm']);alpha=-LD(source['K_cm'])/M
    factor=LD(read(OUT/'normalization.json')['factor']);e0=bg['epsilon'].astype(LD);q0=bg['normalized'].astype(LD)
    gamma=1-1/np.sqrt(2);paths_out={};rows=[]
    for n,a,r in [(128,8,8),(128,4,8),(128,8,4),(64,8,8)]:
        packets,error=inf.emitted(photon,n,SimpleNamespace(T=T));h=T/n;pieces=[]
        for count in np.arange(17)*(n//16):
            value=np.zeros((2,4),LD)
            for j in range(count):
                for stage,offset in enumerate([gamma,1.]):
                    age=(count-j-offset)*h
                    if age>0:value+=kernels[age,a,r]*packets[2*j+stage]
            pieces.append(np.sum(value,axis=1,dtype=LD))
        pieces=np.asarray(pieces,LD)/factor
        compact=np.load(folder/f'gr/wave-{n}-g8.npz')['free_scalar'].astype(LD)/factor
        deps=LD(G)/LD(C)**4/M*pieces[:,1];mass=(alpha+q0)*deps
        charge=(compact+pieces[:,0]+mass)/(1-e0)
        data=dict(t=clock,compact=compact,exterior=pieces[:,0],arrived_energy_erg=pieces[:,1],mass_numerator=mass,charge=charge)
        np.savez_compressed(OUT/f'return-{n}-a{a}-r{r}.npz',**data);paths_out[n,a,r]=data;rows.append(dict(steps=n,angular=a,radial=r,port_relative=error))
    fine=paths_out[128,8,8];controls={}
    for key in ['compact','exterior','arrived_energy_erg','charge']:
        norm=max(np.max(abs(fine[key])),LD('1e-290'))
        controls[key]={name:float(np.max(abs(paths_out[p][key]-fine[key]))/norm) for name,p in [('time',(64,8,8)),('angular',(128,4,8)),('radial',(128,8,4))]}
    original=read(mixed.OUT/'result.json')['endpoint_reciprocal_compact_return']
    result=dict(classification='Counterexample candidate',passed=all(v['time']<.02 and max(v['angular'],v['radial'])<.002 for v in controls.values()),
        controls=controls,rows=rows,physical_transport_return_endpoint=float(fine['charge'][-1]),return_over_applied_mixed_charge=float(fine['charge'][-1]/original),
        actual_mixed_field_applied_to_photons_and_free_material=True,cached_exterior_kernels_reused=True,
        actual_signed_emission_read_at_fixed_operator_infinity=True,additional_exterior_mixed_source_closed=False,
        physical_final_charge_solved=False,full_goal_complete=False)
    write(OUT/'result.json',result);print(json.dumps(result),flush=True);assert result['passed'],result


if __name__=='__main__':
    action=sys.argv[1];sweep=int(sys.argv[2]) if len(sys.argv)>2 else 0;assert action in CAPS
    receipt=OUT/f'{action}-{sweep}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));mixed.inc.native.deadline(CAPS[action])
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                assert sha(OUT/'registered-producer.py' if p==str(Path(__file__)) else p)==h,p
            if action!='repair' and (OUT/'repair-plan.json').exists():
                for p,h in read(OUT/'repair-plan.json')['bindings'].items():assert sha(p)==h,p
            spent=sum(read(p)['seconds'] for p in OUT.rglob('*-receipt.json'));assert spent+CAPS[action]<=TOTAL,(spent,action)
            if sweep>1:assert read(OUT/f'sweep-{sweep-1}/block-result.json')['next_sweep_allowed']
        if action.startswith(('photon_','material_')):
            install();getattr(reuse.coupled,action.split('_')[0])(sweep,action.endswith('pilot'))
        else:globals()[action](sweep) if sweep else globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,sweep=sweep,seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
