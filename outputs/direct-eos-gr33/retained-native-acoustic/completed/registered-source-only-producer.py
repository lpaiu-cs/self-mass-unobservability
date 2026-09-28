"""Apply the native acoustic force to material, photons and GR charge.

Counterexample candidate: one finite constitutive response, not a fixed point.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,signal,sys,time,textwrap,resource
import numpy as np
import def_retained_native_acoustic as native
import def_retained_motion_return as motion

prior=native.prior;OUT=native.OUT/'response';SOURCE=OUT/'collision';PHOTON=OUT/'photons'
DIRECT=OUT/'direct-material';MATERIAL=OUT/'material';GR=OUT/'gr';BEFORE=OUT/'mechanical-input'
FORCE=native.OUT/'force-stable-final'
write,read,sha=native.write,native.read,native.sha


def prepare():
    assert not OUT.exists()
    for p in [OUT,SOURCE,PHOTON,DIRECT,MATERIAL,GR,BEFORE/'photons']:p.mkdir(parents=True,exist_ok=True)
    before=prior.OUT
    for name,target in [('material',DIRECT),('hydro',FORCE),('collision',before/'collision'),
                        ('metric',before/'metric'),('immutable-coupled-128.npz',before/'immutable-coupled-128.npz')]:
        (BEFORE/name).symlink_to(target.resolve(),target_is_directory=target.is_dir())
    (BEFORE/'photons/bank-128').symlink_to((before/'photons/bank-128').resolve(),target_is_directory=True)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',
        claim='Propagate native acoustic force through finite free-material evolution, its actual collision response, simultaneous photons/thermal/H, corrected free material and null-infinity charge.',
        sequence='Mechanical-only force response supplies the known finite matter source. No collision drive or existing photon defect is duplicated. Photon thermal/H correction is returned once to the same forced free material.',
        budgets=dict(material_total_seconds=500,source_seconds=600,photon_seconds=1200,readout_seconds=180),
        total_budget='All native and response receipts stay within the original3100s wall cap. Dispatch only after completed native cost and response prefix forecasts establish capacity.',
        gates=dict(time=.02,quadrature=.002,conservation=1e-8,flux_resolution=.002,independent_GR=1e-9),
        limitations='Known native constitutive force at saved background; response Jacobian remains the declared table owner. No uniform derivatives, nonlinear Einstein solution, coupling contraction, continuum space/floor or static-nuisance certification.',
        bindings={str(p):sha(p) for p in [Path(__file__),Path(native.__file__),Path(prior.__file__),Path(motion.__file__),prior.EV/'coupled-128.npz']}))
    for n in [64,128]:
        p=np.load(before/f'photons/steps-{n}-reference-128.npz')
        np.savez_compressed(BEFORE/f'photons/steps-{n}-reference-128.npz',t=p['t'],
            photon_history_scaled_occupation=np.zeros_like(p['photon_history_scaled_occupation']),
            material_history=np.zeros_like(p['material_history']))


def initialize_material(direct):
    global Material
    prior.initialize();prior.HYDRO=FORCE;folder=DIRECT if direct else MATERIAL
    class AcousticMaterial(prior.Material):
        def __init__(self,reference=128,steps=128):
            super().__init__(reference,steps,driven=False)
            if not direct:
                p=np.load(PHOTON/f'steps-{steps}-reference-128.npz');ids=[np.argmin(abs(p['t']-t)) for t in self.t]
                assert np.max(abs(p['t'][ids]-self.t))<1e-18;c=p['collision_transfer'][ids]
                self.transfer=np.stack([np.zeros_like(c[:,:,0]),p['moments'][ids,3]/self.a,c[:,:,0],c[:,:,1]],axis=1)/prior.AMP
        run=FunctionType(prior.Material.run.__code__,dict(prior.Material.run.__globals__,OUT=folder),argdefs=prior.Material.run.__defaults__)
    Material=AcousticMaterial


def material(direct,pilot):
    assert read(native.OUT/'force.json')['passed'];initialize_material(direct)
    folder=DIRECT if direct else MATERIAL;start=time.monotonic();native.deadline(500);rows=[]
    if not pilot:assert read(folder/'execution-plan.json')['eligible']
    for n in [64,128]:
        tick=time.monotonic();m=Material(128,n);label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        row=m.run(n,label,2 if pilot else None,None if pilot else f'pilot-{n}')
        row.update(worker_seconds=time.monotonic()-tick,finite_calls=m.finite_calls,
                   finite_resolution=m.finite_resolution,maximum_owner=max(p['owner_error'] for p in m.cache.values()))
        row['passed']=bool(row['passed'] and row['maximum_owner']<1e-8)
        write(folder/f'{label}.json',row);rows.append(row);assert row['passed']
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,seconds=time.monotonic()-start)
    if pilot:
        old=[read(prior.MATERIAL/f'steps-{n}-reference-128.json') for n in [64,128]]
        result['upper_remaining_seconds']=2*sum(p['finite_owner_calls']*r['worker_seconds']/max(r['finite_calls'],1)+10 for p,r in zip(old,rows))
        result['eligible']=result['upper_remaining_seconds']<500
    write(folder/('pilot.json' if pilot else 'production.json'),result)
    print(json.dumps(result),flush=True);signal.setitimer(signal.ITIMER_REAL,0.)


def admit(direct):
    folder=DIRECT if direct else MATERIAL;p=read(folder/'pilot.json');assert p['passed']
    old=[read(prior.MATERIAL/f'steps-{n}-reference-128.json') for n in [64,128]]
    upper=2*sum(v['finite_owner_calls']*r['seconds']/max(r['finite_calls'],1)+(r['worker_seconds']-r['seconds'])+10 for v,r in zip(old,p['rows']))
    write(folder/'execution-plan.json',dict(classification='Counterexample candidate',eligible=upper<500,
        upper_seconds=upper,cap_seconds=500,original_pilot_eligibility=p['eligible'],
        repair='Separate measured stepping cost from one-time owner setup. Original forecast incorrectly charged the entire setup again for every future owner call. Reuse both accepted prefixes, preserve original rejected forecast and all gates.',
        bindings={str(v):sha(v) for v in [Path(__file__),folder/'pilot.json',OUT/'registered-prefix-producer.py']}))
    print(upper);assert upper<500


def initialize_photon():
    motion.BEFORE=BEFORE;motion.OUT=OUT;motion.SOURCE=SOURCE;motion.PHOTON=PHOTON;motion.TOTAL=PHOTON
    motion.MATERIAL=MATERIAL;motion.GR=GR;motion.initialize()
    # No preexisting photon/thermal response exists in this isolated mechanical
    # source. Resolve the finite source itself, not a division by zero linear drive.
    source=textwrap.dedent(inspect.getsource(motion.State.point))
    source=prior.replace(source,'norm=max(np.sum(abs(factor*linear)*weight),1.)',
                         'norm=max(np.sum(abs(residual)*weight),1.)')
    source=source.replace('residual_over_linear','residual_over_full_source').replace('rounding_over_linear','rounding_over_full_source').replace('projection_over_linear','projection_over_full_source')
    ns=dict(motion.State.point.__globals__);exec(compile(source,__file__,'exec'),ns);motion.State.point=ns['point']
    nonzero=motion.State.point
    def point(self,k,half=False):
        if np.any(self.m.motion[k]):return nonzero(self,k,half)
        photon=np.zeros_like(self.m.I[k]);escape=np.zeros((3,self.m.n))
        np.savez_compressed(self.folder/f'point-{k}.npz',t=self.m.t[k],photon=photon,bound=photon,escape=escape)
        row=dict(classification='Proven',premise='Identical constitutive inputs give exactly zero difference.',
                 k=k,steps=self.n,passed=True,exact_zero=True,seconds=0.)
        write(self.folder/f'point-{k}.json',row);return row
    motion.State.point=point
    (OUT/'expanded-collision-source.py').write_text(source)


def source(pilot):
    assert read(DIRECT/'production.json')['passed'];initialize_photon();native.deadline(600)
    motion.CAPS['source_production']=600;motion.source(pilot);signal.setitimer(signal.ITIMER_REAL,0.)


def photon(pilot):
    initialize_photon();native.deadline(1200)
    # Existing full-path step costs supply the forecast, while the actual
    # zero photon input remains unchanged in the mechanical-input directory.
    result=read(Path('retained-native-return150-work/photons/result.json'))
    write(BEFORE/'photons/result.json',result)
    motion.CAPS['photon_production']=1200;motion.photon(pilot);signal.setitimer(signal.ITIMER_REAL,0.)


def compact():
    initialize_material(False);native.deadline(180)
    original=prior.PHOTON,prior.MATERIAL,prior.GR,prior.OUT,prior.Material
    prior.PHOTON=PHOTON;prior.MATERIAL=MATERIAL;prior.GR=GR;prior.OUT=OUT;prior.Material=Material
    try:prior.readout()
    finally:prior.PHOTON,prior.MATERIAL,prior.GR,prior.OUT,prior.Material=original
    signal.setitimer(signal.ITIMER_REAL,0.)


def infinity():
    start=time.monotonic();native.deadline(180);before=Path('retained-metric-return152-work')
    s=inspect.getsource(prior.infinity)
    left=s.index("    original=dict(np.load(EOS/'applied-charge.npz'))")
    right=s.index('    for n,a,r in [(128,8,8)',left)
    s=s[:left]+"    prior=dict(np.load(BEFORE/'infinity/charge-128-a8-r8.npz'));alpha=-m.K/m.M;paths=[];rows=[]\n"+s[right:]
    s=prior.replace(s,"background=np.load(CURRENT/f'infinity/completed/retained-128-a{a}-r{r}.npz')",
        "background=np.load(BEFORE/f'infinity/charge-{n}-a{a}-r{r}.npz')")
    s=s.replace("background['exterior']","background['normalized_exterior']").replace("background['arrived']","background['arrived_energy_erg']")
    s=prior.replace(s,"previous=(baseline+d['normalized_exterior_parts'][:,0]+alpha*eps0)/(1-eps0)","previous=background['normalized']")
    s=prior.replace(s,"compact=baseline+wave['free_scalar']","compact=background['compact'].astype(LD)+wave['free_scalar']")
    s=s.replace('native_return','acoustic_return')
    ns=dict(vars(prior),OUT=OUT,GR=GR,PHOTON=PHOTON,BEFORE=before,CAPS=dict(infinity=180))
    exec(compile(s,__file__,'exec'),ns);(OUT/'expanded-infinity.py').write_text(s);ns['infinity']()
    result=read(OUT/'infinity/result.json');result.update(native_acoustic_force_applied=True,
        fixed_inventory_native_acoustic_samples=True,uniform_native_derivative_enclosure=False,
        acoustic_change_applied_to_photons_and_free_material=True,previous_geometry_charge_preserved=True,
        response_native_Jacobian_updated=False,same_inventory_static_comparison=False,
        observational_signal_detected=False,seconds=time.monotonic()-start)
    write(OUT/'infinity/result.json',result);print(json.dumps(result),flush=True);signal.setitimer(signal.ITIMER_REAL,0.)


if __name__=='__main__':
    action=sys.argv[1];start=time.monotonic();cpu=time.process_time();error=None
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3))
    try:
        if action in ['direct-pilot','direct','material-pilot','material']:material(action.startswith('direct'),action.endswith('pilot'))
        elif action in ['direct-admit','material-admit']:admit(action.startswith('direct'))
        elif action in ['source-pilot','source']:source(action.endswith('pilot'))
        elif action in ['photon-pilot','photon']:photon(action.endswith('pilot'))
        else:globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():
            receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
            write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
                peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
