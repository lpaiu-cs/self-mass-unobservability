"""Counterexample candidate: continue the accepted, full-input joint solution.

Only restart/history plumbing changes. Frozen183 equations, fronts and gates
remain the owners. Old128 bytes are explicitly adapted to the equivalent64
prefix; they are never overwritten or treated as a new physical solution.
"""
from pathlib import Path
from types import FunctionType
import gc,json,os,resource,sys,time
import numpy as np
import resolve_full_material_front as prior

OUT=Path('native-front-continuation184-work');OLD=prior.OUT
owner=prior.prior;LD=owner.LD;AMP=owner.AMP
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=20,check=60,coarse=350,fine=550,audit=20)
HISTORY=dict(floor_discard_history='material_floor_discard_history_scaled',
    conserved_history='conserved_material_history',stage_t='joint_stage_times',
    stage_h='joint_stage_weights',stage_states='joint_stage_conserved_scaled',
    stage_native='joint_native_rates_scaled',stage_collision='joint_collision_rates_scaled',
    stage_discard='joint_discard_rates_scaled')


def paths(n,label='pilot'):
    return OUT/f'sweep-1/photons/{label}-{n}.npz'


def state(m):
    return dict(material_floor_discard_scaled=m.floor_discard,
        **{key:np.array(getattr(m,attr)) for attr,key in HISTORY.items()},
        joint_checks_json=np.array(json.dumps(dict(newton=m.newton_iterations,stages=m.stage_log))),
        restart_guide=m.guide_g)


def restore(m,z,row,steps):
    assert row['steps']==row['base_steps']==steps
    flags,_,_=prior.split_flags(m,steps);begin=row['completed_steps']
    assert abs(float(z['t'][-1])-begin*m.t[-1]/steps)<1e-18
    assert len(z['actual_step_edges'])-1==sum(1+flags[:begin])
    for attr,key in HISTORY.items():setattr(m,attr,list(z[key]))
    m.floor_discard=z['material_floor_discard_scaled'].copy()
    checks=json.loads(str(z['joint_checks_json']))
    m.newton_iterations=checks['newton'];m.stage_log=checks['stages']
    m.guide_g=z['restart_guide'].copy()
    assert len(m.stage_t)==len(m.stage_h)==len(m.stage_log)==2*len(m.newton_iterations)
    assert np.array_equal(z['joint_stage_times'],z['accepted_angular_times'])
    assert np.array_equal(z['joint_stage_weights'],z['accepted_angular_quadrature_weights'])
    for attr,key in [('velocity_jet_error','velocity_jet_relative'),('mapping_error','primitive_mapping_relative')]:
        setattr(m,attr,max(getattr(m,attr),row[key]))


def initialize():
    prior.OUT=OUT;prior.initialize();parent=owner.Model.run
    source=(OUT/'expanded-joint-run.py').read_text()
    changes=[
        ('    count=steps if limit is None else limit','    previous={}\n    count=steps if limit is None else limit'),
        ('    if restart is not None:\n        H,_,e,_,_=self.lift(begin*base_h)',
         '    if restart is not None:\n        restore(self,z,previous,steps)\n'
         "        moment=previous['frequency_moment_relative']\n"
         "        if 'restart_x' in z:\n"
         "            x=z['restart_x'].copy();g=z['restart_g'].copy();ledger=z['restart_ledger'].copy();escape=z['restart_escape'].copy()\n"
         "            impulse=z['restart_impulse'].copy();ports=z['restart_ports'].copy();transfer=z['restart_transfer'].copy()\n"
         '        H,_,e,_,_=self.lift(begin*base_h)'),
        ("energy_offset_t=times,conserved_material_history=self.conserved_history,material_floor_discard_scaled=self.floor_discard,material_floor_discard_history_scaled=self.floor_discard_history,",
         'energy_offset_t=times,**state(self),restart_x=x,restart_g=g,restart_ledger=ledger,restart_escape=escape,restart_impulse=impulse,restart_ports=ports,restart_transfer=transfer,'),
        ('transfer=transfer,accepted_angular_times=', 'transfer=transfer,**state(self),accepted_angular_times='),
        ("maximum_extended_stage_residual=max(v['extended_residual'] for v in LINEAR),maximum_refinement_steps=max(v['corrections'] for v in LINEAR),initial_GMRES_stagnations=sum(v['initial_info']!=0 for v in LINEAR)",
         "maximum_extended_stage_residual=max(previous.get('maximum_extended_stage_residual',0.),max((v['extended_residual'] for v in LINEAR),default=0.)),maximum_refinement_steps=max(previous.get('maximum_refinement_steps',0),max((v['corrections'] for v in LINEAR),default=0)),initial_GMRES_stagnations=previous.get('initial_GMRES_stagnations',0)+sum(v['initial_info']!=0 for v in LINEAR)")]
    for a,b in changes:
        assert source.count(a)==1,(a,source.count(a));source=source.replace(a,b)
    ns=dict(parent.__globals__,state=state,restore=restore)
    assert not ns['LINEAR'],'One producer initialization per process'
    exec(compile(source,__file__,'exec'),ns);owner.Model.run=ns['run']
    (OUT/'expanded-continuation-run.py').write_text(source)


def prepare():
    assert not OUT.exists();OUT.mkdir();reuse={}
    assert read(OLD/'time-result.json')['passed'] and read(OLD/'audit-result.json')['passed']
    for s in [0,1]:
        for name in ['photons','material']:(OUT/f'sweep-{s}/{name}').mkdir(parents=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/name for name in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for src in files:
        dst=OUT/src.relative_to(OLD);os.link(src,dst);reuse[str(dst.relative_to(OUT))]=dict(path=str(src),sha256=sha(src))
    for n in [64,128]:
        src=OLD/f'sweep-1/photons/pilot-{n}.npz';p=dict(np.load(src));row=read(src.with_suffix('.json'))
        checks=read((prior.OLD if n==64 else OLD)/'stage-checks.json')
        # Seed only: last pre-floor state gives the same previous closing-stage
        # guide as the uninterrupted owner; physical end state is untouched.
        z=p['joint_stage_conserved_scaled'][-1];b=z[0];s=z[1]
        # Units/kappa are resolved during check with the unchanged constructor.
        p['joint_checks_json']=np.array(json.dumps(checks))
        np.savez_compressed(paths(n,'input'),**p)
        row.update(steps=n,base_steps=n,completed_steps=n//32,
            split_completed_macro_steps=n//32,adapted_from=str(src),adapted_sha256=sha(src),
            physical_step_sequence_unchanged=True)
        write(paths(n,'input').with_suffix('.json'),row)
        files.extend([src,src.with_suffix('.json'),(prior.OLD if n==64 else OLD)/'stage-checks.json'])
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None)
        and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    files += [Path(__file__),OLD/'equivalence.json',OLD/'run-receipt.json',OLD/'time-result.json',OLD/'audit-result.json',owner.drive.FIELDS/'born-g8.npz']
    write(OUT/'reuse.json',reuse)
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='12ffe0846',
        previous_turn='The last continuation only revalidated and reported183; no new authoritative progress. This trial fixes the missing restart histories and applies it to actual accepted states.',
        claim='Continue the accepted full incident photon/B/S/Etilde/H solution without losing its own material stages, floor removal, angular ports or solver history. Then decide the original time gates at the first canonical T/16 output.',
        reuse='No physical prefix is recomputed. Original18364alias is old182128 with exactly equivalent actual stages, explicitly adapted to64completed2macro. Fine128completed4macro remains128. Every original physical array remains unchanged.',
        method='Frozen183time splits, same EOS, background, full incident driver, physical amplitude, equations, floor map, Radau method and original gates. Only restoration, serialization and diagnostic accumulation change.',
        sequence='Zero-step roundtrip both saved paths first, including all old physical arrays and new restart payload. Only after success continue64from2to4macro and128from4to8macro; compare all channels and audit same-solution local ledgers.',
        forecast=dict(coarse_seconds_range=[190,350],fine_seconds_range=[260,550],
            basis='Old coarse-equivalent4steps196.934s and new fine8steps246.395s. Next intervals are unmeasured: coarse1.5x measured plus40s=335.4s; fine1.75x plus40s=471.2s. These are assumptions, not certified bounds. Check actual split counts before dispatch.'),
        gates=dict(time=.02,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,max_Newton_solves=3),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        stop='Any restore, original stage/time/balance/constitutive gate, iteration or cap failure stops this trial. No automatic new clock, second split, further horizon or GR sweep.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def check():
    initialize();rows=[]
    for n in [64,128]:
        m=owner.Model(n);p=dict(np.load(paths(n,'input')));q=p['joint_stage_conserved_scaled'][-1]
        p['restart_guide']=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])
        np.savez_compressed(paths(n,'input'),**p)
        row=m.run(n,f'roundtrip-{n}',n//32,f'input-{n}');out=dict(np.load(paths(n,'roundtrip')))
        src=dict(np.load(OLD/f'sweep-1/photons/pilot-{n}.npz'))
        rebuilt=['delta_packet_scaled_occupation','delta_material','ledger','escape']
        errors={}
        for key,value in src.items():
            if key in ['split_macro_steps','front_cells','front_times']:continue
            if key in rebuilt:
                err=float(np.max(abs(out[key]-value))/max(np.max(abs(value)),LD('1e-290')))
                assert err<4*np.finfo(LD).eps,(key,err);errors[key]=err
            else:assert prior.fingerprint(out[key])==prior.fingerprint(value),key
        assert json.loads(str(out['joint_checks_json']))==json.loads(str(p['joint_checks_json']))
        old=read(paths(n,'input').with_suffix('.json'))
        for key in ['maximum_extended_stage_residual','maximum_refinement_steps','initial_GMRES_stagnations','frequency_moment_relative']:
            assert row[key]==old[key],key
        assert row['new_steps']==row['actual_new_steps']==0 and row['passed']
        flags,_,_=prior.split_flags(m,n);extra=sum(1+flags[n//32:n//16]);assert extra==n//16
        # All continuation uses this serialized roundtrip with exact internal
        # variables. A second zero-step pass must preserve those bits too.
        del m;gc.collect();m=owner.Model(n)
        second=m.run(n,f'identity-{n}',n//32,f'roundtrip-{n}')
        with np.load(paths(n,'identity')) as s:
            for key in out:assert prior.fingerprint(s[key])==prior.fingerprint(out[key]),('Second roundtrip',key)
        rows.append(dict(clock=n,unchanged_physical_prefix=True,zero_physical_steps=True,
            second_serialization_exact=True,legacy_endpoint_roundtrip_relative=errors,new_actual_substeps=int(extra)))
        del m;gc.collect()
    from fractions import Fraction as F
    for k in range(3):assert F(3,4)*F(1,3)**k+F(1,4)==F(1,k+1)
    write(OUT/'check-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        symbolic=dict(classification='Proven',passed=True,scope='Exact rational Radau moments0..2; restart does not change the stage rule.'),
        final_charge_conclusion='unadjudicated'))


def continue_path(n):
    assert read(OUT/'check-result.json')['passed'];initialize()
    pilot,source=owner.clone(owner.pilot,[('initialize();m=Model(64)',f'm=Model({n})'),
        ("m.run(64,'pilot-64',2)",f"m.run({n},'pilot-{n}',{n//16},'roundtrip-{n}')"),
        ("paths(1)[0]/'pilot-64.npz'",f"paths(1)[0]/'pilot-{n}.npz'")])
    (OUT/f'expanded-dispatch-{n}.py').write_text(source);pilot()
    os.rename(OUT/'pilot-result.json',OUT/f'pilot-{n}-result.json')
    os.rename(OUT/'stage-checks.json',OUT/f'stage-checks-{n}.json')


def audit():
    ns=dict(prior.before.compare.__globals__,OUT=OUT)
    compare=FunctionType(prior.before.compare.__code__,ns)
    try:compare()
    finally:
        audit=FunctionType(prior.before.audit.__code__,dict(prior.before.audit.__globals__,OUT=OUT));audit()
        rows=[]
        for n in [64,128]:
            src=dict(np.load(paths(n,'roundtrip')));new=dict(np.load(paths(n)))
            for key in list(HISTORY.values())+['t','moments','radial_ports','photon_history_scaled_occupation','material_history','collision_transfer','accepted_angular_times','accepted_angular_luminosity','accepted_angular_quadrature_weights','actual_step_edges','energy_offset_reference','energy_offset_t']:
                assert prior.fingerprint(new[key][:len(src[key])])==prior.fingerprint(src[key]),key
            row=read(paths(n).with_suffix('.json'));assert row['new_steps']==n//32 and row['actual_new_steps']==n//16
            checks=read(OUT/f'stage-checks-{n}.json');old=json.loads(str(src['joint_checks_json']))
            for key in old:assert checks[key][:len(old[key])]==old[key]
            rows.append(dict(clock=n,all_saved_prefixes_preserved=True,new_macro_steps=row['new_steps'],new_actual_substeps=row['actual_new_steps']))
        write(OUT/'prefix-audit.json',dict(classification='Counterexample candidate',passed=True,rows=rows))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in ['coarse','fine']:continue_path(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
