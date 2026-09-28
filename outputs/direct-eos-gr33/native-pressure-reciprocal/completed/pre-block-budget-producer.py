"""Close the material and finite block tests of the registered175 return.

Counterexample candidate. Consume the completed actual Radau photons once;
reuse the174 material operator and original gates, without another iteration.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,resource,sys,time
import numpy as np
import return_native_pressure_reciprocity as prior
import audit_native_incident_self_gr as audit

run=prior.run;OUT=prior.OUT;paths=prior.paths
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=10,material_pilot=60,material_production=450,residual=30,block=300)


def material_source():
    source=inspect.getsource(run.material)
    old="""            per_call=max(row['seconds']/row['raw_owner_calls'],old['seconds']/old['raw_owner_calls'])
            forecasts.append(max(0,old['raw_owner_calls']-row['raw_owner_calls'])*per_call+max(0,wall-row['seconds'])+5)"""
    new="""            per_step=max(row['seconds']/row['substeps'],old['seconds']/old['substeps'])
            setup=max(old['worker_wall_seconds']-old['seconds'],wall-row['seconds'])
            forecasts.append(max(0,old['substeps']-row['substeps'])*per_step+setup+5)"""
    assert source.count(old)==1;source=source.replace(old,new)
    old='nominal_arithmetic_probe=16.';assert source.count(old)==1
    return source.replace(old,'constitutive_log_step=1e-5,gross_flux_difference_used=False,forward_probe_indicator=None')


def block_source():
    source=inspect.getsource(audit.block)
    def change(a,b):
        nonlocal source
        assert source.count(a)==1,(a,source.count(a));source=source.replace(a,b)
    change("p['collision_transfer'].transpose(0,2,1)/run.AMP", "p['collision_transfer'].astype(LD).transpose(0,2,1)/LD(run.AMP)")
    change('h=m.t[-1]/n;pair(gamma*h);mark=time.monotonic();pair(gamma*h);warm=time.monotonic()-mark',
        "times,_,port=packets(photon/f'steps-{n}-reference-128.npz');pair(times[0]);mark=time.monotonic();pair(times[0]);warm=time.monotonic()-mark")
    change('384*warm', 'stage_count*warm')
    change('        for k in range(n):\n            a,b,mech=pair(k*h+gamma*h);',
        '        for first,last in times.reshape(-1,2):\n            a,b,mech=pair(first);')
    change('a,b,_=pair(k*h+h,mech);', 'a,b,_=pair(last,mech);')
    change("OUT/f'sweep-{sweep}/residual.json'", "OUT/'result.json'")
    change('stage_count=2*n,stage_operator_identical=True', 'stage_count=len(times),angular_port_relative=port,stage_operator_identical=True')
    change('        assert max(list(inputs.values())+source+paired)<.002,rows[-1]', '        print(json.dumps(rows[-1]),flush=True)')
    change('passed=True,rows=rows,', 'passed=maximum<.002,rows=rows,')
    change('original_waveform_change_test_passed=False,original_stop_preserved=True,',
        'actual_Radau_stage_clock_used=True,no_extra_sweep_authorized=True,')
    source+="\n    assert maximum<.002,('Actual finite reciprocal block',maximum)\n"
    return source


def prepare():
    assert read(OUT/'photon-result.json')['passed'] and read(OUT/'photon-result.json')['full_photon_horizon_completed']
    assert read(OUT/'status.json')['state']=='completed' and read(OUT/'run-receipt.json')['error'] is None
    assert read(OUT/'prepare-receipt.json')['seconds']+CAPS['prepare']<prior.CAPS['prepare']
    for p,h in read(OUT/'plan.json')['bindings'].items():
        target=OUT/'initial-producer.py' if Path(p)==Path(prior.__file__) else p
        assert sha(target)==h,p
    assert sha(prior.__file__)==read(OUT/'partition-check-repair.json')['source_sha256']
    files=[Path(__file__),Path(prior.__file__),Path(run.__file__),Path(audit.__file__),
        Path(prior.engine.__file__),Path(prior.engine.face.__file__),OUT/'plan.json',
        OUT/'photon-result.json',OUT/'run-receipt.json',OUT/'partition-check-repair.json']
    files.extend(p/f'steps-{n}-reference-128.npz' for s in [0,1] for p in paths(s) for n in [64,128] if s==0 or p==paths(1)[0])
    write(OUT/'material-block-plan.json',dict(classification='Counterexample candidate',checkpoint='abf148503',
        claim='Apply this completed175 photon collision/impulse history to the same174 analytic material operator, then compare the actual previous and returned lagged inputs and every saved Radau-stage forcing.',
        decision='Only the original0.2percent finite block and paired-E/H criterion permits a GR readout. A failure stops the registered sweep; no automatic third sweep or final-charge sign inheritance.',
        material='Same conservative/primitive map, directional minmod, analytic HLL, EOS partials1e-5 with half/double controls, SSP2, physical CFL and zero additional metric. Preserve accepted4/8macro prefixes. No photon or background recomputation.',
        stages='Use saved actual_step_edges and accepted_angular_times/weights including local bisections. Hold the first-stage mechanical slope at the closing Radau stage exactly as the evolution does; compare both affine sources and unchanged matrices.',
        partition='Compute new M from longdouble material and promoted collision transfer before subtracting or casting. No checker-first binary64 division.',
        forecast='Material uses actual SSP substeps, the larger measured current/old seconds per substep, old full-period late CFL counts, setup and2x margin. Block uses the total SAVED stage count with warmed pair/assembly costs and2x margin; no Krylov solves.',
        budgets=CAPS,total_original_seconds=prior.TOTAL,CPU_threads=1,virtual_GiB=3,
        gates=dict(time=.02,directional=.002,paired=.002,block=.002,conservation=1e-8,owner=1e-8,branch=.01,angular_port=1e-12),
        scope='Finite lagged-input and stage-equation consistency only; not a continuum contraction, uniform derivative bound or full nonlinear GR closure. Final charge remains unadjudicated.',
        final_charge_conclusion='unadjudicated',bindings={str(p):sha(p) for p in files}))
    for label,source in [('material',material_source()),('block',block_source())]:
        compile(source,__file__,'exec');(OUT/f'expanded-reciprocal-{label}.py').write_text(source)
    FunctionType(run.check.__code__,dict(run.check.__globals__,OUT=OUT,paths=paths))()


def material(pilot):
    ns=dict(run.material.__globals__,OUT=OUT,paths=paths,initialize=prior.initialize,CAPS=dict(production=CAPS['material_production']))
    exec(compile(material_source(),__file__,'exec'),ns);ns['material'](pilot)


def residual():
    FunctionType(run.residual.__code__,dict(run.residual.__globals__,OUT=OUT,paths=paths))()


def block():
    assert read(paths(1)[1]/'production.json')['passed'] and read(OUT/'result.json')['passed']
    files=[Path(__file__),OUT/'material-block-plan.json',OUT/'result.json',paths(1)[1]/'production.json']
    files.extend(p/f'steps-{n}-reference-128.npz' for s in [0,1] for p in paths(s) for n in [64,128])
    write(OUT/'block-plan.json',dict(sweep=1,scope=read(OUT/'material-block-plan.json')['scope'],
        bindings={str(p):sha(p) for p in files}))
    count=0
    for n in [64,128]:
        with np.load(paths(1)[0]/f'steps-{n}-reference-128.npz') as p:count+=len(p['accepted_angular_times'])
    owner=SimpleNamespace(initialize=lambda _:prior.initialize(),coupled=run.c,paths=paths,AMP=run.AMP,C=prior.original.C)
    ns=dict(audit.block.__globals__,OUT=OUT,run=owner,packets=run.packets,stage_count=count)
    exec(compile(block_source(),__file__,'exec'),ns)
    try:ns['block'](1)
    finally:
        if (OUT/'block-result.json').exists():
            r=read(OUT/'result.json');b=read(OUT/'block-result.json')
            r.update(passed=b['passed'],reciprocal_block_accepted=b['passed'],maximum_block_defect=b['maximum_block_defect'],
                actual_previous_free_material_input_applied=True,new_free_material_return_completed=True,
                no_extra_sweep_authorized=True,GR_return_completed=False,physical_final_charge_solved=False,
                final_charge_conclusion='unadjudicated',full_goal_complete=False)
            write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'closure-{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));prior.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=prior.TOTAL
        if action!='prepare':
            for p,h in read(OUT/'material-block-plan.json')['bindings'].items():assert sha(p)==h,p
        if action.startswith('material_'):material(action=='material_pilot')
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
