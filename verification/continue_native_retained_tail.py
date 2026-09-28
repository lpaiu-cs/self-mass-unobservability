"""Finish the original fine path, preserving the interrupted phase verbatim."""
from pathlib import Path
from types import FunctionType
import hashlib,io,json,shutil,signal,sys,time

ROOT=Path('outputs/direct-eos-gr33')
OLD=ROOT/'native-retained-tail/supported-temperature';BEFORE=OLD/'evolution'
OUT=ROOT/'native-retained-completion';EV=OUT/'evolution';GR=OUT/'gr'
read=lambda p:json.loads(p.read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def prepare():
    assert not OUT.exists();EV.mkdir(parents=True);GR.mkdir()
    names=['checkpoint-128.npz','history-128.npz','accepted-ports-128.npz','tail-128.json']
    names += [f'state-128-{k}.npz' for k in [0,8,16,24,32]]
    for name in names:shutil.copyfile(BEFORE/name,EV/name)
    # Latest coarse block includes actual late substepping and storage cost.
    late=(601.943791785-380.2656909950001)/16/2
    fixed=read(BEFORE/'cost-64-64.json')['cost']['setup']+43.
    estimate=112*late+fixed;upper=1.5*estimate+20.
    paths=[Path(__file__),Path('verification/evolve_native_retained_tail.py'),
        Path('verification/def_native_updated_gr_return.py'),Path('verification/def_native_refined_thermochemistry.py'),
        Path('verification/def_native_stage_energy_history.py'),OLD/'support/repaired-bank.npz',OLD/'support/spectrum-bank.npz',
        OLD/'checkpoint-control.json',OLD/'production-external.log',BEFORE/'result-64.json',BEFORE/'source-64.npz']
    paths += [BEFORE/name for name in names]
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='6014b35e4',previous_goal_turn='progress',
        claim='Complete the unchanged native retained-material128-step trajectory and apply both actual source histories to the same compact characteristic GR charge operator.',
        decision='Does the original two-clock source comparison pass and does the retained-material change matter at the charge accuracy scale?',
        reuse='Preserve all Phase147 files. Reuse completed64-step results and native banks. Start from the complete16-step fine checkpoint. Recomputing17..32 is necessary because the interrupted32-step checkpoint lacks eight continuation ledgers; compare the preserved32-step physical record before proceeding.',
        change='Only buffered atomic storage and correctly aligned restart histories from Phase147. No equation, input amplitude, grid, native support, horizon or acceptance relaxation.',
        clocks=[64,128],coarse_rerun=False,new_native_states=0,CPU_threads=1,memory_GiB=3,
        budget=dict(fine_seconds=1500,required_prefix_seconds=240,GR_seconds=150,total_scientific_seconds=1650),
        forecast=dict(late_coarse_equivalent_fine_step_seconds=late,setup_and_export_seconds=fixed,
            point_seconds=estimate,upper_seconds=upper,assumption='Same physical-time late substep cost at the fine clock; finer-clock overhead and later support changes remain unmeasured. Use max of this rate and the required16..32 replay rate before the remainder. No completion guarantee.'),
        gates=dict(restart_response=1e-12,energy=1e-8,port=1e-10,baryon=1e-10,time=.02,quadrature=.002,independent=1e-9),
        stop='One fine continuation only. Stop on prefix disagreement, support, conservation, measured remainder forecast or external deadline. No automatic second run, further floor reduction, finer clock or horizon. Preserve any failure.',
        limits='The compact same-initial-operator first variation does not supply updated native collision corrections, null-infinity emission, uniform EOS/derivative, full floor/coupled/nonlinear/static or observation closure.',
        bindings={(p.relative_to(Path.cwd()) if p.is_absolute() else p).as_posix():sha(p) for p in paths}))
    print(json.dumps(dict(prepared=True,forecast_seconds=estimate,upper_seconds=upper,cap=1500)),flush=True)


def numpy_io():
    import numpy as np
    def load(p):return np.load(io.BytesIO(p.read_bytes()))
    return np,load


def evolve():
    start=time.monotonic();plan=read(OUT/'plan.json');cap=plan['budget']['fine_seconds']
    for p,h in plan['bindings'].items():assert sha(p)==h,p
    import evolve_native_retained_tail as run
    np,load=numpy_io();run.support.owner.deadline(start,plan['budget']['required_prefix_seconds'])
    assert int(load(EV/'checkpoint-128.npz')['completed'])==16
    run.EV=EV;run.GR=GR
    run.capture_run=FunctionType(run.capture_run.__code__,dict(run.capture_run.__globals__,OUT=EV))
    try:
        model=run.Capture();prefix_start=time.monotonic();prefix=model.run(128,32,resume=True)
        assert prefix['passed'];checks={}
        for k in [24,32]:
            current=load(EV/f'state-128-{k}.npz');previous=load(BEFORE/f'state-128-{k}.npz');base=load(BEFORE/'state-128-0.npz')
            for key in current.files:
                a,b,z=current[key],previous[key],base[key]
                scale=max(float(np.sum(abs(b-z),dtype=np.longdouble)),1e-100)
                checks[f'state{k}/{key}']=float(np.sum(abs(a-b),dtype=np.longdouble))/scale
        current=load(EV/'history-128.npz');previous=load(BEFORE/'history-128.npz')
        for key in ['moments','t','ports','discard','escape']:
            a,b=current[key],previous[key];assert a.shape==b.shape
            if key=='moments':
                for j in range(a.shape[1]):
                    scale=max(float(np.max(np.sum(abs(b[:,j]),axis=-1))),1e-100)
                    checks[f'moments/{j}']=float(np.max(np.sum(abs(a[:,j]-b[:,j]),axis=-1)))/scale
            else:checks['history/'+key]=float(np.max(abs(a-b)))/max(float(np.max(abs(b))),1e-100)
        a=load(EV/'accepted-ports-128.npz')['angular_luminosity'];b=load(BEFORE/'accepted-ports-128.npz')['angular_luminosity']
        checks['angular_port']=float(np.max(abs(a-b)))/float(np.max(abs(b)))
        prefix_wall=time.monotonic()-prefix_start
        rate=max(prefix_wall/16,plan['forecast']['late_coarse_equivalent_fine_step_seconds'])
        forecast=96*rate+43.;upper=1.5*forecast+20.;remaining=cap-(time.monotonic()-start)
        replay=max(checks.values())<plan['gates']['restart_response']
        eligible=prefix['passed'] and replay and upper<remaining
        record=dict(classification='Counterexample candidate',passed=replay,prefix=prefix,checks=checks,
            prefix_wall_seconds=prefix_wall,elapsed_seconds=time.monotonic()-start,
            remaining_forecast_seconds=forecast,remaining_upper_seconds=upper,remaining_cap_seconds=remaining,eligible=eligible)
        write(OUT/'prefix.json',record);print(json.dumps(record),flush=True)
        assert eligible,('Restart or remainder budget gate failed',record)
        run.support.owner.deadline(start,cap)
        fine=model.run(128,resume=True);assert fine['passed'] and fine['complete']
        run.support.owner.deadline(start,cap)
        result=dict(classification='Counterexample candidate',passed=True,
            paths=[read(BEFORE/'result-64.json'),fine],fine_seconds=time.monotonic()-start,
            coarse_reused=True,prefix_replay_passed=True,both_clocks_completed=True,
            source_paths={64:str(BEFORE/'source-64.npz'),128:str(EV/'source-128.npz')},
            new_GR_charge_available=False,final_charge_solved=False,full_goal_complete=False)
        write(EV/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:
        write(OUT/'evolution-failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


def charge():
    start=time.monotonic();plan=read(OUT/'plan.json');assert read(EV/'result.json')['passed']
    import evolve_native_retained_tail as run
    import verify_native_stage_energy_charge as original
    import sympy as sp
    np,load=numpy_io();cap=plan['budget']['GR_seconds'];run.support.owner.deadline(start,cap)
    paths={64:BEFORE,128:EV};bindings=[Path(__file__),Path(original.__file__),EV/'result.json']
    bindings += [paths[n]/f'source-{n}.npz' for n in [64,128]]
    write(GR/'plan.json',dict(classification='Counterexample candidate',hard_seconds=cap,
        measured_basis='Prior same-operator full three-readout calculation14.386s; latest coarse constructor32.20s plus mounted I/O overhead.150s includes imports, one constructor, three paths, independent direct readout and saved histories.',
        gates=plan['gates'],new_fluid_steps=0,new_native_states=0,
        bindings={str(p):sha(p) for p in bindings}))
    try:
        model=original.gr.Response();waves={};timings=[]
        for n,order in [(128,8),(64,8),(128,4)]:
            mark=time.monotonic();d=dict(load(paths[n]/f'source-{n}.npz'))
            waves[n,order]=original.read(model,d,order)
            run.atomic_npz(GR/f'wave-{n}-g{order}.npz',**waves[n,order])
            timings.append(dict(steps=n,order=order,seconds=time.monotonic()-mark))
            run.support.owner.deadline(start,cap)
        fine=waves[128,8];norm=max(float(np.max(abs(fine['free_scalar']))),1e-300)
        temporal=float(np.max(abs(fine['free_scalar'][::2]-waves[64,8]['free_scalar'])))/norm
        quadrature=float(np.max(abs(fine['free_scalar']-waves[128,4]['free_scalar'])))/norm
        direct,residual=original.independent.direct(model,dict(load(EV/'source-128.npz')),8)
        agreement=float(abs(direct/fine['direct_scalar'][-1]-1));comparisons=[];changes={};baryon=[]
        for n in [64,128]:
            old=load(run.physical.GR/f'wave-{n}-g8.npz');now=waves[n,8]
            delta=now['free_scalar']-old['free_scalar'];changes[n]=delta
            z=load(paths[n]/f'checkpoint-{n}.npz');assert int(z['completed'])==n
            history=load(paths[n]/f'history-{n}.npz');before=load(run.physical.EV/f'history-{n}.npz')
            m=model.model
            def loss(v):return float((v['discard'][-1,2]+m.m.a0*m.cx*v['discard'][-1,0])*m.gas_scale)
            tail=read(paths[n]/f'tail-{n}.json')['rows'];assert len(tail)==n+1
            comparisons.append(dict(steps=n,endpoint=float(now['free_scalar'][-1]),prior_endpoint=float(old['free_scalar'][-1]),
                change=float(delta[-1]),relative_change=float(delta[-1]/old['free_scalar'][-1]),
                discarded_Killing_energy_erg=loss(history),prior_discarded_Killing_energy_erg=loss(before),
                remaining_discard_ratio=loss(history)/loss(before),retained_tail_endpoint_g=tail[-1]['mass_g']))
            source=load(paths[n]/f'source-{n}.npz')
            defect=abs(np.sum(source['baryon_g'][-1],dtype=np.longdouble)+z['discard'][0]*m.gas_scale/run.C**2)
            initial=np.sum(m.flow.initial[0]*m.m.vol,dtype=np.longdouble)*m.gas_scale/run.C**2
            baryon.append(float(defect/initial));assert len(load(paths[n]/f'accepted-ports-{n}.npz')['angular_luminosity'])==n
        delta_norm=max(float(np.max(abs(changes[128]))),1e-300)
        delta_time=float(np.max(abs(changes[128][::2]-changes[64])))/delta_norm
        old_quad=load(run.physical.GR/'wave-128-g4.npz')['free_scalar']
        delta_quad=float(np.max(abs(waves[128,4]['free_scalar']-old_quad-changes[128])))/delta_norm
        A,x,y=sp.symbols('A x y');assert sp.expand(A*(x+y)-A*x-A*y)==0
        passed=temporal<.02 and quadrature<.002 and agreement<1e-9 and max(baryon)<1e-10
        run.support.owner.deadline(start,cap)
        result=dict(classification='Counterexample candidate',passed=bool(passed),comparisons=comparisons,
            time_relative=temporal,quadrature_relative=quadrature,independent_direct_relative=agreement,
            inverse_residual_cm=float(residual),baryon_relative_to_initial=baryon,
            tail_change_time_relative=delta_time,tail_change_quadrature_relative=delta_quad,
            tail_change_resolved_on_original_gates=delta_time<.02 and delta_quad<.002,
            current_actual_sources_applied_to_compact_GR=True,symbolic_linearity_check=True,
            timings=timings,seconds=time.monotonic()-start,new_emission_at_null_infinity=False,
            uniform_EOS_derivative_bound=False,full_floor_feedback_enclosed=False,coupled_fixed_point_verified=False,
            nonlinear_GR=False,final_charge_solved=False,full_goal_complete=False)
        write(GR/'result.json',result);print(json.dumps(result),flush=True)
    except Exception as exc:write(GR/'failure.json',dict(error=repr(exc),seconds=time.monotonic()-start));raise
    finally:signal.setitimer(signal.ITIMER_REAL,0.)


if __name__=='__main__':
    if sys.platform!='win32':
        import resource
        resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,)*2)
    globals()[sys.argv[1]]()
