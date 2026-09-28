"""Return the accepted corrected photon collision history to free matter.

Counterexample candidate. This is the next physical block of the same
solution; a material pass alone does not accept reciprocity or final charge.
"""
from pathlib import Path
import gc,json,resource,shutil,sys,time
import numpy as np
import continue_native_pressure_response as previous

front=previous.prior;original=front.original;c=original.c
OUT=Path('native-pressure-matter174-work');OLD=previous.OUT
read,write,sha=original.read,original.write,original.sha
LD=np.longdouble;AMP=c.AMP
CAPS=dict(prepare=30,check=30,pilot=90,production=450,residual=40)


def paths(sweep):return OUT/f'sweep-{sweep}/photons',OUT/f'sweep-{sweep}/material'


def packets(path):
    """Use accepted Radau times and weights, including every local substep."""
    with np.load(path) as p:
        edges=p['actual_step_edges'];dt=np.diff(edges)
        times=p['accepted_angular_times'];weights=p['accepted_angular_quadrature_weights']
        assert np.all(dt>0) and times.shape==weights.shape==(2*len(dt),)
        assert np.max(abs(times-(edges[:-1,None]+dt[:,None]*front.radau.RK_C).ravel()))<1e-18
        assert np.max(abs(weights-(dt[:,None]*front.radau.RK_B).ravel()))<1e-18
        energy=weights.astype(LD)[:,None]*p['accepted_angular_luminosity'].astype(LD)
        assert energy.shape==(len(times),4)
        angular=np.arange(1,8,2,dtype=LD)/32
        cumulative=np.r_[LD(0),np.cumsum(energy@angular,dtype=LD)]
        ids=np.array([int(np.argmin(abs(edges-t))) for t in p['t']])
        assert np.max(abs(edges[ids]-p['t']))<1e-18
        error=float(np.max(abs(cumulative[2*ids]-p['radial_ports'][:,1,1]))/max(np.sum(abs(energy)@angular),LD('1e-290')))
        assert error<1e-12,error
        return times.copy(),energy,error


def initialize():
    original.OUT=OUT;original.METRIC=OUT/'metric';original.install();original.initialize(1)
    # The exact accepted runner already handles nonzero material feedback.
    # Keep its equation/stages unchanged and only bind the destination owner.
    source=(OLD/'sweep-1/expanded-front-run.py').read_text()
    namespace=dict(c.Response.run.__globals__,stages=front.radau.stages,
        split_flags=front.split_flags,RK_A=front.radau.RK_A,RK_B=front.radau.RK_B,RK_C=front.radau.RK_C)
    exec(compile(source,__file__,'exec'),namespace);c.Response.run=namespace['run']
    (OUT/'sweep-1/expanded-corrected-run.py').write_text(source)


def prepare():
    assert not OUT.exists();result=read(OLD/'result.json')
    assert result['passed'] and result['full_horizon_photon_thermal_H_completed']
    assert read(OLD/'status.json')['state']=='completed' and read(OLD/'run-receipt.json')['error'] is None
    OUT.mkdir();files=[Path(__file__),Path(previous.__file__),Path(front.__file__),Path(original.__file__),
        Path(original.owner.__file__),Path(c.__file__),OLD/'plan.json',OLD/'result.json',OLD/'run-receipt.json',OLD/'sweep-1/expanded-front-run.py']
    for s in [0,1]:
        for p in paths(s):p.mkdir(parents=True)
    for p in (OLD/'sweep-0').rglob('*.npz'):
        q=OUT/p.relative_to(OLD);shutil.copyfile(p,q);assert sha(q)==sha(p);files.append(p)
    for n in [64,128]:
        for suffix in ['.npz','.json']:
            p=OLD/f'sweep-1/photons/steps-{n}-reference-128{suffix}';q=OUT/p.relative_to(OLD)
            shutil.copyfile(p,q);assert sha(q)==sha(p);files.append(p)
    for name in ['normalization.json','photon-conservation-plan.json','sweep-1/photons/result.json']:
        p=OLD/name;shutil.copyfile(p,OUT/name);files.append(p)
    probe=Path('native-mixed-return162-work/material-arithmetic-plan.json');files.append(probe)
    choice=read(probe);assert choice['selected_nominal_arithmetic_probe']==16.
    write(OUT/'material-arithmetic-plan.json',dict(classification='Counterexample candidate',
        selected_factor=choice['selected_factor'],selected_nominal_arithmetic_probe=16.,
        scope='Reuse the existing arithmetic probe on the same constitutive owner; its earlier pass is not evidence for this new input. Enforce the original half/nominal/double and trajectory checks again.',
        bindings={str(probe):sha(probe)}))
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0534c8798',
        claim='Apply the accepted corrected photon collision and impulse histories to the SAME free-material equation, preserving energy/species and actual photon boundary packets.',
        decision='Only a full original-gate material pass permits a reciprocal source audit and subsequent GR return. Final charge remains unadjudicated; no diagnostic work is added to mass.',
        reuse='Copy both completed173photon paths byte-for-byte. No photon replay, new EOS calls for background, new grid, pulse, period or radiation method.',
        time='Same17canonical cumulative collision/impulse inputs, existing SSP2 fluid owner,64/128 macro clocks and physical CFL. Paired4/8macro-step pilots have the same horizon, and are reused in production.',
        gates=dict(time=.02,directional=.002,conservation=1e-8,owner=1e-8,branch=.01,angular_port=1e-12),
        budget=dict(actions=CAPS,total_action_seconds=sum(CAPS.values()),CPU_threads=1,virtual_GiB=3,material_path_pairs=1),
        forecast='Old corrected-probe fine material path took116.415s for9992raw calls. Measure each current pilot and use the larger current/old seconds per raw call, old full-path counts and2x margin; continuation must fit450s. Changed later CFL and arithmetic costs are assumptions.',
        stop='Any original physical/time gate or cap stops the action. Preserve failure and accepted prefix. No automatic probe ladder, extra sweep, grid increase or GR acceptance.',
        final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def check():
    # Exact two-stage Radau moments on unequal intervals (rational arithmetic).
    from fractions import Fraction as F
    edges=[F(0),F(1,4),F(1)];times=[];weights=[]
    for a,b in zip(edges,edges[1:]):
        times.extend([a+(b-a)/3,b]);weights.extend([(b-a)*F(3,4),(b-a)/4])
    for power in range(3):assert sum(w*t**power for w,t in zip(weights,times))==F(1,power+1)
    rows=[]
    for n in [64,128]:
        tt,e,err=packets(paths(1)[0]/f'steps-{n}-reference-128.npz')
        rows.append(dict(base_clock=n,actual_stage_count=len(tt),angular_port_relative=err))
    write(OUT/'packet-check.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        scope='Actual same-solution boundary integral; no GR propagation or final charge yet.'))
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,
        scope='Rational Radau quadrature integrates degrees0..2 exactly on unequal intervals; no full coupled error bound.'))


def material(pilot):
    assert read(OUT/'packet-check.json')['passed'];initialize();folder=paths(1)[1];rows=[];histories=[];forecasts=[];started=time.monotonic()
    if not pilot:assert read(folder/'pilot.json')['eligible']
    for n in [64,128]:
        label=f'pilot-{n}' if pilot else f'steps-{n}-reference-128'
        mark=time.monotonic();m=c.Material(128,n)
        row=m.run(n,label,n//16 if pilot else None,None if pilot else f'pilot-{n}')
        wall=time.monotonic()-mark
        with np.load(folder/f'{label}.npz') as d:
            z=d['delta_scaled'];end=d['t'][-1]
            rates=[m.rhs(end,z,p)[0] for p in [.5,1.,2.]]
            probe=[(np.sum(abs(r-rates[1]),axis=1)/np.maximum(np.sum(abs(rates[1]),axis=1),1.)).astype(float).tolist() for r in [rates[0],rates[2]]]
            with np.load(paths(1)[0]/f'steps-{n}-reference-128.npz') as p:
                clock=p['t'][p['t']<=end+1e-18]
                ids=[int(np.argmin(abs(d['t']-t))) for t in clock];assert np.max(abs(d['t'][ids]-clock))<1e-18
                histories.append(d['history_scaled'][ids].copy())
            if not pilot:
                with np.load(folder/f'pilot-{n}.npz') as old:
                    for key in ['t','history_scaled','ledgers_scaled','norms_scaled','discards_scaled']:
                        assert np.array_equal(d[key][:len(old[key])],old[key]),('Accepted prefix changed',key)
        row.update(worker_wall_seconds=wall,physical_branch_ratio=m.physical_branch_ratio,
            maximum_owner_error=max(v['owner_error'] for v in m.cache.values()),endpoint_probe_half_nominal_double=probe,nominal_arithmetic_probe=16.)
        row['passed']=bool(row['passed'] and row['physical_branch_ratio']<.01 and row['maximum_owner_error']<1e-8 and np.max(probe)<.002)
        write(folder/f'{label}.json',row);rows.append(row);assert row['passed'],row
        if pilot:
            old=read(Path('native-mixed-return162-work/sweep-1/material')/f'steps-{n}-reference-128.json')
            per_call=max(row['seconds']/row['raw_owner_calls'],old['seconds']/old['raw_owner_calls'])
            forecasts.append(max(0,old['raw_owner_calls']-row['raw_owner_calls'])*per_call+max(0,wall-row['seconds'])+5)
        del m;gc.collect()
    errors=c.relative(*histories);result=dict(classification='Counterexample candidate',passed=max(errors)<.02,
        rows=rows,time_comparison=errors,seconds=time.monotonic()-started,
        photon_paths_recomputed=False,final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False)
    if pilot:result.update(upper_remaining_seconds=2*sum(forecasts),eligible=result['passed'] and 2*sum(forecasts)<CAPS['production'])
    write(folder/('pilot.json' if pilot else 'production.json'),result);print(json.dumps(result),flush=True)
    assert result['eligible' if pilot else 'passed'],result


def residual():
    assert read(paths(1)[1]/'production.json')['passed'];rows=[]
    for n in [64,128]:
        with np.load(paths(1)[1]/f'steps-{n}-reference-128.npz') as d,np.load(paths(1)[0]/f'steps-{n}-reference-128.npz') as p:
            ids=[int(np.argmin(abs(d['t']-t))) for t in p['t']];assert np.max(abs(d['t'][ids]-p['t']))<1e-18
            mismatch=c.relative(p['moments'][:,[1,2]]/AMP,d['history_scaled'][ids][:,[2,3]])
            balance=float(np.max(abs(np.sum(d['history_scaled'],axis=2,dtype=LD)+d['discards_scaled']-d['ledgers_scaled'])/np.maximum(d['norms_scaled'],1.)))
        _,_,angular=packets(paths(1)[0]/f'steps-{n}-reference-128.npz');assert balance<1e-8
        rows.append(dict(base_clock=n,photon_material_energy_H_residual=mismatch,material_balance=balance,angular_port=angular))
    write(OUT/'result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        actual_same_photon_history_applied_to_free_material=True,reciprocal_block_accepted=False,
        GR_return_completed=False,final_charge_conclusion='unadjudicated',physical_final_charge_solved=False,full_goal_complete=False))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=sum(CAPS.values())
        if action in ['pilot','production']:material(action=='pilot')
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
