"""Counterexample candidate: finish the original joint64/128 time comparison.

The accepted64 path is immutable. Only the still missing128 prefix is evolved;
the physical operator, stage solver, gates and original clocks stay unchanged.
"""
from pathlib import Path
import json,os,resource,sys,time
import numpy as np
import couple_native_fluid_radau as joint

OUT=Path('native-fluid-time180-work');BEFORE=joint.OUT
read,write,sha=joint.read,joint.write,joint.sha;LD=joint.LD
CAPS=dict(prepare=15,run=450,audit=15)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    coarse=read(BEFORE/'sweep-1/photons/pilot-64.json');admission=read(BEFORE/'fine-admission.json')
    assert coarse['passed'] and not admission['eligible']
    assert not (BEFORE/'sweep-1/photons/pilot-128.npz').exists()
    assert read(BEFORE/'branch-check-result.json')['passed']
    for sweep in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{sweep}/{folder}').mkdir(parents=True)
    sources=list((BEFORE/'sweep-0').rglob('*.npz'))
    sources.extend(BEFORE/name for name in ['normalization.json','photon-conservation-plan.json','check-result.json',
        'sweep-1/photons/pilot-64.npz','sweep-1/photons/pilot-64.json','stages-64.json'])
    reused={}
    for src in sources:
        dst=OUT/src.relative_to(BEFORE);os.link(src,dst);assert sha(src)==sha(dst)
        reused[str(dst.relative_to(OUT))]=dict(source=str(src),sha256=sha(src))
    write(OUT/'reuse.json',reused)
    files=[Path(__file__),BEFORE/'fine-admission.json',BEFORE/'branch64-receipt.json',BEFORE/'branch-check-result.json']+sources
    files.extend(Path(m.__file__) for m in list(sys.modules.values()) if getattr(m,'__file__',None)
        and Path(m.__file__).suffix=='.py' and Path(m.__file__).parent==Path(joint.__file__).parent)
    forecast=2*coarse['stepping_seconds']+coarse['operator_point_seconds']+20
    assert forecast<CAPS['run']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='b00cbabbe',
        previous_turn='Progress: unchanged true joint equations completed the original64prefix after the branch-Jacobian repair. The original480s action budget correctly rejected the still missing128path; no time verdict exists.',
        claim='Determine whether the actual same-solution photon and four-material-variable response survives the original64/128 time comparison before any full-horizon or GR charge claim.',
        decision='All original photon six-channel and material four-channel time differences must be below2percent. A failure stops full-horizon admission and identifies the unresolved component; a pass permits a separately costed continuation of these exact saved states.',
        reuse='Accepted64path and all EOS/background/source banks are reused byte-for-byte. Only the original128 eight-macro-step T/16 path is computed. Existing stage equations, initial-guess history, Radau clocks and front splits are unchanged.',
        reassessment='The former421s estimate exceeded the former210s remainder. That rejection is immutable. Reassess the one still necessary path under a new explicit450s cap, without replaying64 or adding a clock. Saved arrays alone cannot establish its missing time error. A new preconditioner/GPU implementation would require new equation-equivalence evidence before it could replace this bounded calculation.',
        forecast=dict(nominal_seconds=forecast,assumed_range_seconds=[300,450],hard_cap_seconds=CAPS['run'],
            basis='Twice the measured64stepping cost plus measured operator construction and20s reserve. Fine-step Krylov and branch costs are unmeasured; the range is an assumption, not a guarantee.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        gates=dict(time=.02,native_direction=.002,stage=1e-12,physical_stage_moment=1e-13,conservation=1e-8,port=1e-12,max_Newton_solves=3),
        stop='Any original gate, stage iteration or resource cap stops this path. No finer clock, additional split, alternate solver, threshold relaxation or full-horizon dispatch is authorized by this plan.',
        scope='Retained directional response with zero additional metric. Full nonlinear/EOS/spatial/boundary/static/observational closure and final charge remain unproven.',
        final_charge_conclusion='unadjudicated',full_horizon_authorized=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def run():
    joint.OUT=OUT
    joint.pilot((128,))
    # compare_paths saves its failed verdict before raising; preserve it.
    joint.compare_paths()


def audit():
    for rel,row in read(OUT/'reuse.json').items():assert sha(OUT/rel)==sha(row['source'])==row['sha256'],rel
    result=read(OUT/'pilot-result.json');rows=[]
    for n in [64,128]:
        path=OUT/f'sweep-1/photons/pilot-{n}.npz'
        times,energy,port=joint.previous.run.packets(path)
        with np.load(path) as p:
            assert np.array_equal(times,p['joint_stage_times'])
            assert np.array_equal(p['joint_stage_weights'],p['accepted_angular_quadrature_weights'])
            q=p['joint_stage_conserved_scaled'];assert len(q)==len(times)
            native=p['joint_native_rates_scaled'];collision=p['joint_collision_rates_scaled']
            integral=np.sum(p['joint_stage_weights'][:,None,None].astype(LD)*(native+collision),axis=0,dtype=LD)
            # Recover physical values from the archived conserved state, independent
            # of an assumed EOS cx. The endpoint runner gas is already scaled.
            final=p['conserved_material_history'][-1]/joint.AMP
            thermal=p['delta_material'][:,0]/joint.AMP*p['material_energy_units']
            expected=np.column_stack([thermal,final[3],final[0],final[1]])
            error=(np.sum(abs(expected+p['material_floor_discard_scaled']-integral),axis=0,dtype=LD)
                /np.maximum(np.sum(abs(expected),axis=0,dtype=LD),LD('1e-290'))).astype(float).tolist()
            assert max(error)<1e-8,error
            rows.append(dict(clock=n,accepted_stages=len(times),same_solution_material_ledger=error,angular_port_relative=port))
    write(OUT/'audit-result.json',dict(classification='Counterexample candidate',passed=True,rows=rows,
        original_coarse_bytes_unchanged=True,time_comparison_passed=result['passed'],
        final_charge_conclusion='unadjudicated',full_goal_complete=False))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    started=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-started,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
