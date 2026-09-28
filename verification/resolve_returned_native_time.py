"""Counterexample candidate: resolve188's actual returned baryon time error.

First separate the quadrature/base-branch contribution from the evolving-state
contribution, using only the saved Radau trajectories and live native owner.
"""
from pathlib import Path
import gc,json,os,resource,sys,time
import numpy as np
import evolve_same_solution_gr_return as prior

OUT=Path('native-return-time189-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
LD,AMP=prior.LD,prior.AMP
CAPS=dict(prepare=20,diagnose=75)


def prepare():
    assert not OUT.exists();OUT.mkdir();assert not read(OLD/'result.json')['passed']
    files=[];reused={}
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(OLD);os.link(src,dst);reused[str(dst)]=sha(src);files.append(src)
    files += [OLD/f'sweep-1/photons/return-{n}.npz' for n in [64,128]]
    files += [prior.prior.prior.saved(n) for n in [64,128]]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='b4bdec995',
        claim='Identify whether188baryon time failure is dominated by integration of the current native operator/branch history or by the evolved lower state, before modifying the actual coupled evolution.',
        method='Evaluate the coarse Radau dense polynomial at the actual fine stages against the same saved fine base. Integrate with the saved actual weights. Algebraically reconstruct the original endpoint difference including floor. Compare base branches at shared stages, and inspect the local native B/S block only if needed.',
        decision='A dominant operator-history defect requires its temporal representation/branch repair; a dominant evolved-state defect requires a measured stiffness/cadence decision. Do not blindly extend the period or add time grids.',
        gates=dict(saved_native=1e-12,dense_stage=1e-12,decomposition=1e-8),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,new_physical_steps=0,
        forecast='188constructor15s; current native RHS applications were small compared with sparse Jacobians. One constructor plus8fine-stage probes and2shared-stage branch controls expected20..55s, capped75s. No evolution or refinement admitted by this plan.',
        stop='Any anchor, dense-polynomial, algebra or cost failure. Preserve188failure and original2percent gate. No modification to185or188producers.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused))


def gas(m,q):
    return np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])


def dense(m,p,t):
    edges=p['actual_step_edges'];k=int(np.clip(np.searchsorted(edges,t,side='left')-1,0,len(edges)-2))
    h=edges[k+1]-edges[k];theta=LD(t-edges[k])/LD(h)
    initial=np.zeros((m.n,4),LD) if k==0 else gas(m,p['joint_stage_conserved_scaled'][2*k-1])
    initial=initial.copy();initial[~m.material.active(edges[k])]=0
    rates=(p['joint_native_rates_scaled']+p['joint_collision_rates_scaled'])[2*k:2*k+2]/m.units
    b=np.array([LD('1.5')*theta-LD('.75')*theta**2,LD('.75')*theta**2-LD('.5')*theta])
    return initial+LD(h)*np.einsum('j,jnk->nk',b,rates)


def diagnose():
    prior.OUT=OUT;prior.initialize();m=prior.Model(128)
    coarse,fine=[dict(np.load(OLD/f'sweep-1/photons/return-{n}.npz')) for n in [64,128]]
    original_anchor=m.anchor;controls=[];probed=[]
    for p in [coarse,fine]:
        norms=[]
        for j,t in enumerate(p['joint_stage_times']):
            stored=gas(m,p['joint_stage_conserved_scaled'][j]);pred=dense(m,p,float(t))
            error=np.sum(abs(pred-stored)*m.units,axis=0)/np.maximum(np.sum(abs(stored)*m.units,axis=0),LD('1e-290'))
            norms.append(error.astype(float).tolist())
        controls.append(norms)
    assert np.max(controls[0])<1e-12 and np.max(controls[1])<1e-12,('Dense stored stage reconstruction',controls)
    for j,t in enumerate(fine['joint_stage_times']):
        now=float(t);g=dense(m,coarse,now);value=m.native(now,g,details=True)[0]*m.units
        probed.append(value)
    saved=[np.sum(p['joint_stage_weights'][:,None,None].astype(LD)*p['joint_native_rates_scaled'],axis=0,dtype=LD) for p in [coarse,fine]]
    projected=np.sum(fine['joint_stage_weights'][:,None,None].astype(LD)*np.array(probed),axis=0,dtype=LD)
    floor=coarse['material_floor_discard_scaled']-fine['material_floor_discard_scaled']
    delta=coarse['conserved_material_history'][-1,0]/AMP-fine['conserved_material_history'][-1,0]/AMP
    quadrature=saved[0][:,2]-projected[:,2];path=projected[:,2]-saved[1][:,2]
    norm=np.sum(abs(delta),dtype=LD);defect=delta-quadrature-path+floor[:,2]
    common=[];coarse_anchor=dict(np.load(prior.prior.prior.saved(64)))
    for j,t in enumerate(coarse['joint_stage_times']):
        k=int(np.argmin(abs(fine['joint_stage_times']-t)))
        if abs(fine['joint_stage_times'][k]-t)>1e-18:continue
        state=gas(m,fine['joint_stage_conserved_scaled'][k]);m.anchor=coarse_anchor;m.checked_times=set()
        other=m.native(float(t),state)*m.units;m.anchor=original_anchor;m.checked_times=set()
        actual=fine['joint_native_rates_scaled'][k]
        difference=np.sum(abs(other-actual),axis=0)/np.maximum(np.sum(abs(actual),axis=0),LD('1e-290'))
        common.append(dict(time=float(t),same_lower_state_base_branch_relative=difference.astype(float).tolist()))
    result=dict(classification='Counterexample candidate',passed=float(np.sum(abs(defect))/norm)<1e-8,
        scope='Exact saved-history algebra with Radau dense coarse state, evaluated on the fine actual high-state stages. Not an independent evolved correction.',
        dense_stage_relative_max=float(max(np.max(v) for v in controls)),
        baryon_quadrature_and_base_history_L1_over_difference=float(np.sum(abs(quadrature))/norm),
        baryon_evolved_state_L1_over_difference=float(np.sum(abs(path))/norm),
        baryon_floor_L1_over_difference=float(np.sum(abs(floor[:,2]))/norm),
        decomposition_relative=float(np.sum(abs(defect))/norm),common_stage_controls=common,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'decomposition.npz',baryon_difference=delta,quadrature_and_base_history=quadrature,
        evolved_state=path,floor=floor,coarse_dense_on_fine_native=np.array(probed),fine_times=fine['joint_stage_times'],fine_weights=fine['joint_stage_weights'])
    write(OUT/'diagnosis.json',result);print(json.dumps(result),flush=True);assert result['passed']


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            plan=read(OUT/'plan.json')
            for p,h in dict(plan['bindings'],**plan['reused']).items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
