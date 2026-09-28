"""Counterexample candidate: metric return from the complete SAME joint history.

Read an accepted photon/material/GR history and its own actual angular packets.
This is the metric input for a later coupled feedback solve, not that solve.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import json,os,resource,sys,time
import numpy as np
import return_joint_gr_geometry as base
import read_completed_joint_gr as full

OUT=Path('native-complete-metric203-work');GR=full.OUT/'gr'
read,write,sha=base.read,base.write,base.sha
CAPS=dict(prepare=120,metric=600)


def saved(n):return INPUT/f'sweep-1/photons/complete-{n}.npz'


def initialize():
    fn=FunctionType(base.initialize.__code__,dict(base.initialize.__globals__,OUT=OUT));fn()


def prepare():
    result=read(INPUT/'result.json');fields=read(full.OUT/'fields.json')
    assert read(INPUT/'controller-status.json')['state']=='completed'
    assert read(INPUT/'gr-controller-status.json')['state']=='completed'
    assert result['passed'] and result['full_horizon_photon_material_completed']
    assert fields['passed'] and fields['full_declared_input_period']
    assert read(full.OUT/'plan.json')['input_directory']==str(INPUT)
    assert not OUT.exists();OUT.mkdir();(OUT/'metric').mkdir()
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    files=list((full.OUT/'sweep-0').rglob('*.npz'))+[full.OUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for p in files:
        dst=OUT/p.relative_to(full.OUT);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [full.OUT/n for n in ['sources.json','fields.json','plan.json']]+list(GR.glob('*.npz'))
    files += [saved(n) for n in [64,128]]+[INPUT/'result.json']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='2ea3d9f4a',
        claim='Construct the full-period center mass constraint and asymptotically normalized metric return using the accepted SAME photon/material/GR trajectory and its own actual angular emission history.',
        decision='Passing the original metric-time, quadrature, packet-energy and null-ray gates supplies the full-period input for an actual coupled GR-return solve. It neither accepts the old failed5.211566percent returned baryon time comparison nor establishes self-GR/final charge.',
        method='Reuse187exact-center and actual-Radau-packet mapping over all17canonical times. Same531cells,64/128clocks and4/8quadrature. No new physical integration or independent response superposition.',
        boundary='Reuse the actual emitted signed photon packets with their stored stage weights and edges. Retain leading scalar-vacuum lapse and the explicit ADM residual; the missing exterior scalar response and continuous-time metric interpolation remain unresolved.',
        gates=dict(time=.02,quadrature=.002,packet_energy=1e-12,ray=1e-10),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='The old187three-time metric used the same source/packet mapping. This extends that consumer to17times and the longer accepted packet history. Cost is not measured at full horizon; allow10minutes, not a completion-time promise. No automatic extra quadrature or time path.',
        reuse='All accepted physical states, energy and ports are immutable. Only the metric consumer executes. No prior charge, extra diagnostic mass or other solution is added.',
        stop='Any original source, field, mapping, metric/packet/ray/time gate or generous wall cap. Preserve failure before any physical returned evolution.',
        input_directory=str(INPUT),full_horizon_seconds=result['full_horizon_seconds'],
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,
        final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',base.metric.symbolic())


def metric_run():
    # Only the saved-state provider changes. The already checked exact-center
    # and actual-packet algorithms remain owned by187.
    provider=SimpleNamespace(**dict(vars(base.prior),OUT=full.OUT,saved=saved))
    boundary=FunctionType(base.Lapse.boundary.__code__,dict(base.Lapse.boundary.__globals__,prior=provider))
    class Lapse(base.Lapse):pass
    Lapse.boundary=boundary
    fn=FunctionType(base.metric_run.__code__,dict(base.metric_run.__globals__,OUT=OUT,GR=GR,Lapse=Lapse,initialize=initialize));fn()
    fields=dict(np.load(OUT/'metric/metric-128-g8.npz'))
    assert len(fields['t'])==17
    assert abs(fields['t'][-1]-read(OUT/'plan.json')['full_horizon_seconds'])<1e-18
    result=read(OUT/'metric-result.json')
    result.update(full_declared_input_period=True,same_complete_solution_energy_and_packets=True,
        physical_horizon_seconds=float(fields['t'][-1]),new_physical_steps=0,
        GR_input_applied=False,self_GR_return_closed=False,physical_final_charge_solved=False,
        final_charge_conclusion='unadjudicated',full_goal_complete=False)
    write(OUT/'metric-result.json',result)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS
    INPUT=Path(sys.argv[2]) if action=='prepare' else Path(read(OUT/'plan.json')['input_directory'])
    receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));base.prior.run.owner.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        (prepare if action=='prepare' else metric_run)()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
