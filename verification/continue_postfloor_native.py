"""Counterexample candidate: restart the actual solver from its accepted state.

The finite floor map changes the accepted gas. Its pre-map value must not be
the next stage's default branch proposal. Neither the native RHS nor the floor
map changes; every proposed stage still passes the original acceptance gates.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import continue_precise_native as prior

OUT=Path('native-postfloor-continuation200-work');OLD=prior.OUT
before,stable,owner,joint=prior.before,prior.stable,prior.owner,prior.joint
LD=joint.LD;read,write,sha=prior.read,prior.write,prior.sha
KEYS=prior.prior.prior.KEYS
CAPS=dict(prepare=120,check=180,coarse=7200,fine=10800,audit=120)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    failure=read(OLD/'failure-64.json')
    assert read(OLD/'controller-status.json')['state']=='failed' and 'TimeoutError' in failure['error']
    assert failure['actual_accepted_steps']==114
    q=dict(np.load(OLD/'last-accepted-64.npz'))
    cells=np.flatnonzero(np.any(q['g']!=q['restart_guide'],axis=1))
    assert np.array_equal(cells,[260,261]) and not np.any(q['g'][cells])
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    files += [OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    reused={}
    for p in files:dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['last-accepted-64.npz','failed-linear-64.npz','failure-64.json','coarse-receipt.json','linear-64.json','prefix-result.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='bdfb6a4c7',
        claim='Continue the remaining actual coupled period from114accepted substeps with Newton proposals initialized from the accepted post-floor gas.',
        evidence='199accepted the old failed step after native precision repair, then reached its1800s cap inside a later linear correction. Its114state and current linear iterate are saved. The accepted gas is zero at cells260/261, while the stored pre-floor guide at261has B=-5.504e15 and S=-2.269e12. Such a guide does not represent the accepted post-map initial state. This is a proposal-quality hypothesis, not a proven cause of the remaining residual.',
        change='Before each actual stage, use g.copy() as the initial guide, including both time paths. Reuse198native precision and197linear polishing. Do not change equations, floor, fronts, physical source, period, clocks, EOS/photon maps or any acceptance gate.',
        reuse='Restore every114accepted physical state/history/ledger key exactly before changing the proposal. Preserve the interrupted199linear proposal as rejected/unadjudicated evidence; its matrix/RHS must not be reused for a new guide as if identical. No accepted prefix replay.',
        budget_reason='199used1805.43s for one accepted step and eight proposals of the next, with the eighth interrupted.30/45minute whole-path limits were not sufficient for this measured stiff region. Allow2/3hours; costs of later/fine steps remain unmeasured and these are caps, not completion forecasts.',
        iteration_budget='Up to8fresh Newton proposals per stage with the changed initial guide,12linear refinements and restart80/maxiter10. Preserve the earlier8attempted proposals separately. Do not silently accept a capped or failed proposal.',
        persistence='Write the current complete accepted state before every actual stage, so a later interruption cannot discard an accepted substep. This adds compression/IO, not physical integration.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        decision='Exact restart first; then actual coarse, fine and original10channel2percent comparison. Only a passing full pair admits the same-solution complete GR reader.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        stop='Original physical/accuracy failure,8Newton/12linear limits,nonfinite state or2/3hour wall cap. No new grid,period,parameter path or automatic rerun.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'guide-input.json',dict(classification='Counterexample candidate',cells=cells.tolist(),
        accepted=q['g'][cells].astype(float).tolist(),pre_floor_guide=q['restart_guide'][cells].astype(float).tolist(),
        physical_steps_advanced=0,convergence_improvement_demonstrated=False))


def restore_substep(m,previous,moment):
    values=restore_previous(m,previous,moment);moment=values[-1]
    q=dict(np.load(OLD/'last-accepted-64.npz'))
    assert int(q['macro_index'])==61 and int(q['sub_index'])==1
    x,g=q['x'],q['g'];ph=m.moments(x*m.scale)
    total=np.array([ph[0]-np.sum(g[:,1]*m.nu),ph[1]+np.sum(g[:,0]*m.eu)])
    norm=np.array([max(np.sum(abs(x)*m.Nweight),np.sum(abs(g[:,1])*m.nu),1.),max(np.sum(abs(x)*m.Eweight),np.sum(abs(g[:,0])*m.eu),1.)])
    defect=abs(total-q['ledger'])/np.maximum(norm,abs(q['ledger']))
    previous['species_balance_relative']=max(previous['species_balance_relative'],float(defect[0]))
    previous['energy_balance_relative']=max(previous['energy_balance_relative'],float(defect[1]))
    assert max(previous['species_balance_relative'],previous['energy_balance_relative'])<1e-8
    m.floor_discard=q['material_floor_discard_scaled'].copy();m.guide_g=q['restart_guide'].copy()
    for attr,key in before.prior.prior.HISTORY.items():setattr(m,attr,list(q[key]))
    checks=json.loads(str(q['joint_checks_json']));m.newton_iterations=checks['newton'];m.stage_log=checks['stages']
    m.angular_times=list(q['angular_times']);m.angular=list(q['angular'])
    trace=read(OLD/'linear-64.json')
    for group in trace['solves'][:2]:
        calls=trace['calls'][group['begin']:group['end']]
        m.max_iterations=max(m.max_iterations,sum(c['iterations'] for c in calls))
        previous['maximum_refinement_steps']=max(previous['maximum_refinement_steps'],len(calls))
        previous['initial_GMRES_stagnations']+=int(bool(calls) and calls[0]['info']!=0)
    h=q['next_step'];t=q['next_time']-h
    moment=max(moment,*(m.source(t+c*h)[2] for c in joint.C))
    previous['resume_macro']=61;previous['resume_substeps']=1
    return [q[k].copy() if k in KEYS[:7] else list(q[k]) for k in KEYS]+[moment]


def initialize():
    global restore_previous
    init=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT));init(False)
    run=owner.Model.run;restore_previous=run.__globals__['restore_substep']
    owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,restore_substep=restore_substep),argdefs=run.__defaults__)


def check():
    initialize();m=owner.Model(64);row=m.run(64,'resume-check-64',61,restart='interval-15-64')
    q=dict(np.load(OLD/'last-accepted-64.npz'));z=dict(np.load(OUT/'sweep-1/photons/resume-check-64.npz'))
    mapping=dict(x='restart_x',g='restart_g',ledger='restart_ledger',escape='restart_escape',impulse='restart_impulse',ports='restart_ports',transfer='restart_transfer',
        times='t',records='moments',port_history='radial_ports',photon_history='photon_history_scaled_occupation',gas_history='material_history',transfer_history='collision_transfer',stage_weights='accepted_angular_quadrature_weights',actual_edges='actual_step_edges',angular_times='accepted_angular_times',angular='accepted_angular_luminosity')
    fingerprint=before.prior.prior.prior.fingerprint
    for a,b in mapping.items():assert fingerprint(q[a])==fingerprint(z[b]),a
    for key in before.prior.prior.state(m):assert fingerprint(q[key])==fingerprint(z[key]),key
    assert row['actual_new_steps']==row['new_steps']==0 and row['passed']
    from fractions import Fraction as F
    for k in range(3):assert F(3,4)*F(1,3)**k+F(1,4)==F(1,k+1)
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,
        exact_saved_state_and_all_histories=True,new_physical_steps=0,accepted_actual_steps=114,
        remaining_coarse_substeps=5,remaining_fine_substeps=16,
        symbolic=dict(classification='Proven',passed=True,scope='Unchanged Radau moments0..2.')))


def evolve(n):
    factory=FunctionType(prior.factory.__code__,dict(prior.factory.__globals__,OUT=OUT))
    old=SimpleNamespace(right_solver=lambda log:factory(log,n));stable.initialize=initialize
    source=inspect.getsource(before.evolve)
    mark='    def saved_stage(m,t,h,x,g,lus):\n'
    assert source.count(mark)==1;source=source.replace(mark,mark+'        m.guide_g=g.copy()\n')
    begin=source.index('            frame=sys._getframe(1).f_locals')
    end=source.index("            write(OUT/f'failure-{n}.json'",begin)
    save=source[begin:end];source=source[:begin]+source[end:]
    save='\n'.join(line[4:] for line in save.splitlines())+'\n'
    mark='        angular_t=list(m.angular_times);angular=list(m.angular)\n'
    assert source.count(mark)==1;source=source.replace(mark,mark+save)
    ns=dict(before.evolve.__globals__,OUT=OUT,old=old);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-postfloor-evolve.py').write_text(source);ns['evolve'](n)


def audit():
    fn=FunctionType(before.audit.__code__,dict(before.audit.__globals__,OUT=OUT));fn()
    result=read(OUT/'result.json');result.update(maximum_inner_solves=12,maximum_new_Newton_proposals=8,
        previous_interrupted_proposals_preserved=True,native_primitive_precision_preserved=True,initial_guide_from_postfloor_state=True)
    write(OUT/'result.json',result)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
