"""Counterexample candidate: resume saved work with a practical Newton budget."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import complete_stable_full_interval as before

OUT=Path('native-resumed-newton195-work');OLD=before.OUT
stable,owner,joint=before.stable,before.owner,before.joint
read,write,sha=before.read,before.write,before.sha
core_initialize=stable.initialize
CAPS=dict(prepare=25,check=90,coarse=1800,fine=2700,audit=30)
KEYS=['x','g','ledger','escape','impulse','ports','transfer','times','records','port_history','photon_history','gas_history','transfer_history','stage_weights','actual_edges']


def prepare():
    assert not OUT.exists();OUT.mkdir();failure=read(OLD/'failure-64.json')
    assert failure['macro_index']==60 and failure['sub_index']==1 and failure['actual_accepted_steps']==112
    assert 'True native joint Radau equation' in failure['error']
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    reused={}
    for p in inputs:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files=inputs+[OLD/n for n in ['last-accepted-64.npz','rejected-joint-stage.npz','rejected-joint-stage.json','failure-64.json','coarse-receipt.json','linear-64.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='a5d446ae2',
        authorization='2026-09-24user explicitly asks for more generous computation budgets and continuing the experiment to reduce total cost.',
        claim='Continue the actual native coupled solve after the repaired linear systems passed but the three-Newton cap was exhausted.',
        reuse='Restore the exact last accepted substep and all its ledgers from194; do not replay that accepted substep or the preceding15/16history. Reuse the last rejected full Radau pair only as a Newton starting proposal, never as an accepted physical state.',
        resource_change='Allow8Newton proposals per physical step instead of3. The resumed rejected step already used3, so at most5new proposals remain. Keep4inner linear refinements, restart20/maxiter5 and all accuracy gates. Keep the existing time grid, fronts, physical period and two paths.',
        accuracy=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        forecast=dict(coarse_seconds=[300,1200],fine_seconds=[400,1800],
            basis='194spent185.83s on one accepted substep plus3Newton proposals of the next. Saved work is reused. Remaining7/16substeps may require several proposals; allow30/45minute caps to avoid repeated administrative interruption. Later cost remains uncertain.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        decision='Exact zero-step restoration first; then finish coarse, and only after physical acceptance finish fine and the original10-channel2percent time comparison.',
        stop='Original accuracy or constitutive/balance failure,8total Newton limit,4inner refinements,nonfinite state or generous wall cap. Preserve failures; no grid/period extension or additional parameter sweep.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused))


def restore_substep(m,previous,moment):
    q=dict(np.load(OLD/'last-accepted-64.npz'));assert int(q['macro_index'])==60 and int(q['sub_index'])==1
    m.floor_discard=q['material_floor_discard_scaled'].copy();m.guide_g=q['restart_guide'].copy()
    for attr,key in before.prior.prior.HISTORY.items():setattr(m,attr,list(q[key]))
    checks=json.loads(str(q['joint_checks_json']));m.newton_iterations=checks['newton'];m.stage_log=checks['stages']
    # The first new linear solve passed its strict gate, but its exact scalar
    # norm was not serialized. Restore a labelled upper bound, not a measurement.
    previous['maximum_extended_stage_residual']=max(previous['maximum_extended_stage_residual'],1e-14)
    previous['restored_linear_residual_is_upper_bound']=True;m.max_residual=max(m.max_residual,1e-14)
    calls=read(OLD/'linear-64.json')['calls'][:2]
    m.max_iterations=max(m.max_iterations,sum(c['iterations'] for c in calls))
    previous['maximum_refinement_steps']=max(previous['maximum_refinement_steps'],1)
    previous['initial_GMRES_stagnations']+=int(calls[0]['info']!=0)
    h=q['next_step'];t=q['next_time']-h
    moment=max(moment,*(m.source(t+c*h)[2] for c in joint.C))
    m.angular_times=list(q['angular_times']);m.angular=list(q['angular'])
    previous['resume_substeps']=1
    return [q[k].copy() if k in KEYS[:7] else list(q[k]) for k in KEYS]+[moment]


def initialize():
    stable.OUT=OUT;core_initialize()
    source=(OUT/'expanded-continuation-run.py').read_text()
    changes=[('    for k in range(begin,count):',
        '    if steps==64 and begin==60:\n        '+','.join(KEYS+['moment'])+'=restore_substep(self,previous,moment)\n    for k in range(begin,count):'),
        ('for sub in range(parts):',"for sub in range(previous.get('resume_substeps',0) if k==begin else 0,parts):"),
        ('actual_new_steps=len(stage_weights)//2-sum(1+int(v) for v in flags[:begin])',
         "actual_new_steps=len(stage_weights)//2-sum(1+int(v) for v in flags[:begin])-previous.get('resume_substeps',0)")]
    for a,b in changes:assert source.count(a)==1,a;source=source.replace(a,b)
    stage_source=(OUT/'expanded-full-stages.py').read_text()
    replacement="""    seed=getattr(m,'resume_seed',None)
    continuing=seed is not None and t==seed['time']
    limit=5 if continuing else 8
    if continuing:
        guess=seed['solution'].copy();guides=[m.unpack(row)[1].copy() for row in guess.reshape(2,dim)]
        audit=list(seed['equations'])
    for newton in range(limit):"""
    for a,b in [('    for newton in range(3):',replacement),('if newton==2:','if newton==limit-1:')]:
        assert stage_source.count(a)==1,a;stage_source=stage_source.replace(a,b)
    stage=owner.Model.run.__globals__['stages'];ns=dict(stage.__globals__);exec(compile(stage_source,__file__,'exec'),ns)
    run=owner.Model.run;space=dict(run.__globals__,stages=ns['stages'],restore_substep=restore_substep)
    exec(compile(source,__file__,'exec'),space);owner.Model.run=space['run']
    (OUT/'expanded-resumed-run.py').write_text(source);(OUT/'expanded-eight-newton-stage.py').write_text(stage_source)
    init=owner.Model.__init__
    def construct(m,n):
        init(m,n)
        if n==64:
            p=dict(np.load(OLD/'rejected-joint-stage.npz'));p['equations']=read(OLD/'rejected-joint-stage.json')['equations'];m.resume_seed=p
    owner.Model.__init__=construct


def check():
    initialize();m=owner.Model(64);row=m.run(64,'resume-check-64',60,restart='interval-15-64')
    q=dict(np.load(OLD/'last-accepted-64.npz'));z=dict(np.load(OUT/'sweep-1/photons/resume-check-64.npz'))
    mapping=dict(x='restart_x',g='restart_g',ledger='restart_ledger',escape='restart_escape',impulse='restart_impulse',ports='restart_ports',transfer='restart_transfer',
        times='t',records='moments',port_history='radial_ports',photon_history='photon_history_scaled_occupation',gas_history='material_history',transfer_history='collision_transfer',stage_weights='accepted_angular_quadrature_weights',actual_edges='actual_step_edges',angular_times='accepted_angular_times',angular='accepted_angular_luminosity')
    fingerprint=before.prior.prior.prior.fingerprint
    for a,b in mapping.items():assert fingerprint(q[a])==fingerprint(z[b]),a
    for key in before.prior.prior.state(m):assert fingerprint(q[key])==fingerprint(z[key]),key
    assert row['actual_new_steps']==row['new_steps']==0 and row['passed']
    from fractions import Fraction as F
    for k in range(3):assert F(3,4)*F(1,3)**k+F(1,4)==F(1,k+1)
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,new_physical_steps=0,
        all_internal_states_and_ledgers_exact=True,accepted_actual_steps=112,remaining_coarse_substeps=7,
        symbolic=dict(classification='Proven',passed=True,scope='Unchanged exact Radau moments0..2.'),final_charge_conclusion='unadjudicated'))


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];stable.initialize=initialize
            evolve=FunctionType(before.evolve.__code__,dict(before.evolve.__globals__,OUT=OUT));evolve(64 if action=='coarse' else 128)
        elif action=='audit':
            audit=FunctionType(before.audit.__code__,dict(before.audit.__globals__,OUT=OUT));audit()
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
