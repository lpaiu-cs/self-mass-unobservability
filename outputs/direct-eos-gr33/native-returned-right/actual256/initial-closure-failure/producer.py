"""Finish the actual fine return using the existing true-residual solver.

The accepted coarse path and fine canonical15 checkpoint are immutable inputs.
Original physical gates remain authoritative; right preconditioning changes
only the proposal. Reuse the saved failing system instead of solving it twice.
"""
from pathlib import Path
from types import SimpleNamespace
import gc,inspect,os,resource,sys,time
import numpy as np
from scipy.sparse.linalg import gmres as scipy_gmres
import continue_returned_krylov as prior
import right_precondition_full_interval as right

OUT=Path('native-returned-right256-work');OLD=prior.OUT;base=prior.prior.prior
joint=prior.prior.joint
read,write,sha,bind=prior.read,prior.write,prior.sha,prior.bind
CAPS=dict(prepare=600,check=3600,fine=21600,audit=600)


def prepare():
    assert read(OLD/'controller-status.json')['state']=='failed'
    assert read(OLD/'run-64.json')['passed']
    assert read(OLD/'capture-128.json')['actual_steps']==228
    assert 'Four-moment linear residual' in read(OLD/'fine-receipt.json')['error']
    assert not OUT.exists();OUT.mkdir()
    files=list((OLD/'sweep-0').rglob('*.npz'))
    files += [p for d in ['metric','gr'] for p in (OLD/d).iterdir() if p.is_file()]
    files += [p for p in (OLD/'sweep-1/photons').iterdir() if p.is_file() and p.name.endswith(('-64.npz','-64.json','-128.npz','-128.json'))]
    files += [OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json','symbolic.json','metric-result.json','metric-receipt.json','coarse-receipt.json','run-64.json','checks-64.json','recovered-64.npz','recovered-128.npz']]
    for p in files:
        dst=OUT/p.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(p,dst)
    (OUT/'sweep-1/material').mkdir(parents=True,exist_ok=True)
    files += [OLD/n for n in ['fine-receipt.json','controller-status.json','linear-232.json','krylov-expansion.json','last-pair-128.npz','original-limit-232-30.npz']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    p=read(OLD/'plan.json');p.update(classification='Conjectural',
        claim='Complete the actual original119/231 returned pair with unchanged stage and time gates, then read its own charge.',
        failure=read(OLD/'fine-receipt.json'),
        repair='Reuse existing191right preconditioning A P y=b-Ax0 with restart80,maxiter20 and at most12 extended refinements. Original254solver remains first for other systems. Exactly matched saved229system reuses the accepted right-solve proposal. All proposals pass the original long-double vector1e-14 and four-moment1e-13 gates and original nonlinear/constitutive/conservation/time gates.',
        evidence='Fine completed228/231; at229the final left-preconditioned callback was0 while true correction residual1.91676 and final stage-linear residual9.50075e-11. This is an accuracy failure, not a wall or memory stop. More identical left iterations are not an established repair.',
        reuse='Reuse complete coarse119 and fine215 canonical checkpoint. Revalidate only saved fine stages200..215 for omitted native/branch diagnostics. Replay216..228 because their full runtime state and histories were not serialized; require last228pair exactly equal. Do not rerun the primary, metric, coarse or earlier fine215.',
        forecast='Previous fine32remaining steps took2906.52s before failure. Replaying13late steps is assumed25..45min from that measurement; right correction costs are unmeasured. Allow1hour saved-system check and6hours actual fine,16GiB,CPU3, without additional clocks or period.',
        decision='Stop on any original gate, identity mismatch or generous cap; no grid/horizon/tolerance enlargement. A saved-system check is preparatory and feeds the actual continuation automatically.',
        budgets=CAPS,scientific_gates_changed=False,full_goal_complete=False,final_charge_conclusion='unadjudicated')
    p['bindings'].update({str(p):sha(p) for p in dict.fromkeys(files)});write(OUT/'plan.json',p)
    import sympy as s
    a,p,y,b,x=s.symbols('a p y b x');assert s.expand((b-a*(x+p*y))-((b-a*x)-a*p*y))==0
    write(OUT/'right-symbolic.json',dict(classification='Proven',passed=True,scope='Right correction leaves the original linear equation unchanged; no physical or discretization error bound.'))


def right_solve(m,op,P,rhs,guess):
    calls=[];template=right.right_solver(calls)
    inner=bind(template.__globals__['gmres'],joint=SimpleNamespace(gmres=scipy_gmres))
    source=inspect.getsource(joint.solve)
    for a,b in [('restart=20,maxiter=5','restart=80,maxiter=20'),('range(4)','range(12)'),('assert k<3,','assert k<11,')]:
        assert source.count(a)==1;source=source.replace(a,b)
    namespace=dict(joint.solve.__globals__,gmres=inner);exec(compile(source,__file__,'exec'),namespace)
    try:return namespace['solve'](m,op,P,rhs,guess)
    finally:write(OUT/'right-calls.json',dict(classification='Counterexample candidate',calls=calls))


def check():
    Model=bind(base.initialize,OUT=OUT)();m=Model(128)
    saved=dict(np.load(OLD/'original-limit-232-30.npz'));last=dict(np.load(OLD/'last-pair-128.npz'))
    flags=np.load(OLD/'sweep-1/photons/interval-15-128.npz')['split_macro_steps'];clock=[];base_h=m.t[-1]/128
    for k in range(128):
        parts=1+int(flags[k]);step=base_h/parts
        clock.extend((k*base_h+sub*step,step) for sub in range(parts))
    t,h=clock[228];assert t==m.anchor['actual_step_edges'][228]
    x=last['photons'][-1].copy();g=last['gas'][-1].copy();m.guide_g=g.copy();g[~m.material.active(t)]=0.
    assert abs(last['time']+last['step']-t)<1e-18
    lus=[joint.splu(joint.sparse.eye(m.n*m.q,format='csc')-a*h*m.A) for a in [5/12,1/4]]
    captured={}
    class Captured(Exception):pass
    def intercept(m,op,P,rhs,guess):captured.update(op=op,P=P,rhs=rhs,guess=guess);raise Captured()
    stage=bind(Model.run.__globals__['stages'],solve=intercept)
    try:stage(m,t,h,x,g,lus)
    except Captured:pass
    assert captured
    mapping=dict(classification='Counterexample candidate',rhs_exact=np.array_equal(captured['rhs'],saved['rhs']),guess_exact=np.array_equal(captured['guess'],saved['guess']),
        residual_exact=np.array_equal(saved['rhs']-captured['op'].matvec(saved['solution']),saved['residual']),actual_step=229,new_physical_steps=0)
    write(OUT/'system-reconstruction.json',mapping);assert all(mapping[k] for k in ['rhs_exact','guess_exact','residual_exact']),mapping
    answer=right_solve(m,captured['op'],captured['P'],saved['rhs'],saved['solution'])
    residual=saved['rhs']-captured['op'].matvec(answer)
    relative=float(np.linalg.norm(residual)/np.linalg.norm(saved['rhs']));moments=joint.physical_norm(m,residual)/joint.scales(m,saved['rhs'],answer)
    result=dict(classification='Counterexample candidate',passed=relative<1e-14 and max(moments)<1e-13,relative=relative,moments=[float(v) for v in moments],actual_step=229,physical_stage_accepted=False,final_charge_conclusion='unadjudicated')
    assert result['passed'];np.savez_compressed(OUT/'accepted-linear.npz',rhs=saved['rhs'],guess=saved['guess'],solution=answer)
    write(OUT/'right-result.json',result)


def initialize():
    Model=bind(prior.initialize,OUT=OUT)();run=Model.run;stage=run.__globals__['stages'];original=stage.__globals__['solve']
    cached=dict(np.load(OUT/'accepted-linear.npz'))
    def solve(m,op,P,rhs,guess):
        if np.array_equal(rhs,cached['rhs']) and np.array_equal(guess,cached['guess']):
            answer=cached['solution'].copy();residual=rhs-op.matvec(answer)
            relative=float(np.linalg.norm(residual)/np.linalg.norm(rhs));moments=joint.physical_norm(m,residual)/joint.scales(m,rhs,answer)
            assert relative<1e-14 and max(moments)<1e-13
            write(OUT/'actual-cache-reuse.json',dict(classification='Counterexample candidate',rhs_guess_exact=True,relative=relative,moments=[float(v) for v in moments]))
            m.max_residual=max(m.max_residual,relative);return answer
        try:return original(m,op,P,rhs,guess)
        except AssertionError as exc:
            if not str(exc).startswith("('Four-moment linear residual'"):raise
            tb=exc.__traceback__
            while tb and not {'sol','residual','rhs'}.issubset(tb.tb_frame.f_locals):tb=tb.tb_next
            assert tb is not None;failed=tb.tb_frame.f_locals
            np.savez_compressed(OUT/'new-left-failure.npz',rhs=rhs,guess=guess,solution=failed['sol'],residual=failed['residual'])
            write(OUT/'new-left-failure.json',dict(error=repr(exc)))
            return right_solve(m,op,P,rhs,failed['sol'])
    Model.run=bind(run,stages=bind(stage,solve=solve));return Model


def seed_row(Model):
    row=read(base.OLD/'run-128.json');cut=read(OLD/'plan.json')['restart_time_seconds']
    native=[v for v in row['anchor_checks'] if v['time']<=cut+1e-18]
    branches=[v for v in row['branch_checks'] if v['time']<=cut+1e-18]
    m=Model(128);z=np.load(OUT/'sweep-1/photons/interval-15-128.npz');errors=[]
    for t,q,r in zip(z['joint_stage_times'],z['joint_stage_conserved_scaled'],z['joint_native_rates_scaled']):
        if t<=cut+1e-18:continue
        g=np.column_stack([(q[2]-m.kappa*q[0])/m.eu,q[3]/m.nu,q[0]/m.bu,q[1]/m.su])
        actual=m.native(t,g)*m.units
        errors.append(float(np.max(np.sum(abs(actual-r),axis=0)/np.maximum(np.sum(abs(r),axis=0),joint.LD('1e-290')))))
    assert max(errors)<1e-12,errors
    native+=m.anchor_checks;branches+=m.branch_checks
    write(OUT/'prefix-native.json',dict(classification='Counterexample candidate',passed=True,stored_stage_native_relative=errors,physical_steps_recomputed=0))
    row.update(anchor_checks=native,branch_checks=branches,intervals=row['intervals'][:14]+[read(OUT/'sweep-1/photons/interval-15-128.json')]);del m;gc.collect();return row


def fine():
    Model=initialize();row=seed_row(Model)
    source=inspect.getsource(base.evolve)
    changes=[("interval-14-{n}.npz","interval-15-{n}.npz"),("oldrow=read(OLD/f'run-{n}.json')","oldrow=seed"),("oldrow['intervals'][:14]","oldrow['intervals'][:15]"),("restart=f'interval-14-{n}'","restart=f'interval-15-{n}'"),('for j in [15,16]:','for j in [16]:')]
    for a,b in changes:assert source.count(a)==1,(a,source.count(a));source=source.replace(a,b)
    marker="        write(OUT/f'capture-{n}.json',dict(actual_steps=len(moments)//2,new_steps=(len(moments)-count)//2))"
    source=base.base.replace(source,marker,marker+"\n        if len(moments)//2==228: replay(t,h,x,g,pair)")
    def replay(t,h,x,g,pair):
        with np.load(OLD/'last-pair-128.npz') as old:
            for k,v in dict(time=t,step=h,x_initial=x,g_initial=g,photons=[v[0] for v in pair],gas=[v[1] for v in pair]).items():assert np.array_equal(old[k],v),('Accepted228pair changed',k)
        write(OUT/'replay-228.json',dict(classification='Counterexample candidate',passed=True,actual_pair_exact=True))
    namespace=dict(base.evolve.__globals__,OUT=OUT,OLD=OLD,initialize=lambda:Model,seed=row,replay=replay)
    exec(compile(source,__file__,'exec'),namespace);namespace['evolve'](128)


def audit():
    bind(prior.audit,OUT=OUT)();r=read(OUT/'result.json');r.update(original254fine_failure_preserved=True,coarse119_reused=True,fine215_reused=True)
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/f'{action}-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
