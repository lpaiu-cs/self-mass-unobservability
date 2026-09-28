"""Counterexample candidate: balance all acceptance channels in the Krylov solve.

Reuse the accepted substeps and the interrupted linear iterate. Change numerical
coordinates and work limits, never the physical equations or accuracy gates.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
import resume_full_native_newton as prior

OUT=Path('native-balanced-krylov196-work');OLD=prior.OUT;EARLIER=prior.OLD
before,stable,owner,joint=prior.before,prior.stable,prior.owner,prior.joint
LD=joint.LD;read,write,sha=prior.read,prior.write,prior.sha;KEYS=prior.KEYS
CAPS=dict(prepare=30,check=100,coarse=1800,fine=2700,audit=30)


def prepare():
    assert not OUT.exists();OUT.mkdir();failure=read(OLD/'failure-64.json')
    assert failure['actual_accepted_steps']==113 and failure['macro_index']==61 and failure['sub_index']==0
    assert read(OLD/'controller-status.json')['state']=='failed'
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    files += [OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    reused={}
    for p in files:dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['last-accepted-64.npz','failed-linear-64.npz','failure-64.json','coarse-receipt.json','linear-64.json','stage-progress-64.json']]
    files += [EARLIER/'last-accepted-64.npz']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='a48226ed3',
        claim='Complete the remaining actual coupled interval using a Krylov norm that resolves the original vector and all four physical acceptance channels together.',
        evidence='The expanded Newton solve accepted the old failed actual step at5.43115e-13, but the next step exhausted4linear refinements with vector1.56758e-9 and photon/material physical errors0.1162/0.0150. The last inner residual reported0 while its true correction residual was0.347. This is a failure of the original proposal/work budget, not acceptance or a wall timeout.',
        coordinates='For fixed positive diagonal D, solve D*A*P*D^-1*y=D*(b-A*x0), then x=x0+P*D^-1*y. D uses the maximum of the original normalized-vector weight and each physical moment weight divided by its unchanged gate and scale. The final acceptance remains in the original coordinates, with consistent80digit baryon-row arithmetic.',
        work_budget='Per original user resource expansion, use restart80/maxiter10, up to12linear refinements, and the already-authorized8Newton proposals. Preserve all accuracy gates, original two time paths, fronts, source and horizon. The stored failed linear iterate is reused only after exact original RHS/branch-guess identity.',
        reuse='Restore113accepted substeps and their complete physical histories. Only6coarse and16fine substeps remain. Recover energy/species maxima at both extra saved endpoints; do not omit their balance checks on restart.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        forecast=dict(coarse_seconds=[300,1500],fine_seconds=[500,2200],basis='Prior actual continuation371.08s included5new native proposals of the repaired step and4failed linear attempts at the next. Both accepted work and the failed linear iterate are saved. Larger Krylov bases cost more per iteration; no speed guarantee is inferred.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        stop='Original physical/accuracy gate,8Newton/12linear limits,nonfinite scaling or30/45minute wall cap. No new time grid,period or parameter sweep.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as sp
    A,P,D,b,x,y=sp.symbols('A P D b x y',nonzero=True)
    assert sp.expand(D*(A*(x+P*y/D)-b)-(D*A*P*y/D-D*(b-A*x)))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Invertible diagonal row scaling and right preconditioning preserve the linear equation; this is not a conditioning or physical convergence theorem.'))


def restore_substep(m,previous,moment):
    q=dict(np.load(OLD/'last-accepted-64.npz'));assert int(q['macro_index'])==61 and int(q['sub_index'])==0
    # Recompute the omitted diagnostic extrema at BOTH extra accepted endpoints.
    for folder in [EARLIER,OLD]:
        with np.load(folder/'last-accepted-64.npz') as z:
            x,g=z['x'],z['g'];photon=m.moments(x*m.scale)
            total=np.array([photon[0]-np.sum(g[:,1]*m.nu),photon[1]+np.sum(g[:,0]*m.eu)])
            norm=np.array([max(np.sum(abs(x)*m.Nweight),np.sum(abs(g[:,1])*m.nu),1.),max(np.sum(abs(x)*m.Eweight),np.sum(abs(g[:,0])*m.eu),1.)])
            defect=abs(total-z['ledger'])/np.maximum(norm,abs(z['ledger']))
            previous['species_balance_relative']=max(previous['species_balance_relative'],float(defect[0]))
            previous['energy_balance_relative']=max(previous['energy_balance_relative'],float(defect[1]))
            h=z['next_step'];t=z['next_time']-h;moment=max(moment,*(m.source(t+c*h)[2] for c in joint.C))
    assert max(previous['species_balance_relative'],previous['energy_balance_relative'])<1e-8
    m.floor_discard=q['material_floor_discard_scaled'].copy();m.guide_g=q['restart_guide'].copy()
    for attr,key in before.prior.prior.HISTORY.items():setattr(m,attr,list(q[key]))
    checks=json.loads(str(q['joint_checks_json']));m.newton_iterations=checks['newton'];m.stage_log=checks['stages']
    m.angular_times=list(q['angular_times']);m.angular=list(q['angular'])
    previous['maximum_extended_stage_residual']=max(previous['maximum_extended_stage_residual'],1e-14)
    previous['restored_linear_residual_is_upper_bound']=True;m.max_residual=max(m.max_residual,1e-14)
    m.max_iterations=max(m.max_iterations,400);previous['maximum_refinement_steps']=max(previous['maximum_refinement_steps'],3)
    previous['resume_substeps']=0;previous['resume_macro']=61
    return [q[k].copy() if k in KEYS[:7] else list(q[k]) for k in KEYS]+[moment]


def initialize():
    source=inspect.getsource(prior.initialize).split('    init=owner.Model.__init__',1)[0]
    ns=dict(vars(prior),OUT=OUT,OLD=OLD,restore_substep=restore_substep);exec(compile(source,__file__,'exec'),ns);ns['initialize']()
    source=(OUT/'expanded-resumed-run.py').read_text()
    a='=restore_substep(self,previous,moment)\n';assert source.count(a)==1
    source=source.replace(a,a+"        begin=previous['resume_macro'];error=previous['energy_balance_relative'];species_error=previous['species_balance_relative']\n")
    run=owner.Model.run;ns=dict(run.__globals__);exec(compile(source,__file__,'exec'),ns);owner.Model.run=ns['run']
    (OUT/'expanded-resumed-run.py').write_text(source)


def balanced_factory(log,n):
    def solve(m,op,P,rhs,guess):
        branch_guess=guess.copy()
        if n==64 and not log:
            p=dict(np.load(OLD/'failed-linear-64.npz'))
            assert np.array_equal(rhs,p['rhs']) and np.array_equal(guess,p['guess']),'Saved failed system changed'
            guess=p['solution'].copy();write(OUT/'linear-seed-identity.json',dict(classification='Counterexample candidate',passed=True,exact_RHS_and_branch_guess=True,accepted_as_physical_state=False))
        scales=joint.scales(m,rhs,guess);base=LD('1e14')/max(np.linalg.norm(rhs),LD('1e-290'))
        photon=np.maximum(m.Nweight/scales[0],m.Eweight/scales[1])*LD('1e13')
        gas=m.units/np.array([scales[1],scales[0],scales[2],scales[3]])*LD('1e13')
        weights=np.tile(np.maximum(m.pack(photon,gas),base),2);D=np.asarray(weights/np.max(weights),float)
        assert np.all(np.isfinite(D)) and np.all(D>0)
        def krylov(actual,b,**kwargs):
            pre=kwargs.pop('M');initial=kwargs.pop('x0');callback=kwargs.pop('callback');history=[]
            x0=np.zeros_like(b,dtype=LD) if initial is None else np.asarray(initial,LD)
            residual=np.asarray(b,LD)-actual.matvec(x0)
            system=joint.LinearOperator(actual.shape,lambda v:D*actual.matvec(pre.matvec(v/D)),dtype=float)
            def record(v):history.append(float(v));callback(v)
            start=time.monotonic();kwargs.update(restart=80,maxiter=10)
            y,info=joint.gmres(system,np.asarray(D*residual,float),x0=None,callback=record,**kwargs)
            answer=x0+np.asarray(pre.matvec(y/D),LD);r=np.asarray(b,LD)-actual.matvec(answer)
            log.append(dict(info=int(info),seconds=time.monotonic()-start,iterations=len(history),history=history,
                true_relative=float(np.linalg.norm(r)/max(np.linalg.norm(b),1e-290)),weight_min=float(min(D)),weight_max=float(max(D))))
            return answer,info
        source=inspect.getsource(joint.solve)
        for a,b in [('range(4)','range(12)'),('assert k<3,','assert k<11,'),('gmres(op,np.asarray(rhs,float),x0=np.asarray(guess,float)','gmres(op,rhs,x0=guess')]:
            assert source.count(a)==1,a;source=source.replace(a,b)
        ns=dict(joint.solve.__globals__,gmres=krylov);exec(compile(source,__file__,'exec'),ns)
        (OUT/'expanded-balanced-linear.py').write_text(source)
        try:return ns['solve'](m,op,P,rhs,guess)
        except BaseException as exc:
            tb=exc.__traceback__
            while tb and tb.tb_frame.f_code!=ns['solve'].__code__:tb=tb.tb_next
            if tb and 'sol' in tb.tb_frame.f_locals:
                sol=tb.tb_frame.f_locals['sol'];np.savez_compressed(OUT/f'failed-linear-{n}.npz',rhs=rhs,guess=branch_guess,initial_iterate=guess,solution=sol,residual=rhs-op.matvec(sol))
            raise
    return solve


def check():
    source=inspect.getsource(prior.check).replace("m.run(64,'resume-check-64',60,","m.run(64,'resume-check-64',61,").replace('accepted_actual_steps=112','accepted_actual_steps=113').replace('remaining_coarse_substeps=7','remaining_coarse_substeps=6')
    ns=dict(vars(prior),OUT=OUT,OLD=OLD,initialize=initialize);exec(compile(source,__file__,'exec'),ns);ns['check']()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];stable.initialize=initialize;n=64 if action=='coarse' else 128
            old=SimpleNamespace(right_solver=lambda log:balanced_factory(log,n))
            evolve=FunctionType(before.evolve.__code__,dict(before.evolve.__globals__,OUT=OUT,old=old));evolve(n)
        elif action=='audit':
            audit=FunctionType(before.audit.__code__,dict(before.audit.__globals__,OUT=OUT));audit()
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
