"""Counterexample candidate: apply the passing lattice solver to actual evolution."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import continue_precise_momentum as prior
import solve_native_rounding_neighbours as rounding

OUT=Path('native-rounding-continuation218-work');OLD=prior.OUT
old,legacy,owner,joint,LD=prior.old,prior.legacy,prior.owner,prior.joint,prior.LD
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,check=300,coarse=7200,fine=10800,audit=180)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    admission=read(rounding.OUT/'refined-result.json')
    assert admission['linear_passed'] and not admission['actual_stage_passed']
    assert read(OLD/'failure-64.json')['actual_accepted_steps']==116
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    files += [OLD/f'sweep-1/photons/interval-15-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    reused={}
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    for src in files:
        dst=OUT/src.relative_to(OLD);os.link(src,dst);reused[str(dst)]=sha(src)
    files += [OLD/n for n in ['last-accepted-64.npz','failed-linear-64.npz','failure-64.json','coarse-receipt.json','linear-64.json','stage-progress-64.json']]
    files += [rounding.OUT/n for n in ['pair-producer.py','plan.json','result.json','check-receipt.json','refinement-plan.json','refined-result.json','refined-proposal.npz','refine-receipt.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply the demonstrated same-operator numerical repair to the remaining actual coupled solution; only the original nonlinear, physical and paired-time gates can accept it.',
        evidence='210accepted116steps then failed117linear4.093584e-10.215fine correction reduced3.301275e-13;216full collision arithmetic changed nothing.217representable pairs reduced1.668182e-14; triples pass7.980254e-15. The actual nonlinear3.788254e-7still fails: this is only a Newton proposal.',
        method='Reuse210B/S precision operator and202true native B. Extend the existing near-converged corrector with217paired then triple representable-coordinate corrections, verifying the unchanged full linear operator and physical norms. Eligibility1e-8is only a corrector trigger, not an acceptance tolerance.',
        reuse='Restore116accepted states and every history exactly. Seed the first117linear solution from217only after exact equality of the original210RHS and branch guess. Preserve all prior failed states. No accepted physical prefix replay.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,max_Newton=8,max_linear=12,GMRES=dict(restart=80,maxiter=10),
        forecast='210used3889.54s before117linear failure.217pair repair including reconstruction68.66s; triple receipt is bound below. Allow2/3hours for3coarse/16fine physicalsteps. Later/fine speed remains unmeasured; caps are not completion estimates.',
        stop='Any exact restart/seed failure, unchanged physical/accuracy gate,8Newton/12linear limit or wall cap. Save accepted state before each actual step and failed linear proposals. No spatial, temporal or path expansion.',
        decision='Complete coarse/fine and original10channel2percent comparison before full-period GR readout. Prior dense photon recovery,GR source,EOS derivative,self-GR and final charge failures remain.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as sp
    A,E,x,d,b=sp.symbols('A E x d b',commutative=False)
    assert sp.expand(b-A*(x+E*d)-((b-A*x)-A*E*d))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Column corrections preserve the same exact linear equation. No convergence guarantee or physical error bound.'))


def initialize(seed=True):
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    run=owner.Model.run;source=inspect.getsource(old.restore_substep)
    changes=[("int(q['macro_index'])==61 and int(q['sub_index'])==1","int(q['macro_index'])==62 and int(q['sub_index'])==1"),
        ("previous['resume_macro']=61;previous['resume_substeps']=1","previous['resume_macro']=62;previous['resume_substeps']=1")]
    for a,b in changes:assert source.count(a)==1,a;source=source.replace(a,b)
    ns=dict(old.restore_substep.__globals__,OUT=OUT,OLD=OLD,restore_previous=run.__globals__['restore_substep']);exec(compile(source,__file__,'exec'),ns)
    owner.Model.run=FunctionType(run.__code__,dict(run.__globals__,restore_substep=ns['restore_substep']),argdefs=run.__defaults__)
    if seed:
        constructor=owner.Model.__init__
        def construct(m,n):
            constructor(m,n)
            if n==64:
                z=np.load(OLD/'failed-linear-64.npz');p=np.load(OLD/'last-accepted-64.npz')
                m.resume_seed=dict(time=p['next_time'][()],solution=z['guess'].copy(),equations=[])
        owner.Model.__init__=construct


def check():
    source=inspect.getsource(old.check).replace("m.run(64,'resume-check-64',61,","m.run(64,'resume-check-64',62,").replace('accepted_actual_steps=114','accepted_actual_steps=116').replace('remaining_coarse_substeps=5','remaining_coarse_substeps=3')
    ns=dict(old.check.__globals__,OUT=OUT,OLD=OLD,initialize=lambda:initialize(False));exec(compile(source,__file__,'exec'),ns);ns['check']()


def factory(log,n):
    fine=FunctionType(prior.polish_function.__code__,dict(prior.polish_function.__globals__,OUT=OUT))()
    def polish(m,op,rhs,sol):
        sol=fine(m,op,rhs,sol);records=[]
        for triplets in [False,True]:
            if np.linalg.norm(rhs-op.matvec(sol))<LD('1e-14')*max(np.linalg.norm(rhs),LD('1e-290')):break
            sol=rounding.neighbours(m,op,rhs,sol,records,triplets=triplets)
        path=OUT/f'rounding-{n}.json';previous=read(path) if path.exists() else []
        write(path,previous+[records]);return sol
    source=inspect.getsource(legacy.factory).replace('first_system=[True]','first_system=[False]').replace('relative<1e-10','relative<1e-8')
    ns=dict(vars(legacy),OUT=OUT,OLD=OLD,polish=polish);exec(compile(source,__file__,'exec'),ns)
    solve=ns['factory'](log,n);first=[True]
    def seeded(m,op,P,rhs,guess):
        if n==64 and first[0]:
            first[0]=False;z=np.load(OLD/'failed-linear-64.npz');proposal=np.load(rounding.OUT/'refined-proposal.npz')
            assert np.array_equal(rhs,z['rhs']) and np.array_equal(guess,z['guess']),'Stored branch and RHS changed'
            assert np.array_equal(rhs-op.matvec(proposal['solution']),proposal['residual']),'Stored operator changed'
            guess=proposal['solution'].copy()
            write(OUT/'linear-seed-identity.json',dict(classification='Counterexample candidate',passed=True,exact_RHS_and_branch_guess_and_residual=True,accepted_as_physical_state=False))
        return solve(m,op,P,rhs,guess)
    return seeded


def evolve(n):
    FunctionType(prior.evolve.__code__,dict(prior.evolve.__globals__,OUT=OUT,initialize=initialize,factory=factory))(n)


def audit():
    FunctionType(prior.audit.__code__,dict(prior.audit.__globals__,OUT=OUT))()
    r=read(OUT/'result.json');r.update(representable_coordinate_solver_applied=True,final_charge_conclusion='unadjudicated',full_goal_complete=False);write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
