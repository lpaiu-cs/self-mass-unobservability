"""Apply the saved118integer correction to the actual unchanged coupled solve."""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import solve_saved_integer_native as integer

prior=integer.evolution;OUT=Path('native-integer-continuation229-work');OLD=prior.OUT
old,legacy,owner,joint,LD=prior.old,prior.legacy,prior.owner,prior.joint,prior.LD
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,check=300,coarse=7200,fine=10800,audit=180)


def prepare():
    assert read(integer.OUT/'result.json')['linear_passed'];assert read(OLD/'failure-64.json')['actual_accepted_steps']==117
    assert not OUT.exists();OUT.mkdir();files=[]
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/interval-15-{n}{ext}' for n in [64,128] for ext in ['.npz','.json']]
    for src in inputs:
        dst=OUT/src.relative_to(OLD);dst.parent.mkdir(parents=True,exist_ok=True);os.link(src,dst);files += [src,dst]
    (OUT/'sweep-1/material').mkdir(exist_ok=True)
    files += [OLD/n for n in ['last-accepted-64.npz','failed-linear-64.npz','failure-64.json','coarse-receipt.json','linear-64.json','rounding-64.json']]
    files += [integer.OUT/n for n in ['plan.json','result.json','proposal.npz','solve-receipt.json','integer-search.json','full-operator-checks.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Apply228passing integer proposal to actual118and continue the same declared full-period coupled solution at every original gate.',
        evidence='218accepted117then timed out after7202s.228exactly reconstructed its RHS/residual and passed linear5.08025e-15in40.4s. Actual nonlinear5.50329e-10still fails1e-12; the proposal is not an accepted state.',
        method='Restore117state and every existing history without replay. Seed the118linear solve only after exact RHS/guess/residual checks. Reuse210continuous correction, then joint bounded integer-ULP correction on later unresolved linear systems. Full original nonlinear stage, physical, native, constitutive, conservation and paired-time gates remain.',
        gates=read(OLD/'plan.json')['gates'],budgets=CAPS,CPU_threads=1,virtual_GiB=8,max_Newton=8,max_linear=12,GMRES=dict(restart=80,maxiter=10),
        integer_search=dict(max_ULPs=1024,seconds_per_search=90,max_nodes=5000,maximum_improving_passes=4),
        forecast='Saved-system reconstruction plus integer correction40.4s, actual118cost unknown. Allow2hours for2coarse remainingsteps and3hours for16fine. This is a new demonstrated solver correction, not an unchanged timeout retry. No accepted physical prefix replay, new grid/period/path or tolerance change.',
        stop='Any exact restart/proposal identity, original physical/accuracy gate,8Newton/12linear limit or wall cap. Save last accepted state and failed proposals.227stage-GR and223/224pipelines stay unchanged.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(integer.OUT/'symbolic.json'))


def initialize(seed=True):
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    run=owner.Model.run;s=inspect.getsource(old.restore_substep)
    for a,b in [("int(q['macro_index'])==61 and int(q['sub_index'])==1","int(q['macro_index'])==63 and int(q['sub_index'])==0"),
                ("trace['solves'][:2]","trace['solves'][:len(checks['newton'][-1])]"),
                ("previous['resume_macro']=61;previous['resume_substeps']=1","previous['resume_macro']=63;previous['resume_substeps']=0")]:
        assert s.count(a)==1;s=s.replace(a,b)
    ns=dict(old.restore_substep.__globals__,OUT=OUT,OLD=OLD,restore_previous=run.__globals__['restore_substep']);exec(compile(s,__file__,'exec'),ns)
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
    s=inspect.getsource(old.check).replace("m.run(64,'resume-check-64',61,","m.run(64,'resume-check-64',63,").replace('accepted_actual_steps=114','accepted_actual_steps=117').replace('remaining_coarse_substeps=5','remaining_coarse_substeps=2')
    ns=dict(old.check.__globals__,OUT=OUT,OLD=OLD,initialize=lambda:initialize(False));exec(compile(s,__file__,'exec'),ns);ns['check']()


def factory(log,n):
    fine=FunctionType(prior.prior.polish_function.__code__,dict(prior.prior.polish_function.__globals__,OUT=OUT))()
    source=inspect.getsource(prior.rounding.neighbours);a=source.index('        guesses={};count=');b=source.index('        best=sol;',a)
    source=source[:a]+"        order=list(range(len(hard)));guesses=integer_guesses(PH,pr,goal,integer_log)\n"+source[b:]
    def polish(m,op,rhs,sol):
        sol=fine(m,op,rhs,sol);records=[];integer_log=[]
        ns=dict(prior.rounding.neighbours.__globals__,integer_guesses=integer.integer_guesses,integer_log=integer_log)
        exec(compile(source,__file__,'exec'),ns);sol=ns['neighbours'](m,op,rhs,sol,records)
        path=OUT/f'integer-{n}.json';previous=read(path) if path.exists() else []
        write(path,previous+[dict(searches=integer_log,full_operator_checks=records)]);return sol
    s=inspect.getsource(legacy.factory).replace('first_system=[True]','first_system=[False]').replace('relative<1e-10','relative<1e-8')
    ns=dict(vars(legacy),OUT=OUT,OLD=OLD,polish=polish);exec(compile(s,__file__,'exec'),ns);solve=ns['factory'](log,n);first=[True]
    def seeded(m,op,P,rhs,guess):
        if n==64 and first[0]:
            first[0]=False;z=np.load(OLD/'failed-linear-64.npz');proposal=np.load(integer.OUT/'proposal.npz')
            assert np.array_equal(rhs,z['rhs']) and np.array_equal(guess,z['guess'])
            assert np.array_equal(rhs-op.matvec(proposal['solution']),proposal['residual'])
            guess=proposal['solution'].copy();write(OUT/'proposal-identity.json',dict(classification='Counterexample candidate',passed=True,exact_RHS_branch_and_residual=True,physical_state_accepted=False))
        return solve(m,op,P,rhs,guess)
    return seeded


def evolve(n):FunctionType(prior.prior.evolve.__code__,dict(prior.prior.evolve.__globals__,OUT=OUT,initialize=initialize,factory=factory))(n)
def audit():FunctionType(prior.prior.audit.__code__,dict(prior.prior.audit.__globals__,OUT=OUT))()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    os.sched_setaffinity(0,{min(os.sched_getaffinity(0))});resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];evolve(64 if action=='coarse' else 128)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
