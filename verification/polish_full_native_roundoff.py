"""Counterexample candidate: repair representable residuals in the actual solve.

The saved failed iterate has two dominant baryon residuals below one ULP of
their large state components. Minimize the unchanged full residual through
better resolved, coupled material coordinates, then apply every original gate.
"""
from pathlib import Path
from types import FunctionType,SimpleNamespace
import inspect,json,os,resource,sys,time
import numpy as np
from scipy import sparse
import balance_full_native_krylov as prior

OUT=Path('native-roundoff-polish197-work');OLD=prior.OUT
before,stable,owner,joint=prior.before,prior.stable,prior.owner,prior.joint
LD=joint.LD;read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=60,coarse=1800,fine=2700,audit=60)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    assert read(OLD/'controller-status.json')['state']=='failed'
    assert read(OLD/'failure-64.json')['actual_accepted_steps']==113
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    inputs=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    inputs += [OLD/f'sweep-1/photons/interval-15-{n}{suffix}' for n in [64,128] for suffix in ['.npz','.json']]
    reused={}
    for p in inputs:dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    # No physical state changed during196; prove that its restart is still195.
    left=dict(np.load(OLD/'last-accepted-64.npz'));right=dict(np.load(prior.OLD/'last-accepted-64.npz'))
    fingerprint=before.prior.prior.prior.fingerprint
    assert left.keys()==right.keys()
    for k in left:assert fingerprint(left[k])==fingerprint(right[k]),k
    np.savez_compressed(OUT/'saved-residual-input.npz',**dict(np.load(OLD/'failed-linear-64.npz')))
    files=inputs+[OLD/n for n in ['failed-linear-64.npz','last-accepted-64.npz','failure-64.json','coarse-receipt.json','linear-64.json','restart-check.json']]
    files += [prior.OLD/'last-accepted-64.npz',prior.EARLIER/'last-accepted-64.npz']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='e3d677ce3',
        claim='Finish the actual coupled interval by reducing the unchanged full linear residual through representable material coordinates, reusing the last196iterate.',
        evidence='After12balanced corrections the physical errors are below1.81e-19, but the vector7.230025e-14fails1e-14. Its norm3.7854e-5is almost entirely two B rows at cell261; their state ULPs are1.22070e-4and4.88281e-4. Further corrections to those coordinates can round away. This observation is not a proof of the minimum attainable residual.',
        method='At a near-converged failed linear iterate, select up to24 non-B material columns coupled to the dominant B rows. Keep columns whose full operator response to one state ULP is below one quarter of the original residual allowance. Solve their restricted full-row least-squares correction, round into the actual stored array, and re-evaluate the original full vector and physical residuals. At most4polishing passes; require actual residual reduction. Fall back to the existing12balanced corrections only if necessary.',
        invariants='No equation, state precision, physical normalization, gate, time grid, source, front or horizon is changed. A polished iterate is a linear proposal only; the original native nonlinear stage, constitutive, ledger, packet and time gates still decide acceptance.',
        reuse='The196last accepted state is bit-identical to195. Reuse113accepted substeps and the latest failed linear iterate, only after exact RHS/branch-guess identity. Six coarse and16fine actual substeps remain.',
        gates=dict(linear=1e-14,physical_linear=1e-13,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,time=.02),
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='196took554.38s, mostly repeated Krylov attempts after the physical channels converged. Stored-array inspection took2.47s. Up to24local columns cost unmeasured additional operator calls per polishing pass; allow the same generous30/45minute production caps. No new prefix replay or separate expensive system reproduction.',
        stop='Original physical/accuracy failure,8Newton/12linear limits,4polish passes per refinement, nonfinite candidate or30/45minute cap. No automatic new grid or period.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'restart-check.json',dict(classification='Counterexample candidate',passed=True,
        all196accepted_state_keys_exactly_match195=True,prior_exact_restart_verified=read(OLD/'restart-check.json')['passed'],new_physical_steps=0))
    import sympy as sp
    A,E,x,d,b=sp.symbols('A E x d b',commutative=False)
    assert sp.expand(b-A*(x+E*d)-((b-A*x)-A*E*d))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='A correction on selected columns preserves the original full linear equation; no convergence or attainable-precision guarantee.'))


def initialize():
    restore=FunctionType(prior.restore_substep.__code__,dict(prior.restore_substep.__globals__,OUT=OUT))
    init=FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT,restore_substep=restore));init()


def polish(m,op,rhs,sol):
    scale=max(np.linalg.norm(rhs),LD('1e-290'));goal=LD('1e-14')*scale
    cv=inspect.getclosurevars(op._CustomLinearOperator__matvec_impl).nonlocals
    blocks=cv['blocks'];dim=len(rhs)//2;records=[]
    for attempt in range(4):
        residual=rhs-op.matvec(sol);norm=np.linalg.norm(residual)
        if norm<goal:break
        if norm>goal*10000:break
        target=[]
        for stage,row in enumerate(residual.reshape(2,dim)):
            _,gas=m.unpack(row)
            target += [(stage,int(k)) for k in np.flatnonzero(abs(gas[:,2])>goal/8)]
        if not target:break
        columns=sorted({col for stage,cell in target for col,value in blocks[stage][cell] if value and col%4!=2})
        candidates=[]
        for col in columns:
            index=(col//(4*m.n))*dim+m.size+col%(4*m.n)
            quantum=abs(np.spacing(sol[index])) if sol[index] else LD(0)
            score=max((abs(LD(str(dict(blocks[stage][cell]).get(col,0))))*quantum for stage,cell in target),default=LD(0))
            if score<goal/4:candidates.append((float(score),index,quantum))
        candidates.sort();selected=[];images=[];norms=[]
        for _,index,quantum in candidates[:24]:
            unit=np.zeros_like(sol);unit[index]=1;image=op.matvec(unit);cn=np.linalg.norm(image)
            if not np.isfinite(cn) or cn==0 or cn*quantum>=goal/4:continue
            selected.append(index);norms.append(cn);images.append(sparse.csc_matrix(np.asarray(image/cn,float)[:,None]))
        if not images:break
        matrix=sparse.hstack(images,format='csc');rows=np.unique(matrix.nonzero()[0])
        delta,_,rank,_=np.linalg.lstsq(matrix[rows].toarray(),np.asarray(residual[rows],float),rcond=1e-14)
        candidate=sol.copy();candidate[selected]+=np.asarray(delta,LD)/np.asarray(norms,LD)
        actual=rhs-op.matvec(candidate);newnorm=np.linalg.norm(actual)
        moments=joint.physical_norm(m,actual)/joint.scales(m,rhs,candidate)
        accepted=bool(newnorm<norm and max(moments)<1e-13 and np.all(np.isfinite(candidate)))
        records.append(dict(before=float(norm/scale),after=float(newnorm/scale),rank=int(rank),columns=selected,
            rows=len(rows),physical_relative=moments.astype(float).tolist(),accepted=accepted))
        if not accepted:break
        sol=candidate
    path=OUT/'polish.json';old=read(path) if path.exists() else []
    write(path,old+[dict(classification='Counterexample candidate',attempts=records)])
    return sol


def factory(log,n):
    source=inspect.getsource(prior.balanced_factory)
    source=source.replace('    def solve(m,op,P,rhs,guess):','    first_system=[True]\n    def solve(m,op,P,rhs,guess):')
    source=source.replace('if n==64 and not log:', 'if n==64 and first_system[0]:\n            first_system[0]=False')
    source=source.replace("        ns=dict(joint.solve.__globals__,gmres=krylov);exec(compile(source,__file__,'exec'),ns)",'''        source=source.replace('sol,info=gmres(op,rhs,x0=guess,**options);sol=sol.astype(LD)', 'sol=guess.copy();info=0')
        mark='        if relative<1e-14 and max(moments)<1e-13:'
        extra="        if relative>=1e-14 and relative<1e-10 and max(moments)<1e-14:\\n            sol=polish(m,op,rhs,sol);residual=rhs-op.matvec(sol)\\n            relative=float(np.linalg.norm(residual)/max(np.linalg.norm(rhs),1e-290));moments=physical_norm(m,residual)/scales(m,rhs,sol)\\n"
        assert source.count(mark)==1;source=source.replace(mark,extra+mark)
        source=source.replace('delta,_=gmres(', 'delta,code=gmres(').replace('sol+=delta', 'sol+=delta\\n        if k==0:info=code')
        ns=dict(joint.solve.__globals__,gmres=krylov,polish=polish);exec(compile(source,__file__,'exec'),ns)''')
    ns=dict(vars(prior),OUT=OUT,OLD=OLD,polish=polish);exec(compile(source,__file__,'exec'),ns)
    return ns['balanced_factory'](log,n)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        if action in ['coarse','fine']:
            assert read(OUT/'restart-check.json')['passed'];stable.initialize=initialize;n=64 if action=='coarse' else 128
            old=SimpleNamespace(right_solver=lambda log:factory(log,n))
            evolve=FunctionType(before.evolve.__code__,dict(before.evolve.__globals__,OUT=OUT,old=old));evolve(n)
        elif action=='audit':
            audit=FunctionType(before.audit.__code__,dict(before.audit.__globals__,OUT=OUT));audit()
            result=read(OUT/'result.json');result.update(maximum_inner_solves=12,maximum_Newton_proposals=8);write(OUT/'result.json',result)
        else:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
