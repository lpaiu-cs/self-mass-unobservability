"""Roundoff-aware correction of the same linear system, with original gates.

Counterexample candidate: nearby representable hard variables plus the existing
fine-variable correction. No physical output, stage equation or gate is fitted.
"""
from pathlib import Path
from types import FunctionType
from itertools import combinations
import inspect,json,os,resource,sys,time
import numpy as np
from scipy import sparse
import repair_full_momentum_collision as prior

OUT=Path('native-rounding217-work');OLD=prior.OLD
LD,joint,read,write,sha=prior.LD,prior.joint,prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,check=900,refine=900)


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[]
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material']:(OUT/folder).mkdir(parents=True)
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        os.link(src,OUT/src.relative_to(OLD));files.append(src)
    files += [OLD/n for n in ['failed-linear-64.npz','last-accepted-64.npz','failure-64.json']]
    files += [prior.OUT/n for n in ['result.json','repair-result.json','symbolic.json']]+[Path('native-broad-polish215-work/proposal.npz')]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='see_git_checkpoint',
        claim='Determine whether nearby representable gas coordinates solve the original retained117th system; if so, use that original-equation solver in actual continuation.',
        evidence='215fine-column correction leaves3.301e-13.216full collision evaluation changes nothing and photon directions do not add an eligible column at this inactive atmospheric row. The fine correction cannot move large coarse-ULP S/B coordinates.',
        method='Use the SAME fixed linear operator. Form original native B/S column images once. Split fine columns by the existing norm-times-ULP threshold. For up to8hard columns, search pairs within64ULPs of their least-squares centre (centre bounded1024ULPs), eliminate fine corrections by SVD, then evaluate the24best proposals with the full original operator and physical gates. At most4improving passes. These are numerical linear-solver coordinates, not physical parameter trials.',
        gates=dict(linear=1e-14,physical_linear=1e-13,actual_stage=1e-12,physical_stage=1e-13),budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='Constructor/system reconstruction about25s; at most96full residual checks plus small dense pair searches. Allow15minutes; no GMRES or physical step. Preserve failures and never accept a merely physical-moment pass.',
        stop='Exact reconstruction mismatch, original gates or bounded search/cap. No automatic grid/period change or prolonged repeated solve if no useful correction.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    write(OUT/'symbolic.json',read(prior.OUT/'symbolic.json'))


def neighbours(m,op,rhs,sol,log,triplets=False):
    scale=max(np.linalg.norm(rhs),LD('1e-290'));goal=LD('1e-14')*scale;dim=len(rhs)//2
    blocks=inspect.getclosurevars(op._CustomLinearOperator__matvec_impl).nonlocals['blocks_by_component']
    for attempt in range(4):
        residual=rhs-op.matvec(sol);norm=np.linalg.norm(residual)
        if norm<goal:break
        target=[]
        for stage,row in enumerate(residual.reshape(2,dim)):
            _,gas=m.unpack(row)
            target += [(stage,int(k),c) for c in [2,3] for k in np.flatnonzero(abs(gas[:,c])>goal/8)]
        columns=sorted({col for stage,cell,c in target for col,value in blocks[c][stage][cell] if value})
        data=[];images=[]
        for col in columns:
            index=(col//(4*m.n))*dim+m.size+col%(4*m.n);quantum=abs(np.spacing(sol[index])) if sol[index] else LD(0)
            unit=np.zeros_like(sol);unit[index]=1;image=op.matvec(unit);cn=np.linalg.norm(image)
            if cn==0 or not np.isfinite(cn):continue
            data.append((index,quantum,cn));images.append(sparse.csc_matrix(np.asarray(image/cn,float)[:,None]))
        if not images:break
        matrix=sparse.hstack(images,format='csc');rows=np.unique(matrix.nonzero()[0]);M=matrix[rows].toarray();r=np.asarray(residual[rows],float)
        fine=[j for j,(_,q,cn) in enumerate(data) if q*cn<goal/4];hard=[j for j in range(len(data)) if j not in fine]
        if not fine or not hard:break
        F=M[:,fine];u,s,vt=np.linalg.svd(F,full_matrices=False);keep=s>s[0]*1e-14;Q=u[:,keep];pinv=(vt[keep].T/s[keep])@Q.T
        pr=r-Q@(Q.T@r);H=M[:,hard]*np.array([float(data[j][1]*data[j][2]) for j in hard]);PH=H-Q@(Q.T@H)
        sizes=np.linalg.norm(PH,axis=0);order=[int(i) for i in np.argsort(sizes) if sizes[i]>float(goal)/128][:8]
        guesses={};count=3 if triplets else 2;extent=16 if triplets else 64
        grid=np.stack(np.meshgrid(*([np.arange(-extent,extent+1)]*count),indexing='ij'),axis=-1).reshape(-1,count)
        for selection in combinations(order,count):
            C=PH[:,selection];centre=np.linalg.lstsq(C,pr,rcond=1e-14)[0]
            if np.max(abs(centre))>1024:continue
            points=grid+np.rint(centre).astype(int);G=C.T@C;d=C.T@pr
            scores=np.einsum('ni,ij,nj->n',points,G,points)-2*(points@d)+pr@pr
            for k in np.argsort(scores)[:5]:
                shifts=np.zeros(len(hard),dtype=int);shifts[list(selection)]=points[k]
                key=tuple(shifts);guesses[key]=min(guesses.get(key,float('inf')),float(scores[k]))
        best=sol;best_norm=norm;tested=[]
        for key,score in sorted(guesses.items(),key=lambda item:item[1])[:24]:
            candidate=sol.copy();shifts=np.array(key)
            for j,k in enumerate(hard):
                if shifts[j]:index,q,_=data[k];candidate[index]+=LD(int(shifts[j]))*q
            delta=pinv@(r-H@shifts)
            for j,d in zip(fine,delta):index,_,cn=data[j];candidate[index]+=LD(d)/cn
            actual=rhs-op.matvec(candidate);newnorm=np.linalg.norm(actual);moments=joint.physical_norm(m,actual)/joint.scales(m,rhs,candidate)
            passed=bool(newnorm<best_norm and max(moments)<1e-13 and np.all(np.isfinite(candidate)))
            tested.append(dict(indices=[data[hard[j]][0] for j in np.flatnonzero(shifts)],ULPs=shifts[shifts!=0].tolist(),relative=float(newnorm/scale),physical=moments.astype(float).tolist(),improved=passed))
            if passed:best,best_norm=candidate,newnorm
            if best_norm<goal:break
        log.append(dict(pass_index=attempt,before=float(norm/scale),after=float(best_norm/scale),fine_columns=len(fine),hard_columns=len(hard),searched_hard_columns=len(order),projected_hard_norms=sizes.tolist(),tested=tested))
        if best_norm>=norm:break
        sol=best
    return sol


def check(refine=False):
    m,op,rhs,z,p,cs,ss,t,h,v=FunctionType(prior.system.__code__,dict(prior.system.__globals__,OUT=OUT))()
    start=OUT/'proposal.npz' if refine else Path('native-broad-polish215-work/proposal.npz')
    sol=np.load(start)['solution'];logs=[];sol=neighbours(m,op,rhs,sol,logs,triplets=refine);r=rhs-op.matvec(sol);scale=np.linalg.norm(rhs)
    physical=joint.physical_norm(m,r)/joint.scales(m,rhs,sol);rates=[];m.precise_values={};owner=prior.prior.prior
    for j,row in enumerate(sol.reshape(2,-1)):
        x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True);native=owner.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
    defect=(sol.reshape(2,-1)-v-h*(joint.A@np.array(rates))).ravel();defect=owner.precise_defect(m,t,h,v,sol,defect)
    result=dict(classification='Counterexample candidate',linear_relative=float(np.linalg.norm(r)/scale),physical=physical.astype(float).tolist(),actual_stage_relative=float(np.linalg.norm(defect)/scale),
        linear_passed=bool(np.linalg.norm(r)/scale<1e-14 and max(physical)<1e-13),actual_stage_passed=bool(np.linalg.norm(defect)/scale<1e-12 and max(joint.physical_norm(m,defect)/joint.scales(m,rhs,sol))<1e-13),new_physical_steps=0,new_Krylov_iterations=0,final_charge_conclusion='unadjudicated')
    prefix='refined-' if refine else ''
    write(OUT/f'{prefix}search.json',logs);np.savez_compressed(OUT/f'{prefix}proposal.npz',solution=sol,residual=r,actual_defect=defect);write(OUT/f'{prefix}result.json',result);print(json.dumps(result),flush=True)


def refine():
    write(OUT/'refinement-plan.json',dict(classification='Conjectural',
        claim='Apply three-coordinate representable corrections to the same failed117th linear system; passing only admits a Newton proposal, not a physical stage.',
        reason='The paired search reduced3.301e-13 to1.668e-14 in68.66s; it cannot compensate three coarse lattice directions together. This is now close to the unchanged1e-14gate.',
        method='At most8hard columns, triples within16ULPs of the least-squares centre bounded1024ULPs; deduplicate shared proposals; full original operator checks the24best, at most4improving passes. Reuse pair proposal; no prefix replay.',
        caps=CAPS,virtual_GiB=6,CPU_threads=1,forecast='Measured pair search68.66s. Allow900s for small projected triple enumeration; at most96full checks, no new physical steps or Krylov iterations.',
        gates=read(OUT/'plan.json')['gates'],bindings={str(p):sha(p) for p in [Path(__file__),OUT/'pair-producer.py',OUT/'proposal.npz',OUT/'result.json']},final_charge_conclusion='unadjudicated'))
    check(True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                bound=OUT/'pair-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(bound)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
