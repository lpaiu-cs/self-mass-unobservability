"""Evaluate the full native-plus-collision momentum row before rounding.

Counterexample candidate: inspect the saved117th system, keep all original gates.
"""
from pathlib import Path
from types import FunctionType
import inspect,json,os,resource,sys,time
import numpy as np
import continue_precise_momentum as prior

OUT=Path('native-full-momentum216-work');OLD=prior.OUT
LD=prior.LD;joint=prior.joint;precision=prior.prior.precision;mp=precision.mp;hp=precision.scalar
read,write,sha=prior.read,prior.write,prior.sha
CAPS=dict(prepare=180,check=900,repair=900)


def prepare():
    assert not OUT.exists();OUT.mkdir();files=[]
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material']:(OUT/folder).mkdir(parents=True)
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        os.link(src,OUT/src.relative_to(OLD));files.append(src)
    files += [OLD/n for n in ['failed-linear-64.npz','last-accepted-64.npz','failure-64.json','coarse-receipt.json']]
    files += [Path('native-broad-polish215-work')/n for n in ['proposal.npz','result.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='a79af854d',
        claim='Determine whether the remaining momentum residual is actually the floating collision evaluation rather than an unavoidable state lattice. Repair that arithmetic before changing state representation or restarting a long solve.',
        method='Reconstruct210RHS and residual exactly. Reuse215polished proposal. Expand the original linear photon collision, escape and local four-gas momentum coefficients at40/80digits at the observed limiting atmospheric cell261. Combine them with the same native80digit rows and exact binary h/A coefficients before rounding. Independently compare direct high-precision collision contractions; no fitted state or physical output.',
        gates=dict(linear=1e-14,physical_linear=1e-13,precision=1e-25,actual_stage=1e-12),budgets=CAPS,CPU_threads=1,virtual_GiB=6,
        forecast='Same saved-system construction44.4s in215 including polishing. Sparse one-cell collision contraction cost unmeasured; allow15minutes,0GMRES/physical steps. Stop on coefficient mismatch, precision disagreement or cap. Do not widen acceptance or rerun the old long system without a useful correction.',
        decision='If full-row arithmetic explains the defect, apply the consistent same operator to actual Krylov/refinement and independently test true nonlinear acceptance. Otherwise preserve this rejected hypothesis and address representation.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as s
    w,l,S,x,B,g,e,E=s.symbols('w l S x B g e E')
    assert s.expand(w*(-l*x+S*x+B*g)+e*x+E*g-((-w*l+w*S+e)*x+(w*B+E)*g))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Collapse the same linear collision and escape momentum contraction; no physical or EOS error bound.'))


def rows(m,c,cells,digits):
    result={};width=m.q*m.nf
    with mp.workdps(digits):
        for cell in cells:
            assert cell>=m.nb,'Inventory response is zero only in the atmospheric rows used here'
            weight=[hp(e)*hp(mu) for e,mu in zip(m.Eweight[cell].ravel(),np.broadcast_to(m.mu[:,None],(m.q,m.nf)).ravel())]
            den=hp(m.a[cell])*hp(m.su[cell]);ph={cell*width+k:-w*hp(loss)+hp(esc) for k,(w,loss,esc) in enumerate(zip(weight,c['loss'][cell].ravel(),c['esc'][2,cell].ravel()))}
            S=c['S'].tocsr()
            for k,w in enumerate(weight):
                i=cell*width+k
                for j in range(S.indptr[i],S.indptr[i+1]):
                    col=int(S.indices[j]);ph[col]=ph.get(col,mp.mpf(0))+w*hp(S.data[j])
            gas=[]
            for B,Be in [(c['B'],c['Be']),(c['DBS'],c['DBSe'])]:
                for component in range(2):gas.append(-(sum(w*hp(v) for w,v in zip(weight,B[cell,...,component].ravel()))+hp(Be[2,cell,component]))/den)
            result[cell]=(tuple((k,-v/den) for k,v in ph.items() if v),gas)
    return result


def full_operator(m,op,cs,h,cells,digits):
    blocks_by_component=inspect.getclosurevars(op._CustomLinearOperator__matvec_impl).nonlocals['blocks_by_component']
    native=blocks_by_component[3]
    with mp.workdps(digits):
        collision=[rows(m,c,cells,digits) for c in cs];dh=hp(h);da=[[hp(v) for v in row] for row in joint.A]
        coefficients={}
        for i in range(2):
            for cell in cells:
                gas={k:mp.mpf(str(v)) for k,v in native[i][cell]};ph=[]
                for j in range(2):
                    pp,gg=collision[j][cell];factor=-dh*da[i][j]
                    ph.append(tuple((k,factor*v) for k,v in pp))
                    for c,v in enumerate(gg):
                        k=j*4*m.n+4*cell+c;gas[k]=gas.get(k,mp.mpf(0))+factor*v
                coefficients[i,cell]=(tuple(gas.items()),ph)
    def apply(value):
        assert len(blocks_by_component)==2
        result=np.asarray(op.matvec(value),LD).copy();pairs=[m.unpack(row) for row in value.reshape(2,-1)]
        with mp.workdps(digits):
            gas=[hp(v) for _,g in pairs for v in g.ravel()]
            xs=[x.ravel() for x,_ in pairs]
            for (i,cell),(gg,pp) in coefficients.items():
                val=sum(v*gas[k] for k,v in gg)
                for j,terms in enumerate(pp):val+=sum(v*hp(xs[j][k]) for k,v in terms)
                m.unpack(result.reshape(2,-1)[i])[1][cell,3]=LD(str(val))
        return result
    return joint.LinearOperator(op.shape,apply,dtype=float),collision


def system():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)
    m=prior.owner.Model(64);z=dict(np.load(OLD/'failed-linear-64.npz'));p=dict(np.load(OLD/'last-accepted-64.npz'));t,h=p['next_time'][()],p['next_step'][()]
    v=m.pack(p['x'],p['g']);dim=len(v);guides=[m.unpack(row)[1] for row in z['guess'].reshape(2,-1)]
    maps=[m.jacobian(t+c*h,g) for c,g in zip(joint.C,guides)];Js=[r[0] for r in maps];cs=[m.local(t+c*h) for c in joint.C];ss=[m.source(t+c*h) for c in joint.C]
    def L(j,value):
        x,g=m.unpack(value);ph,q,*_=m.collision(cs[j],x,g)
        return m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph,q+(Js[j]@g.ravel()).reshape(m.n,4))
    def mat(value):
        x=value.reshape(2,dim);return (x-h*(joint.A@np.array([L(j,row) for j,row in enumerate(x)]))).ravel()
    raw=joint.LinearOperator((2*dim,)*2,mat,dtype=float);op=prior.stable_operator(m,raw,[],True)
    affine=[b-(J@g.ravel()).reshape(m.n,4) for (J,b),g in zip(maps,guides)]
    src=np.array([m.pack(s[0]/(m.scale*joint.AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+a) for s,c,a in zip(ss,cs,affine)])
    rhs=(np.tile(v,(2,1))+h*(joint.A@src)).ravel();rhs=prior.prior.precise_rhs(m,t,h,v,maps,guides,rhs)
    assert np.array_equal(rhs,z['rhs']) and np.array_equal(rhs-op.matvec(z['solution']),z['residual'])
    return m,op,rhs,z,p,cs,ss,t,h,v


def check():
    m,op,rhs,z,p,cs,ss,t,h,v=system()
    sol=np.load('native-broad-polish215-work/proposal.npz')['solution'];old=rhs-op.matvec(sol);norm=np.linalg.norm(rhs);results=[];vectors=[]
    for digits in [40,80]:
        start=time.monotonic();fixed,collision=full_operator(m,op,cs,h,[261],digits);r=rhs-fixed.matvec(sol);vectors.append(r)
        results.append(dict(digits=digits,relative=float(np.linalg.norm(r)/norm),seconds=time.monotonic()-start))
    # Direct sparse high-precision photon contraction, without collapsed rows.
    direct=[]
    with mp.workdps(80):
        for j,(x,g) in enumerate(m.unpack(row) for row in sol.reshape(2,-1)):
            c=cs[j];cell=261;width=m.q*m.nf;S=c['S'].tocsr();xx=x.ravel();gg=list(map(hp,g[cell]));value=mp.mpf(0)
            for k in range(width):
                q,f=divmod(k,m.nf);i=cell*width+k
                ph=-hp(c['loss'][cell,q,f])*hp(xx[i])+sum(hp(S.data[a])*hp(xx[S.indices[a]]) for a in range(S.indptr[i],S.indptr[i+1]))
                ph+=sum(hp(c['B'][cell,q,f,a])*gg[a]+hp(c['DBS'][cell,q,f,a])*gg[2+a] for a in range(2))
                value+=ph*hp(m.Eweight[cell,q,f])*hp(m.mu[q])+hp(c['esc'][2,cell,q,f])*hp(xx[i])
            value+=sum(hp(c['Be'][2,cell,a])*gg[a]+hp(c['DBSe'][2,cell,a])*gg[2+a] for a in range(2));value/=-hp(m.a[cell])*hp(m.su[cell])
            pp,gc=collision[j][cell];collapsed=sum(v*hp(xx[k]) for k,v in pp)+sum(v*w for v,w in zip(gc,gg))
            error=abs(value-collapsed)/max(abs(value),mp.mpf('1e-290'));assert error<mp.mpf('1e-60');direct.append(float(error))
    change=float(np.linalg.norm(vectors[0]-vectors[1])/norm);assert change<1e-25
    result=dict(classification='Counterexample candidate',exact_saved_RHS_and_residual=True,old_linear_relative=float(np.linalg.norm(old)/norm),full_momentum_evaluations=results,precision_change=change,direct_collision_control=direct,full_arithmetic_change=float(np.linalg.norm(vectors[-1]-old)/norm),
        linear_passed=bool(np.linalg.norm(vectors[-1])/norm<1e-14),new_Krylov_iterations=0,new_physical_steps=0,final_charge_conclusion='unadjudicated',full_goal_complete=False)
    np.savez_compressed(OUT/'comparison.npz',old_residual=old,full_residual=vectors[-1]);write(OUT/'result.json',result);print(json.dumps(result),flush=True)


def photon_columns(m,sol,target,collision,h,goal):
    cells={cell for _,cell,component in target if component==3};dim=len(sol)//2;groups={}
    for j,compiled in enumerate(collision):
        for cell in cells:
            if cell not in compiled:continue
            for k,value in compiled[cell][0]:
                index=j*dim+k;quantum=abs(np.spacing(sol[index])) if sol[index] else LD(0)
                effect=abs(LD(str(value))*h)*max(abs(joint.A[:,j]));score=effect*quantum
                if effect==0 or score>=goal/4:continue
                group=j,cell,(k%(m.q*m.nf))//m.nf
                if group not in groups or effect>groups[group][0]:groups[group]=(effect,float(score),index,quantum)
    return [(score,index,quantum) for _,score,index,quantum in groups.values()]


def coupled_polish(m,op,rhs,sol,collision,h):
    FunctionType(prior.polish_function.__code__,dict(prior.polish_function.__globals__,OUT=OUT))()
    source=(OUT/'expanded-BS-polish.py').read_text()
    changes=[('if norm>goal*10000:break','if norm>goal*1000000:break'),
        ('candidates.sort();selected=[];images=[];norms=[]','candidates += photon_columns(m,sol,target,collision,h,goal)\n        candidates.sort();selected=[];images=[];norms=[]')]
    for a,b in changes:assert source.count(a)==1,a;source=source.replace(a,b)
    ns=dict(prior.legacy.polish.__globals__,OUT=OUT,photon_columns=photon_columns,collision=collision,h=h);exec(compile(source,__file__,'exec'),ns)
    (OUT/'expanded-coupled-polish.py').write_text(source);return ns['polish'](m,op,rhs,sol)


def repair():
    old=read(OUT/'result.json');assert old['full_arithmetic_change']==0 and old['precision_change']==0
    write(OUT/'repair-plan.json',dict(classification='Conjectural',
        rejected_hypothesis='40/80digit full momentum collision assembly and independent direct contraction agree exactly and leave the same3.301e-13residual. Do not describe collision arithmetic as its cause.',
        change='The restricted corrector only searches native gas columns. Add the original photon collision columns that act on the limiting S rows: strongest eligible frequency in each angular bin and stage, filtered by the unchanged coefficient-times-ULP and full-column-times-ULP tests. Keep24totalcolumns,4passes, full original vector and physical improvement tests. No observable or port fitting.',
        claim='Solve the retained linear system at its original gates using its actual coupled variables, then apply the same basis in real continuation if it passes. Any actual nonlinear failure remains a Newton proposal, never an accepted physical state.',
        cap_seconds=900,CPU_threads=1,virtual_GiB=6,new_physical_steps=0,new_Krylov_iterations=0,
        forecast='Original reconstruction and4pass gas polishing44.4s;16extra candidate photon columns add at mostabout32operator evaluations across4passes under24totalcolumns. Allow15minutes; no physical prefix replay.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'collision-producer.py',OUT/'result.json',Path('native-broad-polish215-work/proposal.npz')]},final_charge_conclusion='unadjudicated'))
    m,op,rhs,z,p,cs,ss,t,h,v=system();sol=np.load('native-broad-polish215-work/proposal.npz')['solution']
    fixed,collision=full_operator(m,op,cs,h,[261],80);sol=coupled_polish(m,fixed,rhs,sol,collision,h)
    residual=rhs-fixed.matvec(sol);scale=np.linalg.norm(rhs);physical=joint.physical_norm(m,residual)/joint.scales(m,rhs,sol)
    rates=[];m.precise_values={};dim=len(rhs)//2
    for j,row in enumerate(sol.reshape(2,dim)):
        x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True);native=prior.prior.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
    defect=(sol.reshape(2,dim)-v-h*(joint.A@np.array(rates))).ravel();defect=prior.prior.precise_defect(m,t,h,v,sol,defect)
    result=dict(classification='Counterexample candidate',linear_relative=float(np.linalg.norm(residual)/scale),physical=physical.astype(float).tolist(),actual_stage_relative=float(np.linalg.norm(defect)/scale),
        linear_passed=bool(np.linalg.norm(residual)/scale<1e-14 and max(physical)<1e-13),actual_stage_passed=bool(np.linalg.norm(defect)/scale<1e-12 and max(joint.physical_norm(m,defect)/joint.scales(m,rhs,sol))<1e-13),new_physical_steps=0,new_Krylov_iterations=0,final_charge_conclusion='unadjudicated')
    np.savez_compressed(OUT/'coupled-proposal.npz',solution=sol,residual=residual,actual_defect=defect);write(OUT/'repair-result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists();start=time.monotonic();error=None
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():
                bound=OUT/'collision-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(bound)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
