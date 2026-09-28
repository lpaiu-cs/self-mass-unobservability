"""Counterexample candidate: exact momentum rows in the saved failed system."""
from pathlib import Path
from types import FunctionType
from decimal import Decimal,localcontext
import inspect,json,os,resource,sys,time
import numpy as np
import continue_native_flux_precision as prior

OUT=Path('native-momentum209-work');OLD=prior.OUT
read,write,sha=prior.read,prior.write,prior.sha
owner,joint,LD=prior.owner,prior.joint,prior.LD
CAPS=dict(prepare=180,check=600)


def prepare():
    assert not OUT.exists();OUT.mkdir();assert read(OLD/'controller-status.json')['state']=='failed'
    assert 'Four-moment linear residual' in read(OLD/'coarse-receipt.json')['error']
    files=list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    for folder in ['sweep-0/photons','sweep-0/material','sweep-1/photons','sweep-1/material']:(OUT/folder).mkdir(parents=True,exist_ok=True)
    reused={}
    for p in files:
        dst=OUT/p.relative_to(OLD);os.link(p,dst);reused[str(dst)]=sha(p)
    files += [OLD/n for n in ['failed-linear-64.npz','last-accepted-64.npz','coarse-receipt.json','failure-64.json','linear-64.json','polish.json']]
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None) and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Resolve the actual116th-stage linear residual in saved202data; distinguish momentum-row arithmetic error from unrepresentable or insufficiently converged stored solution before another physical continuation.',
        evidence='202coarse terminated after4181.74s before its7200s cap.115steps accepted; next linear residual1.79440e-14fails1e-14. Saved norm is dominated by momentum cell261:1.45688e-5residual at1.52588e-5state ULP; original norm allowance8.24210e-6. B residuals are about1e-8absolute. Existing polishing targets only B rows.',
        method='Reconstruct the last saved branch Jacobian, RHS and original80digit B operator bit-identically. Assemble S native rows with exact binary h/A/J coefficients and subtract their collision contribution at80digits. Compare40/80digits and the original operator before accepting any arithmetic change. Evaluate actual original nonlinear defect separately.',
        scope='One stored system, no Krylov iteration or accepted physical step. Keep original1e-14linear/1e-13physical and1e-12actual stage gates. A stored linear pass is not physical acceptance.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=6,forecast='Same model and two Jacobian setup previously20..30s; sparse momentum arithmetic is small. Allow10minutes without replaying115accepted steps.',
        stop='RHS/residual/branch mismatch, precision disagreement or wall cap. Preserve202failure and all physical gates.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)},reused=reused,final_charge_conclusion='unadjudicated',full_goal_complete=False))
    import sympy as s
    x,y,j,k,h,a,b,q,r=s.symbols('x y j k h a b q r')
    assert s.expand(x-h*(a*(j*x+q)+b*(k*y+r))-((1-h*a*j)*x-h*b*k*y-h*(a*q+b*r)))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Native momentum preassembly plus original collision is the identical linear equation, not a physical error bound.'))


def momentum_operator(m,op,cs,h,Js,precision=80):
    def exact(v):
        p,q=v.as_integer_ratio();return Decimal(p)/Decimal(q)
    with localcontext() as ctx:
        ctx.prec=precision;dh=exact(h);da=[[exact(v) for v in row] for row in joint.A];blocks=[]
        for i in range(2):
            rows=[]
            for k in range(3,4*m.n,4):
                co={i*4*m.n+k:Decimal(1)}
                for j,J in enumerate(Js):
                    for p in range(J.indptr[k],J.indptr[k+1]):
                        col=j*4*m.n+int(J.indices[p]);co[col]=co.get(col,Decimal(0))-dh*da[i][j]*exact(J.data[p])
                rows.append(tuple(co.items()))
            blocks.append(rows)
    def apply(value):
        result=np.asarray(op.matvec(value),LD).copy();pairs=[m.unpack(row) for row in value.reshape(2,-1)]
        collisions=np.array([m.collision(c,x,g)[1][:,3] for c,(x,g) in zip(cs,pairs)])
        with localcontext() as ctx:
            ctx.prec=precision;gas=[exact(v) for _,g in pairs for v in g.ravel()]
            collision=[[exact(v) for v in row] for row in collisions]
            for i,row in enumerate(result.reshape(2,-1)):
                g=m.unpack(row)[1]
                for cell,co in enumerate(blocks[i]):g[cell,3]=LD(str(sum(v*gas[k] for k,v in co)-dh*sum(da[i][j]*collision[j][cell] for j in range(2))))
        return result
    return joint.LinearOperator(op.shape,apply,dtype=float)


def initialize():
    FunctionType(prior.initialize.__code__,dict(prior.initialize.__globals__,OUT=OUT))(False)


def check():
    initialize();m=owner.Model(64);z=dict(np.load(OLD/'failed-linear-64.npz'));p=dict(np.load(OLD/'last-accepted-64.npz'))
    t,h=p['next_time'][()],p['next_step'][()];v=m.pack(p['x'],p['g']);dim=len(v)
    guides=[m.unpack(row)[1] for row in z['guess'].reshape(2,-1)]
    maps=[m.jacobian(t+c*h,g) for c,g in zip(joint.C,guides)];Js=[q[0] for q in maps]
    cs=[m.local(t+c*h) for c in joint.C];ss=[m.source(t+c*h) for c in joint.C]
    def L(j,value):
        x,g=m.unpack(value);ph,q,*_=m.collision(cs[j],x,g)
        return m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph,q+(Js[j]@g.ravel()).reshape(m.n,4))
    def mat(value):
        x=value.reshape(2,dim);return (x-h*(joint.A@np.array([L(j,row) for j,row in enumerate(x)]))).ravel()
    raw=joint.LinearOperator((2*dim,)*2,mat,dtype=float);op=prior.stable.stable_operator(m,raw,[],True)
    affine=[b-(J@g.ravel()).reshape(m.n,4) for (J,b),g in zip(maps,guides)]
    src=np.array([m.pack(s[0]/(m.scale*prior.joint.AMP)+c['q'],m.gas(c['q'],c['qb'],c['qe'])+a) for s,c,a in zip(ss,cs,affine)])
    rhs=(np.tile(v,(2,1))+h*(joint.A@src)).ravel();rhs=prior.precise_rhs(m,t,h,v,maps,guides,rhs)
    assert np.array_equal(rhs,z['rhs']),'Saved RHS must reproduce bit-identically'
    original=rhs-op.matvec(z['solution']);assert np.array_equal(original,z['residual']),'Saved residual must reproduce bit-identically'
    residuals=[]
    for digits in [40,80]:residuals.append(rhs-momentum_operator(m,op,cs,h,Js,digits).matvec(z['solution']))
    norm=np.linalg.norm(rhs);fixed=residuals[-1];moments=joint.physical_norm(m,fixed)/joint.scales(m,rhs,z['solution'])
    rates=[];m.precise_values={}
    for j,row in enumerate(z['solution'].reshape(2,dim)):
        x,g=m.unpack(row);ph,q,*_=m.collision(cs[j],x,g,True);native=prior.precise_native(m,t+joint.C[j]*h,g,m.native(t+joint.C[j]*h,g,details=True))
        rates.append(m.pack((m.A@x.reshape(m.n*m.q,m.nf)).reshape(x.shape)+ph+ss[j][0]/(m.scale*joint.AMP),q+native[0]))
    defect=(z['solution'].reshape(2,dim)-v-h*(joint.A@np.array(rates))).ravel();defect=prior.precise_defect(m,t,h,v,z['solution'],defect)
    row=dict(classification='Counterexample candidate',exact_saved_rhs_and_residual=True,old_linear_relative=float(np.linalg.norm(original)/norm),precise_momentum_linear_relative=float(np.linalg.norm(fixed)/norm),linear_physical=moments.astype(float).tolist(),
        momentum_arithmetic_effect=float(np.linalg.norm(fixed-original)/norm),precision_difference=float(np.linalg.norm(residuals[0]-fixed)/norm),
        actual_stage_relative=float(np.linalg.norm(defect)/norm),actual_stage_physical=(joint.physical_norm(m,defect)/joint.scales(m,rhs,z['solution'])).astype(float).tolist(),
        linear_gates_passed=bool(np.linalg.norm(fixed)/norm<1e-14 and max(moments)<1e-13),new_Krylov_iterations=0,new_physical_steps=0,physical_stage_accepted=False,final_charge_conclusion='unadjudicated')
    np.savez_compressed(OUT/'residual-comparison.npz',original=original,precise=fixed,actual=defect)
    write(OUT/'result.json',row);print(json.dumps(row),flush=True);assert row['precision_difference']<1e-25,row


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
