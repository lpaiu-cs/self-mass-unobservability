"""Counterexample candidate: evaluate cancelling baryon rows before rounding.

Only the extended residual evaluation changes. The exact binary coefficients,
Krylov budget, physical stage equations, time grid and gates are unchanged.
"""
from pathlib import Path
from types import FunctionType
from decimal import Decimal,localcontext
import inspect,json,os,resource,sys,time
import numpy as np
import right_precondition_full_interval as old

OUT=Path('native-stable-residual193-work');read,write,sha=old.read,old.write,old.sha
owner,joint,LD=old.owner,old.joint,old.LD
CAPS=dict(prepare=20,check=90,consistent_check=90)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    for s in [0,1]:
        for part in ['photons','material']:(OUT/f'sweep-{s}/{part}').mkdir(parents=True)
    files=list((old.OUT/'sweep-0').rglob('*.npz'))+[old.OUT/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]
    reused={}
    for src in files:
        dst=OUT/src.relative_to(old.OUT);os.link(src,dst);reused[str(dst)]=sha(src)
    files += [Path(__file__),Path(old.__file__),Path('native-linear-floor192-work/arithmetic.json'),Path('native-linear-floor192-work/rejected-system.npz')]
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Remove demonstrated cancellation in the actual baryon linear residual, then test the same saved failing system at the unchanged1e-14/1e-13 gates.',
        evidence='The stored rejected solution has an80digit exact-coefficient baryon residual2.20654e-14, while the existing evaluation reports5.00219e-14. Its evaluation error2.93089e-14 alone exceeds the gate. Preassembly in long double does not repair it.',
        method='Keep float64 Krylov proposals and original4refinements. For extended-precision residual calls evaluate each baryon row as one80decimal-digit sum with the exact original binary h,A,J coefficients and actual stored solution, then round the small result once. Other rows use the original operator. No collision source enters the B row.',
        decision='Only a matching original system and all original linear gates admit a separately budgeted remaining-interval continuation. Report the old arithmetic residual too; do not replace physical-stage acceptance with this check.',
        forecast='Original reconstruction24.59s and right-Krylov21.24s; allow90s including selective decimal arithmetic. One system check only.',
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,new_physical_steps=0,
        stop='Any original gate,4refinement limit,source mismatch or90s cap. No automatic iteration,grid,tolerance or horizon increase.',
        bindings={str(p):sha(p) for p in files},reused=reused,final_charge_conclusion='unadjudicated'))
    import sympy as sp
    x,y,j,k,h,a,b=sp.symbols('x y j k h a b')
    assert sp.expand(x-h*(a*j*x+b*k*y)-((1-h*a*j)*x-h*b*k*y))==0
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Baryon row preassembly is the same linear polynomial; this is not a physical error bound.'))


def initialize():
    old.prior.prior.prior.OUT=OUT;old.prior.prior.prior.initialize()


def stable_operator(m,op,audit,consistent=False):
    cv=inspect.getclosurevars(op._CustomLinearOperator__matvec_impl).nonlocals
    h=cv['h'];Js=inspect.getclosurevars(cv['L']).nonlocals['Js'];A=joint.A
    def exact(v):
        p,q=v.as_integer_ratio();return Decimal(p)/Decimal(q)
    with localcontext() as c:
        c.prec=80;dh=exact(h);da=[[exact(a) for a in row] for row in A];blocks=[]
        for i in range(2):
            rows=[]
            for k in range(2,4*m.n,4):
                coefficients={i*4*m.n+k:Decimal(1)}
                for j,J in enumerate(Js):
                    for p in range(J.indptr[k],J.indptr[k+1]):
                        col=j*4*m.n+int(J.indices[p]);coefficients[col]=coefficients.get(col,Decimal(0))-dh*da[i][j]*exact(J.data[p])
                rows.append(tuple(coefficients.items()))
            blocks.append(rows)
    def apply(value):
        result=op.matvec(value)
        if not consistent and value.dtype!=np.dtype(LD):return result
        start=time.monotonic();gas=np.array([m.unpack(row)[1].ravel() for row in value.reshape(2,-1)])
        result=np.asarray(result,LD).copy();values=result.reshape(2,-1)
        with localcontext() as c:
            c.prec=80;dg=[exact(v) for v in gas.ravel()]
            for i,rows in enumerate(blocks):
                _,out=m.unpack(values[i])
                for cell,coefficients in enumerate(rows):out[cell,2]=LD(str(sum(v*dg[k] for k,v in coefficients)))
        audit.append(time.monotonic()-start)
        return result
    return joint.LinearOperator(op.shape,apply,dtype=float)


def check(consistent=False):
    logs=[];arithmetic=[];target=OUT/'consistent' if consistent else OUT
    if consistent:target.mkdir()
    def factory(unused):
        solve=old.right_solver(logs)
        def exact_solve(m,op,P,rhs,guess):
            stable=stable_operator(m,op,arithmetic,consistent)
            try:sol=solve(m,stable,P,rhs,guess)
            except AssertionError as exc:
                tb=exc.__traceback__
                while tb and tb.tb_frame.f_code!=joint.solve.__code__:tb=tb.tb_next
                if tb:
                    sol=tb.tb_frame.f_locals['sol'];np.savez_compressed(target/'rejected-solution.npz',solution=sol,residual=rhs-stable.matvec(sol))
                raise
            stable_res=rhs-stable.matvec(sol);original_res=rhs-op.matvec(sol)
            result=dict(classification='Counterexample candidate',passed=True,
                vector_relative=float(np.linalg.norm(stable_res)/np.linalg.norm(rhs)),
                physical_relative=(joint.physical_norm(m,stable_res)/joint.scales(m,rhs,sol)).astype(float).tolist(),
                original_arithmetic_vector_relative=float(np.linalg.norm(original_res)/np.linalg.norm(rhs)),
                selective_decimal_seconds=sum(arithmetic),exact_system_reconstructed=True,
                new_physical_steps=0,final_charge_conclusion='unadjudicated')
            assert result['vector_relative']<1e-14 and max(result['physical_relative'])<1e-13
            np.savez_compressed(target/'accepted-linear-solution.npz',solution=sol,residual=stable_res)
            write(target/'stable-check-result.json',result);print(json.dumps(result),flush=True)
            raise StopIteration('Linear check passed; no physical step accepted')
        return exact_solve
    probe=FunctionType(old.check.__code__,dict(old.check.__globals__,OUT=target,initialize=initialize,right_solver=factory))
    try:probe()
    except StopIteration:pass
    finally:write(target/'solver-cost.json',dict(classification='Counterexample candidate',calls=logs,decimal_calls=len(arithmetic),decimal_seconds=sum(arithmetic)))


def consistent_check():
    assert read(OUT/'check-receipt.json')['error'] is not None
    write(OUT/'consistent-plan.json',dict(classification='Conjectural',
        failure='Selective evaluation exposed a mismatch between the inaccurate float-input Krylov operator and the accurate long-double residual operator: one correction had true relative residual39.4787. The original four-refinement trial remained rejected.',
        repair='Evaluate the same baryon polynomial accurately on every matvec input, including float64 Krylov vectors. Do not switch the arithmetic operator by input dtype. Preserve all equations, gates and iteration limits.',
        decision='One saved-system consistency repair only; no physical continuation unless all original gates pass.',
        budget_seconds=90,CPU_threads=1,virtual_GiB=4,new_physical_steps=0,
        forecast='Selective pass took73.55s including400inner iterations;13exact evaluations cost0.10s, so consistent evaluation adds about3.1s at that count.90s cap.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'selective-producer.py',OUT/'check-receipt.json',OUT/'solver-cost.json']},
        final_charge_conclusion='unadjudicated'))
    check(True)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();error=None
    try:
        if action!='prepare':
            for p,h in dict(read(OUT/'plan.json')['bindings'],**read(OUT/'plan.json')['reused']).items():
                actual=OUT/'selective-producer.py' if Path(p).resolve()==Path(__file__).resolve() else p
                assert sha(actual)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
