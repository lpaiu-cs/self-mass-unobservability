"""Counterexample candidate: localize the rejected residual, without new steps."""
from pathlib import Path
from types import FunctionType
import inspect,json,resource,sys,time
import numpy as np
import right_precondition_full_interval as old

OUT=Path('native-linear-floor192-work');read,write,sha=old.read,old.write,old.sha
LD=old.LD;CAP=75


def arithmetic():
    from decimal import Decimal,localcontext
    from scipy import sparse
    assert (OUT/'diagnosis.json').exists() and not (OUT/'arithmetic-plan.json').exists()
    write(OUT/'arithmetic-plan.json',dict(classification='Conjectural',
        claim='Resolve whether the dominant baryon residual is inaccurate evaluation or an inaccurate stored iterate. Reconstruct the same operator without Krylov work and evaluate all baryon rows using80decimal digits.',
        decision='Compare original residual with exact-binary-coefficient reference and with the algebraically identical preassembled material block. Accept no physical step or tolerance change.',
        budget_seconds=45,CPU_threads=1,virtual_GiB=4,new_physical_steps=0,new_Krylov_iterations=0,
        forecast='The same model/operator setup previously consumed about28s; sparse gas-row arithmetic is small compared with the1.3million-dimensional system.45s cap.',
        bindings={str(p):sha(p) for p in [Path(__file__),OUT/'rejected-system.npz',OUT/'diagnosis.json']}))
    def factory(log):
        def evaluate(m,op,P,rhs,guess):
            z=dict(np.load(OUT/'rejected-system.npz'));assert np.array_equal(rhs,z['rhs'])
            sol=z['solution'];cv=inspect.getclosurevars(op.matvec.__self__._CustomLinearOperator__matvec_impl).nonlocals
            L=cv['L'];h=cv['h'];A=old.joint.A;Js=inspect.getclosurevars(L).nonlocals['Js']
            gas=np.array([m.unpack(row)[1].ravel() for row in sol.reshape(2,-1)])
            raw=np.array([m.unpack(row)[1].ravel() for row in rhs.reshape(2,-1)])
            direct=np.array([m.unpack(row)[1].ravel() for row in z['residual'].reshape(2,-1)])
            def exact(v):
                p,q=v.as_integer_ratio();return Decimal(p)/Decimal(q)
            refs=[];assembled=[];blocks=[]
            with localcontext() as c:
                c.prec=80;dh=exact(h);da=[[exact(a) for a in row] for row in A]
                dg=[[exact(v) for v in row] for row in gas]
                for i in range(2):
                    vals=[];ri=[];ci=[];rr=[]
                    for k in range(2,4*m.n,4):
                        coefficients={i*4*m.n+k:Decimal(1)}
                        for j,J in enumerate(Js):
                            for p in range(J.indptr[k],J.indptr[k+1]):
                                col=j*4*m.n+int(J.indices[p]);coefficients[col]=coefficients.get(col,Decimal(0))-dh*da[i][j]*exact(J.data[p])
                        value=exact(raw[i,k])-sum(v*dg[col//(4*m.n)][col%(4*m.n)] for col,v in coefficients.items())
                        rr.append(LD(str(value)))
                        for col,v in coefficients.items():ri.append(k//4);ci.append(col);vals.append(LD(str(v)))
                    block=sparse.csr_matrix((np.array(vals,LD),(ri,ci)),shape=(m.n,8*m.n))
                    assembled.append(raw[i,2::4]-block@gas.ravel());refs.append(rr);blocks.append(block)
            exact_res=np.array(refs,LD);assembled=np.array(assembled,LD);den=np.linalg.norm(rhs)
            result=dict(classification='Counterexample candidate',
                original_baryon_relative=float(np.linalg.norm(direct[:,2::4])/den),
                exact_binary_baryon_relative=float(np.linalg.norm(exact_res)/den),
                preassembled_baryon_relative=float(np.linalg.norm(assembled)/den),
                original_evaluation_error=float(np.linalg.norm(direct[:,2::4]-exact_res)/den),
                preassembled_evaluation_error=float(np.linalg.norm(assembled-exact_res)/den),
                gate=1e-14,physical_solution_accepted=False,final_charge_conclusion='unadjudicated')
            np.savez_compressed(OUT/'baryon-reference.npz',exact=exact_res,assembled=assembled)
            for i,J in enumerate(Js):sparse.save_npz(OUT/f'gas-J-{i}.npz',J)
            for i,M in enumerate(blocks):sparse.save_npz(OUT/f'baryon-block-{i}.npz',M)
            write(OUT/'arithmetic.json',result);print(json.dumps(result),flush=True)
            raise StopIteration('Arithmetic complete; no new solve')
        return evaluate
    check=FunctionType(old.check.__code__,dict(old.check.__globals__,OUT=OUT,right_solver=factory))
    try:check()
    except StopIteration:pass


def run():
    assert not OUT.exists();OUT.mkdir()
    previous=read(old.OUT/'check-receipt.json');assert previous['error'] is not None
    files=[Path(__file__),Path(old.__file__),old.OUT/'system-reconstruction.json',old.OUT/'check-receipt.json',old.OUT/'right-calls.json',old.OLD/'failed-linear-system.npz']
    write(OUT/'plan.json',dict(classification='Conjectural',
        claim='Locate the remaining true-vector residual before any further solver change. Repeat the identical saved-system solve once solely because the prior failure did not serialize its accumulated solution.',
        decision='Distinguish incomplete Krylov convergence from floating-point operator/solution cancellation by saving the rejected iterate, per-channel norms and a two-part evaluation of the same linear operator. Neither evaluation is accepted as a physical solution.',
        budget_seconds=CAP,CPU_threads=1,virtual_GiB=4,new_physical_steps=0,
        forecast='Identical reconstruction and solver took49.24s; allow75s including serialization. One diagnostic, no extra iteration/grid/tolerance change.',
        stop='One reproduction,75s cap or reconstructed-system mismatch. No physical continuation admitted.',
        bindings={str(p):sha(p) for p in files},final_charge_conclusion='unadjudicated'))
    def factory(log):
        solve=old.right_solver(log)
        def capture(m,op,P,rhs,guess):
            try:return solve(m,op,P,rhs,guess)
            except AssertionError as exc:
                tb=exc.__traceback__
                while tb and tb.tb_frame.f_code!=old.joint.solve.__code__:tb=tb.tb_next
                assert tb is not None
                sol=tb.tb_frame.f_locals['sol'].copy();residual=rhs-op.matvec(sol)
                high=sol.astype(float).astype(LD);low=sol-high
                split=rhs-op.matvec(high)-op.matvec(low)
                np.savez_compressed(OUT/'rejected-system.npz',rhs=rhs,guess=guess,solution=sol,residual=residual,split_residual=split)
                blocks=[]
                for name,value in [('rhs',rhs),('solution',sol),('residual',residual),('split_residual',split)]:
                    pairs=[m.unpack(v) for v in value.reshape(2,-1)]
                    gas=np.array([g for x,g in pairs]);photons=np.array([x for x,g in pairs])
                    blocks.append(dict(name=name,norm=float(np.linalg.norm(value)),photon_norm=float(np.linalg.norm(photons)),
                        gas_norm=np.sqrt(np.sum(gas**2,axis=(0,1))).astype(float).tolist(),
                        largest_index=int(np.argmax(abs(value))),largest_value=float(value[np.argmax(abs(value))])))
                report=dict(classification='Counterexample candidate',original_error=repr(exc),blocks=blocks,
                    split_vector_relative=float(np.linalg.norm(split)/np.linalg.norm(rhs)),
                    direct_vector_relative=float(np.linalg.norm(residual)/np.linalg.norm(rhs)),
                    operator_split_disagreement=float(np.linalg.norm(split-residual)/np.linalg.norm(rhs)),
                    same_system_replayed=True,new_physical_steps=0,physical_solution_accepted=False,final_charge_conclusion='unadjudicated')
                write(OUT/'diagnosis.json',report);print(json.dumps(report),flush=True)
                raise
        return capture
    check=FunctionType(old.check.__code__,dict(old.check.__globals__,OUT=OUT,right_solver=factory))
    try:check()
    except AssertionError:
        assert (OUT/'diagnosis.json').exists()
    else:raise AssertionError('Prior failure did not reproduce')


if __name__=='__main__':
    action=sys.argv[1] if len(sys.argv)>1 else 'run'
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));old.joint.previous.original.inf.incident.native.deadline(45 if action=='arithmetic' else CAP)
    start=time.monotonic();error=None
    try:globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/('arithmetic-receipt.json' if action=='arithmetic' else 'receipt.json'),dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
