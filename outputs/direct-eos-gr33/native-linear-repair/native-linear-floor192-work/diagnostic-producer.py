"""Counterexample candidate: localize the rejected residual, without new steps."""
from pathlib import Path
from types import FunctionType
import inspect,json,resource,sys,time
import numpy as np
import right_precondition_full_interval as old

OUT=Path('native-linear-floor192-work');read,write,sha=old.read,old.write,old.sha
LD=old.LD;CAP=75


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
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));old.joint.previous.original.inf.incident.native.deadline(CAP)
    start=time.monotonic();error=None
    try:run()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(OUT/'receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
