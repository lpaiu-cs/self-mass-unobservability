"""Run the exterior connection with the public OdeSolution array boundary.

Keep the registered physical equations and failed geometry-check receipt.
SciPy accepts scalar/one-dimensional query times, so restore caller shapes.
"""
from pathlib import Path
import json, resource, sys, time
import numpy as np
import couple_native_exterior_incident as run


def optical(self,r):
    r=np.asarray(r)
    return run.C*self.D*self.vacuum.sol((r.ravel()-self.r0)/(run.C*self.D))[0].reshape(r.shape)


def primitive(self,r):
    r=np.asarray(r)
    return self.vacuum.sol((r.ravel()-self.r0)/(run.C*self.D))[1].reshape(r.shape)/self.r0


if __name__=='__main__':
    action=sys.argv[1];assert action in ['check_repair','pilot','production']
    run.Interaction.optical=optical;run.Interaction.primitive=primitive
    dest=run.OUT/f'{action}-receipt.json';assert not dest.exists()
    cap=run.CAPS['check']-run.read(run.OUT/'check-receipt.json')['seconds'] if action=='check_repair' else run.CAPS[action]
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.previous.incident.native.deadline(cap)
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        for p,h in run.read(run.OUT/'plan.json')['bindings'].items():assert run.sha(p)==h,p
        plan=run.OUT/'array-boundary-repair-plan.json'
        if action=='check_repair':
            assert not plan.exists();failed=run.read(run.OUT/'check-receipt.json')
            assert failed['error'].startswith('ValueError')
            run.write(plan,dict(classification='Counterexample candidate',source_sha256=run.sha(__file__),
                change='Flatten only public SciPy OdeSolution query arrays and reshape its results to the original caller dimensions. The physical equations, solver tolerances and integration gates are unchanged. The original geometry check stopped before any scattering evaluation.',
                original_producer_sha256=run.sha(run.__file__),failure_sha256=run.sha(run.OUT/'check-receipt.json'),
                remaining_check_seconds=cap))
        else:assert run.read(plan)['source_sha256']==run.sha(__file__)
        getattr(run,'check' if action=='check_repair' else action)()
    except Exception as exc:error=repr(exc);raise
    finally:run.write(dest,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=run.sha(__file__)))
