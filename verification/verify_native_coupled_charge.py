"""Independent closed-form ray and normalization checks, no new evolution."""
from types import FunctionType,SimpleNamespace
import json
import time
import numpy as np
import sympy as sp
import def_native_coupled_charge as task


def main():
    out=task.OUT;assert not (out/'exterior-audit.json').exists();start=time.monotonic()
    task.write(out/'exterior-audit-plan.json',dict(classification='Counterexample candidate',seconds=10,native_calls=0,new_fluid_steps=0,
        check='Use flat vacuum with constant boundary occupation, where the continuous-angle retarded energy has an independent cubic closed form. Verify the exact mass-normalization algebra and bound a declared sign-coherent inner baryon debit.',
        gate=1e-9,source_sha256=task.sha(__file__)))
    R=7e9;T=task.prior.END;L=1.5e33
    class Flat:
        def __init__(self,m):pass
        def __call__(self,x):return np.ones_like(x),np.zeros_like(x),np.ones_like(x),np.ones_like(x),np.ones_like(x)
    def metric(r):return np.zeros_like(r),np.ones_like(r),np.ones_like(r),np.ones_like(r),np.zeros_like(r)
    model=SimpleNamespace(m=SimpleNamespace(RJ=R,R=R,rf=np.array([R]),bg=SimpleNamespace(metric=metric)),bulk=SimpleNamespace(edges_mu=np.linspace(-1,1,9)))
    green=SimpleNamespace(Geometry=Flat,polynomial=task.green.polynomial)
    fn=FunctionType(task.rays.__code__,dict(task.rays.__globals__,green=green),argdefs=task.rays.__defaults__)
    t=np.linspace(0,T,17);u=np.linspace(0,T,129);energy,row=fn(model,t,np.full_like(t,L),u)
    k=R/task.C;exact=L*(u*u/(2*k)-u*u*u/(6*k*k));error=float(np.max(abs(energy-exact))/max(exact));assert error<1e-9
    M,K,Q,dm=sp.symbols('M K Q dm',nonzero=True)
    stable=(-Q/M+(-K/M)*dm/M)/(1-dm/M)
    assert sp.simplify(-(K+Q)/(M-dm)+K/M-stable)==0
    d=np.load(out/'source-896-128.npz');w=d['weight'];B=np.asarray(d['baryon_g'].sum(1),float);Bmax=float(max(abs(B)));duration=float(d['t'][-1])
    # This is a conditional envelope for ONE sign-coherent mass debit with
    # total variation <=Bmax, supported in the stored causal radial region.
    # It is not a bound on arbitrary unresolved positive/negative flows.
    debit_bound=task.G*float(d['cx'])/(2*task.C*float(d['M_cm']))*float(max(abs(w)))*duration*Bmax
    direct=np.load(out/'wave-896-128.npz')['normalized_direct'];debit=np.load(out/'wave-896-128.npz')['components'][3]
    result=dict(classification='Counterexample candidate',passed=True,flat_geodesic_arrival_relative=error,
        exact_mass_normalization_algebra=True,conditional_baryon_debit_normalized_bound=debit_bound,
        conditional_baryon_premise='A single-sign compensating mass history of magnitude no greater than the represented net baryon imbalance; spatial support inside the stored causal region and no additional hidden counterflows.',
        bound_over_direct_peak=debit_bound/max(abs(direct)),modeled_debit_peak=float(max(abs(debit))),
        bound_type='Finite stored coefficient maximum and declared source envelope; continuous physical-domain uniformity is not certified.',
        seconds=time.monotonic()-start,full_GR_source_bound=False,final_charge_solved=False)
    task.write(out/'exterior-audit.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
