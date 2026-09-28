"""Replay the failed native partition with one fixed positive spectral measure."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
import def_photon_native_measure as model

OUT=model.OUT;ex=model.ex;h=model.h


def prepare():
    assert not (OUT/'replay-plan.json').exists()
    paths=[Path(__file__),Path(model.__file__),OUT/'plan.json',OUT/'result.json',OUT/'surf1000T0015-bank.npz',
        OUT.parent/'def-photon-partition-audit/result.json',OUT.parent/'def-photon-display-fixed/requests.json']
    ex.write(OUT/'replay-plan.json',dict(classification='Counterexample candidate',checkpoint='48d04e7',
        claim='Use the already frozen native-node opacity field to predict the exact window whose provider 1/3/999-group means disagreed. No new fit or physical query.',
        method='Integrate the same piecewise-constant source-node opacities over arbitrary new boundaries through one cumulative Planck/Rosseland primitive. Compare its 1/3/999 partitions and the prior native one-group mean.',
        gates=dict(native_gray_relative=.001,partition_relative=1e-11,complex_energy_balance=1e-12),
        budget=dict(hard_seconds=30,new_physical_queries=0,new_native_EOS_calls=0,new_stellar_steps=0),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths}))


def cumulative(points,edges,k,T,column):
    """Integral of blackbody weight times k in one fixed spectral step field."""
    weights=model.loss.weights(edges,T)[:,column]
    prefix=np.r_[0,np.cumsum(weights*k)]
    index=np.searchsorted(edges,points,side='right')-1
    index=np.clip(index,0,len(k)-1)
    lo=edges[index]/T;hi=points/T
    x,w=leggauss(24);u=(lo[:,None]+hi[:,None])/2+(hi-lo)[:,None]*x/2
    p=15/np.pi**4*u**3*np.exp(-u)/(-np.expm1(-u))
    if column:p*=u/(4*(-np.expm1(-u)))
    return prefix[index]+(p@w)*(hi-lo)/2*k[index]


def run():
    assert not (OUT/'replay-result.json').exists();start=time.monotonic()
    plan=json.loads((OUT/'replay-plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    bank=np.load(OUT/'surf1000T0015-bank.npz');T=.0015
    request=next(r for r in json.loads((OUT.parent/'def-photon-display-fixed/requests.json').read_text()) if r['name']=='phLine4T0')
    a,b=map(float,[request['fields']['egplow'],request['fields']['egphigh']]);checks=[]
    for count in [1,3,999]:
        boundaries=np.geomspace(a,b,count+1);w=model.loss.weights(boundaries,T)
        P=np.diff(cumulative(boundaries,bank['edges_keV'],bank['absorption'],T,0))
        H=np.diff(cumulative(boundaries,bank['edges_keV'],1/bank['total'],T,1))
        assert np.all(P>0) and np.all(H>0)
        checks.append(dict(groups=count,Rosseland_Planck=[float(w[:,1].sum()/H.sum()),float(P.sum()/w[:,0].sum())]))
    values=np.array([r['Rosseland_Planck'] for r in checks]);partition=float(np.max(abs(values/values[0]-1)))
    old=json.loads((OUT.parent/'def-photon-partition-audit/result.json').read_text())
    comparison=abs(values[0]/old['checks'][0]['Rosseland_Planck']-1)
    # Independent small dense controls for the energy-norm operator bound.
    # The proof is Loewner order plus the resolvent identity, not these samples.
    eps=5e-5;rng=np.random.default_rng(6001);max_ratio=0.
    for _ in range(12):
        capacities=np.exp(rng.uniform(-2,2,8));rates=np.exp(rng.uniform(-3,3,7))
        columns=np.zeros((8,7));columns[0,:]=-np.sqrt(capacities[1:]/capacities[0]);columns[1:,:]=np.eye(7)
        L=(columns*rates)@columns.T
        altered=(columns*(rates*(1+rng.uniform(-eps,eps,7))))@columns.T
        for omega in np.geomspace(.001,1000,17):
            z=1j*omega;I=np.eye(8)
            error=np.linalg.norm(z*(np.linalg.inv(L+z*I)-np.linalg.inv(altered+z*I)),2)
            bound=eps/(2*np.sqrt(1-eps))
            assert error<=bound+1e-10
            max_ratio=max(max_ratio,float(error/bound))
    result=dict(classification='Counterexample candidate',passed=bool(partition<1e-11 and comparison.max()<.001),
        checks=checks,partition_relative=partition,native_one_group_relative=comparison.tolist(),
        original_native_fine_partition_failure_preserved=True,
        interpretation='One fixed source-node field predicts the old one-group mean and retains additivity under all three requested partitions. The original provider fine-group output remains rejected; this does not establish its internal integration algorithm or a continuum line-error bound.',
        seconds=time.monotonic()-start)
    ex.write(OUT/'replay-result.json',result)
    ex.write(OUT/'coupled-rounding-theorem.json',dict(classification='Proven',passed=True,
        assumptions='Fixed positive matter/photon heat capacities and independent relative errors <=eps<1 in positive absorptive rates. These are printed-rate errors, not physical/model, spectral-resolution or EOS errors.',
        generator='In energy-norm coordinates z=(sqrt(C_m)*dT,E_i/sqrt(C_i)), the generator is -L=-sum_i gamma_i*b_i*b_i^T, b_i=(-sqrt(C_i/C_m),e_i). Its conservation nullspace is unchanged by rate perturbations.',
        proof='The perturbed Lp obeys (1-eps)L<=Lp<=(1+eps)L and deltaL=L^(1/2) F L^(1/2), ||F||<=eps. The resolvent identity, L<=Lp/(1-eps), and sup_lambda sqrt(lambda)/|lambda+iw|=1/sqrt(2|w|) give ||iw[(Lp+iw)^(-1)-(L+iw)^(-1)]||<=eps/[2sqrt(1-eps)] for w!=0. The conserved nullspace cancels exactly.',
        eps=eps,uniform_energy_norm_resolvent_bound=eps/(2*np.sqrt(1-eps)),
        finite_dense_controls=204,maximum_control_to_bound_ratio=max_ratio))
    print('REPLAY',result,flush=True);assert result['seconds']<30


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
