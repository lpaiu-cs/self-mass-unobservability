"""A fixed nested spatial contrast preserving the fine weak source exactly."""
from pathlib import Path
from types import SimpleNamespace
import json
import signal
import time
import numpy as np
from scipy.sparse import coo_matrix
import def_gr_energy_modes as modes

OUT=modes.fem.OUT/'nested';OUT.mkdir(exist_ok=True)
assert not (OUT/'result.json').exists()
modes.write(OUT/'plan.json',dict(classification='Counterexample candidate',
    claim='Separate spatial trial-space error from moving quadrature and heat-source debit in the failed every-other-node contrast.',
    method='The same coarse node subset, Galerkin restriction Kc=P.T*Kf*P, Mc=P.T*Mf*P, Fc=P.T*Lf, and exact fine heat lift in physical readouts. No new nodes, no altered energy input, no renewed EOS. Restrict u=q+H rather than separately reinterpolating H.',
    gates=dict(spatial_relative=.02,projection_relative=.02,projection_order=1.5),
    budget=dict(models=1,maximum_basis=512,hard_seconds=60,CPU_threads=1,new_EOS_calls=0),
    forecast='One native assembly4s plus one512-vector projection and three series:8-15s based on22.60s for four prior models. Hard60s.',
    bindings={str(p):modes.digest(p) for p in [Path(__file__),Path(modes.__file__),Path(modes.fem.__file__),modes.OUT/'result.json']}))
signal.alarm(60);start=time.monotonic();fine=modes.fem.Model(modes.evolution.BANK/'fine-bank.npz')
si=fine.surface_index;nodes=np.unique(np.r_[np.arange(0,si+1,2),si,np.arange(si+1,len(fine.grid),2),len(fine.grid)-1])
grid=fine.grid[nodes];cs=int(np.where(nodes==si)[0][0]);nf=cs+1;n=len(grid);size=n+nf-1
indices=np.full((n,2),-1,int);indices[:nf,0]=2*np.arange(nf);indices[:nf,1]=2*np.arange(nf)+1
indices[nf:,1]=2*nf+np.arange(n-nf);indices[-1,1]=-1
row=[];col=[];data=[]
for field in [0,1]:
    valid=fine.indices[:,field]>=0;r=fine.grid[valid];ix=np.clip(np.searchsorted(grid,r,side='right')-1,0,len(grid)-2)
    f=(r-grid[ix])/(grid[ix+1]-grid[ix])
    for delta,values in [(0,1-f),(1,f)]:
        target=indices[ix+delta,field];keep=target>=0
        row.extend(fine.indices[valid,field][keep]);col.extend(target[keep]);data.extend(values[keep])
P=coo_matrix((data,(row,col)),shape=(fine.size,size)).tocsc()
coarse=SimpleNamespace(M=(P.T@fine.M@P).tocsc(),K=(P.T@fine.K@P).tocsc(),load=(P.T@fine.load).tocsc(),
    heat=fine.heat,original=fine.original,size=size)
projection=modes.Projection(coarse);_,Q,orth=projection.basis(512);physical=P@Q
K,skew=modes.energy_matrix(fine,physical)
original_OUT=modes.OUT;modes.OUT=OUT
cases={str(n):modes.series(fine,K[:n,:n],physical[:,:n],'coarse-'+str(n)) for n in [128,256,512]}
modes.OUT=original_OUT
reference=json.loads((modes.OUT/'fine-512.json').read_text());comparisons={}
for field in modes.evolution.FIELDS:
    a,b,c=[np.array([v[field] for v in cases[str(n)]['history']]) for n in [128,256,512]]
    f=np.array([v[field] for v in reference['history']]);norm=max(abs(f).max(),1e-100)
    d1=np.max(abs(a-b))/norm;d2=np.max(abs(b-c))/norm
    comparisons[field]=dict(previous=float(d1),last=float(d2),order=float(np.log2(d1/d2)),spatial=float(np.max(abs(c-f))/norm))
passed=all(v['last']<.02 and v['order']>1.5 and v['spatial']<.02 for v in comparisons.values())
result=dict(classification='Counterexample candidate',passed=passed,comparisons=comparisons,orthogonality=orth,
    stiffness_skew=skew,linear_residual=projection.error,coarse_dofs=size,seconds=time.monotonic()-start,
    original_failure_resolved=False,full_dynamic_charge_solved=False)
modes.write(OUT/'result.json',result);signal.alarm(0);print('NESTED',json.dumps(result),flush=True)
