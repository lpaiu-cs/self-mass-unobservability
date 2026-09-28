from pathlib import Path
import json,numpy as np,sympy as sp
import def_gr_matrix_propagation as m
o=m.OUT;w=o/'weeks';b=m.prior
p=b.Problem(b.prior.OUT/'fine-bank.npz');L,_=m.source_map(p);errors={}
fields=b.FIELDS
cases=[]
for n in [512,1024,2048]:
 d=dict(np.load(w/f'fine-{n}.npz'));speed=d['velocity_native'];scalar=d['scalar_native']
 expected=np.column_stack([np.sqrt((speed*speed)@p.weights),np.sqrt((scalar*scalar)@p.weights)]+[np.sqrt((speed[:,mask]**2)@p.weights[mask]/p.weights[mask].sum()) for mask in p.masks])
 assert np.allclose(expected,d['readouts'],rtol=2e-13,atol=0)
 v=np.interp(p.native,p.bg.grid,p.speed*d['velocity'][:,0])
 f=np.interp(p.native,p.bg.grid,d['response'][:,2]-p.bg.grid*p.bg.nodes['v']*d['response'][:,0])
 e=dict(velocity=float(np.max(abs(v-speed[-1]))/max(abs(speed[-1]).max(),1e-100)),scalar=float(np.max(abs(f-scalar[-1]))/max(abs(scalar[-1]).max(),1e-100)))
 errors[n]=e
 assert max(e.values())<2e-7,e
 cases.append(d['readouts'])
result=json.loads((w/'result.json').read_text());assert result['passed'] is False and result['contrasts']=={}
for j,field in enumerate(fields):
 a,c,d=[v[:,j] for v in cases];norm=abs(d).max();e1=max(abs(a-c))/norm;e2=max(abs(c-d))/norm
 row=result['comparisons'][field]
 assert np.allclose([e1,e2,np.log2(e1/e2)],[row['previous'],row['last'],row['order']],rtol=2e-13,atol=1e-14)
# The exact coordinate identity only uses zero inertia in auxiliary columns.
assert abs(p.D[:,1::4]).sum()==0 and abs(p.D[:,3::4]).sum()==0
model=m.Propagator(p,L);finite=m.Finite(model)
for seed in [17,29]:
 q=np.random.default_rng(seed).normal(size=finite.size)
 full=finite.extend(q)
 # Eulerian subtraction has a known coordinate condition factor; compare
 # against its floating-point bound, not an unscaled two-ulp assertion.
 condition=max(1.,np.max(abs(model.advection*full.reshape(-1,4)[:,0]/model.units[2])))
 assert np.max(abs(finite.restrict(full)-q))<8*np.finfo(float).eps*condition
qz,qp,adv,uz,up=sp.symbols('qz qp adv uz up',nonzero=True)
assert sp.simplify(((qp*up+adv*qz*uz)-adv*qz*uz)/up-qp)==0
z=sp.symbols('z');M,K,F=sp.symbols('M K F')
assert sp.expand((1-z*z*M)*(F/K)+(F/K)*z*z*M-F/K)==0
out=dict(classification='Counterexample candidate',saved_replay_passed=True,endpoint_reconstruction_errors=errors,failed_readout_verdict_reproduced=True,finite_coordinate_identity_passed=True,new_GR_resolvents=0,new_EOS_calls=0,original_failure_resolved=False)
m.write(o/'verification.json',out);print(json.dumps(out),flush=True)
