"""Check patched geometry, unchanged inputs and the finite inversion side."""
from pathlib import Path
import json
import signal
import time
import numpy as np
from scipy.sparse import diags
from scipy.linalg import cholesky_banded
import def_gr_interface_patch as task

go=task.go;OUT=task.OUT;task.verify_plan();task.install();assert not (OUT/'operator-check.json').exists()
signal.alarm(60);start=time.monotonic()
task.write(OUT/'operator-check-plan.json',dict(classification='Counterexample candidate',
    claim='Verify the fixed source-preserving patch, native endpoint readouts and numerical inversion side on the actual accepted finite matrix.',
    budget=dict(hard_seconds=60,new_evolutions=0,new_EOS_calls=0),
    limitation='Floating shifted Cholesky is not an interval SPD proof or continuum stability certificate.',
    bindings={str(p):go.task.digest(p) for p in [Path(__file__),Path(task.__file__),Path(go.__file__),OUT/'p4-2048.npz',OUT/'p4-result.json']}))
p=go.Problem();m=p.model;d=np.load(OUT/'p4-2048.npz');old=np.load(task.PRIOR/'beta1024/p4-2048.npz')
assert abs(m.original.radiation.geometry.tc-0.23080495568542375)<1e-15
for key,value in [('grid',m.grid),('cells',m.cells),('indices',m.indices),('native_radius',m.original.native),('weights',m.original.weights),('masks',m.original.masks)]:
    assert np.array_equal(d[key],value),key
for key in ['native_radius','weights','masks','heat_energy','heat_flux']:assert np.array_equal(d[key],old[key]),key
assert np.all(np.isin(old['cells'],m.cells)) and m.constant_gradient_after==0
cells,ids=task.patch_cells();added=m.cells[~np.isin(m.cells,cells)]
assert len(added)==3*len(ids) and np.all((added>cells[ids.min()])&(added<cells[ids.max()+1]))
source_edge=m.heat.edges[m.heat.face_ids.max()];separation=float(abs(m.cells-source_edge).min())
assert separation<1e-14
native=np.array([p.speed*(m.nativeV[0]@d['qt']),m.nativeV[1]@d['q']]);saved=np.array([d['native_velocity'][-1],d['native_scalar'][-1]])
norm=lambda x:go.transfer.norms(m,x)
endpoint=norm(native-saved)/norm(saved);assert np.max(endpoint)<1e-9
A=(p.K+144*p.M).tocsc();D=diags(1/np.sqrt(A.diagonal()));S=(D@A@D).tocsc();coo=S.tocoo()
width=int(max(abs(coo.row-coo.col)));n=m.size;band=np.zeros((width+1,n))
for j in range(width+1):band[j,:n-j]=S.diagonal(-j)
C=cholesky_banded(band,lower=True);L=diags([C[j,:n-j] for j in range(width+1)],-np.arange(width+1),shape=(n,n)).tocsc()
error=float(np.max(abs((L@L.T-S).data),initial=0));assert error<1e-12 and C[0].min()>0
result=dict(classification='Counterexample candidate',passed=True,added_trial_nodes=len(added),dofs=n,
    physical_horizon_s=float(m.original.radiation.geometry.tc),
    source_support_face_to_trial_node_distance_R=separation,unchanged_native_readouts_and_heat=True,
    endpoint_readout_relative=dict(zip(go.task.FIELDS,endpoint.astype(float))),
    shifted_cholesky_passed=True,shift=144,minimum_factor_diagonal=float(C[0].min()),scaled_factor_residual=error,
    seconds=time.monotonic()-start,new_evolutions=0,
    scope='Actual finite operator and endpoint consistency; no continuum or full physical stellar closure.')
task.write(OUT/'operator-check.json',result);signal.alarm(0);print(json.dumps(result,indent=2))
