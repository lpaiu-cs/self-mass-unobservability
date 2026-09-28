"""Replay actual GR readouts and reconstruct the same constrained endpoint."""
from pathlib import Path
import json
import numpy as np
import def_radial_conduction_reconnect as work

out=work.OUT;task=work.old.task
read=lambda n:json.loads((out/n).read_text())
for p,h in read('plan.json')['bindings'].items():assert work.digest(Path(p))==h,p
result=read('result.json');hist={n:read('heat-'+str(n)+'.json')['history'] for n in [16,32,64]}
rows=[read(n+'.json') for n in ['heat-16','heat-32','heat-64','coefficient-64','outer-64']]
balance=max(r['heat_telescoping'] for r in rows);residual=max(r['linear_residual'] for r in rows)
assert balance<2e-13 and residual<1e-9
for field,expected in result['comparisons'].items():
    a,b,c=[np.array([r[field] for r in hist[n]]) for n in [16,32,64]]
    norm=max(abs(c).max(),1e-100);e1=float(max(abs(a-b[::2]))/norm);e2=float(max(abs(b-c[::2]))/norm)
    assert e1==expected['time_previous'] and e2==expected['time_last']
    assert float(np.log2(e1/e2))==expected['order']
assert not result['passed'] and not result['original_global_gates_passed']
ray=dict(np.load(task.coupled.OUT/'fine-rays.npz'));rad=task.Radiation(ray,out/'fine-bank.npz',False)
bg=task.coupled.Background(rad,2);state=dict(np.load(out/'heat-64.npz'))
assert np.array_equal(bg.grid,state['grid'])
r=bg.grid;q=state['response'];v=state['velocity'];src=rad.heat.source(1.,bg.nodes)
N,a=rad.geometry.metric(r);b=1/a**2;xi=r*q[:,0];A4=np.exp(-8*bg.nodes['phi']**2)
Dm=r*r*b*bg.nodes['v']*q[:,2]-(4*np.pi*r*r*A4*bg.nodes['p']+r*r*b*bg.nodes['v']**2/2)*xi+src[4]
dl=(np.divide(Dm,r,out=np.zeros_like(r),where=r>0)-np.divide(bg.nodes['m']*xi,r*r,out=np.zeros_like(r),where=r>0))/b
rho=np.interp(r,rad.heat.rnative,bg.native['raw'][::-1,0]);geo=task.h.gr.G*.1*rad.geometry.R**2/task.h.gr.C**4
cvT=np.interp(r,rad.heat.rnative,bg.native['thermo'][::-1,3])
logT=bg.nodes['adiabatic_T_rho']*q[:,1]/bg.nodes['gamma']-src[1]/(rho*geo*cvT)
speed=a/N*r*v[:,0]*task.h.gr.C
native=bg.native['radius_cm'][::-1]/(100*rad.geometry.R);weights=bg.native['dm'][::-1];weights/=weights.sum()
RMS=float(np.sqrt(weights@np.interp(native,r,speed)**2))
assert abs(RMS/result['endpoint']['velocity_mass_RMS_m_s']-1)<1e-12
assert all(np.isfinite(x).all() for x in [Dm,dl,logT,speed])
np.savez_compressed(out/'reconstructed-GR-endpoint.npz',grid=r,velocity_m_s=speed,Delta_log_T=logT,Delta_log_rho=q[:,1]/bg.nodes['gamma'],Delta_m_over_R=Dm,Delta_lambda=dl)
summary=dict(classification='Counterexample candidate',saved_replay_passed=True,reproduced_failure=True,
    max_heat_balance=balance,max_linear_residual=residual,
    endpoint_velocity_RMS=RMS,max_abs_Delta_log_T=float(max(abs(logT[r<1]))),
    max_abs_Delta_lambda=float(max(abs(dl[r<1]))),new_time_paths=0,new_native_EOS_calls=0,
    meaning='Replay confirms actual coupled GR endpoint and failed time-order verdict; it does not turn failure into acceptance.')
work.write(out/'verification.json',summary);print(summary)
