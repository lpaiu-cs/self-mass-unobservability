"""Reconstruct temperature independently through canonical physical pressure."""
import json
import signal
import time
import numpy as np
import def_gr_temperature_feedback as task

OUT=task.OUT;assert not (OUT/'temperature-check.json').exists()
signal.alarm(60);start=time.monotonic();task.patch.install();p=task.Problem();m=p.model
saved=np.load(task.patch.OUT/'p4-2048.npz');q=saved['q'];E=saved['heat_energy']
r=m.original.native;data,point=task.go.task.fem.coefficients(m.bg,r)
V,D=m.evaluation(r);value=np.column_stack([v@q for v in V]);derivative=np.column_stack([v@q for v in D])
src=np.column_stack([v@E for v in task.go.task.fem.source_points(m.heat,point)])
A=data[:,:4].reshape(-1,2,2);B=data[:,4:8].reshape(-1,2,2)
g=np.einsum('nij,nj->ni',data[:,16:22].reshape(-1,2,3),src)
P=np.einsum('nij,nj->ni',B,derivative-np.einsum('nij,nj->ni',A,value)-g)
R=data[:,28:32].reshape(-1,2,2);S=data[:,32:36].reshape(-1,2,2)
eta=(P-np.einsum('nij,nj->ni',S,value)/2)[:,0]/R[:,0,0]
rr,loss,J=src.T;density=eta/point['gamma']+rr
gr=task.go.task.fem.base.task.h.gr;geo=gr.G*.1*m.heat.geometry.R**2/gr.C**4
rho=m.heat.d['raw'][::-1,0];cvT=m.heat.d['thermo'][::-1,5]
lag=point['adiabatic_T_rho']*density-loss/(rho*geo*cvT)
ids=np.clip(np.searchsorted(r,r,side='right')-1,0,len(r)-2)
slope=np.diff(m.heat.d['lnT'][::-1])[ids]/np.diff(r)[ids]
euler=lag-r*value[:,0]*slope;reference=p.Tq@q+p.TE@E
error=float(abs(euler-reference).max()/max(abs(euler).max(),1e-100));assert error<1e-12
row=dict(classification='Counterexample candidate',passed=True,
    canonical_pressure_temperature_relative=error,maximum_Eulerian_delta_lnT=float(abs(euler).max()),
    seconds=time.monotonic()-start,new_EOS_calls=0,new_time_paths=0,
    scope='Independent pressure-to-density-to-temperature reconstruction at all5735 original native material samples. Fixed background interpolation, not a physical EOS error enclosure.')
task.write(OUT/'temperature-check.json',row);signal.alarm(0);print(json.dumps(row))
