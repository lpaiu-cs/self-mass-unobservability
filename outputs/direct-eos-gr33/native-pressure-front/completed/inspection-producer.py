"""Counterexample candidate: pressure-mode coefficients on the saved problem."""
from pathlib import Path
import resource,time,json
import numpy as np
import apply_native_direct_radau as run

out=Path('native-pressure-front172-work');assert not out.exists();out.mkdir()
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.original.inf.incident.native.deadline(60)
start=time.monotonic();cpu=time.process_time()
run.write(out/'inspection-plan.json',dict(classification='Counterexample candidate',
    claim='Separate the local gas pressure reaction rate, coefficient variation and pulse arrival before selecting a repair of the actual171pressure mismatch.',
    decision='Use saved states and the same local operator. No physical integration, parameter scan or production. Do not infer full-system eigenvalues from a two-variable gas block.',
    budget=dict(seconds=60,CPU_threads=1,virtual_GiB=3),
    bindings={str(p):run.sha(p) for p in [Path(__file__),Path(run.__file__),run.OUT/'result.json',run.OUT/'source-check.json']}))
run.initialize();m=run.original.c.Response(128);h=m.t[-1]/64;end=m.t[-1]/16
p0=m.point(0)['pressure_map'];p1=m.point(1)['pressure_map'];pdot=(p1-p0)/end
rows=[]
for t in [0.,-float(m.redshift_driver.xc[15]/run.original.C),end]:
    c=m.local(t)
    D=np.stack([m.gas(c['B'][...,j],c['Bb'][...,j],c['Be'][...,j]) for j in range(2)],axis=2)
    for cell in [15,16]:
        a=c['pressure_map'][cell];ratio=a[1]/a[0];T=np.array([[1.,ratio],[0.,1.]])
        td=(pdot[cell,1]-ratio*pdot[cell,0])/a[0]
        transformed=(T@D[cell]+np.array([[0.,td],[0.,0.]]))@np.linalg.inv(T)
        ev=np.linalg.eigvals(D[cell])
        rows.append(dict(t=float(t),cell=cell,gas_matrix_times_h=(h*D[cell]).tolist(),
            gas_eigenvalues_times_h=[[float(v.real*h),float(v.imag*h)] for v in ev],
            pressure_coordinate_matrix_times_h=(h*transformed).tolist(),
            pressure_map=a.tolist(),pressure_map_relative_interval_change=((p1[cell]-p0[cell])/np.maximum(abs(p0[cell]),1e-290)).tolist()))
states=[np.load(run.OUT/f'sweep-1/photons/pilot-{n}.npz') for n in [64,128]]
rates=[];c=m.local(end)
for n,z in zip([64,128],states):
    x=z['delta_packet_scaled_occupation']/(m.scale*run.radau.AMP);g=z['delta_material']/run.radau.AMP
    photon=m.collision(c,x,np.zeros_like(g))[1];gas=m.collision(c,np.zeros_like(x),g)[1]
    contributions=[(np.sum(a*b,axis=1)*m.volume*run.radau.AMP) for a,b in [(c['pressure_map'],photon),(c['pressure_map'],gas),(pdot,g)]]
    rates.append(dict(steps=n,cells=[dict(cell=k,photon_rate=float(contributions[0][k]),gas_rate=float(contributions[1][k]),map_rate=float(contributions[2][k]),
        pressure=float(z['moments'][-1,6,k]),gas_pressure_terms=(c['pressure_map'][k]*g[k]*m.volume[k]*run.radau.AMP).astype(float).tolist()) for k in [15,16]]))
arrivals=-np.asarray(m.redshift_driver.xc,dtype=float)/run.original.C
result=dict(classification='Counterexample candidate',rows=rows,rates=rates,
    deep_cell_arrivals=[dict(cell=k,seconds=float(t)) for k,t in enumerate(arrivals[:m.nb])],
    prefix_end_seconds=float(end),deep_cells=m.nb,
    scope='Actual local gas block and saved endpoint response only. No full coupled spectral bound or unique-error-cause proof.')
run.write(out/'inspection.json',result)
run.write(out/'inspection-receipt.json',dict(action='inspection',seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
    peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=None,source_sha256=run.sha(__file__)))
print(json.dumps(result),flush=True)
