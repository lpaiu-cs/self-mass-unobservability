"""Measure exact banded solves before deciding whether direct evolution fits."""
from pathlib import Path
import json
import time
import signal
import numpy as np
import def_gr_spatial_banded as fast

base=fast.direct;OUT=base.OUT/'banded';assert not OUT.exists();OUT.mkdir();signal.alarm(90)
base.write(OUT/'plan.json',dict(classification='Counterexample candidate',
    reason='General sparse resolvent solves forecast1267s for the fixed comparison and were not launched. Use the already narrow bandwidth with pivoted band LU; do not change any equation, input or cutoff.',
    budget=dict(pilot_models=3,pilot_solves_each=16,hard_seconds=90,new_EOS_calls=0),
    gates=dict(component_transform_relative=1e-8,linear_residual=1e-9),
    bindings={str(p):base.digest(p) for p in [Path(__file__),Path(fast.__file__),Path(base.__file__),base.OUT/'pilot-budget.json']}))
measurements={}
for degree in [4,2,1]:
    start=time.monotonic();p=fast.Problem(degree);setup=time.monotonic()-start;start=time.monotonic()
    samples=np.array([p.transform(6+2j*np.pi*k/4) for k in range(16)]);elapsed=time.monotonic()-start
    reference=np.load(base.OUT/f'pilot-p{degree}.npy');fluid=np.zeros(p.model.size,bool);fluid[p.model.indices[:p.model.surface_index+1,0]]=True
    errors=[float(np.max(abs(samples[:,0,mask]-reference[:,0,mask]))/max(np.max(abs(reference[:,0,mask])),1e-100)) for mask in [fluid,~fluid]]
    assert max(errors)<1e-8,errors
    np.save(OUT/f'pilot-p{degree}.npy',samples)
    measurements[str(degree)]=dict(setup_seconds=setup,solve16_seconds=elapsed,linear_residual=p.error,component_transform_relative=errors,half_bandwidth=p.width)
    del p
forecast=1.35*sum((5 if degree=='4' else 1)*(v['setup_seconds']+2049/16*v['solve16_seconds']+5) for degree,v in measurements.items())
base.write(OUT/'budget.json',dict(classification='Counterexample candidate',measurements=measurements,forecast_seconds=forecast,
    full_workflow_allowed=forecast<300,maximum_seconds=300))
signal.alarm(0);print('BAND',json.dumps(measurements),'FORECAST',forecast,flush=True)
