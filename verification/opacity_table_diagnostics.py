"""Locate the failed native-value match without changing its tolerance."""
import json
import numpy as np
import opacity_tables as t

model=t.Opacity();source=dict(np.load(t.o.OUT/'tables-baseline-captured.npz'))
values=dict(np.load(t.OUT/'evaluation.npz'));indices=np.argsort(values['value_relative_errors'])[-12:][::-1];rows=[]
for i in indices:
    p=source['parameters'][i];z=p[2];grid=model.zgrid
    if z<=grid[0]: rad=model.hydrogen(0,p[1],p[3],p[4])
    elif z>=grid[-1]: rad=model.hydrogen(len(grid)-1,p[1],p[3],p[4])
    else:
        j=int(np.searchsorted(grid,z,side='right'))-1;u=(z-grid[j])/(grid[j+1]-grid[j])
        rad=(1-u)*model.hydrogen(j,p[1],p[3],p[4])+u*model.hydrogen(j+1,p[1],p[3],p[4])
    new_rad=10**rad[0];old_rad=source['inner'][i,8]
    quantized=p.copy();quantized[:5]=np.float32(quantized[:5]).astype(float)
    rows.append(dict(cell=int(i),parameters=p[:5].tolist(),native=source['outputs'][i].tolist(),
        candidate=values['values'][i].tolist(),relative_total_value_change=float(values['value_relative_errors'][i]),
        relative_radiative_value_change=float(new_rad/old_rad-1),
        relative_total_change_at_rounded_inputs=float(model(quantized)[0]/source['outputs'][i,0]-1),
        native_radiative_fraction=float(source['outputs'][i,0]/old_rad)))
t.save('value-diagnostics.json',dict(classification='Counterexample candidate',
    original_value_gate_passed=False,failed_value_cells=int(np.count_nonzero(values['value_relative_errors']>1e-4)),
    rows=rows,root_cause_fully_identified=False))
print(json.dumps(rows[:3],indent=2),flush=True)
