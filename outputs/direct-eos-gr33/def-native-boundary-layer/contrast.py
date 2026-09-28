"""Descriptive mixed time/radial contrast; no new accepted trajectory or gate."""
import json
import numpy as np
import def_native_boundary_layer as run

rows=[];differences=[]
for steps in [64,128]:
    a=np.load(run.prior.OUT/f'wave-{steps}.npz');b=np.load(run.OUT/f'wave-{steps}.npz');assert np.array_equal(a['t'],b['t'])
    differences.append(b['direct_relative']-a['direct_relative'])
    z=np.load(run.OUT/f'coupled-{steps}.npz')
    rows.append(dict(steps=steps,shared_mass_g=float(z['scalar_join_mass']),shared_energy_erg=float(z['scalar_join_energy'])))
scale=float(np.max(abs(b['direct_relative'])))
data=dict(classification='Counterexample candidate',new_fluid_steps=0,clock_identical=True,
    coarse_clock_radial_direct_relative=float(np.max(abs(differences[0]))/scale),
    fine_clock_radial_direct_relative=float(np.max(abs(differences[1]))/scale),
    mixed_difference_over_fine_direct=float(np.max(abs(differences[1]-differences[0]))/scale),
    local_ports=rows,shared_mass_time_relative=abs(rows[1]['shared_mass_g']-rows[0]['shared_mass_g'])/abs(rows[1]['shared_mass_g']),
    interpretation='Paired empirical contrasts only. A small charge comparison does not certify convergence of a tiny shared material port, global space/time error, or the physical charge.')
run.write(run.OUT/'contrast.json',data);print(json.dumps(data))
