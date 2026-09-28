"""Replay the coordinate repair and amplitude-rescale the tiny motion response."""
from pathlib import Path
import json
import resource
import signal
import time
import numpy as np
from scipy.interpolate import CubicHermiteSpline
import def_native_moving_rays as motion

task=motion.task
OUT=task.OUT


def main():
    assert not (OUT/'audit.json').exists();start=time.monotonic()
    result=json.loads((OUT/'motion-result.json').read_text())
    saved=json.loads((OUT/'motion-p4-64.json').read_text())
    previous=result['accounted_total_seconds'];forecast=previous+1.5*(saved['seconds']+saved['setup_seconds']+5)
    assert forecast<180,forecast
    # The source is extremely small. A same-grid amplitude rescaling checks
    # that absolute residual floors did not manufacture its measured response.
    task.write(OUT/'audit-plan.json',dict(classification='Counterexample candidate',
        claim='Verify source-coordinate equivalence, current-nullspace restoration and the tiny additive response under exact linear amplitude rescaling.',
        path='One same p4-64 path with unit maximum surface displacement, restoring physical amplitude only during comparison.',
        rationale='Absolute 1e-100 residual denominators may mask tiny equations. A dimensionless unit driver tests the numerical result without using that floor as accuracy evidence.',
        prior_seconds=previous,forecast_total_seconds=forecast,total_budget_seconds=180,
        gates=dict(source_maps_relative=2e-12,amplitude_relative=1e-7),
        bindings={str(p.relative_to(task.old.ROOT)):task.old.photons.digest(p) for p in [Path(__file__),Path(motion.__file__),Path(task.__file__),OUT/'motion-result.json',OUT/'motion-p4-64.npz',OUT/'p4-64.npz']}))
    signal.alarm(int(180-previous));resource.setrlimit(resource.RLIMIT_AS,(int(5e9),int(5e9)))
    m=task.Model();p=m.bg.sample(m.native)
    rng=np.random.default_rng(96);face=rng.normal(size=m.n).astype(np.longdouble)
    direct=m.source_values(p,face);maps=m.sources(p)
    source_errors=[float(max(abs(a-b@face))/max(max(abs(a)),1e-100)) for a,b in zip(direct,maps)]
    assert max(source_errors)<2e-12,source_errors
    null=m.source_values(p,np.ones(m.n,np.longdouble))[1];assert max(abs(null[2:]))==0
    a=np.load(OUT/'p4-64.npz');b=np.load(OUT/'p4-32.npz')
    amplitude=float(max(abs(a['surface'][:,0])));assert amplitude>0
    hist=CubicHermiteSpline(a['emission_times'],a['surface'][:,0]/amplitude,a['surface'][:,2]/amplitude)
    rays=motion.Rays(m.bg,m.points['r'][m.outside]);f0=float(m.f0[-1]);R=m.bg.R
    m.F0*=0;m.TE0*=0;m.LE0*=0;m.J0*=0
    m.photon_force=lambda t:rays.force(t,hist,f0,R,m.photon_test)
    report=m.evolve(float(a['emission_times'][-1]*m.bg.tc),64,'motion-unit-p4-64')
    c=np.load(OUT/'motion-unit-p4-64.npz');d=np.load(OUT/'motion-p4-64.npz');relative={}
    for f in ['temperature','velocity','scalar','surface','q','v','e','d']:
        relative[f]=float(np.max(abs(c[f]*amplitude-d[f]))/max(np.max(abs(d[f])),1e-100))
    assert max(relative.values())<1e-7,relative
    surface_error=np.max(abs(b['surface']-a['surface'][::2]),axis=0)/np.maximum(np.max(abs(a['surface']),axis=0),1e-100)
    old=np.load(task.prior.OUT/'consistent-p4-64.npz');q=old['q'];v=old['w']-m.H@old['f']
    old_surface=[float((V@z)[0]) for z in [q,v] for V in m.surfaceV]
    reference=np.load(task.prior.OUT/'consistent-p6-64.npz')
    remaining=float(np.max(abs(a['velocity']-reference['velocity']))/np.max(abs(reference['velocity'])))
    reports=[json.loads((OUT/f'{label}.json').read_text()) for label in ['p4-8','p4-32','p4-64','motion-p4-8','motion-p4-32','motion-p4-64','motion-unit-p4-64']]
    for label in ['p4-8','p4-32','p4-64','motion-p4-8','motion-p4-32','motion-p4-64','motion-unit-p4-64']:
        z=np.load(OUT/f'{label}.npz');assert all(np.all(np.isfinite(z[k])) for k in z.files)
        assert all(np.max(abs(z[k][0]))==0 for k in ['temperature','velocity','scalar','surface'])
    # Replay immutable source bindings for both production plans.
    for file in ['plan.json','motion-plan.json','audit-plan.json']:
        for path,digest in json.loads((OUT/file).read_text())['bindings'].items():assert task.old.photons.digest(task.old.ROOT/path)==digest,path
    out=dict(classification='Counterexample candidate',checks_passed=True,symbolic=task.controls(),
        source_map_relative=source_errors,uniform_current_local_debit=0,amplitude_rescale_relative=relative,
        amplitude=amplitude,physical_surface_time_relative=surface_error.tolist(),old_surface_endpoint=old_surface,
        repaired_surface_endpoint=a['surface'][-1].tolist(),saved_p6_velocity_space_relative=remaining,
        original_space_failure_resolved=False,
        max_linear_residual=max(r['max_linear_residual'] for r in reports),
        max_local_deviation_heat_identity=max(r['max_local_deviation_heat_identity'] for r in reports),
        memory_GB=max(r['memory_GB'] for r in reports),seconds=time.monotonic()-start,
        accounted_total_seconds=previous+time.monotonic()-start,
        scope='Unit-driver outputs are a linear arithmetic control, not a unit-displacement physical star. Restore the saved amplitude before physical interpretation. The kinematic L0*zeta contribution is one mixed-order term; other terms of that order remain open.',
        moving_surface_fully_solved=False,final_charge_solved=False,full_goal_complete=False)
    task.write(OUT/'audit.json',out);signal.alarm(0);print(json.dumps(out),flush=True)


if __name__=='__main__':main()
