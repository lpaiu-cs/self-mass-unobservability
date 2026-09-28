"""Verify geometry and each separated charge component on saved results.

Counterexample candidate. Exact binary inputs,100-digit readout arithmetic;
neither input accuracy nor the physical charge is certified to100 digits.
"""
from pathlib import Path
import json, resource, time
import numpy as np
import mpmath as mp
import read_native_incident_infinity as run


def main():
    dest=run.OUT/'normalization-detail.json';assert not dest.exists()
    files=[Path(__file__),run.OUT/'final-result.json',run.OUT/'charge-parts.npz',
           run.OUT/'body-128-a8-r8.npz',run.OUT/'self-128-a8-r8.npz',run.OUT/'direct-g8.npz',
           run.BEFORE/'gr/source-128-reference-128.npz',run.BACKGROUND]
    run.write(run.OUT/'normalization-audit-plan.json',dict(classification='Counterexample candidate',
        claim='Verify each separated physical charge component and coordinate matching, without a direct-field-dominated error norm. Retain background normalization corrections smaller than a longdouble ulp.',
        budget_seconds=30,remaining_from_original_450_seconds=True,new_evolution_steps=0,
        gates=dict(component_relative=1e-12,geometry_relative=1e-14),
        bindings={str(p):run.sha(p) for p in files}))
    run.incident.native.deadline(30);start=time.monotonic();cpu=time.process_time()
    run.prior.initialize();driver=run.incident.Driver(8);m=run.exterior.Exterior()
    source=np.load(run.BEFORE/'gr/source-128-reference-128.npz')
    # Saved source edges are Jordan radii; the actual ray launch and optical
    # incoming data use the Einstein areal radius returned by the same map.
    mapped=driver.geo.metric(source['edges'][-1:]-driver.model.m.RJ)[3][0]*driver.model.m.R
    assert abs(mapped/driver.r0-1)<1e-14
    assert driver.r0==m.r0 and float(source['M_cm'])==m.M==float(driver.bg.z['ADM_mass'])
    assert abs(float(source['K_cm'])/m.K-1)<1e-14 and driver.T==m.T
    geometry=dict(classification='Counterexample candidate',passed=True,
        source_outer_Jordan_radius_cm=float(source['edges'][-1]),ray_outer_Einstein_radius_cm=m.r0,
        mapped_source_outer_radius_cm=float(mapped),mass_cm=m.M,scalar_K_cm=m.K,horizon_seconds=m.T)
    d=np.load(run.OUT/'charge-parts.npz');parts=[np.load(run.OUT/f'{tag}-128-a8-r8.npz') for tag in ['body','self']]
    raw=np.load(run.OUT/'direct-g8.npz')['normalized'];mp.mp.dps=100
    def B(value):
        n,den=value.as_integer_ratio();return mp.mpf(n)/mp.mpf(den)
    rows=[];errors=[]
    for j in range(len(d['t'])):
        a=B(d['alpha0'][()]);e0=B(d['background_epsilon'][j]);q0=B(d['background_normalized'][j])
        de=sum((B(p['epsilon_increment'][j]) for p in parts),mp.mpf(0))
        numerators=[B(raw[j])]+[B(p['compact'][j])+B(p['exterior'][j])+(a+q0)*B(p['epsilon_increment'][j]) for p in parts]
        linear=[v/(1-e0) for v in numerators];rational=[v/(1-e0-de) for v in numerators]
        stored=[B(d[key][j]) for key in ['direct','body','self_gr']]
        errors.append([float(abs(v-w)/max(abs(v),mp.mpf('1e-290'))) for v,w in zip(linear,stored)])
        s0=q0*(1-e0)-a*e0
        ds=B(raw[j])+sum((B(p['compact'][j])+B(p['exterior'][j]) for p in parts),mp.mpf(0))
        exact=(s0+ds+a*(e0+de))/(1-e0-de)-q0
        assert abs(exact-sum(rational))<mp.mpf('1e-110')
        rows.append(dict(time=float(d['t'][j]),linear_parts=[str(v) for v in linear],
            background_normalization_corrections=[str(v*e0/(1-e0)) for v in numerators],
            rational_cross_terms=[str(v*de/(1-e0-de)) for v in linear],rational_parts=[str(v) for v in rational],
            exact_rational_total_increment=str(exact)))
    errors=np.asarray(errors);assert errors.max()<1e-12
    compact=[]
    for folder in [run.BEFORE,run.SELF/'sweep-2']:
        fine=np.load(folder/'gr/wave-128-g8.npz')['free_scalar'];coarse=np.load(folder/'gr/wave-128-g4.npz')['free_scalar']
        compact.append(float(np.max(abs(fine-coarse))/max(np.max(abs(fine)),1e-290)))
    assert max(compact)<.002
    result=dict(classification='Counterexample candidate',passed=True,geometry=geometry,
        per_component_normalization_max_relative=errors.max(axis=0).tolist(),
        inherited_compact_quadrature_relative=compact,exact_binary_inputs=True,arithmetic_digits=100,rows=rows,
        original_first_Born_input_was_applied_to_evolution=True,new_or_higher_Born_input_applied=False,
        raw_flag_clarification='direct-result.first_return_refed_to_fluid=false means this readout did not change or refeed the existing input; it does not retract the Phase155 applied first-Born input.',
        scope='Separated fixed-input readout arithmetic only; extra printed digits are not a physical uncertainty estimate.',
        seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
        peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,full_goal_complete=False)
    run.write(dest,result)
    print(json.dumps({k:v for k,v in result.items() if k!='rows'}),flush=True)


if __name__=='__main__':
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));main()
