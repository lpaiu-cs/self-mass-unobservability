"""Apply the saved radial conduction bank to the original failed GR path.

Counterexample candidate: frozen-source linear fluid/scalar/metric response.
Unknown outer currents and photons remain omitted components, not zero-flux
physical boundary conditions. Keep the old failure and its acceptance gates.
"""
from pathlib import Path
import argparse
import inspect
import json
import resource
import signal
import time
import numpy as np
import def_conduction_radau as old
import def_conduction_resolved_bank as radial

OUT=radial.OUT.parent/'def-radial-conduction-reconnect'
write=radial.ex.write;digest=radial.h.digest


def solver_source():
    # Same actual GR equations, source debit, momentum lift and Radau solver.
    # Add local readouts around BOTH the original and new support endpoints.
    source=inspect.getsource(old.solve)
    anchor='    r=bg.grid;N,a=radiation.geometry.metric(r)'
    source=source.replace(anchor,anchor+'''
    native_faces=bg.native['faces_cm']/(100*radiation.geometry.R)
    width=538000/radiation.geometry.R
    masks=[abs(native-native_faces[i])<=width for i in [4123,1723]]
    assert all(np.any(mask) for mask in masks)
''')
    anchor="return dict(tau=t,velocity_mass_RMS_m_s=float(np.sqrt(weights@(speed*speed))),scalar_mass_RMS=float(np.sqrt(weights@(scalar*scalar)))),q,p"
    assert source.count(anchor)==1
    source=source.replace(anchor,"return dict(tau=t,velocity_mass_RMS_m_s=float(np.sqrt(weights@(speed*speed))),scalar_mass_RMS=float(np.sqrt(weights@(scalar*scalar))),**{name:float(np.sqrt((weights[mask]@speed[mask]**2)/weights[mask].sum())) for name,mask in zip(['old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s'],masks)}),q,p")
    return source


def solve(steps,bank,outer,label):
    source=(OUT/'solver-source.py').read_text()
    ns=dict(vars(old),OUT=OUT);exec(compile(source,str(OUT/'solver-source.py'),'exec'),ns)
    return ns['solve'](steps,bank,False,outer,label)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(old.__file__),Path(old.task.__file__),Path(old.task.coupled.__file__),
        radial.OUT/'partial-bank.npz',radial.OUT/'result.json',radial.state.OUT/'partial-bank.npz',
        old.OUT/'result.json',old.task.coupled.surface.OUT/'background.npz',
        old.task.coupled.thermal.OUT/'coefficients.npz']
    assert json.loads((radial.OUT/'result.json').read_text())['passed']
    assert not json.loads((old.OUT/'result.json').read_text())['passed']
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='af9b712d',
        claim='Apply the already accepted Phase58 radial transport input to the actual Phase56 GR fluid/scalar/metric evolution and retest its velocity failure. The original artificial core cut is crossed using calculated coefficients.',
        decision='Report original global gates and local velocity time convergence at the old and new support boundaries separately. A larger new source must not hide the failed core component. Do not claim a full physical outer boundary, complete photon input or nonlinear feedback.',
        model='Unchanged frozen-source first-order GR equations, native background, microscopic pole startup, energy debit, heat momentum lift, spatial mesh and0.23080495568542375s horizon. Replace only1612 core faces by the saved4012 supported faces. Other source components remain unknown, not certified zero.',
        paths=['heat-16','heat-32','heat-64','coefficient-64','outer-64'],
        coefficient_contrast='Phase58 prior28-anchor eta-mixture bank on identical faces; this previously failed its own local tau gate and is only a sensitivity contrast, not an independent error certificate.',
        gates=dict(time_relative=.02,time_order=1.5,coefficient_relative=.02,outer_relative=.002,linear_residual=1e-9,heat_balance=2e-13),
        budget=dict(paths=5,steps=240,hard_seconds=180,cpu_threads=1,memory_GB=3,new_native_EOS_calls=0,new_collision_states=0,automatic_expansion=False),
        forecast='Prior six same-grid Radau paths107.68s. Five paths forecast80-160s; larger12-pole radial source evaluation is unmeasured. Hard180s. Reuse all input and no new pilot/long stellar run.',
        bindings={str(p):digest(p) for p in paths}))
    for name,p in [('fine',radial.OUT/'partial-bank.npz'),('coarse',radial.state.OUT/'partial-bank.npz')]:
        b=dict(np.load(p));assert np.array_equal(b['known_faces'],np.arange(1723,5735))
        np.savez_compressed(OUT/(name+'-bank.npz'),faces=b['known_faces'],mode_K_SI=b['mode_K_SI'],poles_proper_s=b['poles_proper_s'],unknown_faces=b['unknown_faces'])
    write(OUT/'symbolic.json',old.task.symbolic())
    (OUT/'solver-source.py').write_text(solver_source())


def run():
    assert not (OUT/'result.json').exists();signal.alarm(180)
    resource.setrlimit(resource.RLIMIT_AS,(int(3e9),int(3e9)))
    plan=json.loads((OUT/'plan.json').read_text())
    for p,h in plan['bindings'].items():assert digest(Path(p))==h,p
    begin=time.monotonic();cases={}
    for n in [16,32,64]:cases[str(n)]=solve(n,OUT/'fine-bank.npz',2,'heat-'+str(n))
    cases['coarse']=solve(64,OUT/'coarse-bank.npz',2,'coefficient-64')
    cases['outer']=solve(64,OUT/'fine-bank.npz',3,'outer-64')
    comparisons={}
    for field in ['velocity_mass_RMS_m_s','scalar_mass_RMS','old_interface_velocity_RMS_m_s','new_interface_velocity_RMS_m_s']:
        series=lambda name:np.array([x[field] for x in cases[name]['history']])
        a,b,c=series('16'),series('32'),series('64');norm=max(abs(c).max(),1e-100)
        d1=np.max(abs(a-b[::2]))/norm;d2=np.max(abs(b-c[::2]))/norm
        comparisons[field]=dict(time_previous=float(d1),time_last=float(d2),order=float(np.log2(d1/d2)),
            coefficients=float(np.max(abs(c-series('coarse')))/norm),outer=float(np.max(abs(c-series('outer')))/norm))
    passed=all(x['time_last']<.02 and x['order']>1.5 and x['coefficients']<.02 and x['outer']<.002 for x in comparisons.values())
    passed=passed and max(x['heat_telescoping'] for x in cases.values())<2e-13
    original=json.loads((old.OUT/'result.json').read_text())
    result=dict(classification='Counterexample candidate',actual_GR_fluid_scalar_metric_evolved=True,passed=passed,
        original_global_gates_passed=all(comparisons[k]['time_last']<.02 and comparisons[k]['order']>1.5 and comparisons[k]['coefficients']<.02 and comparisons[k]['outer']<.002 for k in ['velocity_mass_RMS_m_s','scalar_mass_RMS']),
        comparisons=comparisons,previous_failed_comparisons=original['comparisons'],endpoint=cases['64']['history'][-1],
        seconds=time.monotonic()-begin,memory_GB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024/1e9,
        original_cut_crossed_with_computed_input=True,new_outer_support_cut_remains=True,
        whole_star_heat_closed=False,full_temperature_feedback=False,photon_transport=False,full_dynamic_charge_solved=False)
    write(OUT/'result.json',result);signal.alarm(0);print('GR RECONNECTED',json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
