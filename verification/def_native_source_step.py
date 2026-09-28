"""Native-EOS source-step pairs with scalar backreaction and conserved work.

Counterexample candidate: fixed-volume representative-cell source integration.
It is an operator control, NOT a star, scalar exterior, thermal mode, or GR run.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

import def_native_coupling as coupling

ld = np.longdouble
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/direct-eos-gr33/def-native-source-step'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(name, value):
    (OUT/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def prepare():
    assert not OUT.exists()
    coupling.initial.molecular.bindings()
    files = [Path(__file__), Path(coupling.__file__),
        coupling.initial.molecular.OUT/'molecular-state-17-4.npz',
        coupling.initial.molecular.OUT/'manifest.json',
        Path(coupling.initial.molecular.model.LIB)]
    OUT.mkdir()
    save('plan.json', dict(classification='Counterexample candidate',
        checkpoint='3065706d', bindings={str(p.resolve()):digest(p) for p in files},
        cells=[0,1175,1176,2972,5734], steps=[16,32,64], duration_scaled=6.283185307179586,
        beta=-4, background_phi=.001, drive_amplitude=.0002,
        scalar_stiffness_over_rest_energy=16,
        pilot=dict(cell=2972,steps=16),
        source_problem='Fixed Einstein volume, baryon density and nuclear composition. Native molecular entropy is constant during this isolated reversible source substep. A scalar degree of freedom has kinetic energy K*p^2/2, spring K*(phi-phi0)^2/2 and an external constant counterforce balancing the native matter derivative at phi0. An applied sinusoidal force supplies separately measured work. The scalar and actual native matter energy exchange through an exact specific-energy divided difference.',
        method='Implicit midpoint for scalar velocity plus native matter discrete energy gradient. Solve entropy for temperature at each iterate. Paired unforced path and sign-reversed background/drive controls. No post-step energy rescaling.',
        purpose='Validate the reusable native EOS frame map and reciprocal source update before adding it to the spherical constrained evolution. Artificial scalar stiffness and scaled time are controls, not fitted stellar-mode parameters.',
        gates=dict(native_entropy_temperature_residual=2e-12,
            stage_absolute_residual=2e-18, normalized_total_work_defect=1e-9,
            field_time_refinement_order=1.8, native_derivative_relative_difference=1e-5,
            maximum_stage_iterations=8, zero_drive_drift=0., parity_difference=1e-13),
        budget=dict(workers=1,gpu=False,measured_native_seconds=[.02377,.02503,.02172,.04713,.05241],
            estimated_native_calls_maximum=6000, estimated_wall_seconds=[180,600],
            hard_timeout_seconds=600, maximum_run_repeats=1,
            stop='Stop and preserve any native/conservation/refinement failure or time cap; no automatic larger grid, more steps, stronger amplitude or relaxed acceptance.'),
        full_GR_evolution=False, physical_stellar_drive=False, physical_EOS_certified=False,
        transport_included=False, observational_closure=False))
    save('symbolic.json', coupling.symbolic())
    print('PREPARED native reciprocal source-step pairs', flush=True)


def bindings():
    plan = json.loads((OUT/'plan.json').read_text())
    for name, value in plan['bindings'].items():
        assert digest(name) == value, name
    return plan


def derivative(z):
    return -4*z['phi']*z['A']*(z['restJ']+z['uJ']-3*ld(z['raw'][1])/z['rhoJ'])


def curvature(z):
    w = ld(z['raw'][1])/z['rhoJ']
    gamma = (w-ld(z['raw'][9]))/ld(z['raw'][10])
    chi_s = ld(z['raw'][5])+ld(z['raw'][6])*gamma
    h = z['restJ']+z['uJ']-3*w
    return -4*z['A']*h+16*z['phi']**2*z['A']*(h+w*(9*chi_s-12))


def path(provider, data, cell, steps, sign=1, driven=True):
    plan = json.loads((OUT/'plan.json').read_text())
    rho, T, X = data['lnd'][cell], data['lnT'][cell], data['X'][cell]
    phi0 = ld(sign)*ld(str(plan['background_phi']))
    amplitude = ld(sign)*ld(str(plan['drive_amplitude'])) if driven else ld(0)
    z0 = provider.state(rho,T,phi0,X); z = z0
    K = ld(plan['scalar_stiffness_over_rest_energy'])*z0['restJ']
    g0 = derivative(z0); p = ld(0); work = ld(0)
    h = ld(str(plan['duration_scaled']))/steps
    history = [[0.,float(phi0),0.,float(T),0.,0.]]
    maximum_defect, maximum_entropy, maximum_residual, maximum_iterations = 0.,0.,0.,0
    for step in range(steps):
        drive = amplitude*np.sin((ld(step)+ld('.5'))*h)
        dphi = (h*p-h*h/2*(z['phi']-phi0+(derivative(z)-g0)/K-drive))/(1+h*h/4)
        for iteration in range(plan['gates']['maximum_stage_iterations']):
            phi = z['phi']+dphi
            guess = z['logT']+z['isentropic_dlogT_dlogA']*(-2*phi**2-z['logA'])
            trial = provider.isentropic(rho,phi,X,z0['entropy'],guess)
            du = coupling.matter_specific_increment(z,trial)
            gradient = derivative(z) if dphi == 0 else du/dphi
            residual = dphi-h*p+h*h/2*(z['phi']+dphi/2-phi0+(gradient-g0)/K-drive)
            if abs(residual) <= plan['gates']['stage_absolute_residual']:
                break
            dphi -= residual/(1+h*h/4*(1+curvature(trial)/K))
        else:
            raise RuntimeError(('Reciprocal source step failed',cell,steps,step,float(residual)))
        p1 = 2*dphi/h-p
        work += K*drive*dphi
        matter = coupling.matter_specific_increment(z0,trial)
        displacement = trial['phi']-phi0
        scalar = K/2*(p1*p1+displacement*displacement)-g0*displacement
        scale = K*ld(str(plan['drive_amplitude']))**2
        defect = float(abs(matter+scalar-work)/scale)
        maximum_defect = max(maximum_defect,defect)
        maximum_entropy = max(maximum_entropy,float(abs(trial['entropy_temperature_residual'])))
        maximum_residual = max(maximum_residual,float(abs(residual)))
        maximum_iterations = max(maximum_iterations,iteration+1)
        assert defect <= plan['gates']['normalized_total_work_defect'], (cell,step,defect)
        z,p = trial,p1
        history.append([float((step+1)*h),float(z['phi']),float(p),float(z['logT']),float(work),defect])
    array = np.asarray(history)
    label = f'cell-{cell}-steps-{steps}-sign-{sign}-driven-{int(driven)}'
    np.savez_compressed(OUT/(label+'.npz'), history=array)
    return dict(cell=cell,steps=steps,sign=sign,driven=driven,
        final_phi=float(z['phi']), final_velocity=float(p), final_logT=float(z['logT']),
        maximum_field_change=float(np.max(abs(array[:,1]-array[0,1]))),
        maximum_temperature_change=float(np.max(abs(array[:,3]-array[0,3]))),
        maximum_normalized_work_defect=maximum_defect,
        maximum_native_entropy_temperature_residual=maximum_entropy,
        maximum_stage_residual=maximum_residual,maximum_stage_iterations=maximum_iterations,
        history_file=(OUT/(label+'.npz')).relative_to(ROOT).as_posix())


def controls(provider,data):
    rows=[]
    for cell in [0,1175,1176,2972,5734]:
        rho,T,X = data['lnd'][cell],data['lnT'][cell],data['X'][cell]
        zero = provider.state(rho,T,0.,X)
        original = provider.native(2,float(rho),float(T),np.asarray(X,float))
        assert np.array_equal(zero['raw'],original)
        assert zero['P'] == ld(original[1])
        centre = provider.state(rho,T,ld('.001'),X)
        differences=[]
        for step in [ld('.00002'),ld('.00001')]:
            minus=provider.isentropic(rho,centre['phi']-step,X,centre['entropy'],centre['logT'])
            plus=provider.isentropic(rho,centre['phi']+step,X,centre['entropy'],centre['logT'])
            numerical=(coupling.matter_specific_increment(centre,plus)-coupling.matter_specific_increment(centre,minus))/(2*step)
            differences.append(float(abs(numerical/derivative(centre)-1)))
        assert max(differences)<1e-5,(cell,differences)
        rows.append(dict(cell=cell,zero_scalar_native_array_exact=True,isentropic_derivative_relative_differences=differences))
    return rows


def run(pilot=False):
    plan=bindings(); began=time.monotonic()
    output='pilot.json' if pilot else 'result.json'
    assert not (OUT/output).exists()
    data=dict(np.load(coupling.initial.molecular.OUT/'molecular-state-17-4.npz'))
    provider=coupling.Matter()
    if pilot:
        start=provider.calls
        row=path(provider,data,plan['pilot']['cell'],plan['pilot']['steps'])
        result=dict(classification='Counterexample candidate',passed=True,row=row,
            seconds=time.monotonic()-began,native_calls=provider.calls-start)
        save(output,result);print(json.dumps(result,indent=2));return
    pilot_result=json.loads((OUT/'pilot.json').read_text());assert pilot_result['passed']
    control=controls(provider,data); rows=[]
    for cell in plan['cells']:
        for steps in plan['steps']:
            if cell==plan['pilot']['cell'] and steps==plan['pilot']['steps']:
                row=pilot_result['row']
            else: row=path(provider,data,cell,steps)
            rows.append(row)
            save('partial.json',dict(classification='Counterexample candidate',rows=rows))
            assert time.monotonic()-began<plan['budget']['hard_timeout_seconds']
        rows.append(path(provider,data,cell,plan['steps'][0],driven=False))
        print('COMPLETED native coupled source cell',cell,provider.calls,flush=True)
    negative=path(provider,data,2972,64,sign=-1)
    positive=next(r for r in rows if r['cell']==2972 and r['steps']==64 and r['driven'])
    left=np.load(ROOT/positive['history_file'])['history'];right=np.load(ROOT/negative['history_file'])['history']
    parity=float(max(abs(left[:,1:3]+right[:,1:3]).max(),abs(left[:,3]-right[:,3]).max()))
    refinements=[]
    for cell in plan['cells']:
        chosen=[r for r in rows if r['cell']==cell and r['driven']]
        values=np.array([[r['final_phi'],r['final_velocity']] for r in chosen])
        errors=np.max(abs(np.diff(values,axis=0)),axis=1)
        order=float(np.log2(errors[0]/errors[1]))
        refinements.append(dict(cell=cell,endpoint_field_differences=errors.tolist(),order=order,
            passed=bool(errors[1]<errors[0] and order>=plan['gates']['field_time_refinement_order'])))
    zero=all(r['maximum_field_change']==0 and r['maximum_temperature_change']==0 for r in rows if not r['driven'])
    result=dict(classification='Counterexample candidate',completed=True,
        passed=zero and parity<plan['gates']['parity_difference'] and all(r['passed'] for r in refinements),
        native_provider_controls=control,rows=rows,refinements=refinements,
        zero_drive_exactly_stationary=zero,field_sign_temperature_even_parity_difference=parity,
        sign_reversed_path=negative,seconds=time.monotonic()-began,native_calls=provider.calls,
        full_GR_evolution=False,physical_stellar_drive=False,physical_EOS_certified=False,
        transport_included=False,observational_closure=False)
    save(output,result)
    save('manifest.json',dict(sha256={p.relative_to(ROOT).as_posix():digest(p) for p in OUT.iterdir()
        if p.is_file() and p.name!='manifest.json'}))
    assert result['passed'],result
    print('PASS native EOS reciprocal source pairs',json.dumps({k:v for k,v in result.items() if k not in ['rows','native_provider_controls','sign_reversed_path']}),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','pilot','run'])
    command=parser.parse_args().command
    if command=='prepare':prepare()
    else:run(command=='pilot')
