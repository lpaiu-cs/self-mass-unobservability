"""Finite-temperature nonzero-background DEF evolution without HSE projection.

Counterexample candidate: a finite, reflecting material cavity with prescribed
scalar boundary. It is not the free-surface/radiative stellar response.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import time
import numpy as np
import sympy as sp
import def_normalized_charge as q

s,e,ld=q.s,q.e,q.ld
OUT=e.g.OUT/'def-nonstationary-interior'


def initialize(pool,h):
    star,_=s.initialize(pool,-4,ld('.001'),h)
    saved=dict(np.load(q.OUT/'finite.npz'))
    delta=saved['delta'].copy();phi=ld('.001')*saved['normalized_field']
    pi=np.zeros(star.n,dtype=ld);phir=saved['Phi'].copy()
    z=star.fluid(delta,phi,star.native_aux(delta,phi))
    # The old zero-field reference remains an arithmetic reference, not a
    # physical equilibrium. Reconstruct the *instantaneous* metric with p=z.
    del star.previous
    z.update(s.prior.metric(star,z['dEm'],z['R'],pi,phir),Pi=pi,Phi=phir)
    star.previous=z
    z.update(s.metric(star,z['dEm'],z['R'],pi,phir),scalar_residual=0.)
    star.finish_moments(delta,z,z['dEm']);z['at']=np.zeros(star.n,dtype=ld)
    star.previous=z;star.field_guess=phi.copy(),pi.copy(),phir.copy()
    star.initial_guess=delta.copy()
    boundary=phi[-1]+phir[-1]*star.distance[-1]
    star.amplitude=boundary;star.previous_boundary=boundary
    star.equilibrium_pressure_flux=np.zeros(star.n+1,dtype=ld)
    star.equilibrium_gravity=np.zeros(star.n,dtype=ld)
    star.equilibrium_pressure=np.zeros(star.n,dtype=ld)
    baryon=abs(z['B']/star.B0-1).max()
    assert baryon<2e-13,baryon
    return star,delta,z,boundary


def symbolic():
    h,c,div,momentum,at,S,N,nur,E,P,geom,alpha,trace,phir,scale=sp.symbols(
        'h c div momentum at S N nur E P geom alpha trace phir scale')
    force=-at*S-c*N*nur*E+c*N*P*geom+c*N*alpha*trace*phir
    lhs=(momentum+h*c*div-h*force)/scale
    assert sp.expand(lhs*scale-momentum-h*(c*div+at*S+c*N*nur*E-c*N*P*geom-c*N*alpha*trace*phir))==0
    # A nonzero heat momentum is not compatible with a stationary metric.
    r,g,a,v,w,Q=sp.symbols('r g a v w Q',positive=True)
    flux=(w*v+Q*(1+v*v))/(1-v*v)
    mdot=-4*sp.pi*g*c*r*r*N*flux/a
    assert sp.simplify(mdot.subs(v,0)+4*sp.pi*g*c*r*r*N*Q/a)==0
    return dict(classification='Proven',passed=True,
        scope='Removing the three reference momentum terms recovers the unprojected finite-volume radial momentum balance. Nonzero heat flux requires evolving mass at zero material velocity. No continuum error certificate.')


def prepare():
    assert not OUT.exists();OUT.mkdir()
    files=[Path(__file__),Path(s.__file__),Path(s.old.__file__),Path(s.prior.__file__),
        Path(s.m.__file__),Path(s.m.prior.__file__),Path(s.m.method.__file__),
        q.OUT/'finite.npz',q.OUT/'manifest.json',s.OUT/'initial.npz',s.OUT/'initial-manifest.json']
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        checkpoint='7950a546',bindings={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in files},
        symbolic=symbolic(),cells=24,phi_infinity=.001,drive_amplitude=.00001,
        time_steps=[12,24,48],paths=['minus','undriven','plus'],
        duration_crossing_times=2,pulse_duration_crossing_times=1,
        drive='Same finite nonzero initial constrained state; add signed sin^8 pulse to the common scalar surface value. Reflecting material and heat wall remain explicit.',
        primary='Odd paired response (plus-minus)/2 of Jordan temperature, material momentum and interior scalar field; report even part and unforced background drift separately.',
        decision='Can the native finite-temperature, two-heat-carrier, 26-species matter and scalar evolve on the same nonstationary nonzero background without artificial hydrostatic force subtraction? If not, fix that coupled solve before exterior or orbit integration.',
        gates=dict(native_norm=1.,scalar=2e-17,mass_work_identity=2e-15,baryon=1e-9,isotope=1e-9,
            minimum_time_order=.7,maximum_relative_last_difference=.2,minimum_signal_over_time_difference=5),
        budget=dict(workers=4,blas_threads=1,gpu=False,pilot_steps=2,pilot_timeout_seconds=90,
            production_timeout_seconds=600,maximum_steps=252,automatic_expansion=False,
            prior_measurement='Phase40: 4 small early pulse steps took 5.34 s. Nonzero-background and unprojected-force iteration counts unmeasured. Pilot must forecast full work within 600 s before continuation.'),
        scope='Actual nonlinear interior evolution of the declared finite cavity, not stationary transfer, physical atmosphere, vacuum outgoing charge, EOS certification, orbital response or full dynamic-charge completion. The saved GR fine structure is not silently substituted for this conserved coarse Cauchy state.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert e.digest(s.ROOT/rel)==digest,rel
    return p


def path(pool,label,steps,limit=None):
    p=bindings();folder=OUT/f'{label}-{steps}';folder.mkdir(exist_ok=True)
    started=time.monotonic();R=ld(json.loads((s.OUT/'initial-result.json').read_text())['radius_cm'])
    tc=R/e.C;h=2*tc/steps
    star,delta,z,boundary=initialize(pool,h);z0={k:v.copy() for k,v in z.items() if isinstance(v,np.ndarray)}
    if not (folder/'initial.npz').exists():np.savez_compressed(folder/'initial.npz',delta=delta,**z0)
    sign={'minus':-1,'undriven':0,'plus':1}[label];history=[]
    prior_seconds=0.;prior_calls=0
    if (folder/'result.json').exists():
        prior=json.loads((folder/'result.json').read_text());history=prior['history']
        prior_seconds=prior['seconds'];prior_calls=prior['native_calls']
        k=len(history);saved=dict(np.load(folder/f'step-{k:03d}.npz'))
        delta=saved.pop('delta');saved.pop('residual');z.update(saved)
        star.field_guess=tuple(z[key].copy() for key in ['psi','Pi','Phi'])
        t=ld(k)*h
        star.amplitude=boundary+(sign*ld(str(p['drive_amplitude']))*np.sin(ld(str(np.pi))*t/tc)**8 if k<steps//2 else ld(0))
    try:
        for j in range(len(history)+1,(limit or steps)+1):
            previous=(delta,z);star.previous=z;star.initial_guess=delta.copy()
            star.previous_boundary=star.amplitude
            t=ld(j)*h
            pulse=sign*ld(str(p['drive_amplitude']))*np.sin(ld(str(np.pi))*t/tc)**8 if j<steps//2 else ld(0)
            star.amplitude=boundary+pulse
            def log(row):
                row['step']=j
                with (folder/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
            delta,z,value=s.stage(star,previous,log)
            before=previous[1]
            work=z['Hgrad'][-1]*(z['bf'][-1]+before['bf'][-1])/2*star.scalar_area[-1]/(4*np.pi*e.GRAV)*(z['Phi'][-1]+before['Phi'][-1])/2*(star.amplitude-star.previous_boundary)
            dm=(z['dmf'][-1]-before['dmf'][-1])/e.GRAV
            res=np.sum(z['H']*star.volume*star.heat0*value[:,1])
            scale=np.sum(star.volume*(abs(z['Ephi'])+abs(before['Ephi'])+abs(z['dEm']-before['dEm'])))+abs(work)
            identity=float(abs(dm-work-res)/max(scale,ld('1e-100')))
            baryon=float(abs(np.sum((z['B']-z0['B'])*star.volume)/np.sum(z0['B']*star.volume)))
            inventory=float(abs(np.sum((z['dBX']-z0['dBX'])*star.volume[:,None],axis=0)/np.sum(z0['B']*star.volume)).max())
            cone=FunctionType(s.old.two.cones.__code__,dict(s.old.two.cones.__globals__,
                e=SimpleNamespace(C=e.C,TAU=z['tau_cond'])))(z)
            row=dict(step=j,t=float(t),boundary=float(star.amplitude),native_norm=float(np.max(abs(value)/s.ATOL)),
                scalar=z['scalar_residual'],identity=identity,baryon=baryon,isotope=inventory,
                cone_inside=cone['sampled_cone_inside_light_cone'])
            assert identity<p['gates']['mass_work_identity'],row
            assert baryon<1e-9 and inventory<1e-9 and row['cone_inside'],row
            history.append(row)
            np.savez_compressed(folder/f'step-{j:03d}.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
            print('NONSTATIONARY',label,steps,j,row['native_norm'],flush=True)
        result=dict(classification='Counterexample candidate',completed_steps=len(history),steps=steps,
            seconds=prior_seconds+time.monotonic()-started,native_calls=prior_calls+star.pool.evaluations,history=history)
        e.write(folder/'result.json',result);return result
    except Exception as error:
        if hasattr(star,'last_state'):np.savez_compressed(folder/'failed-iterate.npz',delta=star.last_delta,**{k:v for k,v in star.last_state.items() if isinstance(v,np.ndarray)})
        e.write(folder/'failure.json',dict(error=repr(error),seconds=time.monotonic()-started,step=j,native_calls=star.pool.evaluations));raise


def pilot():
    bindings()
    with ProcessPoolExecutor(max_workers=4,initializer=s.old.imported.initial.original.worker_init) as pool:
        row=path(pool,'plus',12,2)
    e.write(OUT/'pilot.json',dict(row,
        projected_seconds=row['seconds']/row['completed_steps']*252,
        forecast_scope='Measured only first two steps; later nonlinear cost unknown. Production retains a hard timeout.'))


def run():
    p=bindings();pilot=json.loads((OUT/'pilot.json').read_text())
    assert pilot['projected_cutoff_seconds']<p['budget']['production_timeout_seconds']
    start=time.monotonic();rows=[]
    with ProcessPoolExecutor(max_workers=4,initializer=s.old.imported.initial.original.worker_init) as pool:
        for steps in p['time_steps']:
            for label in p['paths']:
                rows.append(dict(label=label,**path(pool,label,steps,steps//2)))
    weights=np.load(s.OUT/'geometry.npz')['baryons'];weights=weights/weights.sum()
    def norm(v):return float(np.sqrt(np.sum(weights*abs(v)**2)))
    signals={key:[] for key in ['temperature','velocity','scalar']};even={key:[] for key in signals};drift=[]
    for steps in p['time_steps']:
        states={label:dict(np.load(OUT/f'{label}-{steps}'/f'step-{steps//2:03d}.npz')) for label in p['paths']}
        initial=np.load(OUT/f'plus-{steps}'/'initial.npz')
        def read(z):return dict(temperature=z['delta'][:,1]-z['logA'],velocity=z['delta'][:,2],scalar=z['psi'])
        minus,zero,plus=[read(states[label]) for label in p['paths']]
        for key in signals:
            signals[key].append((plus[key]-minus[key])/2)
            even[key].append(norm((plus[key]+minus[key])/2-zero[key]))
        iz=read(initial)
        drift.append({key:norm(zero[key]-iz[key]) for key in signals})
    comparisons=[]
    for key,values in signals.items():
        amplitudes=[norm(v) for v in values];differences=[norm(values[i+1]-values[i]) for i in [0,1]]
        order=float(np.log2(differences[0]/differences[1]));ratio=amplitudes[-1]/differences[-1]
        comparisons.append(dict(readout=key,odd_rms=amplitudes,even_rms=even[key],time_differences=differences,
            time_order=order,signal_over_last_difference=ratio,
            passed=bool(order>=p['gates']['minimum_time_order'] and ratio>=5)))
    result=dict(classification='Counterexample candidate',completed=True,passed=all(v['passed'] for v in comparisons),
        comparisons=comparisons,undriven_rms_drift=drift,paths=rows,
        total_steps=sum(v['completed_steps'] for v in rows),native_calls=sum(v['native_calls'] for v in rows),
        production_seconds=time.monotonic()-start,full_dynamic_charge_solved=False,free_surface_evolved=False,
        radiative_boundary_connected=False,stationary_transfer_solved=False,symbolic=symbolic())
    e.write(OUT/'result.json',result)
    print('FINAL',json.dumps({k:v for k,v in result.items() if k!='paths'}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','pilot','run'])
    globals()[p.parse_args().action]()
