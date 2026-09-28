"""Measured-budget comparison of the exponential reciprocal coupled method."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import FunctionType, SimpleNamespace
import argparse
import json
import shutil
import time
import numpy as np
import def_exponential_coupled as x

s,e,ld=x.s,x.e,x.ld
OUT=x.OUT/'production'


def prepare():
    assert not OUT.exists();pilot=json.loads((x.OUT/'result.json').read_text());assert pilot['passed']
    expected=pilot['seconds']/2*126;assert expected<600,expected
    template=json.loads((x.previous.OUT/'plan.json').read_text())
    files=[Path(__file__),Path(x.__file__),x.OUT/'plan.json',x.OUT/'result.json',x.OUT/'controls.json']
    files+=list(x.OUT.glob('*.npz'))
    OUT.mkdir()
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(s.ROOT).as_posix():e.digest(p) for p in files},
        time_steps=[12,24,48],paths=['minus','undriven','plus'],amplitude=.00001,
        duration_crossing_times=1,gates=template['gates'],
        method='Native Phase44 matter and unprojected momentum; exponential discrete-gradient scalar with autonomous prescribed-boundary energy reservoir. Full nonlinear scalar/metric iteration, no artificial matter heat correction.',
        comparison='Same initial Cauchy data, grids, pulse, cutoff readouts and gates as the failed Phase44 paths. Reuse the two successful new-method pilot steps. Do not import any old-method step.',
        budget=dict(hard_timeout_seconds=600,expected_wall_seconds=[430,590],measured_projection_seconds=expected,
            native_steps=126,reused_native_steps=2,workers=4,blas_threads=1,automatic_expansion=False),
        limitation='Finite reflecting matter/heat cavity. No physical free surface, outgoing exterior, orbital result or full completion.'))


def bindings():
    p=json.loads((OUT/'plan.json').read_text())
    for rel,digest in p['bindings'].items():assert e.digest(s.ROOT/rel)==digest,rel
    for rel,digest in json.loads((x.OUT/'plan.json').read_text())['bindings'].items():assert e.digest(s.ROOT/rel)==digest,rel
    return p


def run_path(pool,label,steps):
    plan=bindings();folder=OUT/f'{label}-{steps}';assert not folder.exists();folder.mkdir()
    started=time.monotonic();R=ld(json.loads((s.OUT/'initial-result.json').read_text())['radius_cm']);h=2*R/e.C/steps
    amplitude=ld(str(plan['amplitude']))*{'minus':-1,'undriven':0,'plus':1}[label]
    star,delta,z,boundary=x.initialize(pool,h,amplitude);initial=z.copy();history=[];prefix_calls=0;prefix_seconds=0.
    np.savez_compressed(folder/'initial.npz',delta=delta,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
    if label=='plus' and steps==12:
        pilot=json.loads((x.OUT/'result.json').read_text());history=pilot['history'];prefix_seconds=pilot['seconds'];prefix_calls=history[-1]['native_calls']
        for j in [1,2]:shutil.copy2(x.OUT/f'step-{j:03d}.npz',folder/f'step-{j:03d}.npz')
        shutil.copy2(x.OUT/'iterations.jsonl',folder/'iterations.jsonl')
        saved=dict(np.load(folder/'step-002.npz'));delta=saved.pop('delta');saved.pop('residual');z.update(saved)
        star.field_guess=[z[k].copy() for k in ['psi','Pi','Phi']]
        star.amplitude=star.boundary_basis@z['canonical'][star.n:star.n+9]
    try:
        for j in range(len(history)+1,steps//2+1):
            before=z;star.previous=z;star.previous_boundary=star.amplitude;star.initial_guess=delta.copy()
            def log(row):
                row['step']=j
                with (folder/'iterations.jsonl').open('a') as stream:stream.write(json.dumps(row)+'\n')
            delta,z,value=s.stage(star,(delta,z),log)
            dy=z['canonical']-before['canonical'];mid=(z['canonical']+before['canonical'])/2
            work=-R/e.GRAV*(star.clock_matrix@mid)@dy
            dm=(z['dmf'][-1]-before['dmf'][-1])/e.GRAV
            res=np.sum(z['H']*star.volume*star.heat0*value[:,1])
            scale=np.sum(star.volume*(abs(z['Ephi'])+abs(before['Ephi'])+abs(z['dEm']-before['dEm'])))+abs(work)
            identity=float(abs(dm-work-res)/max(scale,ld('1e-100')))
            baryon=float(abs(np.sum((z['B']-initial['B'])*star.volume)/np.sum(initial['B']*star.volume)))
            isotope=float(abs(np.sum((z['dBX']-initial['dBX'])*star.volume[:,None],axis=0)/np.sum(initial['B']*star.volume)).max())
            cone=FunctionType(s.old.two.cones.__code__,dict(s.old.two.cones.__globals__,e=SimpleNamespace(C=e.C,TAU=z['tau_cond'])))(z)
            expected=boundary+amplitude*np.sin(ld(str(np.pi))*2*ld(j)/steps)**8
            row=dict(step=j,native_norm=float(np.max(abs(value)/s.ATOL)),scalar=z['scalar_residual'],
                mass_work_identity=identity,baryon=baryon,isotope=isotope,
                boundary_error=float(abs(star.amplitude-expected)),clock_work_erg=float(work),
                cone_speed=cone['maximum_local_rest_characteristic_speed_over_c'])
            np.savez_compressed(folder/f'step-{j:03d}.npz',delta=delta,residual=value,**{k:v for k,v in z.items() if isinstance(v,np.ndarray)})
            history.append(row)
            assert identity<2e-15 and baryon<1e-9 and isotope<1e-9 and cone['sampled_cone_inside_light_cone'] and row['boundary_error']<2e-19,row
            e.write(OUT/'progress.json',dict(label=label,steps=steps,completed=j,**row))
        row=dict(label=label,steps=steps,completed_steps=len(history),history=history,
            seconds=prefix_seconds+time.monotonic()-started,native_calls=prefix_calls+star.pool.evaluations)
        e.write(folder/'result.json',row);print('COMPLETED',label,steps,row['seconds'],flush=True);return row
    except Exception as error:
        e.write(folder/'failure.json',dict(error=repr(error),history=history,seconds=time.monotonic()-started));raise


def compare(paths):
    plan=bindings();weights=np.load(s.OUT/'geometry.npz')['baryons'];weights/=weights.sum()
    norm=lambda v:float(np.sqrt(np.sum(weights*v*v)))
    signals={key:[] for key in ['temperature','velocity','scalar']};evens={key:[] for key in signals}
    for steps in plan['time_steps']:
        readings=[]
        for label in plan['paths']:
            z=np.load(OUT/f'{label}-{steps}'/f'step-{steps//2:03d}.npz')
            readings.append(dict(temperature=z['delta'][:,1]-z['logA'],velocity=z['delta'][:,2],scalar=z['psi']))
        minus,zero,plus=readings
        for key in signals:
            signals[key].append((plus[key]-minus[key])/2)
            evens[key].append(norm((plus[key]+minus[key])/2-zero[key]))
    rows=[]
    for key,values in signals.items():
        differences=[norm(values[i+1]-values[i]) for i in [0,1]]
        order=float(np.log2(differences[0]/differences[1]));amplitude=norm(values[-1]);ratio=amplitude/differences[-1]
        rows.append(dict(readout=key,odd_rms=list(map(norm,values)),even_rms=evens[key],time_differences=differences,
            time_order=order,signal_over_last_difference=ratio,passed=bool(order>=.7 and ratio>=5)))
    return dict(classification='Counterexample candidate',completed=True,passed=all(r['passed'] for r in rows),
        comparisons=rows,paths=paths,total_steps=sum(r['completed_steps'] for r in paths),
        native_calls=sum(r['native_calls'] for r in paths),full_dynamic_charge_solved=False)


def run():
    plan=bindings();assert not (OUT/'result.json').exists() and not (OUT/'failure.json').exists()
    started=time.monotonic();paths=[]
    try:
        with ProcessPoolExecutor(max_workers=4,initializer=s.old.imported.initial.original.worker_init) as pool:
            for steps in plan['time_steps']:
                for label in plan['paths']:paths.append(run_path(pool,label,steps))
        result=compare(paths);result['seconds']=time.monotonic()-started;e.write(OUT/'result.json',result)
        print('FINAL',json.dumps({k:v for k,v in result.items() if k!='paths'}),flush=True)
    except Exception as error:
        e.write(OUT/'failure.json',dict(error=repr(error),completed_paths=paths,seconds=time.monotonic()-started));raise


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
