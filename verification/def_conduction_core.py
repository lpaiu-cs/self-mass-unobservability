"""Microscopic core pole bank on the saved stellar states.

Convex interpolation preserves the measured positive response. Withheld
states test interpolation; this does not certify omitted collision physics.
"""
from pathlib import Path
import argparse
import json
import time
import resource
import numpy as np
from scipy.interpolate import CubicSpline
import def_electron_screened_fermi as microscopic

ex=microscopic.exchange
screen=microscopic.screen
h=ex.h
OUT=ex.OUT.parent/'def-conduction-core'
ANCHORS=np.linspace(4122,5734,5).astype(int)
WITHHELD=((ANCHORS[:-1]+ANCHORS[1:])//2).astype(int)


def cell(index):
    state=ex.equilibrium(index);maximum=np.sqrt(1-1/(1+2*state['theta']*(state['eta']+40))**2)*(1+1e-12)
    grid=microscopic.nodes(state,maximum)
    values=np.array([screen.polarization(a,state) for a in grid])
    vf=state['xF']/np.sqrt(1+state['xF']**2)
    mids=(grid[1:]+grid[:-1])/2
    selected=np.unique(np.r_[np.linspace(0,255,16).astype(int),np.argsort(abs(mids-vf))[:16]])
    reference=np.array([screen.polarization(a,state) for a in mids[selected]])
    table_error=float(np.max(abs(reference-CubicSpline(grid,values)(mids[selected])))/state['qs2'])
    ex.amplitude=microscopic.leading.amplitude(state,grid,values)
    fine=[];coarse=[]
    for seed in ex.SEEDS:
        a,b,_=ex.bracket(state,17,seed,15);fine.append(a);coarse.append(b)
    EE=np.mean(fine,axis=0);answer=ex.transfer(state,EE,7)
    small=ex.transfer(state,np.mean(coarse,axis=0),7);low=ex.transfer(state,EE,5)
    errors={k:dict(quadrature=float(abs(small[k]/answer[k]-1)),basis=float(abs(low[k]/answer[k]-1))) for k in ['K','tau']}
    np.savez_compressed(OUT/f'cell-{index}.npz',**state,EE=EE,replicate_EE=np.array(fine),phase=grid,polarization=values,**{'fine_'+k:v for k,v in answer.items()})
    return dict(cell=int(index),errors=errors,table_error=table_error,native_ratio=answer['K']/state['native_K'])


def bank(points,anchors):
    d=np.load(ex.model.base.thermal.OUT/'coefficients.npz');native=d['K'][:,1]*1e-5
    rows=[np.load(OUT/f'cell-{i}.npz') for i in anchors]
    a=np.clip(np.searchsorted(anchors,points)-1,0,len(anchors)-2);b=a+1
    f=(points-anchors[a])/(anchors[b]-anchors[a]);assert np.all((f>=0)&(f<=1))
    ratio=np.array([q['fine_K']/q['native_K'] for q in rows])
    K=np.exp(np.interp(points,np.arange(len(native)),np.log(native))+(1-f)*np.log(ratio[a])+f*np.log(ratio[b]))
    poles=np.array([q['fine_poles'] for q in rows]);fraction=np.array([q['fine_weights']/q['fine_poles'] for q in rows]);fraction/=fraction.sum(1)[:,None]
    rates=np.c_[poles[a],poles[b]]
    weights=np.c_[(1-f)[:,None]*fraction[a],f[:,None]*fraction[b]]*K[:,None]
    assert np.min(weights)>=0 and np.min(rates)>0
    return weights,rates


def prepare():
    assert not OUT.exists();OUT.mkdir()
    paths=[Path(__file__),Path(microscopic.__file__),Path(ex.__file__),Path(screen.__file__),Path(microscopic.leading.__file__),
           ex.model.base.thermal.OUT/'coefficients.npz',microscopic.OUT/'result.json']
    paths+=list(microscopic.OUT.glob('cell-*.npz'))
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='cca3717e',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths},
        claim='Construct the state-dependent collision pole bank required by the actual core conduction contribution to the existing free-surface fluid/scalar/metric time response.',
        anchors=ANCHORS.tolist(),withheld=WITHHELD.tolist(),reuse=[4122,4928,5734],
        interpolation='Interpolate the computed/native DC ratio logarithmically and mix the adjacent positive normalized DC pole weights. Native conductivity only supplies the interpolation coordinate trend; there is no fitted normalization. Withheld states are not added to the bank.',
        domain='Original degenerate, almost fully ionized cells 4122..5734 only. Faces joining two core cells define an isolated core transport contribution. Missing outer heat currents are not physically certified zero.',
        gates=dict(quadrature=.03,basis=.02,polarization=1e-5,native_compatibility=.2,withheld_K=.01,withheld_tau=.01),
        budget=dict(new_states=6,events_per_state=4*2**17,hard_seconds=120,workers=1,BLAS_threads=1,native_calls=0,automatic_expansion=False)))
    for i in [4122,4928,5734]:(OUT/f'cell-{i}.npz').write_bytes((microscopic.OUT/f'cell-{i}.npz').read_bytes())
    start=time.monotonic();s=ex.equilibrium(4928);d=np.load(microscopic.OUT/'polarization-4928.npz')
    ex.amplitude=microscopic.leading.amplitude(s,d['phase'],d['polarization']);ex.bracket(s,12,ex.SEEDS[0]);pilot=time.monotonic()-start
    ex.write(OUT/'pilot.json',dict(seconds=pilot,forecast_seconds=pilot*128*10+20,
        assumption='Linear event scaling plus 20 seconds for ten polarization tables, based on prior 5.71 seconds for three states. Stop at 120 seconds; no automatic expansion.'))
    print('PILOT',pilot,'FORECAST',pilot*128*10+20,flush=True)


def run():
    plan=json.loads((OUT/'plan.json').read_text());assert not (OUT/'result.json').exists()
    for p,s in plan['bindings'].items():assert h.digest(h.ROOT/p)==s,p
    assert plan['revised_forecast_seconds']<110
    start=time.monotonic();checks=[]
    for i in sorted(set(plan['anchors']+plan['withheld'])-set(plan['reuse'])):
        row=cell(i);checks.append(row);ex.write(OUT/'progress.json',checks);print('CELL',i,row,flush=True)
        assert time.monotonic()-start<115
    comparisons=[]
    for i in WITHHELD:
        w,r=bank(np.array([i]),ANCHORS);direct=np.load(OUT/f'cell-{i}.npz')
        comparisons.append(dict(cell=int(i),K=float(abs(w.sum()/direct['fine_K']-1)),tau=float(abs(np.sum(w/r)/w.sum()/direct['fine_tau']-1))))
    faces=np.arange(4123,5735);weights,rates=bank(faces-.5,ANCHORS)
    np.savez_compressed(OUT/'bank.npz',faces=faces,mode_K_SI=weights,poles_proper_s=rates)
    coarse_w,coarse_r=bank(faces-.5,ANCHORS[::2]);np.savez_compressed(OUT/'coarse-bank.npz',faces=faces,mode_K_SI=coarse_w,poles_proper_s=coarse_r)
    passed=all(x['table_error']<1e-5 and abs(x['native_ratio']-1)<.2 and all(e['quadrature']<.03 and e['basis']<.02 for e in x['errors'].values()) for x in checks)
    passed=passed and all(x['K']<.01 and x['tau']<.01 for x in comparisons)
    result=dict(classification='Counterexample candidate',passed=passed,checks=checks,withheld=comparisons,
        seconds=time.monotonic()-start,peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        covered_cells=1613,covered_internal_faces=1612,physical_error_certified=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',passed,'WITHHELD',comparisons,'SECONDS',result['seconds'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run']);globals()[p.parse_args().action]()
