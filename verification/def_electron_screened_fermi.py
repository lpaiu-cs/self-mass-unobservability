"""Resolve the Fermi-velocity polarization feature at fixed table size."""
from pathlib import Path
import json
import time
import numpy as np
from scipy.interpolate import CubicSpline
import def_electron_screened_leading as leading

screen=leading.screen
exchange=leading.exchange
OUT=exchange.OUT/'fermi-screening'


def nodes(state,maximum):
    vf=state['xF']/np.sqrt(1+state['xF']**2)
    width=max(.03,3*state['theta']/(vf*(1+state['xF']**2)**1.5))
    lo=np.arctan(-vf/width);span=np.arctan((maximum-vf)/width)-lo
    target=np.linspace(0,1,257);left=np.zeros(257);right=np.full(257,maximum)
    for _ in range(45):
        mid=(left+right)/2;cdf=.5*mid/maximum+.5*(np.arctan((mid-vf)/width)-lo)/span
        left=np.where(cdf<target,mid,left);right=np.where(cdf<target,right,mid)
    result=(left+right)/2;result[0]=0;result[-1]=maximum
    return result


def main():
    assert not OUT.exists();OUT.mkdir();h=screen.h
    previous=leading.OUT;plan=json.loads((previous/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    paths=[Path(__file__),Path(leading.__file__),previous/'plan.json',previous/'preparation.json']
    plan['bindings'].update({p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths})
    plan['preserved_interpolation_failure']='Original uniform 257-node table failed the unchanged 1e-5 withheld criterion at core cell 5734 (1.52718567038592e-5). No leading-screening production ran.'
    plan['table_change']='Keep 257 nodes and the same integration formula; distribute them with a 50/50 uniform and Fermi-centered Cauchy CDF. Width max(0.03,3*T/(v_F*E_F^3)). Do not increase event counts or change physical assumptions or gates.'
    exchange.write(OUT/'plan.json',plan)
    records=[];start=time.monotonic()
    for index in plan['cells']:
        state=exchange.equilibrium(index);old=np.load(previous/f'polarization-{index}.npz')
        grid=nodes(state,float(old['phase'][-1]));values=np.array([screen.polarization(a,state) for a in grid])
        np.savez_compressed(OUT/f'polarization-{index}.npz',phase=grid,polarization=values)
        sample=(np.linspace(0,len(grid)-2,16).astype(int)+.5)*grid[-1]/(len(grid)-1)
        # Also inspect the new nonuniform panel midpoints near the Fermi level.
        vf=state['xF']/np.sqrt(1+state['xF']**2)
        near=np.argsort(abs((grid[1:]+grid[:-1])/2-vf))[:16]
        sample=np.r_[sample,(grid[near]+grid[near+1])/2]
        exact=np.array([screen.polarization(a,state) for a in sample]);estimate=CubicSpline(grid,values)(sample)
        direct=exchange.model.kinetic(np.array([index]),256)['conductivity_SI'][0]
        projected=exchange.transfer(state,np.zeros_like(state['G']),7)['K']
        records.append(dict(cell=index,relative_static_compressibility=float(abs(values[0,0]/state['qs2']-1)),withheld_interpolation_relative=float(np.max(abs(exact-estimate))/state['qs2']),EI_direct_SI_comparison=float(abs(projected/direct-1))))
    state=exchange.equilibrium(plan['cells'][1]);d=np.load(OUT/f"polarization-{state['index']}.npz")
    amp=leading.amplitude(state,d['phase'],d['polarization']);vertex=leading.recovery.controls(state,amp)
    exchange.amplitude=amp;began=time.monotonic();exchange.bracket(state,12,exchange.SEEDS[0]);elapsed=time.monotonic()-began
    exchange.write(OUT/'preparation.json',dict(records=records,vertex_checks=vertex,pilot_seconds=elapsed,
        production_linear_forecast_seconds=elapsed*384,table_seconds=began-start,estimate='Same event count and 120s process cap. 257 table nodes redistributed, not increased.'))
    print('PREPARATION',json.dumps(records),'FORECAST',elapsed*384,flush=True)
    screen.OUT=OUT;screen.make_amplitude=leading.amplitude;screen.run()


if __name__=='__main__':main()
