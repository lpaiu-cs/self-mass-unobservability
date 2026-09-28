"""Native directional thermal/composition response for the new background.

The two durations are differencing probes along the initial reaction vector,
not a frozen-rate stellar evolution. Native enthalpy includes chemical changes.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import json
import time
import numpy as np
import def_free_surface_thermal as thermal

h=thermal.h
OUT=thermal.OUT/'reactive-tangent'


def initialize():
    global data,sources,eos,rest
    data=dict(np.load(thermal.OUT/'coefficients.npz'));sources=dict(np.load(thermal.OUT/'sources.npz'))
    eos=h.molecular.model.EOS();rest=(thermal.g.c.W/thermal.g.c.A-1)*(h.gr.C*100)**2


def cell(i):
    raw=data['raw'][i];X=data['X'][i];lapse=data['A'][i]*data['N'][i]
    cp=data['thermo'][i,5];p=raw[1];rho=raw[0];T=data['lnT'][i]
    h0=np.longdouble(raw[2])+np.longdouble(p)/rho;rows=[];calls=0
    total0=(np.longdouble((X/thermal.g.c.A)@thermal.g.c.W)*(h.gr.C*100)**2+np.longdouble(raw[2]))
    for dt in [8.,4.]:
        x=np.asarray(X.astype(np.longdouble)+np.longdouble(dt*lapse)*sources['dxdt'][i],float)
        assert x.min()>=0 and abs(x.sum()-1)<1e-12,(i,x.min())
        drest=rest.astype(np.longdouble)@(x.astype(np.longdouble)-X.astype(np.longdouble))
        loss=np.longdouble(dt*lapse)*(sources['neutrino'][i]+sources['thermal_neutrino'][i])
        target=h0-drest-loss;lt=float(T+(-drest-loss)/cp)
        budget=max(2.,32*abs(np.spacing(float(target))),float(abs(drest+loss))*1e-8)
        for iteration in range(5):
            a=eos(1,float(np.log(p)),lt,x);calls+=1
            enthalpy=np.longdouble(a[2])+np.longdouble(a[1])/a[0]
            error=enthalpy-target
            if abs(error)<=budget:break
            capacity=a[10]-a[1]/a[0]*a[8]
            lt-=float(error/capacity)
        else:raise RuntimeError(('reactive enthalpy inverse',i,dt,float(error),budget))
        lr=float(np.log(np.longdouble(a[0])/rho))
        de=np.longdouble(rho)*(total0*np.expm1(np.longdouble(lr))+np.exp(np.longdouble(lr))*(drest+np.longdouble(a[2])-raw[2]))
        rows.append([lr/dt,float(de/dt),float((np.longdouble(lt)-T)/dt),float(error),budget,float(drest),float(loss)])
    return i,np.array(rows),calls


def block(indices):return [cell(int(i)) for i in indices]


def run():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    files=[Path(__file__),Path(thermal.__file__),thermal.OUT/'coefficients.npz',thermal.OUT/'sources.npz']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in files},
        claim='Compute the coupled thermal/composition EOS forcing at fixed instantaneous pressure, conserving atomic rest plus native enthalpy minus actual neutrino loss. Includes u(X) without 26 separate EOS Jacobian directions.',
        probe_coordinate_seconds=[8.,4.],
        scope='Directional probes along the initial native reaction vector. They are not finite physical evolution; radiation/heat-flux divergence vanishes initially because the specified initial flux is zero. Transport starts through its flux time derivative.',
        gates=dict(energy_absolute_floor_erg_g=2,energy_ulp=32,energy_relative=1e-8,weighted_forcing_refinement=.01),
        budget=dict(workers=8,pilot_cells=32,hard_timeout_seconds=200,maximum_EOS_calls=57350,automatic_expansion=False)))
    indices=np.unique(np.linspace(0,5734,32).astype(int));rows=[]
    with ProcessPoolExecutor(max_workers=8,initializer=initialize) as pool:
        began=time.monotonic()
        for part in pool.map(block,np.array_split(indices,8)):rows.extend(part)
        seconds=time.monotonic()-began
        # Exclude no costs: this conservative forecast includes pool setup.
        forecast=seconds*5735/len(indices)
        h.write(OUT/'pilot.json',dict(classification='Counterexample candidate',seconds=seconds,projection_seconds=forecast,within_budget=forecast<190))
        np.savez_compressed(OUT/'pilot.npz',indices=[r[0] for r in rows],values=[r[1] for r in rows],calls=[r[2] for r in rows])
        assert forecast<190,('Probe budget',forecast)
        remaining=np.setdiff1d(np.arange(5735),indices)
        for part in pool.map(block,np.array_split(remaining,90)):
            rows.extend(part)
            assert time.monotonic()-start<195
    rows.sort(key=lambda r:r[0]);a=np.array([r[1] for r in rows]);calls=sum(r[2] for r in rows)
    d=np.load(thermal.OUT/'coefficients.npz');dm=d['dm'];rho=d['raw'][:,0]
    # Mass-weighted density and energy-source contrasts; report absolute floors
    # instead of demanding relative resolution where the native change is zero.
    norm=np.array([dm@abs(a[:,1,0]),(dm/rho)@abs(a[:,1,1])])
    error=np.array([dm@abs(a[:,1,0]-a[:,0,0]),(dm/rho)@abs(a[:,1,1]-a[:,0,1])])
    score=error/np.maximum(norm,1e-100)
    np.savez_compressed(OUT/'forcing.npz',rows=a,rho_log_rate=a[:,1,0],energy_density_rate=a[:,1,1],logT_rate=a[:,1,2])
    result=dict(classification='Counterexample candidate',cells=len(rows),native_EOS_calls=calls,seconds=time.monotonic()-start,
        maximum_energy_inverse_score=float(np.max(abs(a[:,:,3])/a[:,:,4])),weighted_forcing_relative_difference=score.tolist(),
        forcing_gate_passed=bool(np.max(score)<.01),maximum_logT_rate=float(abs(a[:,1,2]).max()),
        actual_composition_energy_included=True,physical_evolution=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':run()
