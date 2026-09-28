"""Resume only missing native energy probes with a consistently precise inverse.

Frozen predecessors retain their original defects and budgets. This module
never repeats accepted cells or treats a source-direction probe as evolution.
"""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import time
import numpy as np
import def_reactive_paired_tangent as paired

h=paired.h
OLD=paired.OUT/'budgeted'
OUT=paired.OUT/'energy-consistent'


def inverse(query, lp, x, target, guess, budget):
    """Use the same long-double residual for iteration, selection and acceptance."""
    best=(np.longdouble('inf'),guess,None)
    seen={}
    def sample(t):
        nonlocal best
        if t not in seen:
            a=query(1,lp,t,x)
            defect=np.longdouble(a[2])+np.longdouble(a[1])/a[0]-target
            seen[t]=(defect,t,a)
            if abs(defect)<abs(best[0]):best=seen[t]
        return seen[t]
    t=guess
    for _ in range(15):
        defect,_,a=sample(t)
        if abs(best[0])<=budget:return best
        trial=float(np.longdouble(t)-defect/(a[10]-a[1]/a[0]*a[8]))
        if trial in seen:break
        t=trial
    lo=hi=best[1]
    for _ in range(32):
        lo=float(np.nextafter(lo,-np.inf));hi=float(np.nextafter(hi,np.inf))
        sample(lo);sample(hi)
        if abs(best[0])<=budget:break
    return best


def cell(i):
    d,s=paired.data,paired.sources
    X=d['X'][i];lapse=d['A'][i]*d['N'][i];T=d['lnT'][i]
    lp=float(paired.background[i,2]);calls=0
    def query(mode,p,t,x):
        nonlocal calls
        calls+=1;return paired.eos(mode,p,t,x)
    raw=query(1,lp,T,X);rho=raw[0];p=raw[1]
    cp=raw[10]-p/rho*raw[8]
    h0=np.longdouble(raw[2])+np.longdouble(p)/rho
    total0=np.longdouble((X/paired.old.thermal.g.c.A)@paired.old.thermal.g.c.W)*(h.gr.C*100)**2+np.longdouble(raw[2])
    rows=[]
    for dt in [8.,4.]:
        x=np.asarray(X.astype(np.longdouble)+np.longdouble(dt*lapse)*s['dxdt'][i],float)
        assert x.min()>=0 and abs(x.sum()-1)<1e-12
        drest=paired.rest.astype(np.longdouble)@(x.astype(np.longdouble)-X.astype(np.longdouble))
        loss=np.longdouble(dt*lapse)*(s['neutrino'][i]+s['thermal_neutrino'][i])
        target=h0-drest-loss
        budget=max(2.,32*abs(np.spacing(float(target))),float(abs(drest+loss))*1e-8)
        error,lt,a=inverse(query,lp,x,target,float(np.longdouble(T)+(-drest-loss)/cp),budget)
        assert abs(error)<=budget,(i,dt,float(error),budget)
        lr=float(np.log(np.longdouble(a[0])/rho))
        de=np.longdouble(rho)*(total0*np.expm1(np.longdouble(lr))+np.exp(np.longdouble(lr))*(drest+np.longdouble(a[2])-raw[2]))
        rows.append([lr/dt,float(de/dt),float((np.longdouble(lt)-T)/dt),float(error),budget,float(drest),float(loss)])
    return i,np.array(rows),calls


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    files=[OLD/'reused-pilot.npz',*sorted(OLD.glob('block-*.npz'))]
    records=[dict(np.load(p)) for p in files]
    ids=np.concatenate([q['indices'] for q in records]);assert len(np.unique(ids))==len(ids)
    missing=np.setdiff1d(np.arange(5735),ids)
    assert len(ids)==5483 and len(missing)==252
    inputs=[Path(__file__),Path(paired.__file__),Path(paired.old.__file__),OLD/'plan.json',*files,
            paired.old.thermal.OUT/'coefficients.npz',paired.old.thermal.OUT/'sources.npz']
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b3552307',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in inputs},
        failed_run=dict(accepted_cells=5483,missing_cells=252,reported_cell=3968,probe_seconds=8,
                       residual_erg_g=-60.31103515625,budget_erg_g=32),
        root_fix='The predecessor selected candidates using binary64 enthalpy and a rounded target but accepted using long-double sums and the original target. Use the identical precise residual everywhere; bounded adjacent-float search resolves native representability without weakening the energy gate.',
        unchanged_gates=dict(energy_absolute_floor_erg_g=2,energy_ulp=32,energy_relative=1e-8,weighted_forcing_refinement=.01),
        budget=dict(workers=8,hard_timeout_seconds=90,expected_seconds=[10,60],maximum_missing_cells=252,maximum_adjacent_pairs=32),
        scope='Reuse every successful block. Same 8/4-second directional probes, no physical time evolution, no mesh or gate expansion.'))
    rows=[]
    with ProcessPoolExecutor(max_workers=8,initializer=paired.initialize) as pool:
        futures={pool.submit(cell,int(i)):int(i) for i in missing}
        for future in as_completed(futures):
            row=future.result();rows.append(row)
            np.savez_compressed(OUT/f'cell-{row[0]}.npz',indices=[row[0]],values=[row[1]],calls=[row[2]])
            assert time.monotonic()-start<85
    records.append(dict(indices=np.array([r[0] for r in rows]),values=np.array([r[1] for r in rows]),calls=np.array([r[2] for r in rows])))
    ids=np.concatenate([q['indices'] for q in records]);order=np.argsort(ids)
    assert np.array_equal(ids[order],np.arange(5735))
    a=np.concatenate([q['values'] for q in records])[order]
    d=np.load(paired.old.thermal.OUT/'coefficients.npz');dm=d['dm'];rho=d['raw'][:,0]
    norm=np.array([dm@abs(a[:,1,0]),(dm/rho)@abs(a[:,1,1])])
    error=np.array([dm@abs(a[:,1,0]-a[:,0,0]),(dm/rho)@abs(a[:,1,1]-a[:,0,1])])
    score=error/np.maximum(norm,1e-100)
    np.savez_compressed(OUT/'forcing.npz',rows=a,rho_log_rate=a[:,1,0],energy_density_rate=a[:,1,1],logT_rate=a[:,1,2])
    result=dict(classification='Counterexample candidate',cells=5735,reused_cells=len(ids)-len(missing),
        seconds=time.monotonic()-start,new_native_calls=int(sum(r[2] for r in rows)),
        maximum_energy_inverse_score=float(np.max(abs(a[:,:,3])/a[:,:,4])),
        weighted_forcing_relative_difference=score.tolist(),forcing_gate_passed=bool(np.max(score)<.01),
        actual_composition_energy_included=True,physical_evolution=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
