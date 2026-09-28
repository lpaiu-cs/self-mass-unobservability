"""All-cell native EOS audit of the finite-pressure mechanical background."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
import time
import numpy as np
import def_hydrostatic_background as h

OUT=h.OUT/'absolute-shoot'


def finalize_saved():
    """Finish metadata after the accepted solve's parent-table path error."""
    if (OUT/'result.json').exists():return
    start=time.monotonic();plan=json.loads((OUT/'plan.json').read_text())
    for rel,digest in plan['bindings'].items():assert h.digest(h.ROOT/rel)==digest,rel
    saved=np.load(OUT/'background-0.001.npz');x=saved['parameters']
    progress=json.loads((OUT/'progress-0.001.json').read_text())
    assert np.array_equal(x,progress['history'][-1]['parameters'])
    s=h.Structure(.001);error,_,_,ys=s.branches(x)
    assert abs(error).max()<plan['gates']['interface']
    assert np.array_equal(ys,saved['faces'][0])
    R=ys[0]*s.R;M=ys[1]*s.B;q=.001*s.mu*x[4]
    ext=h.exterior.exact(h.mp.mpf(M/R),h.mp.mpf(q));ADM=M+R*q*q*float(ext[0])
    row=dict(classification='Counterexample candidate',phi_infinity=.001,parameters=x.tolist(),
        interface=error.tolist(),objective_calls=len(progress['history']),radius_m=R,mass_geom_m=M,ADM_geom_m=ADM,
        alpha_over_phi_infinity=s.mu*x[4]*R*float(ext[1])/ADM,
        table_sha256=h.digest(h.OUT/'extended-table.npz'),fixed_inventory=True,
        thermal_stationarity=False,physical_atmosphere=False)
    h.write(OUT/'metadata-recovery.json',dict(classification='Counterexample candidate',
        reason='Accepted background was saved before the report tried to hash extended-table.npz in the child output instead of its unchanged parent.',
        source_sha256=h.digest(Path(__file__)),saved_background_sha256=h.digest(OUT/'background-0.001.npz'),
        new_native_calls=0,new_background_fits=0,replayed_interface=error.tolist(),seconds=time.monotonic()-start))
    h.write(OUT/'result-0.001.json',row)
    h.write(OUT/'result.json',dict(classification='Counterexample candidate',mechanical_matching_passed=True,
        rows=[row],shoot_progress_seconds=progress['seconds'],original_runner_wall_seconds=40.19,
        original_runner_exit_status=1,metadata_recovery_completed=True,native_audit_passed=False,full_dynamic_charge_solved=False))


def initialize():
    global eos,data,table,background
    eos=h.molecular.model.EOS();data,table=h.inputs()
    background=dict(np.load(OUT/'background-0.001.npz'))


def evaluate(indices):
    rows=[];raw=[]
    for i in indices:
        lp=background['states'][i,2];lr,lt,u=background['thermo'][i]
        a=eos(1,float(lp),float(lt),data['X'][i]);raw.append(a)
        rows.append([i,float(np.log(a[0])-lr),float((a[3]-table['reference'][i,3])*np.exp(lt)/a[10]),
                     float((a[2]*1e-4/h.gr.C**2-u)/(data['CX'][i]+u))])
    return np.asarray(rows),np.asarray(raw)


def run():
    target=OUT/'native-audit.json';assert not target.exists()
    finalize_saved()
    plan=json.loads((OUT/'plan.json').read_text());result=json.loads((OUT/'result.json').read_text())
    assert result['mechanical_matching_passed']
    data,table=h.inputs();n=len(data['dm']);b=dict(np.load(OUT/'background-0.001.npz'))
    assert np.array_equal(b['dm'],data['dm']) and np.array_equal(b['X'],data['X'])
    assert np.array_equal(b['entropy'],table['reference'][:,3])
    files=[Path(__file__),OUT/'background-0.001.npz',OUT/'result.json',OUT/'plan.json']
    h.write(OUT/'native-audit-plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in files},
        claim='Evaluate the original native finite-temperature EOS independently at every reconstructed midpoint. Retain the same nuclear and entropy inventories; quantify interpolation disagreement.',
        gates=plan['gates'],budget=dict(workers=4,hard_timeout_seconds=180,maximum_native_calls=n,automatic_expansion=False),
        scope='All saved midpoints, not continuous EOS/space errors, thermal stationarity, atmosphere or dynamic response.'))
    start=time.monotonic();selected=np.unique(np.linspace(0,n-1,32).astype(int));rows=[];raw=[]
    with ProcessPoolExecutor(max_workers=4,initializer=initialize) as pool:
        pilot=time.monotonic()
        for r,a in pool.map(evaluate,np.array_split(selected,4)):rows.extend(r);raw.extend(a)
        seconds=time.monotonic()-pilot;projection=seconds*n/len(selected)
        h.write(OUT/'native-audit-budget.json',dict(classification='Counterexample candidate',pilot_seconds=seconds,
            conservative_projection_seconds=projection,hard_timeout_seconds=180,
            basis='32 distributed cells, includes worker initialization; other states unmeasured. Pilot outputs are reused.'))
        assert projection<170,('Native audit exceeds measured budget',projection)
        remaining=np.setdiff1d(np.arange(n),selected)
        for r,a in pool.map(evaluate,np.array_split(remaining,int(np.ceil(len(remaining)/128)))):
            rows.extend(r);raw.extend(a)
            assert time.monotonic()-start<175,'Native audit wall budget'
    rows=np.asarray(rows);raw=np.asarray(raw);order=np.argsort(rows[:,0]);rows,raw=rows[order],raw[order]
    assert np.array_equal(rows[:,0],np.arange(n))
    maxima=np.max(abs(rows[:,1:]),axis=0)
    passed=maxima[0]<plan['gates']['native_logrho'] and maxima[1]<plan['gates']['native_entropy_over_cv']
    np.savez_compressed(OUT/'native-audit.npz',rows=rows,raw=raw)
    value=dict(classification='Counterexample candidate',passed=bool(passed),cells=n,
        maximum_logrho_error=float(maxima[0]),maximum_entropy_over_cv=float(maxima[1]),
        maximum_energy_relative=float(maxima[2]),native_calls=n,seconds=time.monotonic()-start,
        fixed_baryon_isotope_entropy_inventory=True,physical_EOS_certified=False,thermal_stationarity=False,
        physical_atmosphere=False,continuous_errors_certified=False,full_dynamic_charge_solved=False)
    h.write(target,value);print(json.dumps(value),flush=True)


if __name__=='__main__':run()
