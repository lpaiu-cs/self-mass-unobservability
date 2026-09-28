"""Budgeted full native chemical-energy direction, with durable block reuse."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json
import time
import numpy as np
import def_reactive_direct_gate as direct

old=direct.paired.old
h=old.h
OUT=direct.ns['OUT'].parent/'budgeted'


def block(indices):
    rows=direct.block(indices)
    path=OUT/f'block-{int(indices[0]):04d}.npz'
    np.savez_compressed(path,indices=[r[0] for r in rows],values=[r[1] for r in rows],calls=[r[2] for r in rows])
    h.write(path.with_suffix('.json'),dict(plan_sha256=h.digest(OUT/'plan.json'),output_sha256=h.digest(path),cells=len(rows)))
    return int(indices[0]),len(rows)


def main():
    assert not OUT.exists();OUT.mkdir();start=time.monotonic()
    pilot=OUT.parent/'warm-budget/pilot.npz';measurement=json.loads((pilot.parent/'pilot.json').read_text())
    files=[Path(__file__),Path(direct.__file__),Path(direct.paired.__file__),Path(old.__file__),
        Path(old.thermal.g.s.__file__),old.thermal.OUT/'coefficients.npz',old.thermal.OUT/'sources.npz',pilot]
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='b3552307',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in files},
        purpose='Supply the missing actual composition/internal-energy forcing to the same-inventory free-surface mass/scalar response. A rest-mass-only nuclear heating shortcut omits this term.',
        alternatives='Reuse the full native background and source values; analytic coordinate transform already avoids new thermal derivative evaluations. Two actual reaction directions replace 26 separate composition columns. Old-model chemical derivatives cannot replace the current molecular EOS.',
        budget_reassessment='Preserve the failed 200-second forecasts. The unavoidable paired native chemical calls remain the bottleneck after reusing saved inputs and measuring warm workers. Allocate one bounded run, with no grid, frequency, duration or acceptance expansion.',
        budget=dict(workers=8,BLAS_threads=1,GPU=False,expected_wall_seconds=[340,510],hard_timeout_seconds=600,maximum_runs=1,save_every_block=True),
        measurement=measurement,probe_coordinate_seconds=[8.,4.],
        gates=dict(energy_absolute_floor_erg_g=2,energy_ulp=32,energy_relative=1e-8,weighted_forcing_refinement=.01),
        scope='Initial-source directional EOS probes, not a physical eight-second evolution or a complete nuclear Jacobian. Transport and full time-dependent exterior remain separate.',
        stop='Native energy/input failure or 590-second internal deadline; preserve completed blocks. No automatic restart or extra probes.'))
    p=np.load(pilot);selected=p['indices'].astype(int)
    np.savez_compressed(OUT/'reused-pilot.npz',**dict(p))
    remaining=np.setdiff1d(np.arange(5735),selected);chunks=np.array_split(remaining,90)
    complete=len(selected)
    with ProcessPoolExecutor(max_workers=8,initializer=direct.initialize) as pool:
        futures=[pool.submit(block,c) for c in chunks]
        for future in as_completed(futures):
            _,count=future.result();complete+=count
            h.write(OUT/'progress.json',dict(completed_cells=complete,total_cells=5735,seconds=time.monotonic()-start))
            assert time.monotonic()-start<590
    records=[np.load(OUT/'reused-pilot.npz')]+[np.load(OUT/f'block-{int(c[0]):04d}.npz') for c in chunks]
    ids=np.concatenate([p['indices'] for p in records]);order=np.argsort(ids)
    assert np.array_equal(ids[order],np.arange(5735))
    a=np.concatenate([p['values'] for p in records])[order];calls=int(sum(np.sum(p['calls']) for p in records))
    d=np.load(old.thermal.OUT/'coefficients.npz');dm=d['dm'];rho=d['raw'][:,0]
    norm=np.array([dm@abs(a[:,1,0]),(dm/rho)@abs(a[:,1,1])])
    error=np.array([dm@abs(a[:,1,0]-a[:,0,0]),(dm/rho)@abs(a[:,1,1]-a[:,0,1])]);score=error/np.maximum(norm,1e-100)
    np.savez_compressed(OUT/'forcing.npz',rows=a,rho_log_rate=a[:,1,0],energy_density_rate=a[:,1,1],logT_rate=a[:,1,2])
    result=dict(classification='Counterexample candidate',cells=5735,native_EOS_calls_including_reused_pilot=calls,
        seconds=time.monotonic()-start,maximum_energy_inverse_score=float(np.max(abs(a[:,:,3])/a[:,:,4])),
        weighted_forcing_relative_difference=score.tolist(),forcing_gate_passed=bool(np.max(score)<.01),
        maximum_logT_rate=float(abs(a[:,1,2]).max()),actual_composition_energy_included=True,
        physical_evolution=False,full_dynamic_charge_solved=False)
    h.write(OUT/'result.json',result);print(json.dumps(result),flush=True)


if __name__=='__main__':main()
