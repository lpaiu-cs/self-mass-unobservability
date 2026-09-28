"""Counterexample candidate: isolate the frozen step-42 nonlinear failure.

Compare the existing sparse iteration matrix, its full tangent action and the
native residual at the same saved iterate. Only composition-response refresh
differs between the two arms; equations, histories and time step stay fixed.
"""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json

import numpy as np
import gr_compatible_full_duration as full

m, e, ld = full.cached, full.e, full.ld
OUT = m.BASE/'step42-analysis'


def failed_state(star):
    with np.load(full.OUT/'path-1/failed-iterate.npz') as cp:
        delta = cp['delta'].copy()
        y = star.base+delta
        star.material_cache.clear()
        star.material_cache.update({e.material_key(row):aux.copy() for row,aux in
            zip(zip(y[:,0],y[:,1],y[:,5:]),cp['aux'])})
        return delta, star.evaluate(delta)


def run():
    assert not OUT.exists()
    plan = full.bindings()
    OUT.mkdir()
    inputs = [Path(__file__),full.OUT/'plan.json',full.OUT/'path-1/failure.json',
              *[full.OUT/'path-1'/f for f in ['step-0040.npz','step-0041.npz','failed-iterate.npz']]]
    e.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        question='Separate sparse tangent truncation, lagged composition response and native-versus-tangent mismatch at the identical failed stage.',
        bindings={p.relative_to(e.ROOT).as_posix():e.digest(p) for p in inputs}))
    with ProcessPoolExecutor(max_workers=15,initializer=e.worker_init) as pool:
        star = m.initialize(pool)
        times = full.time_nodes(plan,1)
        older, previous = [full.restored(star,full.OUT/'path-1',k,times[k]) for k in (40,41)]
        delta,z = failed_state(star)
        h = times[42]-times[41]
        coefficients = m.weights(h,times[41]-times[40])
        value,z = m.residual(star,delta,previous,older,h,coefficients)
        norm = float(np.max(abs(value)/m.ATOL))
        assert norm == 18.768935482897138,norm
        cell,field = np.unravel_index(np.argmax(abs(value)/m.ATOL),value.shape)
        star.composition_context = (previous,older,h,coefficients)
        def trial(change):
            candidate = delta.copy()
            candidate[:,:5] += change
            candidate[:,5:] = m.species(star,candidate,previous,older,h,coefficients,z)
            return m.residual(star,candidate,previous,older,h,coefficients)
        zero,zero_z = trial(np.zeros_like(delta[:,:5]))
        report = dict(classification='Counterexample candidate',initial_score=norm,
            cell=int(cell),field=int(field),rho=float(z['rho'][cell]),temperature=float(z['T'][cell]),
            step_seconds=str(h),zero_correction_species_projection_score=float(np.max(abs(zero)/m.ATOL)),
            arms=[])
        print('FAILED ITERATE REPRODUCED',json.dumps(report),flush=True)
        for name,anchor in [('lagged_composition',previous),('refreshed_composition',(delta,z))]:
            star.composition_data = m.composition_response(star,*anchor)
            tangent = m.tangent(star,delta,z)
            matrix = m.jacobian(tangent,delta,previous,older,h,coefficients)
            factor = m.splu(matrix)
            f = np.asarray(value[:,:5]/m.SCALE,float).ravel()
            direction = -factor.solve(f)
            change = direction.reshape(star.n,5).astype(ld)*m.SCALE
            predicted = (matrix@direction).reshape(star.n,5)*m.SCALE
            probe = np.zeros_like(delta)
            probe[:,:5] = change
            plus = m.residual(tangent,delta+probe,previous,older,h,coefficients)[0]
            minus = m.residual(tangent,delta-probe,previous,older,h,coefficients)[0]
            full_action = (plus-minus)/2
            native_plus,_ = trial(change)
            native_minus,_ = trial(-change)
            native_action = (native_plus-native_minus)/2
            tiny,_ = trial(change/128)
            defect = full_action[:,:5]-predicted
            native_defect = native_action[:,:5]-full_action[:,:5]
            item = dict(arm=name,
                linear_solve_defect=float(np.max(abs(matrix@direction+f))),
                sparse_vs_full_tangent_maximum=(np.max(abs(defect),axis=0)/m.ATOL[:5]).astype(float).tolist(),
                tangent_vs_native_maximum=(np.max(abs(native_defect),axis=0)/m.ATOL[:5]).astype(float).tolist(),
                limiting_cell=dict(predicted=float(predicted[cell,field]),full_tangent=float(full_action[cell,field]),
                                   native=float(native_action[cell,field]),residual=float(value[cell,field])),
                trial_scores={'full':float(np.max(abs(native_plus)/m.ATOL)),
                              'one_over_128':float(np.max(abs(tiny)/m.ATOL))})
            report['arms'].append(item)
            np.savez_compressed(OUT/(name+'.npz'),direction=change,composition_coefficients=star.composition_data['coefficients'],
                                composition_inverse=star.composition_data['inverse'],native_action=native_action,
                                full_tangent_action=full_action,predicted=predicted,residual=value)
            e.write(OUT/'result.json',report)
            print('FROZEN ITERATE ARM',json.dumps(item),flush=True)
        assert report['arms'][0]['trial_scores']['full'] > norm
    e.write(OUT/'manifest.json',dict(sha256={p.relative_to(e.ROOT).as_posix():e.digest(p)
        for p in [Path(__file__),*OUT.iterdir()] if p.is_file()}))


if __name__ == '__main__':
    run()
