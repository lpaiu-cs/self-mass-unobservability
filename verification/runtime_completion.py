"""Request 12 external WSL producer. Never modifies archived runtime or REQUEST10."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import shutil

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'outputs/research-completion/runtime12'
RUNTIME = Path.home() / 'work/nutimo_pilot'
RUN = RUNTIME / 'run_request12'
PAR = 'parfile-planetGR-max-bestfit'
TIM = '0337_20211005-sorted-sun5deg-res25microsec-58631_58780_clipped.tim'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    mode = sys.argv[1]
    worker = int(sys.argv[2]) if len(sys.argv)>2 else 0
    workers = int(sys.argv[3]) if len(sys.argv)>3 else 1
    OUT.mkdir(parents=True, exist_ok=True)
    src = RUNTIME / ('nutimo_sepdyn/src' if mode == 'probe' else 'nutimo_request12/src')
    sys.path.insert(0, str(src))
    run = RUN
    if mode == 'jac':
        run = RUNTIME/f'run_request12_jac{worker}'
        if not run.exists(): shutil.copytree(RUN,run)
    os.chdir(run)
    import python_Fittriple_interface as pfi
    assert Path(pfi.__file__).parent == src
    base = np.load(ROOT/'request10_external/baseline_planetGR.npz', allow_pickle=True)
    for key in ['SEPDYN_A', 'SEPDYN_W', 'SEPDYN_PH', 'SEPDYN_TAU']:
        os.environ.pop(key, None)
    t0 = time.perf_counter()
    fit = pfi.PyFittriple(PAR, TIM)
    fit.Compute_lnposterior()
    res = fit.Get_time_residuals().copy()
    error = float(np.max(np.abs(res-base['res'])))
    report = dict(mode=mode, source_sha256=sha(src/'AllTheories3Bodies.cpp'),
                  library_sha256=sha(src/'libFittriplecpp.so'), interface_sha256=sha(pfi.__file__),
                  par_sha256=sha(PAR), tim_sha256=sha(TIM),
                  zero_baseline_max_us=error, baseline_seconds=time.perf_counter()-t0)
    label = f'jac{worker}' if mode=='jac' else mode
    (OUT/f'{label}-preflight.json').write_text(json.dumps(report, indent=2)+'\n')
    print('PREFLIGHT', report, flush=True)
    assert error < 1e-7, 'Live baseline differs from frozen record; no promotion'
    if mode == 'probe':
        return
    if mode == 'jac':
        meta = json.loads((ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())
        names = [str(base['names'][j]) for j in base['fmap']]
        assert names == meta['columns']
        scales = np.asarray(base['scales'])[base['fmap']].astype(float)
        for j, name in enumerate(names):
            if j % workers != worker: continue
            for fraction in ([.5, .25] if j >= 21 else [.5]):
                path = OUT/f'jac_{j:02d}_{fraction}.npz'
                if path.exists():
                    continue
                h = float(meta['abs_steps'][j])*fraction
                outputs = []
                for sign in [1, -1]:
                    shift = np.zeros(len(names)); shift[j] = sign*h/scales[j]
                    fit.Set_fitted_parameter_relativeshifts(shift)
                    fit.Compute_lnposterior()
                    outputs.append(fit.Get_time_residuals().copy())
                rp, rm = outputs
                assert np.isfinite(rp).all() and np.isfinite(rm).all()
                np.savez(path, plus=rp, minus=rm, h=h, dcol=(rp-rm)/(2*h), name=name)
                print('JAC', j, name, fraction, 'elapsed', time.perf_counter()-t0, flush=True)
    elif mode == 'nonlinear':
        data=np.load(OUT/'nonlinear-gap-input.npz')
        delta=data['parameter_delta']
        scales=np.asarray(base['scales'])[base['fmap']].astype(float)
        names=[str(base['names'][j]) for j in base['fmap']]
        params=dict(zip(base['names'],base['params']))
        rejected=[]
        for fraction in [.25,.5,1.]:
            trial={name:params[name]+fraction*delta[j] for j,name in enumerate(names)}
            eccentricity=float(np.hypot(trial['eta_extra1'],trial['kappa_extra1']))
            assert eccentricity>=1, 'Revisit previously observed domain failure'
            rejected.append(dict(fraction=fraction,eccentricity=eccentricity,reason='Outside bound-orbit domain; runtime would clamp eccentricity'))
        (OUT/'nonlinear-domain-rejection.json').write_text(json.dumps(rejected,indent=2)+'\n')
        for fraction in [.001,.003,.01]:
            fit.Set_fitted_parameter_relativeshifts(fraction*delta/scales)
            fit.Compute_lnposterior()
            actual=fit.Get_time_residuals().copy()-res
            np.savez(OUT/f'nonlinear_gap_{fraction}.npz',actual_us=actual,
                     prediction_us=fraction*data['prediction_us'],fraction=fraction,
                     input_sha256=sha(OUT/'nonlinear-gap-input.npz'))
            print('NONLINEAR',fraction,'max_us',float(np.max(np.abs(actual))),flush=True)
    elif mode == 'transient':
        for tau in [2., 52., 500.]:
            for amplitude in [1e-8, 5e-9]:
                path = OUT/f'transient_{tau:g}_{amplitude:g}.npz'
                if path.exists():
                    continue
                outputs = []
                for sign in [1, -1]:
                    os.environ['SEPDYN_A'] = repr(sign*amplitude)
                    os.environ['SEPDYN_TAU'] = repr(tau/24.077445558945514)
                    # Integrate_Allways rereads the environment on each computation.
                    fit.Compute_lnposterior()
                    outputs.append(fit.Get_time_residuals().copy())
                rp, rm = outputs
                dev = float(max(np.max(np.abs(rp-res)), np.max(np.abs(rm-res))))
                np.savez(path, plus=rp, minus=rm, amplitude=amplitude, tau_days=tau,
                         dcol=(rp-rm)/(2*amplitude), maxdev_us=dev)
                assert np.isfinite(rp).all() and np.isfinite(rm).all() and dev < 1000
                assert dev > 1e-8, 'Drive was not applied'
                print('TRANSIENT', tau, amplitude, dev, 'elapsed', time.perf_counter()-t0, flush=True)
        os.environ['SEPDYN_A'] = '0'
        fit.Compute_lnposterior()
        assert np.max(np.abs(fit.Get_time_residuals()-res)) < 1e-7
    else:
        raise ValueError(mode)
    print('DONE', mode, flush=True)


if __name__ == '__main__':
    main()
