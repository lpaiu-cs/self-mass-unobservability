"""Counterexample candidate: isolate nonuniform-mesh acoustic energy injection.

Only the pressure face interpolation changes in the constant-medium control.
This does not modify, certify or repair the completed nonlinear GR trajectories.
"""
from pathlib import Path
import json
import os

import numpy as np
from scipy.linalg import eigvals
import sympy as sp

import gr_time_convergence_cause as control


def symbolic():
    t, p, q, u, v = sp.symbols('t p q u v', real=True)
    face_velocity = (1-t)*u+t*v
    original_pressure = (1-t)*p+t*q
    dual_pressure = t*p+(1-t)*q
    work = lambda pressure: (p-q)*face_velocity+u*(pressure-p)+v*(q-pressure)
    assert sp.expand(work(original_pressure)-(1-2*t)*(p-q)*(u-v)) == 0
    assert sp.expand(work(dual_pressure)) == 0
    assert work(original_pressure).subs({t:sp.Rational(3, 5), p:1, q:0, u:1, v:0}) < 0
    return dict(classification='Proven', original_face_work='A*(1-2*t)*(p_L-p_R)*(v_L-v_R)',
        dual_face_work=0, scope='Exact interior-face algebra for constant-medium linear acoustic energy with the stated cell-volume weights. This is not a theorem for the full GR model.')


def operators(star, cells):
    face_v = np.zeros((star.n+1, cells), dtype=control.ld)
    face_p, face_dual = face_v.copy(), face_v.copy()
    t = (star.rf[1:-1]-star.r[:-1])/np.diff(star.r)
    for j in range(cells):
        value = np.zeros(star.n, dtype=control.ld)
        value[j] = 1
        face_v[:, j] = star.faces(value, odd=True)
        face_p[:, j] = star.faces(value)
        face_dual[1:-1, j] = value[:-1]+(1-t)*(value[1:]-value[:-1])
        face_dual[0, j], face_dual[-1, j] = value[0], value[-1]
    def divergence(face):
        return np.asarray(np.diff(4*np.pi*star.rf[:, None]**2*face, axis=0)[:cells]/star.volume[:cells, None], float)
    D = divergence(face_v)
    source = np.diag(np.asarray(star.area_difference_over_volume[:cells], float))
    G, dual = divergence(face_p)-source, divergence(face_dual)-source
    volume = np.asarray(star.volume[:cells], float)
    adjoint = -D.T*volume[None, :]/volume[:, None]
    error = np.linalg.norm((dual-adjoint)*volume[:, None], np.inf)/np.linalg.norm(adjoint*volume[:, None], np.inf)
    assert error < 1e-12, error
    return D, G, dual, volume, float(error)


def main():
    output = control.SOURCE.parent/'gr-acoustic-adjoint-cause'
    assert not output.exists()
    control.initialize(128)
    sources = [Path(__file__), Path(control.__file__), Path(control.gr.e.__file__),
               Path(control.gr.prior.__file__), control.SOURCE/'time-refinement.json',
               output.parent/'gr-time-cause-central-128/manifest.json',
               output.parent/'gr-time-cause-central-128/exact-time/manifest.json',
               output.parent/'gr-time-cause-central-256/manifest.json']
    output.mkdir()
    plan = dict(classification='Counterexample candidate', sound_speed_cm_s=1.16e8,
        cells=[64, 128, 256],
        assumptions='Constant medium, same actual radial cells, zero perturbations outside each window. Symmetric acoustic fluxes and the original cell-pressure geometric source. Compare only original face pressure weights with their dual weights; no EOS or gravity, no full-GR rerun.',
        sha256={str(p.relative_to(control.ROOT)):control.gr.e.digest(p) for p in sources})
    control.gr.e.write(output/'plan.json', plan)
    rows, spectra = [], {}
    for cells in plan['cells']:
        D, G, dual, volume, error = operators(control.STAR, cells)
        matrices = [np.block([[np.zeros_like(D), -plan['sound_speed_cm_s']*D],
                             [-plan['sound_speed_cm_s']*gradient, np.zeros_like(D)]])
                    for gradient in [G, dual]]
        original, corrected = [eigvals(matrix) for matrix in matrices]
        assert original.real.max() > 10 and abs(corrected.real).max() < 1e-10
        spectra[f'original_{cells}'], spectra[f'dual_{cells}'] = original, corrected
        rows.append(dict(cells=cells, adjoint_identity_relative=error,
            original_maximum_real=float(original.real.max()),
            dual_maximum_absolute_real=float(abs(corrected.real).max()),
            leading_original_modes=[[float(v.real), float(v.imag)] for v in
                                    sorted(original, key=lambda v:v.real, reverse=True)[:4]]))
    core128 = output.parent/'gr-time-cause-central-128'
    core256 = output.parent/'gr-time-cause-central-256'
    with np.load(core128/'trajectories.npz') as a, np.load(core256/'trajectories.npz') as b:
        window_errors = {name: (np.max(abs(a[name]-b[name]), axis=(0, 1))/np.maximum(
            np.max(abs(a[name]), axis=(0, 1)), 1e-300)).tolist() for name in a.files}
    assert max(max(v) for v in window_errors.values()) < 1e-8
    result = dict(classification='Counterexample candidate', symbolic=symbolic(),
        spatial_controls=rows, full_control_128_256_history_relative=window_errors,
        original_GR_verdict_preserved=True, full_GR_repaired=False)
    np.savez_compressed(output/'spectra.npz', **spectra)
    control.gr.e.write(output/'result.json', result)
    plot(output, core128, spectra)
    control.gr.e.write(output/'manifest.json', dict(sha256={p.name:control.gr.e.digest(p)
        for p in output.iterdir() if p.is_file()}))
    print(json.dumps(result), flush=True)


def plot(output, source, spectra):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7))
    with np.load(source/'trajectories.npz') as trajectories:
        for r, color in [(1, '#64748b'), (2, '#d97706'), (4, '#2563eb')]:
            times, actual = [], []
            for step in range(1, 40):
                with np.load(control.SOURCE/f'path-{r}'/f'step-{step*r:04d}.npz') as cp:
                    times.append(float(cp['time_seconds']))
                    actual.append(float(cp['delta'][0, 2]))
            axes[0].plot(times, trajectories[f'path_{r}'][:, 0, 2], color=color, label=f'Fluid control r={r}')
            axes[0].scatter(times, actual, s=12, color=color, marker='x')
    axes[0].set(title='Actual GR (crosses) and fluid control', xlabel='Coordinate time [s]', ylabel='Center velocity / c')
    axes[0].legend(fontsize=8)
    exact = json.loads((source/'exact-time/result.json').read_text())['exact_time_errors']
    for j, name in enumerate(['Density', 'Temperature', 'Velocity']):
        errors = np.array([row['exact_time_error'][j] for row in exact])
        axes[1].loglog([row['refinement'] for row in exact], errors/errors[0], 'o-', label=name)
    axes[1].set(title='Error against exact-time affine control', xlabel='Time refinement r', ylabel='Error / error at r=1')
    axes[1].legend(fontsize=8)
    for name, color, label in [('original_128', '#dc2626', 'Original interpolation'), ('dual_128', '#2563eb', 'Compatible interpolation')]:
        values = spectra[name]
        axes[2].scatter(values.real, values.imag, s=13, alpha=.8, color=color, label=label)
    axes[2].axvline(0, color='#94a3b8', lw=.7)
    axes[2].set(title='Constant-medium spatial operator', xlabel='Growth rate Re(lambda) [1/s]', ylabel='Im(lambda) [rad/s]')
    axes[2].legend(fontsize=8, loc='upper left')
    for ax in axes:
        ax.grid(alpha=.2)
    fig.suptitle('Cause investigation: central GR error reproduced by an unstable spatial acoustic operator', fontsize=13)
    fig.text(.5, .015, 'Counterexample candidate. Finite controls; original full-GR convergence verdict remains FAILED.', ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .035, 1, .93))
    fig.savefig(output/'cause.png', dpi=160)
    plt.close(fig)


if __name__ == '__main__':
    assert set(os.sched_getaffinity(0)) <= set(range(16))
    main()
