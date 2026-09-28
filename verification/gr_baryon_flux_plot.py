"""Portable rendering of the corrected coupled paths and local energy defect."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'outputs/direct-eos-gr33/gr-baryon-flux-evolution'


def export():
    initial = np.load(OUT/'initial.npz')
    data = dict(radius=np.asarray(initial['radius']/initial['radius'][-1], float))
    for n in [4, 8, 16]:
        final = np.load(OUT/f'path-{n}/step-{n:04d}.npz')
        history = json.loads((OUT/f'path-{n}/progress.json').read_text())['rows']
        data[f'delta{n}'] = np.asarray(final['delta'][:, :4], float)
        data[f'history{n}'] = np.array([[r['time_seconds']*1000, r['maximum_changes'][1]] for r in history])
    for label, folder, n in [('heat8', OUT.parent/'gr-flux-coupled-evolution', 8), ('both8', OUT, 8), ('both16', OUT, 16)]:
        data[label] = np.asarray(np.load(folder/f'budget-{n}.npz')['normalized_cell_energy_defect'], float)
    path = OUT/'figure-data.npz'
    assert not path.exists()
    np.savez_compressed(path, **data)


def render():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data = np.load(OUT/'figure-data.npz')
    result = json.loads((OUT/'time-refinement.json').read_text())
    assert result['passed']
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for n, color in zip([4, 8, 16], ['#c78210', '#2ca2b2', '#194d8c']):
        axes[0, 0].plot(*data[f'history{n}'].T, label=f'{n} steps', color=color)
        axes[0, 1].plot(data['radius'], data[f'delta{n}'][:, 1], label=f'{n} steps', color=color)
    axes[0, 0].set(xlabel='Coordinate time (ms)', ylabel=r'$\max |\Delta\ln T|$', title='Actual coupled time path')
    axes[0, 1].set(xlabel='r / outer node radius', ylabel=r'$\Delta\ln T$', xlim=(.995, 1), title='Outer envelope at 0.4113 ms')
    for key, label, color in [('heat8', 'Heat correction only, 8 steps', '#be493c'),
                               ('both8', 'Heat + baryon, 8 steps', '#2ca2b2'),
                               ('both16', 'Heat + baryon, 16 steps', '#194d8c')]:
        axes[1, 0].semilogy(data['radius'], np.maximum(abs(data[key]), 1e-18), label=label, color=color)
    axes[1, 0].set(xlabel='r / outer node radius', ylabel='|Cell energy defect| / initial heat capacity',
                   title='Stable energy increments + shared face exchange')
    errors = np.asarray(result['endpoint_maximum_differences'])
    for j, label in [(1, 'Temperature'), (3, 'Heat flux')]:
        axes[1, 1].plot([0, 1], errors[:, j]/errors[0, j], 'o-', label=f"{label}: order {result['observed_orders'][j]:.2f}")
    axes[1, 1].set(xticks=[0, 1], xticklabels=['4 vs 8', '8 vs 16'], ylabel='Endpoint difference / first difference',
                   title='Finite time refinement, same spatial model', yscale='log')
    for ax in axes.flat:
        ax.grid(alpha=.2)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle('5,735-cell coupled evolution after two spatial consistency fixes\nNumerical candidate; computational outer boundary; no physical EOS certificate', fontsize=11)
    path = OUT/'coupled-conservation.png'
    assert not path.exists()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    files = [Path(__file__), path, OUT/'figure-data.npz', OUT/'time-refinement.json',
             OUT/'budget-8-manifest.json', OUT/'budget-16-manifest.json',
             OUT.parent/'gr-flux-coupled-evolution/budget-8-manifest.json']
    (OUT/'figure-manifest.json').write_text(json.dumps(dict(sha256={p.relative_to(ROOT).as_posix():
          hashlib.sha256(p.read_bytes()).hexdigest() for p in files}), indent=2)+'\n')
    print(path)


if __name__ == '__main__': globals()[sys.argv[1]]()
