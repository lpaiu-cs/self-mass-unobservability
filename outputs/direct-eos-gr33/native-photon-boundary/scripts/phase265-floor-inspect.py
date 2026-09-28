"""Read-only: how close the accepted phase-259 stages came to the 1e-13 material gate."""
import json
from pathlib import Path
import numpy as np
for n in [64, 128]:
    c = json.loads(Path(f'native-short-return259-work/checks-{n}.json').read_text())
    last = np.array([s[-1]['material_relative'] for s in c['newton']])  # accepted proposal per step
    tries = np.array([len(s) for s in c['newton']])
    comp = last.reshape(len(last), 2, 4)
    worst = comp.max(axis=1)  # per step, per component (max over the two stages)
    names = ['Etilde', 'H', 'B', 'S']
    print(f'clock {n}: steps {len(last)}, proposals used: {np.bincount(tries).tolist()}')
    for k, name in enumerate(names):
        v = worst[:, k]
        print(f'  {name}: max {v.max():.4e} (step {int(v.argmax())+1}), median {np.median(v):.3e}, steps >5e-14: {int((v > 5e-14).sum())}, >8e-14: {int((v > 8e-14).sum())}, >9e-14: {int((v > 9e-14).sum())}')
