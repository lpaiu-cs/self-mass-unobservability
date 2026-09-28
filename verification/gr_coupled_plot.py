"""Render portable displays without loading the WSL-only EOS/long-double runtime."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
folder=ROOT/'outputs/direct-eos-gr33/gr-coupled-evolution'
comparison=json.loads((folder/'time-comparison.json').read_text())
assert comparison['passed']
data=np.load(folder/'figure-data.npz')
fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
for row,color in zip(comparison['runs'],['#e69f00','#56b4e9','#0072b2']):
    stem=f"k{row['steps']}_"
    x,delta,history=[data[stem+k] for k in ['radius','delta','history']]
    label=f"{row['steps']} steps"
    axes[0,0].plot(x,delta[:,1],color=color,label=label,lw=1.3)
    axes[0,1].plot(history[:,0],history[:,1],color=color,label=label,lw=1.3)
    axes[1,0].plot(x,delta[:,2]*data['light_cm_per_s'],color=color,label=label,lw=1.3)
    axes[1,1].plot(x,data[stem+'mass_delta'],color=color,label=label,lw=1.3)
axes[0,0].set(xlabel='r / outer node radius',ylabel=r'$\Delta\ln T$',title='Envelope temperature change',xlim=(.85,1))
axes[0,1].set(xlabel='Coordinate time (ms)',ylabel=r'$\max |\Delta\ln T|$',title='Nonlinear time path')
axes[1,0].set(xlabel='r / outer node radius',ylabel='Velocity (cm/s)',title='Fluid response')
axes[1,1].set(xlabel='r / outer node radius',ylabel=r'$\Delta m$ (geometric cm)',title='Enclosed energy-mass response')
for ax in axes.flat:
    ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8)
fig.suptitle('5,735-cell native-EOS coupled evolution\nNumerical candidate with a computational outer boundary',fontsize=12)
path=folder/'coupled-time-path.png';assert not path.exists()
fig.savefig(path,dpi=180);plt.close(fig)
files=[path,folder/'figure-data.npz',folder/'time-comparison.json',Path(__file__)]
(folder/'figure-manifest.json').write_text(json.dumps(dict(sha256={
    p.relative_to(ROOT).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}),indent=2)+'\n')
print('COUPLED FIGURE',path)
