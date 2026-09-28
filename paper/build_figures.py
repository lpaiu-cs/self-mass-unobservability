"""Plot committed interval and coverage results; never runs a timing fit."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def main():
    data = json.loads((ROOT / "request10_external/sep_dynamic/sep_phase_marg_10_8e.json").read_text())
    rows = sorted((float(k.removeprefix("tau_")), v) for k, v in data["anchors"].items())
    tau = [t for t, _ in rows]
    full = [r["u95pm_K10_fullrank"] for _, r in rows]
    trunc = [r["u95pm_K10"] for _, r in rows]
    fisher = [r["u95pm_fisher"] for _, r in rows]
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.15), layout="constrained")
    for values, label, color, marker in [
        (full, "Full nuisance, K=10", "#9b2c2c", "s"),
        (trunc, "Truncated, K=10", "#1f567d", "o"),
        (fisher, "Truncated, K=1", "#555555", "^"),
    ]:
        axes[0].loglog(tau, values, marker=marker, color=color, label=label)
    axes[0].set_ylabel(r"Gaussian envelope $U_\beta$")
    axes[0].legend(fontsize=8, loc="upper left", frameon=False)
    axes[0].set_title("(a) Conditional interval construction", fontsize=10)
    ratio = [f / t for f, t in zip(full, trunc)]
    axes[1].semilogx(tau, ratio, "o-", color="#9b2c2c")
    axes[1].axhline(1, color="#777777", lw=0.8, ls=":")
    axes[1].set_ylabel("Full / truncated (K=10)")
    axes[1].set_ylim(0, 20)
    axes[1].set_title("(b) Nuisance-space dependence", fontsize=10)
    for ax in axes:
        ax.set_xlabel(r"Relaxation time $\tau_\chi$ (day)")
        ax.set_xticks(tau, [str(int(t)) for t in tau])
        ax.grid(alpha=0.2, which="major")
        ax.spines[["top", "right"]].set_visible(False)
    output = ROOT / "paper/figures/conditional-intervals.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    print(f"Wrote {output}")

    folder=ROOT/'outputs/research-completion'
    coverage=json.loads((folder/'coverage-audit.json').read_text())['rows']
    estimated=json.loads((folder/'estimated-covariance-audit.json').read_text())['rows']
    scenarios=['white','extra_fourier_025','extra_fourier_1']
    fig,ax=plt.subplots(figsize=(6.8,3.15),layout='constrained')
    for method,label,color,marker in [
        ('truncated_diag','Truncated, diagonal','#777777','x'),
        ('full_diag','Full, diagonal','#9b2c2c','s'),
        ('full_matched_GLS','Full, oracle GLS','#1f567d','o'),
    ]:
        values=[min(r['U_K1']['fraction'] for r in coverage if r['scenario']==sc and r['method']==method) for sc in scenarios]
        ax.plot(range(3),values,marker=marker,color=color,label=label)
    values=[min(r['U_K1']['fraction'] for r in estimated if r['true_a']==a) for a in [0.,.25,1.]]
    ax.plot(range(3),values,'D--',color='#19734b',label='Full, estimated GLS',markersize=4)
    ax.axhline(.95,ls=':',color='black',lw=.8)
    ax.set(xticks=range(3),xticklabels=['White noise','Extra Fourier RMS 0.25','Extra Fourier RMS 1'],
           ylim=(.5,1.01),ylabel='Minimum pointwise coverage')
    ax.legend(loc='lower left',fontsize=8,frameon=False)
    ax.grid(alpha=.2,axis='y'); ax.spines[['top','right']].set_visible(False)
    output=ROOT/'paper/figures/coverage-validation.pdf'
    fig.savefig(output,metadata={'CreationDate':None,'ModDate':None}); plt.close(fig)
    print(f'Wrote {output}')

    comparison=json.loads((folder/'comparator-audit.json').read_text())['rows']
    phase=json.loads((folder/'phase-refinement.json').read_text())['rows']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.25),layout='constrained')
    names=['P0','P1','P2','P3','P4','even_P4']
    lower=[min(r['relative_information'] for r in comparison if r['comparator']==name) for name in names]
    upper=[max(r['relative_information'] for r in comparison if r['comparator']==name) for name in names]
    for x,(lo,hi) in enumerate(zip(lower,upper)):
        axes[0].plot([x,x],[lo,hi],color='#1f567d',lw=2)
        axes[0].plot([x,x],[lo,hi],'o',color='#1f567d',ms=4)
    axes[0].set(yscale='log',xticks=range(6),xticklabels=['P0','P1','P2','P3','P4','E4'],
                ylim=(1e-7,2),ylabel='Information / instantaneous-only',xlabel='Allowed derivative comparator')
    axes[0].text(.04,.035,'P5: exact zero',transform=axes[0].transAxes,fontsize=8)
    axes[0].set_title('(a) Retained-information range',fontsize=10)
    phase=sorted([r for r in phase if r['a']==0.08316109496333533],key=lambda r:r['tau'])
    lags=[r['tau'] for r in phase]
    axes[1].loglog(lags,[r['refined_U'] for r in phase],'o-',color='#1f567d',label='Local-refined U')
    axes[1].loglog(lags,[r['continuous_upper_bound'] for r in phase],'s--',color='#9b2c2c',label='Analytic upper envelope')
    axes[1].set(xlabel='Relaxation time (day)',ylabel='Conditional unit-drive U')
    axes[1].set_xticks(lags,[str(int(v)) for v in lags])
    axes[1].legend(fontsize=7.5,frameon=False,loc='lower right')
    axes[1].set_title('(b) Three-phase envelope, fixed covariance',fontsize=9)
    for ax in axes:
        ax.grid(alpha=.2,axis='y'); ax.spines[['top','right']].set_visible(False)
    output=ROOT/'paper/figures/comparator-phase-validation.pdf'
    fig.savefig(output,metadata={'CreationDate':None,'ModDate':None}); plt.close(fig)
    print(f'Wrote {output}')


if __name__ == "__main__":
    main()

    folder=ROOT/'outputs/research-completion'
    calibration=json.loads((folder/'simultaneous-validation.json').read_text())['rows']
    live=json.loads((folder/'runtime12-analysis.json').read_text())['transient']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
    c=[r['inclusion'] for r in calibration]
    values=[r['fraction'] for r in c]
    axes[0].errorbar(range(5),values,yerr=[[v-r['lo95'] for v,r in zip(values,c)],
                                        [r['hi95']-v for v,r in zip(values,c)]],fmt='o',color='#1f567d',capsize=3)
    axes[0].axhline(.95,color='#777777',ls=':',lw=1)
    axes[0].set(xticks=range(5),xticklabels=['0','.25','1','4','Shape\nstress'],ylim=(.94,.97),
                xlabel='Generating covariance amplitude',ylabel='Joint-region inclusion')
    axes[0].set_title('(a) Independent validation',fontsize=10)
    rows=[next(f for f in r['fits'] if f['a']==0.08316109496333533) for r in live]
    axes[1].semilogx([r['tau'] for r in live],[r['sigma_beta_ratio'] for r in rows],'o-',color='#19734b')
    axes[1].axhline(1,color='#777777',ls=':',lw=1)
    axes[1].set(xticks=[2,52,500],xticklabels=['2','52','500'],ylim=(.999,1.007),
                xlabel='Relaxation time (day)',ylabel='Standard error with / without transient')
    axes[1].set_title('(b) Dedicated transient response',fontsize=10)
    for ax in axes:
        ax.grid(alpha=.2,axis='y'); ax.spines[['top','right']].set_visible(False)
    output=ROOT/'paper/figures/remaining-levers-validation.pdf'
    fig.savefig(output,metadata={'CreationDate':None,'ModDate':None}); plt.close(fig)
    print(f'Wrote {output}')
