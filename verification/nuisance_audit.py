"""Request 11.1: new read-only audit of frozen timing arrays. No timing engine."""
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "request10_external/scripts"))
import sep_common as sc

OUT = ROOT / "outputs/research-completion"
LAGS = [2., 5., 18., 52., 200., 500.]
CUTS = [1e-2, 1e-3, 1e-4, 1e-6, 0.]


def ndtr(x):
    x=np.asarray(x)
    return .5+.5*np.fromiter((math.erf(v/math.sqrt(2)) for v in x.flat),dtype=float,count=x.size).reshape(x.shape)


def geometry():
    manifest = json.loads((ROOT / "paper/revision-manifest.json").read_text())
    for name, sha in manifest["sha256"].items():
        if name.startswith("request10_external/"):
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha, name
    inp = sc.load_inputs(str(ROOT / "request10_external"))
    meta = json.loads((ROOT / "request10_external/finite_jacobian_v2_meta.json").read_text())
    t, sw = inp['t'], inp['sw']
    names = meta['columns'] + ['offset']
    columns = [sw[:,None]*inp['J'], sw[:,None]]
    fourier = []
    for j in range(1,31):
        phase = 2*np.pi*j*(t-t.min())/np.ptp(t)
        fourier += [np.cos(phase), np.sin(phase)]
        names += [f'fourier_{j}_cos', f'fourier_{j}_sin']
    columns += [sw[:,None]*np.column_stack(fourier), (sw*inp['jD'])[:,None]]
    names += ['static_SEP_guard']
    raw = np.column_stack(columns)
    norms = np.linalg.norm(raw,axis=0)
    assert np.all(np.isfinite(raw)) and np.all(norms>0)
    b = raw/norms
    u,s,vh = np.linalg.svd(b,full_matrices=False)
    qt = sc.build_nuisance(inp)['Q']
    qf = sc.build_nuisance(inp,sv_cut=0,keep_guard=False)['Q']
    ud,sd,_ = np.linalg.svd(sc.proj_out(qt,qf),full_matrices=False)
    d = ud[:,sd>0.5]
    c6,order = sc.projected_sixcols(inp,np.empty((inp['N'],0)))
    return dict(inp=inp, names=names, norms=norms, b=b, u=u, s=s, vh=vh,
                qt=qt,qf=qf,d=d,c6=c6,order=order)


def intervals(b,sigma):
    """Vector version of the archived Gaussian-mass definition, not new coverage."""
    b,sigma = np.broadcast_arrays(np.asarray(b),np.asarray(sigma))
    assert np.all(sigma>0)
    lo,hi = np.zeros_like(b),np.abs(b)+10*sigma
    for _ in range(60):
        mid=(lo+hi)/2
        mass=ndtr((mid-b)/sigma)-ndtr((-mid-b)/sigma)
        lo,hi=np.where(mass<.95,mid,lo),np.where(mass>=.95,mid,hi)
    return (lo+hi)/2


def stencil_grid(tau,toffs,oms):
    pairs=[sc.template_stencil(tau,float(t),oms) for t in toffs]
    return np.array([p[0] for p in pairs]),np.array([p[1] for p in pairs])


def fit_grid(gram,score,wc,wb,scale):
    a=np.einsum('ki,ij,kj->k',wc,gram,wc)
    b=np.einsum('ki,ij,kj->k',wc,gram,wb)
    c=np.einsum('ki,ij,kj->k',wb,gram,wb)
    information=c-b*b/a
    assert np.all(a>0) and np.all(information>0)
    beta=(wb@score-b/a*(wc@score))/information
    return beta,scale/np.sqrt(information)


def main():
    for b,s in [(0.,1.),(3.,.2),(-2.,4.)]:
        assert np.isclose(intervals(b,s),sc.u95_of(b,s),rtol=1e-13)
    g=geometry(); inp=g['inp']; n=inp['N']; qt,qf,d=g['qt'],g['qf'],g['d']
    q_qr,_=np.linalg.qr(g['b'],mode='reduced')
    diagnostics={
        'rank_full':qf.shape[1], 'rank_truncated':qt.shape[1], 'omitted_dimension':d.shape[1],
        'minimum_relative_singular_value':float(g['s'][-1]/g['s'][0]),
        'machine_rank':int(np.linalg.matrix_rank(g['b'])),
        'qr_svd_subspace_residual':float(np.linalg.norm(sc.proj_out(qf,q_qr))),
        'orthogonality_full':float(np.linalg.norm(qf.T@qf-np.eye(qf.shape[1]))),
        'orthogonality_truncated':float(np.linalg.norm(qt.T@qt-np.eye(qt.shape[1]))),
        'guard_residual_truncated':float(np.linalg.norm(sc.proj_out(qt,g['b'][:,-1]))),
    }
    assert diagnostics['machine_rank']==90,diagnostics
    assert max(diagnostics[k] for k in diagnostics if 'residual' in k or 'orthogonality' in k)<1e-7,diagnostics
    assert d.shape[1]==qf.shape[1]-qt.shape[1]
    loss=np.sum(sc.proj_out(qt,g['u'])**2,axis=0)
    mode_rows=[]
    for j,sv in enumerate(g['s']):
        vector=g['vh'][j]
        top=np.argsort(np.abs(vector))[-5:][::-1]
        mode_rows.append(dict(mode=j+1,relative_sv=float(sv/g['s'][0]),lost_fraction=float(loss[j]),
            timing_share=float(np.sum(vector[:28]**2)),fourier_share=float(np.sum(vector[29:89]**2)),
            offset_share=float(vector[28]**2),guard_share=float(vector[89]**2),
            top_loadings=[dict(name=g['names'][k],coefficient=float(vector[k])) for k in top]))
    y=inp['sw']*inp['res0']; toffs=np.arange(0.,inp['P_out'],inp['P_in']/24)
    stencils={tau:stencil_grid(tau,toffs,inp['OMS']) for tau in LAGS}
    stored=json.loads((ROOT/'request10_external/sep_dynamic/sep_phase_marg_10_8e.json').read_text())
    rows=[]
    for cut in CUTS:
        nu=sc.build_nuisance(inp,sv_cut=cut,keep_guard=cut>0)
        q=nu['Q']; cp=sc.proj_out(q,g['c6']); yp=sc.proj_out(q,y)
        gram=cp.T@cp; score=cp.T@yp; scale,_=sc.noise_scale(q,inp['sw'],inp['res0'])
        for tau in LAGS:
            beta,sigma=fit_grid(gram,score,*stencils[tau],scale)
            for factor in [1.,10.]:
                upper=intervals(beta,factor*sigma); j=int(np.argmax(upper))
                row=dict(cut=cut,rank=q.shape[1],tau=tau,K=factor,noise_scale=scale,
                    U=float(upper[j]),toff=float(toffs[j]),beta=float(beta[j]),sigma=float(sigma[j]))
                rows.append(row)
                if factor==10 and tau!=500 and cut in (0.,1e-3):
                    key='u95pm_K10_fullrank' if cut==0 else 'u95pm_K10'
                    assert np.isclose(row['U'],stored['anchors'][sc.anchor_key(tau)][key],rtol=1e-5),(row,key)
        print('cut',cut,'rank',q.shape[1],'scale',scale,flush=True)
    cp=sc.proj_out(qt,g['c6']); yp=sc.proj_out(qt,y)
    points=[(2.,stored['anchors']['tau_2']['worst_toff_K10'])]
    points += [(r['tau'],r['toff']) for r in rows if r['K']==10 and r['cut'] in (0.,1e-3)]
    leakage=[]; priors=[]
    for tau,toff in sorted(set(points)):
        wc,wb=sc.template_stencil(tau,toff,inp['OMS']); x=cp@np.column_stack([wc,wb])
        l=x@np.linalg.inv(x.T@x)[:,1]
        sensitivity=float(np.linalg.norm(d.T@l)/np.linalg.norm(l))
        leakage.append(dict(tau=tau,toff=toff,unit_noise_sigma=float(np.linalg.norm(l)),
            omitted_unit_bias_sigma=sensitivity,bias_sigma_by_R={str(r):r*sensitivity for r in [1,3,10,30]},
            fitted_omitted_residual_norm=float(np.linalg.norm(d.T@y))))
        if toff==stored['anchors']['tau_2']['worst_toff_K10']:
            for a in [0.,1.,3.,10.,float('inf')]:
                fraction=1. if np.isinf(a) else a*a/(1+a*a)
                px=x-fraction*d@(d.T@x); py=yp-fraction*d@(d.T@yp)
                inverse=np.linalg.inv(x.T@px); estimate=inverse@(x.T@py)
                priors.append(dict(tau=tau,toff=toff,prior_width='infinity' if np.isinf(a) else a,
                    beta=float(estimate[1]),unit_noise_conditional_sigma=float(np.sqrt(inverse[1,1]))))
                if a==0: assert np.isclose(np.sqrt(inverse[1,1]),np.linalg.norm(l),rtol=1e-9)
                if np.isinf(a):
                    xf=sc.proj_out(qf,g['c6']@np.column_stack([wc,wb]))
                    assert np.isclose(inverse[1,1],np.linalg.inv(xf.T@xf)[1,1],rtol=1e-7)
    result=dict(status='Imported from prior work',interpretation='New audit of frozen linearized arrays, not a timing refit or physical prior.',
        diagnostics=diagnostics,singular_modes=mode_rows,grid_origins=len(toffs),intervals=rows,
        leakage=leakage,artificial_prior_sensitivity=priors,
        input_sha256={p:sha for p,sha in json.loads((ROOT/'paper/revision-manifest.json').read_text())['sha256'].items() if p.startswith('request10_external/')})
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'nuisance-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    with (OUT/'nuisance-intervals.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps(diagnostics,indent=2))
    print('PASS: archived anchor reproduction, projection geometry, interval implementation, prior endpoints.')


if __name__=='__main__':
    main()
