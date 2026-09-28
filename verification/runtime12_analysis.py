"""Adjudicate live response, derivative convergence, and nonlinear displacement gates."""
import json
import math
import numpy as np

import comparator_audit as comp
import nuisance_audit as audit
from physical_drive_completion import drive, stencil

LIVE=audit.OUT/'runtime12'


def fit_columns(x,y,dof):
    norms=np.linalg.norm(x,axis=0)
    q,r=np.linalg.qr(x/norms,mode='reduced')
    assert min(abs(np.diag(r)))>1e-10
    inv=np.linalg.inv(r)
    coef=np.linalg.solve(r,q.T@y)/norms
    scale=np.linalg.norm(y-q@(q.T@y))/math.sqrt(dof)
    sigma=scale*np.sqrt(np.sum(inv*inv,axis=1))/norms
    return coef,sigma


def main():
    g=comp.prepare(); sw=g['inp']['sw']; n=len(sw); d=drive()
    transient=[]
    for tau in [2.,52.,500.]:
        full=np.load(LIVE/f'transient_{tau:g}_1e-08.npz')['dcol']
        half=np.load(LIVE/f'transient_{tau:g}_5e-09.npz')['dcol']
        v=audit.sc.proj_out(g['qf'],sw*half)
        delta=audit.sc.proj_out(g['qf'],sw*(half-full))
        change=float(np.linalg.norm(delta)/np.linalg.norm(v))
        fits=[]
        for a in comp.AMPLITUDES:
            m=comp.metric(g,a); h=g['h']
            y=g['yp']+h@((m['root']-1)*g['y'])
            vt=v+h@((m['root']-1)*(h.T@v))
            x=h@(m['c']@stencil(g,d,tau))
            coef0,sig0=fit_columns(x,y,n-90-2)
            coef,sig=fit_columns(np.column_stack([x,vt]),y,n-90-3)
            fits.append(dict(a=a,beta=float(coef[1]),sigma_beta=float(sig[1]),
                             initial_coupling=float(coef[2]),sigma_initial_coupling=float(sig[2]),
                             sigma_beta_ratio=float(sig[1]/sig0[1]),
                             conditional_U=float(audit.intervals(coef[1],sig[1]))))
        transient.append(dict(tau=tau,weighted_derivative_change=float(np.linalg.norm(sw*(half-full))/np.linalg.norm(sw*half)),
                              nuisance_projected_derivative_change=change,convergence_5percent_pass=change<.05,
                              fits=fits))
    meta=json.loads((audit.ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())
    jnew=np.column_stack([np.load(LIVE/f'jac_{j:02d}_0.5.npz')['dcol'] for j in range(28)])
    derivative=[]
    for j,name in enumerate(meta['columns']):
        old=g['inp']['J'][:,j]; new=jnew[:,j]
        row=dict(name=name,half_vs_archived_weighted_relative=float(np.linalg.norm(sw*(new-old))/np.linalg.norm(sw*new)))
        if j>=21:
            quarter=np.load(LIVE/f'jac_{j:02d}_0.25.npz')['dcol']
            row['quarter_vs_half_weighted_relative']=float(np.linalg.norm(sw*(quarter-new))/np.linalg.norm(sw*quarter))
        derivative.append(row)
    bnew=g['b'].copy(); bnew[:,:28]=sw[:,None]*jnew/g['norms'][None,:28]
    q,s,_=np.linalg.svd(bnew,full_matrices=False)
    assert len(s)==90 and s[-1]>s[0]*1e-12
    maxsin=float(np.linalg.norm(audit.sc.proj_out(q,g['qf']),ord=2))
    perturbation=float(np.linalg.norm(bnew-g['b'],ord=2))
    residual=audit.sc.proj_out(q,g['inp']['sw']*g['inp']['res0'])
    projection=[]
    for tau in audit.LAGS:
        raw=g['c6']@stencil(g,d,tau)
        xold=audit.sc.proj_out(g['qf'],raw); xnew=audit.sc.proj_out(q,raw)
        co,so=fit_columns(xold,g['yp'],n-90-2); cn,sn=fit_columns(xnew,residual,n-90-2)
        projection.append(dict(tau=tau,new_sigma_over_archived=float(sn[1]/so[1]),
                               new_beta=float(cn[1]),new_conditional_U=float(audit.intervals(cn[1],sn[1]))))
    nonlinear=[]
    for fraction in [.001,.003,.01]:
        z=np.load(LIVE/f'nonlinear_gap_{fraction}.npz')
        assert str(z['input_sha256'])==audit.hashlib.sha256((LIVE/'nonlinear-gap-input.npz').read_bytes()).hexdigest()
        actual=z['actual_us']; predicted=z['prediction_us']
        error=sw*(actual-predicted)
        relative=float(np.linalg.norm(error)/np.linalg.norm(sw*predicted))
        nonlinear.append(dict(fraction=fraction,weighted_nonlinear_error_ratio=relative,
                              residual_error_unit_noise_norm=float(np.linalg.norm(error)),
                              remaining_error_after_nuisance_norm=float(np.linalg.norm(audit.sc.proj_out(g['qf'],error))),
                              max_actual_us=float(np.max(np.abs(actual))),linear_5percent_pass=relative<.05))
    result=dict(status='Imported from prior work',scope='New live external evaluations; registered finite amplitudes, steps and displacement only, no full nonlinear noise/pulse fit or EOS matching',
                transient=transient,derivative_rows=derivative,
                halfstep_basis=dict(rank=90,maximum_principal_sine=maxsin,normalized_matrix_change_norm=perturbation,
                                    archived_smin=float(g['s'][-1]),rigorous_derivative_error_bound=None),
                physical_diagonal_projection_comparison=projection,nonlinear_gap=nonlinear,
                original_displacement_rejection=json.loads((LIVE/'nonlinear-domain-rejection.json').read_text()))
    (audit.OUT/'runtime12-analysis.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('transient convergence',[(r['tau'],r['nuisance_projected_derivative_change']) for r in transient])
    print('derivative max',max(r['half_vs_archived_weighted_relative'] for r in derivative),'principal sine',maxsin)
    print('nonlinear',nonlinear)


if __name__=='__main__': main()
