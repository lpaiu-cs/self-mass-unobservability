"""Request 11.2b: finite-grid REML with co-fitted signal and noise covariance."""
import json
import math

import numpy as np

import coverage_audit as coverage
import nuisance_audit as audit

SEED=2026090903
GRID=np.r_[0.,np.geomspace(.01,4.,100)]


def reml_scan(x,y,rest,eigenvalues,dof):
    if y.ndim==1: y=y[:,None]
    weights=1/(1+GRID[:,None]**2*eigenvalues[None,:])
    cc=weights@(x[:,0]**2); cb=weights@(x[:,0]*x[:,1]); bb=weights@(x[:,1]**2)
    information=bb-cb*cb/cc
    assert np.all(information>0)
    score_c=(weights*x[:,0])@y; score_b=(weights*x[:,1])@y
    bhat=(score_b-cb[:,None]/cc[:,None]*score_c)/information[:,None]
    chat=(score_c-cb[:,None]*bhat)/cc[:,None]
    residual=weights@(y*y)+np.asarray(rest)-chat*score_c-bhat*score_b
    assert np.all(residual>0)
    sigma2=residual/dof
    objective=(np.log1p(GRID[:,None]**2*eigenvalues).sum(axis=1)+np.log(cc*information))[:,None]+dof*np.log(sigma2)
    j=np.argmin(objective,axis=0); col=np.arange(y.shape[1])
    result=dict(a=GRID[j],beta=bhat[j,col],sigma_beta=np.sqrt(sigma2[j,col]/information[j]),
                white_scale=np.sqrt(sigma2[j,col]),index=j)
    return result,objective


def main():
    g=audit.geometry(); basis,c6,red,dim_d,rest_dof=coverage.compressed_model(g)
    h=basis[:,dim_d:]; c=c6[dim_d:]; lred=red[dim_d:]
    ev,rotation=np.linalg.eigh(lred@lred.T); ev=np.maximum(ev,0)
    c=rotation.T@c; noise_modes=rotation.T@lred
    rng=np.random.default_rng(SEED); nsim=coverage.NSIM
    white=rng.normal(size=(len(ev),nsim)); rz=rng.normal(size=(lred.shape[1],nsim))
    rest=rng.chisquare(rest_dof,size=nsim)
    # The omitted D noise is fitted by the full model; the remaining rest dof is unchanged.
    dof=g['inp']['N']-g['qf'].shape[1]-2
    y_full=audit.sc.proj_out(g['qf'],g['inp']['sw']*g['inp']['res0'])
    y_data=rotation.T@(h.T@y_full)
    data_rest=float(y_full@y_full-y_data@y_data)
    assert data_rest>0
    previous=json.loads((audit.OUT/'nuisance-audit.json').read_text())
    points=sorted(set([(tau,0.) for tau in audit.LAGS]+[(r['tau'],r['toff']) for r in previous['intervals'] if r['cut'] in [0.,1e-3] and r['K']==10]))
    rows=[]; fits=[]; boundary=[]; worst=1.
    for tau,toff in points:
        wc,wb=audit.sc.template_stencil(tau,toff,g['inp']['OMS'])
        x=c@np.column_stack([wc,wb]); _,sf=coverage.beta_weight(x,np.eye(len(ev)))
        fitted,objective=reml_scan(x,y_data,data_rest,ev,dof)
        direct_x=audit.sc.proj_out(g['qf'],g['c6']@np.column_stack([wc,wb]))
        direct_coef=np.linalg.solve(direct_x.T@direct_x,direct_x.T@y_full)
        direct_rss=np.linalg.norm(y_full-direct_x@direct_coef)**2
        zero_objective=np.linalg.slogdet(x.T@x)[1]+dof*math.log(direct_rss/dof)
        assert np.isclose(objective[0,0],zero_objective,rtol=1e-10)
        fits.append(dict(tau=tau,toff=toff,a=float(fitted['a'][0]),white_scale=float(fitted['white_scale'][0]),
                         beta=float(fitted['beta'][0]),sigma_beta=float(fitted['sigma_beta'][0]),
                         U_K1=float(audit.intervals(fitted['beta'][0],fitted['sigma_beta'][0]))))
        for amplitude in [0.,.25,1.]:
            noise=white+amplitude*(noise_modes@rz)
            est,obj=reml_scan(x,noise,rest,ev,dof)
            # Jointly refitting the signal makes the covariance objective translation-invariant.
            check,shifted=reml_scan(x,noise[:,:3]+50*sf*x[:,1,None],rest[:3],ev,dof)
            assert np.allclose(obj[:,:3],shifted,rtol=1e-9,atol=1e-7)
            assert np.allclose(check['beta'],est['beta'][:3]+50*sf,rtol=1e-7,atol=1e-15)
            boundary.append(dict(tau=tau,toff=toff,true_a=amplitude,median_a=float(np.median(est['a'])),
                                 zero_fraction=float(np.mean(est['index']==0)),upper_fraction=float(np.mean(est['index']==len(GRID)-1))))
            for snr in [0.,-2.,2.,-5.,5.,-20.,20.,-50.,50.]:
                truth=snr*sf; estimate=est['beta']+truth
                u=coverage.wilson(int(coverage.mass_coverage(truth,estimate,est['sigma_beta']).sum()),nsim)
                signed=coverage.wilson(int((np.abs(estimate-truth)<=coverage.Z95*est['sigma_beta']).sum()),nsim)
                rows.append(dict(tau=tau,toff=toff,true_a=amplitude,snr_full=snr,beta_true=truth,U_K1=u,signed_K1=signed))
                if snr!=0: worst=min(worst,u['fraction'])
        print('REML completed',tau,round(toff,4),flush=True)
    threshold=.95-5*math.sqrt(.95*.05/nsim)
    result=dict(seed=SEED,nsim=nsim,a_grid=GRID.tolist(),profile_dof=dof,scope='Estimated covariance within the prespecified fixed Fourier family; not an astrophysical noise validation.',
                min_U_coverage=worst,screen_threshold=threshold,screen_pass=worst>=threshold,rows=rows,amplitude_estimates=boundary,stored_residual_fits=fits)
    (audit.OUT/'estimated-covariance-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Min estimated-covariance U coverage:',worst,'screen:',result['screen_pass'])
    print('Upper-grid boundary fraction:',max(r['upper_fraction'] for r in boundary))


if __name__=='__main__':
    main()
