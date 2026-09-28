"""Six-coefficient confidence regions: phase/lag selection follows joint inference."""
import json
import math
import sys
import numpy as np

import comparator_audit as comp
import coverage_audit as coverage
import estimated_covariance_audit as estimated
import nuisance_audit as audit
from physical_drive_completion import drive, stencil

NSIM = 8192
CHI6 = 12.591587243743977
CALPATH = audit.OUT/'simultaneous-calibration.json'


def setup():
    g=comp.prepare(); q,r=np.linalg.qr(g['c'],mode='reduced')
    assert np.linalg.norm(q@r-g['c'])/np.linalg.norm(g['c'])<1e-12
    g.update(q6=q,r6=r,dof=g['inp']['N']-90-6)
    return g


def scan(g,y,rest):
    if y.ndim==1: y=y[:,None]
    best=np.full(y.shape[1],np.inf)
    out=dict(statistic=np.empty(y.shape[1]),a=np.empty(y.shape[1]),scale=np.empty(y.shape[1]),
             coefficient=np.empty((6,y.shape[1])),index=np.empty(y.shape[1],dtype=int))
    matrices=[]
    for j,a in enumerate(estimated.GRID):
        p=1/(1+a*a*g['eigenvalues']); qp=g['q6'].T*p
        gram=qp@g['q6']; score=qp@y
        coefficient=np.linalg.solve(gram,score)
        fitted=np.sum(score*coefficient,axis=0)
        rss=p@(y*y)+rest-fitted
        assert np.all(rss>0)
        s2=rss/g['dof']
        objective=np.log1p(a*a*g['eigenvalues']).sum()+np.linalg.slogdet(gram)[1]+g['dof']*np.log(s2)
        keep=objective<best; best[keep]=objective[keep]
        out['statistic'][keep]=fitted[keep]/s2[keep]
        out['a'][keep]=a; out['scale'][keep]=np.sqrt(s2[keep]); out['index'][keep]=j
        out['coefficient'][:,keep]=coefficient[:,keep]
        matrices.append(gram)
    out['grams']=matrices
    return out


def draws(g,seed,a,misspecified=False):
    rng=np.random.default_rng(seed)
    n=len(g['eigenvalues'])
    if misspecified:
        t=g['inp']['t']; red=[]
        for j in range(31,61):
            angle=2*np.pi*j*(t-t.min())/np.ptp(t)
            red.extend([np.cos(angle)/j,np.sin(angle)/j])
        red=g['inp']['sw'][:,None]*np.column_stack(red)
        red*=math.sqrt(len(t))/np.linalg.norm(red)
        projected=g['h'].T@red
        y=rng.normal(size=(n,NSIM))+a*projected@rng.normal(size=(60,NSIM))
    else:
        y=np.sqrt(1+a*a*g['eigenvalues'])[:,None]*rng.normal(size=(n,NSIM))
    rest=rng.chisquare(g['inp']['N']-90-n,size=NSIM)
    return y,rest


def envelope_stat(g,y):
    p=1/(1+16*g['eigenvalues']); qp=g['q6'].T*p
    gram=qp@g['q6']; score=qp@y
    return np.sum(score*np.linalg.solve(gram,score),axis=0)


def physical_projection(g,fit,threshold):
    d=drive(); j=int(fit['index'][0]); gram=fit['grams'][j]
    mean=fit['coefficient'][:,0]; scale=float(fit['scale'][0]); rows=[]
    for tau in audit.LAGS:
        x=g['r6']@stencil(g,d,tau)
        normal=x.T@gram@x; coeff=np.linalg.solve(normal,x.T@gram@mean)
        delta=mean-x@coeff; qmin=float(delta@gram@delta/scale**2)
        half=math.sqrt(max(0,threshold-qmin)*scale**2*np.linalg.inv(normal)[1,1])
        rows.append(dict(tau=tau,minimum_statistic=qmin,empty=qmin>threshold,
                         beta=float(coeff[1]),joint_region_abs_beta_upper=None if qmin>threshold else float(abs(coeff[1])+half)))
    return rows


def main():
    mode=sys.argv[1]; assert mode in ['calibrate','validate']
    seed=2026090912 if mode=='calibrate' else 2026090913
    g=setup(); rows=[]; quantiles=[]
    threshold=CHI6 if mode=='calibrate' else json.loads(CALPATH.read_text())['threshold']
    scenarios=[(a,False) for a in [0.,.25,1.,4.]]
    if mode=='validate': scenarios.append((1.,True))
    for i,(a,miss) in enumerate(scenarios):
        y,rest=draws(g,seed+100*i,a,miss); out=scan(g,y,rest)
        # Refitting all six signals makes the covariance selection invariant to any signal in their span.
        shift=np.arange(1.,7.)/3
        shifted=scan(g,y[:,:3]+g['q6']@shift[:,None],rest[:3])
        assert np.array_equal(shifted['index'],out['index'][:3])
        assert np.allclose(shifted['coefficient'],out['coefficient'][:,:3]+shift[:,None],rtol=1e-8,atol=1e-8)
        order=math.ceil((NSIM+1)*.95)-1
        quantile=float(np.sort(out['statistic'])[order]); quantiles.append(quantile)
        rows.append(dict(a=a,misspecified=miss,order_statistic_95=quantile,
                         inclusion=coverage.wilson(np.count_nonzero(out['statistic']<=threshold),NSIM),
                         null_false_positive=coverage.wilson(np.count_nonzero(out['statistic']>threshold),NSIM),
                         covariance_upper_boundary=float(np.mean(out['index']==len(estimated.GRID)-1)),
                         covariance_envelope_known_scale_inclusion=coverage.wilson(np.count_nonzero(envelope_stat(g,y)<=CHI6),NSIM)))
        print(mode,a,miss,rows[-1]['inclusion']['fraction'],flush=True)
    result=dict(status='Imported from prior work',seed=seed,nsim=NSIM,rows=rows,
                scope='Six-dimensional region on frozen linear arrays, fixed Fourier family; continuous phase/lag projection inherits region inclusion when true mean is in its span.',
                threshold=max(CHI6,*quantiles) if mode=='calibrate' else threshold,
                nominal_chi6=CHI6,signal_translation_control=True)
    if mode=='validate':
        data=scan(g,g['y'],g['y_rest'])
        result['data']=dict(a=float(data['a'][0]),scale=float(data['scale'][0]),
                            omnibus_null_statistic=float(data['statistic'][0]),
                            exceeds_threshold=bool(data['statistic'][0]>threshold),
                            physical_lag_sections=physical_projection(g,data,threshold))
        result['calibration_sha256']=audit.hashlib.sha256(CALPATH.read_bytes()).hexdigest()
    path=CALPATH if mode=='calibrate' else audit.OUT/'simultaneous-validation.json'
    path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('threshold',result['threshold'])


if __name__=='__main__':
    main()
