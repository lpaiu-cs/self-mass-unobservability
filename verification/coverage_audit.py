"""Request 11.2: conditional linear-model coverage; no new timing integration."""
import json
import math
from statistics import NormalDist

import numpy as np

import nuisance_audit as audit

NSIM=8192
SEED=2026090902
Z95=NormalDist().inv_cdf(.975)
SCENARIOS=[('white',0.,0.),('omitted_R3',3.,0.),('omitted_R10',10.,0.),
           ('omitted_R30',30.,0.),('extra_fourier_025',0.,.25),('extra_fourier_1',0.,1.)]


def wilson(hits,n):
    p=hits/n; den=1+Z95**2/n
    center=(p+Z95**2/(2*n))/den
    half=Z95*math.sqrt(p*(1-p)/n+Z95**2/(4*n*n))/den
    return dict(hits=int(hits),n=n,fraction=p,lo95=max(0.,center-half),hi95=min(1.,center+half))


def mass_coverage(truth,estimate,sigma):
    if abs(truth)<=Z95*np.min(sigma):
        return np.ones_like(estimate,dtype=bool)
    mass=audit.ndtr((abs(truth)-estimate)/sigma)-audit.ndtr((-abs(truth)-estimate)/sigma)
    return mass<=.95+1e-13


def beta_weight(x,precision):
    c,b=x[:,0],x[:,1]
    residual=b-c*(c@precision@b)/(c@precision@c)
    information=residual@precision@residual
    assert information>0
    l=precision@residual/information
    assert np.allclose(l@x,[0,1],atol=1e-7)
    return l,1/math.sqrt(information)


def compressed_model(g):
    inp=g['inp']; t=inp['t']; modes=[]
    for j in range(31,61):
        phase=2*np.pi*j*(t-t.min())/np.ptp(t)
        modes.extend([np.cos(phase)/j**2,np.sin(phase)/j**2])
    red=inp['sw'][:,None]*np.column_stack(modes)
    red*=math.sqrt(inp['N'])/np.linalg.norm(red)
    candidates=np.column_stack([g['c6'],red])
    residual=audit.sc.proj_out(g['qf'],candidates)
    normalized=residual/np.linalg.norm(residual,axis=0)
    h,s,_=np.linalg.svd(normalized,full_matrices=False)
    h=h[:,s>s[0]*1e-12]
    basis=np.column_stack([g['d'],h]); dim_d=g['d'].shape[1]
    coordinates=basis.T@candidates
    target=audit.sc.proj_out(g['qt'],candidates)
    relative=float(np.linalg.norm(target-basis@coordinates)/np.linalg.norm(target))
    orth=float(np.linalg.norm(basis.T@basis-np.eye(basis.shape[1])))
    assert max(relative,orth)<1e-7,(relative,orth)
    rest_dof=inp['N']-g['qt'].shape[1]-basis.shape[1]
    assert rest_dof>0
    # Direct-vs-compressed checks include a deterministic mean and residual scale.
    rng=np.random.default_rng(SEED-1)
    v=rng.normal(size=inp['N'])
    vp=audit.sc.proj_out(g['qt'],v)
    z=basis.T@vp; rest=vp-basis@z
    assert np.isclose(vp@vp,z@z+rest@rest,rtol=1e-10)
    wc,wb=audit.sc.template_stencil(2.,0.,inp['OMS'])
    w=np.column_stack([wc,wb]); x=coordinates[:,:6]@w
    for start,q in [(0,g['qt']),(dim_d,g['qf'])]:
        xdirect=audit.sc.proj_out(q,g['c6']@w)
        c,b=xdirect[:,0],xdirect[:,1]; res=b-c*(c@b)/(c@c)
        ldirect=res/(res@res)
        l,_=beta_weight(x[start:],np.eye(basis.shape[1]-start))
        assert np.isclose(ldirect@v,l@z[start:],rtol=1e-6,atol=1e-17)
        direct=audit.sc.proj_out(q,v)
        assert np.isclose(direct@direct,z[start:]@z[start:]+rest@rest,rtol=1e-9)
    return basis,coordinates[:,:6],coordinates[:,6:],dim_d,rest_dof,


def envelope(truth,mean,noise,rest,qrank,start,c6,stencils,n_toas,n=512):
    yy=noise[start:,:n]+mean[start:,None]
    scale=np.sqrt((np.sum(yy*yy,axis=0)+rest[:n])/(n_toas-qrank))
    c=c6[start:]; gram=c.T@c; wc,wb=stencils
    a=np.einsum('ki,ij,kj->k',wc,gram,wc)
    b=np.einsum('ki,ij,kj->k',wc,gram,wb)
    info=np.einsum('ki,ij,kj->k',wb,gram,wb)-b*b/a
    assert np.all(info>0)
    weights=(wb-wc*(b/a)[:,None])/info[:,None]
    score=c.T@yy; sigma=1/np.sqrt(info)
    results={}
    for k in [1.,10.]:
        covered=np.max(Z95*k*sigma)*scale>=abs(truth)
        for lo in range(0,len(wc),128):
            remaining=np.flatnonzero(~covered)
            if not len(remaining): break
            sl=slice(lo,min(lo+128,len(wc)))
            estimates=weights[sl]@score[:,remaining]
            sig=k*sigma[sl,None]*scale[None,remaining]
            masses=audit.ndtr((abs(truth)-estimates)/sig)-audit.ndtr((-abs(truth)-estimates)/sig)
            covered[remaining]=np.any(masses<=.95+1e-13,axis=0)
        results[str(k)]=wilson(int(covered.sum()),n)
    return results


def main():
    g=audit.geometry(); basis,c6,red,dim_d,rest_dof=compressed_model(g)
    n=g['inp']['N']; m=basis.shape[1]; rng=np.random.default_rng(SEED)
    white=rng.normal(size=(m,NSIM)); red_z=rng.normal(size=(red.shape[1],NSIM))
    rest=rng.chisquare(rest_dof,size=NSIM)
    old=json.loads((audit.OUT/'nuisance-audit.json').read_text())
    points=sorted(set([(t,0.) for t in audit.LAGS]+[(r['tau'],r['toff']) for r in old['intervals'] if r['cut'] in [0.,1e-3] and r['K']==10]))
    rows=[]; envelopes=[]; controls=[]
    origins=np.arange(0.,g['inp']['P_out'],g['inp']['P_in']/24)
    scenario_cache={}
    for name,radius,amplitude in SCENARIOS:
        noise=white+amplitude*(red@red_z)
        entries=[]
        for method,start in [('truncated_diag',0),('full_diag',dim_d),('full_matched_GLS',dim_d)]:
            cov=np.eye(m-start)
            if method=='full_matched_GLS': cov+=amplitude**2*(red[start:]@red[start:].T)
            precision=np.linalg.inv(cov)
            pn=precision@noise[start:]
            quadratic=np.sum(noise[start:]*pn,axis=0)+rest
            entries.append((method,start,precision,pn,quadratic))
        scenario_cache[name]=(noise,entries)
    for tau,toff in points:
        wc,wb=audit.sc.template_stencil(tau,toff,g['inp']['OMS']); x=c6@np.column_stack([wc,wb])
        lt,_=beta_weight(x,np.eye(m)); lf,sf=beta_weight(x[dim_d:],np.eye(m-dim_d))
        direction=lt[:dim_d]/np.linalg.norm(lt[:dim_d])
        full_origin=next(r['toff'] for r in old['intervals'] if r['cut']==0 and r['K']==10 and r['tau']==tau)
        for name,radius,amplitude in SCENARIOS:
            noise,entries=scenario_cache[name]
            for snr in [0.,-2.,2.,-5.,5.,-20.,20.,-50.,50.]:
                truth=snr*sf; mean=truth*x[:,1].copy()
                mean[:dim_d]-=(1 if snr>=0 else -1)*radius*direction
                for method,start,precision,pn,quadratic in entries:
                    l,sigma=beta_weight(x[start:],precision); mu=mean[start:]
                    estimate=l@noise[start:]+l@mu
                    rank=g['qt'].shape[1] if start==0 else g['qf'].shape[1]
                    scale=np.sqrt((quadratic+2*mu@pn+mu@precision@mu)/(n-rank))
                    row=dict(tau=tau,toff=toff,scenario=name,method=method,snr_full=snr,beta_true=truth,
                             sigma_assumed_unit=sigma,bias_in_assumed_sigma=float((l@mu-truth)/sigma),
                             error_sd_in_assumed_sigma=float(np.std(estimate,ddof=1)/sigma),
                             median_scale=float(np.median(scale)))
                    for k in [1.,10.]:
                        sig=k*sigma*scale
                        row[f'U_K{k:g}']=wilson(int(mass_coverage(truth,estimate,sig).sum()),NSIM)
                        row[f'signed_K{k:g}']=wilson(int((np.abs(estimate-truth)<=Z95*sig).sum()),NSIM)
                    known=wilson(int(mass_coverage(truth,estimate,np.full(NSIM,sigma)).sum()),NSIM)
                    row['U_known_scale_K1']=known; rows.append(row)
                    if method=='full_matched_GLS' and snr!=0:
                        controls.append(known['fraction'])
                if tau in [2.,200.] and toff==full_origin and name in ['white','omitted_R30','extra_fourier_1'] and snr==50:
                    stencils=audit.stencil_grid(tau,origins,g['inp']['OMS'])
                    for method,start,rank in [('truncated_diag',0,g['qt'].shape[1]),('full_diag',dim_d,g['qf'].shape[1])]:
                        value=envelope(truth,mean,noise,rest,rank,start,c6,stencils,n)
                        envelopes.append(dict(tau=tau,toff=toff,scenario=name,method=method,beta_true=truth,coverage=value))
        print('completed tau',tau,'origin',round(toff,4),flush=True)
    threshold=.95-5*math.sqrt(.95*.05/NSIM)
    result=dict(seed=SEED,nsim=NSIM,scope='Conditional coverage on frozen linear arrays; known-covariance GLS is an oracle control, not fitted astrophysical noise.',
                compressed_dimension=m,omitted_dimension=dim_d,residual_chi2_dof=rest_dof,
                gaussian_control_min_coverage=min(controls),gaussian_control_screen_threshold=threshold,
                gaussian_control_screen_pass=min(controls)>=threshold,rows=rows,registered_grid_envelopes=envelopes)
    (audit.OUT/'coverage-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Gaussian positive-control min:',min(controls),'screen threshold:',threshold)
    print('Wrote',len(rows),'pointwise rows and',len(envelopes),'full-grid envelope rows')
    assert result['gaussian_control_screen_pass'], 'Known-Gaussian control failed; inspect without retuning.'


if __name__=='__main__':
    main()
