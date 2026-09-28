"""Request 11.5: expanded phases, missing transient response and linear gap steps."""
import hashlib
import json
import math
from statistics import NormalDist
import numpy as np
import sympy as sp

import nuisance_audit as audit
import comparator_audit as comparator

Z95=NormalDist().inv_cdf(.975)


def phase_rows(g):
    om=np.array([g['inp']['OMS'][k] for k in ['in','out','dif']])
    legacy=np.arange(0.,g['inp']['P_out'],g['inp']['P_in']/24)[:,None]*om
    pair=np.array(np.meshgrid(np.arange(72)*2*np.pi/72,np.arange(72)*2*np.pi/72,indexing='ij')).reshape(2,-1).T
    two=np.column_stack([pair,pair[:,0]-pair[:,1]])
    three=np.array(np.meshgrid(*([np.arange(24)*2*np.pi/24]*3),indexing='ij')).reshape(3,-1).T
    grids=[('legacy_origins',legacy),('two_longitudes',two),('three_phases',three)]
    rows=[]
    prior=json.loads((audit.OUT/'nuisance-audit.json').read_text())
    for amplitude in comparator.AMPLITUDES:
        m=comparator.metric(g,amplitude)
        scale=math.sqrt((m['y']@m['y']+g['y_rest'])/(g['inp']['N']-90))
        for tau in audit.LAGS:
            v=comparator.weights(om,[0.,0.,0.],[0]); wb=comparator.pole(om,[0.,0.,0.],tau)
            distance=np.linalg.norm(comparator.residual(wb,v))
            information_floor=m['s'][-1]**2*distance**2
            upper=(m['signal_data_norm']+Z95*scale)/math.sqrt(information_floor)
            for name,phi in grids:
                f=np.exp(1j*phi); h=f/(1+1j*tau*om[None,:])
                wc=np.stack([f.real,f.imag],axis=2).reshape(-1,6)
                wp=np.stack([h.real,h.imag],axis=2).reshape(-1,6)
                beta,sigma=audit.fit_grid(m['gram'],m['score'],wc,wp,scale)
                information=(scale/sigma)**2
                assert min(information)>=information_floor*(1-1e-8)
                u=audit.intervals(beta,sigma); j=int(np.argmax(u))
                assert max(u)<=upper*(1+1e-8)
                # Direct SVD fit at the observed extremum checks the six-dimensional Gram shortcut.
                x=m['c']@np.column_stack([wc[j],wp[j]])
                fit=np.linalg.lstsq(x,m['y'],rcond=None)[0]
                rb=comparator.residual(x[:,1],x[:,0,None]); info=rb@rb
                assert np.isclose(beta[j],fit[1],rtol=1e-6,atol=1e-16)
                assert np.isclose(sigma[j],scale/math.sqrt(info),rtol=1e-7)
                if amplitude==0 and name=='legacy_origins':
                    old=next(r['U'] for r in prior['intervals'] if r['tau']==tau and r['cut']==0 and r['K']==1)
                    assert np.isclose(max(u),old,rtol=1e-7)
                rows.append(dict(a=amplitude,tau=tau,domain=name,grid_size=len(phi),
                                 grid_max_U=float(u[j]),maximizing_phases=phi[j].tolist(),beta=float(beta[j]),sigma=float(sigma[j]),
                                 continuous_all_phase_upper_bound=upper,information_lower_bound=information_floor,
                                 data_scale=scale,signal_space_data_norm=m['signal_data_norm']))
    return rows


def transient_rows(g):
    t=g['inp']['t']; x=t-min(t)
    om=np.array([g['inp']['OMS'][k] for k in ['in','out','dif']])
    periodic=np.column_stack([np.ones_like(t)]+[fun(w*t) for w in om for fun in [np.cos,np.sin]])
    q,_=np.linalg.qr(periodic,mode='reduced')
    rows=[]
    for tau in audit.LAGS:
        h=np.exp(-x/tau); hp=h-q@(q.T@h)
        fraction=float(np.linalg.norm(hp)/np.linalg.norm(h))
        assert fraction>1e-3
        rows.append(dict(tau=tau,transient_outside_known_drive_span=fraction,
                         final_over_initial=math.exp(-np.ptp(t)/tau),
                         time_to_one_percent_days=tau*math.log(100),
                         time_to_one_permille_days=tau*math.log(1000)))
    # A causal local differential operator annihilates every recorded periodic/static input.
    z,tau,w1,w2,w3=sp.symbols('z tau w1 w2 w3',positive=True)
    annihilator=z*(z*z+w1*w1)*(z*z+w2*w2)*(z*z+w3*w3)
    assert annihilator.subs(z,0)==0
    assert all(sp.simplify(annihilator.subs(z,sgn*sp.I*w))==0 for w in [w1,w2,w3] for sgn in [-1,1])
    transient=sp.factor(annihilator.subs(z,-1/tau))
    assert sp.simplify(transient+(1/tau)*(1/tau**2+w1**2)*(1/tau**2+w2**2)*(1/tau**2+w3**2))==0
    return dict(rows=rows,causal_annihilator=str(annihilator),transient_response=str(transient),
                conclusion='Six harmonic responses plus a static response do not determine the exponential response without an additional operator restriction or transient column.')


def gap_rows(g):
    t=g['inp']['t']; order=np.argsort(t); sorted_t=t[order]
    gaps=np.flatnonzero(np.diff(sorted_t)>1.)
    assert len(gaps)==565
    frozen=np.load(audit.ROOT/'request10_external/baseline_planetGR.npz',allow_pickle=True)
    p=dict(zip(frozen['names'],frozen['params']))
    cycle_us=86400e6/p['spinfreq']
    rows=[]
    for amplitude in [0.,1.]:
        m=comparator.metric(g,amplitude); q6=g['h']@m['u']
        data=g['yp']+g['h']@((m['root']-1)*g['y'])
        data-=q6@(q6.T@data)
        for lo in range(0,len(gaps),64):
            indices=gaps[lo:lo+64]
            cuts=(sorted_t[indices]+sorted_t[indices+1])/2
            raw=cycle_us*g['inp']['sw'][:,None]*(t[:,None]>cuts[None,:])
            pj=audit.sc.proj_out(g['qf'],raw); hc=g['h'].T@pj
            whitened=pj+g['h']@((m['root']-1)[:,None]*hc)
            coefficient=q6.T@whitened
            remain=whitened-q6@coefficient
            norm2=np.sum(remain**2,axis=0); cross=data@remain
            compressed=np.sum(pj**2,axis=0)+np.sum(((m['root']**2-1)[:,None])*hc**2,axis=0)-np.sum(coefficient**2,axis=0)
            assert np.allclose(norm2,compressed,rtol=1e-7,atol=1e-4)
            for j,index in enumerate(indices):
                rows.append(dict(a=amplitude,after_sorted_index=int(index),left_day=float(sorted_t[index]),right_day=float(sorted_t[index+1]),
                                 gap_days=float(sorted_t[index+1]-sorted_t[index]),remaining_unit_noise_norm=math.sqrt(norm2[j]),
                                 delta_chi2_plus=float(norm2[j]+2*cross[j]),delta_chi2_minus=float(norm2[j]-2*cross[j])))
    return dict(gap_threshold_days=1.,gaps=565,cycle_microseconds=cycle_us,rows=rows,
                scope='Idealized linear permanent single-cycle steps only; full nuisance and all six free harmonics fitted. Not arbitrary multiple slips or nonlinear pulse-number certification.')


def numerical_audit(g):
    old=audit.ROOT/'request10_external/finite_jacobian.npy'
    j1=np.load(old); j2=g['inp']['J']
    assert j1.shape==j2.shape
    b1=g['b'].copy(); b1[:,:28]=g['inp']['sw'][:,None]*j1/g['norms'][:28]
    q1,s1,_=np.linalg.svd(b1,full_matrices=False)
    cosines=np.linalg.svd(q1.T@g['qf'],compute_uv=False)
    q2,_=np.linalg.qr(g['b'],mode='reduced')
    discrepancy=float(np.linalg.norm(audit.sc.proj_out(g['qf'],q2)))
    assert discrepancy<1e-7
    meta=json.loads((audit.ROOT/'request10_external/finite_jacobian_v2_meta.json').read_text())
    oldmeta=json.loads((audit.ROOT/'request10_external/finite_jacobian_meta.json').read_text())
    assert oldmeta['columns']==meta['columns'], 'Derivative columns must have the same order before comparison.'
    differences=[]
    for name in meta['v2_replaced_columns']:
        j=meta['columns'].index(name)
        rel=np.linalg.norm(g['inp']['sw']*(j1[:,j]-j2[:,j]))/np.linalg.norm(g['inp']['sw']*j2[:,j])
        differences.append(dict(column=name,v1_v2_weighted_relative_difference=float(rel),
                                stored_halfstep_deviation=meta['v2_replaced_columns'][name]['halfstep_lin_dev']))
    return dict(qr_svd_subspace_residual=discrepancy,v1_min_relative_singular_value=float(s1[-1]/s1[0]),
                v1_v2_max_principal_sine=math.sqrt(max(0.,1-min(cosines)**2)),differences=differences,
                v1_sha256=hashlib.sha256(old.read_bytes()).hexdigest(),
                full_derivative_error_certificate=False,
                scope='v1 uses superseded coarse planet steps and is diagnostic only. Stored half-step scalars do not identify their projected error vectors or bound all timing derivatives.')


def main():
    g=comparator.prepare()
    phases=phase_rows(g); transients=transient_rows(g); gaps=gap_rows(g); numeric=numerical_audit(g)
    result=dict(status='Imported from prior work',scope='Fixed linear arrays, declared phase grids and analytic bounds, missing transient response, idealized single-cycle gaps; no timing engine',
                phase_rows=phases,transients=transients,gaps=gaps,numerical=numeric)
    (audit.OUT/'phase-state-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for a in comparator.AMPLITUDES:
        print('Covariance a=',a)
        for tau in [2.,200.,500.]:
            r=[v for v in phases if v['a']==a and v['tau']==tau]
            print('phase U',tau,[(v['domain'],v['grid_max_U']) for v in r],'analytic bound',r[0]['continuous_all_phase_upper_bound'])
    for a in [0.,1.]:
        r=[v for v in gaps['rows'] if v['a']==a]
        print('gap min norm / min delta chi2',a,min(v['remaining_unit_noise_norm'] for v in r),min(min(v['delta_chi2_plus'],v['delta_chi2_minus']) for v in r))
    print('Numerical',numeric)
    print('PASS: phase bound and direct fits; legacy reproduction; causal transient obstruction; all 565 single-step gaps; numerical diagnostics')


if __name__=='__main__':
    main()
