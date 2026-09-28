"""Request 11.4: bounded comparators on frozen response arrays, no timing engine."""
import json
import math
import numpy as np
import sympy as sp

import nuisance_audit as audit
import coverage_audit as coverage

AMPLITUDES=[0.,0.08316109496333533,1.]


def prepare():
    g=audit.geometry()
    basis,c6,red,nd,_=coverage.compressed_model(g)
    h=basis[:,nd:]; eig,rot=np.linalg.eigh(red[nd:]@red[nd:].T)
    eig=np.maximum(eig,0)
    yp=audit.sc.proj_out(g['qf'],g['inp']['sw']*g['inp']['res0'])
    g.update(h=h@rot, eigenvalues=eig, c=rot.T@c6[nd:],y=rot.T@(h.T@yp),yp=yp)
    g['y_rest']=float(yp@yp-g['y']@g['y'])
    assert g['y_rest']>0
    return g


def metric(g,amplitude):
    root=1/np.sqrt(1+amplitude**2*g['eigenvalues'])
    c=root[:,None]*g['c']; y=root*g['y']
    u,s,vh=np.linalg.svd(c,full_matrices=False)
    assert s[-1]>s[0]*1e-12
    return dict(c=c,y=y,root=root,u=u,s=s,vh=vh,gram=c.T@c,score=c.T@y,
                signal_data_norm=float(np.linalg.norm(u.T@y)))


def weights(omegas,phases,powers):
    n=np.asarray(powers); angle=np.asarray(phases)[:,None]+np.pi*n[None,:]/2
    amplitude=(omegas[:,None]/max(omegas))**n[None,:]
    return np.stack([amplitude*np.cos(angle),amplitude*np.sin(angle)],axis=1).reshape(6,len(n))


def pole(omegas,phases,tau):
    h=np.exp(1j*np.asarray(phases))/(1+1j*omegas*tau)
    return np.column_stack([h.real,h.imag]).ravel()


def residual(target,columns):
    normalized=columns/np.linalg.norm(columns,axis=0)
    u,s,_=np.linalg.svd(normalized,full_matrices=False)
    assert s[-1]>s[0]*1e-12
    return target-u@(u.T@target)


def algebra_checks():
    t=sp.symbols('t',real=True)
    f=sp.sin(t)+sp.cos(2*t)
    for n in [1,3,5]:
        assert sp.integrate(sp.expand_trig(f*sp.diff(f,t,n)),(t,0,2*sp.pi))==0
    tau,lo,hi=sp.symbols('tau lo hi',positive=True)
    ratio=(1+hi**2*tau**2)/(1+lo**2*tau**2)
    assert sp.simplify(sp.diff(ratio,tau)-2*tau*(hi**2-lo**2)/(1+lo**2*tau**2)**2)==0
    # Reciprocal gradient-system transfer, including noncommuting positive matrices.
    damping=np.array([[2.,.2],[.2,1.]]); stiffness=np.array([[2.,1.],[1.,3.]])
    b=np.array([1.,.3]); ev,vec=np.linalg.eigh(damping)
    invroot=(vec/np.sqrt(ev))@vec.T
    rates,modes=np.linalg.eigh(invroot@stiffness@invroot)
    residues=(modes.T@invroot@b)**2/rates
    assert min(rates)>0 and min(residues)>=0
    for w in [.01,1.,10.]:
        direct=b@np.linalg.solve(stiffness+1j*w*damping,b)
        spectral=np.sum(residues/(1+1j*w/rates))
        assert abs(direct-spectral)<1e-13
    return ['odd_conservative_terms','monotone_quadrature_ratio','reciprocal_positive_relaxation_spectrum']


def main():
    checks=algebra_checks(); g=prepare()
    om=np.array([g['inp']['OMS'][k] for k in ['in','out','dif']])
    matching=json.loads((audit.OUT/'physical-matching.json').read_text())
    phases=[('zero',[0.,0.,0.]),('independent',[.7,1.1,2.2]),
            ('physical_phase_only',[matching['leading_physical_phase_radians'][k] for k in ['in','out','dif']])]
    gap=10*max(om); low,high=1,0
    fast_ratio=(1+(om[high]/gap)**2)/(1+(om[low]/gap)**2)
    rows=[]; witnesses=[]; diagnostics=[]
    for amplitude in AMPLITUDES:
        m=metric(g,amplitude); c=m['c']
        direct=audit.sc.proj_out(g['qf'],g['c6'])
        direct+=g['h']@((m['root']-1)[:,None]*(g['h'].T@direct))
        assert np.allclose(direct.T@direct,m['gram'],rtol=1e-7,atol=1e-3)
        diagnostics.append(dict(a=amplitude,map_singular_values=m['s'].tolist(),condition=float(m['s'][0]/m['s'][-1])))
        for tau in audit.LAGS:
            for label,phi in phases:
                wp=pole(om,phi,tau); target=c@wp
                base=np.linalg.norm(residual(target,c@weights(om,phi,[0])))**2
                for name,powers in [(f'P{n}',list(range(n+1))) for n in range(6)]+[('even_P4',[0,2,4])]:
                    v=weights(om,phi,powers)
                    rp=residual(target,c@v); information=float(rp@rp)
                    rc=residual(wp,v); coefficient_distance=float(np.linalg.norm(rc))
                    lower=float(m['s'][-1]**2*(rc@rc))
                    q,_=np.linalg.qr(c@v,mode='reduced')
                    qr=target-q@(q.T@target)
                    assert np.linalg.norm(rp-qr)<1e-7*np.linalg.norm(target)
                    if name=='P5':
                        assert np.linalg.norm(rp)/np.linalg.norm(target)<1e-9
                    else:
                        assert information>=lower*(1-1e-6) and information>0
                    rows.append(dict(a=amplitude,tau=tau,phase=label,comparator=name,
                                     relative_information=information/base,unit_sigma_beta=None if name=='P5' else 1/math.sqrt(information),
                                     continuous_phase_information_lower_bound=0. if name=='P5' else lower,
                                     coefficient_residual_norm=coefficient_distance,exact_collapse=name=='P5'))
                # S=-Im(H)/omega; dephase the fitted complex carrier coefficients first.
                v=np.zeros(6)
                for k,factor in [(low,1.),(high,-fast_ratio)]:
                    v[2*k:2*k+2]=factor*np.array([math.sin(phi[k]),-math.cos(phi[k])])/om[k]
                violation=float(v@wp)
                analytic=tau/(1+(om[low]*tau)**2)-fast_ratio*tau/(1+(om[high]*tau)**2)
                assert np.isclose(violation,analytic,rtol=1e-12) and violation>0
                assert abs(v@weights(om,phi,[0])[:,0])<1e-12
                sigma=float(np.linalg.norm((m['vh']@v)/m['s']))
                witnesses.append(dict(a=amplitude,tau=tau,phase=label,fast_gap_rate_per_day=gap,
                                      slow_positive_beta_witness=violation,unit_noise_witness_sigma=sigma,
                                      minimum_distance_per_positive_beta=violation/sigma))
    # The unprojected coefficient distance is unchanged by a block-orthogonal phase rotation.
    for tau in audit.LAGS:
        for name in ['P0','P1','P2','P3','P4','even_P4']:
            values=[r['coefficient_residual_norm'] for r in rows if r['tau']==tau and r['comparator']==name]
            assert np.ptp(values)<1e-8*max(values)
    checks+=['direct_whitened_metric','comparator_QR_SVD','degree_five_collapse','continuous_phase_information_bound','quadrature_witness']
    result=dict(status='Imported from prior work',scope='Conditional frozen-array identifiability; known phases and hypothetical fast relaxation gap; no detection or EOS constraint',
                checks=checks,diagnostics=diagnostics,rows=rows,positive_fast_spectrum_witnesses=witnesses,
                gap_is_physically_matched=False)
    (audit.OUT/'comparator-audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('PASS:',', '.join(checks))
    for name in ['P0','P1','P2','P3','P4','P5','even_P4']:
        values=[r['relative_information'] for r in rows if r['comparator']==name]
        print(name,'relative information min/max',min(values),max(values))
    print('Map conditions:',[r['condition'] for r in diagnostics])


if __name__=='__main__':
    main()
