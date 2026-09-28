"""Corrected leading coplanar drive, explicit normalization, and fixed-array fits."""
import json
import math
import numpy as np

import comparator_audit as comp
import estimated_covariance_audit as est
import nuisance_audit as audit


def drive():
    z = np.load(audit.ROOT/'request10_external/baseline_planetGR.npz', allow_pickle=True)
    p = dict(zip(z['names'], z['params']))
    gm = 6.67408e-11*1.9884754153381438e30
    a_b = math.hypot(p['absini_o'], p['abcosi_o'])
    inc_o = math.acos(p['abcosi_o']/a_b)
    a_p = p['apsini_i']/math.sin(inc_o+math.radians(p['delta_i']))
    q = p['masspar_p']; a_i = a_p*(1+1/q)
    m_b = a_i**3*(2*math.pi/(p['period_i']*86400))**2/gm
    m_p = m_b/(1+q); m_i = q*m_p
    mass_function = 4*math.pi**2*a_b**3/(p['period_o']*86400)**2/gm
    m_o = .4
    for _ in range(100):
        m_o = (mass_function*(m_o+m_b)**2)**(1/3)
    a_o = a_b*(m_o+m_b)/m_o
    f = m_i/m_b
    e = np.array([math.hypot(p['eta_p'], p['kappa_p']), math.hypot(p['eta_b'], p['kappa_b'])])
    varpi = np.array([math.atan2(p['kappa_p'], p['eta_p']), math.atan2(p['kappa_b'], p['eta_b'])])
    potential_factors = gm/299792458.**2*np.array([m_i/a_i, m_o/a_o])
    amplitudes = np.array([potential_factors[0]*e[0], potential_factors[1]*e[1], potential_factors[1]*f*a_i/a_o])
    ustar = amplitudes.sum()
    matching = json.loads((audit.OUT/'physical-matching.json').read_text())
    phi = np.array([matching['leading_physical_phase_radians'][k] for k in ['in','out','dif']])
    return dict(masses_solar=[m_p,m_i,m_o], semi_major_m=[a_i,a_o], f=f,
                eccentricities=e.tolist(), varpi=varpi.tolist(),
                potential_factors=potential_factors.tolist(), amplitudes=amplitudes.tolist(),
                Ustar=float(ustar), normalized_amplitudes=(amplitudes/ustar).tolist(), phases=phi.tolist(),
                mass_o_parameter_difference=float(m_o-p['mass_o']))


def stencil(g, d, tau):
    om = np.array([g['inp']['OMS'][k] for k in ['in','out','dif']])
    scale = np.repeat(d['normalized_amplitudes'], 2)
    wc = scale*comp.weights(om, d['phases'], [0])[:,0]
    wb = scale*comp.pole(om, d['phases'], tau)
    return np.column_stack([wc,wb])


def orbital_grid(d, n):
    mean = np.arange(n)*2*np.pi/n
    vectors = []
    for e, angle, a in zip(d['eccentricities'], d['varpi'], d['semi_major_m']):
        eccentric = mean.copy()
        for _ in range(8):
            eccentric -= (eccentric-e*np.sin(eccentric)-mean)/(1-e*np.cos(eccentric))
        assert np.max(np.abs(eccentric-e*np.sin(eccentric)-mean)) < 2e-14
        xy = a*np.column_stack([np.cos(eccentric)-e, math.sqrt(1-e*e)*np.sin(eccentric)])
        rotation = np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])
        vectors.append(xy@rotation.T)
    ri, ro = vectors
    gm_c2 = 6.67408e-11*1.9884754153381438e30/299792458.**2
    exact = gm_c2*(d['masses_solar'][1]/np.linalg.norm(ri,axis=1)[:,None]
                  +d['masses_solar'][2]/np.linalg.norm(ro[None,:,:]+d['f']*ri[:,None,:],axis=2))
    exact -= exact.mean()
    amp = d['amplitudes']; vp = d['varpi']
    leading = amp[0]*np.cos(mean)[:,None]+amp[1]*np.cos(mean)[None,:]-amp[2]*np.cos(mean[:,None]-mean[None,:]+vp[0]-vp[1])
    error = exact-leading
    return dict(n=n, omitted_input_rms_over_leading_rms=float(np.linalg.norm(error)/np.linalg.norm(leading)),
                sampled_max_omitted_over_Ustar=float(np.max(np.abs(error))/d['Ustar']),
                exact_rms_over_Ustar=float(np.sqrt(np.mean(exact*exact))/d['Ustar']))


def main():
    d = drive(); g = comp.prepare(); rows=[]
    assert np.isclose(sum(d['normalized_amplitudes']),1)
    for tau in audit.LAGS:
        w = stencil(g,d,tau); x=g['c']@w
        fitted,_=est.reml_scan(x,g['y'],g['y_rest'],g['eigenvalues'],g['inp']['N']-90-2)
        rows.append(dict(tau=tau,a=float(fitted['a'][0]),beta=float(fitted['beta'][0]),
                         sigma=float(fitted['sigma_beta'][0]),
                         conditional_pointwise_U=float(audit.intervals(fitted['beta'][0],fitted['sigma_beta'][0]))))
    result=dict(status='Imported from prior work',scope='New computation on frozen arrays; leading coplanar Newtonian drive, no full scalar-tensor force or EOS matching',
                drive=d,torus_checks=[orbital_grid(d,n) for n in [128,256]],fits=rows,
                controls=['runtime parameter_set=6 mass/length convention','Kepler equation residual','Ustar normalization'],
                omitted_timing_error_bound=None)
    (audit.OUT/'corrected-physical-drive.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
