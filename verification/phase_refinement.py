"""Registered local follow-through; never represents a global optimization certificate."""
import json
import math
import numpy as np
import nuisance_audit as audit
import comparator_audit as comparator


def main():
    g=comparator.prepare(); om=np.array([g['inp']['OMS'][k] for k in ['in','out','dif']])
    previous=json.loads((audit.OUT/'phase-state-audit.json').read_text())['phase_rows']
    offsets=np.array(np.meshgrid([-1.,0.,1.],[-1.,0.,1.],[-1.,0.,1.],indexing='ij')).reshape(3,-1).T
    result=[]
    for amplitude in comparator.AMPLITUDES:
        m=comparator.metric(g,amplitude)
        scale=math.sqrt((m['y']@m['y']+g['y_rest'])/(g['inp']['N']-90))
        for tau in audit.LAGS:
            seeds=[r for r in previous if r['a']==amplitude and r['tau']==tau]
            def values(phi):
                f=np.exp(1j*phi); h=f/(1+1j*tau*om[None,:])
                wc=np.stack([f.real,f.imag],axis=2).reshape(-1,6)
                wp=np.stack([h.real,h.imag],axis=2).reshape(-1,6)
                beta,sigma=audit.fit_grid(m['gram'],m['score'],wc,wp,scale)
                return audit.intervals(beta,sigma)
            runs=[]
            for seed in seeds:
                point=np.array(seed['maximizing_phases']); best=seed['grid_max_U']; caps=[]
                for level in range(17):
                    before=best; step=(math.pi/12)*2**(-level)
                    for move in range(24):
                        trials=(point+step*offsets)%(2*math.pi); u=values(trials); j=int(np.argmax(u))
                        if u[j]<=best*(1+1e-13): break
                        point,best=trials[j],float(u[j])
                    else:
                        caps.append(level)
                assert best>=seed['grid_max_U']*(1-1e-12)
                runs.append(dict(seed_domain=seed['domain'],U=best,phases=point.tolist(),
                                 last_level_relative_improvement=(best-before)/best,move_cap_levels=caps))
            best=max(runs,key=lambda r:r['U']); upper=seeds[0]['continuous_all_phase_upper_bound']
            assert best['U']<=upper*(1+1e-8)
            legacy=next(r['grid_max_U'] for r in seeds if r['domain']=='legacy_origins')
            result.append(dict(a=amplitude,tau=tau,legacy_U=legacy,refined_U=best['U'],
                               refined_over_legacy=best['U']/legacy,continuous_upper_bound=upper,runs=runs,
                               global_maximum_certified=False))
            print('refined',amplitude,tau,best['U'],'ratio',best['U']/legacy,flush=True)
    out=dict(scope='Post-coarse-result registered local refinement; all seeds retained, analytic upper envelope unchanged',
             plan='REQUEST11_5B_PHASE_REFINEMENT_PLAN.md',rows=result)
    (audit.OUT/'phase-refinement.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('move-cap hits',sum(len(s['move_cap_levels']) for r in result for s in r['runs']))
    print('max last-level relative improvement',max(s['last_level_relative_improvement'] for r in result for s in r['runs']))


if __name__=='__main__':
    main()
