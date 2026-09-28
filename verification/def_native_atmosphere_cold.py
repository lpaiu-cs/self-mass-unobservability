"""Continue the saved native adiabat; keep all earlier runs immutable."""
from pathlib import Path
from types import FunctionType
import argparse
import json
import time
import numpy as np
import def_native_atmosphere_density as warm

h=warm.h
OUT=warm.parent.OUT/'cold-continuation'


def table():
    assert not OUT.exists();OUT.mkdir()
    data,tab=h.inputs();X=data['X'][0];s0=tab['reference'][0,3]
    prior=json.loads((warm.OUT/'table.json').read_text());rows=prior['rows'].copy()
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),Path(warm.__file__),warm.OUT/'table.json',warm.OUT/'result.json',warm.parent.BACKGROUND/'background-0.001.npz',warm.parent.BACKGROUND/'lapse.npz']},
        claim='Extend the same native entropy/composition toward its dilute cold limit, then reuse the same full DEF atmosphere equations. No EOS substitution or extrapolation.',
        decision='Locate the remaining material boundary and measure the enthalpy left to the native molecular ground state; finite-temperature endpoints remain finite pressure.',
        density_step=.5,maximum_new_roots=64,stop_temperature_K=.1,
        budget=dict(hard_timeout_seconds=60,expected_seconds=[2,15],workers=1,new_evolution_steps=0,automatic_expansion=False),
        entropy_energy_gate=1,stop='First native failure or entropy-gate failure, 50 seconds, 64 new roots, or T <= 0.1 K. Preserve every accepted row.'))
    eos=h.molecular.model.EOS();lt=rows[-1]['lnT'];lr0=rows[-1]['logrho'];started=time.monotonic();calls=0
    def sample(lr,t):
        nonlocal calls
        a=eos(2,float(lr),float(t),X);calls+=1
        score=abs(np.exp(t)*(a[3]-s0))/max(2.,32*np.spacing(abs(a[2]+a[1]/a[0])))
        return score,a,t
    try:
        for k in range(1,65):
            old=np.array(rows[-1]['raw']);slope=(old[1]/old[0]-old[9])/old[10]
            lr=lr0-.5*k;lt-=.5*slope;best=None
            for _ in range(18):
                trial=sample(lr,lt)
                if best is None or trial[0]<best[0]:best=trial
                if trial[0]<.25:break
                lt-=np.clip(np.exp(lt)*(trial[1][3]-s0)/trial[1][10],-.15,.15)
            score,a,lt=best
            assert score<=1,('Original entropy gate',k,score)
            gamma=a[5]+a[6]*(a[1]/a[0]-a[9])/a[10]
            assert abs(gamma/a[4]-1)<1e-8
            rows.append(dict(logrho=float(np.log(a[0])),lnT=float(lt),logP=float(np.log(a[1])),raw=a.tolist(),gamma1=float(gamma),entropy_score=float(score)))
            h.write(OUT/'progress.json',dict(classification='Counterexample candidate',rows=rows,new_roots=k,native_calls=calls,seconds=time.monotonic()-started))
            if np.exp(lt)<=.1:break
            assert time.monotonic()-started<50,'Cold continuation wall budget'
    except Exception as exc:
        h.write(OUT/'failure.json',dict(classification='Counterexample candidate',error=repr(exc),requested_logrho=lr,requested_lnT=float(lt),accepted_new_roots=len(rows)-len(prior['rows']),native_calls=calls,seconds=time.monotonic()-started))
        raise
    h.write(OUT/'table.json',dict(classification='Counterexample candidate',passed=True,rows=rows,new_roots=len(rows)-len(prior['rows']),native_calls=calls,seconds=time.monotonic()-started,reached_cold_threshold=bool(np.exp(lt)<=.1)))
    print('COLD',len(rows),'T',np.exp(lt),'P',a[1],'ground_enthalpy_erg_g',a[2]+a[11]+a[1]/a[0],flush=True)


def run():
    FunctionType(warm.run.__code__,dict(warm.run.__globals__,OUT=OUT,__file__=__file__))()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['table','run']);globals()[p.parse_args().action]()
