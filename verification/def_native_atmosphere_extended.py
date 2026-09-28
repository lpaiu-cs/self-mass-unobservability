"""Same molecular EOS, saved extended arithmetic in density coordinates."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
import mpmath as mp
from types import FunctionType
import def_native_atmosphere_density as warm
import gr_molecular_precision as precision

h=warm.h
OUT=warm.parent.OUT/'extended-cold'


def table():
    assert not OUT.exists();OUT.mkdir()
    mp.mp.dps=65;data,tab=h.inputs();X=data['X'][0];s0=mp.mpf(float(tab['reference'][0,3]))
    old=json.loads((warm.OUT/'table.json').read_text());rows=old['rows'].copy()
    h.write(OUT/'plan.json',dict(classification='Counterexample candidate',
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in [Path(__file__),Path(warm.__file__),warm.OUT/'table.json',warm.OUT/'result.json',precision.OUT/'manifest.json',warm.parent.OUT/'cold-continuation/precision-probe.json']},
        correction='Reuse the saved matching 33-digit native EOS. Only arithmetic changes; no coefficients, energy gate, molecular or plasma model changes.',
        gates=dict(entropy_energy_score=1,overlap_logrho=1e-12,overlap_logP=1e-10,overlap_entropy_relative=1e-10),
        budget=dict(hard_timeout_seconds=90,maximum_new_roots=24,new_evolution_steps=0,automatic_expansion=False),
        target_temperature_K=100,density_step=.5,
        stop='First native failure, entropy gate failure, 80 seconds, 24 new roots or T <= 100 K. Preserve accepted rows and all original binary64 failures.'))
    eos=precision.EOS();ym=(X/precision.g.c.A)@eos.mapping;cx=float(ym@eos.weights);eps=ym/cx
    seed=np.zeros(24);seed[2]=1;seed/=seed@eos.weights
    calls=0;start=time.monotonic()
    def sample(lr,lt):
        nonlocal calls
        eos.raw(0,mp.mpf(-20),mp.mpf(float(np.log(1e6))),seed,1.)
        a=eos.raw(2,lr+mp.log(cx),lt,eps,cx);calls+=1
        budget=max(2.,32*np.spacing(abs(float(a[2]+a[1]/a[0]))))
        return abs(mp.exp(lt)*(a[3]-s0))/budget,a
    lr0=mp.mpf(rows[-1]['logrho']);lt=mp.mpf(rows[-1]['lnT']);_,control=sample(lr0,lt)
    prev=np.array(rows[-1]['raw']);overlap=dict(logrho=float(abs(mp.log(control[0])-lr0)),logP=float(abs(mp.log(control[1])-rows[-1]['logP'])),entropy_relative=float(abs(control[3]/prev[3]-1)))
    assert overlap['logrho']<1e-12 and overlap['logP']<1e-10 and overlap['entropy_relative']<1e-10,overlap
    h.write(OUT/'overlap.json',dict(classification='Counterexample candidate',passed=True,**overlap))
    a=control
    try:
        for k in range(1,25):
            lr=lr0-mp.mpf('.5')*k;lt-=mp.mpf('.5')*(a[1]/a[0]-a[9])/a[10]
            for _ in range(18):
                score,a=sample(lr,lt)
                if score<mp.mpf('.25'):break
                lt-=max(mp.mpf('-.15'),min(mp.mpf('.15'),mp.exp(lt)*(a[3]-s0)/a[10]))
            assert score<=1,('Original entropy gate',k,str(score))
            gamma=a[5]+a[6]*(a[1]/a[0]-a[9])/a[10]
            assert abs(gamma/a[4]-1)<mp.mpf('1e-8')
            rows.append(dict(logrho=float(mp.log(a[0])),lnT=float(lt),logP=float(mp.log(a[1])),raw=[float(v) for v in a],raw_extended=[mp.nstr(v,40) for v in a],gamma1=float(gamma),entropy_score=float(score)))
            h.write(OUT/'progress.json',dict(classification='Counterexample candidate',rows=rows,new_roots=k,native_calls=calls,seconds=time.monotonic()-start))
            if mp.exp(lt)<=100:break
            assert time.monotonic()-start<80,'Extended continuation wall budget'
    except Exception as exc:
        h.write(OUT/'failure.json',dict(classification='Counterexample candidate',error=repr(exc),accepted_new_roots=len(rows)-len(old['rows']),native_calls=calls,seconds=time.monotonic()-start))
        raise
    h.write(OUT/'table.json',dict(classification='Counterexample candidate',passed=True,rows=rows,new_roots=len(rows)-len(old['rows']),native_calls=calls,seconds=time.monotonic()-start,reached_cold_threshold=bool(mp.exp(lt)<=100)))
    print('EXTENDED',len(rows),'T',float(mp.exp(lt)),'P',float(a[1]),'seconds',time.monotonic()-start,flush=True)


def run():
    FunctionType(warm.run.__code__,dict(warm.run.__globals__,OUT=OUT,__file__=__file__))()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['table','run']);globals()[p.parse_args().action]()
