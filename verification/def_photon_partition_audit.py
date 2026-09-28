"""Check native group additivity and the conditional loss-moment theorem.

No opacity rescaling and no new stellar evolution. Failed source partitions
remain failed even if the abstract moment inequality is correct.
"""
from pathlib import Path
from types import FunctionType
import argparse
import json
import signal
import time
import numpy as np
import sympy as sp
import def_photon_display_fixed as fixed

native=fixed.native; ex=fixed.ex; h=fixed.h
loss=fixed.previous.previous
OUT=fixed.OUT.parent/'def-photon-partition-audit'


def proof():
    x,f,c,r,a=sp.symbols('x f c r a',positive=True)
    gap=(x-f)*(x-c)**2/((x+r)*(f+r)*(c+r)**2)
    polynomial=sp.cancel(gap+1/(x+r))
    assert sp.denom(polynomial).has(x) is False
    assert sp.degree(polynomial,x)==2
    assert sp.factor(polynomial-1/(x+r)-gap)==0
    assert sp.factor(polynomial.subs(x,f)-1/(f+r))==0
    assert sp.factor(polynomial.subs(x,c)-1/(c+r))==0
    assert sp.factor(sp.diff(polynomial,x).subs(x,c)+1/(c+r)**2)==0
    assert loss.symbolic()['passed']
    # Independent finite measures test the actual complex resolvent, including
    # near-constant and extremely separated opacities, without table helpers.
    rng=np.random.default_rng(5901); worst=0.; tested=0
    for n in [1,2,7,31]:
        for _ in range(25):
            k=np.exp(rng.uniform(-12,12,n)); mu=rng.dirichlet(np.ones(n))
            floor=float(k.min()/2); H=float(mu@(1/k)); M=float(mu@k)
            aa=1/H; cc=(M-floor)/(1-floor*H)
            omega=np.r_[0.,np.geomspace(k.min()*1e-4,k.max()*1e4,81)]
            direct=np.sum(mu[:,None]/(k[:,None]+1j*omega),axis=0)
            exact_error=abs(1/(aa+1j*omega)-direct)
            envelope=2*np.sqrt(2)*(M-aa)*omega/((floor+omega)*(aa+omega)*(cc+omega))
            scale=128*np.finfo(float).eps*H
            assert np.all(exact_error<=np.maximum(envelope,0)+scale)
            if M-aa>1e-10*M:
                worst=max(worst,float(np.max(exact_error/(envelope+scale))))
            for p in [1,2]:
                maximum=float(loss.maxima(floor,np.array([aa]),np.array([cc]),p)[0])
                sampled=omega**p/((floor+omega)*(aa+omega)*(cc+omega))
                assert sampled.max()<=maximum*(1+1e-12)
            tested+=len(omega)
    return dict(classification='Proven',passed=True,
        majorant='P2(x)-1/(x+r)=(x-f)(x-c)^2/[(x+r)(f+r)(c+r)^2]>=0 for x>=f; P2 is quadratic, touches f and is tangent at c. Its expectation is fixed by the first two reweighted moments.',
        scope='The exact conditional moment inequality is proved algebraically. Floating-point controls are separate finite checks, not a certified interval implementation or proof that incompatible native partitions share an opacity measure.',
        numerical_checks=tested,maximum_error_to_bound_ratio=worst)


def prepare():
    assert not OUT.exists();OUT.mkdir()
    old=next(r for r in json.loads((fixed.OUT/'requests.json').read_text()) if r['name']=='phLine4T0')
    rows=[]
    for count in [2,4]:
        name=f'phPart{count}'
        rows.append(dict(old,name=name,boundaries=count,fields=dict(old['fields'],mixname=name,ngpengs=str(count))))
    ex.write(OUT/'requests.json',rows)
    (OUT/'lanl-tops-form.html').write_bytes((fixed.OUT/'lanl-tops-form.html').read_bytes())
    paths=[Path(__file__),Path(loss.__file__),OUT/'requests.json',fixed.OUT/'phLine4T0-table.txt',fixed.OUT/'phLine4T0-generated-results-request.json',fixed.OUT/'result.json']
    ex.write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='0f6b09ed',
        claim='Determine whether identical native opacity inputs produce additive harmonic and arithmetic integrals under nested partitions of the same line window.',
        decision='If the 1/3/999 group integrals disagree beyond the original 1e-3 mean gate, reject fine native group moments as an identified common loss measure. Do not normalize them or increase resolution automatically.',
        boundaries=[2,4,1000],reuse_999_groups=True,gates=dict(mean_relative=.001,input_printing_score=1),
        budget=dict(queries=2,forecast_seconds=5,hard_seconds=45,per_query_seconds=20,CPU_workers=1,new_stellar_steps=0,new_native_EOS_calls=0),
        bindings={p.relative_to(h.ROOT).as_posix():h.digest(p) for p in paths}))
    ex.write(OUT/'symbolic.json',proof())


def fetch():
    plan=json.loads((OUT/'plan.json').read_text())
    for p,sha in plan['bindings'].items():assert h.digest(h.ROOT/p)==sha,p
    start=time.monotonic();results=[]
    for row in json.loads((OUT/'requests.json').read_text()):
        assert not list(OUT.glob(row['name']+'-*'))
        def display(form):
            v=native.retrieval.fields(form);ex.write(OUT/(row['name']+'-generated-results-request.json'),v)
            return dict(v,egplow=f"{float(v['egplow'])*(1-1e-4):.17g}",egphigh=f"{float(v['egphigh'])*(1+1e-4):.17g}")
        env=dict(vars(native.retrieval),OUT=OUT,fields=display)
        env['save']=lambda name,value:ex.write(OUT/name,value)
        fn=FunctionType(native.retrieval.fetch.__code__,env,argdefs=native.retrieval.fetch.__defaults__)
        began=time.monotonic();signal.alarm(20)
        try:result=fn(row)
        except Exception as exc:
            ex.write(OUT/'failure.json',dict(error=repr(exc),name=row['name']));raise
        finally:signal.alarm(0)
        result['seconds']=time.monotonic()-began;results.append(result)
        ex.write(OUT/'progress.json',results);print('QUERY',row['name'],result,flush=True)
        assert time.monotonic()-start<45
    ex.write(OUT/'retrieval.json',dict(seconds=time.monotonic()-start,queries=len(results)))


def analyze():
    assert not (OUT/'result.json').exists()
    rows=json.loads((OUT/'requests.json').read_text())
    original=next(r for r in json.loads((fixed.OUT/'requests.json').read_text()) if r['name']=='phLine4T0')
    header=native.reader.namespace()['header'];checks=[]
    for folder,row in [(OUT,r) for r in rows]+[(fixed.OUT,original)]:
        path=folder/(row['name']+'-table.txt')
        common=header(path,dict(row,fields=dict(row['fields'],datype='gray')))
        assert common['input_passed']
        lines=[s.strip() for s in path.read_text().splitlines() if s.strip()]
        n=int(next(s for s in lines if s.startswith('Photon grid')).split()[-2]);assert n==row['boundaries']-1
        j=next(i for i,s in enumerate(lines) if s.startswith('Energy') and 'density =' in s)
        tokens=[s.split() for s in lines[j+1:j+1+n]];g=np.array(tokens,float)
        assert np.isfinite(g).all() and np.all(g>0) and np.all(g[:,1:]!=1e10)
        e=np.geomspace(float(row['fields']['egplow']),float(row['fields']['egphigh']),n+1)
        assert max(native.audit.score(native.audit.F(float(v)),t[0]) for v,t in zip(e[:-1],tokens))<=1
        w=loss.weights(e,.0015)
        means=[float(w[:,1].sum()/sum(w[:,1]/g[:,1])),float(w[:,0]@g[:,2]/w[:,0].sum())]
        checks.append(dict(groups=n,Rosseland_Planck=means,gray=common['means']))
    assert all(c['gray']==checks[0]['gray'] for c in checks)
    values=np.array([c['Rosseland_Planck'] for c in checks]);diff=abs(values/values[0]-1)
    result=dict(classification='Counterexample candidate',passed=bool(diff.max()<.001),checks=checks,
        relative_to_one_group=diff.tolist(),gray_headers_identical=True,
        conclusion='Nested native group means are additive within the gate.' if diff.max()<.001 else 'Native group moments fail partition additivity at fixed physical inputs. These partitions cannot all be exact moments of one common opacity measure within the gate. The internal provider integration algorithm is not established by this test.',
        local_loss_response_accepted=False,whole_star_heat_closed=False,full_dynamic_charge_solved=False)
    ex.write(OUT/'result.json',result);print('RESULT',result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['prepare','fetch','analyze']);globals()[parser.parse_args().action]()
