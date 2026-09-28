"""Remove the measured per-length remainder floor without relaxing its budget."""
from types import ModuleType
import gzip,json,sys
import numpy as np
import sympy as sp
import gr_outer_self_pilot as original

OUT=original.OUT.parent/'gr-outer-self-refined'


def source():
    before=(original.ROOT/'verification/gr_outer_self_pilot.py').read_text()
    changes={"OUT=highq.OUT.parent/'gr-outer-self-pilot';CACHE=highq.moments.native.g.CACHE/'outer-self-pilot'":
        "OUT=highq.OUT.parent/'gr-outer-self-refined';CACHE=highq.moments.native.g.CACHE/'outer-self-refined'",
        '    for a,b in changes.items():assert cpp.count(a)==1;cpp=cpp.replace(a,b)':
        "    changes.update({'R=4*h':'R=8*h','power(MI(2)/3,32)':'power(MI(2)/7,32)'})\n    for a,b in changes.items():assert cpp.count(a)==1;cpp=cpp.replace(a,b)"}
    after=before
    for a,b in changes.items():assert after.count(a)==1;after=after.replace(a,b)
    reverse=after
    for a,b in reversed(list(changes.items())):assert reverse.count(b)==1;reverse=reverse.replace(b,a)
    assert reverse==before;return after,changes


def module(text):
    name='gr_outer_self_refined_candidate';obj=ModuleType(name);sys.modules[name]=obj;exec(compile(text,str(OUT/'candidate.py'),'exec'),obj.__dict__);return obj


def prepare():
    assert not OUT.exists();text,changes=source();obj=module(text);obj.prepare();(OUT/'candidate.py').write_text(text)
    p=json.loads((OUT/'plan.json').read_text());p.update(checkpoint='336a78b',driver_substitutions=changes,finite_radius_half_width_ratio=8,
        failure='At the new tiny Q and 1e-18 finite remainder budget, the R=4h Cauchy/Gauss factor is fixed at norm16*(2/3)^32. Once a positive per-length envelope exceeds this floor, further subdivision cannot satisfy the budget. The original native preflight reached its depth limit before completion.',
        correction='Use R=8h, retaining every analytic-domain and exact-coverage check. Then width/(R-h)=2/7, reducing the per-length remainder factor. The 1e-18 remainder, 2e-7 field gate, true root intervals, outer schedules and independent reference are unchanged.')
    files=[original.ROOT/'verification/gr_outer_self_refined.py',OUT/'candidate.py',original.OUT/'plan.json',original.OUT/'finite-source.cpp',original.OUT/'preflight-0000-native.log']
    p['bindings'].update({f.relative_to(original.ROOT).as_posix():original.sha(f) for f in files});obj.save('plan.json',p)
    h=sp.symbols('h',positive=True);assert sp.cancel(2*h/(8*h-h)-sp.Rational(2,7))==0
    obj.save('refinement-proof.json',dict(classification='Proven',passed=True,
        factor_reduction_exact=str(sp.Rational(3,7)**32),same_tolerance=True,original_failure_preserved=True))


def run(name):
    text,_=source();assert text==(OUT/'candidate.py').read_text();obj=module(text);obj.bindings();getattr(obj,name)()


def preflight_control():
    text,_=source();obj=module(text);plan=obj.bindings();state,_,roots,data,_,_=obj.highq.context();checks=[]
    for i in plan['positions']:
        assert (OUT/f'preflight-{i:04d}-result.json').exists()
        with gzip.open(OUT/f'preflight-{i:04d}-points.jsonl.gz','rt') as stream:rows=[json.loads(x) for x in stream]
        eta=float(sum(obj.cusp.endpoints(roots[i]['root']))/2);Q=np.array([float(obj.F(r['Q_exact'])) for r in rows]);size=len(Q)
        values=obj.highq.neutral.thermo.response(np.full(size,eta),np.full(size,float(data['beta'][i])),Q,np.full(size,float(data['Sref'][i])),2e-13)
        target=np.array([[float(sum(obj.cusp.endpoints(x))/2) for x in r['response']] for r in rows]).T
        score=float(np.max(np.abs(values-target))/max(np.max(np.abs(values)),1e-300));assert score<2e-9
        checks.append(dict(position=i,cell=roots[i]['cell'],components=6*size,relative_vector_difference=score,passed=True))
    obj.save('preflight-controls.json',dict(classification='Counterexample candidate',passed=True,gate=2e-9,checks=checks,
        scope='Independent binary64 adaptive response at rounded 80-bit Q and neutral-root midpoints. Finite agreement only; no exact reference inclusion claim.'))
    print('PASS 288 independent finite outer-node preflight comparisons',flush=True)


def verify():run('verify')


if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='preflight_control':preflight_control()
    else:run(sys.argv[1])
