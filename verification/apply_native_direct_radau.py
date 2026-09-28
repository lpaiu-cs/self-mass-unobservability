"""Counterexample candidate: apply the physical source at actual Radau stages.

No affine lift or fitted forcing: preserve the direct stiff stage balance.
"""
from pathlib import Path
import json,resource,shutil,sys,time
import numpy as np
import solve_native_radau_transfer as radau

OUT=Path('native-direct-radau171-work');OLD=Path('native-known-forcing170-work')
original=radau.prior;read,write,sha=radau.read,radau.write,radau.sha
CAPS=dict(prepare=30,check=60,pilot=180,finish=30);TOTAL=300


def initialize():
    assert not (OUT/'lift-plan.json').exists() and not (OUT/'primitive-plan.json').exists()
    radau.OUT=OUT;radau.initialize()


def prepare():
    assert not OUT.exists() and not read(OLD/'result.json')['passed'];OUT.mkdir()
    files=[Path(__file__),Path(radau.__file__),Path(original.__file__),OLD/'result.json',OLD/'stiff-limit.json']
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    for p in (OLD/'sweep-0').rglob('*.npz'):
        shutil.copyfile(p,OUT/p.relative_to(OLD));files.append(p)
    for name in ['normalization.json','photon-conservation-plan.json']:
        shutil.copyfile(OLD/name,OUT/name);files.append(OLD/name)
    write(OUT/'plan.json',dict(classification='Counterexample candidate',checkpoint='d40b5c21fa57d904f0ada4194862f1842fc8f1fa',
        claim='Apply the unchanged face-redshift correction source directly to the actual coupled Radau stages, removing the affine/source-fitting stiff-limit defect.',
        rationale='170proved on the scalar stiff test that exact source primitives alter the leading stage balance, and full LH fitting can leave a spurious nonzero stiff limit. This does not identify every transport error. Test the simplest existing stage-consistent source path before any clock/order increase.',
        method='Reuse168simultaneous two-stage Radau without167lift-plan/primitive-plan. Existing zero-lift owner supplies H=0. Use actual source S(tj), collision, gas, ghost/work and actual angular samples with original3/4,1/4 weights. No moment-fitted rate or separate angular correction.',
        decision='Unchanged4/8-step pair must pass all six2percent time gates and physical checks before any full coupled production/charge readout. The final-charge conclusion remains unadjudicated until that same corrected coupled solution is completed.',
        inputs='Same531cells8angles152frequencies,17background knots,64/128 clocks and0.214652ms horizon; zero initial correction, previous material motion and extra metric.',
        gates=dict(time=.02,energy_H=1e-8,stage=1e-12,physical_stage_moment=1e-13,port=1e-12,source=1e-12),
        budget=dict(actions=CAPS,total_action_seconds=TOTAL,CPU_threads=1,virtual_GiB=3),
        forecast='168with exact H cost85.80s;170with moment fits142.44s. Removing all known-source moment integrations should cost60-110s subject to contention;180s hard cap. Krylov changes remain unmeasured. One short pair only.',
        stop='Any original gate or cap: preserve failure and stop. No automatic refinement, new scheme, period or production.',
        bindings={str(p):sha(p) for p in dict.fromkeys(files)}))
    shutil.copyfile(OLD/'stiff-limit.json',OUT/'symbolic.json')


def check():
    initialize();m=original.c.Response(128);rows=[]
    assert not np.any(m.motion) and not np.any(m.energy_offset)
    for t in [m.t[-1]/32,m.t[-1]*.31]:
        assert all(not np.any(v) for v in m.lift(t)[:4]);a={}
        for factor in [0.,1.,-1.,2.]:
            m.drive_scale=factor;s,l,e=m.source(t);expected,ports,err=original.correction(m,t)
            assert np.array_equal(s,expected) and np.array_equal(l,ports) and max(e,err)<1e-12
            a[factor]=(s,l)
            c=m.local(t);assert all(not np.any(c[k]) for k in ['q','qb','qe','mechanical'])
            number,energy=m.moments(s);norm=[np.sum(abs(s)*m.weights,dtype=np.longdouble),np.sum(abs(s)*m.weights*m.E,dtype=np.longdouble)]
            balance=float(np.max(abs(np.array([number+l[0],energy+l[1]-l[2]]))/np.maximum(norm,1e-290)))
            assert balance<1e-12,(t,factor,balance)
        assert all(not np.any(v) for v in a[0.])
        odd=max(float(np.max(abs(a[k][j]-k*a[1.][j]))/max(np.max(abs(a[1.][j])),1e-290)) for k in [-1.,2.] for j in [0,1])
        assert odd<1e-12;rows.append(dict(time=float(t),source_is_exact_original_correction=True,zero_lift=True,zero_exact=True,odd_scaling_relative=odd))
    h=m.t[-1]/64;rates=(-m.A.diagonal()).reshape(m.n,m.q)
    result=dict(classification='Counterexample candidate',passed=True,rows=rows,
        stream_diagonal_h_range=[float(np.min(h*rates)),float(np.max(h*rates))],
        near_pressure_error_cells={str(k):[float(np.min(h*rates[k])),float(np.max(h*rates[k]))] for k in [15,16]},
        scope='Source/ledger/zero-lift checks and actual stream diagonal scale; diagonals are not full coupled spectral eigenvalues or a unique-cause proof.')
    write(OUT/'source-check.json',result);print(json.dumps(result),flush=True)


def pilot():
    assert read(OUT/'source-check.json')['passed'];radau.OUT=OUT;radau.pilot()


def finish():
    radau.OUT=OUT;radau.finish();r=read(OUT/'result.json');r.pop('export_repair',None)
    r['final_charge_conclusion']='unadjudicated on the corrected coupled solution'
    r['source_representation']='Actual direct physical source at Radau stages; no affine primitive or moment-fitted source.'
    write(OUT/'result.json',r)


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
            assert sum(read(p)['seconds'] for p in OUT.glob('*-receipt.json'))+CAPS[action]<=TOTAL
        globals()[action]()
    except Exception as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(action=action,seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,error=error,source_sha256=sha(__file__)))
