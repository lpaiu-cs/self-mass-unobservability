"""Counterexample candidate: include native material fronts in local time splits.

Keep the physical equation. Reuse the old128solution only after checking that
its actual stage sequence and both nominal-clock RHS owners are identical.
"""
from pathlib import Path
from types import FunctionType
import gc,hashlib,json,os,resource,sys,time
import numpy as np
import compare_full_incident_fluid_time as before

prior=before.prior;OUT=Path('native-material-front183-work');OLD=before.OUT
read,write,sha=before.read,before.write,before.sha;LD=prior.LD
CAPS=dict(prepare=15,check=60,run=450,audit=15)


def split_flags(m,steps):
    reference_h=m.t[-1]/64;h=m.t[-1]/steps
    arrivals=-np.asarray(m.redshift_driver.xc,float)/prior.C
    starts=np.r_[arrivals,arrivals+m.redshift_driver.D]
    flags=np.array([any(t<(k+1)*h and t+2*reference_h>k*h and 0<=t<m.t[-1] for t in starts) for k in range(steps)])
    return flags,np.arange(m.n),starts


def initialize():
    prior.OUT=OUT;prior.initialize();owner=prior.Model.run
    prior.Model.run=FunctionType(owner.__code__,dict(owner.__globals__,split_flags=split_flags),argdefs=owner.__defaults__)


def prepare():
    assert not OUT.exists();OUT.mkdir();assert not read(OLD/'time-result.json')['passed']
    reused={};files=[]
    for s in [0,1]:
        for folder in ['photons','material']:(OUT/f'sweep-{s}/{folder}').mkdir(parents=True)
    for src in list((OLD/'sweep-0').rglob('*.npz'))+[OLD/n for n in ['normalization.json','photon-conservation-plan.json','check-result.json']]:
        dst=OUT/src.relative_to(OLD);os.link(src,dst);reused[str(dst.relative_to(OUT))]=dict(path=str(src),sha256=sha(src));files.append(src)
    # Read-only alias for the equivalent new64split prefix. Preserve old file
    # metadata and physical bytes; it is not yet a continuation checkpoint.
    for suffix in ['.npz','.json']:
        src=OLD/f'sweep-1/photons/pilot-128{suffix}';dst=OUT/f'sweep-1/photons/pilot-64{suffix}'
        os.link(src,dst);reused[str(dst.relative_to(OUT))]=dict(path=str(src),sha256=sha(src));files.append(src)
    write(OUT/'reuse.json',reused)
    row=read(OLD/'pilot-result.json');receipt=read(OLD/'retry-receipt.json')
    forecast=2*row['stepping_seconds']+row['operator_point_seconds']+receipt['seconds']-row['seconds']+20
    assert forecast<CAPS['run']
    files += [Path(v.__file__) for v in list(sys.modules.values()) if getattr(v,'__file__',None)
        and Path(v.__file__).parent.name=='verification' and Path(v.__file__).suffix=='.py']
    files += [Path(__file__),OLD/'time-result.json',OLD/'time-localization.json',OLD/'front-cascade.json',prior.drive.FIELDS/'born-g8.npz']
    write(OUT/'plan.json',dict(classification='Conjectural',checkpoint='ea88eb835',
        claim='Repair the omitted native-material front in time splitting, then test the actual same full-input photon/material solution against original2percent gates.',
        evidence='182Btime5.89185percent fails; its floor contribution is1.14e-9 of the difference. Native Etilde transport also dominates its error. The photon-slow-only rule misses14us and4us fronts; a no-fit exact-pulse cascade predicts6.77percent and5.89percent paired error, reduced below1percent by one split. This is mechanistic evidence, not a full-system guarantee.',
        correction='Use all existing native cell primary leading/trailing arrival times with the same window2*T/64. Bisect each overlapping macro interval exactly once. No selection or fit by output error, no extra cell, EOS, waveform, background or stage method. Stiff photon cells do not justify dropping independent native material fronts.',
        budget_reassessment='Original182failure stays failed. Before dispatch, verify the modified64prefix is the same actual sequence as saved182128 and the two nominal-clock RHS owners are identical. Reuse that path; run only modified128 over T/32, eight actual substeps. Effective step is T/256, disclosed explicitly; this is a newly justified material-front repair, not automatic grid escalation.',
        gates=dict(time=.02,stage=1e-12,physical_stage=1e-13,constitutive=.002,conservation=1e-8,port=1e-12,max_Newton_solves=3),
        reuse='All saved physical arrays remain byte-identical. The alias pilot-64 retains original128metadata and must not be used as an ordinary restart without an explicit metadata/history adapter.',
        forecast=dict(seconds=forecast,assumed_range_seconds=[350,450],hard_cap_seconds=450,
            basis='Twice measured four-substep196.934s stepping plus point/setup/I/O and20s reserve. New fine-step branches remain an assumption.'),
        budgets=CAPS,CPU_threads=1,virtual_GiB=4,
        stop='Any equivalence, original physical/time gate, iteration or cap fails this trial. No second bisection, finer base clock, new method, automatic full horizon or GR run.',
        scope='Retained directional full declared incident response. Previous native gross-recovery derivative failure and physical EOS, nonlinear, boundary, static and observational gaps remain. Final charge is unadjudicated.',
        final_charge_conclusion='unadjudicated',full_goal_complete=False,bindings={str(p):sha(p) for p in dict.fromkeys(files)}))


def fingerprint(value):
    if hasattr(value,'tocsr'):
        a=value.tocsr();return [a.shape,*[fingerprint(v) for v in [a.data,a.indices,a.indptr]]]
    a=np.ascontiguousarray(value)
    if a.dtype==np.dtype(LD):
        # x86 extended precision has six padding bytes; hash its80 value bits.
        assert sys.byteorder=='little' and a.dtype.itemsize==16 and np.finfo(LD).nmant==63
        data=a.view(np.uint8).reshape(-1,16)[:,:10].tobytes()
    else:data=a.tobytes()
    return [str(a.dtype),list(a.shape),hashlib.sha256(data).hexdigest()]


def check():
    initialize();p=dict(np.load(OLD/'sweep-1/photons/pilot-128.npz'));owners=[];sequences=[]
    for n in [64,128]:
        m=prior.Model(n);flags,cells,fronts=split_flags(m,n);h=m.t[-1]/n;edges=[0.];tt=[];ww=[]
        for k in range(n//32):
            parts=2 if flags[k] else 1
            for sub in range(parts):
                t=k*h+sub*h/parts;hh=h/parts
                tt.extend(t+hh*prior.joint.C);ww.extend(hh*prior.joint.B);edges.append(t+hh)
        sequences.append(dict(clock=n,edges=edges,stage_times=list(tt),weights=list(ww),split_macros=flags[:n//32].tolist()))
        if n==64:
            for a,key in [(edges,'actual_step_edges'),(tt,'joint_stage_times'),(ww,'joint_stage_weights')]:assert np.array_equal(a,p[key]),key
        fixed={k:fingerprint(getattr(m,k)) for k in ['t','I','A','Nweight','Eweight','scale','eu','nu','bu','su','kappa']}
        stages=[]
        for t,z in zip(p['joint_stage_times'],p['joint_stage_conserved_scaled']):
            g=np.column_stack([(z[2]-m.kappa*z[0])/m.eu,z[3]/m.nu,z[0]/m.bu,z[1]/m.su])
            c=m.local(float(t));s=m.source(float(t));native=m.native(float(t),g)
            stages.append(dict(local={k:fingerprint(v) for k,v in c.items()},source=[fingerprint(v) for v in s],native=fingerprint(native),
                pressure=fingerprint(m.pressure(float(t),g)),active=fingerprint(m.material.active(float(t)))))
        owners.append(dict(fixed=fixed,stages=stages));del m;gc.collect()
    assert owners[0]==owners[1],'Nominal clock changes physical RHS'
    write(OUT/'equivalence.json',dict(classification='Counterexample candidate',passed=True,
        owners_identical=True,original128_is_new64_step_sequence=True,all_stored_physical_bytes_reused=True,
        prefix_only=True,ordinary_restart_ready=False,sequences=sequences,owner_fingerprints=owners[0]))
    # Exact Radau weights on the changed intervals still integrate degree0..2.
    from fractions import Fraction as F
    for n in [4,8]:
        for power in range(3):assert sum((F(3,4)*((F(k)+F(1,3))/n)**power+F(1,4)*(F(k+1)/n)**power)/n for k in range(n))==F(1,power+1)
    write(OUT/'symbolic.json',dict(classification='Proven',passed=True,scope='Exact rational Radau quadrature moments0..2 on the candidate physical step intervals. No full-system error bound.'))


def run():
    assert read(OUT/'equivalence.json')['passed'];initialize()
    pilot,source=prior.clone(prior.pilot,[('initialize();m=Model(64)','m=Model(128)'),
        ("m.run(64,'pilot-64',2)","m.run(128,'pilot-128',4)"),("paths(1)[0]/'pilot-64.npz'","paths(1)[0]/'pilot-128.npz'")])
    (OUT/'expanded-dispatch.py').write_text(source);pilot()
    compare=FunctionType(before.compare.__code__,dict(before.compare.__globals__,OUT=OUT));compare()


def audit():
    fn=FunctionType(before.audit.__code__,dict(before.audit.__globals__,OUT=OUT));fn()


if __name__=='__main__':
    action=sys.argv[1];assert action in CAPS;receipt=OUT/f'{action}-receipt.json';assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));prior.joint.previous.original.inf.incident.native.deadline(CAPS[action])
    start=time.monotonic();cpu=time.process_time();error=None
    try:
        if action!='prepare':
            for p,h in read(OUT/'plan.json')['bindings'].items():assert sha(p)==h,p
        globals()[action]()
    except BaseException as exc:error=repr(exc);raise
    finally:
        if OUT.exists():write(receipt,dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=sha(__file__)))
