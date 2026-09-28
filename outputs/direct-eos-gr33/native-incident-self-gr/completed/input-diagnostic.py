from pathlib import Path
import resource,time,json
import numpy as np
import solve_native_incident_self_gr as run
out=run.OUT;start=time.monotonic();cpu=time.process_time();error=None
resource.setrlimit(resource.RLIMIT_AS,(3*1024**3,3*1024**3));run.coupled.base.drive.native.deadline(45)
try:
    assert sum(run.read(p)['seconds'] for p in out.rglob('*-receipt.json'))+45<4600
    run.initialize(1);m=run.coupled.Material(128,64);m.driver=run.SavedMetric(128)
    source=run.coupled.base.old.aligned.run_source
    ns=dict(m.run_owner.__func__.__globals__);exec(compile(source,__file__,'exec'),ns)
    m.run_owner=ns['run'].__get__(m,type(m))
    row=m.run(64,'common-metric-ssp2-64',4)
    actual=np.load(run.paths(1)[1]/'common-metric-ssp2-64.npz')['delta_scaled']
    fine=np.load(out/'sweep-1/material-analytic/steps-128-reference-128.npz')['history_scaled'][8]
    relative=(np.sum(abs(actual-fine),axis=1)/np.maximum(np.sum(abs(fine),axis=1),1.)).astype(float).tolist()
    result=dict(classification='Counterexample candidate',row=row,first_interval_common_metric_time_difference=relative,
        original_metric_two_path_first_interval_difference=[.07579788928909155,.2729043696230876,.46008682247390575,.0032998064409323785],
        changed='Only prescribe the SAME fine generated GR field to both clocks. Original SSP2 and separate64/128 photon collision transfers retained. Fine completed trajectory reused; one4-step coarse prefix.',
        passed=row['passed'] and max(relative)<.02,production_authorized=False)
    run.write(out/'input-diagnostic.json',result);print(json.dumps(result),flush=True)
except Exception as exc:error=repr(exc);raise
finally:run.write(out/'input-diagnostic-receipt.json',dict(seconds=time.monotonic()-start,CPU_seconds=time.process_time()-cpu,
    peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,error=error,source_sha256=run.sha(__file__)))
