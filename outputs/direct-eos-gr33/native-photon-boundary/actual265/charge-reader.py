"""Use the unchanged full-period charge reader on the repaired actual return."""
from pathlib import Path
import resource,sys,time
import read_full_return_charge as reader

ROOT=Path('native-photon-boundary-ld2-charge265-work');OUT=ROOT/'full'
ACTUAL=Path('native-photon-boundary-ld2-265-work')
read,write,sha,bind=reader.read,reader.write,reader.sha,reader.bind
endpoint,dense,charge=reader.endpoint,reader.dense,reader.charge
reader.ACTUAL=ACTUAL;reader.OUT=OUT
endpoint.INPUT=ACTUAL;endpoint.OUT=OUT/'endpoint'
dense.INPUT=endpoint.OUT;dense.OUT=OUT/'dense'
charge.INPUT=dense.OUT;charge.OUT=OUT/'charge'
charge.Response.setup=bind(endpoint.returned.Response.setup,INPUT=dense.OUT,bridge=charge.bridge)
reader.endpoint_prepare=bind(reader.endpoint_prepare,INPUT=ACTUAL,OUT=endpoint.OUT)
reader.charge_prepare=bind(reader.charge_prepare,INPUT=dense.OUT,OUT=charge.OUT,ACTUAL=ACTUAL)


if __name__=='__main__':
    action=sys.argv[1];assert action in reader.CAPS and action!='regression'
    assert read(ACTUAL/'controller-status.json')['state']=='completed'
    result=read(ACTUAL/'result.json');assert result['passed'] and result['scientific_gates_changed'] is False
    original=Path(reader.__file__)
    assert read(reader.ROOT/'regression-receipt.json')['source_sha256']==sha(original)
    assert read(reader.ROOT/'regression.json')['passed']
    OUT.mkdir(parents=True,exist_ok=True);part,verb=action.split('_',1)
    folder=dict(endpoint=endpoint.OUT,dense=dense.OUT,charge=charge.OUT)[part]
    receipt=folder/(verb+'-receipt.json');assert not receipt.exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    reader.base.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(reader.CAPS[action])
    start=time.monotonic();error=None
    try:
        reader.execute(action)
        if verb=='prepare':
            p=read(folder/'plan.json');p.update(actual_solution_directory=str(ACTUAL),
                input_adapter='Only input/output directories change. The original251short end-to-end229array/charge proof and every numerical gate remain applicable to the same245/246/247reader; full257input has its own completed physical audit.',
                input_adapter_sha256=sha(__file__),claim='Read the actual265full-period complete-photon-boundary return with its own captured matter/photon/port and applied metric histories.')
            p['bindings'].update({str(Path(__file__)):sha(__file__),str(original):sha(original),str(ACTUAL/'result.json'):sha(ACTUAL/'result.json')})
            write(folder/'plan.json',p)
    except BaseException as exc:error=repr(exc);raise
    finally:
        if folder.exists():write(receipt,dict(seconds=time.monotonic()-start,error=error,source_sha256=sha(__file__),original_reader_sha256=sha(original),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
