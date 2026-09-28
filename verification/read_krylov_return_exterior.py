"""Apply the unchanged frozen-exterior reader to the repaired actual return."""
from pathlib import Path
import inspect,resource,sys,time
import read_returned_exterior as prior
import verify_returned_exterior as verify

ROOT=Path('native-krylov-exterior254-work')


if __name__=='__main__':
    assert sys.argv[1]=='full'
    original=prior.ROOT
    assert prior.read(original/'common/audit.json')['passed']
    assert prior.read(original/'check-receipt.json')['source_sha256']==prior.sha(prior.__file__)
    assert prior.self_check()['passed']
    source=inspect.getsource(prior.run)
    for a,b in [('native-full-return249-work/full','native-returned-krylov254-work'),
                ('native-full-charge251-work/full/charge','native-krylov-charge254-work/full/charge')]:
        assert source.count(a)==1;source=source.replace(a,b)
    prior.ROOT=ROOT;namespace=dict(prior.run.__globals__)
    exec(compile(source,__file__,'exec'),namespace)
    ROOT.mkdir(exist_ok=True);assert not (ROOT/'full-receipt.json').exists()
    resource.setrlimit(resource.RLIMIT_AS,(16*1024**3,)*2)
    prior.common.prior.prior.joint.previous.original.inf.incident.native.deadline(1800)
    start=time.monotonic();error=None
    try:
        namespace['run']('full')
        p=prior.read(ROOT/'full/plan.json')
        p['bindings'].update({str(Path(__file__)):prior.sha(__file__)})
        p['input_adapter']='Only254actual and charge paths replace249/251. Identical252reader, initial exterior, normalization, component controls and scope limits.'
        prior.write(ROOT/'full/plan.json',p)
    except BaseException as exc:error=repr(exc);raise
    finally:
        prior.write(ROOT/'full-receipt.json',dict(seconds=time.monotonic()-start,error=error,source_sha256=prior.sha(prior.__file__),
            input_adapter_sha256=prior.sha(__file__),peak_RSS_bytes=1024*resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
    verify.audit('full')
