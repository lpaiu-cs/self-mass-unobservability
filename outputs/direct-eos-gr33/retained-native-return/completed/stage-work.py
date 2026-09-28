"""Move only this task's intermediate I/O; preserve every numerical array."""
from pathlib import Path
import hashlib,json,signal,time
src=Path('outputs/direct-eos-gr33/retained-native-return')
dst=Path('/home/lpaiu/work/native-retained-tail-runtime/retained-native-return150-work')
read=lambda p:json.loads(p.read_text())
sha=lambda b:hashlib.sha256(b).hexdigest()
write=lambda p,v:p.write_text(json.dumps(v,indent=2)+'\n')
start=time.monotonic()
def timeout(*_):raise TimeoutError('120s staging cap')
signal.signal(signal.SIGALRM,timeout);signal.alarm(120)
assert dst.parent.resolve()==Path('/home/lpaiu/work/native-retained-tail-runtime') and not dst.exists()
dst.mkdir();copies={}
try:
    for p in sorted(src.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(src);data=p.read_bytes();target=dst/rel;target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(data);h=sha(data);assert sha(target.read_bytes())==h;copies[str(rel)]=h
    p=Path('outputs/direct-eos-gr33/native-retained-completion/evolution/coupled-128.npz');input_path=p;data=p.read_bytes()
    h=sha(data);assert h==read(src/'plan.json')['bindings'][p.as_posix()]
    target=dst/'immutable-coupled-128.npz';target.write_bytes(data);assert sha(target.read_bytes())==h
    old=read(src/'execution-plan.json');spent=sum(read(src/n)['seconds'] for n in ['collision-receipt.json','collision_resume-receipt.json','hydro_pilot-receipt.json','bank-receipt.json'])
    plan=dict(old,first_attempt_seconds=spent+time.monotonic()-start,remaining_budgets=dict(bank=120,photon_pilot=90,photon_production=900,hydro_pilot_resume=150,hydro_production=800,material_pilot=60,material_production=780,readout=180),
        source_sha256=sha(Path('verification/def_retained_native_return.py').read_bytes()),
        repair='Retain three timed-out attempts and all completed native/coefficient states. Use ext4 intermediate files and a SHA-verified immutable coupled-history alias. Reuse coefficient knots0..3; recheck sampled jets at knot0 without replacing its file. Save accepted angular arrays in both checkpoint and final NPZ. No photon physical step has yet run. Native parallelism, numerical gates and aggregate2880s wall/charged CPU ceiling remain as registered.',
        staging_seconds=time.monotonic()-start,staging_cap_seconds=120,previous_execution_plan_sha256=sha((src/'second-execution-plan.json').read_bytes()))
    for p in [src/'executed-parallel-plan.py',src/'second-execution-plan.json',src/'bank-failure.json',src/'bank-receipt.json',*sorted((src/'photons/bank-128').glob('point-*.npz'))]:
        plan['bindings'][p.as_posix()]=sha(p.read_bytes())
    write(dst/'execution-plan.json',plan)
    result=dict(passed=True,copies=copies,immutable_input=dict(path=str(input_path),sha256=h),seconds=time.monotonic()-start,
        first_attempt_seconds=plan['first_attempt_seconds'],source_sha256=plan['source_sha256'],work_directory=str(dst))
    write(dst/'staging-receipt.json',result);write(src/'staging-receipt.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='copies'}),flush=True)
finally:signal.alarm(0)
