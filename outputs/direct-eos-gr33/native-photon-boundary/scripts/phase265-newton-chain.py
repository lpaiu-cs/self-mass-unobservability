"""Retry after the first attempt's Newton stop: build, check, evolve and read the photon-boundary pair with12Newton proposals."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

geometric=Path('native-geometric-clock265-status.json');status=Path('native-photon-boundary-newton265-chain.json')
producer=Path('verification/apply_photon_geometric_boundary_newton.py');readout=Path('verification/read_photon_boundary_exterior.py')
launcher=Path('.phase251-reader-launch.py');follow=Path('.phase265-newton-followthrough.py')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH='/home/lpaiu/work/nutimo_pilot/request13_deps:verification')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
bound={str(p):sha(p) for p in [Path(__file__),producer,readout,launcher,follow]}
start=time.monotonic();done=[]
def step(name,command,cpu):
    for p,h in bound.items():assert sha(p)==h,p
    write(status,dict(state='running',step=name,completed=done,seconds=time.monotonic()-start,bindings=bound))
    with open(f'.phase265-newton-chain-{name}.stdout.log','w') as o,open(f'.phase265-newton-chain-{name}.stderr.log','w') as e:
        code=subprocess.run(['taskset','-c',str(cpu),sys.executable,*command],stdout=o,stderr=e,env=env).returncode
    done.append(dict(step=name,returncode=code,seconds=time.monotonic()-start));assert code==0,(name,code)
try:
    while True:
        s=read(geometric)
        if s['state']=='completed':break
        assert s['state']!='failed',s
        write(status,dict(state='waiting',step='geometric_boundary',pending_items=s.get('pending_items'),seconds=time.monotonic()-start,bindings=bound))
        time.sleep(30)
    assert read('native-geometric-clock265-work/result.json')['passed']
    step('prepare',[str(launcher),str(producer),'prepare'],6)
    step('check',[str(launcher),str(producer),'check'],6)
    step('followthrough',[str(follow)],2)
    step('complete_readout',[str(readout),'full'],2)
    result=read('native-photon-boundary-newton-exterior265-work/complete/result.json')
    write(status,dict(state='completed',completed=done,seconds=time.monotonic()-start,
        final_total_clock128=result['final_total_clock128'],final_sign_negative=result['final_sign_negative'],readout_passed=result['passed']))
except BaseException as exc:
    write(status,dict(state='failed',completed=done,error=repr(exc),seconds=time.monotonic()-start));raise
