"""Prepare, check (both clocks in parallel), finish and read the photon-boundary pair with the long-double fallback."""
from pathlib import Path
import hashlib,json,os,subprocess,sys,time

status=Path('native-photon-boundary-ld265-chain.json')
producer=Path('verification/apply_photon_geometric_boundary_ld.py');readout=Path('verification/read_photon_boundary_exterior.py')
solver=Path('verification/solve_long_double_fgmres.py');launcher=Path('.phase251-reader-launch.py');follow=Path('.phase265-ld-followthrough.py')
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONPATH='/home/lpaiu/work/nutimo_pilot/request13_deps:verification')
read=lambda p:json.loads(Path(p).read_text())
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v):
    tmp=Path(str(p)+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');os.replace(tmp,p)
bound={str(p):sha(p) for p in [Path(__file__),producer,readout,solver,launcher,follow]}
start=time.monotonic();done=[]
def launch(name,command,cpu):
    for p,h in bound.items():assert sha(p)==h,p
    o=open(f'.phase265-ld-chain-{name}.stdout.log','w');e=open(f'.phase265-ld-chain-{name}.stderr.log','w')
    return subprocess.Popen(['taskset','-c',str(cpu),sys.executable,*command],stdout=o,stderr=e,env=env)
def steps(group):
    write(status,dict(state='running',steps=[g[0] for g in group],completed=done,seconds=time.monotonic()-start,bindings=bound))
    children=[(name,launch(name,command,cpu)) for name,command,cpu in group]
    for name,child in children:
        code=child.wait();done.append(dict(step=name,returncode=code,seconds=time.monotonic()-start))
    for name,child in children:assert child.returncode==0,(name,child.returncode)
try:
    steps([('prepare',[str(launcher),str(producer),'prepare'],6)])
    steps([('boundary',[str(launcher),str(producer),'boundary'],6)])
    steps([('check64',[str(launcher),str(producer),'check64'],4),('check128',[str(launcher),str(producer),'check128'],6)])
    steps([('followthrough',[str(follow)],2)])
    steps([('complete_readout',[str(readout),'full'],2)])
    result=read('native-photon-boundary-ld-exterior265-work/complete/result.json')
    write(status,dict(state='completed',completed=done,seconds=time.monotonic()-start,bindings=bound,
        final_total_clock128=result['final_total_clock128'],final_sign_negative=result['final_sign_negative'],readout_passed=result['passed']))
except BaseException as exc:
    write(status,dict(state='failed',completed=done,error=repr(exc),seconds=time.monotonic()-start,bindings=bound));raise
