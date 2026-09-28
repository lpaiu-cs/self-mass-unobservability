import sys,json,signal,time,traceback
import numpy as np
sys.path.insert(0,'verification')
import def_native_two_way_atmosphere as task
p=task.OUT
plan=dict(classification='Counterexample candidate',claim='Locate the first failed conservative primitive using only saved complete step104 plus at most5 original macrosteps.',decision='Distinguish an EOS support exit from an inadmissible trial time step or inverse solver bug before any further production.',wall_seconds=45,expected_seconds=[18,35],forecast='Measured896 production361.89s for108 completed steps; five replay steps plus setup is expected18-35s, late stiffness remains uncertain.',maximum_macrosteps=5,extra_native_states=0,unchanged='Same grids, horizon, equations, gates and original frozen producer. No full path replay. No final charge inference.',source_sha256=task.sha(task.__file__))
task.write(p/'root-replay-plan.json',plan);signal.signal(signal.SIGALRM,task.optical.timeout);signal.alarm(45);start=time.monotonic();m=task.Coupled(896,8);f=m.flow;d=np.load(p/'cells-896-steps-128.npz');h=task.END/128
U=d['snapshot_U'][-1].copy();I=d['snapshot_I'][-1].copy();xb=d['snapshot_bulk_I'][-1]/m.bulk.scale;theta=d['snapshot_theta'][-1].copy();eta=d['snapshot_eta'][-1].copy();u=m.bulk.eos.gas(theta,eta)[1]
j0=round(float(d['snapshot_t'][-1])/h);assert j0==104
original=f.primitive.__globals__['brentq'];failures=[]
def traced(fn,a,b,**kw):
    try:return original(fn,a,b,**kw)
    except ValueError:
        cl=dict(zip(fn.__code__.co_freevars,[v.cell_contents for v in fn.__closure__]));i=int(cl['i']);bad=cl['U'];rho=cl['rho']
        row=dict(cell=i,x_cm=float(f.base.x[i]),U=bad[:,i].tolist(),rho=float(rho[i]),floor=f.eos.floor,low_T=float(np.exp(a)),high_T=float(np.exp(b)),f_low=float(fn(a)),f_high=float(fn(b)),time_step=j+1)
        failures.append(row);task.write(p/'root-failure-detail.json',row);np.savez_compressed(p/'root-failure-state.npz',U=bad,seed=f.seed);print(json.dumps(row),flush=True);raise
f.primitive.__globals__['brentq']=traced
completed=0;failure=None
try:
    f.primitive(U)
    for j in range(j0,j0+5):
        U,I,*_=m.local(U,I,j*h,h/2)
        xb,I,u,theta,eta,*_=m.radiate(xb,I,u,theta,eta,h)
        U,I,*_=m.local(U,I,j*h+h/2,h/2);f.primitive(U);completed+=1
        print(json.dumps(dict(step=j+1,seconds=time.monotonic()-start)),flush=True)
except Exception as e:failure=repr(e);traceback.print_exc()
task.write(p/'root-replay.json',dict(classification='Counterexample candidate',failure=failure,completed_replayed_macros=completed,seconds=time.monotonic()-start,failures=failures,native_constructor_calls=m.spectrum.native.ion.calls));signal.alarm(0)
