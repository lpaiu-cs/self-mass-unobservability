import sys,json,signal,time,traceback
import numpy as np
sys.path.insert(0,'verification')
import def_native_two_way_atmosphere as task
p=task.OUT;signal.signal(signal.SIGALRM,task.optical.timeout);signal.alarm(20);start=time.monotonic();m=task.Coupled(896,8);f=m.flow;d=np.load(p/'cells-896-steps-128.npz');h=task.END/128
U=d['U'].copy();I=d['I'].copy();original=f.primitive.__globals__['brentq'];failures=[]
def traced(fn,a,b,**kw):
    try:return original(fn,a,b,**kw)
    except ValueError:
        cl=dict(zip(fn.__code__.co_freevars,[v.cell_contents for v in fn.__closure__]));i=int(cl['i']);bad=cl['U']
        np.savez_compressed(p/'root-failure-state.npz',U=bad,seed=f.seed)
        row=dict(cell=i,x_cm=float(f.base.x[i]),U=bad[:,i].tolist(),floor=f.eos.floor,low_T=float(np.exp(a)),high_T=float(np.exp(b)),f_low=float(fn(a)),f_high=float(fn(b)),time_step=109)
        failures.append(row);task.write(p/'root-failure-detail.json',row);print(json.dumps(row),flush=True);raise
f.primitive.__globals__['brentq']=traced;failure=None
try:U,I,*_=m.local(U,I,108*h,h/2)
except Exception as e:failure=repr(e);traceback.print_exc()
task.write(p/'terminal-halfstep-diagnostic.json',dict(classification='Counterexample candidate',failure=failure,seconds=time.monotonic()-start,failures=failures,native_constructor_calls=m.spectrum.native.ion.calls,scope='Diagnosis only. The failed terminal radiation may include partial local updates; this is not a resumed accepted trajectory.',prior_diagnostic_error='The exact replay reproduced the109th macrostep failure, but diagnostic closure extraction assumed a nonexistent rho capture. This attempt saves U before reading other fields.'));signal.alarm(0)
