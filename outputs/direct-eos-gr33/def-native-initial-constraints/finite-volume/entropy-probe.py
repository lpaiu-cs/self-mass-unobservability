import json,time
import numpy as np
import def_native_initial_constraints as t
start=time.monotonic();d=t.FiniteVolumeData().model.bulk.d;n=t.previous.chem.old.Native(cap=20)
z=np.load(t.OUT/'finite-volume/balanced-initial-state.npz');j=1;t.previous.chem.setup(n,d,j)
lt0=np.log(d['T'][j]);x=float(z['density_log_ratio'][j]);y=d['y0'][j];target=d['raw'][j,3]
base=n.state(0.,lt0,y)['raw'];rows=[];lt=lt0
for i in range(10):
    raw=n.state(x,lt,y)['raw'];delta=(raw[3]-target)*np.exp(lt)/raw[10]
    rows.append(dict(i=i,logT=lt,delta=float(delta),relative_to_ref=float((raw[3]-base[3])*np.exp(lt)/raw[10])))
    lt-=delta
t.write(t.OUT/'finite-volume/entropy-probe.json',dict(native_calls=n.ion.calls,reference_entropy_difference=float((base[3]-target)*np.exp(lt0)/base[10]),rows=rows,seconds=time.monotonic()-start))
print(json.dumps(rows))
