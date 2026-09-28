"""Compare the failed background-energy map to exact-input decimal arithmetic."""
from decimal import Decimal,localcontext
from pathlib import Path
from types import FunctionType
import json,resource,time
import numpy as np
import read_complete_radau_history as r

out=Path('native-complete-radau224-check-work');r.OUT=out;start=time.monotonic()
resource.setrlimit(resource.RLIMIT_AS,(6*1024**3,6*1024**3));r.endpoint.evolution.joint.previous.original.inf.incident.native.deadline(300)
FunctionType(r.endpoint.initialize.__code__,dict(r.endpoint.initialize.__globals__,OUT=out))()
m=r.base.run.owner.Model(64);model=r.base.gr.Response();d=dict(np.load(out/'gr/source-64.npz'));material=m.material
p=np.load(out/'input/material-64.npz');step=14;u=r.LD(5)/6;t=p['actual_step_edges'][step]+u*(p['actual_step_edges'][step+1]-p['actual_step_edges'][step])
cf=model.coeff(material.rE);E=cf['Eg']*r.C**4/r.base.gr.base.G*material.V;P=cf['Pg']*r.C**4/r.base.gr.base.G*material.V
k=int(np.clip(np.searchsorted(m.t,t,side='right')-1,0,15));w=(t-m.t[k])/(m.t[k+1]-m.t[k]);q0=material.point(k)['Q'];q1=material.point(k+1)['Q']
Q=(1-w)*q0+w*q1;old=E+P-(Q[2]+material.rest*Q[0])/material.a
stable=(1-w)*(E+P-(q0[2]+material.rest*q0[0])/material.a)+w*(E+P-(q1[2]+material.rest*q1[0])/material.a)
def dec(v):
    a,b=np.longdouble(v).as_integer_ratio();return Decimal(a)/Decimal(b)
rows=[];refs=[]
for precision in [50,80]:
    with localcontext() as ctx:
        ctx.prec=precision;ww=dec(w);rest=dec(material.rest);ref=[]
        for i in range(m.n):
            val=dec(E[i])+dec(P[i])-((1-ww)*(dec(q0[2,i])+rest*dec(q0[0,i]))+ww*(dec(q1[2,i])+rest*dec(q1[0,i])))/dec(material.a[i])
            ref.append(r.LD(str(val)))
    refs.append(np.array(ref))
vol=(3*m.geometry(float(t))[0][0]+m.geometry(float(t))[0][2])*r.AMP
state=r.prior.prior.evaluate(d['state_coeff_gas_nonrest_energy_erg'][:,step],u)
reference=state+refs[-1]*vol;norm=max(np.sum(abs(reference)),r.LD('1e-290'))
row=dict(classification='Counterexample candidate',step=step,theta=float(u),time=float(t),
    old_error_over_readout=float(np.sum(abs((old-refs[-1])*vol))/norm),stable_error_over_readout=float(np.sum(abs((stable-refs[-1])*vol))/norm),
    decimal50_80_exact=bool(np.array_equal(*refs)),max_contrast_old_abs=float(np.max(abs(old-refs[-1]))),max_contrast_stable_abs=float(np.max(abs(stable-refs[-1]))),
    seconds=time.monotonic()-start,source_sha256=r.sha(__file__))
r.write(out/'energy-probe.json',row);print(json.dumps(row));assert row['decimal50_80_exact']
