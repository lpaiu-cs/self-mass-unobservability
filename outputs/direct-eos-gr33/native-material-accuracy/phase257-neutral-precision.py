from pathlib import Path
from decimal import Decimal,localcontext
import json,time,numpy as np
p=Path('native-returned-right256-work/sweep-1/photons/return-128.npz');z=np.load(p)
step=219;h=np.longdouble(4)*z['joint_stage_weights'][2*step+1]
y=z['joint_stage_conserved_scaled'][2*step:2*step+2,3];x=z['joint_stage_conserved_scaled'][2*step-1,3];unit=z['material_neutral_units']
n=z['joint_native_rates_scaled'][2*step:2*step+2,:,1];c=z['joint_collision_rates_scaled'][2*step:2*step+2,:,1]
weights=np.array([np.longdouble('.75'),np.longdouble('.25')]);rates=(n+c)/unit
raw=(y[1]/unit-x/unit-h*np.einsum('j,jn->n',weights,rates))*unit
physical=y[1]-x-h*np.einsum('j,jn->n',weights,n+c)
def dec(v):
 a,b=v.as_integer_ratio();return Decimal(a)/Decimal(b)
start=time.monotonic()
with localcontext() as ctx:
 ctx.prec=80;dh=dec(h);dw=[dec(v) for v in weights]
 exact=np.array([np.longdouble(str(dec(y[1,i])-dec(x[i])-dh*sum(dw[j]*(dec(n[j,i])+dec(c[j,i])) for j in range(2)))) for i in range(len(x))])
norm=np.sum(abs(y[1]),dtype=np.longdouble)
result=dict(classification='Counterexample candidate',actual_step=220,gas_component='H',original_reader_relative=float(np.sum(abs(raw))/norm),direct_conserved_relative=float(np.sum(abs(physical))/norm),exact80_conserved_relative=float(np.sum(abs(exact))/norm),arithmetic_effect=float(np.sum(abs(exact-raw))/norm),precision_seconds=time.monotonic()-start,initial_floor_not_applied_in_this_diagnostic=True,source_gate=1e-12,new_physical_steps=0)
print(json.dumps(result));Path('phase257-neutral-precision.json').write_text(json.dumps(result,indent=2)+'\n')
